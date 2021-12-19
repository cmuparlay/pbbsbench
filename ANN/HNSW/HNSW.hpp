#ifndef _HNSW_HPP
#define _HNSW_HPP

#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <unordered_map>
#include <queue>
#include <iterator>
#include <type_traits>
#include <limits>
// #include "parallelize.h"
#define DEBUG_OUTPUT 1

namespace ANN{

enum class type_metric{
	L2, ANGULAR, DOT
};

struct point{
	double x, y;
};

template<typename T>
struct node{
	uint32_t id;
	uint32_t level;
	T data;
	std::vector<node*> *neighbors;
};

struct dist{
	double d;
	node<point> *u;
/*
	dist(double d_, node<point> *u_) :
		d(d_), u(u_)
	{
	}

	dist& operator=(const dist &other)	// To rvalue
	{
		d = other.d;
		u = other.u;
		return *this;
	}
*/
};

struct nearest{
	constexpr bool operator()(const dist &lhs, const dist &rhs) const{
		return lhs.d>rhs.d;
	}
};

struct farthest{
	constexpr bool operator()(const dist &lhs, const dist &rhs) const{
		return lhs.d<rhs.d;
	}
};


template<typename T>
inline double distance(const node<T> &u, const node<T> &v)
{
	const auto &pu=u.data, &pv=v.data;
	return (pu.x-pv.x)*(pu.x-pv.x)+(pu.y-pv.y)*(pu.y-pv.y);
}


template<typename T>
inline auto neighbourhood(const node<T> &u, uint32_t level)
	-> std::vector<node<T>*>&
{
	return u.neighbors[level];
}

// `set_neighbourhood` will consume `vNewConn`
template<typename T>
void set_neighbourhood(node<T> &u, uint32_t level, std::vector<node<T>*>& vNewConn)
{
	u.neighbors[level] = std::move(vNewConn);
}

template<typename T>
void add_connection(std::vector<node<T>*> &neighbors, node<T> &u, uint32_t level)
{
	for(auto pv : neighbors)
	{
		assert(&u!=pv);
		pv->neighbors[level].push_back(&u);
		u.neighbors[level].push_back(pv);
	}
}

template<typename T, template<typename> class Allocator=std::allocator>
class HNSW
{
public:
	template<typename Iter>
	HNSW(Iter begin, Iter end);
	std::vector<T*> search(const T &q, uint32_t k, uint32_t ef);
private:
	typedef uint32_t type_index;

	// uint32_t cnt_level;
	// std::vector<double> probability; // To purge
	// std::vector<uint32_t> cnt_neighbor;
	node<T> *entrance; // To init
	// auto m, max_m0, m_L; // To init
	uint32_t m_l = 15;
	uint32_t m = 1000;
	// uint32_t level_max = 30; // To init
	uint32_t ef_construction = 120;
	uint32_t n = 0;
	Allocator<node<T>> allocator;
	std::vector<node<T>*> node_pool;
/*
	void set_probability(uint32_t step, double multiplier_level, double eps_prob)
	{
		uint32_t n = 0;
		// traverse the levels
		for(uint32_t i=0;; ++i)
		{
			const auto p = exp(-i/multiplier_level)*(1-exp(-1/multiplier_level));
			if(p<esp) break;
			probability.push_back(p);
			n += step;
			cnt_neighbor.push_back(n);
		}
	}
*/
	node<T>* insert(const T &q, uint32_t id);

	void select_neighbors_simple_impl(const point &u, 
		std::priority_queue<dist,std::vector<dist>,farthest> &C, uint32_t M)
	{
		/*
		list res;
		for(uint32_t i=0; i<M; ++i)
		{
			res.insert(C.pop_front());
		}
		return res;
		*/
		(void)u;
		while(C.size()>M) C.pop();
	}

	auto select_neighbors_simple(const point &u, 
		std::priority_queue<dist,std::vector<dist>,farthest> C, uint32_t M)
	{
		// The parameter C is intended to be copy constructed
		select_neighbors_simple_impl(u, C, M);
		return C;
	}

	// To optimize
	auto select_neighbors(const point &u, 
		const std::priority_queue<dist,std::vector<dist>,farthest> &C, uint32_t M,
		uint32_t level, bool extendCandidate=false, bool keepPrunedConnections=false)
	{
		(void)level, (void)extendCandidate, (void)keepPrunedConnections;
		return select_neighbors_simple(u,C,M);
	}

	uint32_t get_level_random()
	{
		static thread_local int32_t anchor;
		// uint32_t esp;
		// asm volatile("movl %0, %%esp":"=a"(esp));
		static thread_local std::mt19937 gen{anchor};
		static thread_local std::uniform_real_distribution<> dis(std::numeric_limits<double>::min(), 1.0);
		const uint32_t res = uint32_t(-log(dis(gen))*m_l);
		return res;
	}

	auto search_layer(const node<T> &u, const std::vector<node<T>*> &eps, uint32_t ef, uint32_t l_c) const; // To static
};


template<typename T, template<typename> class Allocator>
template<typename Iter>
HNSW<T,Allocator>::HNSW(Iter begin, Iter end)
{
	static_assert(std::is_same_v<typename std::iterator_traits<Iter>::value_type, T>);
	static_assert(std::is_base_of_v<
		std::random_access_iterator_tag, typename std::iterator_traits<Iter>::iterator_category>);

	if(begin==end) return;

	const auto level_ep = get_level_random();
	entrance = allocator.allocate(1);
	new(entrance) node<T>{0, level_ep, *begin, new std::vector<node<T>*>[level_ep+1]/*anything else*/};
	#if DEBUG_OUTPUT
		fprintf(stderr, "[%u] at lv.%u (%.2f,%.2f)**\n", 0, level_ep, begin->x, begin->y);
	#endif
	node_pool.push_back(entrance);
	
	n = std::distance(begin, end);
	for(uint32_t i=1; i<n; ++i)
	{
		auto *p = insert(*(begin+i), i);
		node_pool.push_back(p);
	}

	#if DEBUG_OUTPUT
		for(const auto *pu : node_pool)
		{
			fprintf(stderr, "[%u] (%.2f,%.2f)\n", pu->id, pu->data.x, pu->data.y);
			for(int32_t l=pu->level; l>=0; --l)
			{
				fprintf(stderr, "\tlv. %d:", l);
				for(const auto *k : pu->neighbors[l])
					fprintf(stderr, " %u", k->id);
				fputs("\n", stderr);
			}
		}
	#endif
}

template<typename T, template<typename> class Allocator>
node<T>* HNSW<T,Allocator>::insert(const T &q, uint32_t id)
{
	std::vector<node<T>*> eps = {entrance};
	const auto level_ep = entrance->level;
	const auto level_u = get_level_random();
	auto *const pu = allocator.allocate(1);		// To add pointer manager
	auto &u = *new(pu) node<T>{id, level_u, q, new std::vector<node<T>*>[level_u+1]};
	#if DEBUG_OUTPUT
		fprintf(stderr, "[%u] at lv.%u (%.2f,%.2f)\n", id, level_u, q.x, q.y);
	#endif

	for(uint32_t l=level_ep; l>level_u; --l)
	{
		const auto res = search_layer(u, eps, 1, l); // To optimize
		eps[0] = res.top().u;
	}

	for(int32_t l_c=std::min(level_u,level_ep); l_c>=0; --l_c) // TODO: fix the type
	{
		auto res = search_layer(u, eps, ef_construction, l_c);
		auto neighbors_queue = select_neighbors(q, res, m, l_c);
		std::vector<node<T>*> neighbors;
		while(neighbors_queue.size()>0)
		{
			neighbors.push_back(neighbors_queue.top().u);
			neighbors_queue.pop();
		}
		add_connection(neighbors, u, l_c);
		for(auto *pv: neighbors)	// To be const
		{
			auto &vConn = neighbourhood(*pv, l_c);
			if(vConn.size()>m)
			{
				// if l_c==0 then M_max = M_max0	// To check
				/*
				std::priority_queue<dist,std::vector<dist>,farthest> dist_v;
				for(const auto *pw : vConn)
					dist_v.emplace(distance(*pw,*pv),*pw);
				auto vNewConn = select_neighbors(pv->data, dist_v, m, l_c);
				set_neighbourhood(*pv, l_c, vNewConn);
				*/
				std::sort(vConn.begin(), vConn.end(), [=](const node<T> *lhs, const node<T> *rhs){
					return distance(*lhs,*pv)<distance(*rhs,*pv);
				});
				vConn.resize(m);
			}
		}
		eps.clear();
		while(res.size()>0)
		{
			eps.push_back(res.top().u);
			res.pop();
		}
	}

	if(level_u>level_ep)
	{
		entrance = &u;
		#if DEBUG_OUTPUT
			fprintf(stderr, "[%u]**\n", u.id);
		#endif
	}
	return pu;
}

template<typename T, template<typename> class Allocator>
auto HNSW<T,Allocator>::search_layer(const node<T> &u, const std::vector<node<T>*> &eps, uint32_t ef, uint32_t l_c) const
{
	std::vector<bool> visited(n);
	std::priority_queue<dist,std::vector<dist>,nearest> C;
	std::priority_queue<dist,std::vector<dist>,farthest> W;

	for(auto *ep : eps)
	{
		visited[ep->id] = true;
		const auto d = distance(u,*ep);
		C.push({d,ep});
		W.push({d,ep});
	}

	while(C.size()>0)
	{
		const auto &c = *C.top().u; C.pop();
		const auto &f = *W.top().u;
		if(distance(c,u)>distance(f,u))
			break;
		for(auto *pv: neighbourhood(c, l_c))
		{
			if(visited[pv->id]) continue;
			visited[pv->id] = true;
			const auto &f = *W.top().u;
			if(distance(*pv,u)<distance(f,u)||W.size()<ef)
			{
				const auto d = distance(u,*pv);
				C.push({d,pv});
				W.push({d,pv});
				if(W.size()>ef) W.pop();
			}
		}
	}
	return W;
}

template<typename T, template<typename> class Allocator>
std::vector<T*> HNSW<T,Allocator>::search(const T &q, uint32_t k, uint32_t ef)
{
	node<T> u{0, 0, q, nullptr}; // To optimize
	std::priority_queue<dist,std::vector<dist>,farthest> W;
	auto *ep = entrance;
	for(int l_c=entrance->level; l_c>0; --l_c) // TODO: fix the type
	{
		W = search_layer(u, std::vector{ep}, 1, l_c);
		ep = W.top().u;
	}
	W = search_layer(u, std::vector{ep}, ef, 0);
	W = select_neighbors(q, W, k, 0);
	std::vector<T*> res;
	while(W.size()>0)
	{
		res.push_back(&W.top().u->data);
		W.pop();
	}
	return res;
}

} // namespace HNSW

#endif // _HNSW_HPP