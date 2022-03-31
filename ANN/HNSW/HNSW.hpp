#ifndef _HNSW_HPP
#define _HNSW_HPP

#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <random>
#include <memory>
#include <vector>
#include <unordered_map>
#include <queue>
#include <iterator>
#include <type_traits>
#include <limits>
// #include "parallelize.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#define DEBUG_OUTPUT 0

namespace ANN{

enum class type_metric{
	L2, ANGULAR, DOT
};

struct point{
	double x, y;
};

template<typename U, template<typename> class Allocator=std::allocator>
class HNSW
{
	using T = typename U::type_point;
public:
	template<typename Iter>
	HNSW(Iter begin, Iter end, uint32_t dim, float m_l=1, uint32_t m=100, uint32_t ef_construction=50, float alpha=5, float batch_base=2);
	std::vector<T*> search(const T &q, uint32_t k, uint32_t ef);
private:
	typedef uint32_t type_index;

	struct node{
		// uint32_t id;
		uint32_t level;
		T data;
		std::vector<node*> *neighbors;
	};

	struct dist{
		double d;
		node *u;
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


	// uint32_t cnt_level;
	// std::vector<double> probability; // To purge
	// std::vector<uint32_t> cnt_neighbor;
	node *entrance; // To init
	// auto m, max_m0, m_L; // To init
	uint32_t dim;
	float m_l;
	uint32_t m;
	// uint32_t level_max = 30; // To init
	uint32_t ef_construction;
	float alpha;
	uint32_t n;
	Allocator<node> allocator;
	std::vector<node*> node_pool;
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
	static auto neighbourhood(const node &u, uint32_t level)
		-> std::vector<node*>&
	{
		return u.neighbors[level];
	}

	// `set_neighbourhood` will consume `vNewConn`
	static void set_neighbourhood(node &u, uint32_t level, std::vector<node*>& vNewConn)
	{
		u.neighbors[level] = std::move(vNewConn);
	}

	static void add_connection(std::vector<node*> &neighbors, node &u, uint32_t level)
	{
		for(auto pv : neighbors)
		{
			assert(&u!=pv);
			pv->neighbors[level].push_back(&u);
			u.neighbors[level].push_back(pv);
		}
	}

	// node* insert(const T &q, uint32_t id);
	template<typename Iter>
	void insert(Iter begin, Iter end);

	void select_neighbors_simple_impl(const T &u, 
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

	auto select_neighbors_simple(const T &u, 
		std::priority_queue<dist,std::vector<dist>,farthest> C, uint32_t M)
	{
		// The parameter C is intended to be copy constructed
		select_neighbors_simple_impl(u, C, M);
		return C;
	}

	// To optimize
	auto select_neighbors_heuristic(const T &u, 
		const std::priority_queue<dist,std::vector<dist>,farthest> &C, uint32_t M,
		uint32_t level, bool extendCandidate, bool keepPrunedConnections)
	{
		(void)level, (void)extendCandidate;

		std::priority_queue<dist,std::vector<dist>,farthest> C_cp=C, W_d;
		/*
		if(extendCandidate)
		{
			while(C_cp.size())
			{
				const auto e = C_cp.top();
				C_cp.pop();
				for(const auto &e_adj : neighbourhood(e.u,level))
				{
					if(!W.find(e_adj))
						W.push(e_adj);
				}
		}
		*/
		std::vector<node*> R;
		std::priority_queue<dist,std::vector<dist>,nearest> W;
		while(C_cp.size())
		{
			W.push(C_cp.top());
			C_cp.pop();
		}
		//auto &W = C_cp;

		while(W.size()>0 && R.size()<M)
		{
			const auto e = W.top();
			W.pop();
			const auto d_q = e.d;

			bool is_good = true;
			for(const auto &r : R)
			{
				const auto d_r = U::distance(e.u->data, r->data, dim);
				if(d_r>d_q*alpha)
				{
					is_good = false;
					break;
				}
			}

			if(is_good)
				R.push_back(e.u);
			else
				W_d.push(e);
		}

		std::priority_queue<dist,std::vector<dist>,farthest> res;
		for(const auto &r : R)
		{
			res.push({U::distance(u,r->data,dim), r});
		}
		if(keepPrunedConnections)
		{
			while(W_d.size()>0 && res.size()<M)
				res.push(W_d.top()), W_d.pop();
		}
		return res;
	}

	auto select_neighbors(const T &u, 
		const std::priority_queue<dist,std::vector<dist>,farthest> &C, uint32_t M,
		uint32_t level, bool extendCandidate=false, bool keepPrunedConnections=false)
	{
		/*
		(void)level, (void)extendCandidate, (void)keepPrunedConnections;
		return select_neighbors_simple(u,C,M);
		*/
		return select_neighbors_heuristic(u, C, M, level, extendCandidate, keepPrunedConnections);
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

	auto search_layer(const node &u, const std::vector<node*> &eps, uint32_t ef, uint32_t l_c) const; // To static
	auto get_threshold_m(uint32_t level){
		return level==0? m*2: m;
	}
};


template<typename T, template<typename> class Allocator>
template<typename Iter>
HNSW<T,Allocator>::HNSW(Iter begin, Iter end, uint32_t dim_, float m_l_, uint32_t m_, uint32_t ef_construction_, float alpha_, float batch_base)
	: dim(dim_), m_l(m_l_), m(m_), ef_construction(ef_construction_), alpha(alpha_), n(std::distance(begin,end))
{
	static_assert(std::is_same_v<typename std::iterator_traits<Iter>::value_type, T>);
	static_assert(std::is_base_of_v<
		std::random_access_iterator_tag, typename std::iterator_traits<Iter>::iterator_category>);

	if(n==0) return;

	const auto level_ep = get_level_random();
	entrance = allocator.allocate(1);
	new(entrance) node{level_ep, *begin, new std::vector<node*>[level_ep+1]/*anything else*/};
	#if DEBUG_OUTPUT
		fprintf(stderr, "[%u] at lv.%u (%.2f,%.2f)**\n", 0, level_ep, begin->x, begin->y);
	#endif
	node_pool.push_back(entrance);

	uint32_t batch_begin=0, batch_end=1;
	while(batch_end<n)
	{
		batch_begin = batch_end;
		batch_end = std::min(n, (uint32_t)std::ceil(batch_begin*batch_base));
		// batch_end = batch_begin+1;

		insert(begin+batch_begin, begin+batch_end);

	}
/*
	for(uint32_t i=1; i<n; ++i)
	{
		auto *p = insert(*(begin+i), i);
		node_pool.push_back(p);
	}
*/
	#if DEBUG_OUTPUT
		for(const auto *pu : node_pool)
		{
			fprintf(stderr, "[%u] (%.2f,%.2f)\n", U::get_id(pu->data), pu->data.x, pu->data.y);
			for(int32_t l=pu->level; l>=0; --l)
			{
				fprintf(stderr, "\tlv. %d:", l);
				for(const auto *k : pu->neighbors[l])
					fprintf(stderr, " %u", U::get_id(k->data));
				fputs("\n", stderr);
			}
		}
	#endif
}

template<typename U, template<typename> class Allocator>
// typename HNSW<U,Allocator>::node* HNSW<U,Allocator>::insert(const T &q, uint32_t id)
template<typename Iter>
void HNSW<U,Allocator>::insert(Iter begin, Iter end)
{
	const auto level_ep = entrance->level;
	const auto size_batch = std::distance(begin,end);
	auto node_new = std::make_unique<node*[]>(size_batch);
	auto eps = std::make_unique<std::vector<node*>[]>(size_batch);

	// first, query the nearest point as the starting point for each node to insert
	parlay::parallel_for(0, size_batch, [&](uint32_t i){
		const T &q = *(begin+i);
		auto &eps_u = eps[i]; 
		eps_u.push_back(entrance);
		const auto level_u = get_level_random();
		auto *const pu = allocator.allocate(1);		// TODO: add pointer manager

		auto &u = *new(pu) node{level_u, q, new std::vector<node*>[level_u+1]};
		#if DEBUG_OUTPUT
			fprintf(stderr, "[%u] at lv.%u (%.2f,%.2f)\n", U::get_id(q), level_u, q.x, q.y);
		#endif
		node_new[i] = pu;

		for(uint32_t l=level_ep; l>level_u; --l)
		{
			const auto res = search_layer(u, eps_u, 1, l); // TODO: optimize
			eps_u[0] = res.top().u;
		}
	});

	// then we process them layer by layer (from high to low)
	for(int32_t l_c=level_ep; l_c>=0; --l_c) // TODO: fix the type
	{
		parlay::sequence<parlay::sequence<std::pair<node*,node*>>> edge_add(size_batch);

		parlay::parallel_for(0, size_batch, [&](uint32_t i){
			auto &u = *node_new[i];
			if((uint32_t)l_c>u.level) return;

			auto &eps_u = eps[i];
			auto res = search_layer(u, eps_u, ef_construction, l_c);
			auto neighbors_queue = select_neighbors(u.data, res, m, l_c);
			// move the content from `neighbors_queue` to `u.neighbors[l_c]`
			auto &nbh_u = neighbourhood(u, l_c);
			auto &edge_u = edge_add[i];
			nbh_u.resize(neighbors_queue.size());
			edge_u.resize(neighbors_queue.size());
			for(uint32_t j=0; neighbors_queue.size()>0; ++j)
			{
				auto *pv = neighbors_queue.top().u;
				neighbors_queue.pop();
				nbh_u[j] = pv;
				edge_u[j] = std::make_pair(pv, &u);
			}

			eps_u.clear();
			while(res.size()>0)
			{
				eps_u.push_back(res.top().u); // TODO: optimize
				res.pop();
			}
		});

		// now we add edges in the other direction
		auto edge_add_flatten = parlay::flatten(edge_add);
		auto edge_add_grouped = parlay::group_by_key(edge_add_flatten);

		parlay::parallel_for(0, edge_add_grouped.size(), [&](size_t j){
			node *pv = edge_add_grouped[j].first;
			const auto &nbh_v_add = edge_add_grouped[j].second;
			auto &nbh_v = neighbourhood(*pv,l_c);
			const uint32_t size_nbh_total = nbh_v.size()+nbh_v_add.size();

			if(size_nbh_total>get_threshold_m(l_c))
			{
				auto dist_nbh = std::make_unique<dist[]>(size_nbh_total);
				for(size_t k=0; k<nbh_v.size(); ++k)
					dist_nbh[k] = dist{U::distance(nbh_v[k]->data,pv->data,dim), nbh_v[k]};
				for(size_t k=0; k<nbh_v_add.size(); ++k)
					dist_nbh[k+nbh_v.size()] = dist{U::distance(nbh_v_add[k]->data,pv->data,dim), nbh_v_add[k]};

				std::sort(dist_nbh.get(), dist_nbh.get()+size_nbh_total, farthest());

				nbh_v.resize(m);
				for(size_t k=0; k<m; ++k)
					nbh_v[k] = dist_nbh[k].u;
			}
			else nbh_v.insert(nbh_v.end(),nbh_v_add.begin(), nbh_v_add.end());
		});
	}

	// finally, update the entrance
	auto *node_highest = *std::max_element(
		node_new.get(), node_new.get()+size_batch, [](const node *u, const node *v){
			return u->level < v->level;
	});
	if(node_highest->level>level_ep)
	{
		entrance = node_highest;
		// #if DEBUG_OUTPUT
			fprintf(stderr, "[%u]** at lev %u\n", U::get_id(entrance->data), entrance->level);
		// #endif
	}

	// and add new nodes to the pool
	node_pool.insert(node_pool.end(), node_new.get(), node_new.get()+size_batch);
}

template<typename U, template<typename> class Allocator>
auto HNSW<U,Allocator>::search_layer(const node &u, const std::vector<node*> &eps, uint32_t ef, uint32_t l_c) const
{
	std::vector<bool> visited(n);
	std::priority_queue<dist,std::vector<dist>,nearest> C;
	std::priority_queue<dist,std::vector<dist>,farthest> W;

	for(auto *ep : eps)
	{
		visited[U::get_id(ep->data)] = true;
		const auto d = U::distance(u.data,ep->data,dim);
		C.push({d,ep});
		W.push({d,ep});
	}

	while(C.size()>0)
	{
		const auto &c = *C.top().u; C.pop();
		const auto &f = *W.top().u;
		if(U::distance(c.data,u.data,dim)>U::distance(f.data,u.data,dim))
			break;
		for(auto *pv: neighbourhood(c, l_c))
		{
			if(visited[U::get_id(pv->data)]) continue;
			visited[U::get_id(pv->data)] = true;
			const auto &f = *W.top().u;
			if(U::distance(pv->data,u.data,dim)<U::distance(f.data,u.data,dim)||W.size()<ef)
			{
				const auto d = U::distance(u.data,pv->data,dim);
				C.push({d,pv});
				W.push({d,pv});
				if(W.size()>ef) W.pop();
			}
		}
	}
	return W;
}

template<typename U, template<typename> class Allocator>
std::vector<typename HNSW<U,Allocator>::T*> HNSW<U,Allocator>::search(const T &q, uint32_t k, uint32_t ef)
{
	node u{0, q, nullptr}; // To optimize
	std::priority_queue<dist,std::vector<dist>,farthest> W;
	auto *ep = entrance;
	for(int l_c=entrance->level; l_c>0; --l_c) // TODO: fix the type
	{
		W = search_layer(u, std::vector{ep}, 1, l_c);
		ep = W.top().u;
	}
	W = search_layer(u, std::vector{ep}, ef, 0);
	W = select_neighbors_simple(q, W, k);
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
