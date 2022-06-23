#ifndef _HNSW_HPP
#define _HNSW_HPP

#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <random>
#include <memory>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <queue>
#include <set>
#include <iterator>
#include <type_traits>
#include <limits>
#include <thread>
// #include "parallelize.h"
#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/delayed_sequence.h>
#include <parlay/random.h>
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
	/*
		Construct from the vectors [begin, end).
		std::iterator_trait<Iter>::value_type ought to be convertible to T
		dim: 				vector dimension
		m_l: 				control the # of levels (larger m_l leads to more layer)
		m: 					max degree
		ef_construction:	beam size during the construction
		alpha:				parameter of the heuristic (similar to the one in vamana)
		batch_base: 		growth rate of the batch size (discarded because of two passes)
	*/
	template<typename Iter>
	HNSW(Iter begin, Iter end, uint32_t dim, float m_l=1, uint32_t m=100, uint32_t ef_construction=50, float alpha=5, float batch_base=2, bool do_fixing=false);

	/*
		Construct from the saved model
		getter(i) returns the actual data (convertible to type T) of the vector with id i
	*/
	template<typename G>
	HNSW(const std::string &filename_model, G getter);

	std::vector<std::pair<uint32_t,double>> search(const T &q, uint32_t k, uint32_t ef);
	std::vector<std::tuple<uint32_t,uint32_t,double>> search_ex(const T &q, uint32_t k, uint32_t ef);
	// save the current model to a file
	void save(const std::string &filename_model) const;
public:
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
	};

	struct dist_ex : dist
	{
		uint32_t depth;
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

	struct cmp_id{
		constexpr bool operator()(const dist &lhs, const dist &rhs) const{
			return U::get_id(lhs.u->data)<U::get_id(rhs.u->data);
		}
	};

	std::vector<node*> entrance; // To init
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
	void insert(Iter begin, Iter end, bool from_blank);

	template<typename Queue>
	void select_neighbors_simple_impl(const T &u, Queue &C, uint32_t M)
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
		std::vector<typename Queue::value_type> tie;
		double dist_tie = 1e20;
		while(C.size()>M)
		{
			const auto &t = C.top();
			if(t.d+1e-6<dist_tie) // t.d<dist_tie
			{
				dist_tie = t.d;
				tie.clear();
			}
			if(fabs(dist_tie-t.d)<1e-6) // t.d==dist_tie
				tie.push_back(t);
			C.pop();
		}
		if(fabs(dist_tie-C.top().d)<1e-6) // C.top().d==dist_tie
			while(!tie.empty())
			{
			//	C.push({dist_tie,tie.back()});
				C.push(tie.back());
				tie.pop_back();
			}
	}

	template<typename Queue>
	auto select_neighbors_simple(const T &u, Queue C, uint32_t M)
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
		(void)extendCandidate;

		std::priority_queue<dist,std::vector<dist>,farthest> C_cp=C, W_d;
		std::set<dist,cmp_id> W_tmp;
		while(!C_cp.empty())
		{
			auto &e = C_cp.top();
			W_tmp.insert(e);
			if(extendCandidate)
			{
				for(auto *e_adj : neighbourhood(*e.u,level))
				{
					if(e_adj==nullptr) continue;
					if(W_tmp.find(dist{0,e_adj})==W_tmp.end())
						W_tmp.insert(dist{U::distance(e_adj->data,u,dim),e_adj});
				}
			}
			C_cp.pop();
		}

		std::priority_queue<dist,std::vector<dist>,nearest> W;
		for(auto &e : W_tmp)
			W.push(e);
		W_tmp.clear();

		std::vector<node*> R;
		while(W.size()>0 && R.size()<M)
		{
			const auto e = W.top();
			W.pop();
			const auto d_q = e.d;

			bool is_good = true;
			for(const auto &r : R)
			{
				const auto d_r = U::distance(e.u->data, r->data, dim);
				//if(d_r*(level+1)>d_q*alpha*(entrance->level+1))
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
		uint32_t level, bool extendCandidate=false, bool keepPrunedConnections=true)
	{
		/*
		(void)level, (void)extendCandidate, (void)keepPrunedConnections;
		return select_neighbors_simple(u,C,M);
		*/
		return select_neighbors_heuristic(u, C, M, level, extendCandidate, keepPrunedConnections);
	}

	uint32_t get_level_random()
	{
		// static thread_local int32_t anchor;
		// uint32_t esp;
		// asm volatile("movl %0, %%esp":"=a"(esp));
		static thread_local std::hash<std::thread::id> h;
		static thread_local std::mt19937 gen{h(std::this_thread::get_id())};
		static thread_local std::uniform_real_distribution<> dis(std::numeric_limits<double>::min(), 1.0);
		const uint32_t res = uint32_t(-log(dis(gen))*m_l);
		return res;
	}

	auto search_layer(const node &u, const std::vector<node*> &eps, uint32_t ef, uint32_t l_c) const; // To static
	auto search_layer_ex(const node &u, const std::vector<node*> &eps, uint32_t ef, uint32_t l_c) const; // To static
	auto get_threshold_m(uint32_t level){
		return level==0? m*2: m;
	}

	void fix_edge()
	{
		fprintf(stderr, "Start fixing edges...\n");

		for(int32_t l_c=entrance[0]->level; l_c>=0; --l_c)
		{
			parlay::sequence<parlay::sequence<std::pair<node*,node*>>> edge_add(n);

			parlay::parallel_for(0, n, [&](uint32_t i){
				auto &u = *node_pool[i];
				if(l_c>u.level) return;

				auto &edge_v = edge_add[i];
				edge_v.clear();
				for(auto *pv : neighbourhood(u,l_c))
				{
					const auto &nbh_v = neighbourhood(*pv,l_c);
					if(std::find_if(nbh_v.begin(),nbh_v.end(),[&](const node *pu_extant){
						return pu_extant==&u;
					})==nbh_v.end())
						edge_v.emplace_back(pv, &u);
				}
			});

			auto edge_add_flatten = parlay::flatten(edge_add);
			auto edge_add_grouped = parlay::group_by_key(edge_add_flatten);

			parlay::parallel_for(0, edge_add_grouped.size(), [&](size_t j){
				node *pv = edge_add_grouped[j].first;
				auto &nbh_v = neighbourhood(*pv,l_c);
				auto &nbh_v_add = edge_add_grouped[j].second;

				nbh_v.insert(nbh_v.end(),nbh_v_add.begin(), nbh_v_add.end());
			});
		}
	}
};


template<typename U, template<typename> class Allocator>
template<typename G>
HNSW<U,Allocator>::HNSW(const std::string &filename_model, G getter)
{
	std::ifstream model(filename_model, std::ios::binary);
	if(!model.is_open())
		throw std::runtime_error("Failed to open the model");

	auto read = [&](auto &data, auto ...args){
		auto read_impl = [&](auto &f, auto &data, auto ...args){
			using T = std::remove_reference_t<decltype(data)>;
			if constexpr(std::is_pointer_v<std::decay_t<T>>)
			{
				auto read_array = [&](auto &data, size_t size, auto ...args){
					for(size_t i=0; i<size; ++i)
						f(f, data[i], args...);
				};
				// use the array extent as the size
				if constexpr(sizeof...(args)==0 && std::is_array_v<T>)
				{
					read_array(data, std::extent_v<T>);
				}
				else
				{
					static_assert(sizeof...(args), "size was not provided");
					read_array(data, args...);
				}
			}
			else
			{
				static_assert(std::is_standard_layout_v<T>);
				model.read((char*)&data, sizeof(data));
			}
		};
		read_impl(read_impl, data, args...);
	};

	char model_type[5] = {'\000'};
	read(model_type, 4);
	if(strcmp(model_type,"HNSW"))
		throw std::runtime_error("Wrong type of model");
	uint32_t version;
	read(version);
	if(version>1)
		throw std::runtime_error("Unsupported version");

	size_t code_U, size_node;
	read(code_U);
	read(size_node);
	if((typeid(U).hash_code()^sizeof(U))!=code_U)
		throw std::runtime_error("Inconsistent type `U`");
	if(sizeof(node)!=size_node)
		throw std::runtime_error("Inconsistent type `node`");

	// read parameter configuration
	read(dim);
	read(m_l);
	read(m);
	read(ef_construction);
	read(alpha);
	read(n);
	puts("Configuration loaded");
	printf("dim = %u\n", dim);
	printf("m_l = %f\n", m_l);
	printf("m = %u\n", m);
	printf("efc = %u\n", ef_construction);
	printf("alpha = %f\n", alpha);
	printf("n = %u\n", n);
	// read indices
	std::unordered_map<uint32_t,node*> addr;
	node_pool.reserve(n);
	for(uint32_t i=0; i<n; ++i)
	{
		auto *u = new node;
		read(u->level);
		uint32_t id_u;
		read(id_u);
		u->data = getter(id_u);
		addr[id_u] = u;
		node_pool.push_back(u);
	}
	for(node *u : node_pool)
	{
		u->neighbors = new std::vector<node*>[u->level+1];
		for(uint32_t l=0; l<=u->level; ++l)
		{
			size_t size;
			read(size);
			auto &nbh_u = u->neighbors[l];
			nbh_u.reserve(size);
			for(size_t i=0; i<size; ++i)
			{
				uint32_t id_v;
				read(id_v);
				nbh_u.push_back(addr.at(id_v));
			}
		}
	}
	// read entrances
	size_t size;
	read(size);
	entrance.reserve(size);
	for(size_t i=0; i<size; ++i)
	{
		uint32_t id_u;
		read(id_u);
		entrance.push_back(addr.at(id_u));
	}
}

template<typename T, template<typename> class Allocator>
template<typename Iter>
HNSW<T,Allocator>::HNSW(Iter begin, Iter end, uint32_t dim_, float m_l_, uint32_t m_, uint32_t ef_construction_, float alpha_, float batch_base, bool do_fixing)
	: dim(dim_), m_l(m_l_), m(m_), ef_construction(ef_construction_), alpha(alpha_), n(std::distance(begin,end))
{
	static_assert(std::is_same_v<typename std::iterator_traits<Iter>::value_type, T>);
	static_assert(std::is_base_of_v<
		std::random_access_iterator_tag, typename std::iterator_traits<Iter>::iterator_category>);

	if(n==0) return;

	auto perm = parlay::random_permutation<uint32_t>(n, 1206);
	auto rand_seq = parlay::delayed_seq<T>(n, [&](uint32_t i){
		return *(begin+perm[i]);
	});

	const auto level_ep = get_level_random();
	node *entrance_init = allocator.allocate(1);
	new(entrance_init) node{level_ep, *rand_seq.begin(), new std::vector<node*>[level_ep+1]/*anything else*/};
	#if DEBUG_OUTPUT
		fprintf(stderr, "[%u] at lv.%u (%.2f,%.2f)**\n", 0, level_ep, rand_seq.begin()->x, rand_seq.begin()->y);
	#endif
	node_pool.push_back(entrance_init);
	entrance.push_back(entrance_init);

	uint32_t batch_begin=0, batch_end=1;
	float progress = 0.0;
	while(batch_end<n)
	{
		batch_begin = batch_end;
		batch_end = std::min(n, (uint32_t)std::ceil(batch_begin*batch_base));
		/*
		if(batch_end>batch_begin+100)
			batch_end = batch_begin+100;
		*/
		// batch_end = batch_begin+1;

		insert(rand_seq.begin()+batch_begin, rand_seq.begin()+batch_end, true);
		insert(rand_seq.begin()+batch_begin, rand_seq.begin()+batch_end, false);

		if(batch_end>n*(progress+0.1))
		{
			progress = float(batch_end)/n;
			fprintf(stderr, "Done: %.2f\n", progress);
		}
	}

	if(do_fixing) fix_edge();

	#if DEBUG_OUTPUT
		for(const auto *pu : node_pool)
		{
			fprintf(stderr, "[%u] (%.2f,%.2f)\n", U::get_id(pu->data), pu->data[0], pu->data[1]);
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
template<typename Iter>
void HNSW<U,Allocator>::insert(Iter begin, Iter end, bool from_blank)
{
	const auto level_ep = entrance[0]->level;
	const auto size_batch = std::distance(begin,end);
	auto node_new = std::make_unique<node*[]>(size_batch);
	auto nbh_new = std::make_unique<std::vector<node*>[]>(size_batch);
	auto eps = std::make_unique<std::vector<node*>[]>(size_batch);
	//const float factor_m = from_blank? 0.5: 1;
	const auto factor_m = 1;

	// first, query the nearest point as the starting point for each node to insert
	if(from_blank)
	{
	parlay::parallel_for(0, size_batch, [&](uint32_t i){
		const T &q = *(begin+i);
		const auto level_u = get_level_random();
		auto *const pu = allocator.allocate(1);		// TODO: add pointer manager

		new(pu) node{level_u, q, new std::vector<node*>[level_u+1]};
		#if DEBUG_OUTPUT
			fprintf(stderr, "[%u] at lv.%u (%.2f,%.2f)\n", U::get_id(q), level_u, q.x, q.y);
		#endif
		node_new[i] = pu;
	});
	}
	else
	{
	parlay::parallel_for(0, size_batch, [&](uint32_t i){
		node_new[i] = *(node_pool.end()-size_batch+i);
	});
	}

	// TODO: merge ops
	parlay::parallel_for(0, size_batch, [&](uint32_t i){
		auto &u = *node_new[i];
		const auto level_u = u.level;
		auto &eps_u = eps[i]; 
		// eps_u.push_back(entrance);
		eps_u = entrance;
		for(uint32_t l=level_ep; l>level_u; --l)
		{
			const auto res = search_layer(u, eps_u, 1, l); // TODO: optimize
			eps_u.clear();
			eps_u.push_back(res.top().u);
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
			auto neighbors_queue = select_neighbors(u.data, res, get_threshold_m(l_c)*factor_m, l_c);
			// move the content from `neighbors_queue` to `u.neighbors[l_c]`
			auto &nbh_u = nbh_new[i];
			auto &edge_u = edge_add[i];
			nbh_u.clear();
			edge_u.clear();
			nbh_u.reserve(neighbors_queue.size());
			edge_u.reserve(neighbors_queue.size());
			for(uint32_t j=0; neighbors_queue.size()>0; ++j)
			{
				auto *pv = neighbors_queue.top().u;
				neighbors_queue.pop();
				// nbh_u[j] = pv;
				// edge_u[j] = std::make_pair(pv, &u);
				nbh_u.push_back(pv);
				edge_u.emplace_back(pv, &u);
			}

			eps_u.clear();
			while(res.size()>0)
			{
				eps_u.push_back(res.top().u); // TODO: optimize
				res.pop();
			}
		});

		parlay::parallel_for(0, size_batch, [&](uint32_t i){
			auto &u = *node_new[i];
			if((uint32_t)l_c<=u.level)
				neighbourhood(u,l_c) = std::move(nbh_new[i]);
		});

		// now we add edges in the other direction
		auto edge_add_flatten = parlay::flatten(edge_add);
		auto edge_add_grouped = parlay::group_by_key(edge_add_flatten);

		parlay::parallel_for(0, edge_add_grouped.size(), [&](size_t j){
			node *pv = edge_add_grouped[j].first;
			auto &nbh_v = neighbourhood(*pv,l_c);
			auto &nbh_v_add = edge_add_grouped[j].second;

			for(auto it=nbh_v_add.begin(); it!=nbh_v_add.end();)
			{
				bool is_extant = *it==pv||std::find_if(nbh_v.begin(), nbh_v.end(), [&](const node *pu_extant){
					return *it==pu_extant;
				})!=nbh_v.end();
				it = is_extant? nbh_v_add.erase(it): std::next(it);
			}

			const uint32_t size_nbh_total = nbh_v.size()+nbh_v_add.size();

			const auto m_s = get_threshold_m(l_c)*factor_m;
			if(size_nbh_total>m_s)
			{
				auto dist_nbh = std::make_unique<dist[]>(size_nbh_total);
				for(size_t k=0; k<nbh_v.size(); ++k)
					dist_nbh[k] = dist{U::distance(nbh_v[k]->data,pv->data,dim), nbh_v[k]};
				for(size_t k=0; k<nbh_v_add.size(); ++k)
					dist_nbh[k+nbh_v.size()] = dist{U::distance(nbh_v_add[k]->data,pv->data,dim), nbh_v_add[k]};

				std::sort(dist_nbh.get(), dist_nbh.get()+size_nbh_total, farthest());

				nbh_v.resize(m_s);
				for(size_t k=0; k<m_s; ++k)
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
		entrance.clear();
		entrance.push_back(node_highest);
	#if DEBUG_OUTPUT
		fprintf(stderr, "[%u]** at lev %u\n", U::get_id(entrance->data), entrance->level);
	#endif
	}
	else if(node_highest->level==level_ep)
		entrance.push_back(node_highest);

	// and add new nodes to the pool
	if(from_blank)
	node_pool.insert(node_pool.end(), node_new.get(), node_new.get()+size_batch);
}

template<typename U, template<typename> class Allocator>
auto HNSW<U,Allocator>::search_layer(const node &u, const std::vector<node*> &eps, uint32_t ef, uint32_t l_c) const
{
	auto W_ex = search_layer_ex(u, eps, ef, l_c);
	std::priority_queue<dist,std::vector<dist>,farthest> W;
	while(!W_ex.empty())
	{
		W.push(W_ex.top());
		W_ex.pop();
	}
	return W;
}

template<typename U, template<typename> class Allocator>
auto HNSW<U,Allocator>::search_layer_ex(const node &u, const std::vector<node*> &eps, uint32_t ef, uint32_t l_c) const
{
	std::vector<bool> visited(n);
	std::priority_queue<dist_ex,std::vector<dist_ex>,nearest> C;
	std::priority_queue<dist_ex,std::vector<dist_ex>,farthest> W;

	for(auto *ep : eps)
	{
		visited[U::get_id(ep->data)] = true;
		const auto d = U::distance(u.data,ep->data,dim);
		C.push({d,ep,1});
		W.push({d,ep,1});
	}

	while(C.size()>0)
	{
		const auto &dc = C.top().depth;
		const auto &c = *C.top().u;
		const auto &f = *W.top().u;
		C.pop();
		if(U::distance(c.data,u.data,dim)>U::distance(f.data,u.data,dim))
			break;
		for(auto *pv: neighbourhood(c, l_c))
		{
			if(visited[U::get_id(pv->data)]) continue;
			visited[U::get_id(pv->data)] = true;
			const auto &f = *W.top().u;
			if(W.size()<ef||U::distance(pv->data,u.data,dim)<U::distance(f.data,u.data,dim))
			{
				const auto d = U::distance(u.data,pv->data,dim);
				C.push({d,pv,dc+1});
				W.push({d,pv,dc+1});
				if(W.size()>ef) W.pop();
			}
		}
	}
	return W;
}


template<typename U, template<typename> class Allocator>
std::vector<std::pair<uint32_t,double>> HNSW<U,Allocator>::search(const T &q, uint32_t k, uint32_t ef)
{
	auto res_ex = search_ex(q,k,ef);
	std::vector<std::pair<uint32_t,double>> res;
	res.reserve(res_ex.size());
	for(const auto &e : res_ex)
		res.emplace_back(std::get<0>(e), std::get<2>(e));

	return res;
}

template<typename U, template<typename> class Allocator>
std::vector<std::tuple<uint32_t,uint32_t,double>> HNSW<U,Allocator>::search_ex(const T &q, uint32_t k, uint32_t ef)
{
	node u{0, q, nullptr}; // To optimize
	std::priority_queue<dist,std::vector<dist>,farthest> W;
	auto eps = entrance;
	for(int l_c=entrance[0]->level; l_c>0; --l_c) // TODO: fix the type
	{
		W = search_layer(u, eps, 1, l_c);
		eps.clear();
		eps.push_back(W.top().u);
		/*
		while(!W.empty())
		{
			eps.push_back(W.top().u);
			W.pop();
		}
		*/
	}
	auto W_ex = search_layer_ex(u, eps, ef, 0);
	W_ex = select_neighbors_simple(q, W_ex, k);
	std::vector<std::tuple<uint32_t,uint32_t,double>> res;
	while(W_ex.size()>0)
	{
		res.push_back({U::get_id(W_ex.top().u->data), W_ex.top().depth, W_ex.top().d});
		W_ex.pop();
	}
	return res;
}

template<typename U, template<typename> class Allocator>
void HNSW<U,Allocator>::save(const std::string &filename_model) const
{
	std::ofstream model(filename_model, std::ios::binary);
	if(!model.is_open())
		throw std::runtime_error("Failed to create the model");

	const auto write = [&](const auto &data, auto ...args){
		auto write_impl = [&](auto &f, const auto &data, auto ...args){
			using T = std::remove_reference_t<decltype(data)>;
			if constexpr(std::is_pointer_v<std::decay_t<T>>)
			{
				auto write_array = [&](const auto &data, size_t size, auto ...args){
					for(size_t i=0; i<size; ++i)
						f(f, data[i], args...);
				};
				// use the array extent as the size
				if constexpr(sizeof...(args)==0 && std::is_array_v<T>)
				{
					write_array(data, std::extent_v<T>);
				}
				else
				{
					static_assert(sizeof...(args), "size was not provided");
					write_array(data, args...);
				}
			}
			else
			{
				static_assert(std::is_standard_layout_v<T>);
				model.write((const char*)&data, sizeof(data));
			}
		};
		write_impl(write_impl, data, args...);
	};
	// write header (version number, type info, etc)
	write("HNSW", 4);
	write(uint32_t(1));
	write(typeid(U).hash_code()^sizeof(U));
	write(sizeof(node));
	// write parameter configuration
	write(dim);
	write(m_l);
	write(m);
	write(ef_construction);
	write(alpha);
	write(n);
	// write indices
	for(const auto *u : node_pool)
	{
		write(u->level);
		write(U::get_id(u->data));
	}
	for(const auto *u : node_pool)
	{
		for(uint32_t l=0; l<=u->level; ++l)
		{
			write(u->neighbors[l].size());
			for(const auto *v : u->neighbors[l])
				write(U::get_id(v->data));
		}
	}
	// write entrances
	write(entrance.size());
	for(const auto *u : entrance)
		write(U::get_id(u->data));
}

} // namespace HNSW

#endif // _HNSW_HPP

