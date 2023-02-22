#pragma once
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/format.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
template<class mymesh>
class loop {


	typedef typename mymesh::LCC_3					LCC_3;
	typedef typename LCC_3::Dart_handle                                 Dart_handle;
	typedef typename LCC_3::Point                                       Point;
	using C0 = typename mymesh::C0;
	using C1 = typename mymesh::C1;
	using C2 = typename mymesh::C2;
	using C3 = typename mymesh::C3;
public:
	/*readmesh() {};*/
	loop(mymesh& out) :mesh(out.mesh()),
		VC0(out.VC0), VC1(out.VC1), VC2(out.VC2), VC3(out.VC3)
	{


	};
	loop(mymesh& out, std::string outfile) :mesh(out.mesh()),
		VC0(out.VC0), VC1(out.VC1), VC2(out.VC2), VC3(out.VC3), output_file(outfile)
	{


	};
	loop(mymesh& out, int rn, std::string outfile = "") :mesh(out.mesh()), reduce_number(rn),
		VC0(out.VC0), VC1(out.VC1), VC2(out.VC2), VC3(out.VC3), output_file(outfile)
	{


	};
	bool  get_weighted_loop();
	void hex_unweighted_loop();

protected:
	LCC_3& mesh;

	std::vector<C0>& VC0;
	std::vector<C1>& VC1;
	std::vector<C2>& VC2;
	std::vector<C3>& VC3;
	std::string output_file;
	int reduce_number = 0;
};
template <class T, class C>
class Cycle
{
public:
	T* head()
	{
		return m_cycle.empty() ? NULL : *(m_cycle.rbegin());
	}

	void add(T* pw)
	{
		if (!(pw->generator))
			return;
		auto [it, flag] = m_cycle.insert(pw);
		if (!flag) {
			m_cycle.erase(it);
		}


	};

	bool empty() { return m_cycle.empty(); };

	void print()
	{
		for (auto viter = m_cycle.begin(); viter != m_cycle.end(); viter++)
		{
			T* pv = *viter;
			std::cout << pv->id << " ";
		}
		std::cout << std::endl;
	}

protected:
	std::set<T*, C> m_cycle;
};

template<class C>
struct classcomp {
	bool operator() (const C* lhs, const C* rhs) const
	{
		return (lhs->id) < (rhs->id);
	}
};
template<typename E>
struct finish_city :
	public boost::base_visitor<finish_city<typename E>>
{

	finish_city(std::vector<int>& temp, std::vector<E>& temp2) :
		pre(temp),
		preedge(temp2)
	{
		//std::cerr << "pre size is " << pre.size() << '\n';
		//std::cerr << "temp size is "<<temp.size() << '\n';
	}
	typedef boost::on_tree_edge event_filter;
	template <class Edge, class Graph>
	inline void operator()(Edge e, Graph& g) {
		pre[e.m_target] = e.m_source;
		preedge[e.m_target] = e;
		/*std::cerr << boost::format("%d >= %d\n") % e.m_target % pre.size();*/
	   /* if (e.m_target >= pre.size()) {
			std::cerr << boost::format("%d >= %d\n")% e.m_target% pre.size();
		}*/
		/*std::cerr << boost::format("(%d -> %d) ") % e.m_source % e.m_target;*/
	}
	std::vector<E>& preedge;
	std::vector<int>& pre;
};


struct distance_city :
	public boost::base_visitor<distance_city>
{

	distance_city(std::vector<int>& temp) :
		pre(temp)

	{

	}
	typedef boost::on_tree_edge event_filter;
	template <class Edge, class Graph>
	inline void operator()(Edge e, Graph& g) {
		pre[e.m_target] = pre[e.m_source] + 1;

		/*std::cerr << boost::format("%d >= %d\n") % e.m_target % pre.size();*/
	   /* if (e.m_target >= pre.size()) {
			std::cerr << boost::format("%d >= %d\n")% e.m_target% pre.size();
		}*/
		/*std::cerr << boost::format("(%d -> %d) ") % e.m_source % e.m_target;*/
	}

	std::vector<int>& pre;
};
template<class mymesh>
inline bool loop<mymesh>::get_weighted_loop()
{
	using namespace boost;
	typedef property< edge_weight_t, double, property< edge_name_t, int> > Flow;
	typedef adjacency_list< vecS, vecS, undirectedS, no_property,
		Flow >
		Graph;
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	typedef graph_traits< Graph >::edge_descriptor Edge;
	typedef std::pair< int, int > E;
	double mmax = -1;
	int start = 0;
	for (int i = 0; i < VC0.size(); i++) {
		Point p = mesh.point(VC0[i].da);
		if (p.y() > mmax) {
			mmax = p.y();
			start = i;
		}
	}

	VC0[start].sharp = true;
	std::vector< Edge > spanning_tree;
	std::vector<int> depth(VC0.size(), 0);
	std::vector<E> edge_array;
	std::vector<Flow> weights;
	Graph subg;
	for (int i = 0; i < VC1.size(); i++) {
		auto pE = VC1[i].da;
		auto pE1 = mesh.beta<1>(pE);
		int id1 = mesh.info<0>(pE);
		int id2 = mesh.info<0>(pE1);

		edge_array.push_back(E(id1, id2));
		double aa = std::sqrt(CGAL::squared_distance(mesh.point(pE), mesh.point(pE1)));
		weights.push_back(Flow(aa, i));
	}


	Graph g(edge_array.begin(), edge_array.end(), weights.begin(), VC0.size());
	if (num_edges(g) != VC1.size()) {
		std::cerr << "num_edges(g) != VC1.size()\n";
	}

	{
		/*std::vector< int > p(num_vertices(g));*/
		std::vector<int> p(num_vertices(g));
		std::vector<double> dis(num_vertices(g));
		/*boost::associative_property_map< std::unordered_map<int, int> >
			address_map(p);*/
		dijkstra_shortest_paths(g, start, predecessor_map(p.data()).distance_map(dis.data()));
		if (abs(dis[start]) > 0.0000001) {
			std::cerr << "dis[start]!=0\n";
		}
		int count = 0;
		for (int i = 0; i < p.size(); i++) {
			count += ((i == p[i]) ? 1 : 0);
		}
		if (count != 1)
			std::cerr << "count!=1\n";

		std::unordered_set<E, boost::hash<E>> eset;
		for (auto& it : edge_array) {
			eset.insert({ it.first,it.second });
		}
		Edge e;
		bool flag;
		count = 0;
		for (int i = 0; i < num_vertices(g); i++) {
			if (i == p[i])	count++;
			if (eset.find({ i,p[i] }) != eset.end()) {
				//flag:prevent  the case i==p[i]
				tie(e, flag) = edge(i, p[i], g);
				if (flag)
					spanning_tree.push_back(e);
			}
			else {
				tie(e, flag) = edge(p[i], i, g);
				if (flag)
					spanning_tree.push_back(e);
			}
		}
		if (count != 1) {
			std::cerr << "i==p[i] !=1\n";
		}
		if (spanning_tree.size() + 1 != VC0.size()) {
			std::cerr << "spanning_tree.size() + 1 != m_pMesh->vertices().size()";
		}



		if (spanning_tree.size() + 1 != g.m_vertices.size()) {
			std::cerr << "spanning_tree.size() + 1 != g.m_vertices.size()\n";
		}
		if (g.m_vertices.size() != VC0.size()) {
			std::cerr << "g.m_vertices.size() != m_pMesh->vertices().size()\n";
		}


		{
			edge_array.clear();
			weights.clear();
			auto ww = get(edge_weight, g);
			auto name = get(edge_name, g);
			for (auto& it : spanning_tree) {
				edge_array.push_back(E(it.m_source, it.m_target));

				weights.push_back(Flow(ww[it], name[it]));
			}
			int N = VC0.size();

			subg = Graph(edge_array.begin(), edge_array.end(), weights.begin(), N);
		}


		property_map< Graph, edge_name_t >::type name = get(edge_name, g);

		std::unordered_set<int> gset;
		std::unordered_map<int, double> e_dis_map;

		/*std::set<std::pair<int, int>> gset;*/
		for (auto& it : spanning_tree)
			gset.insert(name[it]);
		if (spanning_tree.size() != gset.size()) {
			std::cout << "spanning_tree.size() != gset.size()\n";
		}
		count = 0;
		for (auto& ei : spanning_tree) {
			bool got = (gset.find(name[ei]) == gset.end());
			if (got)	count++;
		}
		if (count > 0)
			std::cerr << "spanning_tree check error\n";
		if (num_edges(g) != VC1.size()) {
			std::cout << "num_edges(g) != VC1.size()\n";
		}
		count = 0;
		int count1 = 0;
		{
			/*{
			std::ofstream vc1file("E:\\hjs\\matlabcode\\hjs1\\file1\\vcn.data");
			vc1file << format("%d %d %d %d\n") % VC0.size() % VC1.size() % VC2.size() % VC3.size();

		}
		{
			std::ofstream vc1file("E:\\hjs\\matlabcode\\hjs1\\file1\\vc1.data");

			for (int i = 0; i < VC1.size(); i++) {
				auto pE = VC1[i].da;
				auto pE1 = mesh.beta<1>(pE);
				int id1 = mesh.info<0>(pE);
				int id2 = mesh.info<0>(pE1);
				vc1file << format("%d %d\n") % id1% id2;
			}
		}
		{
			std::ofstream vc1file("E:\\hjs\\matlabcode\\hjs1\\file1\\edge.data");
			for (auto[ei, ei_end] = edges(g); ei != ei_end; ++ei) {

				vc1file << format("%d %d\n") % ei->m_source% ei->m_target;

			}
		}

		{
			std::ofstream vc1file("E:\\hjs\\matlabcode\\hjs1\\file1\\span.data");
			for (auto& ei : spanning_tree) {
				vc1file << format("%d %d\n") % ei.m_source%ei.m_target;

			}
		}*/
		/*{
			std::ofstream vc1file("E:\\hjs\\matlabcode\\hjs1\\file1\\gset.data");
			for (auto& ei : gset) {
				auto[a, b] = ei;
				vc1file << format("%d %d\n") % a%b;

			}
		}*/
		}
		Graph generator_tree;
		{
			edge_array.clear();
			weights.clear();
			auto ww = get(edge_all, g);
			for (auto [ei, ei_end] = edges(g); ei != ei_end; ++ei) {
				bool got = (gset.find(name[*ei]) == gset.end());
				if (got) {
					count++;
					edge_array.push_back(E(ei->m_source, ei->m_target));

					weights.push_back(ww[*ei]);
				}
				else count1++;
			}
			int N = VC0.size();

			generator_tree = Graph(edge_array.begin(), edge_array.end(), weights.begin(), N);
		}

		if (count + count1 != VC1.size()) {
			std::cerr << "count + count1 != VC1.size()\n";
		}
		if (count + spanning_tree.size() != VC1.size()) {
			std::cerr << "edges(g),gset,VC1 check error\n";

		}

		{
			std::vector<int> parent(num_vertices(g));
			std::vector<bool> visited(num_vertices(g), 0);
			name = get(edge_name, generator_tree);

			std::function<int(int)> find = [&](int x) {
				if (x != parent[x])
					parent[x] = find(parent[x]);
				return parent[x];
			};
			std::function<void(int)>   tarjan = [&](int u) {
				parent[u] = u;
				visited[u] = 1;
				for (auto [ei, en] = adjacent_vertices(u, subg); ei != en; ei++) {
					int t1 = *ei;
					if (!visited[t1]) {
						tarjan(t1);
						parent[t1] = u;
					}
				}
				for (auto [ei, en] = out_edges(u, generator_tree); ei != en; ei++) {
					if (visited[ei->m_target]) {
						int a = name[*ei];
						name[*ei] = find(ei->m_target);
						//std::cerr << format("%d %d\n")%a%name[*ei];

					}
				}



			};
			for (int i = 0; i < num_vertices(g); i++) {
				parent[i] = i;
			}
			tarjan(start);
		}


		int vid = 1;
		count = 0;
		graph_traits<Graph>::edge_iterator ei, ei_end;
		for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
		{
			C1& pE = VC1[name[*ei]];
			bool got = (gset.find(name[*ei]) == gset.end());
			if (got) {
				//means it does not belong to tree.
				pE.generator = true;

			}
			else {
				pE.id = vid++;
				/*			pE->generator() = false;*/
				count++;
			}
		}
		std::multimap<double, C1*> mmap;
		auto gname = get(edge_name, generator_tree);
		for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
		{
			C1& pE = VC1[name[*ei]];
			int s = ei->m_source, t = ei->m_target;
			auto [e, f] = edge(s, t, generator_tree);
			bool got = (gset.find(name[*ei]) == gset.end());
			if (got) {
				mmap.insert({ dis[s] + dis[t] - 2 * dis[gname[e]],&pE });


			}

		}
		for (auto& [it, CC] : mmap) {
			CC->id = vid++;
		}
		if (count != spanning_tree.size()) {
			std::cerr << "count != spanning_tree.size()\n";
		}

		count = 0;
		for (auto& it : VC1) {
			count += (it.generator == false ? 1 : 0);
		}
		if (count + 1 != VC0.size()) {
			std::cerr << "count + 1 != m_pMesh->vertices().size()\n";
		}
	}






	{


		int count = 0;
		struct C2s {
			int id;
			C2* cc;
		};

		std::vector<C2*> fvec;
		fvec.reserve(VC2.size());
		{
			std::vector<C2s> tvec;
			tvec.reserve(VC2.size());
			int mmax = -1;
			for (auto& it : VC2) {
				mmax = -1;
				Dart_handle dc1 = it.da;
				mmax = std::max((VC1.data() + mesh.info<1>(dc1))->id, mmax);
				dc1 = mesh.beta<1>(dc1);
				mmax = std::max((VC1.data() + mesh.info<1>(dc1))->id, mmax);
				dc1 = mesh.beta<1>(dc1);
				mmax = std::max((VC1.data() + mesh.info<1>(dc1))->id, mmax);
				tvec.push_back({ mmax,&it });
			}
			std::sort(tvec.begin(), tvec.end(), [](C2s& c1, C2s& c2)
				{return c1.id < c2.id; });
			for (auto& it : tvec)
				fvec.push_back(it.cc);
		}
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

		//shuffle(fvec.begin(), fvec.end(), std::default_random_engine(seed));


		int tcount = 0;
		for (auto& fiter : fvec)
		{
			/*std::cout << format("%d/%d\n")% facecount++%m_pMesh->faces().size();*/
			C2* pF = fiter;

			//std::cout << format("%d %d\n") % tcount++%fvec.size();
			//insert your code here



			Cycle<C1, classcomp<C1>> ecycle;
			C2* tpF = pF;

			while (1) {
				Dart_handle dc1 = tpF->da;
				ecycle.add(VC1.data() + mesh.info<1>(dc1));
				dc1 = mesh.beta<1>(dc1);
				ecycle.add(VC1.data() + mesh.info<1>(dc1));
				dc1 = mesh.beta<1>(dc1);
				ecycle.add(VC1.data() + mesh.info<1>(dc1));



				C1* pE = ecycle.head();
				if (pE == NULL) {
					count++;
					pF->generator = true;
					break;
				}
				if (!(pE->generator)) {
					if (pE->id != 0)
						std::cerr << "!(pE->generator()) and PE is not loop\n";
					else
						std::cerr << "!(pE->generator())\n";
				}
				if (pE->pair == -1) {
					pE->pair = pF - VC2.data();
					break;
				}
				tpF = VC2.data() + pE->pair;
			}

		}


		count = 0;

		for (auto& it : VC1) {
			count += (it.generator) && (it.pair == -1) ? 1 : 0;
		}
		int g3 = VC0.size() - VC1.size() + VC2.size() - VC3.size();
		if (1 - g3 != count) {
			std::cerr << format("1 - g3(%d) != count(%d)\n") % g3 % count;
		}
		/*assert(count == 2);*/
		std::cout << format("loop num:%d\n") % count;
		if (count != 2) {
			/*std::cerr << "count != 2\n";*/
		}

	}

	{

		class mycount {
		public:
			mycount(std::vector<C1*>& ff) :i1(0), i2(0) {
				for (auto& it : ff)
					eset.insert(it);

			};
			void count(C1* it) {
				if (eset.find(it) != eset.end()) {
					i1++;
				}
				else {
					i2++;
				}
			}
			void reset() {
				i1 = 0;
				i2 = 0;
			}
			void print() {
				std::cerr << format("total edges:%5d,inter edges:%5d,boundary edges:%5d\n") % (i1 + i2) % i2 % i1;
			}

		private:
			std::unordered_set<C1*> eset;
			int i1, i2;

		};
		/*std::vector<C1*> ff;
		{
			for (auto a1 = mesh.one_dart_per_cell<2>().begin(), a2 = mesh.one_dart_per_cell<2>().end(); a1 != a2; a1++) {

				if (!mesh.is_free<3>(a1))
					continue;


				for (auto dc1 = mesh.one_dart_per_incident_cell<1, 2>(a1).begin(),
					dc2 = mesh.one_dart_per_incident_cell<1, 2>(a1).end();
					dc1 != dc2; dc1++) {

					ff.push_back(VC1.data() + mesh.attribute<1>(dc1)->info());
				}


			}
		}
		mycount cc(ff);*/

		int N = VC0.size();
		std::vector<int> pre(N, 0);
		std::vector<Edge> preedge(N, Edge());
		if (subg.m_edges.size() + 1 != subg.m_vertices.size()) {
			std::cerr << "subg is wrong\n";
		}


		std::vector<C1*> evec;
		std::vector<std::vector<Edge>> loop;


		for (auto& it : VC1) {
			if ((it.generator) && (it.pair == -1)) {
				evec.push_back(&it);
			}
		}
		int loopn = -1;

		for (auto& it : evec) {
			/*if (m_boundary_edges.find(it) == m_boundary_edges.end()) {
				std::cerr << "m_boundary_edges.find(it) == m_boundary_edges.end()\n";
			}*/
			loopn++;


			int id1 = mesh.info<0>(it->da);
			int id2 = mesh.info<0>(mesh.beta<1>(it->da));


			auto v1 = vertex(id1, subg);
			auto v2 = vertex(id2, subg);

			pre[id1] = id1;
			preedge[id1] = Edge();



			/*, visitor(make_bfs_visitor(fc))*/
			auto vis = visitor(make_bfs_visitor(
				finish_city<Edge>(pre, preedge)));
			breadth_first_search(subg, id1, vis);
			//check
			{
				int count = 0;
				for (int i = 0; i < N; i++) {
					/*  std::cout << pre[i] << ' ';*/
					if (pre[i] == i)
						count++;
				}
				if (count != 1) {
					std::cerr << format("pre is wrong %d\n") % count;
				}
				count = 0;
				for (int i = 0; i < N; i++) {

					if (preedge[i] == Edge())
						count++;
				}
				if (count != 1) {
					std::cerr << format("preedge is wrong %d\n") % count;
				}

			}
			int i1, i2 = id2;

			std::list<Edge> tel;
			std::vector<Edge> el;
			std::vector<Edge> e2;



			Edge ee = preedge[i2];
			auto name = get(edge_name, g);


			std::vector<Edge> vec;
			//std::vector<int> veccount(g.m_vertices.size(), 0);
			/*cc.reset();*/
			while (ee != Edge()) {
				tel.push_back(edge(target(ee, subg), source(ee, subg), g).first);
				i2 = ee.m_source;
				ee = preedge[i2];
			}
			//veccount do not calculate endpoint of it;
			//simplify loop
			{

				std::vector<int> pre(num_vertices(g));
				auto reduce_weight_loop = [&pre](std::vector<Edge>& el, int n, Graph& g)->std::vector<Edge> {
					const double INF = std::numeric_limits<double>::infinity();
					Vertex s = source(el[0], g), t = source(el[n], g);
					auto weight_m = get(edge_weight, g);
					double sum = 0;
					for (int i = 0; i < n; i++) {
						sum += weight_m[el[i]];
					}


					typedef std::pair<double, int> iPair;

					std::priority_queue< iPair, std::vector <iPair>, std::greater<iPair> > pq;
					std::vector<double> dist(num_vertices(g), INF);
					pq.push(std::make_pair(0, s));
					dist[s] = 0;
					pre[s] = s;
					std::vector<bool> f(num_vertices(g), false);

					/* Looping till priority queue becomes empty (or all
					distances are not finalized) */
					while (!pq.empty())
					{
						// The first vertex in pair is the minimum distance
						// vertex, extract it from priority queue.
						// vertex label is stored in second of pair (it
						// has to be done this way to keep the vertices
						// sorted distance (distance must be first item
						// in pair)
						int u = pq.top().second;
						pq.pop();
						f[u] = true;
						if (u == t)
							break;
						// 'i' is used to get all adjacent vertices of a vertex
						for (auto [ei, en] = out_edges(u, g); ei != en; ei++) {
							int v = ei->m_target;
							double weight = weight_m[*ei];
							if (f[v] == false && dist[v] > dist[u] + weight)
							{
								// Updating distance of v
								dist[v] = dist[u] + weight;
								pre[v] = u;
								pq.push(std::make_pair(dist[v], v));
							}
						}

					}




					// invoke variant 2 of Dijkstra's algorithm
					//std::cout << format("<s:%d t:%d>\n") % s % t;


					if (abs(dist[t] - sum) < 0.00001) {
						return std::vector(el);
					}
					/*for (auto& it : pre)std::cout << it << ' ';*/
					//std::cout << '\n';
					std::list<Edge> el1;
					//std::cout << format("dis %.6f %.6f\n")%sum%dis[t];
					while (pre[t] != t) {
						el1.push_front(edge(pre[t], t, g).first);

						t = pre[t];
					}

					//ccc(el, g);
					std::vector<Edge> el2(el.size() - n + el1.size());
					std::copy(el.begin() + n, el.end(), el2.begin() + el1.size());
					//ccc(el2, g);
					std::copy(el1.begin(), el1.end(), el2.begin());

					//ccc(el2, g);

					return el2;
				};

				auto reduce_unweight_loop = [&pre](std::vector<Edge>& el, int n, Graph& g)->std::vector<Edge> {

					Vertex s = source(el[0], g), t = source(el[n], g);
					std::vector<int> dist(num_vertices(g), 0);
					int sum = n;


					std::vector<bool> visited(num_vertices(g), 0);


					// Create a queue for BFS
					std::list<int> queue;

					// Mark the current node as visited and enqueue it
					visited[s] = true;
					dist[s] = 0;
					pre[s] = s;
					queue.push_back(s);



					while (!queue.empty())
					{
						// Dequeue a vertex from queue and print it
						int u = queue.front();
						if (u == t)
							break;
						queue.pop_front();

						// Get all adjacent vertices of the dequeued
						// vertex s. If a adjacent has not been visited, 
						// then mark it visited and enqueue it
						for (auto [vi, vn] = adjacent_vertices(u, g); vi != vn; vi++)
						{

							if (!visited[*vi])
							{
								visited[*vi] = true;
								pre[*vi] = u;
								dist[*vi] = dist[u] + 1;
								queue.push_back(*vi);
							}
						}
					}


					// invoke variant 2 of Dijkstra's algorithm
					//std::cout << format("<s:%d t:%d>\n") % s % t;


					if (dist[t] == sum) {
						return std::vector(el);
					}
					/*for (auto& it : pre)std::cout << it << ' ';*/
					//std::cout << '\n';
					std::list<Edge> el1;
					//std::cout << format("dis %.6f %.6f\n")%sum%dis[t];
					while (pre[t] != t) {
						el1.push_front(edge(pre[t], t, g).first);

						t = pre[t];
					}

					//ccc(el, g);
					std::vector<Edge> el2(el.size() - n + el1.size());
					std::copy(el.begin() + n, el.end(), el2.begin() + el1.size());
					//ccc(el2, g);
					std::copy(el1.begin(), el1.end(), el2.begin());

					//ccc(el2, g);

					return el2;
				};
				tel.push_back(edge(id1, id2, g).first);
				el.resize(tel.size());
				std::copy(tel.begin(), tel.end(), el.begin());

				tel.clear();
				srand(time(NULL));

				int loop = 3000;


				auto ccc = [](auto& el, Graph& g) {
					typename property_map<Graph, vertex_index_t>::type
						index = get(vertex_index, g);
					std::cout << "out-edges: ";
					for (auto& it : el) {
						Vertex src = source(it, g), targ = target(it, g);
						std::cout << "(" << index[src] << ","
							<< index[targ] << ") ";
					};

					std::cout << std::endl;
				};
				while (loop--) {

					//std::cout << "-----start----\n";
					int ell = el.size();
					if (ell <= 3)
						break;
					int length = std::max((ell / 10), 2);

					int iSecret = rand() % ell;
					//std::cout << format("length:%d iSecret%d\n") % length%iSecret;
					e2.resize(el.size());


					std::copy(el.begin() + iSecret, el.end(), e2.begin());

					std::copy(el.begin(), el.begin() + iSecret, e2.begin() + ell - iSecret);

					//ccc(e2, g);
					el = reduce_weight_loop(e2, length, g);
					//std::cout << format("length:%d->%d\n") % e2.size()%el.size();
					//flag = (el.size() == e2.size()) ? 0 : 1;
					//ccc(el, g);
					//std::cout << "----end----\n\n";

				}


			}
			for (auto& it : el) {
				//VC1.data() + name[it]
				vec.push_back(it);
				//veccount[it.m_source]++;
				//veccount[it.m_target]++;
			}

			//std::cout << format("pass %d times\n") % veccount[start];
			std::cout << format("one loop finished\n");
			/*cc.count(it);

			cc.print();*/
			loop.push_back(vec);
		}
		std::vector<double>dis;
		auto name = get(edge_name, g);
		auto weight = get(edge_weight, g);
		std::vector<std::list<int>> loop_c1(VC1.size());
		dis.reserve(loop.size());
		int loopcount = 0;
		for (auto& vec : loop) {
			double sum = 0;
			for (Edge& it : vec) {
				loop_c1[name[it]].push_back(loopcount);
				sum += weight[it];
			}
			loopcount++;
			dis.push_back(sum);





		}
		std::vector<int>index_dis;
		index_dis.reserve(dis.size());
		std::vector<bool>index_ind(dis.size(), 1);


		for (auto& vec : loop_c1) {
			if (vec.size() < 2)
				continue;
			for (auto& it : vec)
				index_ind[it] = 0;
		}

		for (int i = 0; i < dis.size(); i++) {
			/*if(index_ind[i])*/
			index_dis.push_back(i);

		}
		std::cout << format("loop:%d real loop%d\n") % dis.size() % index_dis.size();
		std::sort(index_dis.begin(), index_dis.end(),
			[&dis](int a, int b) {return dis[a] > dis[b]; });
		for (auto& it : index_dis) {
			std::cout << format("%.3f\n") % dis[it];
		}
		/*{

	for (auto& ei : gset) {
		auto[a, b] = ei;
		vc1file << format("%d %d\n") % a%b;

	}
}*/
		if (output_file != "") {
			system(("mkdir " + output_file).c_str());
		}


		for (int i = 0; i < std::max(std::size_t(0), index_dis.size()); i++) {
			std::ofstream vc1file;
			if (output_file == "") {
				vc1file = std::ofstream((format("E:\\hjs\\matlabcode\\hjs1\\file1\\loop%d.data") % i).str());
			}
			else {
				vc1file = std::ofstream(output_file + (format("\\loop%d.data") % i).str());

			}
			vc1file << format("%1%\n") % dis[index_dis[i]];
			auto& vec = loop[index_dis[i]];
			/*cc.reset();*/
			for (auto& it : vec) {
				VC1[name[it]].sharp = true;
				vc1file << format("%1%\n") % mesh.point(VC0[it.m_source].da);

				/*	cc.count(it);*/
			}
			/*cc.print();*/
		}

	}
	return true;
}

template<class mymesh>
inline void loop<mymesh>::hex_unweighted_loop()
{
	//connected component must be 1
	using namespace boost;
	typedef property<edge_name_t, int> Flow;
	typedef adjacency_list< vecS, vecS, undirectedS, no_property,
		Flow >
		Graph;
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	typedef graph_traits< Graph >::edge_descriptor Edge;
	typedef std::pair< int, int > E;
	double mmax = -1;
	srand(time(NULL));
	int start = rand() % VC0.size();


	VC0[start].sharp = true;
	std::vector< Edge > spanning_tree;
	spanning_tree.reserve(VC1.size());
	std::vector<int> dis;
	std::vector<E> edge_array;
	edge_array.reserve(VC1.size());
	std::vector<Flow> weights;
	edge_array.reserve(VC1.size());
	Graph generator_tree;
	Graph subg;
	for (int i = 0; i < VC1.size(); i++) {
		auto pE = VC1[i].da;
		auto pE1 = mesh.beta<1>(pE);
		int id1 = mesh.info<0>(pE);
		int id2 = mesh.info<0>(pE1);

		edge_array.push_back(E(id1, id2));

		weights.push_back(Flow(i));
	}


	Graph g(edge_array.begin(), edge_array.end(), weights.begin(), VC0.size());
	auto name_g = get(edge_name, g);
	std::vector<int> dis_lca(num_edges(g), 0);
	std::vector<int> is_tree(num_edges(g), 0);
	//tree:1 generator:0 real generator:2
	if (num_edges(g) != VC1.size()) {
		std::cerr << "num_edges(g) != VC1.size()\n";
	}

	{
		/*std::vector< int > p(num_vertices(g));*/
		struct bfs_time_visitor : public default_bfs_visitor
		{



			bfs_time_visitor(std::vector<int>& opre, std::vector<int>& odis) :
				pre(opre), dis(odis) {}
			void tree_edge(Edge e, const Graph& g)const
			{
				pre[e.m_target] = e.m_source;
				dis[e.m_target] = dis[e.m_source] + 1;
			}
			std::vector<int>& dis;
			std::vector<int>& pre;

		};
		std::vector<int> p(num_vertices(g));
		dis = std::vector<int>(num_vertices(g));
		p[start] = start;
		dis[start] = 0;
		breadth_first_search(g, start, visitor(bfs_time_visitor(p, dis)));
		//prim_minimum_spanning_tree(g, p.data());		
		//std::ofstream hhh1("E:\\hjs\\c++\\hh1\\pro2\\build\\Release\\pre_g.data");

		int count = 0;
		for (int i = 0; i < p.size(); i++) {
			//hhh1 << format("%d %d\n") % i%p[i];
			count += ((i == p[i]) ? 1 : 0);
		}
		if (count != 1)
			std::cerr << "count!=1\n";
		for (int i = 0; i < num_vertices(g); i++) {
			auto [aa, bb] = edge(i, p[i], g);
			if (bb)
				is_tree[name_g[aa]] = 1;
		}
		count = 0;
		for (int i = 0; i < VC1.size(); i++) {
			count += is_tree[i];
		}
		if (count + 1 != VC0.size())
			std::cerr << "count += is_tree[i];\n";

		Edge e;
		bool flag;
		count = 0;
		for (int i = 0; i < num_vertices(g); i++) {

			if (i == p[i]) { count++; continue; }
			tie(e, flag) = edge(i, p[i], g);
			assert(flag);
			if (flag)
				spanning_tree.push_back(e);

		}
		if (count != 1) {
			std::cerr << "i==p[i] !=1\n";
		}
		if (spanning_tree.size() + 1 != VC0.size()) {
			std::cerr << "spanning_tree.size() + 1 != m_pMesh->vertices().size()";
		}



		if (spanning_tree.size() + 1 != g.m_vertices.size()) {
			std::cerr << "spanning_tree.size() + 1 != g.m_vertices.size()\n";
		}
		if (g.m_vertices.size() != VC0.size()) {
			std::cerr << "g.m_vertices.size() != m_pMesh->vertices().size()\n";
		}


		{
			edge_array.clear();
			weights.clear();

			auto name = get(edge_name, g);
			for (auto& it : spanning_tree) {
				edge_array.push_back(E(it.m_source, it.m_target));

				weights.push_back(Flow(name[it]));
			}
			int N = VC0.size();

			subg = Graph(edge_array.begin(), edge_array.end(), weights.begin(), N);
		}


		property_map< Graph, edge_name_t >::type name = get(edge_name, g);


		count = 0;
		int count1 = 0;


		{
			edge_array.clear();
			weights.clear();
			auto ww = get(edge_all, g);
			for (auto [ei, ei_end] = edges(g); ei != ei_end; ++ei) {

				if (is_tree[name_g[*ei]] == 0) {
					count++;
					edge_array.push_back(E(ei->m_source, ei->m_target));

					weights.push_back(ww[*ei]);
				}
				else count1++;
			}
			int N = VC0.size();

			generator_tree = Graph(edge_array.begin(), edge_array.end(), weights.begin(), N);
		}

		if (count + count1 != VC1.size()) {
			std::cerr << "count + count1 != VC1.size()\n";
		}
		if (count + spanning_tree.size() != VC1.size()) {
			std::cerr << "edges(g),gset,VC1 check error\n";

		}

		{
			std::vector<int> parent(num_vertices(g));
			std::vector<bool> visited(num_vertices(g), 0);
			name = get(edge_name, generator_tree);

			std::function<int(int)> find = [&](int x) {
				if (x != parent[x])
					parent[x] = find(parent[x]);
				return parent[x];
			};
			std::function<void(int)>   tarjan = [&](int u) {
				parent[u] = u;
				visited[u] = 1;
				for (auto [ei, en] = adjacent_vertices(u, subg); ei != en; ei++) {
					int t1 = *ei;
					if (!visited[t1]) {
						tarjan(t1);
						parent[t1] = u;
					}
				}
				for (auto [ei, en] = out_edges(u, generator_tree); ei != en; ei++) {
					if (visited[ei->m_target]) {
						int a = name[*ei];
						name[*ei] = find(ei->m_target);
						//std::cerr << format("%d %d\n")%a%name[*ei];

					}
				}



			};
			for (int i = 0; i < num_vertices(g); i++) {
				parent[i] = i;
			}
			tarjan(start);
		}


		int vid = 1;
		count = 0;
		std::vector<C1*> idC1(VC1.size() + 5, nullptr);
		graph_traits<Graph>::edge_iterator ei, ei_end;
		for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
		{
			C1& pE = VC1[name[*ei]];

			if (is_tree[name_g[*ei]] == 0) {
				pE.generator = true;

			}
			else {
				pE.id = vid;
				idC1[vid++] = &pE;
				//std::cout << format("tree vid:%d\n")%(vid-1);
				/*			pE->generator() = false;*/
				count++;
			}
		}
		std::multimap<int, C1*> mmap;
		auto gname = get(edge_name, generator_tree);
		//std::list<int> ttt;
		std::vector<int> tp(num_vertices(g));
		std::vector<int> tdis(num_vertices(g));
		//std::ofstream hhh2("E:\\hjs\\c++\\hh1\\pro2\\build\\Release\\td_tdis.data");

		for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
		{
			C1& pE = VC1[name[*ei]];
			int s = ei->m_source, t = ei->m_target;
			auto [e, f] = edge(s, t, generator_tree);


			if (is_tree[name_g[*ei]] == 0) {

				tdis[s] = 0;
				//breadth_first_search(subg, s, visitor(bfs_time_visitor(tp, tdis)));
				//int td = dis[s] + dis[t] - 2 * dis[gname[e]];

				//we should add d(s,t) in this term for weighted_grapg
				int td1 = dis[s] + dis[t];
				int td2 = td1 - 2 * dis[gname[e]];
				dis_lca[name_g[*ei]] = td1;
				/*			assert(td == tdis[t]);*/
						/*	hhh2 << format("%d %d\n") % td%tdis[t];*/
				mmap.insert({ td1,&pE });
				//ttt.push_back(s);
				//ttt.push_back(dis[s]);
				//ttt.push_back(t);
				//ttt.push_back(dis[t]);
				//ttt.push_back(gname[e]);
				//ttt.push_back(dis[gname[e]]);
			}

		}
		//std::vector<int> pF(ttt.begin(), ttt.end());
		//Matrix<int, Dynamic, Dynamic, RowMajor> F = Map<Matrix<int, Dynamic, Dynamic, RowMajor>>(pF.data(), pF.size() / 6, 6);
		//std::ofstream hhh("E:\\hjs\\c++\\hh1\\pro2\\build\\Release\\st_dis.data");
		//hhh << F;

		for (auto& [it, CC] : mmap) {
			CC->id = vid;
			idC1[vid++] = CC;
			//std::cout << format("non_tree vid:%d\n") % (vid - 1);
		}
		if (count != spanning_tree.size()) {
			std::cerr << "count != spanning_tree.size()\n";
		}

		count = 0;
		for (auto& it : VC1) {
			count += (it.generator == false ? 1 : 0);
		}
		if (count + 1 != VC0.size()) {
			std::cerr << "count + 1 != m_pMesh->vertices().size()\n";
		}


		count = 0;
		struct C2s {
			int id;
			C2* cc;
		};

		std::vector<C2*> fvec;
		fvec.reserve(VC2.size());
		{
			std::vector<C2s> tvec;
			tvec.reserve(VC2.size());
			int mmax = -1;
			for (auto& it : VC2) {
				mmax = -1;
				Dart_handle dc1 = it.da;
				mmax = std::max(VC1[mesh.info<1>(dc1)].id, mmax);
				dc1 = mesh.beta<1>(dc1);
				mmax = std::max(VC1[mesh.info<1>(dc1)].id, mmax);
				dc1 = mesh.beta<1>(dc1);
				mmax = std::max(VC1[mesh.info<1>(dc1)].id, mmax);
				dc1 = mesh.beta<1>(dc1);
				mmax = std::max(VC1[mesh.info<1>(dc1)].id, mmax);
				tvec.push_back({ mmax,&it });
			}
			std::sort(tvec.begin(), tvec.end(), [](C2s& c1, C2s& c2)
				{return c1.id < c2.id; });
			for (auto& it : tvec)
				fvec.push_back(it.cc);
		}
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

		//shuffle(fvec.begin(), fvec.end(), std::default_random_engine(seed));


		std::vector<std::list<int>> vl(VC2.size(), std::list<int>());
		std::vector<int> vlind(VC1.size() + 5, -1);
		for (int i = 0; i < VC2.size(); i++) {
			C2* tpF = fvec[i];
			auto& it2 = vl[i];
			Dart_handle dc1 = tpF->da;
			if ((VC1.data() + mesh.info<1>(dc1))->generator)
				it2.push_back(VC1[mesh.info<1>(dc1)].id);
			dc1 = mesh.beta<1>(dc1);
			if ((VC1.data() + mesh.info<1>(dc1))->generator)
				it2.push_back(VC1[mesh.info<1>(dc1)].id);
			dc1 = mesh.beta<1>(dc1);
			if ((VC1.data() + mesh.info<1>(dc1))->generator)
				it2.push_back(VC1[mesh.info<1>(dc1)].id);
			dc1 = mesh.beta<1>(dc1);
			if ((VC1.data() + mesh.info<1>(dc1))->generator)
				it2.push_back(VC1[mesh.info<1>(dc1)].id);
			it2.sort(std::greater<int>());
			int i1 = it2.size();
			/*it2.unique();
			assert(it2.size()==i1);*/
			
		}
		
		
		std::list<int>::iterator it1, it2;
		//Merge from begin of (list1) to (it1) of (it).
		//(it) will change, while (list1) not.
		//Return the first element after (it1) in (it).
		auto merge = []( std::list<int>& it,std::list<int>::iterator it1, std::list<int>& list1 )
			-> std::list<int>::iterator {
			bool edit = 1;
			std::list<int>::iterator ri;
			std::list<int>::iterator it2 = list1.begin();
			while (it1 != it.end() && it2 != list1.end()) {
				if (*it1 < *it2) {
					it.insert(it1, *it2);
					ri = (edit) ? std::prev(it1) : ri;
					edit = 0;
					it2++;
				}
				else if (*it1 > *it2) {
					ri = (edit) ? it1 : ri;
					edit = 0;
					it1++;
					
				}
				else {
					it2++;
					it1 = it.erase(it1);
					ri = (edit)?it1:ri;
				}
			}
			if (it1 == it.end()) {
				if (it2 != list1.end()) {
					it.insert(it1, *it2);
					ri = (edit) ? std::prev(it1) : ri;
					it2++;
					while (it2 != list1.end()) {
						it.insert(it1, *it2);
						it2++;
					}
				}
				
				
				
			}
			return ri;
		};
		int i3 = -1, i1,i2;
		for (auto& it : vl) {
			i3++;
			i1 = vlind[it.front()];
			it1 = it.begin();
			while (it1 != it.end()) {
				i1 = vlind[*it1];
				if (i1 >= 0)
					it1 = merge(it, it1,vl[i1]);
				else
					it1++;
			
				
			}
			if (!std::empty(it)) {
				vlind[it.front()] = i3;
				idC1[it.front()]->pair = 1;
			}


		}
		for (auto& [it, CC] : mmap) {
			assert(it > 0);
			i1 = vlind[CC->id];
			if (i1 < 0) {
				continue;
			}
			std::list<int>& it = vl[i1];
			it1 = std::next(it.begin());
			while (it1 != it.end()) {
				i1 = vlind[*it1];
				if (i1 >= 0)
					it1 = merge(it, it1, vl[i1]);
				else
					it1++;
			}
		}
		i2 = -1;  
		for (auto& it : VC1) {

			i2++;
			if (it.generator) {
				flag = 0;
				i1 = vlind[it.id];
				if (i1 >= 0) {
					if (vl[i1].size() > 1)
						flag = 1;
				}
				else
					flag = 1;
				if (flag==0)
					dis_lca[i2]=0;
					
			}
			//dis_lca[]
			//count += (it.generator) && (it.pair == -1) ? 1 : 0;
		}
		
	

	}





	{





		int count = 0;

		for (auto& it : VC1) {
			count += (it.generator) && (it.pair == -1) ? 1 : 0;
		}
		int g3 = VC0.size() - VC1.size() + VC2.size() - VC3.size();
		if (1 - g3 != count) {
			std::cerr << format("1 - g3(%d) != count(%d)\n") % g3 % count;
		}
		/*assert(count == 2);*/
		std::cout << format("loop num:%d\n") % count;
		if (count != 2) {
			/*std::cerr << "count != 2\n";*/
		}

	}

	{

		class mycount {
		public:
			mycount(std::vector<C1*>& ff) :i1(0), i2(0) {
				for (auto& it : ff)
					eset.insert(it);

			};
			void count(C1* it) {
				if (eset.find(it) != eset.end()) {
					i1++;
				}
				else {
					i2++;
				}
			}
			void reset() {
				i1 = 0;
				i2 = 0;
			}
			void print() {
				std::cerr << format("total edges:%5d,inter edges:%5d,boundary edges:%5d\n") % (i1 + i2) % i2 % i1;
			}

		private:
			std::unordered_set<C1*> eset;
			int i1, i2;

		};


		int N = VC0.size();
		std::vector<int> pre(N, 0);
		std::vector<Edge> preedge(N, Edge());
		std::vector<bool> visited(N, 0);
		if (subg.m_edges.size() + 1 != subg.m_vertices.size()) {
			std::cerr << "subg is wrong\n";
		}


		std::vector<C1*> evec;


		int count = 0, tcount = 0;
		for (auto& it : VC1) {
			if ((it.generator) && (it.pair == -1)) {
				evec.push_back(&it);
				is_tree[count] = 2;
				tcount++;
			}
			count++;
		}
		std::vector<std::list<Edge>> loop(tcount);
		Graph increg = subg;


		std::vector<int> tvec(num_edges(g));
		for (int i = 0; i < tvec.size(); i++)
			tvec[i] = i;
		std::sort(tvec.begin(), tvec.end(), [&dis_lca, &is_tree](const int& a, const int& b) {

			return std::make_pair(dis_lca[a], 2 - is_tree[a]) <
				std::make_pair(dis_lca[b], 2 - is_tree[b]);
			});

		



		int loopcount = 0;
		//is_tree=0,1,2
		for (int i = 0; i < num_edges(g); i++) {

			int ni = tvec[i];
			int aa = is_tree[ni];
			if (is_tree[ni] == 1) {

				continue;
			}
			auto it = VC1.data() + ni;

			const int id1 = mesh.info<0>(it->da);
			const int id2 = mesh.info<0>(mesh.beta<1>(it->da));

			if (is_tree[ni] == 0) {
				add_edge(id1, id2, increg);
				continue;
			}
			auto& tel = loop[loopcount++];



			auto v1 = vertex(id1, subg);
			auto v2 = vertex(id2, subg);
			
			auto breadth_first_search = [&pre, &preedge,&visited](int s,int t, Graph& g) {
			


				

				std::queue<int> record;
				// Create a queue for BFS
				std::list<int> queue;

				// Mark the current node as visited and enqueue it
				visited[s] = true;
				record.push(s);
				pre[s] = s;
				preedge[s] = Edge();
				queue.push_back(s);


				int targ=-1;
				while (!queue.empty())
				{
					// Dequeue a vertex from queue and print it
					int u = queue.front();
					
					queue.pop_front();

					// Get all adjacent vertices of the dequeued
					// vertex s. If a adjacent has not been visited, 
					// then mark it visited and enqueue it
					for (auto& [ei, en] = out_edges(u, g);ei != en; ++ei)
					{
						 targ = target(*ei, g);
						if (!visited[targ])
						{
							visited[targ] = true;
							record.push(targ);
							pre[targ] = u;
							preedge[targ] = *ei;
							if (targ == t)
								break;
							queue.push_back(targ);
						}
					}
					if (targ == t)
						break;
				}
				while (!record.empty())
				{
					visited[record.front()] = 0;
					record.pop();
				}

				


				
			};

			



	
			breadth_first_search(id1, id2, increg);

			int i1, i2 = id2;





			Edge ee = preedge[i2];
			auto name = get(edge_name, g);



			std::vector<Edge> vec;
			auto gname = get(edge_name, generator_tree);
			//id2->id1 loop in tree
			while (ee != Edge()) {
				int s = ee.m_source, t = ee.m_target, eit;
				i2 = ee.m_source;
				tel.push_back(edge(t, i2, g).first);
				ee = preedge[i2];

			}
			//veccount do not calculate endpoint of it;
			//simplify loop
			tel.push_back(edge(id1, id2, g).first);



			//std::cout << format("pass %d times\n") % veccount[start];
			std::cout << format("one loop finished\n");
			/*cc.count(it);

			cc.print();*/
			/*loop.push_back(tel);*/
			add_edge(id1, id2, increg);
		}
		std::vector<int>dis;
		auto name = get(edge_name, g);
		
		dis.reserve(loop.size());
		loopcount = 0;
		//count repeat loop
		for (auto& vec : loop) {

			
			loopcount++;
			dis.push_back(vec.size());





		}
		std::vector<int>index_dis;
		index_dis.reserve(dis.size());





		for (int i = 0; i < dis.size(); i++) {
			/*if(index_ind[i])*/
			index_dis.push_back(i);

		}
		std::cout << format("loop:%d real loop%d\n") % dis.size() % index_dis.size();
		//sort loop by length
		std::sort(index_dis.begin(), index_dis.end(),
			[&dis](int a, int b) {return dis[a] > dis[b]; });
		for (auto& it : index_dis) {
			std::cout << format("%3d ") % dis[it];
		}
		std::cout <<"\n";
		/*{
	for (auto& ei : gset) {
		auto[a, b] = ei;
		vc1file << format("%d %d\n") % a%b;

	}
}*/
		if (output_file != "") {
			system(("mkdir " + output_file).c_str());
		}


		for (int i = 0; i < std::max(std::size_t(0), index_dis.size()); i++) {
			std::ofstream vc1file;
			if (output_file == "") {
				vc1file = std::ofstream((format("loop%d.data") % i).str());
			}
			else {
				vc1file = std::ofstream(output_file + (format("\\loop%d.data") % i).str());

			}
			vc1file << format("%1%\n") % dis[index_dis[i]];
			auto& vec = loop[index_dis[i]];
			/*cc.reset();*/
			for (auto& it : vec) {
				VC1[name[it]].sharp = true;
				vc1file << format("%1%\n") % mesh.point(VC0[it.m_source].da);

				/*	cc.count(it);*/
			}
			/*cc.print();*/
		}

	}

}
