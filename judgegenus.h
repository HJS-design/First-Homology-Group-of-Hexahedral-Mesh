#pragma once
#include <boost/format.hpp>
#include "mesh.h"
#include "geometry.h"
#include "readmesh.cpp"
#include <boost/functional/hash.hpp>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <unordered_set>
#include <tuple>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
using namespace Eigen;
class judgegenus {
	typedef typename mymesh<geometry> tmesh;

public:

	judgegenus() {};

	void draw();
	//this count is not true
	void read();
	//use dart tu count
	void count_use_dart(std::string filename);
	void viewloop();
protected:
	std::vector<std::pair<int, std::array<int, 3>>> plist;
	std::vector<std::array<int, 4>> vlist;
	std::vector<int> g2vec;
	std::vector<int> g3vec;
	igl::opengl::glfw::Viewer viewer;
	Matrix<double, Dynamic, Dynamic, RowMajor> V;
	Matrix<int, Dynamic, Dynamic, RowMajor> F;
	Matrix<int, Dynamic, Dynamic, RowMajor> E;
	int bias_i = -1;
	int countbias_i=0;
};


inline void judgegenus::draw()
{


	tmesh mesh;
	readmesh hh(mesh);

	//hh.read(plist, std::vector<std::array<int, 4>>(vlist.begin(), vlist.begin() + 12 * (bias_i + 1)));
	hh.read_manifold_info(plist, vlist);

	hh.get_view_matrix(V, F, E);
}

void judgegenus::read()
{
	std::ifstream hh("hh.graph");
	//std::istream& hh = std::cin;
	int i1, i2, i3, i4, i5, p_all = 0, vn = 0, pn = 0, en = 0, fn = 0, subp = 0, sube = 0, subf = 0;
	double d1, d2, d3;
	using namespace boost;
	using namespace Eigen;
	typedef std::tuple<Vector3i, std::array<int, 6>> cellinfo;
	hh >> p_all;
	g2vec.reserve(p_all);
	g3vec.reserve(p_all);

	Array3i bias;
	Array3i sizenum;

	class edge {
	public:
		edge(unsigned int oa, unsigned int ob) {
			if (oa < ob) {
				a = oa, b = ob;
			}
			else {
				b = oa, a = ob;
			}
		}
		unsigned int a, b;
	};

	struct hash_name {
		size_t operator()(const edge& p) const {
			return h_pair(std::make_pair(p.a, p.b));
		}
		boost::hash<std::pair<int, int>> h_pair;
	};

	struct equal {
		bool operator()(const edge& p1, const edge& p2) const {
			return (std::tie(p1.a, p1.b) == std::tie(p2.a, p2.b));
		}
	};

	std::unordered_map<int, cellinfo> pmap;

	{
		std::unordered_map<edge, char, hash_name, equal> nemap;
		std::unordered_set<edge, hash_name, equal> teset;

		std::unordered_set<edge, hash_name, equal> fset;



		std::unordered_map<unsigned int, char> npmap;
		std::unordered_set<unsigned int> tpset;
		auto getindex = [&](Array3i v) {
			assert((v > 0).all());
			assert((v < sizenum).all());
			return (v * bias).sum(); };
		auto addedge = [&](int a, int b) {
			auto& [point_a, neigh_a] = pmap[a];
			auto& [point_b, neigh_b] = pmap[b];
			ArrayXXi facep(4, 3);
			Vector3i np = point_a + point_b;
			if (point_a(0) != point_b(0)) {
				facep.row(0) = np + Vector3i{ 0, -1, 1 };
				facep.row(1) = np + Vector3i{ 0, 1, 1 };
				facep.row(2) = np + Vector3i{ 0, 1, -1 };
				facep.row(3) = np + Vector3i{ 0, -1, -1 };

			}
			else if (point_a(1) != point_b(1)) {
				facep.row(0) = np + Vector3i{ -1, 0, 1 };
				facep.row(1) = np + Vector3i{ 1, 0, 1 };
				facep.row(2) = np + Vector3i{ 1, 0, -1 };
				facep.row(3) = np + Vector3i{ -1, 0, -1 };
			}
			else {
				facep.row(0) = np + Vector3i{ -1, 1, 0 };
				facep.row(1) = np + Vector3i{ 1, 1, 0 };
				facep.row(2) = np + Vector3i{ 1, -1, 0 };
				facep.row(3) = np + Vector3i{ -1, -1, 0 };
			}
			std::array<int, 4> index;
			for (int i = 0; i < 4; i++) {
				index[i] = getindex(facep.row(i));
			}
			for (int i = 0; i < 4; i++) {
				unsigned int id1 = index[i];

				edge e1(id1, index[(i + 1) % 4]);
				auto [it, f1] = teset.insert(e1);

				auto [it1, f] = nemap.insert({ e1 ,1 });
				if (!f) {
					if (f1)en--;
					(it1->second)++;
					if (it1->second == 8)sube++;
				}

				{
					auto [it, f1] = tpset.insert(id1);
					auto [it1, f] = npmap.insert({ id1,1 });
					if (!f) {
						if (f1)pn--;
						(it1->second)++;
						if (it1->second == 24)subp++;
					}

				}
			}
		};
		std::vector<int> pvec;
		pvec.reserve(p_all);
		//read file, get bias
		{
			std::array<int, 6> tt;
			ArrayXXi mm(3, p_all);

			ArrayXXi t2(3, 9);
			for (int i = 0; i < p_all; i++) {
				hh >> i1 >> i2 >> i3 >> i4;


				pvec.push_back(i1);
				for (auto& it : tt) hh >> it;
				Vector3i vv(i2, i3, i4);

				mm.col(i) = vv;
				pmap.insert({ i1, {vv, tt} });

			}

			Array3i mmmin = mm.rowwise().minCoeff();
			Array3i mmmax = mm.rowwise().maxCoeff();
			for (auto& [id, cinfo] : pmap) {
				auto& [vv, tt] = cinfo;
				vv = vv.array() - mmmin + 1;
			};
			mmmax = (mmmax - mmmin + 1) * 2 + 3;
			sizenum = mmmax;
			/*      std::cout << mmmax <<"\n----------\n";*/
			bias = Array3i(1, mmmax(0), mmmax(0) * mmmax(1));
			//std::vector<int> hh;
			//for (int i = 0; i < mmmax(2); i++)
			//	for (int j = 0; j < mmmax(1); j++)
			//		for (int k = 0; k < mmmax(0); k++)
			//			hh.push_back(getindex(Array3i(k, j, i)));
			//for (int i = 0; i < hh.size() - 1; i++)
			//	assert(hh[i] + 1 == hh[i + 1]);

		}
		auto g3 = format("-<%d>--p:%d e:%d f:%d v:%d g3:%d----\n");
		auto g2 = format("-    --p:%d e:%d f:%d g2:%d----\n\n");
		//plist,vlist to mesh
		{
			std::unordered_map<int, Array3i> t1;
			ArrayXXi t2(3, 9);
			for (auto& v : pvec) {

				auto& cinfo = pmap[v];

				Vector3i vv = std::get<0>(cinfo);
				{
					auto addedge = [&](Vector3i point_a, Vector3i point_b) {

						ArrayXXi facep(4, 3);
						Vector3i np = point_a + point_b;
						if (point_a(0) != point_b(0)) {
							facep.row(0) = np + Vector3i{ 0, -1, 1 };
							facep.row(1) = np + Vector3i{ 0, 1, 1 };
							facep.row(2) = np + Vector3i{ 0, 1, -1 };
							facep.row(3) = np + Vector3i{ 0, -1, -1 };

						}
						else if (point_a(1) != point_b(1)) {
							facep.row(0) = np + Vector3i{ -1, 0, 1 };
							facep.row(1) = np + Vector3i{ 1, 0, 1 };
							facep.row(2) = np + Vector3i{ 1, 0, -1 };
							facep.row(3) = np + Vector3i{ -1, 0, -1 };
						}
						else {
							facep.row(0) = np + Vector3i{ -1, 1, 0 };
							facep.row(1) = np + Vector3i{ 1, 1, 0 };
							facep.row(2) = np + Vector3i{ 1, -1, 0 };
							facep.row(3) = np + Vector3i{ -1, -1, 0 };
						}
						point_a *= 2;
						vlist.push_back({ getindex(point_a),getindex(facep.row(0)),getindex(facep.row(1)),getindex(facep.row(2)) });
						vlist.push_back({ getindex(point_a),getindex(facep.row(2)),getindex(facep.row(3)),getindex(facep.row(0)) });

					};
					t2.col(0) = vv + Vector3i{ 1, 0, 0 };
					t2.col(1) = vv + Vector3i{ -1, 0, 0 };
					t2.col(2) = vv + Vector3i{ 0,  1, 0 };
					t2.col(3) = vv + Vector3i{ 0, -1, 0 };
					t2.col(4) = vv + Vector3i{ 0,0,  1 };
					t2.col(5) = vv + Vector3i{ 0,0, -1 };
					for (int i = 0; i < 6; i++)
						addedge(vv, t2.col(i));

				}
				{
					vv = vv * 2;
					t2.col(0) = vv;
					t2.col(1) = vv + Vector3i{ -1, -1, -1 };
					t2.col(2) = vv + Vector3i{ -1, -1,  1 };
					t2.col(3) = vv + Vector3i{ -1,  1,  1 };
					t2.col(4) = vv + Vector3i{ -1,  1, -1 };
					t2.col(5) = vv + Vector3i{ 1, -1, -1 };
					t2.col(6) = vv + Vector3i{ 1, -1,  1 };
					t2.col(7) = vv + Vector3i{ 1,  1,  1 };
					t2.col(8) = vv + Vector3i{ 1,  1, -1 };
					/*std::cout << t2 << '\n';*/
					for (int i = 0; i < 9; i++)
						t1.insert({ getindex(t2.col(i)),t2.col(i) });
				}
			}
			for (auto& [id, point] : t1) {
				plist.push_back({ id,{point(0),point(1),point(2)} });
			}
		}



		int count = 0;
		for (auto& v : pvec) {
			auto& cinfo = pmap[v];
			tpset.clear();
			teset.clear();
			fn += 6;
			pn += 8;
			en += 12;
			vn++;
			auto& [aa, info] = cinfo;
			for (auto& it : info) {
				if (it == -1)
					continue;
				auto [a, b] = fset.insert(edge(v, it));

				if (!b) fn--, subf++;
				addedge(it, v);
			}
			int ng2 = pn - en + fn - subp + sube - subf;
			int ng3 = pn - en + fn - vn;
			g2vec.push_back(ng2);
			g3vec.push_back(ng3);
			std::cout << g3 % (++count) % pn % en % fn % vn % (ng3);
			std::cout << g2 % (pn - subp) % (en - sube) % (fn - subf) % (ng2);
		}
	}

	//for (auto& it : g2vec)
	//	std::cout << it << ' ';
	std::cout << "\n\n";
	for (auto& it : g3vec)
		std::cout << it << ' ';
}

inline void judgegenus::count_use_dart(std::string filename)
{
	std::ifstream hh(filename);
	//std::istream& hh = std::cin;
	int i1, i2, i3, i4, i5, p_all = 0, vn = 0, pn = 0, en = 0, fn = 0, subp = 0, sube = 0, subf = 0;
	double d1, d2, d3;
	using namespace boost;
	using namespace Eigen;
	typedef std::tuple<Vector3i, std::array<int, 6>> cellinfo;
	hh >> p_all;

	g3vec.reserve(p_all);

	Array3i bias;
	Array3i sizenum;

	class edge {
	public:
		edge(unsigned int oa, unsigned int ob) {
			if (oa < ob) {
				a = oa, b = ob;
			}
			else {
				b = oa, a = ob;
			}
		}
		unsigned int a, b;
	};

	struct hash_name {
		size_t operator()(const edge& p) const {
			return h_pair(std::make_pair(p.a, p.b));
		}
		boost::hash<std::pair<int, int>> h_pair;
	};

	struct equal {
		bool operator()(const edge& p1, const edge& p2) const {
			return (std::tie(p1.a, p1.b) == std::tie(p2.a, p2.b));
		}
	};

	std::unordered_map<int, cellinfo> pmap;

	{

		auto getindex = [&](Array3i v) {
			assert((v > 0).all());
			assert((v < sizenum).all());
			return (v * bias).sum(); };

		std::vector<int> pvec;
		pvec.reserve(p_all);
		//read file, get bias
		{
			std::array<int, 6> tt;
			ArrayXXi mm(3, p_all);

			ArrayXXi t2(3, 9);
			for (int i = 0; i < p_all; i++) {
				hh >> i1 >> i2 >> i3 >> i4;


				pvec.push_back(i1);
				for (auto& it : tt) hh >> it;
				Vector3i vv(i2, i3, i4);

				mm.col(i) = vv;
				pmap.insert({ i1, {vv, tt} });

			}

			Array3i mmmin = mm.rowwise().minCoeff();
			Array3i mmmax = mm.rowwise().maxCoeff();
			for (auto& [id, cinfo] : pmap) {
				auto& [vv, tt] = cinfo;
				vv = vv.array() - mmmin + 1;
			};
			mmmax = (mmmax - mmmin + 1) * 2 + 3;
			sizenum = mmmax;
			/*      std::cout << mmmax <<"\n----------\n";*/
			bias = Array3i(1, mmmax(0), mmmax(0) * mmmax(1));
			//std::vector<int> hh;
			//for (int i = 0; i < mmmax(2); i++)
			//	for (int j = 0; j < mmmax(1); j++)
			//		for (int k = 0; k < mmmax(0); k++)
			//			hh.push_back(getindex(Array3i(k, j, i)));
			//for (int i = 0; i < hh.size() - 1; i++)
			//	assert(hh[i] + 1 == hh[i + 1]);

		}
		auto g3 = format("-<%d>--p:%d e:%d f:%d v:%d g3:%d----\n");
		auto g2 = format("-    --p:%d e:%d f:%d g2:%d----\n\n");
		//plist,vlist to mesh
		{
			std::unordered_map<int, Array3i> t1;
			ArrayXXi t2(3, 9);
			for (auto& v : pvec) {

				auto& cinfo = pmap[v];

				Vector3i vv = std::get<0>(cinfo);
				{
					auto addedge = [&](Vector3i point_a, Vector3i point_b) {

						ArrayXXi facep(4, 3);
						Vector3i np = point_a + point_b;
						if (point_a(0) != point_b(0)) {
							facep.row(0) = np + Vector3i{ 0, -1, 1 };
							facep.row(1) = np + Vector3i{ 0, 1, 1 };
							facep.row(2) = np + Vector3i{ 0, 1, -1 };
							facep.row(3) = np + Vector3i{ 0, -1, -1 };

						}
						else if (point_a(1) != point_b(1)) {
							facep.row(0) = np + Vector3i{ -1, 0, 1 };
							facep.row(1) = np + Vector3i{ 1, 0, 1 };
							facep.row(2) = np + Vector3i{ 1, 0, -1 };
							facep.row(3) = np + Vector3i{ -1, 0, -1 };
						}
						else {
							facep.row(0) = np + Vector3i{ -1, 1, 0 };
							facep.row(1) = np + Vector3i{ 1, 1, 0 };
							facep.row(2) = np + Vector3i{ 1, -1, 0 };
							facep.row(3) = np + Vector3i{ -1, -1, 0 };
						}
						point_a *= 2;
						vlist.push_back({ getindex(point_a),getindex(facep.row(0)),getindex(facep.row(1)),getindex(facep.row(2)) });
						vlist.push_back({ getindex(point_a),getindex(facep.row(2)),getindex(facep.row(3)),getindex(facep.row(0)) });

					};
					t2.col(0) = vv + Vector3i{ 1, 0, 0 };
					t2.col(1) = vv + Vector3i{ -1, 0, 0 };
					t2.col(2) = vv + Vector3i{ 0,  1, 0 };
					t2.col(3) = vv + Vector3i{ 0, -1, 0 };
					t2.col(4) = vv + Vector3i{ 0,0,  1 };
					t2.col(5) = vv + Vector3i{ 0,0, -1 };
					for (int i = 0; i < 6; i++)
						addedge(vv, t2.col(i));

				}
				{
					vv = vv * 2;
					t2.col(0) = vv;
					t2.col(1) = vv + Vector3i{ -1, -1, -1 };
					t2.col(2) = vv + Vector3i{ -1, -1,  1 };
					t2.col(3) = vv + Vector3i{ -1,  1,  1 };
					t2.col(4) = vv + Vector3i{ -1,  1, -1 };
					t2.col(5) = vv + Vector3i{ 1, -1, -1 };
					t2.col(6) = vv + Vector3i{ 1, -1,  1 };
					t2.col(7) = vv + Vector3i{ 1,  1,  1 };
					t2.col(8) = vv + Vector3i{ 1,  1, -1 };
					/*std::cout << t2 << '\n';*/
					for (int i = 0; i < 9; i++)
						t1.insert({ getindex(t2.col(i)),t2.col(i) });
				}
			}
			for (auto& [id, point] : t1) {
				plist.push_back({ id,{point(0),point(1),point(2)} });
			}
		}

		{
			std::ofstream log("debug.log");
			std::vector<std::array<int, 4>> trace_res;
			{
				tmesh mmesh;
				readmesh hh(mmesh);
				hh.read_with_count(plist, vlist, trace_res);
			}
			{
				tmesh mmesh;
				readmesh hh(mmesh);
				hh.read_non_manifold_info(plist, vlist);
				loop hh1(mmesh);
				drawdart mydraw(mmesh);
				mydraw.draw();
			}
			
			int count = 0;
			for (int i = 11; i < trace_res.size(); i += 12) {

				
				

				auto& res = trace_res[i];
	
				int ng3 = res[0] - res[1] + res[2] - res[3];
				log << g3 % (++count) % res[0] % res[1] % res[2] % res[3] % (ng3);
				//std::cout << g3 % (count) % res[0] % res[1] % res[2] % res[3] % (ng3);
				g3vec.push_back(ng3);
			}
		}
		assert(g3vec.size()==p_all);
		/*std::ofstream log("hh.log");
		for (auto& it : g3vec)
			log << it << ' ';*/


		//std::ifstream log("hh.log");
		//int it;
		//while (log >> it) {
		//	g3vec.push_back(it);
		//}
	}


}

inline void judgegenus::viewloop()
{
	auto key_down = [this](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)->bool
	{

		
		if (key == '1')
		{
			bool f = 0;
			for (int i = (bias_i == -1 ? 1 : bias_i + 1); i < g3vec.size(); i++) {
				if (g3vec[i] != g3vec[i - 1]) {
					bias_i = i;
					f = 1;
					countbias_i++;
					break;
				}
			}
			std::cout << boost::format("<%d>\n") % countbias_i;
			if (!f)
				return false;
			/*bias_i=24;*/
			draw();
			// Clear should be called before drawing the mesh
			viewer.data().clear();
			// Draw_mesh creates or updates the vertices and faces of the displayed mesh.
			// If a mesh is already displayed, draw_mesh returns an error if the given V and
			// F have size different than the current ones
			viewer.data().set_mesh(V, F);
			viewer.core().align_camera_center(V, F);
		};


		return false;
	};
	viewer.callback_key_down = key_down;
	viewer.launch();
}
