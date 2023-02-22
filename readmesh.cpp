#pragma once
#include "readmesh.h"


template<typename eproperty>
class bfs_time_visitor : public boost::default_bfs_visitor
{


public:
	bfs_time_visitor(std::vector<char>& out, eproperty& outname, std::vector<char>& outori) :
		cluster(out), f(1), name(outname), ori(outori) {}
	template < typename Edge, typename Graph >
	void tree_edge(Edge e, const Graph& g) const
	{
		int u = source(e, g), v = target(e, g), wu, wv;

		boost::tie(wu, wv) = name[e];
		wu %= 2; wv %= 2;
		ori[v] = ori[u] ^ (wu == wv);
		cluster[v] = 3 - cluster[u];
	}
	template < typename Edge, typename Graph >
	void non_tree_edge(Edge e, const Graph& g)
	{
		if (cluster[target(e, g)] == cluster[source(e, g)])
			f = 0;

	}
	eproperty& name;
	bool f;
	std::vector<char>& cluster;
	std::vector<char>& ori;
};
template<class LCC_3>
void readmesh<LCC_3>::disp(Dart_handle d)
{
	auto d2 = mesh.beta<1>(d);
	std::cout << boost::format("<%1%>,(%2%)-><%3%>,(%4%)\n") % mesh.info<0>(d) % mesh.point(d)
		% mesh.info<0>(d2) % mesh.point(d2);
}
template<class LCC_3>
void readmesh<LCC_3>::disp_loop(Dart_handle d1, Dart_handle d2)
{
	assert(mesh.info<0>(d1) == mesh.info<0>(mesh.beta<1>(d2)));
	d1 = mesh.beta<1>(d1); d2 = mesh.beta<0>(d2);
	assert(mesh.info<0>(d1) == mesh.info<0>(mesh.beta<1>(d2)));
	d1 = mesh.beta<1>(d1); d2 = mesh.beta<0>(d2);
	assert(mesh.info<0>(d1) == mesh.info<0>(mesh.beta<1>(d2)));

}

template<class mymesh>
inline void readmesh<mymesh>::disp_one_cell()
{
	std::vector<unsigned int> cells(4);
	for (unsigned int i = 0; i <= 3; ++i)
	{
		cells[i] = i;
	}


	std::vector<unsigned int> res = mesh.count_cells(cells);

	/*std::vector<unsigned int> res2(4, 0);
	{
		const int i = 0;
		for (auto a1 = mesh.one_dart_per_cell<i>().begin(), a2 = mesh.one_dart_per_cell<i>().end();
			a1 != a2; a1++) {
			res2[i]++;
		}
	}

	{
		const int i = 1;
		for (auto a1 = mesh.one_dart_per_cell<i>().begin(), a2 = mesh.one_dart_per_cell<i>().end();
			a1 != a2; a1++) {
			res2[i]++;
		}
	}
	{
		const int i = 2;
		for (auto a1 = mesh.one_dart_per_cell<i>().begin(), a2 = mesh.one_dart_per_cell<i>().end();
			a1 != a2; a1++) {
			res2[i]++;
		}
	}
	{
		const int i = 3;
		for (auto a1 = mesh.one_dart_per_cell<i>().begin(), a2 = mesh.one_dart_per_cell<i>().end();
			a1 != a2; a1++) {
			res2[i]++;
		}
	}
	assert(VC0.size() == res2[0]);
	assert(VC1.size() == res2[1]);
	assert(VC2.size() == res2[2]);
	assert(VC3.size() == res2[3]);*/

	assert(VC0.size() == res[0]);
	assert(VC1.size() == res[1]);
	assert(VC2.size() == res[2]);
	assert(VC3.size() == res[3]);

	{
		auto& VC = VC0;
		const int dimension = 0;
		for (int i = 0; i < VC.size(); i++) {
			auto& it = VC[i].da;
			for (auto dc1 = mesh.darts_of_cell<dimension>(it).begin(), dc2 = mesh.darts_of_cell<dimension>(it).end();
				dc1 != dc2; dc1++) {
				assert(mesh.info<dimension>(dc1) == i);
			}
		};
	}
	{
		auto& VC = VC1;
		const int dimension = 1;
		for (int i = 0; i < VC.size(); i++) {
			auto& it = VC[i].da;
			for (auto dc1 = mesh.darts_of_cell<dimension>(it).begin(), dc2 = mesh.darts_of_cell<dimension>(it).end();
				dc1 != dc2; dc1++) {
				assert(mesh.info<dimension>(dc1) == i);
			}
		};
	}
	{
		auto& VC = VC2;
		const int dimension = 2;
		for (int i = 0; i < VC.size(); i++) {
			auto& it = VC[i].da;
			for (auto dc1 = mesh.darts_of_cell<dimension>(it).begin(), dc2 = mesh.darts_of_cell<dimension>(it).end();
				dc1 != dc2; dc1++) {
				auto d = mesh.info<dimension>(dc1);
				assert(d == i);
			}
		};
	}
	{
		auto& VC = VC3;
		const int dimension = 3;
		for (int i = 0; i < VC.size(); i++) {
			auto& it = VC[i].da;
			for (auto dc1 = mesh.darts_of_cell<dimension>(it).begin(), dc2 = mesh.darts_of_cell<dimension>(it).end();
				dc1 != dc2; dc1++) {
				assert(mesh.info<dimension>(dc1) == i);
			}
		};
	}

}

template<class LCC_3>
inline void readmesh<LCC_3>::read_t_file(std::string filename)
{
	std::ifstream hh(filename);
	int i1, i2, i3, i4, i5;
	double d1, d2, d3;


	std::list<std::pair<int, std::array<double, 3>>> plist;
	std::list<std::array<int, 4>> vlist;
	while (hh.good()) {
		std::string s1;
		hh >> s1;
		if (s1 == "Vertex") {
			hh >> i1 >> d1 >> d2 >> d3;
			plist.push_back({ i1,{d1,d2,d3} });

		}
		else if (s1 == "Tet") {

			hh >> i5 >> i1 >> i2 >> i3 >> i4;
			vlist.push_back({ i1,i2,i3,i4 });

		}
	}
	read_manifold_info(plist, vlist);







}
struct Sum_functor_min
{


	template<class Cell_attribute>
	void operator()(Cell_attribute& ca1, Cell_attribute& ca2)
	{
		//std::cout << format("%1% %2% %3%\n") % (typeid(ca1.info()).name()) % (ca1.info()) % (ca2.info());


		ca1.info() = std::min(ca1.info(), ca2.info());
		/*std::cout << "i am in\n";*/





	}


};
template<class mymesh>
void readmesh<mymesh>::read_hex_file_non_manifold_info(std::string name)
{
	std::string filename(name);
	typedef LCC_3::Vector vec;
	std::list<std::array<int, 7>> vlist;
	std::unordered_map<int, std::array<Dart_handle, 6>> darta;
	

	mesh.onmerge_functor<0>() = Sum_functor_min();
	mesh.onmerge_functor<1>() = Sum_functor_min();
	mesh.onmerge_functor<2>() = Sum_functor_min();
	mesh.onmerge_functor<3>() = Sum_functor_min();
	{
		std::array<int, 7> tt;
		int npoint;
		std::ifstream file(filename);
		file >> npoint;
		for (int i = 0; i < npoint; i++) {
			for (auto&it : tt)
				file >> it;
			vlist.push_back(tt);
		}
	}
	{
		Point p, q;
		Dart_handle d, d1, d2;
		for (auto&it : vlist) {
			p = Point(2 * it[1], 2 * it[2], 2 * it[3]);
			d = mesh.make_hexahedron(p + vec(-1, -1, -1), p + vec(1, -1, -1),
				p + vec(1, 1, -1), p + vec(-1, 1, -1),
				p + vec(-1, 1, 1), p + vec(-1, -1, 1),
				p + vec(1, -1, 1), p + vec(1, 1, 1));
			d1 = mesh.beta<1>(d);
			d2 = mesh.beta<2>(d);
			darta.insert({ it[0],{mesh.beta<1,2>(d1),mesh.beta<0,0,2>(d2),
				mesh.beta<2>(d1), (d2),
				(d), mesh.beta<0,2>(d)} });
		}
	}
	{
		std::array<int, 4> res = { 0,0,0,0 };
		for (auto it = mesh.attributes<0>().begin(), itend = mesh.attributes<0>().end(); it != itend; ++it)
		{
			it->info() = res[0]++;
		}
		for (auto it = mesh.darts().begin(), de = mesh.darts().end(); it != de; it++) {

			{
				const int i = 1;
				if (mesh.attribute<i>(it) == NULL)
					mesh.set_attribute<i>(it, mesh.create_attribute<i>(res[i]++));
			}
			{
				const int i = 2;
				if (mesh.attribute<i>(it) == NULL)
					mesh.set_attribute<i>(it, mesh.create_attribute<i>(res[i]++));
			}
			{
				const int i = 3;
				if (mesh.attribute<i>(it) == NULL)
					mesh.set_attribute<i>(it, mesh.create_attribute<i>(res[i]++));
			}

		}
	}

	{
		int i1, i2;
		for (auto&it : vlist) {
			for (int i = 0; i < 3; i++) {
				i1 = it[i + 4];
				if (i1 != -1) {
					mesh.sew<3>(darta[it[0]][i], darta[i1][i + 3]);
				}

			}
		}
	}
	{
		std::array<int, 4> res = { 0,0,0,0 };


		{
			const int i = 0;
			for (auto it = mesh.attributes<i>().begin(), itend = mesh.attributes<i>().end(); it != itend; ++it)
			{
				it->info() = res[i]++;
			}
		}

		{
			const int i = 1;
			for (auto it = mesh.attributes<i>().begin(), itend = mesh.attributes<i>().end(); it != itend; ++it)
			{
				it->info() = res[i]++;
			}
		}
		{
			const int i = 2;
			for (auto it = mesh.attributes<i>().begin(), itend = mesh.attributes<i>().end(); it != itend; ++it)
			{
				it->info() = res[i]++;
			}
		}
		{
			const int i = 3;
			for (auto it = mesh.attributes<i>().begin(), itend = mesh.attributes<i>().end(); it != itend; ++it)
			{
				it->info() = res[i]++;
			}
		}

		std::cout << boost::format("count:%d %d %d %d genus:%d\n")
			% res[0] % res[1] % res[2] % res[3] % (res[0] - res[1] + res[2] - res[3]);

		VC0.resize(res[0]);
		VC1.resize(res[1]);
		VC2.resize(res[2]);
		VC3.resize(res[3]);

		for (auto it(mesh.darts().begin()),
			itend(mesh.darts().end()); it != itend; ++it)
		{
			VC0[mesh.info<0>(it)].da = it;
			VC1[mesh.info<1>(it)].da = it;
			VC2[mesh.info<2>(it)].da = it;
			VC3[mesh.info<3>(it)].da = it;
		}
	}

}



template<class mymesh>
inline void readmesh<mymesh>::get_view_matrix(Matrix<double, Dynamic, Dynamic, RowMajor>& V, Matrix<int, Dynamic, Dynamic, RowMajor>& F, Matrix<int, Dynamic, Dynamic, RowMajor>& E)
{



	std::array<int, 4> res = { VC0.size(),VC1.size() ,VC2.size() ,VC3.size() };


	{
		std::vector<double> pV;
		pV.reserve(res[0] * 3);
		std::vector<int> pF;
		pF.reserve(res[3] * 4);
		std::vector<std::array<double, 3>> tpV(res[0]);
		for (auto a1 = mesh.one_dart_per_cell<0>().begin(), a2 = mesh.one_dart_per_cell<0>().end();
			a1 != a2; a1++) {
			auto& p = mesh.point(a1);

			tpV[mesh.attribute<0>(a1)->info()] = { p.x(),p.y(),p.z() };
		}
		for (auto& it : tpV) {
			pV.push_back(it[0]);
			pV.push_back(it[1]);
			pV.push_back(it[2]);
		}

		int count = 0, cc1 = 0;
		for (auto a1 = mesh.one_dart_per_cell<2>().begin(), a2 = mesh.one_dart_per_cell<2>().end(); a1 != a2; a1++) {
			cc1++;
			if (!mesh.is_free<3>(a1))
				continue;

			count = 0;
			for (auto dc1 = mesh.one_dart_per_incident_cell<0, 2>(a1).begin(),
				dc2 = mesh.one_dart_per_incident_cell<0, 2>(a1).end();
				dc1 != dc2; dc1++) {
				count++;
				pF.push_back(mesh.attribute<0>(dc1)->info());
			}
			assert(count == 3);

		}
		/*assert(cc1 == 7);*/
		V = Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(pV.data(), res[0], 3);
		F = Map<Matrix<int, Dynamic, Dynamic, RowMajor>>(pF.data(), pF.size() / 3, 3);
		pF.clear();
		for (auto& it : VC1) {
			if (it.sharp) {
				pF.push_back(mesh.info<0>(it.da));
				auto tt = mesh.beta<1>(it.da);
				pF.push_back(mesh.info<0>(tt));
			}

		}
		E = Map<Matrix<int, Dynamic, Dynamic, RowMajor>>(pF.data(), pF.size() / 2, 2);

	}

}





template<class mymesh>
template<typename pclass, typename vclass>
void readmesh<mymesh>::read_graph_file(std::string filename, pclass & plist, vclass & vlist)
{
	//Graph file is the coordinate of CT image. 
	//We first convert min coordinate to [1,1,1], then all origin coordinate *2, 
	//then dual vertex is the other integer vertex.
	std::ifstream hh(filename);
	//std::istream& hh = std::cin;
	int i1, i2, i3, i4, i5, p_all = 0, vn = 0, pn = 0, en = 0, fn = 0, subp = 0, sube = 0, subf = 0;
	double d1, d2, d3;
	using namespace boost;
	using namespace Eigen;
	typedef std::tuple<Vector3i, std::array<int, 6>> cellinfo;
	hh >> p_all;



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
			for (auto&[id, cinfo] : pmap) {
				auto&[vv, tt] = cinfo;
				vv = vv.array() - mmmin + 1;
			};
			mmmax = (mmmax - mmmin + 1) * 2 + 3;
			sizenum = mmmax;
			/*      std::cout << mmmax <<"\n----------\n";*/
			bias = Array3i(1, mmmax(0), mmmax(0) * mmmax(1));


		}

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
			for (auto&[id, point] : t1) {
				plist.push_back({ id,{point(0),point(1),point(2)} });
			}
		}
		hh.close();
	}
}

template<class mymesh>
template<typename pclass, typename vclass>
void readmesh<mymesh>::read_graph_file_with_space(std::string filename, pclass & plist, vclass & vlist)
{
	//Graph file is the coordinate of CT image. 
	//We first convert min coordinate to [1,1,1], then all origin coordinate *2, 
	//then dual vertex is the other integer vertex.
	std::ifstream hh(filename);
	//std::istream& hh = std::cin;
	int i1, i2, i3, i4, i5, p_all = 0, vn = 0, pn = 0, en = 0, fn = 0, subp = 0, sube = 0, subf = 0;
	double d1, d2, d3;
	using namespace boost;
	using namespace Eigen;
	typedef std::tuple<Vector3i, std::array<int, 6>> cellinfo;
	Array3d space;
	hh >> p_all;
	hh >> space(0) >> space(1) >> space(2);
	/*std::cout << space;*/


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
	std::unordered_map<int, Array3d > opmap;
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
				Array3d vv1(i2, i3, i4);
				mm.col(i) = vv;
				pmap.insert({ i1, {vv, tt} });
				opmap.insert({ i1,vv1 });
			}

			Array3i mmmin = mm.rowwise().minCoeff();
			Array3i mmmax = mm.rowwise().maxCoeff();
			for (auto&[id, cinfo] : pmap) {
				auto&[vv, tt] = cinfo;
				vv = vv.array() - mmmin + 1;
			};
			mmmax = (mmmax - mmmin + 1) * 2 + 3;
			sizenum = mmmax;
			/*      std::cout << mmmax <<"\n----------\n";*/
			bias = Array3i(1, mmmax(0), mmmax(0) * mmmax(1));


		}

		{
			std::unordered_map<int, Array3d> t1;
			ArrayXXi t2(3, 9);
			ArrayXXd t3(3, 9);
			for (auto& v : pvec) {

				auto& cinfo = pmap[v];

				Vector3i vv = std::get<0>(cinfo);
				Array3d vv1 = opmap[v];
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


					vv1 = vv1 * space;
					t3.col(0) = vv1;
					t3.col(1) = vv1 + Array3d{ -0.5, -0.5, -0.5 }*space;
					t3.col(2) = vv1 + Array3d{ -0.5, -0.5,  0.5 }*space;
					t3.col(3) = vv1 + Array3d{ -0.5,  0.5,  0.5 }*space;
					t3.col(4) = vv1 + Array3d{ -0.5,  0.5, -0.5 }*space;
					t3.col(5) = vv1 + Array3d{ 0.5, -0.5, -0.5 }*space;
					t3.col(6) = vv1 + Array3d{ 0.5, -0.5,  0.5 }*space;
					t3.col(7) = vv1 + Array3d{ 0.5,  0.5,  0.5 }*space;
					t3.col(8) = vv1 + Array3d{ 0.5,  0.5, -0.5 }*space;
					for (int i = 0; i < 9; i++)
						t1.insert({ getindex(t2.col(i)),t3.col(i) });
				}
			}
			for (auto&[id, point] : t1) {
				plist.push_back({ id,{point(0),point(1),point(2)} });
			}
		}
		hh.close();
	}
}



template<class mymesh>
template <typename pclass, typename vclass>
inline void readmesh<mymesh>::read_manifold_info(pclass plist, vclass vlist)
{

	int vn_all = 0, pn_all = 0;

	using namespace boost;
	typedef std::array<int, 9> cellinfo;







	{
		int i1, i2, i3, i4, i5;
		double d1, d2, d3;
		std::unordered_map<int, int> count;
		std::array<int, 4> tt;
		std::unordered_map<int, Point > pmap;

		for (auto& it : plist) {
			auto&[i1, t1] = it;
			auto&[d1, d2, d3] = t1;
			pmap.insert({ i1,Point(d1, d2, d3) });
		}
		for (auto& it : vlist) {

			tt = it;
			for (auto& it : tt) {
				auto[t1, f] = count.insert({ it,pn_all });
				if (f) {
					pn_all++;
					pvec.push_back(pmap[it]);
				}
				it = t1->second;
			}


			std::sort(tt.begin(), tt.end());
			cell.push_back({ tt[0], tt[1], tt[2],3  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[1], tt[3],2  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[2], tt[3],1  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[1], tt[2], tt[3],0  ,vn_all++,tt[0], tt[1], tt[2] ,tt[3] });

		}

		std::cout << format("vertex are:%d\n") % count.size();

	}
	std::cout << boost::format("pn:%d vn:%d\n") % pn_all % vn_all;
	VC0 = std::vector<C0>(pn_all);
	boost::hash<std::tuple<int, int, int>> hash_tupple;
	auto hash_key = [&](int a)->std::size_t { return hash_tupple(std::make_tuple(cell[a][0], cell[a][1], cell[a][2])); };

	auto equal_key = [&](int a, int b) {
		return(std::make_tuple(cell[a][0], cell[a][1], cell[a][2]) ==
			std::make_tuple(cell[b][0], cell[b][1], cell[b][2]));
	};

	std::unordered_set<int, std::function<std::size_t(int)>,
		std::function<bool(int, int)>> uset_f(pn_all, hash_key, equal_key);


	std::vector<char> ori(vn_all, 0);
	{
		typedef std::pair< int, int > E;
		std::vector<E> edge_array;
		std::vector<E> weight;
		typedef adjacency_list< vecS, vecS, directedS, no_property, property<edge_color_t, std::pair<int, int>> > graph;
		typedef graph_traits<graph> GraphTraits;
		typename GraphTraits::edge_descriptor e;
		for (int i = 0; i < cell.size(); i++) {
			auto f = uset_f.insert(i);
			if (f.second)
				continue;
			cellinfo& c1 = cell[*(f.first)], &c2 = cell[i];
			edge_array.push_back({ c1[4],c2[4] });
			weight.push_back({ c1[3],c2[3] });
			edge_array.push_back({ c2[4],c1[4] });
			weight.push_back({ c2[3],c1[3] });

		}
		uset_f.clear();
		graph g(edge_array.begin(), edge_array.end(), weight.begin(), vn_all);
		auto name = get(edge_color, g);
		std::vector<char> cluster(vn_all, 0);
		bfs_time_visitor<decltype(name)> vis(cluster, name, ori);
		for (int i = 0; i < vn_all; i++) {
			if (cluster[i])
				continue;
			cluster[i] = 1;
			breadth_first_search(g, i, visitor(vis));
		}
		if (vis.f == 0)
			std::cerr << "vis.f == 0\n";
		for (auto& it : cluster) {
			if (it == 0)
				std::cerr << "it==0;\n";
		}

	}

	{
		std::vector<bool> is_vol_create(vn_all, 0);
		std::vector < std::array<Dart_handle, 4>> dartv(vn_all);

		std::unordered_map<std::tuple<int, int>, int, boost::hash<std::tuple<int, int>>> m1;
		std::unordered_map<std::tuple<int, int>, int, boost::hash<std::tuple<int, int>>>::iterator it;
		bool f;
		int en = 0, fn = 0, vn = 0;

		auto create_cell = [&, this](cellinfo& arr, int curvol) ->auto{
			if (ori[arr[4]]) {

				auto a1 = mesh.create_vertex_attribute(pvec[arr[5]]);
				a1->info() = arr[5];
				auto a2 = mesh.create_vertex_attribute(pvec[arr[6]]);
				a2->info() = arr[6];
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				a3->info() = arr[7];
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				a4->info() = arr[8];
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);
				is_vol_create[curvol] = 1;
				dartv[curvol] = { mesh.beta<1, 2,1>(d3), mesh.beta<0, 2,0>(d3), mesh.beta<2,1>(d3), mesh.beta<0>(d3) };



				return d3;

			}
			else {
				auto a1 = mesh.create_vertex_attribute(pvec[arr[6]]);
				a1->info() = arr[6];
				auto a2 = mesh.create_vertex_attribute(pvec[arr[5]]);
				a2->info() = arr[5];
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				a3->info() = arr[7];
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				a4->info() = arr[8];
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);

				is_vol_create[curvol] = 1;

				dartv[curvol] = { mesh.beta<0, 2,0>(d3), mesh.beta<1, 2,1>(d3), mesh.beta<2,0>(d3), mesh.beta<1>(d3) };

				return d3;
			}
		};
		std::vector<std::tuple<int, int>> info_1 = { {5,6},{6,7},{5,7},{5,8},{6,8},{7,8} };
		std::array<int, 4> info_0_1 = { 5,6,7,8 };
		std::array<int, 4> info_0_2 = { 6,5,7,8 };
		auto modify_info_01 = [&, this](auto d3, auto& arr) {
			if (ori[arr[4]]) {


				{
					//assign every edge to a info,
					//may remove m1.insert()
					std::array<decltype(d3), 6> dlist = { d3,mesh.beta<1>(d3),mesh.beta<0>(d3)
				,mesh.beta<2,1>(d3) ,mesh.beta<1,2,1>(d3) ,mesh.beta<0,2,1>(d3) };
					for (int i = 0; i < 6; i++) {
						tie(it, f) = m1.insert({ std::make_tuple(arr[std::get<0>(info_1[i])],arr[std::get<1>(info_1[i])]),en });
						if (f) {

							auto p = mesh.attribute<1>(dlist[i]);
							if (p == NULL) {
								p = mesh.create_attribute<1>();
								mesh.set_attribute<1>(dlist[i], p);
							}

							p->info() = en;
							C1 cc;
							cc.da = dlist[i];
							VC1.push_back(cc);
							en++;
						}
					}
				}
				{
					//assign every vertex to a info,
					const std::array<decltype(d3), 4> dlist = { d3,mesh.beta<1>(d3),mesh.beta<0>(d3)
				,mesh.beta<2,0>(d3) };

					for (int i = 0; i < 4; i++)
						VC0[arr[info_0_1[i]]].da = dlist[i];
				}


			}
			else {
				{
					std::array<decltype(d3), 6> dlist = { d3,mesh.beta<0>(d3),mesh.beta<1>(d3)
									,mesh.beta<2,0>(d3) ,mesh.beta<0,2,0>(d3) ,mesh.beta<1,2,0>(d3) };
					for (int i = 0; i < 6; i++) {
						tie(it, f) = m1.insert({ std::make_tuple(arr[std::get<0>(info_1[i])],arr[std::get<1>(info_1[i])]),en });
						if (f) {
							auto p = mesh.attribute<1>(dlist[i]);
							if (p == NULL) {
								p = mesh.create_attribute<1>();
								mesh.set_attribute<1>(dlist[i], p);
							}

							p->info() = en;
							C1 cc;
							cc.da = dlist[i];
							VC1.push_back(cc);
							en++;
						}
					}
				}

				{
					const std::array<decltype(d3), 4> dlist = { d3,mesh.beta<1>(d3),mesh.beta<0>(d3)
				,mesh.beta<2,0>(d3) };

					for (int i = 0; i < 4; i++)
						VC0[arr[info_0_2[i]]].da = dlist[i];
				}

			}


		};
		auto modify_info_2 = [&, this](auto d) {
			auto p = mesh.attribute<2>(d);
			if (p == NULL) {
				p = mesh.create_attribute<2>();
				mesh.set_attribute<2>(d, p);
			}

			p->info() = fn;
			C2 cc;
			cc.da = d;
			VC2.push_back(cc);
			//std::cout << boost::format("%d %d %d\n")%(VC2.size()-1)%fn%(p->info());
			fn++;
		};
		auto modify_info_3 = [&, this](auto d3) {

			auto p = mesh.attribute<3>(d3);
			if (p == NULL) {
				p = mesh.create_attribute<3>();
				mesh.set_attribute<3>(d3, p);
			}

			p->info() = vn;
			C3 cc;
			cc.da = d3;
			VC3.push_back(cc);
			vn++;
		};
		std::array<Dart_handle, 4> dl;
		std::array<int, 4> il;

		int curvol;
		Dart_handle d1, d2, d3;
		mesh.set_automatic_attributes_management(false);


		for (int i = 0; i < cell.size(); i++) {
			/*if ((i % 48) == 0) {
				std::vector<unsigned int> cells(1, 0);

				std::vector<unsigned int> res = mesh.count_cells(cells);
				int ii1 = i / 48;
				std::cout << format("<%d %d>\n") % (ii1) % (res[0] - ii1);
			}*/
			auto aa = uset_f.insert(i);


			if (aa.second) {
				auto& arr = cell[i];
				curvol = arr[4];

				if (!(is_vol_create[curvol])) {
					auto d3 = create_cell(arr, curvol);
					modify_info_01(d3, arr);
					modify_info_3(d3);

				}

				d1 = dartv[curvol][arr[3]];
				modify_info_2(d1);
				continue;
			}






			auto& arr = cell[i];
			curvol = arr[4];
			if (!(is_vol_create[curvol])) {

				auto d3 = create_cell(arr, curvol);
				modify_info_01(d3, arr);
				modify_info_3(d3);
			}

			d2 = dartv[curvol][arr[3]];

			auto id_d11 = *(aa.first);
			d1 = dartv[cell[id_d11][4]][cell[id_d11][3]];
			//std::cout << format("%d %d - %d %d\n") % curvol % arr[3] % cell[id_d11][4] % cell[id_d11][3];
			disp_loop(d1, d2);

			assert(mesh.is_sewable<3>(d1, d2));
			mesh.sew<3>(d1, d2);






		}


	}

	//display_map_and_attributes<0>();
	mesh.set_automatic_attributes_management(true);
	/*display_map_and_attributes<0>();*/
	disp_one_cell();
	mesh.display_characteristics(std::cout);
	//std::cout << format("\n%d %d %d %d\n")%countcell[0]% countcell[1] % countcell[2]% countcell[3];
	assert(mesh.is_valid());
	std::cout << format("\nis valid: %d\n") % mesh.is_valid();
}
template<class mymesh>
template<typename pclass, typename vclass>
void readmesh<mymesh>::read_without_info(pclass plist, vclass vlist)
{
	//number of vulume,number of point
	int vn_all = 0, pn_all = 0;



	//mesh.set_automatic_attributes_management(true);
	using namespace boost;
	typedef std::array<int, 9> cellinfo;







	{
		int i1, i2, i3, i4, i5;
		double d1, d2, d3;
		std::unordered_map<int, int> count;
		std::array<int, 4> tt;
		std::unordered_map<int, Point > pmap;

		for (auto& it : plist) {
			auto&[i1, t1] = it;
			auto&[d1, d2, d3] = t1;
			pmap.insert({ i1,Point(d1, d2, d3) });
		}
		for (auto& it : vlist) {

			tt = it;
			for (auto& it : tt) {
				auto[t1, f] = count.insert({ it,pn_all });
				if (f) {
					pn_all++;
					pvec.push_back(pmap[it]);
				}
				it = t1->second;
			}


			std::sort(tt.begin(), tt.end());
			cell.push_back({ tt[0], tt[1], tt[2],3  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[1], tt[3],2  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[2], tt[3],1  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[1], tt[2], tt[3],0  ,vn_all++,tt[0], tt[1], tt[2] ,tt[3] });

		}


	}

	boost::hash<std::tuple<int, int, int>> hash_tupple;
	auto hash_key = [&](int a)->std::size_t { return hash_tupple(std::make_tuple(cell[a][0], cell[a][1], cell[a][2])); };

	auto equal_key = [&](int a, int b) {
		return(std::make_tuple(cell[a][0], cell[a][1], cell[a][2]) ==
			std::make_tuple(cell[b][0], cell[b][1], cell[b][2]));
	};

	std::unordered_set<int, std::function<std::size_t(int)>,
		std::function<bool(int, int)>> uset_f(pn_all, hash_key, equal_key);


	std::vector<char> ori(vn_all, 0);
	{
		typedef std::pair< int, int > E;
		std::vector<E> edge_array;
		std::vector<E> weight;
		typedef adjacency_list< vecS, vecS, directedS, no_property, property<edge_color_t, std::pair<int, int>> > graph;
		typedef graph_traits<graph> GraphTraits;
		typename GraphTraits::edge_descriptor e;
		for (int i = 0; i < cell.size(); i++) {
			auto f = uset_f.insert(i);
			if (f.second)
				continue;
			cellinfo& c1 = cell[*(f.first)], &c2 = cell[i];
			edge_array.push_back({ c1[4],c2[4] });
			weight.push_back({ c1[3],c2[3] });
			edge_array.push_back({ c2[4],c1[4] });
			weight.push_back({ c2[3],c1[3] });

		}
		uset_f.clear();
		graph g(edge_array.begin(), edge_array.end(), weight.begin(), vn_all);
		auto name = get(edge_color, g);
		std::vector<char> cluster(vn_all, 0);
		bfs_time_visitor<decltype(name)> vis(cluster, name, ori);
		for (int i = 0; i < vn_all; i++) {
			if (cluster[i])
				continue;
			cluster[i] = 1;
			breadth_first_search(g, i, visitor(vis));
		}
		if (vis.f == 0)
			std::cerr << "vis.f == 0\n";
		for (auto& it : cluster) {
			if (it == 0)
				std::cerr << "it==0;\n";
		}

	}

	{
		std::vector<bool> is_vol_create(vn_all, 0);
		std::vector < std::array<Dart_handle, 4>> dartv(vn_all);

		//std::unordered_map<std::tuple<int, int>, int, boost::hash<std::tuple<int, int>>>::iterator it;
		bool f;
		int en = 0, fn = 0, vn = 0;

		auto create_cell = [&, this](cellinfo& arr, int curvol) ->auto{
			if (ori[arr[4]]) {

				auto a1 = mesh.create_vertex_attribute(pvec[arr[5]]);

				auto a2 = mesh.create_vertex_attribute(pvec[arr[6]]);
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);
				is_vol_create[curvol] = 1;
				dartv[curvol] = { mesh.beta<1, 2,1>(d3), mesh.beta<0, 2,0>(d3), mesh.beta<2,1>(d3), mesh.beta<0>(d3) };



				return d3;

			}
			else {
				auto a1 = mesh.create_vertex_attribute(pvec[arr[6]]);
				auto a2 = mesh.create_vertex_attribute(pvec[arr[5]]);
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);

				is_vol_create[curvol] = 1;

				dartv[curvol] = { mesh.beta<0, 2,0>(d3), mesh.beta<1, 2,1>(d3), mesh.beta<2,0>(d3), mesh.beta<1>(d3) };

				return d3;
			}
		};


		std::array<Dart_handle, 4> dl;
		std::array<int, 4> il;

		int curvol;
		Dart_handle d1, d2, d3;




		for (int i = 0; i < cell.size(); i++) {
			//std::cout << outfmt % i%cell.size();
			auto aa = uset_f.insert(i);


			if (aa.second) {
				auto& arr = cell[i];
				curvol = arr[4];

				if (!(is_vol_create[curvol])) {
					auto d3 = create_cell(arr, curvol);



				}

				d1 = dartv[curvol][arr[3]];


			}
			else {
				auto& arr = cell[i];
				curvol = arr[4];
				if (!(is_vol_create[curvol])) {

					auto d3 = create_cell(arr, curvol);


				}

				d2 = dartv[curvol][arr[3]];

				auto id_d11 = *(aa.first);
				d1 = dartv[cell[id_d11][4]][cell[id_d11][3]];
				//std::cout << format("%d %d - %d %d\n") % curvol % arr[3] % cell[id_d11][4] % cell[id_d11][3];




				assert(mesh.is_sewable<3>(d1, d2));
				mesh.sew<3>(d1, d2);
			}












		}


	}

	assert(trace_res.size() == vlist.size());

}
template<class mymesh>
template<typename pclass, typename vclass>
void readmesh<mymesh>::read_non_manifold_info(pclass plist, vclass vlist)
{
	//number of vulume,number of point
	int vn_all = 0, pn_all = 0;
	std::array<int, 4> res = { 0,0,0,0 };
	std::array<int, 4> subres = res;


	mesh.onmerge_functor<0>() = Sum_functor(subres[0]);
	mesh.onmerge_functor<1>() = Sum_functor(subres[1]);
	mesh.onmerge_functor<2>() = Sum_functor(subres[2]);
	mesh.onmerge_functor<3>() = Sum_functor(subres[3]);
	//mesh.set_automatic_attributes_management(true);
	using namespace boost;
	typedef std::array<int, 9> cellinfo;







	{
		int i1, i2, i3, i4, i5;
		double d1, d2, d3;
		std::unordered_map<int, int> count;
		std::array<int, 4> tt;
		std::unordered_map<int, Point > pmap;

		for (auto& it : plist) {
			auto&[i1, t1] = it;
			auto&[d1, d2, d3] = t1;
			pmap.insert({ i1,Point(d1, d2, d3) });
		}
		for (auto& it : vlist) {

			tt = it;
			for (auto& it : tt) {
				auto[t1, f] = count.insert({ it,pn_all });
				if (f) {
					pn_all++;
					pvec.push_back(pmap[it]);
				}
				it = t1->second;
			}


			std::sort(tt.begin(), tt.end());
			cell.push_back({ tt[0], tt[1], tt[2],3  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[1], tt[3],2  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[2], tt[3],1  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[1], tt[2], tt[3],0  ,vn_all++,tt[0], tt[1], tt[2] ,tt[3] });

		}


	}

	boost::hash<std::tuple<int, int, int>> hash_tupple;
	auto hash_key = [&](int a)->std::size_t { return hash_tupple(std::make_tuple(cell[a][0], cell[a][1], cell[a][2])); };

	auto equal_key = [&](int a, int b) {
		return(std::make_tuple(cell[a][0], cell[a][1], cell[a][2]) ==
			std::make_tuple(cell[b][0], cell[b][1], cell[b][2]));
	};

	std::unordered_set<int, std::function<std::size_t(int)>,
		std::function<bool(int, int)>> uset_f(pn_all, hash_key, equal_key);


	std::vector<char> ori(vn_all, 0);
	{
		typedef std::pair< int, int > E;
		std::vector<E> edge_array;
		std::vector<E> weight;
		typedef adjacency_list< vecS, vecS, directedS, no_property, property<edge_color_t, std::pair<int, int>> > graph;
		typedef graph_traits<graph> GraphTraits;
		typename GraphTraits::edge_descriptor e;
		for (int i = 0; i < cell.size(); i++) {
			auto f = uset_f.insert(i);
			if (f.second)
				continue;
			cellinfo& c1 = cell[*(f.first)], &c2 = cell[i];
			edge_array.push_back({ c1[4],c2[4] });
			weight.push_back({ c1[3],c2[3] });
			edge_array.push_back({ c2[4],c1[4] });
			weight.push_back({ c2[3],c1[3] });

		}
		uset_f.clear();
		graph g(edge_array.begin(), edge_array.end(), weight.begin(), vn_all);
		auto name = get(edge_color, g);
		std::vector<char> cluster(vn_all, 0);
		bfs_time_visitor<decltype(name)> vis(cluster, name, ori);
		for (int i = 0; i < vn_all; i++) {
			if (cluster[i])
				continue;
			cluster[i] = 1;
			breadth_first_search(g, i, visitor(vis));
		}
		if (vis.f == 0)
			std::cerr << "vis.f == 0\n";
		for (auto& it : cluster) {
			if (it == 0)
				std::cerr << "it==0;\n";
		}

	}

	{
		std::vector<bool> is_vol_create(vn_all, 0);
		std::vector < std::array<Dart_handle, 4>> dartv(vn_all);

		std::unordered_map<std::tuple<int, int>, int, boost::hash<std::tuple<int, int>>> m1;
		std::unordered_map<std::tuple<int, int>, int, boost::hash<std::tuple<int, int>>>::iterator it;
		bool f;
		int en = 0, fn = 0, vn = 0;

		auto create_cell = [&, this](cellinfo& arr, int curvol) ->auto{
			if (ori[arr[4]]) {

				auto a1 = mesh.create_vertex_attribute(pvec[arr[5]]);
				a1->info() = res[0]++;
				auto a2 = mesh.create_vertex_attribute(pvec[arr[6]]);
				a2->info() = res[0]++;
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				a3->info() = res[0]++;
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				a4->info() = res[0]++;
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);
				is_vol_create[curvol] = 1;
				dartv[curvol] = { mesh.beta<1, 2,1>(d3), mesh.beta<0, 2,0>(d3), mesh.beta<2,1>(d3), mesh.beta<0>(d3) };



				return d3;

			}
			else {
				auto a1 = mesh.create_vertex_attribute(pvec[arr[6]]);
				a1->info() = res[0]++;
				auto a2 = mesh.create_vertex_attribute(pvec[arr[5]]);
				a2->info() = res[0]++;
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				a3->info() = res[0]++;
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				a4->info() = res[0]++;
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);

				is_vol_create[curvol] = 1;

				dartv[curvol] = { mesh.beta<0, 2,0>(d3), mesh.beta<1, 2,1>(d3), mesh.beta<2,0>(d3), mesh.beta<1>(d3) };

				return d3;
			}
		};
		auto modify_info = [&, this](Dart_handle adh) {

			for (auto it = mesh.darts_of_orbit<1, 2>(adh).begin(),
				a2 = mesh.darts_of_orbit<1, 2>(adh).end(); it != a2; it++) {
					{
						const int i = 1;
						if (mesh.attribute<i>(it) == NULL)
							mesh.set_attribute<i>(it, mesh.create_attribute<i>(res[i]++));
					}
					{
						const int i = 2;
						if (mesh.attribute<i>(it) == NULL)
							mesh.set_attribute<i>(it, mesh.create_attribute<i>(res[i]++));
					}
					{
						const int i = 3;
						if (mesh.attribute<i>(it) == NULL)
							mesh.set_attribute<i>(it, mesh.create_attribute<i>(res[i]++));
					}

			}

		};

		std::array<Dart_handle, 4> dl;
		std::array<int, 4> il;

		int curvol;
		Dart_handle d1, d2, d3;

		auto outfmt = format("%d/%d\n");


		for (int i = 0; i < cell.size(); i++) {
			//std::cout << outfmt % i%cell.size();
			auto aa = uset_f.insert(i);


			if (aa.second) {
				auto& arr = cell[i];
				curvol = arr[4];

				if (!(is_vol_create[curvol])) {
					auto d3 = create_cell(arr, curvol);
					modify_info(d3);


				}

				d1 = dartv[curvol][arr[3]];


			}
			else {
				auto& arr = cell[i];
				curvol = arr[4];
				if (!(is_vol_create[curvol])) {

					auto d3 = create_cell(arr, curvol);
					modify_info(d3);

				}

				d2 = dartv[curvol][arr[3]];

				auto id_d11 = *(aa.first);
				d1 = dartv[cell[id_d11][4]][cell[id_d11][3]];
				//std::cout << format("%d %d - %d %d\n") % curvol % arr[3] % cell[id_d11][4] % cell[id_d11][3];




				assert(mesh.is_sewable<3>(d1, d2));
				mesh.sew<3>(d1, d2);
			}














		}


		{
			std::array<int, 4> res = { 0,0,0,0 };


			{
				const int i = 0;
				for (auto it = mesh.attributes<i>().begin(), itend = mesh.attributes<i>().end(); it != itend; ++it)
				{
					it->info() = res[i]++;
				}
			}

			{
				const int i = 1;
				for (auto it = mesh.attributes<i>().begin(), itend = mesh.attributes<i>().end(); it != itend; ++it)
				{
					it->info() = res[i]++;
				}
			}
			{
				const int i = 2;
				for (auto it = mesh.attributes<i>().begin(), itend = mesh.attributes<i>().end(); it != itend; ++it)
				{
					it->info() = res[i]++;
				}
			}
			{
				const int i = 3;
				for (auto it = mesh.attributes<i>().begin(), itend = mesh.attributes<i>().end(); it != itend; ++it)
				{
					it->info() = res[i]++;
				}
			}
			VC0.resize(res[0]);
			VC1.resize(res[1]);
			VC2.resize(res[2]);
			VC3.resize(res[3]);
			std::cout << format("count:%d %d %d %d genus:%d\n")
				% res[0] % res[1] % res[2] % res[3] % (res[0] - res[1] + res[2] - res[3]);
			for (auto it(mesh.darts().begin()),
				itend(mesh.darts().end()); it != itend; ++it)
			{
				VC0[mesh.info<0>(it)].da = it;
				VC1[mesh.info<1>(it)].da = it;
				VC2[mesh.info<2>(it)].da = it;
				VC3[mesh.info<3>(it)].da = it;
			}

		}




	}
	disp_one_cell();
	mesh.display_characteristics(std::cout);
	//std::cout << format("\n%d %d %d %d\n")%countcell[0]% countcell[1] % countcell[2]% countcell[3];
	assert(mesh.is_valid());
	/*std::cout << format("\nis valid: %d\n") % mesh.is_valid();*/
	std::cout << format("\n");

}
class Sum_functor
{
public:
	Sum_functor(int& out) :count(out) {};

	template<class Cell_attribute>
	void operator()(Cell_attribute& ca1, Cell_attribute& ca2)
	{
		//std::cout << format("%1% %2% %3%\n") % (typeid(ca1.info()).name()) % (ca1.info()) % (ca2.info());
		auto& i1 = ca1.info(), &i2 = ca2.info();
		if (i1 != i2)
			count++;
		ca1.info() = std::min(ca1.info(), ca2.info());
		/*std::cout << "i am in\n";*/





	}

protected:
	int& count;
};
template<class mymesh>
template<typename pclass, typename vclass>
inline void readmesh<mymesh>::read_with_count(pclass plist, vclass vlist,
	std::vector<std::array<int, 4>>& trace_res)
{
	//number of vulume,number of point
	int vn_all = 0, pn_all = 0;
	std::array<int, 4> res = { 0,0,0,0 };
	std::array<int, 4> subres = res;
	trace_res.reserve(vlist.size());

	mesh.onmerge_functor<0>() = Sum_functor(subres[0]);
	mesh.onmerge_functor<1>() = Sum_functor(subres[1]);
	mesh.onmerge_functor<2>() = Sum_functor(subres[2]);
	mesh.onmerge_functor<3>() = Sum_functor(subres[3]);
	//mesh.set_automatic_attributes_management(true);
	using namespace boost;
	typedef std::array<int, 9> cellinfo;







	{
		int i1, i2, i3, i4, i5;
		double d1, d2, d3;
		std::unordered_map<int, int> count;
		std::array<int, 4> tt;
		std::unordered_map<int, Point > pmap;

		for (auto& it : plist) {
			auto&[i1, t1] = it;
			auto&[d1, d2, d3] = t1;
			pmap.insert({ i1,Point(d1, d2, d3) });
		}
		for (auto& it : vlist) {

			tt = it;
			for (auto& it : tt) {
				auto[t1, f] = count.insert({ it,pn_all });
				if (f) {
					pn_all++;
					pvec.push_back(pmap[it]);
				}
				it = t1->second;
			}


			std::sort(tt.begin(), tt.end());
			cell.push_back({ tt[0], tt[1], tt[2],3  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[1], tt[3],2  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[2], tt[3],1  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[1], tt[2], tt[3],0  ,vn_all++,tt[0], tt[1], tt[2] ,tt[3] });

		}


	}
	VC0 = std::vector<C0>(pn_all);
	boost::hash<std::tuple<int, int, int>> hash_tupple;
	auto hash_key = [&](int a)->std::size_t { return hash_tupple(std::make_tuple(cell[a][0], cell[a][1], cell[a][2])); };

	auto equal_key = [&](int a, int b) {
		return(std::make_tuple(cell[a][0], cell[a][1], cell[a][2]) ==
			std::make_tuple(cell[b][0], cell[b][1], cell[b][2]));
	};

	std::unordered_set<int, std::function<std::size_t(int)>,
		std::function<bool(int, int)>> uset_f(pn_all, hash_key, equal_key);


	std::vector<char> ori(vn_all, 0);
	{
		typedef std::pair< int, int > E;
		std::vector<E> edge_array;
		std::vector<E> weight;
		typedef adjacency_list< vecS, vecS, directedS, no_property, property<edge_color_t, std::pair<int, int>> > graph;
		typedef graph_traits<graph> GraphTraits;
		typename GraphTraits::edge_descriptor e;
		for (int i = 0; i < cell.size(); i++) {
			auto f = uset_f.insert(i);
			if (f.second)
				continue;
			cellinfo& c1 = cell[*(f.first)], &c2 = cell[i];
			edge_array.push_back({ c1[4],c2[4] });
			weight.push_back({ c1[3],c2[3] });
			edge_array.push_back({ c2[4],c1[4] });
			weight.push_back({ c2[3],c1[3] });

		}
		uset_f.clear();
		graph g(edge_array.begin(), edge_array.end(), weight.begin(), vn_all);
		auto name = get(edge_color, g);
		std::vector<char> cluster(vn_all, 0);
		bfs_time_visitor<decltype(name)> vis(cluster, name, ori);
		for (int i = 0; i < vn_all; i++) {
			if (cluster[i])
				continue;
			cluster[i] = 1;
			breadth_first_search(g, i, visitor(vis));
		}
		if (vis.f == 0)
			std::cerr << "vis.f == 0\n";
		for (auto& it : cluster) {
			if (it == 0)
				std::cerr << "it==0;\n";
		}

	}

	{
		std::vector<bool> is_vol_create(vn_all, 0);
		std::vector < std::array<Dart_handle, 4>> dartv(vn_all);

		std::unordered_map<std::tuple<int, int>, int, boost::hash<std::tuple<int, int>>> m1;
		std::unordered_map<std::tuple<int, int>, int, boost::hash<std::tuple<int, int>>>::iterator it;
		bool f;
		int en = 0, fn = 0, vn = 0;

		auto create_cell = [&, this](cellinfo& arr, int curvol) ->auto{
			if (ori[arr[4]]) {

				auto a1 = mesh.create_vertex_attribute(pvec[arr[5]]);
				a1->info() = res[0]++;
				auto a2 = mesh.create_vertex_attribute(pvec[arr[6]]);
				a2->info() = res[0]++;
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				a3->info() = res[0]++;
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				a4->info() = res[0]++;
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);
				is_vol_create[curvol] = 1;
				dartv[curvol] = { mesh.beta<1, 2,1>(d3), mesh.beta<0, 2,0>(d3), mesh.beta<2,1>(d3), mesh.beta<0>(d3) };



				return d3;

			}
			else {
				auto a1 = mesh.create_vertex_attribute(pvec[arr[6]]);
				a1->info() = res[0]++;
				auto a2 = mesh.create_vertex_attribute(pvec[arr[5]]);
				a2->info() = res[0]++;
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				a3->info() = res[0]++;
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				a4->info() = res[0]++;
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);

				is_vol_create[curvol] = 1;

				dartv[curvol] = { mesh.beta<0, 2,0>(d3), mesh.beta<1, 2,1>(d3), mesh.beta<2,0>(d3), mesh.beta<1>(d3) };

				return d3;
			}
		};
		auto modify_info = [&, this](Dart_handle adh) {

			for (auto it = mesh.darts_of_orbit<1, 2>(adh).begin(),
				a2 = mesh.darts_of_orbit<1, 2>(adh).end(); it != a2; it++) {
					{
						const int i = 1;
						if (mesh.attribute<i>(it) == NULL)
							mesh.set_attribute<i>(it, mesh.create_attribute<i>(res[i]++));
					}
					{
						const int i = 2;
						if (mesh.attribute<i>(it) == NULL)
							mesh.set_attribute<i>(it, mesh.create_attribute<i>(res[i]++));
					}
					{
						const int i = 3;
						if (mesh.attribute<i>(it) == NULL)
							mesh.set_attribute<i>(it, mesh.create_attribute<i>(res[i]++));
					}

			}

		};

		std::array<Dart_handle, 4> dl;
		std::array<int, 4> il;

		int curvol;
		Dart_handle d1, d2, d3;

		auto outfmt = format("%d/%d\n");

		int count = 0;
		for (int i = 0; i < cell.size(); i++) {
			//std::cout << outfmt % i%cell.size();
			auto aa = uset_f.insert(i);


			if (aa.second) {
				auto& arr = cell[i];
				curvol = arr[4];

				if (!(is_vol_create[curvol])) {
					auto d3 = create_cell(arr, curvol);
					modify_info(d3);


				}

				d1 = dartv[curvol][arr[3]];


			}
			else {
				auto& arr = cell[i];
				curvol = arr[4];
				if (!(is_vol_create[curvol])) {

					auto d3 = create_cell(arr, curvol);
					modify_info(d3);

				}

				d2 = dartv[curvol][arr[3]];

				auto id_d11 = *(aa.first);
				d1 = dartv[cell[id_d11][4]][cell[id_d11][3]];
				//std::cout << format("%d %d - %d %d\n") % curvol % arr[3] % cell[id_d11][4] % cell[id_d11][3];




				assert(mesh.is_sewable<3>(d1, d2));
				mesh.sew<3>(d1, d2);
			}


			count++;
			if (count == 4) {

				trace_res.push_back({ res[0] - subres[0],
					res[1] - subres[1],res[2] - subres[2],res[3] - subres[3] });
				count = 0;
			}









		}


	}

	assert(trace_res.size() == vlist.size());

	//std::cout << format("\n%d %d %d %d\n")%countcell[0]% countcell[1] % countcell[2]% countcell[3];
}

template<class mymesh>
template<typename pclass, typename vclass>
void readmesh<mymesh>::count_genus(pclass plist, vclass vlist)
{
	int vn_all = 0, pn_all = 0;

	using namespace boost;
	typedef std::array<int, 9> cellinfo;







	{
		int i1, i2, i3, i4, i5;
		double d1, d2, d3;
		std::unordered_map<int, int> count;
		std::array<int, 4> tt;
		std::unordered_map<int, Point > pmap;

		for (auto& it : plist) {
			auto&[i1, t1] = it;
			auto&[d1, d2, d3] = t1;
			pmap.insert({ i1,Point(d1, d2, d3) });
		}
		for (auto& it : vlist) {

			tt = it;
			for (auto& it : tt) {
				auto[t1, f] = count.insert({ it,pn_all });
				if (f) {
					pn_all++;
					pvec.push_back(pmap[it]);
				}
				it = t1->second;
			}


			std::sort(tt.begin(), tt.end());
			cell.push_back({ tt[0], tt[1], tt[2],3  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[1], tt[3],2  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[0], tt[2], tt[3],1  ,vn_all,tt[0], tt[1], tt[2] ,tt[3] });
			cell.push_back({ tt[1], tt[2], tt[3],0  ,vn_all++,tt[0], tt[1], tt[2] ,tt[3] });

		}


	}
	VC0 = std::vector<C0>(pn_all);
	boost::hash<std::tuple<int, int, int>> hash_tupple;
	auto hash_key = [&](int a)->std::size_t { return hash_tupple(std::make_tuple(cell[a][0], cell[a][1], cell[a][2])); };

	auto equal_key = [&](int a, int b) {
		return(std::make_tuple(cell[a][0], cell[a][1], cell[a][2]) ==
			std::make_tuple(cell[b][0], cell[b][1], cell[b][2]));
	};

	std::unordered_set<int, std::function<std::size_t(int)>,
		std::function<bool(int, int)>> uset_f(pn_all, hash_key, equal_key);


	std::vector<char> ori(vn_all, 0);
	{
		typedef std::pair< int, int > E;
		std::vector<E> edge_array;
		std::vector<E> weight;
		typedef adjacency_list< vecS, vecS, directedS, no_property, property<edge_color_t, std::pair<int, int>> > graph;
		typedef graph_traits<graph> GraphTraits;
		typename GraphTraits::edge_descriptor e;
		for (int i = 0; i < cell.size(); i++) {
			auto f = uset_f.insert(i);
			if (f.second)
				continue;
			cellinfo& c1 = cell[*(f.first)], &c2 = cell[i];
			edge_array.push_back({ c1[4],c2[4] });
			weight.push_back({ c1[3],c2[3] });
			edge_array.push_back({ c2[4],c1[4] });
			weight.push_back({ c2[3],c1[3] });

		}
		uset_f.clear();
		graph g(edge_array.begin(), edge_array.end(), weight.begin(), vn_all);
		auto name = get(edge_color, g);
		std::vector<char> cluster(vn_all, 0);
		bfs_time_visitor<decltype(name)> vis(cluster, name, ori);
		for (int i = 0; i < vn_all; i++) {
			if (cluster[i])
				continue;
			cluster[i] = 1;
			breadth_first_search(g, i, visitor(vis));
		}
		if (vis.f == 0)
			std::cerr << "vis.f == 0\n";
		for (auto& it : cluster) {
			if (it == 0)
				std::cerr << "it==0;\n";
		}

	}

	{
		std::vector<bool> is_vol_create(vn_all, 0);
		std::vector < std::array<Dart_handle, 4>> dartv(vn_all);

		std::unordered_map<std::tuple<int, int>, int, boost::hash<std::tuple<int, int>>> m1;
		std::unordered_map<std::tuple<int, int>, int, boost::hash<std::tuple<int, int>>>::iterator it;
		bool f;
		int en = 0, fn = 0, vn = 0;

		auto create_cell = [&, this](cellinfo& arr, int curvol) ->auto{
			if (ori[arr[4]]) {

				auto a1 = mesh.create_vertex_attribute(pvec[arr[5]]);
				auto a2 = mesh.create_vertex_attribute(pvec[arr[6]]);
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);
				is_vol_create[curvol] = 1;
				dartv[curvol] = { mesh.beta<1, 2,1>(d3), mesh.beta<0, 2,0>(d3), mesh.beta<2,1>(d3), mesh.beta<0>(d3) };



				return d3;

			}
			else {
				auto a1 = mesh.create_vertex_attribute(pvec[arr[6]]);
				auto a2 = mesh.create_vertex_attribute(pvec[arr[5]]);
				auto a3 = mesh.create_vertex_attribute(pvec[arr[7]]);
				auto a4 = mesh.create_vertex_attribute(pvec[arr[8]]);
				auto d3 = mesh.make_tetrahedron(a1, a2, a3, a4);

				is_vol_create[curvol] = 1;

				dartv[curvol] = { mesh.beta<0, 2,0>(d3), mesh.beta<1, 2,1>(d3), mesh.beta<2,0>(d3), mesh.beta<1>(d3) };

				return d3;
			}
		};


		std::array<Dart_handle, 4> dl;
		std::array<int, 4> il;

		int curvol;
		Dart_handle d1, d2, d3;
		mesh.set_automatic_attributes_management(false);
		std::array<size_type, 4> amark;
		for (int i = 0; i < 4; i++)
			amark[i] = mesh.get_new_mark();

		for (int i = 0; i < cell.size(); i++) {

			auto aa = uset_f.insert(i);


			if (aa.second) {
				auto& arr = cell[i];
				curvol = arr[4];

				if (!(is_vol_create[curvol])) {
					auto d3 = create_cell(arr, curvol);


				}

				d1 = dartv[curvol][arr[3]];

				continue;
			}






			auto& arr = cell[i];
			curvol = arr[4];
			if (!(is_vol_create[curvol])) {

				auto d3 = create_cell(arr, curvol);

		}

			d2 = dartv[curvol][arr[3]];

			auto id_d11 = *(aa.first);
			d1 = dartv[cell[id_d11][4]][cell[id_d11][3]];
			//std::cout << format("%d %d - %d %d\n") % curvol % arr[3] % cell[id_d11][4] % cell[id_d11][3];

#ifdef _DEBUG
			disp_loop(d1, d2);
#endif // _DEBUG


			assert(mesh.is_sewable<3>(d1, d2));
			mesh.sew<3>(d1, d2);






	}


}





}



template<class LCC_3>
template<int N>
inline void readmesh<LCC_3>::display_map_and_attributes()
{
	/*auto ss = boost::format("index: %1% (%2%)\n");
	for (auto it = mesh.attributes<N>().begin(), itend = mesh.attributes<N>().end();
		it != itend; ++it)
	{
		std::cout << ss % mesh.info_of_attribute<N>(it) % mesh.point_of_vertex_attribute(it);
	}
	std::cout << boost::format("is valid: %d\n\n") % mesh.is_valid();*/


	/*std::cout << "--display_map_and_attributes--s-\n";
	for (auto a1 = mesh.one_dart_per_cell<0>().begin(), a2 = mesh.one_dart_per_cell<0>().end(); a1 != a2; a1++) {
		auto ff1 = boost::format("<%1%>(%2%) ");
		for (auto dc1 = mesh.darts_of_cell<0>(a1).begin(), dc2 = mesh.darts_of_cell<0>(a1).end();
			dc1 != dc2; dc1++) {
			std::cout << ff1 % mesh.info<0>(dc1) % mesh.point(dc1);
		}
		std::cout << "\n";
	}
	std::cout << "--display_map_and_attributes--e-\n";*/

	std::cout << "--display_map_and_attributes--s-\n";
	for (auto a1 = mesh.one_dart_per_cell<0>().begin(), a2 = mesh.one_dart_per_cell<0>().end(); a1 != a2; a1++) {
		auto ff1 = boost::format("<%1%>");
		for (auto dc1 = mesh.darts_of_cell<0>(a1).begin(), dc2 = mesh.darts_of_cell<0>(a1).end();
			dc1 != dc2; dc1++) {
			std::cout << ff1 % mesh.info<0>(dc1);
		}
		std::cout << "\n";
	}
	std::cout << "--display_map_and_attributes--e-\n";
}
