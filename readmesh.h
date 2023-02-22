

#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <cstdlib>
using namespace Eigen;
template<class mymesh>
class readmesh {


	typedef typename mymesh::LCC_3					LCC_3;
	typedef typename LCC_3::Dart_handle                                 Dart_handle;
	typedef typename LCC_3::Point                                       Point;
	using C0 = typename mymesh::C0;
	using C1 = typename mymesh::C1;
	using C2 = typename mymesh::C2;
	using C3 = typename mymesh::C3;
public:
	/*readmesh() {};*/
	readmesh(mymesh& out) :mesh(out.mesh()),
		VC0(out.VC0), VC1(out.VC1), VC2(out.VC2), VC3(out.VC3)
	{

	};

	template<int N>
	void display_map_and_attributes();

	void disp(Dart_handle d);
	void disp_loop(Dart_handle d, Dart_handle d2);
	void disp_one_cell();
	void read_t_file(std::string filename);
	template <typename pclass, typename vclass>
	void read_graph_file(std::string filename,pclass &plist, vclass &vlist);
	template <typename pclass, typename vclass>
	void read_graph_file_with_space(std::string filename, pclass &plist, vclass &vlist);

	void read_hex_file_non_manifold_info(std::string filename);
	template <typename pclass, typename vclass>
	void read_manifold_info(pclass plist, vclass vlist);
	template <typename pclass, typename vclass>
	void read_without_info(pclass plist, vclass vlist);
	template <typename pclass, typename vclass>
	void read_non_manifold_info(pclass plist, vclass vlist);
	//my count is not true for un_manifold mesh in read() 
	template <typename pclass, typename vclass>
	void read_with_count(pclass plist, vclass vlist, std::vector<std::array<int, 4>>& trace_res);
	
	template <typename pclass, typename vclass>
	void count_genus(pclass plist, vclass vlist);
	void get_view_matrix(Matrix<double, Dynamic, Dynamic, RowMajor>& V,
		Matrix<int, Dynamic, Dynamic, RowMajor>& F,

		Matrix<int, Dynamic, Dynamic, RowMajor>& E);
protected:
	LCC_3& mesh;

	std::vector<C0>& VC0;
	std::vector<C1>& VC1;
	std::vector<C2>& VC2;
	std::vector<C3>& VC3;
	std::vector<std::array<int, 9>>cell;
	std::vector<Point> pvec;
};

