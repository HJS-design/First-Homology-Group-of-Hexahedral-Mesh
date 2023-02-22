#pragma once

#include <igl/opengl/glfw/Viewer.h>
template<class mymesh>
class drawdart {



	typedef typename mymesh::LCC_3					LCC_3;
	typedef typename LCC_3::Dart_handle                                 Dart_handle;
	typedef typename LCC_3::Point                                       Point;
	using C0 = typename mymesh::C0;
	using C1 = typename mymesh::C1;
	using C2 = typename mymesh::C2;
	using C3 = typename mymesh::C3;
public:
	/*readmesh() {};*/
	drawdart(mymesh& out) :mesh(out.mesh()),
		VC0(out.VC0), VC1(out.VC1), VC2(out.VC2), VC3(out.VC3)
	{

	};
	void draw();

protected:
	LCC_3& mesh;
	std::vector<C0>& VC0;
	std::vector<C1>& VC1;
	std::vector<C2>& VC2;
	std::vector<C3>& VC3;
};

template<typename LCC_3>
inline void drawdart<LCC_3>::draw()
{
	std::cout << "---draw---\n";
	std::vector<unsigned int> cells = { 0,3 };

	std::array<int, 4> res = { VC0.size(),VC1.size() ,VC2.size() ,VC3.size() };
	using namespace Eigen;
	Matrix<double, Dynamic, Dynamic, RowMajor> V;
	Matrix<int, Dynamic, Dynamic, RowMajor> F;

	Matrix<int, Dynamic, Dynamic, RowMajor> E;
	Matrix<double, Dynamic, Dynamic, RowMajor> myPoint;
	
	{
		std::vector<double> pV;
		pV.reserve(res[0] * 3);
		std::vector<int> pF;
		pF.reserve(res[3] * 4);
		std::vector<std::array<double, 3>> tpV(res[0]);
		for (int i = 0; i < VC0.size(); i++) {
			auto& p = mesh.point(VC0[i].da);
			tpV[i] = { p.x(),p.y(),p.z() };
		}
		
		for (auto& it : tpV) {
			pV.push_back(it[0]);
			pV.push_back(it[1]);
			pV.push_back(it[2]);
		}

		int count = 0, cc1 = 0,i1,i2;
		std::list<int> tpF;
		for (int i = 0; i < VC2.size(); i++) {
			auto a1 = VC2[i].da;
			if (!mesh.is_free<3>(a1))
				continue;

			/*count = 0;*/
			tpF.clear();
			for (auto dc1 = mesh.one_dart_per_incident_cell<0, 2>(a1).begin(),
				dc2 = mesh.one_dart_per_incident_cell<0, 2>(a1).end();
				dc1 != dc2; dc1++) {
				/*	count++;*/
				tpF.push_back(mesh.info<0>(dc1));
				
			}
			auto it = tpF.begin();
			i1 = *(it++);
			i2 = *(it++);
			while (it != tpF.end()) {
				pF.push_back(i1);
				pF.push_back(i2);
				i2 = *(it++);
				pF.push_back(i2);
			}
			
		}
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
		pV.clear();
		for (auto& it : VC0) {
			if (it.sharp) {
				
				auto& p = mesh.point(it.da);
				pV.push_back(p.x());
				pV.push_back(p.y());
				pV.push_back(p.z());
			
			}

		}
		myPoint = Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(pV.data(), pV.size() / 3, 3);
	}





	//std::cout << boost::format("%d %d %d %d\n") % V.rows() % V.cols() % F.rows() % F.cols();
	/*std::cout << V;
	std::cout << "\n";
	std::cout << F;*/
	auto C =
		(V.rowwise() - V.colwise().minCoeff()).array().rowwise() /
		(V.colwise().maxCoeff() - V.colwise().minCoeff()).array();
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	//std::ofstream vv("E:\\hjs\\matlabcode\\hjs1\\file1\\V.data");
	//vv << V;
	//std::ofstream ff("E:\\hjs\\matlabcode\\hjs1\\file1\\F.data");
	//ff << F;
	viewer.data().set_colors(C);
	if(E.size()>0)
		viewer.data().set_edges(V,E, Matrix<double, 1, 3>(1,0,0));
	if (myPoint.size() > 0)
		viewer.data().set_points(myPoint, Matrix<double, 1, 3>(0, 1, 0));
	viewer.data().line_width=3;
	//viewer.core().background_color= Vector4f(1,1,1,1);
	// Launch the viewer
	viewer.launch();
	
}
