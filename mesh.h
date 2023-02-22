#pragma once
#include <CGAL/Generalized_map.h>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>

template<class geometry>
class mymesh{
	
public:

	struct Myitem
	{
		template<class Refs>
		struct Dart_wrapper
		{
			
			typedef CGAL::Cell_attribute_with_point< Refs, int, CGAL::Tag_true>
				V_attribute;
			typedef CGAL::Cell_attribute<Refs, int, CGAL::Tag_true> E_attrib;
			typedef std::tuple<V_attribute, E_attrib, E_attrib, E_attrib> Attributes;
		};
	};
	typedef CGAL::Linear_cell_complex_traits
		<3, CGAL::Exact_predicates_inexact_constructions_kernel> Traits;
	typedef CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, Traits, Myitem> LCC_3;
	typedef typename LCC_3::Dart_handle                                 Dart_handle;

	using C0 = typename geometry::template C0<Dart_handle>;
	using C1 = typename geometry::template C1<Dart_handle>;
	using C2 = typename geometry::template C2<Dart_handle>;
	using C3 = typename geometry::template C3<Dart_handle>;
	


	LCC_3& mesh() { return _mesh; }
	std::vector<C0> VC0;
	std::vector<C1> VC1;
	std::vector<C2> VC2;
	std::vector<C3> VC3;
protected:
	LCC_3 _mesh;

};