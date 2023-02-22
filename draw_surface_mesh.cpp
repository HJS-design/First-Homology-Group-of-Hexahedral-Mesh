
#include "mesh.h"
#include "geometry.h"
#include "readmesh.cpp"
#include "drawdart.h"

#include "loop.h"
#include "judgegenus.h"
#include<boost/algorithm/string.hpp>

int main(int argc, char *argv[]) {




	typedef typename mymesh<geometry> tmesh;
	tmesh mmesh;
	readmesh hh(mmesh);
	hh.read_hex_file_non_manifold_info(boost::lexical_cast<std::string>(argv[1]));
	loop hh1(mmesh,boost::lexical_cast<int>(argv[2]), boost::lexical_cast<std::string>(argv[3]));
	hh1.hex_unweighted_loop();
	drawdart mydraw(mmesh);
	mydraw.draw();




}

