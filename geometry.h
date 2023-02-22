#pragma once

struct geometry {
	template<class dart>
	struct C0 {
		dart da;
		bool sharp = false;
	};
	template<class dart>
	struct C1 {
		dart da;
		int id = 0;
		bool generator = false;
		int pair = -1;
		bool sharp = false;

	};
	template<class dart>
	struct C2 {
		dart da; int id = 0; bool generator = false;
	};
	template<class dart>
	struct C3 { dart da; };
};