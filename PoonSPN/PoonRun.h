/*  Arvid Frydenlund 2014-05-06

A reconstruction of Poon's "Run" class

This class acts as the main.
*/

#ifndef POONRUN_H
#define POONRUN_H

#include <string>
#include <vector>

#include "PoonParameter.h"

using namespace std;


class PoonRun{

public:
	// domains
	const string DOM_OLIVETTI_ = "O";
	const string DOM_CALTECH_ = "C";

	// directory
	const string expDir_ = "/projects/dm/2/hoifung/projects/dspn/release";

	const string oliveDataDir_ = expDir_ + "/data/olivetti";
	const string oliveRstDir_ = expDir_ + "/results/olivetti/completions";
	const string oliveMdlDir_ = expDir_ + "/results/olivetti/models";

	const string calDataDir_ = expDir_ + "/data/caltech";
	const string calRstDir_ = expDir_ + "/results/caltech/completions";
	const string calMdlDir_ = expDir_ + "/results/caltech/models";

	/* don't know, think it finds all files with ".raw.rescale"
	static FilenameFilter caltechImgNameFilter_ = new FilenameFilter(){
		@Override
		public boolean accept(File dir, String name) {
			if (name.indexOf(".raw.rescale")>0) return true;
			return false;
		}
	};

	*/

	PoonRun(int argc, char* argv[]); //main

private:
	//For running the different data sets
	void runOlivetti(PoonParameter& Params);
	void runCaltech(PoonParameter& Params);

	

};




#endif