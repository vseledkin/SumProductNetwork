/*  Arvid Frydenlund 2014-05-06

A reconstruction of Poon's "Run" class

This class acts as the main.
*/

#ifndef POONRUN_H
#define POONRUN_H

#include <string>
#include <vector>

#include "PoonParameter.h"



class PoonRun{

public:
	// domains
	const std::string DOM_OLIVETTI_ = "O";
	const std::string DOM_CALTECH_ = "C";

	// directory
	const std::string expDir_ = "/projects/dm/2/hoifung/projects/dspn/release";

	const std::string oliveDataDir_ = expDir_ + "/data/olivetti";
	const std::string oliveRstDir_ = expDir_ + "/results/olivetti/completions";
	const std::string oliveMdlDir_ = expDir_ + "/results/olivetti/models";

	const std::string calDataDir_ = expDir_ + "/data/caltech";
	const std::string calRstDir_ = expDir_ + "/results/caltech/completions";
	const std::string calMdlDir_ = expDir_ + "/results/caltech/models";

	PoonRun(int argc, char* argv[]); //main

private:
	//For running the different data sets
	void runOlivetti(std::shared_ptr<PoonParameter> Params);
	void runCaltech(std::shared_ptr<PoonParameter> Params);

};




#endif