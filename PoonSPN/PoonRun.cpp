
#include "PoonRun.h"
#include "PoonParameter.h"
#include "PoonDataset.h"
#include "PoonGenerativeLearning.h"
#include "Utils/ImagePGM.h"

#include <sstream> 

using namespace std;

PoonRun::PoonRun(int argc, char* argv[]){

	SharedParams<PoonParameter> spp;
	std::shared_ptr<PoonParameter> Params = spp.instance();

	if(!Params->processArgs(argc, argv)){
		return;  //failed to process args
	}

	if (Params->domain_ == this->DOM_OLIVETTI_) {
		runOlivetti(Params);
	}
	else if (Params->domain_ == this->DOM_CALTECH_) {
		runCaltech(Params);
	}
	else{
		cout << "Did not run, check the value of the domain" << endl;
	}
}


void PoonRun::runCaltech(std::shared_ptr<PoonParameter>Params){
	//TODO LATER ..or never which ever
}


void PoonRun::runOlivetti(std::shared_ptr<PoonParameter> params){
	
	bool debug = false;

	//reading in data set
	PoonDataset data(params);
	data.loadOlivetti();

	if (debug){
		//write out images to confirm they work
		int i = 0;
		for (auto& inst : data.getTrain()){
			stringstream s;
			s << data.rstDir_ << "/" << i << "_test_read_in.pgm";
			ImagePGM image(inst.vals_);
			image.rescaleImageValues();
			image.writeImage(s.str());
			i++;
		}
	}
	

	
	// learn
	PoonGenerativeLearning l;
	l.learn(data.getTrain(), params);
	PoonSPN dspn = l.getDSPN();


	// complete
	//ImageCompletion.completeLeft(dspn, data.getTest(), "olive", oliveRstDir_);
	//ImageCompletion.completeBottom(dspn, data.getTest(), "olive", oliveRstDir_);
	
	cout << "Done, runOlivetti!" << endl;
}


