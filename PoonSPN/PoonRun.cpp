
/*
package eval;

import spn.SPN;

import spn.GenerativeLearning;

import java.io.*;
import java.util.*;

import common.*;

import mpi.MPI;
*/

#include "PoonRun.h"
#include "PoonParameter.h"
#include "PoonDataset.h"

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
	//TODO LATER
}


void PoonRun::runOlivetti(std::shared_ptr<PoonParameter> Params){
	
	bool debug = true;

	//reading in data set
	PoonDataset data(Params);
	data.loadOlivetti();

	if (debug){
		//write out images to confirm they work
	}
	

	/*
	// learn
	GenerativeLearning l = new GenerativeLearning();
	l.learn(data.getTrain());
	SPN dspn = l.getDSPN();
	if (MyMPI.myOffset_ == 0) {
		dspn.saveDSPN(oliveMdlDir_ + "/olive");
	}

	// complete
	ImageCompletion.completeLeft(dspn, data.getTest(), "olive", oliveRstDir_);
	ImageCompletion.completeBottom(dspn, data.getTest(), "olive", oliveRstDir_);
	*/
	cout << "Done, runOlivetti!" << endl;
}


