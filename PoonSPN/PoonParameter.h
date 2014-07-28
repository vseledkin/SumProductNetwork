/*
Like Poon's "Parameter" class

Stores variables and process arguments

This should be a singleton.  The reasoning is that it will be inconviente to make fast changes if I have to pass parmaters or a param object around all the time.
I should pass it into the main class though and try to pass it through as much as possible.
*/

#ifndef POONPARAMETER_H
#define POONPARAMETER_H

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace std;


class PoonParameter {
public:
	int maxIter_;
	double thresholdLLHChg_;
	int batch_size_ = 50;
	double sparsePrior_;

	// SPN
	int numSumPerRegion_;
	int inputDim1_;
	int inputDim2_;
	int baseResolution_;
	double smoothSumCnt_;
	int numComponentsPerVar_;

	// Eval
	int maxTestSize_;
	string domain_;
	int numSlavePerClass_;
	int numSlaveGrp_;

	//Set default values
	PoonParameter(){
		int maxIter_ = 30;
		double thresholdLLHChg_ = 0.1;
		int batch_size_ = 50;
		double sparsePrior_ = 1;

		// SPN
		int numSumPerRegion_ = 20;
		int inputDim1_ = 64;
		int inputDim2_ = 64;
		int baseResolution_ = 4;
		double smoothSumCnt_ = 0.01;
		int numComponentsPerVar_ = 4;

		// Eval
		int maxTestSize_ = 50;
		string domain_ = "";
		int numSlavePerClass_ = 50;
		int numSlaveGrp_ = 1;
	};

	//Poon's proc and procArgs

	//return true if valid, false if not
	bool processArgs(int argc, char* argv[]){
		//parse the arguments
		
		vector<string> args(argv + 1, argv + argc);
		for (unsigned int i = 0; i < args.size(); i++){
			if (args[i] == "-d") {
				this->domain_ = args[++i];
			}
			else if (args[i] == "-ncv") {
				this->numComponentsPerVar_ = atoi(args[++i].c_str());
			}
			else if (args[i] == "-nsr") {
				this->numSumPerRegion_ = atoi(args[++i].c_str());
			}
			else if (args[i] == "-sp") {
				this->sparsePrior_ = atof(args[++i].c_str());
			}
			else if (args[i] == "-br") {
				this->baseResolution_ = atoi(args[++i].c_str());
			}
			else if (args[i] == "-ct") {
				this->thresholdLLHChg_ = atof(args[++i].c_str());
			}
			else if (args[i] == "-bs") {
				this->batch_size_ = atoi(args[++i].c_str());
			}
			else if (args[i] == "-ns") {
				this->numSlavePerClass_ = atoi(args[++i].c_str());
			}
			else if (args[i] == "-nsg") {
				this->numSlaveGrp_ = atoi(args[++i].c_str());
			}
		}

		//safty check
		if (this->domain_ == "") {
			cout << "\n\nOptions: [-d <domain>]\n[-sp <sparsePrior>]\n[-br <baseResolution>]\n[-ncv <numComponentsPerVar>]\n[-nsr <numSumPerRegion>]\n[-ct <convergencyThrehold>]\n[-bs <batchSize>]\n[-ns <numSlavePerCat>]\n[-nsg <numSlaveGrp>]\n";
			return false;
		}
		
		return true;
	};

};

#endif