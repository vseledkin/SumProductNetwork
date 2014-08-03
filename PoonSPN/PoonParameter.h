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
#include <memory>
#include <iostream>
#include <stdlib.h>
#include <random>


class PoonParameter {
public:
	//random



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
	std::string domain_;
	int numSlavePerClass_;
	int numSlaveGrp_;


	//MPI values
	// buffer
	static int buf_idx_;
	static const int buf_size_ = 10000000;
	static const int buf_size_double_ = 100;
	static std::vector<int> buf_int_;  //instantiation in cpp
	static std::vector<double> buf_double_;  //instantiation in cpp
	static std::string buf_char_;

	// MPI util
	/*
	static double recvDouble(int src, int tag) {
		MPI.COMM_WORLD.Recv(MyMPI.buf_double_, 0, 1, MPI.DOUBLE, src, tag);
		return MyMPI.buf_double_[0];
	}
	static void sendDouble(int dest, int tag, double d) {
		MyMPI.buf_double_[0] = d;
		MPI.COMM_WORLD.Send(MyMPI.buf_double_, 0, 1, MPI.DOUBLE, dest, tag);
	}
	static char recvChar(int src, int tag) {
		MPI.COMM_WORLD.Recv(MyMPI.buf_char_, 0, 1, MPI.CHAR, src, tag);
		return MyMPI.buf_char_[0];
	}
	static void sendChar(int dest, int tag, char c) {
		MyMPI.buf_char_[0] = c;
		MPI.COMM_WORLD.Send(MyMPI.buf_char_, 0, 1, MPI.CHAR, dest, tag);
	}*/

	//Set default values
	PoonParameter(){
		maxIter_ = 30;
		thresholdLLHChg_ = 0.1;
		batch_size_ = 50;
		sparsePrior_ = 1;

		// SPN
		numSumPerRegion_ = 20;
		inputDim1_ = 64;
		inputDim2_ = 64;
		baseResolution_ = 4;
		smoothSumCnt_ = 0.01;
		numComponentsPerVar_ = 4;

		// Eval
		maxTestSize_ = 50;
		domain_ = "";
		numSlavePerClass_ = 50;
		numSlaveGrp_ = 1;

		//MPI
		buf_idx_ = 0;
	};

	//Poon's proc and procArgs

	//return true if valid, false if not
	bool processArgs(int argc, char* argv[]){
		//parse the arguments
		
		std::vector<std::string> args(argv + 1, argv + argc);
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
			std::cout << "\n\nOptions: [-d <domain>]\n[-sp <sparsePrior>]\n[-br <baseResolution>]\n[-ncv <numComponentsPerVar>]\n[-nsr <numSumPerRegion>]\n[-ct <convergencyThrehold>]\n[-bs <batchSize>]\n[-ns <numSlavePerCat>]\n[-nsg <numSlaveGrp>]\n";
			return false;
		}
		
		return true;
	};

};

//a class to turn the paramaters into a singletonesque class
//the instance will stay alive as long as there is an instance of it alive somwhere.
//This alows one to use it as singleton if they desire.
template <typename PoonParameter>
class SharedParams {
public:
	static std::shared_ptr<PoonParameter> instance()
	{
		auto ptr = s_instance.lock();

		if (!ptr) {
			ptr = std::make_shared<PoonParameter>();
			s_instance = ptr;
		}

		return ptr;
	}

private:
	static std::weak_ptr<Target> s_instance;
};

template <typename PoonParameter> std::weak_ptr<PoonParameter> SharedPrams<PoonParameter>::s_instance;

#endif