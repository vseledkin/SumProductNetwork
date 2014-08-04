#ifndef POONDATASET_H
#define POONDATASET_H

#include <vector>
#include <string>
#include <set>
#include "PoonInstance.h"
#include "PoonParameter.h"


//Future note, this class should be abstract with a virtual load class
class PoonDataset{
public:

	std::vector<std::vector<int> > tmp_;  // buffer for proc imgs 
	std::shared_ptr<PoonParameter> params;

	// data
	std::string expDir_ = "../..";
	std::string olivettiRawFileName_ = expDir_ + "/data/olivetti/olivetti.raw";
	std::string dataDir_ = expDir_ + "/data/olivetti";
	std::string rstDir_ = expDir_ + "/results/olivetti/completions";
	std::string mdlDir_ = expDir_ + "/results/olivetti/models";
	int RESCALE_LEN_ = 100;

	
private:
	//note that these should be wrapped into their own class so I can overload get instance
	std::vector<PoonInstance> train_;
	std::vector<PoonInstance> test_;

public:
	PoonDataset(std::shared_ptr<PoonParameter> p){
		this->params = p;
		expDir_ = params->inputPath_;

		olivettiRawFileName_ = expDir_ + "/data/olivetti/olivetti.raw";
		dataDir_ = expDir_ + "/data/olivetti";
		rstDir_ = expDir_ + "/results/olivetti/completions";
		mdlDir_ = expDir_ + "/results/olivetti/models";
	};

	
	std::vector<PoonInstance>& getTrain(){ return train_; };
	std::vector<PoonInstance>& getTest(){ return test_; };

	// dataset
	void loadOlivetti();

private:
	
	void setInstance(std::vector<std::vector<double> >& buf, PoonInstance& inst);

	//silly function, should have randomized options and be more robust
	void genTestIdx(int maxSize, int testSize, std::set<int>& tis) {
		for (int i = maxSize - testSize; i<maxSize; i++) {
			tis.insert(i); if (tis.size() == testSize) break;
		}
	};

	void readOlivettiInstance(std::vector<std::vector<double> >& faces, int pi, PoonInstance& inst);
};




#endif