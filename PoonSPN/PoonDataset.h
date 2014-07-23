#ifndef POONDATASET_H
#define POONDATASET_H

#include <vector>
#include <string>
#include <set>
#include "PoonInstance.h"
#include "PoonParameter.h"

using namespace std;

//Future note, this class should be abstract with a virtual load class
class PoonDataset{
public:

	vector<vector<int> > tmp_;  // buffer for proc imgs 
	PoonParameter params;

	// data
	string expDir_ = "../..";
	string olivettiRawFileName_ = expDir_ + "/data/olivetti/olivetti.raw";
	string calDataDir_ = expDir_ + "/data/caltech";
	string calRstDir_ = expDir_ + "/results/caltech/completions";
	string calMdlDir_ = expDir_ + "/results/caltech/models";
	int RESCALE_LEN_ = 100;
	/*
	static FilenameFilter calFileNameFilter_ = new FilenameFilter(){
		public boolean accept(File dir, String name) {
			if (name.indexOf(".raw.rescale")<0) return false; else return true;
		}
	};*/
	
private:
	vector<PoonInstance> train_;
	vector<PoonInstance> test_;

public:
	PoonDataset(PoonParameter& p){
		this->params = p;
	};

	vector<PoonInstance>& getTrain(){ return train_; };
	vector<PoonInstance>& getTest(){ return test_; };

	// dataset
	void loadOlivetti();

private:
	
	void setInstance(vector<vector<double> >& buf, PoonInstance inst);

	//silly function, should have randomized options and be more robust
	void genTestIdx(int maxSize, int testSize, set<int>& tis) {
		for (int i = maxSize - testSize; i<maxSize; i++) {
			tis.insert(i); if (tis.size() == testSize) break;
		}
	};

	void readOlivettiInstance(vector<vector<double> >& faces, int pi, PoonInstance& inst);
};




#endif