#include <vector>
#include <map>
#include <string>
#include <cmath>

#include "PoonProdNode.h"
#include "PoonSumNode.h"
#include "PoonParameter.h"

class PoonRegion
{
public:
	int id_;
	int a1_, a2_, b1_, b2_;
	int a_, b_;	// a=a2-a1, b=b2-b1	
	int interval_;	// for coarse resolution	

	// for pixel region only: gaussian units
	std::vector<double> means_;
	std::vector<double> vars_;
	std::vector<double> cnts_;
	double ttlCnt_ = 0;

	// data structure for a parse
	std::map<int, int> inst_type_;
	std::map<int, string> inst_decomp_;
	std::map<string, PoonProdNode> decomp_prod_);

	// each region is alloted a set of sum nodes
	std::vector<PoonSumNode> types_;

	// for MAP computation
	int defMapTypeIdx_;
	std::vector<string> mapDecomps_;
	double defMapSumPrb_;
	double defMapProdPrb_;

	//
	double invar_ = sqrt(20);

	static std::map<int, PoonRegion> id_regions_;

	// NOTE: dimension limited by range of 32-bit integer to around 215
	// for larger dimension, use long for id, or switch to string
	static int getRegionId(int a1, int a2, int b1, int b2, PoonParameter& params) {
		int id = ((a1*params.inputDim1_ + a2 - 1)*params.inputDim2_ + b1)*params.inputDim2_ + b2 - 1;
		if (id_regions_.find(id) == id_regions_.end())
			id_regions_[id] = PoonRegion(id, a1, a2, b1, b2, params);
		return id;
	}


	PoonRegion(int id, int a1, int a2, int b1, int b2, PoonParameter& params);




	~PoonRegion();
};

