#ifndef POONREGION_H
#define POONREGION_H

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>
#include <iostream>

#include "PoonProdNode.h"
#include "PoonSumNode.h"
#include "PoonParameter.h"
#include "PoonInstance.h"

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
	std::map<int, std::string> inst_decomp_;
	std::map<std::string, std::shared_ptr<PoonProdNode> > decomp_prod_;

	// each region is alloted a set of sum nodes
	std::vector<std::shared_ptr<PoonSumNode>> types_;

	// for MAP computation
	int defMapTypeIdx_;
	std::vector<std::string> mapDecomps_;
	double defMapSumPrb_;
	double defMapProdPrb_;

	//
	double  invar_;

	static std::map<int, PoonRegion> id_regions_;



	PoonRegion(){}; //default construc
	
	
	PoonRegion& PoonRegion::operator=(PoonRegion& other){  //copy operator
		if (this != &other){
			id_ = other.id_;
			a1_ = other.a1_;
			a2_ = other.a2_;
			b1_ = other.b1_;
			b2_ = other.b2_;
			a_ = other.a_;
			b_ = other.b_;	// a=a2-a1, b=b2-b1	
			interval_ = other.interval_;	// for coarse resolution	

			// for pixel region only: gaussian units
			means_ = other.means_;
			vars_ = other.vars_;
			cnts_ = other.cnts_;
			ttlCnt_ = other.ttlCnt_;

			// data structure for a parse
			inst_type_ = other.inst_type_;
			inst_decomp_ = other.inst_decomp_;
			decomp_prod_ = other.decomp_prod_;

			// each region is alloted a set of sum nodes
			types_ = other.types_;

			// for MAP computation
			defMapTypeIdx_ = other.defMapTypeIdx_;
			mapDecomps_ = other.mapDecomps_;
			defMapSumPrb_ = other.defMapSumPrb_;
			defMapProdPrb_ = other.defMapProdPrb_;

			//
			invar_ = other.invar_;
		}

		return *this;
	};


	PoonRegion& PoonRegion::operator=(PoonRegion&& other){  //move operator, c++11
		if (this != &other){
			id_ = other.id_;
			a1_ = other.a1_;
			a2_ = other.a2_;
			b1_ = other.b1_;
			b2_ = other.b2_;
			a_ = other.a_;
			b_ = other.b_;	// a=a2-a1, b=b2-b1	
			interval_ = other.interval_;	// for coarse resolution	

			// for pixel region only: gaussian units
			means_ = other.means_;
			vars_ = other.vars_;
			cnts_ = other.cnts_;
			ttlCnt_ = other.ttlCnt_;

			// data structure for a parse
			inst_type_ = other.inst_type_;
			inst_decomp_ = other.inst_decomp_;
			decomp_prod_ = other.decomp_prod_;

			// each region is alloted a set of sum nodes
			types_ = other.types_;

			// for MAP computation
			defMapTypeIdx_ = other.defMapTypeIdx_;
			mapDecomps_ = other.mapDecomps_;
			defMapSumPrb_ = other.defMapSumPrb_;
			defMapProdPrb_ = other.defMapProdPrb_;

			//
			invar_ = other.invar_;
		}

		return *this;
	};

	PoonRegion(int id, int a1, int a2, int b1, int b2, std::shared_ptr<PoonParameter> params);

	// NOTE: dimension limited by range of 32-bit integer to around 215
	// for larger dimension, use long for id, or switch to string
	static int getRegionId(int a1, int a2, int b1, int b2, std::shared_ptr<PoonParameter> params) {
		int id = ((a1*params->inputDim1_ + a2 - 1) * params->inputDim2_ + b1) * params->inputDim2_ + b2 - 1;
		if (id_regions_.find(id) == id_regions_.end()){
			id_regions_[id] = PoonRegion(id, a1, a2, b1, b2, params);
			//std::cout << "new region" << std::endl;
		}
		return id;
	}

	static PoonRegion& getRegion(int id, std::shared_ptr<PoonParameter> params) {
		
		if (id_regions_.find(id) == id_regions_.end()) {
			std::cout << "not found" << std::endl;
			int b2 = id % params->inputDim2_ + 1;
			int x = id / params->inputDim2_;
			int b1 = x % params->inputDim2_;
			x = x / params->inputDim2_;
			int a2 = x % params->inputDim1_ + 1;
			int a1 = x / params->inputDim1_;
			return getRegion(getRegionId(a1, a2, b1, b2, params), params);  //this might now work, check, weird recussion

		}else{
			//std::cout << "found" << std::endl;
			return id_regions_.at(id);
		}
	}

	int getId() { return id_; }

	std::string myStr() {
		std::ostringstream s;
		s << "< " << a1_ << ", " << a2_ << ", " << b1_ << ", " << b2_ << " >";
		return s.str();
	}

	// initialization
	void resetTypes(int numTypes, std::shared_ptr<PoonParameter> param);

	void setTypes(int numTypes, std::shared_ptr<PoonParameter> params);

	// set value for input layer	
	void setBase(double val) {
		setBaseGauss(val);
	}

	double cmpGauss(double v, double mean) {
		double m = mean - v;
		return -(m*m / 2);
	}

	void setBaseGauss(double v);

	void setBaseForSumOut(std::shared_ptr<PoonParameter> params);

	// compute MAP state at inference time
	void inferMAP(int instIdx, PoonInstance& inst);

	// compute MAP state at learning time: could tap a previous unused node
	void inferMAPForLearning(int instIdx, PoonInstance& inst, std::shared_ptr<PoonParameter> params);

	// downward trace-back step
	void setCurrParseToMAP(int instIdx, std::shared_ptr<PoonParameter> params);

	// clear an existing parse for incremental EM 
	void clearCurrParse(int instIdx, std::shared_ptr<PoonParameter> params);

	// clear parse from other slaves
	void clearCurrParseFromBuf(int chosenType, int ri1, int ri2, int ti1, int ti2, std::shared_ptr<PoonParameter>params);

	void setCurrParseFromBuf(int chosenType, int ri1, int ri2, int ti1, int ti2, std::shared_ptr<PoonParameter> params);


};

#endif