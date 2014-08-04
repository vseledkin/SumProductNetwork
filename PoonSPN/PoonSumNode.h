#ifndef POONSUMNODE_H
#define POONSUMNODE_H

#include "PoonNode.h"
#include "PoonParameter.h"

#include <map>
#include <string>
#include <memory>

class PoonSumNode : public PoonNode{
public:
	std::map<std::string, std::shared_ptr<PoonNode>> chds_;
	std::map<std::string, double> chdCnts_;
	double cnt_;

	PoonSumNode(std::shared_ptr<PoonParameter> param) {
		cnt_ = param->smoothSumCnt_;
	}

	PoonSumNode() {  //need a blank constuctor because the SPN has a sum node as the root; Annoying!
		SharedParams<PoonParameter> spp;
		auto params = spp.instance();
		cnt_ = params->smoothSumCnt_;
	};

	/*
	PoonSumNode(PoonSumNode&& other){ //move constructor'
		chds_(other.chds_.begin(), other.chds_.end());
		chdCnts_(other.chdCnts_);

	};*/

	void eval();

	void passDerivative();

	double getChdCnt(std::string di) {
		return chdCnts_.at(di);
	}

	void setChdCnt(std::string di, double cnt){
		chdCnts_[di] = cnt;
	}

	void removeChdOnly(std::string decompIdx, double cnt);

	void addChdOnly(std::string decompIdx, double cnt, std::shared_ptr<PoonNode> n);

};

#endif