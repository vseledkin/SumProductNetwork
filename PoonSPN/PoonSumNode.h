#pragma once
#include "PoonNode.h"
#include "PoonParameter.h"
#include <map>
#include <string>

class PoonSumNode : public PoonNode
{
public:
	std::map<string, PoonNode> chds_;
	std::map<string, double> chdCnts_;
	double cnt_;

	PoonSumNode(std::shared_ptr<PoonParameter> param) {
		cnt_ = param->smoothSumCnt_;
	}

	PoonSumNode() {  //need a blank constuctor because the SPN has a sum node as the root; Annoying!
		SharedParams<PoonParameter> spp;
		auto params = spp.instance();
		cnt_ = params->smoothSumCnt_;
	};

	void eval();

	void passDerivative();

	double getChdCnt(string di) {
		return chdCnts_.at(di);
	}

	void setChdCnt(string di, double cnt){
		chdCnts_[di] = cnt;
	}

	void removeChdOnly(string decompIdx, double cnt);

	void PoonSumNode::addChdOnly(string decompIdx, double cnt, PoonNode& n);

	void PoonSumNode::removeChdOnly(string decompIdx, double cnt);
};

