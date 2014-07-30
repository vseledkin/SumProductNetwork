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

	PoonSumNode(PoonParameter& param) {
		cnt_ = param.smoothSumCnt_;
	}

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

