
#include "PoonSumNode.h"
#include <cmath>
#include "Utils/Utils.h"

using namespace std;

void PoonSumNode::eval() {
	double v = 0;
	string maxi;
	double maxl = 0;
	for (auto& kv : chds_) { //c++11 key value pairs
		string i = kv.first;
		double l = kv.second->getLogVal();
		if (l != this->ZERO_LOGVAL_){
			if (maxi.empty() || maxl < l) {
				maxi = i; 
				maxl = l;
			}
		}
	}
	if (maxi.empty()){
		this->setVal(this->ZERO_LOGVAL_);
		return; 
	}
	for (auto& kv : chds_) {
		string i = kv.first;
		if (chdCnts_.find(i) == chdCnts_.end())
			continue;
		double l = l = kv.second->getLogVal();
		if (l == this->ZERO_LOGVAL_) 
			continue;
		v += getChdCnt(i)*exp(l - maxl);
	}
	this->setVal(log(v / cnt_) + maxl);
}

void PoonSumNode::passDerivative() {
	if (logDerivative_ == this->ZERO_LOGVAL_) 
		return;

	for (auto& kv : chds_) {
		string di = kv.first;
		auto n = kv.second;
		double l = logDerivative_ + log(getChdCnt(di) / cnt_);
		if (n->getLogDerVal() == this->ZERO_LOGVAL_)
			n->setLogDerVal(l);
		else
			n->setLogDerVal(addLog(l, n->getLogVal()));
	}
}

void PoonSumNode::addChdOnly(string decompIdx, double cnt, shared_ptr<PoonNode> n) {
	if (chds_.find(decompIdx) == chds_.end()) {
		chds_[decompIdx] = n;
	}
	if (chdCnts_.find(decompIdx) == chdCnts_.end()) {
		chdCnts_[decompIdx] = cnt;
	}
	else 
		chdCnts_[decompIdx] = cnt + getChdCnt(decompIdx);
	cnt_ += cnt;
}

void PoonSumNode::removeChdOnly(string decompIdx, double cnt) {
	double cc = getChdCnt(decompIdx);
	cc -= cnt;
	if (cc == 0) { //double test if zero, bad
		chds_.erase(decompIdx);
		chdCnts_.erase(decompIdx);
	}
	else 
		chdCnts_[decompIdx] = cc;
	cnt_ -= cnt;
}