
#include "PoonProdNode.h"
#include "Utils/Utils.h"

using namespace std;


void PoonProdNode::passDerivative() {

	if (this->logDerivative_ == this->ZERO_LOGVAL_) return;

	if (this->logval_ == this->ZERO_LOGVAL_) {
		int cnt = 0;
		for (auto n : chds_){
			if (n->getLogVal() == this->ZERO_LOGVAL_) {
				cnt++; 
				if (cnt > 1)
					return; 
			}
		}
	}

	for (auto n : chds_) {
		if (n->getLogVal() == this->ZERO_LOGVAL_) {
			double l = 0;
			for (auto m : chds_) {
				if (m->getLogVal() != this->ZERO_LOGVAL_) 
					l += m->getLogVal();
			}
			l += logDerivative_;
			if (n->getLogVal() == this->ZERO_LOGVAL_)
				n->setVal(l);
			else 
				n->setVal(addLog(n->getLogVal(), l));
		}
		else if (logval_ != this->ZERO_LOGVAL_) {
			double l = this->logDerivative_ + logval_ - n->getLogVal();
			if (n->getLogDerVal() == this->ZERO_LOGVAL_)
				n->setLogDerVal(l);
			else
				n->setLogDerVal(addLog(n->getLogDerVal(), l));
		}
	}
}

void PoonProdNode::eval() {
	logval_ = 0;
	for (auto n : chds_) {
		double v = n->getLogVal();
		if (v == this->ZERO_LOGVAL_) {
			logval_ = this->ZERO_LOGVAL_; 
			return; 
		}
		logval_ += v;
	}
}

void PoonProdNode::addChd(shared_ptr<PoonNode> n) {
	this->chds_.push_back(n);
}