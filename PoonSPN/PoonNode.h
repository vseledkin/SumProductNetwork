#ifndef POONNODE_H
#define POONNODE_H

#include <limits>  
#include <math.h> 

class PoonNode{
public:
	const double ZERO_LOGVAL_ = std::numeric_limits<double>::min(); //negative infinity
protected:
	double logval_;
	double logDerivative_;

public:
	PoonNode(){
		logDerivative_ = ZERO_LOGVAL_;
	};

	double getLogVal() { return logval_; }; //naming is bad getLogVal vs setVal
	double getLogDerVal() { return logDerivative_; }; //naming is bad getLogVal vs setVal
	void setVal(double v) { 
		if (v == 0){
			logval_ = ZERO_LOGVAL_;
		}
		else{
			logval_ = log(v);
		}
	};
	void setLogDerVal(double v){ logDerivative_ = v; };

	// evaluate root
	virtual void eval() = 0;

	// propagate derivative to children
	virtual void passDerivative() = 0;

	~PoonNode(){};
};

#endif