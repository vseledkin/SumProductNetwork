

#include <limits>  
#include <math.h> 

//an abstract node class
class PoonNode
{
public:
	const double ZERO_LOGVAL_ = std::numeric_limits<double>::min(); //negative infinity
protected:
	double logval_;
	double logDerivative_;

public:
	PoonNode(){
		logDerivative_ = ZERO_LOGVAL_;
	};
	~PoonNode();

	double getLogVal() { return logval_; }; //naming is bad getLogVal vs setVal
	void setVal(double v) { 
		if (v == 0){
			logval_ = ZERO_LOGVAL_;
		}
		else{
			logval_ = log(v);
		}
	};

	// evaluate root
	virtual void eval() = 0;

	// propagate derivative to children
	virtual void passDerivative() = 0;
};

