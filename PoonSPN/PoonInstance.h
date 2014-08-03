#ifndef POONINSTANCE_H
#define POONINSTANCE_H

#include <vector>

// represent an image instance
class PoonInstance
{
public:
	std::vector<std::vector<double> >  vals_;
	double mean_;
	double std_;

	PoonInstance(){
		this->mean_ = 0; 
		this->std_ = 1;
	};
};



#endif

