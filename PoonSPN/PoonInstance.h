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

	PoonInstance(int dim1, int dim2){
		this->mean_ = 0; 
		this->std_ = 1;
		vals_.resize(dim1, std::vector<double>(dim2, 0.0)); //initialize
	};
};



#endif

