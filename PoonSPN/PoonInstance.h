#ifndef POONINSTANCE_H
#define POONINSTANCE_H

#include <vector>
using namespace std;

// represent an image instance
class PoonInstance
{
public:
	vector<vector<double> >  vals_;
	double mean_;
	double std_;

	PoonInstance(){this->mean_ = 0; this->std_ = 1;};
};



#endif

