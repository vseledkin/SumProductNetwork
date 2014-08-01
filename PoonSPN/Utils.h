
#include <iostream>
#include <cstdio>
#include <ctime>
#include <math.h> 

#include "PoonInstance.h"


//old poon funtions
int getIntVal(PoonInstance& inst, double p) {
	return (int)(p*inst.std_ + inst.mean_);
};

static double addLog(double l1, double l2) {
	if (l1>l2) {
		return l1 + log(1.0 + exp(l2 - l1));
	}
	else {
		return l2 + log(1.0 + exp(l1 - l2));
	}
};


