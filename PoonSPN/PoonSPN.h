#pragma once

#include <vector>
#include <string>
#include <set>
#include "PoonInstance.h"
#include "PoonParameter.h"

using namespace std;


class PoonSPN{
private:
	PoonParameter params;

	bool isRecordingUpdate_ = true;	// record update in clearparse/setcurrparse
	vector<PoonInstance> trainingSet_;  //note that these should be wrapped into their own class so I can overload get instance

	// completion
	bool completeByMarginal_ = true; // complete pixel by marginal

	// root
	PoonSumNode root_;
	PoonRegion rootRegion_;

	// coarser resolution for larger regions
	int coarseDim1_, coarseDim2_;


public:
	PoonSPN(PoonParameter& params);
	~PoonSPN();



	// ----------------------------------------------------------
	// Bottom
	// ----------------------------------------------------------
	void PoonSPN::completeBottomImg(PoonInstance& inst);

	void PoonSPN::cmpMAPBottomHalf(PoonInstance& inst);

	// compute marginal by differentiation; see Darwiche-03 for details 
	void PoonSPN::cmpMAPBottomHalfMarginal(PoonInstance& inst);


	void PoonSPN::inferMAPBottomHalf(int ii, PoonInstance& inst);

	void PoonSPN::setMAPBottomToBuf(int instIdx, PoonInstance& inst);

	// ----------------------------------------------------------
	// Left
	// ----------------------------------------------------------
	void PoonSPN::completeLeftImg(PoonInstance& inst);

	void PoonSPN::cmpMAPLeftHalf(PoonInstance& inst);

	// compute marginal by differentiation; see Darwiche-03 for details
	void PoonSPN::cmpMAPLeftHalfMarginal(PoonInstance& inst);

	double PoonSPN::cmpMarginal(PoonRegion& r);

	void PoonSPN::setMAPLeftToBuf(int instIdx, PoonInstance& inst);

	// ----------------------------------------------------------
	// Learning
	// ----------------------------------------------------------
	void PoonSPN::init();

	// init: set mean/variance by equal quantiles from training for each pixel
	void PoonSPN::initUnitRegion(PoonRegion& r);

	void PoonSPN::clearUnusedInSPN();

	// ----------------------------------------------------------
	// Computation
	// ----------------------------------------------------------
	// derivative
	void PoonSPN::cmpDerivative();

	void PoonSPN::cmpDerivative(PoonRegion& r);

	void PoonSPN::initDerivative(PoonRegion& r);

	void PoonSPN::initDerivative();


	// evaluation: upward pass
	void PoonSPN::eval();

	void PoonSPN::eval(PoonRegion& r);

	// compute MAP
	void PoonSPN::inferMAPLeftHalf(int ii, PoonInstance& inst);
	void PoonSPN::inferMAPForLearning(int ii, PoonInstance& inst);

	// clear/set parses
	void PoonSPN::clearCurrParse(int ii);

	void PoonSPN::setCurrParseToMAP(int ii);

	void PoonSPN::setCurrParseFromBuf();
	void PoonSPN::clearCurrParseFromBuf();

	void PoonSPN::sendUpdate(int dest);

	void PoonSPN::recvUpdate(int src);

	// compute log probability
	double PoonSPN::llh(PoonInstance& inst);

	// set dspn input
	void PoonSPN::setInput(PoonInstance& inst);
	void PoonSPN::setInputOccludeLeftHalf(PoonInstance& inst);
	void PoonSPN::setInputOccludeBottomHalf(PoonInstance& inst);
	// -------------------------------------------------------------- //
	// load/save
	// -------------------------------------------------------------- //
	void PoonSPN::saveDSPN(string mdlName);
	void PoonSPN::loadDSPN(string mdlName, PoonSPN& spn);
	void PoonSPN::addChd(PoonRegion r, PoonSumNode n, string di, double cc);
	// -------------------------------------------------------------- //
	// utils
	// -------------------------------------------------------------- //
	void printParams();
}






};
