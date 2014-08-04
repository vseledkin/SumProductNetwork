#include "PoonSPN.h"
#include "Utils/Utils.h"
#include "PoonDecomposition.h"

#include <algorithm>  
#include <set>

using namespace std;

bool PoonSPN::isRecordingUpdate_ = true;

// ----------------------------------------------------------
// Bottom
// ----------------------------------------------------------
void PoonSPN::completeBottomImg(PoonInstance& inst) {
	//Utils.logTimeMS("before complete bottom half");
	if (completeByMarginal_) {
		cmpMAPBottomHalfMarginal(inst);
		//Utils.logTimeMS("Complete bottom by Marginal");
	}
	else {
		cmpMAPBottomHalf(inst);
		//Utils.logTimeMS("Complete bottom by MPE");
	}
};

void PoonSPN::cmpMAPBottomHalf(PoonInstance& inst) {
	isRecordingUpdate_ = false; // inference now; no need to update count
	int cmpIdx = -1;	// temp idx for cmpMAP				
	inferMAPBottomHalf(cmpIdx, inst);
	setCurrParseToMAP(cmpIdx);
	setMAPBottomToBuf(cmpIdx, inst);
	clearCurrParse(cmpIdx);
	isRecordingUpdate_ = true;
};

// compute marginal by differentiation; see Darwiche-03 for details 
void PoonSPN::cmpMAPBottomHalfMarginal(PoonInstance& inst) {
	setInputOccludeBottomHalf(inst);
	eval();
	cmpDerivative();

	for (int i = 0; i < params->inputDim1_ / 2; i++) {
		for (int j = 0; j < params->inputDim2_; j++)
			params->buf_int_[params->buf_idx_++] = getIntVal(inst, inst.vals_[i][j]);
	}
	for (int i = params->inputDim1_ / 2; i < params->inputDim1_; i++) {
		for (int j = 0; j < params->inputDim2_; j++) {
			int ri = PoonRegion::getRegionId(i, i + 1, j, j + 1, params);
			PoonRegion& r = PoonRegion::getRegion(ri, params);
			double p = cmpMarginal(r);
			params->buf_int_[params->buf_idx_++] = getIntVal(inst, p);//(int)(p*255);
		}
	}
};


void PoonSPN::inferMAPBottomHalf(int ii, PoonInstance& inst) {
	setInputOccludeBottomHalf(inst);

	// finePoonRegion 
	for (int ca = 0; ca < coarseDim1_; ca++){
		for (int cb = 0; cb < coarseDim2_; cb++){
			for (int a = 1; a <= params->baseResolution_; a++){
				for (int b = 1; b <= params->baseResolution_; b++)
				{
					if (a == 1 && b == 1) continue;
					for (int a1 = ca*params->baseResolution_; a1 <= (ca + 1)*params->baseResolution_ - a; a1++) {
						int a2 = a1 + a;
						for (int b1 = cb*params->baseResolution_; b1 <= (cb + 1)*params->baseResolution_ - b; b1++) {
							int b2 = b1 + b;
							int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
							PoonRegion& r = PoonRegion::getRegion(ri, params);
							r.inferMAP(ii, inst);
						}
					}
				}
			}
		}
	}

	// coarsePoonRegion
	for (int ca = 1; ca <= coarseDim1_; ca++){
		for (int cb = 1; cb <= coarseDim2_; cb++) {
			if (ca == 1 && cb == 1) continue;	// taken care of below in fine

			for (int a1 = 0; a1 <= params->inputDim1_ - ca*params->baseResolution_; a1 += params->baseResolution_) {
				int a2 = a1 + ca*params->baseResolution_;
				for (int b1 = 0; b1 <= params->inputDim2_ - cb*params->baseResolution_; b1 += params->baseResolution_) {
					int b2 = b1 + cb*params->baseResolution_;
					int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
					PoonRegion& r = PoonRegion::getRegion(ri, params);
					r.inferMAP(ii, inst);
				}
			}
		}
	}
}

void PoonSPN::setMAPBottomToBuf(int instIdx, PoonInstance& inst) {
	//
	for (int a1 = 0; a1 <= params->inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= params->inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			if (a1 < params->inputDim2_ / 2) {
				params->buf_int_[params->buf_idx_++] = getIntVal(inst, inst.vals_[a1][b1]);
			}
			else {
				int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
				PoonRegion& r = PoonRegion::getRegion(ri, params);
				int vi = r.inst_type_.at(instIdx);
				params->buf_int_[params->buf_idx_++] = getIntVal(inst, r.means_[vi]);
			}
		}
	}
};

// ----------------------------------------------------------
// Left
// ----------------------------------------------------------
void PoonSPN::completeLeftImg(PoonInstance& inst) {
	//Utils.logTimeMS("before complete left half");
	if (completeByMarginal_) {
		cmpMAPLeftHalfMarginal(inst);
		//Utils.logTimeMS("Complete left by Marginal");
	}
	else {
		cmpMAPLeftHalf(inst);
		//Utils.logTimeMS("Complete left by MAP");
	}
};

void PoonSPN::cmpMAPLeftHalf(PoonInstance& inst) {
	isRecordingUpdate_ = false;
	int cmpIdx = -1;	// temp idx for cmpMAP
	inferMAPLeftHalf(cmpIdx, inst);
	setCurrParseToMAP(cmpIdx);
	setMAPLeftToBuf(cmpIdx, inst);
	clearCurrParse(cmpIdx);
	isRecordingUpdate_ = true;
};

// compute marginal by differentiation; see Darwiche-03 for details
void PoonSPN::cmpMAPLeftHalfMarginal(PoonInstance& inst) {
	setInputOccludeLeftHalf(inst);
	eval();
	cmpDerivative();

	for (int i = 0; i < params->inputDim1_; i++) {
		for (int j = 0; j < params->inputDim2_ / 2; j++) {
			int ri = PoonRegion::getRegionId(i, i + 1, j, j + 1, params);
			PoonRegion& r = PoonRegion::getRegion(ri, params);
			double p = cmpMarginal(r);
			params->buf_int_[params->buf_idx_++] = getIntVal(inst, p);
		}
		for (int j = params->inputDim2_ / 2; j < params->inputDim2_; j++)
			params->buf_int_[params->buf_idx_++] = getIntVal(inst, inst.vals_[i][j]);
	}
};

double PoonSPN::cmpMarginal(PoonRegion& r) {
	double t = 0, d = 0;
	double md = 100;

	for (int i = 0; i < r.types_.size(); i++) {
		auto n = r.types_.at(i);
		if (n->getLogDerVal() == n->ZERO_LOGVAL_)
			continue;
		if (md == 100 || n->getLogDerVal() > md)
			md = n->getLogDerVal();
	}
	for (int i = 0; i < r.types_.size(); i++) {
		auto n = r.types_.at(i);
		if (n->getLogDerVal() == n->ZERO_LOGVAL_)
			continue;
		double p = exp(n->getLogDerVal() - md);
		d += r.means_[i] * p; t += p;
	}
	d /= t;
	return d;
};

void PoonSPN::setMAPLeftToBuf(int instIdx, PoonInstance& inst) {
	for (int a1 = 0; a1 <= params->inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= params->inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			if (b1 >= params->inputDim2_ / 2) {
				params->buf_int_[params->buf_idx_++] = getIntVal(inst, inst.vals_[a1][b1]);
			}
			else {
				int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
				PoonRegion& r = PoonRegion::getRegion(ri, params);
				int vi = r.inst_type_.at(instIdx);
				params->buf_int_[params->buf_idx_++] = getIntVal(inst, r.means_[vi]);
			}
		}
	}
};

// ----------------------------------------------------------
// Learning
// ----------------------------------------------------------
void PoonSPN::init() {
	// coarsePoonRegion
	for (int ca = 1; ca <= coarseDim1_; ca++){
		for (int cb = 1; cb <= coarseDim2_; cb++) {
			if (ca == 1 && cb == 1) continue;	// taken care of below in fine

			for (int a1 = 0; a1 <= params->inputDim1_ - ca*params->baseResolution_; a1 += params->baseResolution_) {
				int a2 = a1 + ca*params->baseResolution_;
				for (int b1 = 0; b1 <= params->inputDim2_ - cb*params->baseResolution_; b1 += params->baseResolution_) {
					int b2 = b1 + cb*params->baseResolution_;

					// coarsePoonRegions
					int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
					PoonRegion& r = PoonRegion::getRegion(ri, params);
					if (ca == coarseDim1_ && cb == coarseDim2_) {
						r.resetTypes(1, params);	// one sum node as root						
						rootRegion_ = r;
						root_ = r.types_.at(0);
					}
					else r.resetTypes(params->numSumPerRegion_, params);
				}
			}
		}
	}

	// finePoonRegion 
	for (int ca = 0; ca < coarseDim1_; ca++){
		for (int cb = 0; cb < coarseDim2_; cb++){
			for (int a = 1; a <= params->baseResolution_; a++){
				for (int b = 1; b <= params->baseResolution_; b++){
					for (int a1 = ca*params->baseResolution_; a1 <= (ca + 1)*params->baseResolution_ - a; a1++) {
						int a2 = a1 + a;
						for (int b1 = cb*params->baseResolution_; b1 <= (cb + 1)*params->baseResolution_ - b; b1++) {
							int b2 = b1 + b;
							int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
							PoonRegion r = PoonRegion::getRegion(ri, params);
							if (a == 1 && b == 1) {
								initUnitRegion(r);
							}
							else r.resetTypes(params->numSumPerRegion_, params);
						}
					}
				}
			}
		}
	}
};

// init: set mean/variance by equal quantiles from training for each pixel
//why is this not part of a constructor for region?
void PoonSPN::initUnitRegion(PoonRegion& r) {
	r.resetTypes(params->numComponentsPerVar_, params);

	r.means_.clear();
	r.means_.resize(params->numComponentsPerVar_);
	r.vars_.clear();
	r.vars_.resize(params->numComponentsPerVar_);
	r.cnts_.clear();
	r.cnts_.resize(params->numComponentsPerVar_);


	int ttlCnt = trainingSet_.size();
	int cnt = (int)ceil(ttlCnt*1.0 / params->numComponentsPerVar_);

	vector<double> vals(ttlCnt, 0.0);
	for (int ii = 0; ii < trainingSet_.size(); ii++) {
		vals[ii] = trainingSet_.at(ii).vals_[r.a1_][r.b1_];
	}

	sort(vals.begin(), vals.end());
	for (int bi = 0; bi < params->numComponentsPerVar_; bi++) {
		int ac = 0;
		for (int ii = bi*cnt; ii < (bi + 1)*cnt && ii < ttlCnt; ii++, ac++) {
			r.means_[bi] += vals[ii];
			r.vars_[bi] += vals[ii] * vals[ii];
		}
		r.means_[bi] /= ac;
		r.vars_[bi] /= ac; r.vars_[bi] -= r.means_[bi] * r.means_[bi];
		r.cnts_[bi] = ac;
	}
	r.ttlCnt_ = ttlCnt;
};

void PoonSPN::clearUnusedInSPN() {
	// coarse
	for (int ca = 1; ca <= coarseDim1_; ca++){
		for (int cb = 1; cb <= coarseDim2_; cb++) {
			if (ca == 1 && cb == 1) 
				continue;	// taken care of below in fine
			for (int a1 = 0; a1 <= params->inputDim1_ - ca*params->baseResolution_; a1 += params->baseResolution_) {
				int a2 = a1 + ca*params->baseResolution_;
				for (int b1 = 0; b1 <= params->inputDim2_ - cb*params->baseResolution_; b1 += params->baseResolution_) {
					int b2 = b1 + cb*params->baseResolution_;
					int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
					PoonRegion& r = PoonRegion::getRegion(ri, params);

					set<string> decomps;
					for (auto n : r.types_) {
						if (n->chds_.size() > 0) {
							double tc = 0;
							for (auto& kv : n->chdCnts_) {
								tc += n->chdCnts_.at(kv.first);  // could just be kv
								decomps.insert(kv.first);
							}
						}
					}

					// clear dead decomp_prod
					set<string> deadDecomps;
					for (auto& kv : r.decomp_prod_) {
						if (decomps.find(kv.first) == decomps.end()) {
							deadDecomps.insert(kv.first);
							continue;
						}
					}

					for (string di : deadDecomps) {
						r.decomp_prod_.erase(di);
						PoonDecomposition::remove(di);
					}
				}
			}
		}
	}

	// finePoonRegion
	for (int ca = 0; ca < coarseDim1_; ca++){
		for (int cb = 0; cb < coarseDim2_; cb++){
			for (int a = 1; a <= params->baseResolution_; a++){
				for (int b = 1; b <= params->baseResolution_; b++)
				{
					for (int a1 = ca*params->baseResolution_; a1 <= (ca + 1)*params->baseResolution_ - a; a1++) {
						int a2 = a1 + a;
						for (int b1 = cb*params->baseResolution_; b1 <= (cb + 1)*params->baseResolution_ - b; b1++) {
							int b2 = b1 + b;
							int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
							PoonRegion& r = PoonRegion::getRegion(ri, params);

							// clear dead decomp_prod
							set<string> decomps;
							for (auto n : r.types_) {
								if (n->chds_.size() > 0) {
									for (auto& kv : n->chdCnts_) {
										decomps.insert(kv.first);
									}
								}
							}
							set<string> deadDecomps;
							for (auto& kv : r.decomp_prod_) {
								if (decomps.find(kv.first) == decomps.end()) {
									deadDecomps.insert(kv.first);
									continue;
								}
							}

							for (string di : deadDecomps) {
								r.decomp_prod_.erase(di);
								PoonDecomposition::remove(di);
							}
						}
					}
				}
			}
		}
	}
};

// ----------------------------------------------------------
// Computation
// ----------------------------------------------------------
// derivative
void PoonSPN::cmpDerivative() {
	initDerivative();

	root_->setLogDerVal(0.0);
	root_->passDerivative();
	for (auto& kv : root_->chds_) {
		auto n = root_->chds_.at(kv.first);
		n->passDerivative();
	}

	// coarsePoonRegion
	for (int ca = coarseDim1_; ca >= 1; ca--){
		for (int cb = coarseDim2_; cb >= 1; cb--) {
			if (ca == 1 && cb == 1) continue;	// taken care of below in fine

			for (int a1 = 0; a1 <= params->inputDim1_ - ca*params->baseResolution_; a1 += params->baseResolution_) {
				int a2 = a1 + ca*params->baseResolution_;
				for (int b1 = 0; b1 <= params->inputDim2_ - cb*params->baseResolution_; b1 += params->baseResolution_) {
					int b2 = b1 + cb*params->baseResolution_;

					// coarsePoonRegions
					int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
					PoonRegion& r = PoonRegion::getRegion(ri, params);
					cmpDerivative(r);
				}
			}
		}
	}


	// finePoonRegion 
	for (int ca = coarseDim1_ - 1; ca >= 0; ca--){
		for (int cb = coarseDim2_ - 1; cb >= 0; cb--){
			for (int a = params->baseResolution_; a >= 1; a--){
				for (int b = params->baseResolution_; b >= 1; b--)
				{
					if (a == 1 && b == 1) continue;	// take care in setInput
					for (int a1 = ca*params->baseResolution_; a1 <= (ca + 1)*params->baseResolution_ - a; a1++) {
						int a2 = a1 + a;
						for (int b1 = cb*params->baseResolution_; b1 <= (cb + 1)*params->baseResolution_ - b; b1++) {
							int b2 = b1 + b;
							int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
							PoonRegion& r = PoonRegion::getRegion(ri, params);
							cmpDerivative(r);
						}
					}
				}
			}
		}
	}
};

void PoonSPN::cmpDerivative(PoonRegion& r) {
	for (int i = 0; i < r.types_.size(); i++) {
		r.types_.at(i)->passDerivative();
	}
	for (auto& kv : r.decomp_prod_) {
		auto n = r.decomp_prod_.at(kv.first);
		n->passDerivative();
	}
};

void PoonSPN::initDerivative(PoonRegion& r) {
	for (auto& kv : r.decomp_prod_) {
		auto n = r.decomp_prod_.at(kv.first);
		n->setLogDerVal(n->ZERO_LOGVAL_);
	}
	for (auto n : r.types_) {
		n->setLogDerVal(n->ZERO_LOGVAL_);
	}
};

void PoonSPN::initDerivative() {
	for (int ca = coarseDim1_; ca >= 1; ca--){
		for (int cb = coarseDim2_; cb >= 1; cb--) {
			if (ca == 1 && cb == 1) continue;	// taken care of below in fine

			for (int a1 = 0; a1 <= params->inputDim1_ - ca*params->baseResolution_; a1 += params->baseResolution_) {
				int a2 = a1 + ca*params->baseResolution_;
				for (int b1 = 0; b1 <= params->inputDim2_ - cb*params->baseResolution_; b1 += params->baseResolution_) {
					int b2 = b1 + cb*params->baseResolution_;

					// coarsePoonRegions
					int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
					PoonRegion r = PoonRegion::getRegion(ri, params);
					initDerivative(r);
				}
			}
		}
	}

	// finePoonRegion 
	for (int ca = coarseDim1_ - 1; ca >= 0; ca--){
		for (int cb = coarseDim2_ - 1; cb >= 0; cb--){
			for (int a = params->baseResolution_; a >= 1; a--){
				for (int b = params->baseResolution_; b >= 1; b--)
				{
					for (int a1 = ca*params->baseResolution_; a1 <= (ca + 1)*params->baseResolution_ - a; a1++) {
						int a2 = a1 + a;
						for (int b1 = cb*params->baseResolution_; b1 <= (cb + 1)*params->baseResolution_ - b; b1++) {
							int b2 = b1 + b;
							int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
							PoonRegion& r = PoonRegion::getRegion(ri, params);
							initDerivative(r);
						}
					}
				}
			}
		}
	}
};


// evaluation: upward pass
void PoonSPN::eval() {
	// finePoonRegion 
	for (int ca = 0; ca < coarseDim1_; ca++){
		for (int cb = 0; cb < coarseDim2_; cb++){
			for (int a = 1; a <= params->baseResolution_; a++)
			for (int b = 1; b <= params->baseResolution_; b++)
			{
				if (a == 1 && b == 1) continue;	// take care in setInput
				for (int a1 = ca*params->baseResolution_; a1 <= (ca + 1)*params->baseResolution_ - a; a1++) {
					int a2 = a1 + a;
					for (int b1 = cb*params->baseResolution_; b1 <= (cb + 1)*params->baseResolution_ - b; b1++) {
						int b2 = b1 + b;
						int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
						PoonRegion& r = PoonRegion::getRegion(ri, params);
						eval(r);
					}
				}
			}
		}
	}

	// coarsePoonRegion
	for (int ca = 1; ca <= coarseDim1_; ca++){
		for (int cb = 1; cb <= coarseDim2_; cb++) {
			if (ca == 1 && cb == 1)
				continue;

			for (int a1 = 0; a1 <= params->inputDim1_ - ca*params->baseResolution_; a1 += params->baseResolution_) {
				int a2 = a1 + ca*params->baseResolution_;
				for (int b1 = 0; b1 <= params->inputDim2_ - cb*params->baseResolution_; b1 += params->baseResolution_) {
					int b2 = b1 + cb*params->baseResolution_;

					// coarsePoonRegions
					int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
					PoonRegion r = PoonRegion::getRegion(ri, params);
					eval(r);
				}
			}
		}
	}
};

void PoonSPN::eval(PoonRegion& r) {
	for (auto& kv : r.decomp_prod_) {
		auto n = r.decomp_prod_.at(kv.first);
		n->eval();
	}
	for (auto n : r.types_) {
		if (n->chds_.size() > 0)
			n->eval();
		else
			n->setVal(n->ZERO_LOGVAL_);
	}
};

// compute MAP
void PoonSPN::inferMAPLeftHalf(int ii, PoonInstance& inst) {
	setInputOccludeLeftHalf(inst);

	// finePoonRegion 
	for (int ca = 0; ca < coarseDim1_; ca++){
		for (int cb = 0; cb < coarseDim2_; cb++){
			for (int a = 1; a <= params->baseResolution_; a++){
				for (int b = 1; b <= params->baseResolution_; b++){
					{
						if (a == 1 && b == 1) continue;
						for (int a1 = ca*params->baseResolution_; a1 <= (ca + 1)*params->baseResolution_ - a; a1++) {
							int a2 = a1 + a;
							for (int b1 = cb*params->baseResolution_; b1 <= (cb + 1)*params->baseResolution_ - b; b1++) {
								int b2 = b1 + b;
								int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
								PoonRegion r = PoonRegion::getRegion(ri, params);
								r.inferMAP(ii, inst);
							}
						}
					}
				}
			}
		}
	}

	// coarsePoonRegion
	for (int ca = 1; ca <= coarseDim1_; ca++){
		for (int cb = 1; cb <= coarseDim2_; cb++) {
			if (ca == 1 && cb == 1) 
				continue;	// taken care of below in fine

			for (int a1 = 0; a1 <= params->inputDim1_ - ca*params->baseResolution_; a1 += params->baseResolution_) {
				int a2 = a1 + ca*params->baseResolution_;
				for (int b1 = 0; b1 <= params->inputDim2_ - cb*params->baseResolution_; b1 += params->baseResolution_) {
					int b2 = b1 + cb*params->baseResolution_;
					int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
					PoonRegion& r = PoonRegion::getRegion(ri, params);
					r.inferMAP(ii, inst);
				}
			}
		}
	}
};

void PoonSPN::inferMAPForLearning(int ii, PoonInstance& inst) {
	setInput(inst);

	// finePoonRegion 
	for (int ca = 0; ca < coarseDim1_; ca++){
		for (int cb = 0; cb < coarseDim2_; cb++){
			for (int a = 1; a <= params->baseResolution_; a++){
				for (int b = 1; b <= params->baseResolution_; b++){
					if (a == 1 && b == 1)
						continue;
					for (int a1 = ca*params->baseResolution_; a1 <= (ca + 1)*params->baseResolution_ - a; a1++) {
						int a2 = a1 + a;
						for (int b1 = cb*params->baseResolution_; b1 <= (cb + 1)*params->baseResolution_ - b; b1++) {
							int b2 = b1 + b;
							int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
							PoonRegion& r = PoonRegion::getRegion(ri, params);
							r.inferMAPForLearning(ii, inst, params);
						}
					}
				}
			}
		}
	}

	// coarsePoonRegion
	for (int ca = 1; ca <= coarseDim1_; ca++){
		for (int cb = 1; cb <= coarseDim2_; cb++) {
			if (ca == 1 && cb == 1) continue;	// taken care of below in fine

			for (int a1 = 0; a1 <= params->inputDim1_ - ca*params->baseResolution_; a1 += params->baseResolution_) {
				int a2 = a1 + ca*params->baseResolution_;
				for (int b1 = 0; b1 <= params->inputDim2_ - cb*params->baseResolution_; b1 += params->baseResolution_) {
					int b2 = b1 + cb*params->baseResolution_;
					int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
					PoonRegion& r = PoonRegion::getRegion(ri, params);
					r.inferMAPForLearning(ii, inst, params);
				}
			}
		}
	}
};

// clear/set parses
void PoonSPN::clearCurrParse(int ii) {
	rootRegion_.clearCurrParse(ii, params);
};

void PoonSPN::setCurrParseToMAP(int ii) {
	rootRegion_.setCurrParseToMAP(ii, params);
};

void PoonSPN::setCurrParseFromBuf() {
	// --- update format: instId,PoonRegionId, type, decomp(rid1, rid2, type1, type2) - in buf_inf				
	int k = 0;
	while (k < params->buf_idx_) {
		int ri = params->buf_int_[k++];
		int chosenType = params->buf_int_[k++];
		int ri1 = params->buf_int_[k++];
		int ri2 = params->buf_int_[k++];
		int ti1 = params->buf_int_[k++];
		int ti2 = params->buf_int_[k++];
		PoonRegion& r = PoonRegion::getRegion(ri, params);
		r.setCurrParseFromBuf(chosenType, ri1, ri2, ti1, ti2, params);
	}
};
void PoonSPN::clearCurrParseFromBuf() {
	// --- update format: instId,PoonRegionId, type, decomp(rid1, rid2, type1, type2) - in buf_inf				
	int k = 0;
	while (k < params->buf_idx_) {
		int ri = params->buf_int_[k++];
		int chosenType = params->buf_int_[k++];
		int ri1 = params->buf_int_[k++];
		int ri2 = params->buf_int_[k++];
		int ti1 = params->buf_int_[k++];
		int ti2 = params->buf_int_[k++];
		PoonRegion& r = PoonRegion::getRegion(ri, params);
		r.clearCurrParseFromBuf(chosenType, ri1, ri2, ti1, ti2, params);
	}
};

void PoonSPN::sendUpdate(int dest) {
	/*
	if (params->buf_idx_ >= params->buf_size_) 
		cout << "ERR: buffer overflow to " << dest << endl;
	MPI.COMM_WORLD.Send(params->buf_int_, 0, params->buf_idx_, MPI.INT, dest, 0);
	*/
};

void PoonSPN::recvUpdate(int src) {
	/*
	Status status = MPI.COMM_WORLD.Recv(params->buf_int_, params->buf_idx_, params->buf_size_, MPI.INT, src, 0);
	params->buf_idx_ += status.count;
	if (params->buf_idx_ >= params->buf_size_)
		cout << "ERR: buffer overflow from " << src << endl;
	*/
};

// compute log probability
double PoonSPN::llh(PoonInstance& inst) {
	setInput(inst);
	eval();
	return root_->getLogVal();
};

// set dspn input
void PoonSPN::setInput(PoonInstance& inst) {
	for (int a1 = 0; a1 <= params->inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= params->inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
			PoonRegion& r = PoonRegion::getRegion(ri, params);
			r.setBase(inst.vals_[a1][b1]);
		}
	}
};

void PoonSPN::setInputOccludeLeftHalf(PoonInstance& inst) {
	for (int a1 = 0; a1 <= params->inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= params->inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
			PoonRegion& r = PoonRegion::getRegion(ri, params);
			if (b1 < params->inputDim2_ / 2) //r.setBase(0,0);	// log 1,1
				r.setBaseForSumOut(params);
			else
				r.setBase(inst.vals_[a1][b1]);
		}
	}
};

void PoonSPN::setInputOccludeBottomHalf(PoonInstance& inst) {
	for (int a1 = 0; a1 <= params->inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= params->inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			int ri = PoonRegion::getRegionId(a1, a2, b1, b2, params);
			PoonRegion& r = PoonRegion::getRegion(ri, params);
			if (a1 >= params->inputDim1_ / 2) //r.setBase(0,0);	// log 1,1
				r.setBaseForSumOut(params);
			else
				r.setBase(inst.vals_[a1][b1]);
		}
	}
};

// -------------------------------------------------------------- //
// load/save
// -------------------------------------------------------------- //
//void PoonSPN::saveDSPN(string mdlName){}
//void PoonSPN::loadDSPN(string mdlName, PoonSPN& spn){}



void PoonSPN::addChd(PoonRegion r, PoonSumNode n, string di, double cc) {
	n.setChdCnt(di, cc);
	
	if (r.decomp_prod_.find(di) == r.decomp_prod_.end()) {


		PoonDecomposition d = PoonDecomposition::getDecomposition(di);
		auto np = make_shared<PoonProdNode>();
		r.decomp_prod_[di] = np;
		PoonRegion& r1 = PoonRegion::getRegion(d.regionId1_, params);
		PoonRegion& r2 = PoonRegion::getRegion(d.regionId2_, params);

		np->addChd(r1.types_.at(d.typeId1_));
		np->addChd(r2.types_.at(d.typeId2_));
		n.chds_[di] = np;
	}
	else{
		auto np = r.decomp_prod_.at(di);
		n.chds_[di] = np;
	}

};

// -------------------------------------------------------------- //
// utils
// -------------------------------------------------------------- //
void PoonSPN::printParams() {

};


