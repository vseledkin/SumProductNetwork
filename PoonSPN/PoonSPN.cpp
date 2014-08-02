#include "stdafx.h"
#include "PoonSPN.h"


PoonSPN::PoonSPN(PoonParameter& params){
	this->params = params;
	this->coarseDim1_ = params.inputDim1_ / params.baseResolution_;  //int division, will floor
	this->coarseDim2_ = params.inputDim2_ / params.baseResolution_;
	this->isRecordingUpdate_ = true;
}

PoonSPN::~PoonSPN()
{
}


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
}

void PoonSPN::cmpMAPBottomHalf(PoonInstance& inst) {
	isRecordingUpdate_ = false; // inference now; no need to update count
	int cmpIdx = -1;	// temp idx for cmpMAP				
	inferMAPBottomHalf(cmpIdx, inst);
	setCurrParseToMAP(cmpIdx);
	setMAPBottomToBuf(cmpIdx, inst);
	clearCurrParse(cmpIdx);
	isRecordingUpdate_ = true;
}

// compute marginal by differentiation; see Darwiche-03 for details 
void PoonSPN::cmpMAPBottomHalfMarginal(PoonInstance& inst) {
	setInputOccludeBottomHalf(inst);
	eval();
	cmpDerivative();

	for (int i = 0; i < Parameter.inputDim1_ / 2; i++) {
		for (int j = 0; j < Parameter.inputDim2_; j++)
			MyMPI.buf_int_[MyMPI.buf_idx_++] = Utils.getIntVal(inst, inst.vals_[i][j]);
	}
	for (int i = Parameter.inputDim1_ / 2; i < Parameter.inputDim1_; i++) {
		for (int j = 0; j < Parameter.inputDim2_; j++) {
			int ri = Region.getRegionId(i, i + 1, j, j + 1);
			Region r = Region.getRegion(ri);
			double p = cmpMarginal(r);
			MyMPI.buf_int_[MyMPI.buf_idx_++] = Utils.getIntVal(inst, p);//(int)(p*255);
		}
	}
}


void PoonSPN::inferMAPBottomHalf(int ii, PoonInstance& inst) {
	setInputOccludeBottomHalf(inst);

	// fine region 
	for (int ca = 0; ca < coarseDim1_; ca++)
	for (int cb = 0; cb < coarseDim2_; cb++)
	for (int a = 1; a <= Parameter.baseResolution_; a++)
	for (int b = 1; b <= Parameter.baseResolution_; b++)
	{
		if (a == 1 && b == 1) continue;
		for (int a1 = ca*Parameter.baseResolution_; a1 <= (ca + 1)*Parameter.baseResolution_ - a; a1++) {
			int a2 = a1 + a;
			for (int b1 = cb*Parameter.baseResolution_; b1 <= (cb + 1)*Parameter.baseResolution_ - b; b1++) {
				int b2 = b1 + b;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				r.inferMAP(ii, inst);
			}
		}
	}

	// coarse region
	for (int ca = 1; ca <= coarseDim1_; ca++)
	for (int cb = 1; cb <= coarseDim2_; cb++) {
		if (ca == 1 && cb == 1) continue;	// taken care of below in fine

		for (int a1 = 0; a1 <= Parameter.inputDim1_ - ca*Parameter.baseResolution_; a1 += Parameter.baseResolution_) {
			int a2 = a1 + ca*Parameter.baseResolution_;
			for (int b1 = 0; b1 <= Parameter.inputDim2_ - cb*Parameter.baseResolution_; b1 += Parameter.baseResolution_) {
				int b2 = b1 + cb*Parameter.baseResolution_;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				r.inferMAP(ii, inst);
			}
		}
	}
}

void PoonSPN::setMAPBottomToBuf(int instIdx, PoonInstance& inst) {
	//
	for (int a1 = 0; a1 <= Parameter.inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= Parameter.inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			if (a1 < Parameter.inputDim2_ / 2) {
				MyMPI.buf_int_[MyMPI.buf_idx_++] = Utils.getIntVal(inst, inst.vals_[a1][b1]);
			}
			else {
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				int vi = r.inst_type_.get(instIdx);
				MyMPI.buf_int_[MyMPI.buf_idx_++] = Utils.getIntVal(inst, r.means_[vi]);
			}
		}
	}
}

// ----------------------------------------------------------
// Left
// ----------------------------------------------------------
void PoonSPN::completeLeftImg(PoonInstance& inst) {
	Utils.logTimeMS("before complete left half");
	if (completeByMarginal_) {
		cmpMAPLeftHalfMarginal(inst);
		Utils.logTimeMS("Complete left by Marginal");
	}
	else {
		cmpMAPLeftHalf(inst);
		Utils.logTimeMS("Complete left by MAP");
	}
}

void PoonSPN::cmpMAPLeftHalf(PoonInstance& inst) {
	isRecordingUpdate_ = false;
	int cmpIdx = -1;	// temp idx for cmpMAP
	inferMAPLeftHalf(cmpIdx, inst);
	setCurrParseToMAP(cmpIdx);
	setMAPLeftToBuf(cmpIdx, inst);
	clearCurrParse(cmpIdx);
	isRecordingUpdate_ = true;
}

// compute marginal by differentiation; see Darwiche-03 for details
void PoonSPN::cmpMAPLeftHalfMarginal(PoonInstance& inst) {
	setInputOccludeLeftHalf(inst);
	eval();
	cmpDerivative();

	for (int i = 0; i < Parameter.inputDim1_; i++) {
		for (int j = 0; j < Parameter.inputDim2_ / 2; j++) {
			int ri = Region.getRegionId(i, i + 1, j, j + 1);
			Region r = Region.getRegion(ri);
			double p = cmpMarginal(r);
			MyMPI.buf_int_[MyMPI.buf_idx_++] = Utils.getIntVal(inst, p);
		}
		for (int j = Parameter.inputDim2_ / 2; j < Parameter.inputDim2_; j++)
			MyMPI.buf_int_[MyMPI.buf_idx_++] = Utils.getIntVal(inst, inst.vals_[i][j]);
	}
}

double PoonSPN::cmpMarginal(PoonRegion& r) {
	double t = 0, d = 0;
	double md = 100;

	for (int i = 0; i<r.types_.size(); i++) {
		SumNode n = r.types_.get(i);
		if (n.logDerivative_ == Node.ZERO_LOGVAL_) continue;
		if (md == 100 || n.logDerivative_>md) md = n.logDerivative_;
	}
	for (int i = 0; i < r.types_.size(); i++) {
		SumNode n = r.types_.get(i);
		if (n.logDerivative_ == Node.ZERO_LOGVAL_) continue;
		double p = Math.exp(n.logDerivative_ - md);
		d += r.means_[i] * p; t += p;
	}
	d /= t;
	return d;
}

void PoonSPN::setMAPLeftToBuf(int instIdx, PoonInstance& inst) {
	for (int a1 = 0; a1 <= Parameter.inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= Parameter.inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			if (b1 >= Parameter.inputDim2_ / 2) {
				MyMPI.buf_int_[MyMPI.buf_idx_++] = Utils.getIntVal(inst, inst.vals_[a1][b1]);
			}
			else {
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				int vi = r.inst_type_.get(instIdx);
				MyMPI.buf_int_[MyMPI.buf_idx_++] = Utils.getIntVal(inst, r.means_[vi]);
			}
		}
	}
}

// ----------------------------------------------------------
// Learning
// ----------------------------------------------------------
void PoonSPN::init() {
	// coarse region
	for (int ca = 1; ca <= coarseDim1_; ca++)
	for (int cb = 1; cb <= coarseDim2_; cb++) {
		if (ca == 1 && cb == 1) continue;	// taken care of below in fine

		for (int a1 = 0; a1 <= Parameter.inputDim1_ - ca*Parameter.baseResolution_; a1 += Parameter.baseResolution_) {
			int a2 = a1 + ca*Parameter.baseResolution_;
			for (int b1 = 0; b1 <= Parameter.inputDim2_ - cb*Parameter.baseResolution_; b1 += Parameter.baseResolution_) {
				int b2 = b1 + cb*Parameter.baseResolution_;

				// coarse regions
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				if (ca == coarseDim1_ && cb == coarseDim2_) {
					r.resetTypes(1);	// one sum node as root						
					rootRegion_ = r;
					root_ = r.types_.get(0);
				}
				else r.resetTypes(Parameter.numSumPerRegion_);
			}
		}
	}

	// fine region 
	for (int ca = 0; ca < coarseDim1_; ca++)
	for (int cb = 0; cb < coarseDim2_; cb++)
	for (int a = 1; a <= Parameter.baseResolution_; a++)
	for (int b = 1; b <= Parameter.baseResolution_; b++)
	{
		for (int a1 = ca*Parameter.baseResolution_; a1 <= (ca + 1)*Parameter.baseResolution_ - a; a1++) {
			int a2 = a1 + a;
			for (int b1 = cb*Parameter.baseResolution_; b1 <= (cb + 1)*Parameter.baseResolution_ - b; b1++) {
				int b2 = b1 + b;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				if (a == 1 && b == 1) {
					initUnitRegion(r);
				}
				else r.resetTypes(Parameter.numSumPerRegion_);
			}
		}
	}
}

// init: set mean/variance by equal quantiles from training for each pixel
void PoonSPN::initUnitRegion(PoonRegion& r) {
	r.resetTypes(Parameter.numComponentsPerVar_);

	r.means_ = new double[Parameter.numComponentsPerVar_];
	r.vars_ = new double[Parameter.numComponentsPerVar_];
	r.cnts_ = new double[Parameter.numComponentsPerVar_];

	int ttlCnt = trainingSet_.size();
	int cnt = (int)Math.ceil(ttlCnt*1.0 / Parameter.numComponentsPerVar_);

	double[] vals = new double[ttlCnt];
	for (int ii = 0; ii < trainingSet_.size(); ii++) {
		vals[ii] = trainingSet_.get(ii).vals_[r.a1_][r.b1_];
	}
	Arrays.sort(vals);
	for (int bi = 0; bi < Parameter.numComponentsPerVar_; bi++) {
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
}

void PoonSPN::clearUnusedInSPN() {
	// coarse
	for (int ca = 1; ca <= coarseDim1_; ca++)
	for (int cb = 1; cb <= coarseDim2_; cb++) {
		if (ca == 1 && cb == 1) continue;	// taken care of below in fine
		for (int a1 = 0; a1 <= Parameter.inputDim1_ - ca*Parameter.baseResolution_; a1 += Parameter.baseResolution_) {
			int a2 = a1 + ca*Parameter.baseResolution_;
			for (int b1 = 0; b1 <= Parameter.inputDim2_ - cb*Parameter.baseResolution_; b1 += Parameter.baseResolution_) {
				int b2 = b1 + cb*Parameter.baseResolution_;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);

				Set<String> decomps = new TreeSet<String>();
				for (SumNode n : r.types_) {
					if (n.chds_.size() > 0) {
						double tc = 0;
						for (String ci : n.chdCnts_.keySet()) {
							tc += n.chdCnts_.get(ci);
							decomps.add(ci);
						}
					}
				}

				// clear dead decomp_prod
				Set<String> deadDecomps = new HashSet<String>();
				for (String di : r.decomp_prod_.keySet()) {
					if (!decomps.contains(di)) {
						deadDecomps.add(di);
						continue;
					}
				}

				for (String di : deadDecomps) {
					r.decomp_prod_.remove(di);
					Decomposition.remove(di);
				}
			}
		}
	}

	// fine region
	for (int ca = 0; ca < coarseDim1_; ca++)
	for (int cb = 0; cb < coarseDim2_; cb++)
	for (int a = 1; a <= Parameter.baseResolution_; a++)
	for (int b = 1; b <= Parameter.baseResolution_; b++)
	{
		for (int a1 = ca*Parameter.baseResolution_; a1 <= (ca + 1)*Parameter.baseResolution_ - a; a1++) {
			int a2 = a1 + a;
			for (int b1 = cb*Parameter.baseResolution_; b1 <= (cb + 1)*Parameter.baseResolution_ - b; b1++) {
				int b2 = b1 + b;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);

				// clear dead decomp_prod
				Set<String> decomps = new TreeSet<String>();
				for (SumNode n : r.types_) {
					if (n.chds_.size() > 0) {
						for (String ci : n.chdCnts_.keySet()) {
							decomps.add(ci);
						}
					}
				}
				Set<String> deadDecomps = new HashSet<String>();
				for (String di : r.decomp_prod_.keySet()) {
					if (!decomps.contains(di)) {
						deadDecomps.add(di);
						continue;
					}
				}

				for (String di : deadDecomps) {
					r.decomp_prod_.remove(di);
					Decomposition.remove(di);
				}
			}
		}
	}
}

// ----------------------------------------------------------
// Computation
// ----------------------------------------------------------
// derivative
void PoonSPN::cmpDerivative() {
	initDerivative();

	root_.logDerivative_ = 0;
	root_.passDerivative();
	for (String di : root_.chds_.keySet()) {
		Node n = root_.chds_.get(di);
		n.passDerivative();
	}

	// coarse region
	for (int ca = coarseDim1_; ca >= 1; ca--)
	for (int cb = coarseDim2_; cb >= 1; cb--) {
		if (ca == 1 && cb == 1) continue;	// taken care of below in fine

		for (int a1 = 0; a1 <= Parameter.inputDim1_ - ca*Parameter.baseResolution_; a1 += Parameter.baseResolution_) {
			int a2 = a1 + ca*Parameter.baseResolution_;
			for (int b1 = 0; b1 <= Parameter.inputDim2_ - cb*Parameter.baseResolution_; b1 += Parameter.baseResolution_) {
				int b2 = b1 + cb*Parameter.baseResolution_;

				// coarse regions
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				cmpDerivative(r);
			}
		}
	}


	// fine region 
	for (int ca = coarseDim1_ - 1; ca >= 0; ca--)
	for (int cb = coarseDim2_ - 1; cb >= 0; cb--)
	for (int a = Parameter.baseResolution_; a >= 1; a--)
	for (int b = Parameter.baseResolution_; b >= 1; b--)
	{
		if (a == 1 && b == 1) continue;	// take care in setInput
		for (int a1 = ca*Parameter.baseResolution_; a1 <= (ca + 1)*Parameter.baseResolution_ - a; a1++) {
			int a2 = a1 + a;
			for (int b1 = cb*Parameter.baseResolution_; b1 <= (cb + 1)*Parameter.baseResolution_ - b; b1++) {
				int b2 = b1 + b;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				cmpDerivative(r);
			}
		}
	}
}

void PoonSPN::cmpDerivative(PoonRegion& r) {
	for (int i = 0; i < r.types_.size(); i++) {
		r.types_.get(i).passDerivative();
	}
	for (String di : r.decomp_prod_.keySet()) {
		Node n = r.decomp_prod_.get(di);
		n.passDerivative();
	}
}

void PoonSPN::initDerivative(PoonRegion& r) {
	for (String di : r.decomp_prod_.keySet()) {
		ProdNode n = r.decomp_prod_.get(di);
		n.logDerivative_ = Node.ZERO_LOGVAL_;
	}
	for (SumNode n : r.types_) {
		n.logDerivative_ = Node.ZERO_LOGVAL_;
	}
}

void PoonSPN::initDerivative() {
	for (int ca = coarseDim1_; ca >= 1; ca--)
	for (int cb = coarseDim2_; cb >= 1; cb--) {
		if (ca == 1 && cb == 1) continue;	// taken care of below in fine

		for (int a1 = 0; a1 <= Parameter.inputDim1_ - ca*Parameter.baseResolution_; a1 += Parameter.baseResolution_) {
			int a2 = a1 + ca*Parameter.baseResolution_;
			for (int b1 = 0; b1 <= Parameter.inputDim2_ - cb*Parameter.baseResolution_; b1 += Parameter.baseResolution_) {
				int b2 = b1 + cb*Parameter.baseResolution_;

				// coarse regions
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				initDerivative(r);
			}
		}
	}

	// fine region 
	for (int ca = coarseDim1_ - 1; ca >= 0; ca--)
	for (int cb = coarseDim2_ - 1; cb >= 0; cb--)
	for (int a = Parameter.baseResolution_; a >= 1; a--)
	for (int b = Parameter.baseResolution_; b >= 1; b--)
	{
		for (int a1 = ca*Parameter.baseResolution_; a1 <= (ca + 1)*Parameter.baseResolution_ - a; a1++) {
			int a2 = a1 + a;
			for (int b1 = cb*Parameter.baseResolution_; b1 <= (cb + 1)*Parameter.baseResolution_ - b; b1++) {
				int b2 = b1 + b;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				initDerivative(r);
			}
		}
	}
}


// evaluation: upward pass
void PoonSPN::eval() {
	// fine region 
	for (int ca = 0; ca < coarseDim1_; ca++)
	for (int cb = 0; cb < coarseDim2_; cb++)
	for (int a = 1; a <= Parameter.baseResolution_; a++)
	for (int b = 1; b <= Parameter.baseResolution_; b++)
	{
		if (a == 1 && b == 1) continue;	// take care in setInput
		for (int a1 = ca*Parameter.baseResolution_; a1 <= (ca + 1)*Parameter.baseResolution_ - a; a1++) {
			int a2 = a1 + a;
			for (int b1 = cb*Parameter.baseResolution_; b1 <= (cb + 1)*Parameter.baseResolution_ - b; b1++) {
				int b2 = b1 + b;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				eval(r);
			}
		}
	}

	// coarse region
	for (int ca = 1; ca <= coarseDim1_; ca++)
	for (int cb = 1; cb <= coarseDim2_; cb++) {
		if (ca == 1 && cb == 1) continue;

		for (int a1 = 0; a1 <= Parameter.inputDim1_ - ca*Parameter.baseResolution_; a1 += Parameter.baseResolution_) {
			int a2 = a1 + ca*Parameter.baseResolution_;
			for (int b1 = 0; b1 <= Parameter.inputDim2_ - cb*Parameter.baseResolution_; b1 += Parameter.baseResolution_) {
				int b2 = b1 + cb*Parameter.baseResolution_;

				// coarse regions
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				eval(r);
			}
		}
	}
}

void PoonSPN::eval(PoonRegion& r) {
	for (String di : r.decomp_prod_.keySet()) {
		ProdNode n = r.decomp_prod_.get(di);
		n.eval();
	}
	for (SumNode n : r.types_) {
		if (n.chds_.size() > 0) n.eval();
		else n.logval_ = Node.ZERO_LOGVAL_;
	}
}

// compute MAP
void PoonSPN::inferMAPLeftHalf(int ii, PoonInstance& inst) {
	setInputOccludeLeftHalf(inst);

	// fine region 
	for (int ca = 0; ca < coarseDim1_; ca++)
	for (int cb = 0; cb < coarseDim2_; cb++)
	for (int a = 1; a <= Parameter.baseResolution_; a++)
	for (int b = 1; b <= Parameter.baseResolution_; b++)
	{
		if (a == 1 && b == 1) continue;
		for (int a1 = ca*Parameter.baseResolution_; a1 <= (ca + 1)*Parameter.baseResolution_ - a; a1++) {
			int a2 = a1 + a;
			for (int b1 = cb*Parameter.baseResolution_; b1 <= (cb + 1)*Parameter.baseResolution_ - b; b1++) {
				int b2 = b1 + b;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				r.inferMAP(ii, inst);
			}
		}
	}

	// coarse region
	for (int ca = 1; ca <= coarseDim1_; ca++)
	for (int cb = 1; cb <= coarseDim2_; cb++) {
		if (ca == 1 && cb == 1) continue;	// taken care of below in fine

		for (int a1 = 0; a1 <= Parameter.inputDim1_ - ca*Parameter.baseResolution_; a1 += Parameter.baseResolution_) {
			int a2 = a1 + ca*Parameter.baseResolution_;
			for (int b1 = 0; b1 <= Parameter.inputDim2_ - cb*Parameter.baseResolution_; b1 += Parameter.baseResolution_) {
				int b2 = b1 + cb*Parameter.baseResolution_;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				r.inferMAP(ii, inst);
			}
		}
	}
}

void PoonSPN::inferMAPForLearning(int ii, PoonInstance& inst) {
	setInput(inst);

	// fine region 
	for (int ca = 0; ca < coarseDim1_; ca++)
	for (int cb = 0; cb < coarseDim2_; cb++)
	for (int a = 1; a <= Parameter.baseResolution_; a++)
	for (int b = 1; b <= Parameter.baseResolution_; b++)
	{
		if (a == 1 && b == 1) continue;
		for (int a1 = ca*Parameter.baseResolution_; a1 <= (ca + 1)*Parameter.baseResolution_ - a; a1++) {
			int a2 = a1 + a;
			for (int b1 = cb*Parameter.baseResolution_; b1 <= (cb + 1)*Parameter.baseResolution_ - b; b1++) {
				int b2 = b1 + b;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				r.inferMAPForLearning(ii, inst);
			}
		}
	}

	// coarse region
	for (int ca = 1; ca <= coarseDim1_; ca++)
	for (int cb = 1; cb <= coarseDim2_; cb++) {
		if (ca == 1 && cb == 1) continue;	// taken care of below in fine

		for (int a1 = 0; a1 <= Parameter.inputDim1_ - ca*Parameter.baseResolution_; a1 += Parameter.baseResolution_) {
			int a2 = a1 + ca*Parameter.baseResolution_;
			for (int b1 = 0; b1 <= Parameter.inputDim2_ - cb*Parameter.baseResolution_; b1 += Parameter.baseResolution_) {
				int b2 = b1 + cb*Parameter.baseResolution_;
				int ri = Region.getRegionId(a1, a2, b1, b2);
				Region r = Region.getRegion(ri);
				r.inferMAPForLearning(ii, inst);
			}
		}
	}
}

// clear/set parses
void PoonSPN::clearCurrParse(int ii) {
	rootRegion_.clearCurrParse(ii);
}

void PoonSPN::setCurrParseToMAP(int ii) {
	rootRegion_.setCurrParseToMAP(ii);
}

void PoonSPN::setCurrParseFromBuf() {
	// --- update format: instId, regionId, type, decomp(rid1, rid2, type1, type2) - in buf_inf				
	int k = 0;
	while (k < MyMPI.buf_idx_) {
		int ri = MyMPI.buf_int_[k++];
		int chosenType = MyMPI.buf_int_[k++];
		int ri1 = MyMPI.buf_int_[k++];
		int ri2 = MyMPI.buf_int_[k++];
		int ti1 = MyMPI.buf_int_[k++];
		int ti2 = MyMPI.buf_int_[k++];
		Region r = Region.getRegion(ri);
		r.setCurrParseFromBuf(chosenType, ri1, ri2, ti1, ti2);
	}
}
void PoonSPN::clearCurrParseFromBuf() {
	// --- update format: instId, regionId, type, decomp(rid1, rid2, type1, type2) - in buf_inf				
	int k = 0;
	while (k < MyMPI.buf_idx_) {
		int ri = MyMPI.buf_int_[k++];
		int chosenType = MyMPI.buf_int_[k++];
		int ri1 = MyMPI.buf_int_[k++];
		int ri2 = MyMPI.buf_int_[k++];
		int ti1 = MyMPI.buf_int_[k++];
		int ti2 = MyMPI.buf_int_[k++];
		Region r = Region.getRegion(ri);
		r.clearCurrParseFromBuf(chosenType, ri1, ri2, ti1, ti2);
	}
}

void PoonSPN::sendUpdate(int dest) {
	if (MyMPI.buf_idx_ >= MyMPI.buf_size_) Utils.println("ERR: buffer overflow to " + dest);
	MPI.COMM_WORLD.Send(MyMPI.buf_int_, 0, MyMPI.buf_idx_, MPI.INT, dest, 0);
}

void PoonSPN::recvUpdate(int src) {
	Status status = MPI.COMM_WORLD.Recv(MyMPI.buf_int_, MyMPI.buf_idx_, MyMPI.buf_size_, MPI.INT, src, 0);
	MyMPI.buf_idx_ += status.count;
	if (MyMPI.buf_idx_ >= MyMPI.buf_size_) Utils.println("ERR: buffer overflow from " + src);
}

// compute log probability
double PoonSPN::llh(PoonInstance& inst) {
	setInput(inst);
	eval();
	return root_.getLogVal();
}

// set dspn input
void PoonSPN::setInput(PoonInstance& inst) {
	for (int a1 = 0; a1 <= Parameter.inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= Parameter.inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			int ri = Region.getRegionId(a1, a2, b1, b2);
			Region r = Region.getRegion(ri);
			r.setBase(inst.vals_[a1][b1]);
		}
	}
}

void PoonSPN::setInputOccludeLeftHalf(PoonInstance& inst) {
	for (int a1 = 0; a1 <= Parameter.inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= Parameter.inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			int ri = Region.getRegionId(a1, a2, b1, b2);
			Region r = Region.getRegion(ri);
			if (b1 < Parameter.inputDim2_ / 2) //r.setBase(0,0);	// log 1,1
				r.setBaseForSumOut();
			else
				r.setBase(inst.vals_[a1][b1]);
		}
	}
}

void PoonSPN::setInputOccludeBottomHalf(PoonInstance& inst) {
	for (int a1 = 0; a1 <= Parameter.inputDim1_ - 1; a1++) {
		int a2 = a1 + 1;
		for (int b1 = 0; b1 <= Parameter.inputDim2_ - 1; b1++) {
			int b2 = b1 + 1;
			int ri = Region.getRegionId(a1, a2, b1, b2);
			Region r = Region.getRegion(ri);
			if (a1 >= Parameter.inputDim1_ / 2) //r.setBase(0,0);	// log 1,1
				r.setBaseForSumOut();
			else
				r.setBase(inst.vals_[a1][b1]);
		}
	}
}

// -------------------------------------------------------------- //
// load/save
// -------------------------------------------------------------- //
void PoonSPN::saveDSPN(string mdlName){}
void PoonSPN::loadDSPN(string mdlName, PoonSPN& spn){}



void PoonSPN::addChd(PoonRegion r, PoonSumNode n, string di, double cc) {
	n.setChdCnt(di, cc);
	ProdNode np = r.decomp_prod_.get(di);
	if (np == null) {
		Decomposition d = Decomposition.getDecomposition(di);
		np = new ProdNode();
		r.decomp_prod_.put(di, np);
		Region r1 = Region.getRegion(d.regionId1_);
		Region r2 = Region.getRegion(d.regionId2_);

		np.addChd(r1.types_.get(d.typeId1_));
		np.addChd(r2.types_.get(d.typeId2_));
	}
	n.chds_.put(di, np);
}

// -------------------------------------------------------------- //
// utils
// -------------------------------------------------------------- //
void PoonSPN::printParams() {

}
}


