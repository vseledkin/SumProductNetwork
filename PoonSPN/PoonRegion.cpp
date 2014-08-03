
#include "PoonRegion.h"
#include "PoonDecomposition.h"
#include "PoonSPN.h"

#include <random>

using namespace std;

// 
PoonRegion::PoonRegion(int id, int a1, int a2, int b1, int b2, std::shared_ptr<PoonParameter> params) {
	id_ = id;
	a1_ = a1;	a2_ = a2;
	b1_ = b1;	b2_ = b2;
	a_ = a2_ - a1_; b_ = b2_ - b1_;

	if (a_  > params->baseResolution_ || b_ > params->baseResolution_) {
		if (a_ % params->baseResolution_ != 0 || b_ % params->baseResolution_ != 0) {
			cout << "ERR: base_res= " << params->baseResolution_ << " " << a1 << ", " << a2 << ", " << b1 << ", " << b2;
			exit(-1);
		}
	}
	if (a_ <= params->baseResolution_ && b_ <= params->baseResolution_) 
		interval_ = 1; 
	else 
		interval_ = params->baseResolution_;
}


// initialization
void PoonRegion::resetTypes(int numTypes, std::shared_ptr<PoonParameter> params) {
	// clean up
	types_.clear();
	inst_type_.clear();
	inst_decomp_.clear();
	decomp_prod_.clear();
	mapDecomps_.clear();
	for (int i = 0; i< numTypes; i++) {
		types_.push_back(PoonSumNode(params));
	}
}

void PoonRegion::setTypes(int numTypes, std::shared_ptr<PoonParameter> params) {
	if (numTypes<types_.size()) { 
		resetTypes(numTypes, params); 
		return; 
	}
	int nn = numTypes - types_.size();
	for (int i = 0; i<nn; i++) {
		types_.push_back(PoonSumNode(params));
	}
}


void PoonRegion::setBaseGauss(double v) {
	defMapTypeIdx_ = -1;
	double mp = 0;
	for (int i = 0; i < types_.size(); i++) {
		PoonSumNode& n = types_.at(i);
		n.setVal(cmpGauss(v, means_[i]));
		if (defMapTypeIdx_ == -1 || n.getLogVal() >mp) {
			defMapTypeIdx_ = i;
			mp = n.getLogVal();
		}
	}
}

void PoonRegion::setBaseForSumOut(std::shared_ptr<PoonParameter> params) {
	defMapTypeIdx_ = -1;
	for (int i = 0; i< params->numComponentsPerVar_; i++) {
		PoonSumNode& n = types_.at(i);
		n.setVal(0);
	}
}

// compute MAP state at inference time
void PoonRegion::inferMAP(int instIdx, PoonInstance& inst) {

	// compute prod values
	for (auto& kv : decomp_prod_) {
		PoonProdNode& n = decomp_prod_.at(kv.first);
		n.eval();
	}

	// evaluate children for sum nodes 
	for (int ti = 0; ti < types_.size(); ti++) {
		if (types_.at(ti).chds_.size() == 0) 
			continue;
		PoonSumNode& n = types_.at(ti);
		n.eval();

		double maxChdPrb = 0;
		vector<string> mapDecompOpt;
		for (auto& kv : n.chds_) {
			auto di = kv.first;
			PoonNode& c = n.chds_.at(di);
			double m = (c.getLogVal() == c.ZERO_LOGVAL_) ? c.ZERO_LOGVAL_ : c.getLogVal + log(n.getChdCnt(di));

			if (mapDecompOpt.empty() || m > maxChdPrb) {
				mapDecompOpt.clear();
				maxChdPrb = m;
			}
			if (m == maxChdPrb) {
				mapDecompOpt.push_back(di);
			}
		}

		// randomly break tie if there is more than one maximum value
		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(0, mapDecompOpt.size());

		mapDecomps_[ti] = mapDecompOpt.at(distribution(generator));
		mapDecompOpt.clear();
	}
}

// compute MAP state at learning time: could tap a previous unused node
void PoonRegion::inferMAPForLearning(int instIdx, PoonInstance& inst, std::shared_ptr<PoonParameter> params) {
	vector<string> defMapDecompOpts;
	string defMapDecomp;

	defMapTypeIdx_ = -1;
	defMapSumPrb_ = 100;
	defMapProdPrb_ = 100;
	defMapDecompOpts.clear();

	// sum: choose a previous unused node
	vector<int> blanks;
	for (int i = 0; i<types_.size(); i++) {
		PoonSumNode& n = types_.at(i);
		if (n.chds_.size() == 0) {
			blanks.push_back(i);
		}
	}
	int chosenBlankIdx = -1;
	if (!blanks.empty()) { //empty or just a safty check?
		if (blanks.size()>1) {

			std::default_random_engine generator;
			std::uniform_int_distribution<int> distribution(0, blanks.size());

			int ci = distribution(generator);
			chosenBlankIdx = blanks.at(ci);
		}
		else 
			chosenBlankIdx = blanks.at(0);
		blanks.clear();
	}

	// find MAP decomposition
	for (int i = a1_ + interval_; i<a2_; i += interval_) {
		int ri1 = getRegionId(a1_, i, b1_, b2_, params);
		int ri2 = getRegionId(i, a2_, b1_, b2_, params);
		PoonRegion r1 = getRegion(ri1, params);
		PoonRegion r2 = getRegion(ri2, params);
		PoonSumNode n1 = r1.types_.at(r1.defMapTypeIdx_);
		PoonSumNode n2 = r2.types_.at(r2.defMapTypeIdx_);
		double lp;

		if (n1.getLogVal() == n1.ZERO_LOGVAL_ || n2.getLogVal() == n2.ZERO_LOGVAL_)
			lp = n1.ZERO_LOGVAL_;
		else
			lp = n1.getLogVal() + n2.getLogVal();

		if (defMapDecompOpts.empty() || lp>defMapProdPrb_) {
			defMapProdPrb_ = lp;
			defMapDecompOpts.clear();
		}
		if (lp == defMapProdPrb_) {
			string di = PoonDecomposition::getIdStr(ri1, ri2, r1.defMapTypeIdx_, r2.defMapTypeIdx_);
			defMapDecompOpts.push_back(di);
		}
	}
	for (int i = b1_ + interval_; i<b2_; i += interval_) {
		int ri1 = getRegionId(a1_, a2_, b1_, i, params);
		int ri2 = getRegionId(a1_, a2_, i, b2_, params);
		PoonRegion r1 = getRegion(ri1, params);
		PoonRegion r2 = getRegion(ri2, params);

		PoonSumNode n1 = r1.types_.at(r1.defMapTypeIdx_);
		PoonSumNode n2 = r2.types_.at(r2.defMapTypeIdx_);
		double lp;
		if (n1.getLogVal() == n1.ZERO_LOGVAL_ || n2.getLogVal() == n2.ZERO_LOGVAL_)
			lp = n1.ZERO_LOGVAL_;
		else
			lp = n1.getLogVal() + n2.getLogVal();

		if (defMapDecompOpts.empty() || lp > defMapProdPrb_) {
			defMapProdPrb_ = lp;
			defMapDecompOpts.clear();
		}
		if (lp == defMapProdPrb_) {
			string di = PoonDecomposition::getIdStr(ri1, ri2, r1.defMapTypeIdx_, r2.defMapTypeIdx_);
			defMapDecompOpts.push_back(di);
		}
	}

	// random break ties for a previously unused node
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, blanks.size());

	defMapDecomp = defMapDecompOpts.at(distribution(generator));
	defMapDecompOpts.clear();

	// evaluate product nodes
	for (auto& kv : decomp_prod_) {
		PoonProdNode n = decomp_prod_.at(kv.first);  //could just kv.second.eval();
		n.eval();
	}

	// evaluate existing sum nodes and children
	vector<int> mapTypes;
	for (int ti = 0; ti < types_.size(); ti++) {
		if (types_.at(ti).chds_.size() == 0) 
			continue;
		PoonSumNode n = types_.at(ti);
		n.eval();

		double maxSumPrb = 0;
		vector<string> mapDecompOpt;

		for (auto& kv : n.chds_) {
			string di = kv.first;
			PoonNode& c = n.chds_.at(di);
			double l = n.getLogVal() + log(n.cnt_);
			double m = c.getLogVal();
			double nl;

			if (l>m) {
				nl = l + log(1 + exp(m - l));
			}
			else {
				nl = m + log(1 + exp(l - m));
			}

			if (mapDecompOpt.empty() || nl>maxSumPrb) {
				mapDecompOpt.clear();
				maxSumPrb = nl;
			}
			if (nl == maxSumPrb) {
				mapDecompOpt.push_back(di);
			}
		}

		if (n.chds_.find(defMapDecomp) == n.chds_.end()) {
			double nl = defMapProdPrb_;
			if (n.getLogVal() != n.ZERO_LOGVAL_) {
				nl = log(n.cnt_) + n.getLogVal();

				if (defMapProdPrb_>nl) {
					nl = defMapProdPrb_ + log(1 + exp(nl - defMapProdPrb_));
				}
				else {
					nl = nl + log(exp(defMapProdPrb_ - nl) + 1);
				}
			}

			nl -= params->sparsePrior_;

			if (mapDecompOpt.empty() || nl>maxSumPrb) {
				mapDecompOpt.clear();
				maxSumPrb = nl;
				mapDecompOpt.push_back(defMapDecomp);
			}
		}

		n.setVal(maxSumPrb - log(n.cnt_ + 1));

		// randomly break tie
		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(0, mapDecompOpt.size());


		mapDecomps_[ti] = mapDecompOpt.at(distribution(generator));
		mapDecompOpt.clear();

		if (mapTypes.empty() || n.getLogVal() > defMapSumPrb_) {
			defMapSumPrb_ = n.getLogVal();
			mapTypes.clear();
		}
		if (n.getLogVal() == defMapSumPrb_) {
			mapTypes.push_back(ti);
		}
	}

	if (chosenBlankIdx >= 0) {
		PoonSumNode& n = types_.at(chosenBlankIdx);
		n.setVal(defMapProdPrb_ - log(n.cnt_ + 1) - params->sparsePrior_);
		mapDecomps_[chosenBlankIdx] = defMapDecomp;

		if (mapTypes.empty() || n.getLogVal() > defMapSumPrb_) {
			defMapSumPrb_ = n.getLogVal();
			mapTypes.clear(); mapTypes.push_back(chosenBlankIdx);
		}
	}

	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, mapTypes.size());

	defMapTypeIdx_ = mapTypes.at(distribution(generator));
	mapTypes.clear();
}

// downward trace-back step
void PoonRegion::setCurrParseToMAP(int instIdx, std::shared_ptr<PoonParameter> params) {
	if (a_ == 1 && b_ == 1) {
		return;
	}

	// type node
	if (types_.size() == 1) 
		inst_type_[instIdx] = 0;	// only one choice

	int chosenType = inst_type_.at(instIdx);
	string di = mapDecomps_[chosenType];

	inst_decomp_[instIdx]  = di;
	PoonDecomposition d = PoonDecomposition::getDecomposition(di);
	PoonRegion r1 = PoonRegion::getRegion(d.regionId1_, params);
	PoonRegion r2 = PoonRegion::getRegion(d.regionId2_, params);

	r1.inst_type_[instIdx] = d.typeId1_;
	r2.inst_type_[instIdx] = d.typeId2_;

	// record update if slave
	if (PoonSPN::isRecordingUpdate_) {
		params->buf_int_[params->buf_idx_++] = id_;
		params->buf_int_[params->buf_idx_++] = chosenType;
		params->buf_int_[params->buf_idx_++] = d.regionId1_;
		params->buf_int_[params->buf_idx_++] = d.regionId2_;
		params->buf_int_[params->buf_idx_++] = d.typeId1_;
		params->buf_int_[params->buf_idx_++] = d.typeId2_;
	}

	// if product node not created, create it now
	if (decomp_prod_.find(di) == decomp_prod_.end()) {
		PoonProdNode np;
		decomp_prod_[di] = np;
		np.addChd(r1.types_.at(d.typeId1_));
		np.addChd(r2.types_.at(d.typeId2_));
	}

	r1.setCurrParseToMAP(instIdx, params);
	r2.setCurrParseToMAP(instIdx, params);
}

// clear an existing parse for incremental EM 
void PoonRegion::clearCurrParse(int instIdx, std::shared_ptr<PoonParameter> params) {
	if (inst_type_.find(instIdx) == inst_type_.end())
		return;
	if (a_ == 1 && b_ == 1) 
		return;
	int cti = inst_type_.at(instIdx);
	string di = inst_decomp_.at(instIdx);

	inst_type_.erase(instIdx);
	inst_decomp_.erase(instIdx);
	PoonDecomposition d = PoonDecomposition::getDecomposition(di);

	// record update if slave
	if (PoonSPN::isRecordingUpdate_) {  //!MyMPI.isClassMaster_
		params->buf_int_[params->buf_idx_++] = id_;
		params->buf_int_[params->buf_idx_++] = cti;
		params->buf_int_[params->buf_idx_++] = d.regionId1_;
		params->buf_int_[params->buf_idx_++] = d.regionId2_;
		params->buf_int_[params->buf_idx_++] = d.typeId1_;
		params->buf_int_[params->buf_idx_++] = d.typeId2_;
	}

	PoonRegion r1 = PoonRegion::getRegion(d.regionId1_, params);
	r1.clearCurrParse(instIdx, params);
	PoonRegion r2 = PoonRegion::getRegion(d.regionId2_, params);
	r2.clearCurrParse(instIdx, params);
}

// clear parse from other slaves
void PoonRegion::clearCurrParseFromBuf(int chosenType, int ri1, int ri2, int ti1, int ti2, std::shared_ptr<PoonParameter> params) {
	if (a_ == 1 && b_ == 1) 
		return;

	string di = PoonDecomposition::getIdStr(ri1, ri2, ti1, ti2);
	PoonSumNode n = types_.at(chosenType);
	n.removeChdOnly(di, 1);
}

void PoonRegion::setCurrParseFromBuf(int chosenType, int ri1, int ri2, int ti1, int ti2, std::shared_ptr<PoonParameter> params) {
	if (a_ == 1 && b_ == 1) 
		return;

	string di = PoonDecomposition::getIdStr(ri1, ri2, ti1, ti2);
	PoonSumNode& n = types_.at(chosenType);
	PoonDecomposition& d = PoonDecomposition::getDecomposition(di);

	// if prodnode not created, create it now
	if (decomp_prod_.find(di) == decomp_prod_.end()) {
		PoonProdNode np;
		decomp_prod_[di] =  np;
		PoonRegion& r1 = PoonRegion::getRegion(d.regionId1_, params);
		PoonRegion& r2 = PoonRegion::getRegion(d.regionId2_, params);

		np.addChd(r1.types_.at(d.typeId1_));
		np.addChd(r2.types_.at(d.typeId2_));
		n.addChdOnly(di, 1, np);
	}
	else{
		PoonProdNode& np = decomp_prod_.at(di);
		n.addChdOnly(di, 1, np);
	}
}
