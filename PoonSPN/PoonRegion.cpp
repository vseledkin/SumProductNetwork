
#include "PoonRegion.h"
#include <random>



// 
PoonRegion::PoonRegion(int id, int a1, int a2, int b1, int b2, PoonParameter& params) {
	id_ = id;
	a1_ = a1;	a2_ = a2;
	b1_ = b1;	b2_ = b2;
	a_ = a2_ - a1_; b_ = b2_ - b1_;

	if (a_  >params.baseResolution_ || b_ > params.baseResolution_) {
		if (a_ % params.baseResolution_ != 0 || b_ % params.baseResolution_ != 0) {
			cout << "ERR: base_res= " << params.baseResolution_ << " " << a1 << ", " << a2 << ", " << b1 << ", " << b2;
			exit(-1);
		}
	}
	if (a_ <= params.baseResolution_ && b_ <= params.baseResolution_) 
		interval_ = 1; 
	else 
		interval_ = params.baseResolution_;
}


// initialization
void PoonRegion::resetTypes(int numTypes, PoonParameter& params) {
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

void PoonRegion::setTypes(int numTypes, PoonParameter& params) {
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

void PoonRegion::setBaseForSumOut(PoonParameter& params) {
	defMapTypeIdx_ = -1;
	for (int i = 0; i< params.numComponentsPerVar_; i++) {
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
void PoonRegion::inferMAPForLearning(int instIdx, PoonInstance& inst, PoonParameter& params) {
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
			string di = Decomposition.getIdStr(ri1, ri2, r1.defMapTypeIdx_, r2.defMapTypeIdx_);
			defMapDecompOpts.add(di);
		}
	}
	for (int i = b1_ + interval_; i<b2_; i += interval_) {
		int ri1 = Region.getRegionId(a1_, a2_, b1_, i);
		int ri2 = Region.getRegionId(a1_, a2_, i, b2_);
		Region r1 = Region.getRegion(ri1);
		Region r2 = Region.getRegion(ri2);

		SumNode n1 = r1.types_.get(r1.defMapTypeIdx_);
		SumNode n2 = r2.types_.get(r2.defMapTypeIdx_);
		double lp;
		if (n1.logval_ == Node.ZERO_LOGVAL_ || n2.logval_ == Node.ZERO_LOGVAL_)
			lp = Node.ZERO_LOGVAL_;
		else
			lp = n1.logval_ + n2.logval_;

		if (defMapDecompOpts.isEmpty() || lp>defMapProdPrb_) {
			defMapProdPrb_ = lp;
			defMapDecompOpts.clear();
		}
		if (lp == defMapProdPrb_) {
			String di = Decomposition.getIdStr(ri1, ri2, r1.defMapTypeIdx_, r2.defMapTypeIdx_);
			defMapDecompOpts.add(di);
		}
	}

	// random break ties for a previously unused node
	defMapDecomp = defMapDecompOpts.get(Utils.random_.nextInt(defMapDecompOpts.size()));
	defMapDecompOpts.clear();

	// evaluate product nodes
	for (String di : decomp_prod_.keySet()) {
		ProdNode n = decomp_prod_.get(di);
		n.eval();
	}

	// evaluate existing sum nodes and children
	ArrayList<Integer> mapTypes = new ArrayList<Integer>();
	for (int ti = 0; ti<types_.size(); ti++) {
		if (types_.get(ti).chds_.size() == 0) continue;
		SumNode n = types_.get(ti);
		n.eval();

		double maxSumPrb = 0;
		ArrayList<String> mapDecompOpt = new ArrayList<String>();

		for (String di : n.chds_.keySet()) {
			Node c = n.chds_.get(di);
			double l = n.logval_ + Math.log(n.cnt_);
			double m = c.logval_;
			double nl;

			if (l>m) {
				nl = l + Math.log(1 + Math.exp(m - l));
			}
			else {
				nl = m + Math.log(1 + Math.exp(l - m));
			}

			if (mapDecompOpt.isEmpty() || nl>maxSumPrb) {
				mapDecompOpt.clear();
				maxSumPrb = nl;
			}
			if (nl == maxSumPrb) {
				mapDecompOpt.add(di);
			}
		}

		if (!n.chds_.containsKey(defMapDecomp)) {
			double nl = defMapProdPrb_;
			if (n.logval_ != Node.ZERO_LOGVAL_) {
				nl = Math.log(n.cnt_) + n.logval_;

				if (defMapProdPrb_>nl) {
					nl = defMapProdPrb_ + Math.log(1 + Math.exp(nl - defMapProdPrb_));
				}
				else {
					nl = nl + Math.log(Math.exp(defMapProdPrb_ - nl) + 1);
				}
			}
			nl -= Parameter.sparsePrior_;

			if (mapDecompOpt.isEmpty() || nl>maxSumPrb) {
				mapDecompOpt.clear();
				maxSumPrb = nl;
				mapDecompOpt.add(defMapDecomp);
			}
		}

		n.logval_ = maxSumPrb - Math.log(n.cnt_ + 1);

		// randomly break tie
		mapDecomps_[ti] = mapDecompOpt.get(Utils.random_.nextInt(mapDecompOpt.size()));
		mapDecompOpt.clear();

		if (mapTypes.isEmpty() || n.logval_>defMapSumPrb_) {
			defMapSumPrb_ = n.logval_;
			mapTypes.clear();
		}
		if (n.logval_ == defMapSumPrb_) {
			mapTypes.add(ti);
		}
	}

	if (chosenBlankIdx >= 0) {
		SumNode n = types_.get(chosenBlankIdx);
		n.logval_ = defMapProdPrb_ - Math.log(n.cnt_ + 1) - Parameter.sparsePrior_;
		mapDecomps_[chosenBlankIdx] = defMapDecomp;

		if (mapTypes.isEmpty() || n.logval_>defMapSumPrb_) {
			defMapSumPrb_ = n.logval_;
			mapTypes.clear(); mapTypes.add(chosenBlankIdx);
		}
	}

	defMapTypeIdx_ = mapTypes.get(Utils.random_.nextInt(mapTypes.size()));
	mapTypes.clear();
}

// downward trace-back step
void PoonRegion::setCurrParseToMAP(int instIdx) {
	if (a_ == 1 && b_ == 1) {
		return;
	}

	// type node
	if (types_.size() == 1) inst_type_.put(instIdx, 0);	// only one choice

	int chosenType = inst_type_.get(instIdx);
	String di = mapDecomps_[chosenType];

	inst_decomp_.put(instIdx, di);
	Decomposition d = Decomposition.getDecomposition(di);
	Region r1 = Region.getRegion(d.regionId1_);
	Region r2 = Region.getRegion(d.regionId2_);

	r1.inst_type_.put(instIdx, d.typeId1_);
	r2.inst_type_.put(instIdx, d.typeId2_);

	// record update if slave
	if (!MyMPI.isClassMaster_ && SPN.isRecordingUpdate_) {
		MyMPI.buf_int_[MyMPI.buf_idx_++] = id_;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = chosenType;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = d.regionId1_;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = d.regionId2_;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = d.typeId1_;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = d.typeId2_;
	}

	// if product node not created, create it now
	ProdNode np = decomp_prod_.get(di);
	if (np == null) {
		np = new ProdNode();
		decomp_prod_.put(di, np);
		np.addChd(r1.types_.get(d.typeId1_));
		np.addChd(r2.types_.get(d.typeId2_));
	}

	r1.setCurrParseToMAP(instIdx);
	r2.setCurrParseToMAP(instIdx);
}

// clear an existing parse for incremental EM 
void PoonRegion::clearCurrParse(int instIdx) {
	if (!inst_type_.containsKey(instIdx)) return;
	if (a_ == 1 && b_ == 1) return;
	int cti = inst_type_.get(instIdx);
	String di = inst_decomp_.get(instIdx);

	inst_type_.remove(instIdx);
	inst_decomp_.remove(instIdx);
	Decomposition d = Decomposition.getDecomposition(di);

	// record update if slave
	if (!MyMPI.isClassMaster_ && SPN.isRecordingUpdate_) {
		MyMPI.buf_int_[MyMPI.buf_idx_++] = id_;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = cti;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = d.regionId1_;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = d.regionId2_;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = d.typeId1_;
		MyMPI.buf_int_[MyMPI.buf_idx_++] = d.typeId2_;
	}

	Region r1 = Region.getRegion(d.regionId1_);
	r1.clearCurrParse(instIdx);
	Region r2 = Region.getRegion(d.regionId2_);
	r2.clearCurrParse(instIdx);
}

// clear parse from other slaves
void PoonRegion::clearCurrParseFromBuf(int chosenType, int ri1, int ri2, int ti1, int ti2) {
	if (a_ == 1 && b_ == 1) return;

	String di = Decomposition.getIdStr(ri1, ri2, ti1, ti2);
	SumNode n = types_.get(chosenType);
	n.removeChdOnly(di, 1);
}

void PoonRegion::setCurrParseFromBuf(int chosenType, int ri1, int ri2, int ti1, int ti2) {
	if (a_ == 1 && b_ == 1) return;

	String di = Decomposition.getIdStr(ri1, ri2, ti1, ti2);
	SumNode n = types_.get(chosenType);
	Decomposition d = Decomposition.getDecomposition(di);

	// if prodnode not created, create it now
	ProdNode np = decomp_prod_.get(di);
	if (np == null) {
		np = new ProdNode();
		decomp_prod_.put(di, np);
		Region r1 = Region.getRegion(d.regionId1_);
		Region r2 = Region.getRegion(d.regionId2_);

		np.addChd(r1.types_.get(d.typeId1_));
		np.addChd(r2.types_.get(d.typeId2_));
	}
	n.addChdOnly(di, 1, np);
}
