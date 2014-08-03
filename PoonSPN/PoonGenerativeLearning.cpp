#include "PoonGenerativeLearning.h"


using namespace std;





void PoonGenerativeLearning::learnHardEM(vector<PoonInstance> train){

	//long startTime = System.currentTimeMillis();
	bool isLog = false;

	spn_.printParams();
	spn_.init();
	//Utils.logTimeMS("init");

	// -------------------------------------------------- //
	// -- process each mini-batch, find map, update weights
	// -------------------------------------------------- //		
	double ollh = 0;
	double origPrior = params->sparsePrior_;

	int numInstPerSlave = (int)ceil(params->batch_size_*1.0 / params->numSlavePerClass_);

	for (int iter = 1; iter <= params->maxIter_; iter++) {
		// anneal prior
		if (iter <= 10) params->sparsePrior_ = origPrior*iter / 10;

		for (int bi = 0; bi<train.size(); bi += params->batch_size_) {
			// master: aggregate update and pass on
			if (MyMPI.isClassMaster_) {
				// recv clear-parse
				MyMPI.buf_idx_ = 0;
				for (int i = 0; i<params->numSlavePerClass_; i++) {
					if (i*numInstPerSlave<params->batch_size_ && bi + i*numInstPerSlave<train.size()) {
						if (isLog) Utils.println("recv clear update from " + (MyMPI.mySlave_ + i));
						spn_.recvUpdate(MyMPI.mySlave_ + i);
					}
				}

				for (int i = 0; i<params->numSlavePerClass_; i++) {
					if (isLog) Utils.println("send clear update to " + (MyMPI.mySlave_ + i));
					spn_.sendUpdate(MyMPI.mySlave_ + i);
				}

				// recv parse from slaves
				MyMPI.buf_idx_ = 0;
				for (int i = 0; i<params->numSlavePerClass_; i++) {
					if (i*numInstPerSlave<params->batch_size_ && bi + i*numInstPerSlave<train.size()) {
						if (isLog) Utils.println("recv parse update from " + (MyMPI.mySlave_ + i));
						spn_.recvUpdate(MyMPI.mySlave_ + i);
					}
				}
				for (int i = 0; i<params->numSlavePerClass_; i++) {
					if (isLog) Utils.println("send parse update to " + (MyMPI.mySlave_ + i));
					spn_.sendUpdate(MyMPI.mySlave_ + i);
				}
			}
			// slave
			else {
				int k = MyMPI.myOffset_;

				if (k*numInstPerSlave<params->batch_size_ && bi + k*numInstPerSlave<train.size()) {
					MyMPI.buf_idx_ = 0;
					for (int i = k*numInstPerSlave; i<(k + 1)*numInstPerSlave && bi + i<train.size(); i++) {
						// map -> update cnt						
						spn_.clearCurrParse(bi + i);
					}
					if (isLog) Utils.println("send clear update to " + MyMPI.masterRank_);
					spn_.sendUpdate(MyMPI.masterRank_);
				}

				if (isLog) Utils.println("recv clear update from " + MyMPI.masterRank_);

				MyMPI.buf_idx_ = 0;
				spn_.recvUpdate(MyMPI.masterRank_);
				spn_.clearCurrParseFromBuf();

				if (isLog) Utils.logTimeMS("clear parse");

				if (k*numInstPerSlave<params->batch_size_ && bi + k*numInstPerSlave<train.size()) {
					MyMPI.buf_idx_ = 0;
					for (int i = k*numInstPerSlave; i<(k + 1)*numInstPerSlave && bi + i<train.size(); i++) {
						// map -> update cnt						
						spn_.inferMAPForLearning(bi + i, train.get(bi + i));
						spn_.setCurrParseToMAP(bi + i);
					}
					if (isLog) Utils.println("send parse update to " + MyMPI.masterRank_);
					spn_.sendUpdate(MyMPI.masterRank_);
				}

				if (isLog) Utils.logTimeMS("cmp map ...");

				MyMPI.buf_idx_ = 0;
				spn_.recvUpdate(MyMPI.masterRank_);
				spn_.setCurrParseFromBuf();

				if (isLog) Utils.logTimeMS("update weight");
			}
		}
		Utils.logTimeMS("finish iter " + iter);

		if (!MyMPI.isClassMaster_) spn_.clearUnusedInSPN();

		// convergence test
		if (MyMPI.isClassMaster_) {
			double llh = 0;
			for (int i = 0; i<params->numSlavePerClass_; i++) {
				int rank = i + MyMPI.mySlave_;
				llh += recvllh(rank);
			}
			llh /= train.size();
			Utils.logTimeMS("[iter=" + iter + "] llh=" + llh + " ollh=" + ollh);
			if (iter == 1) ollh = llh;
			else {
				double dllh = Math.abs(llh - ollh);
				ollh = llh;
				if (dllh<params->thresholdLLHChg_) {
					Utils.println("\tllh converged");

					// send msg break
					for (int k = 0; k<params->numSlavePerClass_; k++)
						sendMsgBreak(k + MyMPI.mySlave_);
					break;
				}
			}

			// sendMsg OK
			for (int k = 0; k<params->numSlavePerClass_; k++)
				sendMsgOK(k + MyMPI.mySlave_);
		}
		else {
			// slave: compute llh, send to class master
			double llh = 0;
			int size = (int)Math.ceil(train.size()*1.0 / params->numSlavePerClass_);
			for (int i = MyMPI.myOffset_*size; i<(MyMPI.myOffset_ + 1)*size && i<train.size(); i++) {
				Instance inst = train.get(i);
				llh += spn_.llh(inst);
			}
			sendllh(MyMPI.masterRank_, llh);

			// recvMsg: continue?
			char msg = recvMsg(MyMPI.masterRank_);
			if (msg == 'B') break;
		}

		if (isLog) Utils.logTimeMS("done with convergence test");
	}

	// time
	long time = System.currentTimeMillis() - startTime;
	Utils.println("Total learning time: " + time);
}

// ----------------------------------------------------- //
// Utils
// ----------------------------------------------------- //	
void PoonGenerativeLearning::sendMsgBreak(int dest) {
	MyMPI.sendChar(dest, 1, 'B');
}
void PoonGenerativeLearning::sendMsgOK(int dest) {
	MyMPI.sendChar(dest, 1, 'O');
}
char PoonGenerativeLearning::recvMsg(int src) {
	return MyMPI.recvChar(src, 1);
}

void PoonGenerativeLearning::sendllh(int dest, double d) {
	int tag = 1;
	MyMPI.sendDouble(dest, tag, d);
}
double PoonGenerativeLearning::recvllh(int src) {
	int tag = 1;
	return MyMPI.recvDouble(src, tag);
}