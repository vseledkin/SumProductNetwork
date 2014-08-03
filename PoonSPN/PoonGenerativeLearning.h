#include <vector>

#include "PoonSPN.h"
#include "PoonInstance.h"
#include "PoonParameter.h"



class PoonGenerativeLearning
{
public:
	PoonSPN spn_;

	std::shared_ptr<PoonParameter> params;

	PoonSPN& getDSPN() { return spn_; };


	PoonGenerativeLearning(){};

	void learn(std::vector<PoonInstance> train, std::shared_ptr<PoonParameter> params){
		this->params = params;
		spn_ = PoonSPN(train, params);

		learnHardEM(train);
	}
	//void saveModel(string mdlFileName){spn_.saveDSPN(mdlFileName);}
	void learnHardEM(std::vector<PoonInstance> train);

	//MPI functions
	void sendMsgBreak(int dest);
	void sendMsgOK(int dest);
	char recvMsg(int src);
	void sendllh(int dest, double d);
	double recvllh(int src);
};

