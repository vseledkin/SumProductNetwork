

#include "PoonDataset.h"
#include <fstream>

using namespace std;

//does normalization and copies buffer to inst
//changed to read in doubles from the start

//For future, normalization should be a option for the data set
//that way you can get the means and stds from the training set and use those
void PoonDataset::setInstance(vector<vector<double> >& buf, PoonInstance& inst){
	double tf = 0;
	double varf = 0;
	int cf = 0;

	for (int i = 0; i < params->inputDim1_; i++){
		for (int j = 0; j < params->inputDim2_; j++) {
			tf += buf[i][j];
			varf += buf[i][j] * buf[i][j]; cf++;
		}
	}
	tf /= cf; 
	varf /= cf; 
	inst.mean_ = tf; 
	inst.std_ = sqrt(varf - tf*tf);

	for (int i = 0; i < params->inputDim1_; i++){
		for (int j = 0; j < params->inputDim2_; j++) {
			inst.vals_[i][j] = (buf[i][j] - inst.mean_) / inst.std_;
		}
	}
}


// --------------------------------------------------------- // 
// Olivetti
// --------------------------------------------------------- // 

void PoonDataset::loadOlivetti(){
	set<int> tis;
	this->genTestIdx(400, params->maxTestSize_, tis);


	ifstream inFile(this->olivettiRawFileName_);

	if (inFile.fail()){
		cout << "Error in opening the olivettiRawFileName_: " << this->olivettiRawFileName_ << endl;
		return;
	}

	vector<vector<double> > faces;

	int imageSize = 4096; //hardcoded like the original
	int numImages = 400;

	faces.resize(numImages, vector<double>(imageSize, 0.0)); //initialize

	//i is the image instace
	//j is the image in a row
	for (int i = 0; i < imageSize; ++i) {
		for (int j = 0; j < numImages; ++j) {
			inFile >> faces[imageSize][numImages];
		}
	}

	inFile.close();

	//make test and training set
	for (int pi = 0; pi<imageSize; pi++) {
		PoonInstance inst;
		readOlivettiInstance(faces, pi, inst);
		if (tis.find(pi) != tis.end()){
			this->test_.push_back(inst);
		}
		else{
			this->train_.push_back(inst);
		}
	}
}

//takes an image as a row and makes it into a 2d matrix
//normalizes image in the set instance funtion
void PoonDataset::readOlivettiInstance(vector<vector<double> >& faces, int pi, PoonInstance& inst){
	vector<vector<double> > buff;
	buff.resize(params->inputDim1_, vector<double>(params->inputDim2_, 0.0)); //initialize

	int k = 0;
	for (int i = 0; i < params->inputDim1_; i++){
		for (int j = 0; j < params->inputDim2_; j++) {
			buff[i][j] = (double)faces[k][pi];
			k++;
		}
	}
	setInstance(buff, inst);
}

