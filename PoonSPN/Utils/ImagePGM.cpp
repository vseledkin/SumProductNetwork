#include "ImagePGM.h"
#include <iostream>

using namespace std;

ImagePGM::ImagePGM(const vector<vector<double> >& buff){
	if(buff.size() < 1)
		return;

	//copy in image
	this->image.resize(buff.size(), vector<int>(buff[0].size(), 0)); //initialize

	for (unsigned int i = 0; i < buff.size(); i++) {
		for (unsigned int j = 0; j < buff[0].size(); j++) {
			//cout << buff[i][j] << endl;
			image[i][j] = (int)(buff[i][j] * 127.0);
		}
	}
	this->height = buff.size();
	this->width = buff[0].size();
};


void ImagePGM::writeImage(const string pathName) const{


	FILE * grey = fopen(pathName.c_str(), "w");

	fprintf(grey, "P2\n%d %d\n255\n", width, height);

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			fprintf(grey, "%d ", image[i][j]);
		}
		fprintf(grey, "\n");
	}

	fclose(grey);

};


void ImagePGM::rescaleImageValues(){
	for (int h = 0; h < height; h++){
		for (int w = 0; w < width; w++){
			image[h][w] += 127;
			if (image[h][w] < 0) 
				image[h][w] = 0;
			if (image[h][w] > 225)
				image[h][w] = 255;
		}
	}
};

