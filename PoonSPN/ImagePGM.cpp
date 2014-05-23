#include "ImagePGM.h"


ImagePGM::ImagePGM(vector<vector <double> >& buff){
	if(buff.size() < 1)
		return;

	//copy in image
	this->image.resize(buff.size(), vector<int>(buff[0].size(), 0.0)); //initialize

	for (unsigned int i = 0; i < buff.size(); i++) {
		for (unsigned int j = 0; j < buff[0].size(); j++) {
			image[i][j] = (int)buff[i][j];
		}
	}
	this->height = buff.size();
	this->width = buff[0].size();
}


void ImagePGM::writeImage(string pathName){

	FILE * grey = fopen(pathName.c_str(), "w");

	fprintf(grey, "P2\n%d %d\n255\n", width, height);

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			fprintf(grey, "%d ", image[i][j]);
		}
		fprintf(grey, "\n");
	}

	fclose(grey);

}
