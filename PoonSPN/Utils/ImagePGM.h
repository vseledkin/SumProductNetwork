#ifndef IMAGEPGM_H
#define IMAGEPGM_H

#include <vector>
#include <string>

/*
A simple class for writing out images in PGM format.
This is mostly useful if you don't want to rely on an
extra image library like libpng.

This only works for greyscale images right now.
*/
class ImagePGM
{
private:
	unsigned int height, width; //store these so that sub images are easier to work with
	std::vector<std::vector <int> > image;
public:

	//assumes buff is zero meaned
	ImagePGM(const std::vector<std::vector <double> >& buff);

	void writeImage(const std::string pathName) const;
	int getWidth() const { return width; };
	int getHeight() const { return height; };

	//takes a zero meaned image to 128 mean
	//mins and maxes at 0 and 255
	void rescaleImageValues();
};

#endif

