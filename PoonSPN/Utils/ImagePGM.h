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
	ImagePGM(std::vector<std::vector <double> >& buff);

	void writeImage(std::string pathName);
	int getWidth(){ return width;};
	int getHeight(){ return height; };
};

#endif

