#include "PoonParameter.h"

std::vector<int> PoonParameter::buf_int_(PoonParameter::buf_size_, 0);  //instantiation of ststic members
std::vector<double> PoonParameter::buf_double_(PoonParameter::buf_size_double_, 0.0);