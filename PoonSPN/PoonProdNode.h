

#include<vector>

#include "PoonNode.h"

class PoonProdNode : public PoonNode
{
public:
	std::vector<PoonNode> chds_;

	PoonProdNode();
	~PoonProdNode();

	void PoonProdNode::passDerivative();

	void eval();

	void addChd(PoonNode& n);
	


};

