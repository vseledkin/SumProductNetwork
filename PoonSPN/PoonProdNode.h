

#include<vector>
#include <memory>

#include "PoonNode.h"

class PoonProdNode : public PoonNode
{
public:
	std::vector<std::shared_ptr<PoonNode>> chds_;

	void PoonProdNode::passDerivative();

	void eval();

	void addChd(std::shared_ptr<PoonNode> n);
	


};

