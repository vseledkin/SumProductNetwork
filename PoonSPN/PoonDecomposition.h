#ifndef POONDECOMPOSITION_H
#define POONDECOMPOSITION_H

#include <map>
#include <string>
#include <sstream>

class PoonDecomposition
{
public:

	static std::map<std::string, PoonDecomposition> id_decomp_;  // .cpp has instantiation of id_decomp_

	std::string id_;
	int regionId1_, regionId2_, typeId1_, typeId2_;

	PoonDecomposition(std::string id, int regionId1, int regionId2, int typeId1, int typeId2) {
		id_ = id;
		regionId1_ = regionId1;	
		regionId2_ = regionId2;
		typeId1_ = typeId1;	
		typeId2_ = typeId2;
	}

	PoonDecomposition& PoonDecomposition::operator=(PoonDecomposition&& other){ //c++11 move operator
		//id_decomp_ = other.id_decomp_;
		id_ = other.id_;
		regionId1_ = other.regionId1_;
		regionId2_ = other.regionId2_;
		typeId1_ = other.typeId1_;
		typeId2_ = other.typeId2_;
	}

	//unused function?
	/*
	static PoonDecomposition& getDecomposition(int regionId1, int regionId2, int typeId1, int typeId2) {
		std::string id = getIdStr(regionId1, regionId2, typeId1, typeId2);
		
		if (id_decomp_.find(id) == id_decomp_.end()){
			//id_decomp_[id] = d; //java code adds null to the map? why
		}else
			return  id_decomp_.at(id);
	}*/

	static PoonDecomposition& getDecomposition(std::string id) {
		if (id_decomp_.find(id) == id_decomp_.end()) {
			if (id.empty()) {
				// blank; doesn't matter, blank is no longer static
				return PoonDecomposition("", -1, -1, -1, -1);
			}
			//String[] ts = id.split(" ");
			std::stringstream stream(id);

			int regionId1, regionId2, typeId1, typeId2;
			stream >> regionId1;
			stream >> regionId2;
			stream >> typeId1;
			stream >> typeId2;
			return PoonDecomposition(id, regionId1, regionId2, typeId1, typeId2);
		}
		else{
			return id_decomp_.at(id);
		}
	};

	std::string getId() { return id_; };

	static void remove(std::string id) {
		id_decomp_.erase(id);
	}
	static std::string getIdStr(int regionId1, int regionId2, int typeId1, int typeId2) {
		std::stringstream id;
		id << regionId1 << " " << regionId2 << " " << typeId1 << " " << typeId2;
		return id.str();
	}
};

#endif