#include "BBNode.h"

BBNode::BBNode(BBNode* parent) : placedNfs_(parent->placedNfs()),
xOnlyNfs_(parent->partialNfs(CoordinateType::X)),
yOnlyNfs_(parent->partialNfs(CoordinateType::Y)),
xCoords_(parent->coordMap(CoordinateType::X)),
yCoords_(parent->coordMap(CoordinateType::Y)) {

	lb_ = parent->lb();
	ub_ = numeric_limits<float>::max();
	fathomed_ = false;
}

BBNode::BBNode() {
	lb_ = 0;
	ub_ = numeric_limits<float>::max();
	fathomed_ = false;
}

BBNode::~BBNode()
{
	
}

vector<Barrier> BBNode::placedNfs() {
	return placedNfs_;
}

bool BBNode::fathomed() const {
	return fathomed_;
}

void BBNode::setFathomed(bool fathomed) {
	fathomed_ = fathomed;
}

vector<Barrier> BBNode::partialNfs(CoordinateType coordinateType)
{
	switch (coordinateType)
	{
	case CoordinateType::X:
		return xOnlyNfs_;
	case CoordinateType::Y:
		return yOnlyNfs_;
	}
}

map<int, float> BBNode::coordMap(CoordinateType coordinateType) {
	
	switch (coordinateType)
	{
	case CoordinateType::X:
		return xCoords_;
	case CoordinateType::Y:
		return yCoords_;
	}
}

float BBNode::lb() const
{
	return lb_;
}

void BBNode::setLb(float lb) {
	lb_ = lb;
}

float BBNode::ub() const
{
	return ub_;
}

void BBNode::setUb(float ub) {
	ub_ = ub;
}

pair<bool, map<int, float>::iterator> BBNode::exists(int nfId, CoordinateType coordinateType) {
	switch (coordinateType)
	{
	case CoordinateType::X:
	{
		auto it = xCoords_.find(nfId);
		return  pair<bool, map<int, float>::iterator>(it != xCoords_.end(), it);
	}
	case CoordinateType::Y:
	{
		auto it = yCoords_.find(nfId);
		return  pair<bool, map<int, float>::iterator>(it != yCoords_.end(), it);
	}
	}
}

pair<bool, Barrier> BBNode::add(Barrier nf, CoordinateType coordinateType) {

	Barrier retNf{};
	
	bool success = false;
	int nfId = nf.id();
	float xCoord = nf.getMinX();
	float yCoord = nf.getMaxY();

	switch (coordinateType)
	{
	case CoordinateType::X:
	{
		xOnlyNfs_.push_back(nf);
		xCoords_.insert(pair<int, float>(nfId, xCoord));
		auto yExists = exists(nfId, CoordinateType::Y);
		success = yExists.first;
		yCoord = success ? (yExists.second)->second : yCoord;
		break;
	}

	case CoordinateType::Y:
		yOnlyNfs_.push_back(nf);
		yCoords_.insert(pair<int, float>(nfId, yCoord));
		auto xExists = exists(nfId, CoordinateType::X);
		success = xExists.first;
		xCoord = success ? (xExists.second)->second : xCoord;
		break;
	}

	// If a complete NF coordinate was found, success becomes true
	if (success) {
		retNf.setBarrier(xCoord, yCoord, nf.getWidth(), nf.getHeight(), nf.getCongestion());
		retNf.setId(nfId);
		placedNfs_.push_back(retNf);
	}

	return pair<bool, Barrier>(success, retNf);
}
