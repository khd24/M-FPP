#pragma once

#include "structures_1.h"
#include "structures_2.h"
#include <vector>
#include <algorithm>
#include <set>
#include <map>

using namespace std;

class BBNode
{
public:

	BBNode(BBNode* parent);
	BBNode();
	~BBNode();

	float lb() const;
	void setLb(float lb);

	float ub() const;
	void setUb(float ub);

	vector<Barrier> placedNfs();
	vector<Barrier> partialNfs(CoordinateType coordinateType);
	map<int, float> coordMap(CoordinateType coordinateType);
	
	pair<bool, map<int, float>::iterator> exists(int nfId, CoordinateType coordinateType);

	pair<bool, Barrier> add(Barrier nf, CoordinateType coordinateType);

	bool fathomed() const;

	void setFathomed(bool fathomed);

private:

	float lb_, ub_;

	// nfs that have both x and y coordinates determined
	vector<Barrier> placedNfs_;

	// nfid and order pair
	vector<Barrier> xOnlyNfs_, yOnlyNfs_;

	// mappings of nfids and coordinates
	map<int, float> xCoords_, yCoords_;

	bool fathomed_;
};

class BBNodeComparer {
public:
	inline bool operator()(const BBNode& lhs, const BBNode& rhs) { return float_gt(lhs.lb(), rhs.lb()); };
};
