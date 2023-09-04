#pragma once

#include<vector>
#include "BBNode.h"
#include "functions_1.h"

using namespace std;
class CandidateGenerator
{
public:

	CandidateGenerator(vector<Barrier> bList, vector<Barrier> nfList, vector<IOPoint> ioList, float lWidth, float lHeight);
	~CandidateGenerator();

	vector<float> getCandidates(CoordinateType coordinateType, vector<Barrier> placedNfs);

	vector<float> getCandidatesPolytomic(CoordinateType coordinateType, vector<Barrier> nfs, int nfOrder);

	vector<float> getPotentialCoordinates(CoordinateType coordinateType, Barrier nf, float coordinate);

	vector<Barrier> getPotentialPolytomicCoordinates(CoordinateType coordinateType, vector<Barrier> nfs, int newNfId, int nfOrder);

	bool isFeasible(vector<Barrier> nfs, Barrier nf);

private:
	float lWidth_, lHeight_;
	vector<float> xGridlines_, yGridlines_;
	vector<Barrier> bList_, nfList_;
};

