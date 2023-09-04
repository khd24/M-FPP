#pragma once
#include "structures_2.h"
#include "functions_4.h"
#include "BBNode.h"
#include "functions_3.h"
#include "functions_lb.h"
#include "CandidateGenerator.h"

using namespace std;
class CandidateEvaluator
{
public:

	CandidateEvaluator(float lWidth, float lHeight, vector<Barrier> bList, vector<Barrier> nfList, vector<IOPoint> ioList, vector<int> branchingOrder, float* efefFlow, float* efnfFlow, float* nfnfFlow);
	~CandidateEvaluator();

	long long int feasibleCandidates();

	void evaluateUb(ObjectiveFunction* gblOptimal, BBNode& bbNode);

	void evaluateLb(BBNode& bbNode, CandidateGenerator candidateGenerator);

private:

	float lWidth_, lHeight_;
	float* efefFlow_, * efnfFlow_, * nfnfFlow_;
	vector<Barrier> bList_, nfList_;
	vector<IOPoint> ioList_;
	vector<int> branchingOrder_;

	int efCount_, nfCount_;
	long long int feasibleCandidates_;

	void updateFlowMatrices(float* efef2, float* efnf2, float* nfnf2, std::vector<int> placed, std::vector<int> to_be_placed);

	pair<bool, float> computeFiniteLB(vector<Barrier> placed, vector<Barrier> toBePlaced, CandidateGenerator candidateGenerator);
};

