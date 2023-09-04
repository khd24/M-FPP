#pragma once

#include <vector>
#include <queue>
#include <stack>
#include <set>
#include <omp.h>
#include "functions_1.h"
#include "functions_4.h"
#include "BBNode.h"
#include "structures_1.h"
#include "structures_2.h"
#include "BBQueue.h"
#include "functions_lb.h"
#include "CandidateGenerator.h"
#include <map>
#include <set>
#include "CandidateEvaluator.h"
#include <iomanip>

using namespace std;

class BBController
{
public:

	BBController(float lWidth, float lHeight, vector<Barrier> bList, vector<Barrier> nfList, vector<IOPoint> ioList, vector<int> branchingOrder, float* efefFlow, float* efnfFlow, float* nfnfFlow);
	~BBController();

	void run(ObjectiveFunction* gblOptimal);

	long long int processed();
	long long int remaining();
	float bestUb();
	float bestLb();

	void setMaxGap(float maxGap);
	void setMaxTime(int maxTime);
	void setBranchingStrategy(BranchingStrategy strategy);

private:
	BBQueue bbQueue_;
	CandidateGenerator candidateGenerator_;
	CandidateEvaluator candidateEvaluator_;

	int nfCount_;
	float lWidth_, lHeight_;
	float bestLb_, bestUb_;

	double startTime_;

	float maxGap_;
	int maxTime_;
	bool enableLb_;

	vector<int> branchingOrder_;

	bool branch(ObjectiveFunction* gblOptimal, BBNode bbNode, CoordinateType coordType);

	void bound(ObjectiveFunction* gblOptimal, BBNode& bbNode);

	void rootInitialize(ObjectiveFunction* gblOptimal);

	void progress(bool header);

	bool stop();

	float gap(float lb);

	const char separator = ' ';
	const int numWidth = 16;
	template<typename T> void printElement(T t, const int& width);
};

