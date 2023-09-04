#include "BBController.h"

BBController::BBController(float lWidth, float lHeight, vector<Barrier> bList, vector<Barrier> nfList, vector<IOPoint> ioList, vector<int> branchingOrder, float* efefFlow, float* efnfFlow, float* nfnfFlow) 
	: candidateGenerator_(bList, nfList, ioList, lWidth, lHeight),
	candidateEvaluator_(lWidth, lHeight, bList, nfList, ioList, branchingOrder, efefFlow, efnfFlow, nfnfFlow) {
	
	lWidth_ = lWidth;
	lHeight_ = lHeight;
	branchingOrder_ = branchingOrder;

	bestLb_ = 0.0f;
	bestUb_ = numeric_limits<float>::max();

	nfCount_ = nfList.size();

	maxGap_ = -1.0;
	maxTime_ = numeric_limits<int>::max();
	enableLb_ = false;

	startTime_ = omp_get_wtime();
}

BBController::~BBController() {

}

void BBController::run(ObjectiveFunction* gblOptimal) {
	
	progress(true);

	rootInitialize(gblOptimal);
	
	while (!bbQueue_.empty()) {		
		BBNode node = bbQueue_.pop();

		if (!node.fathomed() && node.placedNfs().size() < nfCount_ && node.lb() < bestUb_ && gap(node.lb()) > maxGap_) {
			
			bestLb_ = node.lb();

			auto coordType = node.partialNfs(CoordinateType::X).size() != node.partialNfs(CoordinateType::Y).size()
				? CoordinateType::X
				: CoordinateType::Y;
			
			if (branch(gblOptimal, node, coordType)) {
				break;
			}
		}

		progress(false);
	}

	gblOptimal->candidatecount = processed();
	gblOptimal->bestbound = bestLb_;
}

long long int BBController::processed() {
	return bbQueue_.processed();
}

long long int  BBController::remaining() {
	return bbQueue_.size();
}

float BBController::bestUb() {
	return bestUb_;
}

float BBController::bestLb() {
	return bestLb_;
}

void BBController::setMaxGap(float maxGap)
{
	enableLb_ = !float_lt(maxGap, 0);

	if (float_lt(maxGap, 0)) maxGap_ = 0.0001;
	else if (float_gt(maxGap, 1)) maxGap_ = 1.0;
	else maxGap_ = maxGap;
}

void BBController::setMaxTime(int maxTime)
{
	if (maxTime > 0) maxTime_ = maxTime; // at least 1 second.
}

void BBController::setBranchingStrategy(BranchingStrategy strategy)
{
	bbQueue_.setStrategy(strategy);
}

bool BBController::branch(ObjectiveFunction* gblOptimal, BBNode parentNode, CoordinateType coordType) {

	for (int i = 0; i < nfCount_; i++) {

		int nfId = branchingOrder_[i];

		if (parentNode.exists(nfId, coordType).first) continue;

		auto nfs = candidateGenerator_.getPotentialPolytomicCoordinates(coordType, parentNode.partialNfs(coordType), nfId, i);

#pragma omp parallel for
		for (int i = 0; i < nfs.size(); i++) {

			BBNode newNode(parentNode);
			auto addCoordinate = newNode.add(nfs[i], coordType);

			/*if (stop()) {
				return true;
			}*/

			bool fathomed = addCoordinate.first && !candidateGenerator_.isFeasible(parentNode.placedNfs(), addCoordinate.second);

			newNode.setFathomed(fathomed);

			if (!fathomed) {
				bound(gblOptimal, newNode);
			}

#pragma omp critical
			bbQueue_.push(newNode);
		}
	}

	return stop();
}

void BBController::bound(ObjectiveFunction* gblOptimal, BBNode& bbNode) {
	
	candidateEvaluator_.evaluateUb(gblOptimal, bbNode);
	
	if (enableLb_) {
		candidateEvaluator_.evaluateLb(bbNode, candidateGenerator_);
	}

	if (bbNode.ub() < bestUb_) {
		bestUb_ = gblOptimal->value;
	}
}

void BBController::rootInitialize(ObjectiveFunction* gblOptimal) {
	BBNode root;	
	bound(gblOptimal, root);
	bbQueue_.push(root);
}

void BBController::progress(bool header) {
	if (header) {
		
		printElement("NNodes", numWidth);
		printElement("NRem", numWidth);
		printElement("NFeas", numWidth);
		printElement("ObjVal", numWidth);
		printElement("BestBound", numWidth);
		printElement("ETime(s)", numWidth);
		cout << endl;

		return;
	}

	auto curTime = omp_get_wtime();
	
	if (processed() % 10000 == 0) {
		printElement(processed(), numWidth);
		printElement(remaining(), numWidth);
		printElement(candidateEvaluator_.feasibleCandidates(), numWidth);
		printElement(bestUb_, numWidth);
		printElement(bestLb_, numWidth);
		printElement((int)(curTime - startTime_), numWidth);
		cout << endl;
	}
}

bool BBController::stop()
{
	auto curTime = omp_get_wtime();

	return curTime - startTime_ > maxTime_;
}

float BBController::gap(float lb)
{
	if (candidateEvaluator_.feasibleCandidates() > 0) {
		return (bestUb_ - lb) / lb;
	}
	
	return numeric_limits<float>::max();
}

template<typename T> void BBController::printElement(T t, const int& width)
{
	cout << left << setw(width) << setfill(separator) << t;
}
