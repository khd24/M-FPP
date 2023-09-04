#include "CandidateEvaluator.h"

CandidateEvaluator::CandidateEvaluator(float lWidth, float lHeight, vector<Barrier> bList, vector<Barrier> nfList, vector<IOPoint> ioList, vector<int> branchingOrder, float* efefFlow, float* efnfFlow, float* nfnfFlow)
{
	lWidth_ = lWidth;
	lHeight_ = lHeight;
	efefFlow_ = efefFlow;
	efnfFlow_ = efnfFlow;
	nfnfFlow_ = nfnfFlow;

	bList_ = bList;
	nfList_ = nfList;
	ioList_ = ioList;
	branchingOrder_ = branchingOrder;

	feasibleCandidates_ = 0;
	efCount_ = bList.size();
	nfCount_ = nfList.size();
}

CandidateEvaluator::~CandidateEvaluator()
{

}

void CandidateEvaluator::evaluateUb(ObjectiveFunction* gblOptimal, BBNode& bbNode) {

	auto placedNfs = bbNode.placedNfs();

	if (placedNfs.size() != nfCount_) return;

	float* coordinateTuple = new float[2 * (long long)nfCount_];

	for (auto it = placedNfs.begin(); it != placedNfs.end(); ++it) {
		coordinateTuple[2 * it->id()] = it->getMinX();
		coordinateTuple[2 * it->id() + 1] = it->getMaxY();
	}

	evaluateCandidates2(gblOptimal, branchingOrder_.data(), efefFlow_, efnfFlow_, nfnfFlow_, coordinateTuple, bList_, nfList_, ioList_, lWidth_, lHeight_, false);

	gblOptimal->feasibleCandidates = ++feasibleCandidates_;

	bbNode.setUb(gblOptimal->value);

	delete[] coordinateTuple;
}

void CandidateEvaluator::evaluateLb(BBNode& bbNode, CandidateGenerator candidateGenerator) {

	auto placed = bbNode.placedNfs();

	if (placed.size() == nfCount_) {
		bbNode.setLb(bbNode.ub());
		return;
	}

	set<int> placedIds;
	vector<Barrier> toBePlaced;

	for (auto it = placed.begin(); it != placed.end(); ++it) {
		placedIds.insert(it->id());
	}

	auto xcoordMap = bbNode.coordMap(CoordinateType::X);
	auto ycoordMap = bbNode.coordMap(CoordinateType::Y);

	for (int i = 0; i < nfCount_; i++) {
		Barrier nf = nfList_[i];
		nf.setId(i);
		if (placedIds.find(i) == placedIds.end()) {

			auto itx = xcoordMap.find(i);
			auto ity = ycoordMap.find(i);

			if (itx != xcoordMap.end()) {
				nf.setX(itx->second);
			}

			if (ity != ycoordMap.end()) {
				nf.setY(ity->second);
			}

			toBePlaced.push_back(nf);
		}
	}

	auto lbPair = computeFiniteLB(placed, toBePlaced, candidateGenerator);

	if(lbPair.first) bbNode.setLb(lbPair.second);

	bbNode.setFathomed(!lbPair.first);
}

long long int CandidateEvaluator::feasibleCandidates() {
	return feasibleCandidates_;
}

pair<bool, float> CandidateEvaluator::computeFiniteLB(vector<Barrier> placed, vector<Barrier> toBePlaced, CandidateGenerator candidateGenerator) {

	std::vector<Barrier> nflist(placed);

	///////////////////////////////////////////////// Candidate Generation //////////////////////////////////////////////

	float efefLb = 0;
	float efnfLb = 0;
	float nfnfLb = 0;

	bool success = true;

	for (auto it = toBePlaced.begin(); it != toBePlaced.end(); ++it) {

		LBObjectiveFunction optimal;

		bool isFeasible = false;

		Barrier nf = *it;

		auto xcandidates = float_eq(it->getMinX(), -1) ? candidateGenerator.getCandidates(CoordinateType::X, placed) : vector<float>{ it->getMinX() };
		auto ycandidates = float_eq(it->getMaxY(), -1) ? candidateGenerator.getCandidates(CoordinateType::Y, placed) : vector<float>{ it->getMaxY() };

		for (auto itx = xcandidates.begin(); itx != xcandidates.end(); ++itx) {
			auto potentialX = float_eq(it->getMinX(), -1) ? candidateGenerator.getPotentialCoordinates(CoordinateType::X, nf, *itx) : vector<float>{ it->getMinX() };
			for (auto itx2 = potentialX.begin(); itx2 != potentialX.end(); ++itx2) {

				for (auto ity = ycandidates.begin(); ity != ycandidates.end(); ++ity) {
					auto potentialY = float_eq(it->getMaxY(), -1) ? candidateGenerator.getPotentialCoordinates(CoordinateType::Y, nf, *ity) : vector<float>{ it->getMaxY() };
					for (auto ity2 = potentialY.begin(); ity2 != potentialY.end(); ++ity2) {

						nf.setX(*itx2);
						nf.setY(*ity2);

						if (candidateGenerator.isFeasible(placed, nf)) {

							nflist.push_back(nf);

							//Evaluate
							evaluateCandidates3(&optimal, efefFlow_, efnfFlow_, nfnfFlow_, bList_, nflist, ioList_, lWidth_, lHeight_, false, nfCount_);

							nflist.pop_back();

							isFeasible = true;
						}
					}
				}
			}
		}

		if (isFeasible) {
			efefLb = max<float>(efefLb, optimal.efefobj + optimal.efnfobjPlaced + optimal.nfnfobjPlaced);
			efnfLb += optimal.efnfobjToBePlaced + optimal.nfnfobjToBePlaced;
		}

		else
			success = false;
	}

	// Naive NFNF interaction with assumption that io is located at top left corner.
	for (int i = 0; i < toBePlaced.size(); i++) {
		int nfid1 = toBePlaced[i].id();
		float width1 = toBePlaced[i].getWidth();
		float height1 = toBePlaced[i].getHeight();
		
		for (int j = i + 1; j < toBePlaced.size(); j++) {
			auto nfid2 = toBePlaced[j].id();
			float width2 = toBePlaced[j].getWidth();
			float height2 = toBePlaced[j].getHeight();

			auto flow = nfnfFlow_[nfid1 * nfCount_ + nfid2];
			auto distance = min<float>(min<float>(width1, width2), min<float>(height1, height2));

			nfnfLb += flow * distance;
		}
	}

	float lb = success ? efefLb + efnfLb + nfnfLb : numeric_limits<float>::max();

	return pair<bool, float>(success, lb);
}
