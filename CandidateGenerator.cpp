#include "CandidateGenerator.h"

CandidateGenerator::CandidateGenerator(vector<Barrier> bList, vector<Barrier> nfList, vector<IOPoint> ioList, float lWidth, float lHeight)
{
	lWidth_ = lWidth;
	lHeight_ = lHeight;
	bList_ = bList;
	nfList_ = nfList;
	getGridlinecoord(bList, ioList, xGridlines_, yGridlines_, lWidth, lHeight);
}

CandidateGenerator::~CandidateGenerator()
{
}

vector<float> CandidateGenerator::getCandidates(CoordinateType coordinateType, vector<Barrier> placedNfs) {

	vector<float> candidates;

	for (auto it = placedNfs.begin(); it != placedNfs.end(); ++it) {
		switch (coordinateType)
		{
		case CoordinateType::X:
			candidates.push_back(it->getMinX());
			candidates.push_back(it->getMaxX());
			break;

		case CoordinateType::Y:
			candidates.push_back(it->getMinY());
			candidates.push_back(it->getMaxY());
			break;
		}
	}

	switch (coordinateType)
	{
	case CoordinateType::X:
		return vectorUnion(candidates, xGridlines_);

	case CoordinateType::Y:
		return vectorUnion(candidates, yGridlines_);
	}
}

vector<float> CandidateGenerator::getCandidatesPolytomic(CoordinateType coordinateType, vector<Barrier> nfs, int nfOrder) {

	int order = -1;
	vector<float> candidates;

	switch (coordinateType)
	{
	case CoordinateType::X:
		candidates.insert(candidates.end(), xGridlines_.begin(), xGridlines_.end());
		break;

	case CoordinateType::Y:
		candidates.insert(candidates.end(), yGridlines_.begin(), yGridlines_.end());
		break;
	}

	for (auto it = nfs.begin(); it != nfs.end(); ++it) {

		vector<float> newCandidates;

		switch (coordinateType)
		{
		case CoordinateType::X:
			if (!binary_search(candidates.begin(), candidates.end(), it->getMinX()))
				newCandidates.push_back(it->getMinX());

			if (!binary_search(candidates.begin(), candidates.end(), it->getMaxX()))
				newCandidates.push_back(it->getMaxX());

			break;

		case CoordinateType::Y:
			if (!binary_search(candidates.begin(), candidates.end(), it->getMinY()))
				newCandidates.push_back(it->getMinY());

			if (!binary_search(candidates.begin(), candidates.end(), it->getMaxY()))
				newCandidates.push_back(it->getMaxY());

			break;
		}

		if (nfOrder > it->order())
			candidates = vectorUnion(candidates, newCandidates);
		else
			candidates = newCandidates;
	}

	return candidates;
}

vector<float> CandidateGenerator::getPotentialCoordinates(CoordinateType coordinateType, Barrier nf, float coordinate) {

	vector<float> newCoordinates;

	if (!float_lt(coordinate, 0)) {

		switch (coordinateType)
		{
		case CoordinateType::X:
			if (!float_gt(coordinate, lWidth_)) {
				if (!float_lt(coordinate - nf.getWidth(), 0))
					newCoordinates.push_back(coordinate - nf.getWidth());

				if (!float_gt(coordinate + nf.getWidth(), lWidth_))
					newCoordinates.push_back(coordinate);
			}
			break;

		case CoordinateType::Y:
			if (!float_gt(coordinate, lHeight_)) {

				if (!float_lt(coordinate - nf.getHeight(), 0))
					newCoordinates.push_back(coordinate);

				if (!float_gt(coordinate + nf.getHeight(), lHeight_))
					newCoordinates.push_back(coordinate + nf.getHeight());
			}
			break;
		}
	}

	return newCoordinates;
}

vector<Barrier> CandidateGenerator::getPotentialPolytomicCoordinates(CoordinateType coordinateType, vector<Barrier> nfs, int newNfId, int nfOrder) {

	vector<Barrier> retVector;
	Barrier newNf = nfList_[newNfId];

	newNf.setOrder(nfOrder);
	newNf.setId(newNfId);

	auto candidates = getCandidatesPolytomic(coordinateType, nfs, nfOrder);

	for (auto it = candidates.begin(); it != candidates.end(); ++it) {
		auto potentialCoordinates = getPotentialCoordinates(coordinateType, newNf, *it);
		for (auto it2 = potentialCoordinates.begin(); it2 != potentialCoordinates.end(); ++it2) {

			switch (coordinateType)
			{
			case CoordinateType::X:
				newNf.setX(*it2);
				break;
			case CoordinateType::Y:
				newNf.setY(*it2);
				break;
			default:
				break;
			}
			
			retVector.push_back(newNf);
		}
	}

	return retVector;
}

bool CandidateGenerator::isFeasible(vector<Barrier> nfs, Barrier nf) {

	if (float_lt(nf.getMinX(), 0) || float_gt(nf.getMaxX(), lWidth_))
		return false;

	if (float_lt(nf.getMinY(), 0) || float_gt(nf.getMaxY(), lHeight_))
		return false;

	return !barrierOverlap(bList_, nf) && !barrierOverlap(nfs, nf);
}