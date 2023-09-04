#include "CandidateRepository.h"

void CandidateRepository::GenerateNFCandidates() 
{
	vector<float> vgridlinecoord;
	vector<float> hgridlinecoord;

	getGridlinecoord(blist, iolist, vgridlinecoord, hgridlinecoord, lwidth, lheight);

#ifdef _DEBUG
	PrintVector(permutations, nfcount, "permutations");
#endif // DEBUG

	for (auto it = permutations.begin(); it != permutations.end(); ++it) {

		auto perm = *it;

		float *xcoord = new float[nfcount];
		float *ycoord = new float[nfcount];

		GenerateCoordinates(0, perm, xcoord, vgridlinecoord, CoordinateType::X);

		GenerateCoordinates(0, perm, ycoord, hgridlinecoord, CoordinateType::Y);

		delete[] xcoord;
		delete[] ycoord;
	}
}

// Constructor
CandidateRepository::CandidateRepository(int _nfcount, float _lwidth, float _lheight, vector<Barrier> _blist, vector<Barrier> _nflist, vector<IOPoint> _iolist, vector<int*> _permutations)
{
	nfcount = _nfcount;
	lwidth = _lwidth;
	lheight = _lheight;
	blist = _blist;
	nflist = _nflist;
	iolist = _iolist;
	permutations = _permutations;

	candidateCount = 0;
	feasibleCandidateCount = 0;
	xCount = 0;
	yCount = 0;
}

// Destructor
CandidateRepository::~CandidateRepository() {
	Clear<string, float>(x_coordinates);
	Clear<string, float>(y_coordinates);
}

int CandidateRepository::GetCandidateCount() { return candidateCount; }

int CandidateRepository::GetFeasibleCandidateCount() { return feasibleCandidateCount; }

template<typename T1, typename T2>
void CandidateRepository::Clear(map<T1, vector<T2*>> std_map)
{
	for (auto it = std_map.begin(); it != std_map.end(); ++it) {
		Clear<T2>(it->second);
	}

	std_map.clear();
}

template <typename T>
void CandidateRepository::Clear(vector<T*> std_vector)
{
	for (auto it = std_vector.begin(); it != std_vector.end(); ++it) {
		delete[] * it;
	}

	std_vector.clear();
}

template <typename T>
void CandidateRepository::PrintVector(vector<T*> std_vector, int len, string name) {
	
	cout << name << ":" << endl;
	
	for (auto it = std_vector.begin(); it != std_vector.end(); ++it) {
		
		for (int i = 0; i < len - 1; ++i) {
			cout << (*it)[i] << ", ";
		}
		
		cout << (*it)[len - 1] << endl;
	}
}

// This function recursively places each NF from the given permutation and records the new coordinates. At the end of the depth, a tuple is stored. 
void CandidateRepository::GenerateCoordinates(int depth, int *perm, float *coord, vector<float> gridlines, CoordinateType coordinateEnum) {

	if (depth == nfcount) {

		float *coord_to_save = new float[nfcount];
		
		copy(coord, coord + nfcount, coord_to_save);

		if (coordinateEnum == CoordinateType::X) {
			bool added = AddCoordinate(x_coordinates, coord_to_save);
			if (added) ++xCount;
		}

		else if (coordinateEnum == CoordinateType::Y) {
			bool added = AddCoordinate(y_coordinates, coord_to_save);
			if(added) ++yCount;
		}

		return;
	}

	for (auto it = gridlines.begin(); it != gridlines.end(); ++it) {

		auto nfid = perm[depth];

		auto nf = nflist[nfid];

		auto coordinate = *it;

		vector<float> pontentialCoordinates;

		if(float_lt(coordinate, 0))
			continue;

		if (coordinateEnum == CoordinateType::X) {

			if (float_gt(coordinate, lwidth))
				continue;
			
			if(!float_gt(coordinate + nf.getWidth(), lwidth))
				pontentialCoordinates.push_back(coordinate);

			if(!float_lt(coordinate - nf.getWidth(), 0))
				pontentialCoordinates.push_back(coordinate - nf.getWidth());
		}

		else if (coordinateEnum == CoordinateType::Y) {

			if (float_gt(coordinate, lheight))
				continue;

			if (!float_lt(coordinate - nf.getHeight(), 0))
				pontentialCoordinates.push_back(coordinate);

			if (!float_gt(coordinate + nf.getHeight(), lheight))
				pontentialCoordinates.push_back(coordinate + nf.getHeight());
		}
		
		for (auto it2 = pontentialCoordinates.begin(); it2 != pontentialCoordinates.end(); ++it2) {

			coord[nfid] = *it2;

			vector<float> temp_gridlines(gridlines);

			temp_gridlines.push_back(*it2);

			if (coordinateEnum == CoordinateType::X)
				temp_gridlines.push_back(*it2 + nf.getWidth());

			else if (coordinateEnum == CoordinateType::Y)
				temp_gridlines.push_back(*it2 - nf.getHeight());

			removeDuplicates(temp_gridlines);

			GenerateCoordinates(depth + 1, perm, coord, temp_gridlines, coordinateEnum);
		}
	}
}

// This function returns a string hashcode composed of integer components of a list of floating point numbers.
string CandidateRepository::GetHashCode(float *coord, int len) {

	ostringstream s;

	for (int i = 0; i < len - 1; i++) {
		s << ((int)coord[i]) << "_";
	}

	s << ((int)coord[len - 1]);

	return s.str();
}

// This function adds a coordinate to the given map if the coordinate doesn't already exist.
bool CandidateRepository::AddCoordinate(map<string, vector<float*>>& coordinateMap, float* coordinate) {

	string code = GetHashCode(coordinate, nfcount);

	auto it = coordinateMap.find(code);

	if (it == coordinateMap.end()) {
		coordinateMap.insert(pair < string, vector<float*>>(code, vector<float*>()));
	}

	if (!Exists(coordinate, coordinateMap[code], nfcount)) {
		coordinateMap[code].push_back(coordinate);
		return true;
	}

	return false;
}

// This function returns true if the given coordinate tuple exists in the vector.
bool CandidateRepository::Exists(float *coordinate, vector<float*> coordinates, int length) {

	bool ret = false;

	for (auto it = coordinates.begin(); it != coordinates.end(); ++it) {
		
		ret = true;

		auto coordinateToCompare = *it;

		for (int i = 0; i < length; i++) {
			
			if (!float_eq(coordinate[i], coordinateToCompare[i])) {

				ret = false;

				break;
			}
		}
	}

	return ret;
}

void CandidateRepository::EvaluateNFCandidates(ObjectiveFunction *optimal, float* efefflow, float* efnfflow, float* nfnfflow) {

	auto startTime = omp_get_wtime();

	int step = 1;
	int displayNext = 0;
	int percent = 0;
	long long int totalCandidateCount = xCount * yCount;

	ProgressDisplay(displayNext, 0, totalCandidateCount, 0, optimal->value);

	displayNext = step;

	for (auto it1 = x_coordinates.begin(); it1 != x_coordinates.end(); ++it1) {
		for (auto it2 = (it1->second).begin(); it2 != (it1->second).end(); ++it2) {
			for (auto it3 = y_coordinates.begin(); it3 != y_coordinates.end(); ++it3) {
				for (auto it4 = (it3->second).begin(); it4 != (it3->second).end(); ++it4) {

					auto coordinateTuple = new float[2 * nfcount];

					for (int i = 0; i < nfcount; ++i) {
						coordinateTuple[2 * i] = (*it2)[i];
						coordinateTuple[2 * i + 1] = (*it4)[i];
					}

					optimal->candidatecount = ++candidateCount;

					if (IsFeasible(coordinateTuple)) {
						evaluateCandidates2(optimal, permutations[0], efefflow, efnfflow, nfnfflow, coordinateTuple, blist, nflist, iolist, lwidth, lheight, false);
						optimal->feasibleCandidates = ++feasibleCandidateCount;
					}

					int elapsedTime = (int) (omp_get_wtime() - startTime);

					ProgressDisplay(displayNext, step, totalCandidateCount, elapsedTime, optimal->value);

					delete[] coordinateTuple;

				}
			}
		}
	}
}

// Test function
void CandidateRepository::ProgressDisplay(int &displayNext, int step, long long int totalCandidateCount, int elapsedTime, float incObjVal)
{
	// Formatted progress indicator
	int percent = (int) (100 * ((double) candidateCount / totalCandidateCount));
	if (percent >= displayNext)
	{
		cout << "\r" << "[" << std::string(percent / 5, '=') << std::string(100 / 5 - percent / 5, ' ') << "] ";
		cout << percent << "%" << " [Candidates " << candidateCount << " of " << totalCandidateCount << "] " << "[nFeasible: " << feasibleCandidateCount << "] " << "[IncObjval: " << incObjVal <<"]" << " [" << elapsedTime << " s]";
		std::cout.flush();
		displayNext += step;
	}
}

bool CandidateRepository::IsFeasible(float *coordinateTuple) {

	vector<Barrier> blistNew(blist);

	for (int i = 0; i < nfcount; i++) {
		Barrier nf;
		nf.setBarrier(coordinateTuple[2 * i], coordinateTuple[2 * i + 1], nflist[i].getWidth(), nflist[i].getHeight(), nflist[i].getCongestion());

		if (float_lt(nf.getMinX(), 0) || float_gt(nf.getMaxX(), lwidth))
			return false;

		if (float_lt(nf.getMinY(), 0) || float_gt(nf.getMaxY(), lheight))
			return false;

		for (auto it = blistNew.begin(); it != blistNew.end(); ++it) {
			
			auto ef = *it;

			if (BarrierOverlap(nf, ef)) {
				
				return false;
			}
		}

		blistNew.push_back(nf);
	}

	return true;
}

bool CandidateRepository::BarrierOverlap(Barrier b1, Barrier b2) {

	if(!float_lt(b1.getMinX(), b2.getMaxX()) || !float_lt(b2.getMinX(), b1.getMaxX()))
		return false;

	if (!float_lt(b1.getMinY(), b2.getMaxY()) || !float_lt(b2.getMinY(), b1.getMaxY()))
		return false;

	return true;
}

void CandidateRepository::InputCandidate() {
	//(283.031738, 373.087494, 163.618912, 263.159058, 159.817291, 353.500702)
	//(45.1103, 169.68)(39.5851, 169.68)(39.5851, 169.68)

	float xCoord[3] = { 45.1103, 39.5851, 39.5851 };
	float yCoord[3] = { 169.68, 169.68, 169.68 };
	float coordTuple[6] = { 45.1103, 169.68, 39.5851, 169.68, 39.5851, 169.68 };

	auto xHash = GetHashCode(xCoord, 3);
	auto yHash = GetHashCode(yCoord, 3);

	auto it1 = x_coordinates.find(xHash);
	auto it3 = y_coordinates.find(yHash);

	if(it1 != x_coordinates.end()) {
		for (auto it2 = (it1->second).begin(); it2 != (it1->second).end(); ++it2) {
			if (it3 != y_coordinates.end()) {
				for (auto it4 = (it3->second).begin(); it4 != (it3->second).end(); ++it4) {

					float coordinateTuple[6];

					for (int i = 0; i < nfcount; ++i) {
						coordinateTuple[2 * i] = (*it2)[i];
						coordinateTuple[2 * i + 1] = (*it4)[i];
					}

					++candidateCount;

					if (IsFeasible(coordinateTuple))
						++feasibleCandidateCount;

					cout << candidateCount << "\t" << feasibleCandidateCount << endl;
				}
			}
		}
	}
}

