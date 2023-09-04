#ifndef _CANDIDATE_REPOSITORY
#define _CANDIDATE_REPOSITORY

#include <map>
#include <vector>
#include <string>
#include <omp.h>

#include "structures_1.h"
#include "functions_1.h"
#include "functions_4.h"

using namespace std;

class CandidateRepository
{
public:
	CandidateRepository(int _nfcount, float _lwidth, float _lheight, vector<Barrier> _blist, vector<Barrier> _nflist, vector<IOPoint> _iolist, vector<int*> _permutations);
	~CandidateRepository();

	void GenerateNFCandidates();

	void EvaluateNFCandidates(ObjectiveFunction *optimal, float *efefflow, float* efnfflow, float* nfnfflow);

	int GetCandidateCount();
	int GetFeasibleCandidateCount();


private:

	int nfcount;
	float lwidth, lheight;
	map<string, vector<float*>> x_coordinates;
	map<string, vector<float*>> y_coordinates;

	long long int candidateCount, feasibleCandidateCount, xCount, yCount;

	bool AddCoordinate(map<string, vector<float*>>& coordinateMap, float* coordinate);

	void GenerateCoordinates(int depth, int *perm, float *xcoords, std::vector<float> v_gridlines, CoordinateType coordinateEnum);

	string GetHashCode(float *coord, int len);

	bool Exists(float* coordinate, vector<float*> coordinates, int length);

	void ProgressDisplay(int &displayNext, int step, long long int totalCandidateCount, int elapsedTime, float incObjVal);

	bool IsFeasible(float *coordinateTuple);

	bool BarrierOverlap(Barrier b1, Barrier b2);

	vector<Barrier> blist;
	vector<Barrier> nflist;
	vector<IOPoint> iolist;

	vector<int*> permutations;

	template <typename T>
	void Clear(vector<T*> std_vector);

	template<typename T1, typename T2>
	void Clear(map<T1, vector<T2*>> std_map);

	template<typename T>
	void PrintVector(vector<T*> std_vector, int len, string name);

	void InputCandidate();
};

#endif
