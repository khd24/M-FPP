/*
 * functions_alg.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat
 */

#include "functions_alg.h"

 // Optimal Procedure + Heuristic 5
void executeCornerOptimalHeuristic_allPerm(Layout *layout, vector<int*> permutations, float* times, ObjectiveFunction *optimal) {

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;
	int iocount = layout->iocount;

	float lwidth = layout->width;
	float lheight = layout->height;
	vector<Barrier> barrierlist(layout->barrierlist, layout->barrierlist + efcount);
	vector<IOPoint> iolist(layout->iolist, layout->iolist + iocount);
	vector<Barrier> nflist(layout->nflist, layout->nflist + nfcount);

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;

	///////////////////////////////////////////////// Candidate Generation //////////////////////////////////////////////

	cout << "CORNER_HEU_ALL_PERM "; //(total perm: " << permutations.size() << ")" << endl;

	float time1 = (float)omp_get_wtime();

	for (int k = 0; k < permutations.size(); k++) {

		bool inf_flag = false;
		float *candidate = new float[2 * nfcount];
		generateCandidates(optimal, permutations.at(k), 0, nfcount, inf_flag, candidate, barrierlist, nflist, iolist, efefflow, efnfflow, nfnfflow, lwidth, lheight);
		delete[] candidate;
	}

	float time2 = (float)omp_get_wtime();

	//	cout << "cgen done! count: " << candidatelist.size() << " time: (" << time2 - time1 << " s) ";

	times[0] = time2 - time1;

	//////////////////////////////////////////////// Candidate Evaluation /////////////////////////////

//	cout << "C_EVAL (total cnds: " << candidatelist.size() << ")" << endl;

//	(*optimal).candidatecount = candidatelist.size();

/*	for (int i = 0; i < candidatelist.size(); i++) {

		float* coord = candidatelist.at(i);
		int *perm = permlist.at(i);

		evaluateCandidates(optimal, perm, efefflow, efnfflow, nfnfflow, coord, barrierlist, nflist, iolist, lwidth, lheight, false);

		delete[] candidatelist.at(i);
	}*/

	//////////////////////////////////////////////////////////////////////////////////////////////////

	float time3 = (float)omp_get_wtime();

	cout << "ceval done! (opt: " << (*optimal).value << " cnd: " << (*optimal).candidatecount << " time: " << time3 - time1 << " s)" << endl;

	times[1] = time3 - time2;
	times[2] = time3 - time1;

}

// Heuristic 1 + Heuristic 5
void executeCornerOptimalHeuristic_givenPerm(Layout *layout, int* permutation, float* times, ObjectiveFunction *optimal) {

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;
	int iocount = layout->iocount;

	float lwidth = layout->width;
	float lheight = layout->height;

	vector<Barrier> barrierlist(layout->barrierlist, layout->barrierlist + efcount);
	vector<IOPoint> iolist(layout->iolist, layout->iolist + iocount);
	vector<Barrier> nflist(layout->nflist, layout->nflist + nfcount);

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;

	///////////////////////////////////////////////// Candidate Generation //////////////////////////////////////////////

	cout << "CORNER_HEU_GIVEN_PERM ";

	float time1 = (float)omp_get_wtime();

	bool inf_flag = false;
	float *candidate = new float[2 * nfcount];
	generateCandidates(optimal, permutation, 0, nfcount, inf_flag, candidate, barrierlist, nflist, iolist, efefflow, efnfflow, nfnfflow, lwidth, lheight);
	delete[] candidate;


	float time2 = (float)omp_get_wtime();

	//	cout << "cgen done! count: " << (*optimal).candidatecount << " time: (" << time2 - time1 << " s) ";


	times[0] = time2 - time1;

	//////////////////////////////////////////////// Candidate Evaluation /////////////////////////////

//	cout << "C_EVAL (total cnds: " << candidatelist.size() << ")" << endl;

//	(*optimal).candidatecount = candidatelist.size();

/*
	for (int i = 0; i < candidatelist.size(); i++) {

		float* coord = candidatelist.at(i);

		int *perm = permlist.at(i);

		evaluateCandidates(optimal, permutation , efefflow, efnfflow, nfnfflow, coord, barrierlist, nflist, iolist, lwidth, lheight, false);

		delete[] candidatelist.at(i);
	}
*/

//////////////////////////////////////////////////////////////////////////////////////////////////

	float time3 = (float)omp_get_wtime();

	cout << "ceval done! (opt: " << (*optimal).value << " cnd: " << (*optimal).candidatecount << " time: " << time3 - time1 << " s)" << endl;

	times[1] = time3 - time2;
	times[2] = time3 - time1;

}

// Heuristic 2 + Heuristic 5
void executeClass1(Layout *layout, vector<int*> &permutations, float* times, ObjectiveFunction *optimal, int bucketsize) {

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;
	int iocount = layout->iocount;

	float lwidth = layout->width;
	float lheight = layout->height;

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;

	if (permutations.size() > 1)
		cout << "HEURISTIC_2_ALL_PERM_BUCKETSIZE_" << bucketsize << " "; //(total perm: " << permutations.size() << ")" << endl;
	else
		cout << "HEURISTIC_3_SINGLE_PERM_BUCKETSIZE_" << bucketsize << " "; //(total perm: " << permutations.size() << ")" << endl;

	///////////////////////////////////////////////// Candidate Generation //////////////////////////////////////////////

	float time1 = (float)omp_get_wtime();

	for (vector<int*>::iterator it1 = permutations.begin(); it1 != permutations.end(); ++it1) {

		ObjectiveFunction lcl_opt1;
		lcl_opt1.ptr = new float[2 * nfcount];
		fill(lcl_opt1.ptr, lcl_opt1.ptr + 2 * nfcount, -1.0);
		lcl_opt1.value = 0;
		lcl_opt1.candidatecount = 0;
		lcl_opt1.perm = new int[nfcount];

		int M_by_m = nfcount / bucketsize;
		int extra_m = nfcount % bucketsize;
		int extra_size = (extra_m > 0) ? 1 : 0;

		vector<Barrier> barrierlist(layout->barrierlist, layout->barrierlist + efcount);
		vector<IOPoint> iolist(layout->iolist, layout->iolist + iocount);
		vector<Barrier> nflist(layout->nflist, layout->nflist + nfcount);

		vector<int> placed;

		for (int i = 0; i < M_by_m + extra_size; i++) {

			vector<int> to_be_placed;
			vector<Barrier> newnflist;

			int newbucketsize = i < M_by_m ? bucketsize : extra_m;

			///////////// Permute intra-bucket NFs ///////////////
			vector<int*> red_perm;
			int *a = new int[newbucketsize];
			for (int j = 0; j < newbucketsize; j++) {
				a[j] = j;
				int nfid = (*it1)[i * bucketsize + j];
				Barrier nf = nflist.at(nfid);
				newnflist.push_back(nf);
				to_be_placed.push_back(nfid);
			}

			red_perm.push_back(a);

			/////////////////////////////////////////////////////

			float *new_efefflowmatrix = new float[(efcount + placed.size()) * (efcount + placed.size())];
			float *new_efnfflowmatrix = new float[(efcount + placed.size()) * to_be_placed.size()];
			float *new_nfnfflowmatrix = new float[to_be_placed.size() * to_be_placed.size()];

			updateFlowMatrices(new_efefflowmatrix, new_efnfflowmatrix, new_nfnfflowmatrix, efefflow, efnfflow, nfnfflow, efcount, nfcount, placed, to_be_placed);

			ObjectiveFunction lcl_opt2;
			lcl_opt2.ptr = new float[2 * newbucketsize];
			fill(lcl_opt2.ptr, lcl_opt2.ptr + 2 * newbucketsize, -1.0);
			lcl_opt2.value = FLT_MAX;
			lcl_opt2.candidatecount = 0;
			lcl_opt2.perm = new int[newbucketsize];

			for (vector<int*>::iterator it2 = red_perm.begin(); it2 != red_perm.end(); ++it2) {

				bool inf_flag = false;
				float *candidate = new float[2 * newbucketsize];
				generateCandidates(&lcl_opt2, *it2, 0, newbucketsize, inf_flag, candidate, barrierlist, newnflist, iolist, new_efefflowmatrix, new_efnfflowmatrix, new_nfnfflowmatrix, lwidth, lheight);
				delete[] candidate;
				delete[](*it2);
				(*it2) = 0;
			}

			if (lcl_opt2.candidatecount > 0) {
				for (int k = 0; k < to_be_placed.size(); k++) {
					int nfid = to_be_placed.at(k);
					lcl_opt1.ptr[2 * nfid] = lcl_opt2.ptr[2 * k];
					lcl_opt1.ptr[2 * nfid + 1] = lcl_opt2.ptr[2 * k + 1];
					lcl_opt1.perm[i * bucketsize + k] = to_be_placed.at(lcl_opt2.perm[k]);
				}
			}


			lcl_opt1.value = lcl_opt2.candidatecount > 0 ? lcl_opt2.value : FLT_MAX;

			lcl_opt1.candidatecount = lcl_opt2.candidatecount > 0 ? lcl_opt1.candidatecount + lcl_opt2.candidatecount : 0;

			delete[] lcl_opt2.ptr;
			delete[] lcl_opt2.perm;
			lcl_opt2.ptr = 0;
			lcl_opt2.perm = 0;

			delete[] new_efefflowmatrix;
			delete[] new_efnfflowmatrix;
			delete[] new_nfnfflowmatrix;

			if (lcl_opt2.candidatecount == 0)
				break;

			for (int k = 0; k < newnflist.size(); k++) {
				int nfid = to_be_placed.at(k);
				Barrier newb;
				newb.setBarrier(lcl_opt1.ptr[2 * nfid], lcl_opt1.ptr[2 * nfid + 1], newnflist.at(k).getWidth(), newnflist.at(k).getHeight(), newnflist.at(k).getCongestion());
				barrierlist.push_back(newb);

				IOPoint newio;
				newio.x = lcl_opt1.ptr[2 * nfid];
				newio.y = lcl_opt1.ptr[2 * nfid + 1];
				newio.isNF = false;
				iolist.push_back(newio);

				placed.push_back(to_be_placed.at(k));
			}

		}

		if (lcl_opt1.candidatecount > 0) {
			optimal->candidatecount += lcl_opt1.candidatecount;

			if (lcl_opt1.value < optimal->value) {
				optimal->value = lcl_opt1.value;
				copy(lcl_opt1.ptr, lcl_opt1.ptr + 2 * nfcount, optimal->ptr);
				copy(lcl_opt1.perm, lcl_opt1.perm + nfcount, optimal->perm);

			}
		}

		delete[] lcl_opt1.ptr;
		delete[] lcl_opt1.perm;
		lcl_opt1.ptr = 0;
		lcl_opt1.perm = 0;
	}

	float time2 = (float)omp_get_wtime();

	//	cout << "cgen done! count: " << candidatelist.size() << " time: (" << time2 - time1 << " s) ";

	times[0] = time2 - time1;

	//////////////////////////////////////////////// Candidate Evaluation /////////////////////////////

	//	cout << "C_EVAL (total cnds: " << candidatelist.size() << ")" << endl;

	//	(*optimal).candidatecount = candidatelist.size();

	/*	for (int i = 0; i < candidatelist.size(); i++) {

	float* coord = candidatelist.at(i);
	int *perm = permlist.at(i);

	evaluateCandidates(optimal, perm, efefflow, efnfflow, nfnfflow, coord, barrierlist, nflist, iolist, lwidth, lheight, false);

	delete[] candidatelist.at(i);
	}*/

	//////////////////////////////////////////////////////////////////////////////////////////////////

	float time3 = (float)omp_get_wtime();

	cout << "ceval done! (opt: " << (*optimal).value << " cnd: " << (*optimal).candidatecount << " time: " << time3 - time1 << " s)" << endl;

	times[1] = time3 - time2;
	times[2] = time3 - time1;

}

// Heuristic 2 + Heuristic 5
void executeClass2(Layout *layout, vector<int*> &permutations, float* times, ObjectiveFunction *optimal, int bucketsize) {

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;
	int iocount = layout->iocount;

	float lwidth = layout->width;
	float lheight = layout->height;

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;

	if (permutations.size() > 1)
		cout << "HEURISTIC_4_ALL_PERM_BUCKETSIZE_" << bucketsize << " "; //(total perm: " << permutations.size() << ")" << endl;
	else
		cout << "HEURISTIC_5_SINGLE_PERM_BUCKETSIZE_" << bucketsize << " "; //(total perm: " << permutations.size() << ")" << endl;

	getUniquePermutations(permutations, nfcount, bucketsize);

	///////////////////////////////////////////////// Candidate Generation //////////////////////////////////////////////

	float time1 = (float)omp_get_wtime();

	for (vector<int*>::iterator it1 = permutations.begin(); it1 != permutations.end(); ++it1) {

		ObjectiveFunction lcl_opt1;
		lcl_opt1.ptr = new float[2 * nfcount];
		fill(lcl_opt1.ptr, lcl_opt1.ptr + 2 * nfcount, -1.0);
		lcl_opt1.value = 0;
		lcl_opt1.candidatecount = 0;
		lcl_opt1.perm = new int[nfcount];

		int M_by_m = nfcount / bucketsize;
		int extra_m = nfcount % bucketsize;
		int extra_size = (extra_m > 0) ? 1 : 0;

		vector<Barrier> barrierlist(layout->barrierlist, layout->barrierlist + efcount);
		vector<IOPoint> iolist(layout->iolist, layout->iolist + iocount);
		vector<Barrier> nflist(layout->nflist, layout->nflist + nfcount);

		vector<int> placed;

		for (int i = 0; i < M_by_m + extra_size; i++) {

			vector<int> to_be_placed;
			vector<Barrier> newnflist;

			int newbucketsize = i < M_by_m ? bucketsize : extra_m;

			///////////// Permute intra-bucket NFs ///////////////
			vector<int*> red_perm;
			int *a = new int[newbucketsize];
			for (int j = 0; j < newbucketsize; j++) {
				a[j] = j;
				int nfid = (*it1)[i * bucketsize + j];
				Barrier nf = nflist.at(nfid);
				newnflist.push_back(nf);
				to_be_placed.push_back(nfid);
			}

			permute(red_perm, a, 0, newbucketsize);
			delete[] a;

			/////////////////////////////////////////////////////

			float *new_efefflowmatrix = new float[(efcount + placed.size()) * (efcount + placed.size())];
			float *new_efnfflowmatrix = new float[(efcount + placed.size()) * to_be_placed.size()];
			float *new_nfnfflowmatrix = new float[to_be_placed.size() * to_be_placed.size()];

			updateFlowMatrices(new_efefflowmatrix, new_efnfflowmatrix, new_nfnfflowmatrix, efefflow, efnfflow, nfnfflow, efcount, nfcount, placed, to_be_placed);

			ObjectiveFunction lcl_opt2;
			lcl_opt2.ptr = new float[2 * newbucketsize];
			fill(lcl_opt2.ptr, lcl_opt2.ptr + 2 * newbucketsize, -1.0);
			lcl_opt2.value = FLT_MAX;
			lcl_opt2.candidatecount = 0;
			lcl_opt2.perm = new int[newbucketsize];

			for (vector<int*>::iterator it2 = red_perm.begin(); it2 != red_perm.end(); ++it2) {

				bool inf_flag = false;
				float *candidate = new float[2 * newbucketsize];
				generateCandidates(&lcl_opt2, *it2, 0, newbucketsize, inf_flag, candidate, barrierlist, newnflist, iolist, new_efefflowmatrix, new_efnfflowmatrix, new_nfnfflowmatrix, lwidth, lheight);
				delete[] candidate;
				delete[](*it2);
				(*it2) = 0;
			}

			if (lcl_opt2.candidatecount > 0) {
				for (int k = 0; k < to_be_placed.size(); k++) {
					int nfid = to_be_placed.at(k);
					lcl_opt1.ptr[2 * nfid] = lcl_opt2.ptr[2 * k];
					lcl_opt1.ptr[2 * nfid + 1] = lcl_opt2.ptr[2 * k + 1];
					lcl_opt1.perm[i * bucketsize + k] = to_be_placed.at(lcl_opt2.perm[k]);
				}
			}

			lcl_opt1.value = lcl_opt2.candidatecount > 0 ? lcl_opt2.value : FLT_MAX;

			lcl_opt1.candidatecount = lcl_opt2.candidatecount > 0 ? lcl_opt1.candidatecount + lcl_opt2.candidatecount : 0;

			delete[] lcl_opt2.ptr;
			delete[] lcl_opt2.perm;
			lcl_opt2.ptr = 0;
			lcl_opt2.perm = 0;

			delete[] new_efefflowmatrix;
			delete[] new_efnfflowmatrix;
			delete[] new_nfnfflowmatrix;

			if (lcl_opt2.candidatecount == 0)
				break;

			for (int k = 0; k < newnflist.size(); k++) {
				int nfid = to_be_placed.at(k);
				Barrier newb;
				newb.setBarrier(lcl_opt1.ptr[2 * nfid], lcl_opt1.ptr[2 * nfid + 1], newnflist.at(k).getWidth(), newnflist.at(k).getHeight(), newnflist.at(k).getCongestion());
				barrierlist.push_back(newb);

				IOPoint newio;
				newio.x = lcl_opt1.ptr[2 * nfid];
				newio.y = lcl_opt1.ptr[2 * nfid + 1];
				newio.isNF = false;
				iolist.push_back(newio);

				placed.push_back(to_be_placed.at(k));
			}

		}

		if (lcl_opt1.candidatecount > 0) {
			optimal->candidatecount += lcl_opt1.candidatecount;

			if (lcl_opt1.value < optimal->value) {
				optimal->value = lcl_opt1.value;
				copy(lcl_opt1.ptr, lcl_opt1.ptr + 2 * nfcount, optimal->ptr);
				copy(lcl_opt1.perm, lcl_opt1.perm + nfcount, optimal->perm);

			}
		}

		delete[] lcl_opt1.ptr;
		delete[] lcl_opt1.perm;
		lcl_opt1.ptr = 0;
		lcl_opt1.perm = 0;
	}

	float time2 = (float)omp_get_wtime();

	//	cout << "cgen done! count: " << candidatelist.size() << " time: (" << time2 - time1 << " s) ";

	times[0] = time2 - time1;

	//////////////////////////////////////////////// Candidate Evaluation /////////////////////////////

	//	cout << "C_EVAL (total cnds: " << candidatelist.size() << ")" << endl;

	//	(*optimal).candidatecount = candidatelist.size();

	/*	for (int i = 0; i < candidatelist.size(); i++) {

	float* coord = candidatelist.at(i);
	int *perm = permlist.at(i);

	evaluateCandidates(optimal, perm, efefflow, efnfflow, nfnfflow, coord, barrierlist, nflist, iolist, lwidth, lheight, false);

	delete[] candidatelist.at(i);
	}*/

	//////////////////////////////////////////////////////////////////////////////////////////////////

	float time3 = (float)omp_get_wtime();

	cout << "ceval done! (opt: " << (*optimal).value << " cnd: " << (*optimal).candidatecount << " time: " << time3 - time1 << " s)" << endl;

	times[1] = time3 - time2;
	times[2] = time3 - time1;

}

// Repair Heuristic
void executeRepair(Layout *layout, float* times, ObjectiveFunction *optimal, int bucketsize, int max_itn) {

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;
	int iocount = layout->iocount;

	float lwidth = layout->width;
	float lheight = layout->height;

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;

	set<int> nfids;
	for (int i = 0; i < nfcount; i++)
		nfids.insert(i);

	cout << "REPAIR_HEURISTIC_BUCKETSIZE_" << bucketsize << endl; //(total perm: " << permutations.size() << ")" << endl;

	vector<vector<int>> combinations;
	vector<int> combination;
	int count = 0;
	getUniqueCombinations(combinations, combination, nfcount, bucketsize, 0, 0, count);

	///////////////////////////////////////////////// Candidate Generation //////////////////////////////////////////////

	float time1 = (float)omp_get_wtime();

	set<int> visited_combinations;

	int itn_count = 0;

	while (visited_combinations.size() < min<int>(max_itn, combinations.size())) {

		vector<Barrier> barrierlist(layout->barrierlist, layout->barrierlist + efcount);
		vector<IOPoint> iolist(layout->iolist, layout->iolist + iocount);
		vector<Barrier> nflist(layout->nflist, layout->nflist + nfcount);

		random_shuffle(combinations.begin() + visited_combinations.size(), combinations.end());
		vector<int> comb = combinations.at(visited_combinations.size());
		int comb_id = comb.at(bucketsize);

		cout << "REPAIR ITN " << ++itn_count << "\t<";
		for (int i = 0; i < bucketsize - 1; i++)
			cout << comb.at(i) + 1 << ", ";
		cout << comb.at(bucketsize - 1) + 1 << ">\t";

		set<int> placed_nfids(nfids);

		vector<int> to_be_placed;
		vector<Barrier> newnflist;

		vector<int*> red_perm;
		int *a = new int[bucketsize];
		for (int j = 0; j < bucketsize; j++) {
			a[j] = j;
			int nfid = comb.at(j);
			Barrier nf = nflist.at(nfid);
			newnflist.push_back(nf);
			to_be_placed.push_back(nfid);
			placed_nfids.erase(nfid);
		}

		permute(red_perm, a, 0, bucketsize);
		delete[] a;

		vector<int> placed(placed_nfids.begin(), placed_nfids.end());

		for (int k = 0; k < placed.size(); k++) {
			int nfid = placed.at(k);
			Barrier newb;
			newb.setBarrier(optimal->ptr[2 * nfid], optimal->ptr[2 * nfid + 1], nflist.at(nfid).getWidth(), nflist.at(nfid).getHeight(), nflist.at(nfid).getCongestion());
			barrierlist.push_back(newb);

			IOPoint newio;
			newio.x = optimal->ptr[2 * nfid];
			newio.y = optimal->ptr[2 * nfid + 1];
			newio.isNF = false;
			iolist.push_back(newio);
		}

		float *new_efefflowmatrix = new float[(efcount + placed.size()) * (efcount + placed.size())];
		float *new_efnfflowmatrix = new float[(efcount + placed.size()) * to_be_placed.size()];
		float *new_nfnfflowmatrix = new float[to_be_placed.size() * to_be_placed.size()];

		updateFlowMatrices(new_efefflowmatrix, new_efnfflowmatrix, new_nfnfflowmatrix, efefflow, efnfflow, nfnfflow, efcount, nfcount, placed, to_be_placed);

		ObjectiveFunction lcl_opt2;
		lcl_opt2.ptr = new float[2 * bucketsize];
		fill(lcl_opt2.ptr, lcl_opt2.ptr + 2 * bucketsize, -1.0);
		lcl_opt2.value = FLT_MAX;
		lcl_opt2.candidatecount = 0;
		lcl_opt2.perm = new int[bucketsize];

		for (vector<int*>::iterator it2 = red_perm.begin(); it2 != red_perm.end(); ++it2) {

			bool inf_flag = false;
			float *candidate = new float[2 * bucketsize];
			generateCandidates(&lcl_opt2, *it2, 0, bucketsize, inf_flag, candidate, barrierlist, newnflist, iolist, new_efefflowmatrix, new_efnfflowmatrix, new_nfnfflowmatrix, lwidth, lheight);
			delete[] candidate;
			delete[](*it2);
			(*it2) = 0;
		}

		if (lcl_opt2.value < optimal->value) {
			optimal->value = lcl_opt2.value;
			for (int k = 0; k < to_be_placed.size(); k++) {
				int nfid = to_be_placed.at(k);
				optimal->ptr[2 * nfid] = lcl_opt2.ptr[2 * k];
				optimal->ptr[2 * nfid + 1] = lcl_opt2.ptr[2 * k + 1];
			}

			//			visited_combinations.clear();
		}

		//		else 
		visited_combinations.insert(comb_id);

		float time3 = (float)omp_get_wtime();

		cout << optimal->value << "\t" << time3 - time1 << " s" << endl;

		delete[] lcl_opt2.ptr;
		delete[] lcl_opt2.perm;
		lcl_opt2.ptr = 0;
		lcl_opt2.perm = 0;

		delete[] new_efefflowmatrix;
		delete[] new_efnfflowmatrix;
		delete[] new_nfnfflowmatrix;

	}

	float time4 = (float)omp_get_wtime();

	times[1] += time4 - time1;
	times[2] += time4 - time1;

}

void generatePermutation_ALDEP(Layout *layout, int* permutation) {

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;


	//////////////////////////////////////////////// Permutation Generation /////////////////////////////

	vector<NFPriority> placed;
	vector<NFPriority> not_placed;

	for (int j = 0; j < nfcount; j++) {

		NFPriority priority;
		priority.id = j;
		priority.interaction = 0;

		for (int i = 0; i < efcount; i++)
			priority.interaction += efnfflow[i * nfcount + j];

		not_placed.push_back(priority);

	}

	NFPriorityCompare nfPriorityCompare;

	while (placed.size() != nfcount) {
		sort(not_placed.begin(), not_placed.begin() + not_placed.size(), nfPriorityCompare);
		NFPriority priority = not_placed.back();
		not_placed.pop_back();

		placed.push_back(priority);

		for (int i = 0; i < not_placed.size(); i++) {
			int nfid1 = priority.id;
			int nfid2 = not_placed.at(i).id;

			float extra_flow = nfnfflow[nfid1 * nfcount + nfid2];
			not_placed.at(i).interaction += extra_flow;
		}
	}

	for (int i = 0; i < nfcount; i++)
		permutation[i] = placed.at(i).id;

}

void lb_finite(Layout *layout, vector<int*> permutations, float* times, ObjectiveFunction *optimal, bool is_infinitesimal, bool rectilinear_approx) {

	cout << "QSAP_LB ";

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;
	int iocount = layout->iocount;

	float lwidth = layout->width;
	float lheight = layout->height;

	vector<float*> candidatelist;
	vector<Barrier> barrierlist(layout->barrierlist, layout->barrierlist + efcount);
	vector<IOPoint> iolist(layout->iolist, layout->iolist + iocount);
	vector<Barrier> nflist(layout->nflist, layout->nflist + nfcount);

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;
	float *efefflow2 = 0;
	float *efnfflow2 = 0;
	float nfnfflow2 = 0;

	///////////////////////////////////////////////// Candidate Generation //////////////////////////////////////////////


	float time1 = (float)omp_get_wtime();
	optimal->value = 0;
	optimal->candidatecount = 0;

	vector<int> placed;

	vector<float> vgridlines, hgridlines;
	getGridlinecoord(barrierlist, iolist, vgridlines, hgridlines, lwidth, lheight);

	int _hgridline_count = hgridlines.size();
	int _vgridline_count = vgridlines.size();
	int ycell_count = _hgridline_count - 1;
	int xcell_count = _vgridline_count - 1;

	Node *network = new Node[_hgridline_count * _vgridline_count];
	Cell *cells = new Cell[(xcell_count)*(ycell_count)]; //additional space for NF cells

	createNetwork(network, hgridlines, vgridlines);
	createCells(cells, network, hgridlines, vgridlines);
	getBarrierCellOverlap(network, cells, barrierlist, hgridlines, vgridlines);
	getIONodeOverlap(network, iolist, hgridlines, vgridlines);

	int M = nfcount;
	int N = (_hgridline_count * _vgridline_count);

	float *F = nfnfflow;

	float *D = new float[(_hgridline_count)*(_vgridline_count)*(_hgridline_count)*(_vgridline_count)];
	allPairSP(D, network, _hgridline_count, _vgridline_count);

	float *b = new float[nfcount * _hgridline_count * _vgridline_count]; // Fixed-cost matrix
	constructXCosts(b, efnfflow, D, iolist, iocount, nfcount, _hgridline_count, _vgridline_count);

	float efinteraction = 0;

	getEFInteraction(efinteraction, iolist, efefflow, D, N);

	float objval = 0;

	//solveQSAP(objval, b, F, D, M, N);

	optimal->value = efinteraction + objval;

	float time2 = (float)omp_get_wtime();

	times[0] = time2 - time1;

	cout << "done! (lb: " << optimal->value << " time: " << time2 - time1 << " s)" << endl;

	times[1] = time2 - time1;
	times[2] = time2 - time1;

	delete[] network;
	delete[] cells;

	delete[] b;
	delete[] D;

}

void updateFlowMatrices(float *efef2, float *efnf2, float *nfnf2, float *efef, float* efnf, float * nfnf, int efcount, int nfcount, vector<int> placed, vector<int> to_be_placed) {

	int newefcount = efcount + placed.size();
	int newnfcount = to_be_placed.size();

	// Update EF-EF flow matrix
	for (int i = 0; i < efcount; i++) {
		for (int j = 0; j < efcount; j++) {
			efef2[i * newefcount + j] = efef[i * efcount + j];
		}

		for (int j = 0; j < placed.size(); j++) {
			int jj = j + efcount;
			int nfid = placed.at(j);
			efef2[i * newefcount + jj] = efnf[i * nfcount + nfid];
			efef2[jj * newefcount + i] = efnf[i * nfcount + nfid];
		}
	}

	for (int i = 0; i < placed.size(); i++) {
		for (int j = 0; j < placed.size(); j++) {
			int ii = i + efcount;
			int jj = j + efcount;

			int nfid1 = placed.at(i);
			int nfid2 = placed.at(j);

			efef2[ii * newefcount + jj] = nfnf[nfid1 * nfcount + nfid2];
		}
	}


	// Update EF-NF flow matrix
	for (int i = 0; i < efcount; i++) {
		for (int j = 0; j < newnfcount; j++) {
			int nfid = to_be_placed.at(j);
			efnf2[i * newnfcount + j] = efnf[i * nfcount + nfid];
		}
	}

	for (int i = 0; i < placed.size(); i++) {
		int ii = i + efcount;
		int nfid1 = placed.at(i);

		for (int j = 0; j < newnfcount; j++) {
			int nfid2 = to_be_placed.at(j);
			efnf2[ii * newnfcount + j] = nfnf[nfid1 * nfcount + nfid2];
		}
	}

	// Update NF-NF flow matrix
	for (int i = 0; i < newnfcount; i++) {
		for (int j = 0; j < newnfcount; j++) {
			int nfid1 = to_be_placed.at(i);
			int nfid2 = to_be_placed.at(j);

			nfnf2[i * newnfcount + j] = nfnf[nfid1 * nfcount + nfid2];
		}
	}
}

/*
 *	This function is used for testing the corrected procedure for M Facility Placement
 */
 // Optimal Procedure
void executeComprehensiveOptimal(Layout *layout, vector<int*> permutations, float* times, ObjectiveFunction *optimal) {

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;
	int iocount = layout->iocount;

	float lwidth = layout->width;
	float lheight = layout->height;
	vector<Barrier> barrierlist(layout->barrierlist, layout->barrierlist + efcount);
	vector<IOPoint> iolist(layout->iolist, layout->iolist + iocount);
	vector<Barrier> nflist(layout->nflist, layout->nflist + nfcount);

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;

	///////////////////////////////////////////////// Candidate Generation //////////////////////////////////////////////

	cout << "OPT_ALL_PERM "; //(total perm: " << permutations.size() << ")" << endl;

	float time1 = (float)omp_get_wtime();

	vector<float> vgridlinecoord;
	vector<float> hgridlinecoord;

	getGridlinecoord(barrierlist, iolist, vgridlinecoord, hgridlinecoord, lwidth, lheight);

	vector<set<NFCandidate, CompareByPosition>> masterlist;
	vector<set<float, FloatCompare>> original_qx;
	vector<set<float, FloatCompare>> original_qy;

	char* alignment = new char[nfcount];

	for (int m = 0; m < nfcount; m++) {

		set<NFCandidate, CompareByPosition> candidatelist;
		set<float, FloatCompare> qx;
		set<float, FloatCompare> qy;

		Barrier nf = nflist.at(m);
		float nf_width = nf.getWidth();
		float nf_height = nf.getHeight();

		vector<float> nfcoord;

		getSingleNFCoordinates(nfcoord, barrierlist, iolist, nf_width, nf_height, lwidth, lheight);

		int candidatecount = nfcoord.size() / 2;

		for (int i = 0; i < candidatecount; i++) {
			float x = nfcoord.at(2 * i);
			float y = nfcoord.at(2 * i + 1);
			qx.insert(x);
			qy.insert(y);

			NFCandidate cnd;
			cnd.x = x;
			cnd.y = y;

			candidatelist.insert(cnd);

		}

		masterlist.push_back(candidatelist);
		original_qx.push_back(qx);
		original_qy.push_back(qy);
	}


	for (int k = 0; k < permutations.size(); k++) {

		bool inf_flag = false;
		float *candidate = new float[2 * nfcount];
		generateCandidates3(optimal, permutations.at(k), 0, nfcount, inf_flag, candidate, barrierlist, nflist, iolist, efefflow, efnfflow, nfnfflow, lwidth, lheight, masterlist, original_qx, original_qy);
		delete[] candidate;
	}

	float time2 = (float)omp_get_wtime();

	times[0] = time2 - time1;

	//////////////////////////////////////////////// Candidate Evaluation /////////////////////////////

//	cout << "C_EVAL (total cnds: " << candidatelist.size() << ")" << endl;

//	(*optimal).candidatecount = candidatelist.size();

/*	for (int i = 0; i < candidatelist.size(); i++) {

		float* coord = candidatelist.at(i);
		int *perm = permlist.at(i);

		evaluateCandidates(optimal, perm, efefflow, efnfflow, nfnfflow, coord, barrierlist, nflist, iolist, lwidth, lheight, false);

		delete[] candidatelist.at(i);
	}*/

	//////////////////////////////////////////////////////////////////////////////////////////////////

	float time3 = (float)omp_get_wtime();

	cout << "ceval done! (opt: " << (*optimal).value << " cnd: " << (*optimal).candidatecount << " time: " << time3 - time1 << " s)" << endl;

	times[1] = time3 - time2;
	times[2] = time3 - time1;

}

/*
 *	This is the true comprehensive optimal procedure written for IIE Transactions paper.
 */
void executeTrueComprehensiveOptimal(Layout *layout, vector<int*> permutations, float* times, ObjectiveFunction *optimal, float maxGap, int maxTime) {

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;
	int iocount = layout->iocount;

	float lwidth = layout->width;
	float lheight = layout->height;
	vector<Barrier> barrierlist(layout->barrierlist, layout->barrierlist + efcount);
	vector<IOPoint> iolist(layout->iolist, layout->iolist + iocount);
	vector<Barrier> nflist(layout->nflist, layout->nflist + nfcount);

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;

	vector<int> branchingOrder(permutations.data()[0], permutations.data()[0] + nfcount);

	BBController bbController(lwidth, lheight, barrierlist, nflist, iolist, branchingOrder, efefflow, efnfflow, nfnfflow);
	bbController.setBranchingStrategy(BranchingStrategy::BstFS);
	bbController.setMaxGap(maxGap);
	bbController.setMaxTime(maxTime);

	///////////////////////////////////////////////// Candidate Generation //////////////////////////////////////////////

	cout << "OPT_ALL_PERM" << endl; //(total perm: " << permutations.size() << ")" << endl;

	float time1 = (float)omp_get_wtime();

	//candidateRepository.GenerateNFCandidates();
	//candidateRepository.EvaluateNFCandidates(optimal, efefflow, efnfflow, nfnfflow);

	bbController.run(optimal);

	float time2 = (float)omp_get_wtime();

	times[2] = time2 - time1;

	cout << "ceval done! (opt: " << (*optimal).value << " bound: " << optimal->bestbound << " nodes: " << optimal->candidatecount << " feasible: " << optimal->feasibleCandidates << " time: " << time2 - time1 << " s)\n" << endl;
}

void getUniqueCombinations(vector<vector<int>> &combinations, vector<int> &combination, int M, int m, int depth, int ii, int &count) {

	if (depth == m) {
		vector<int> newcombination(combination.begin(), combination.end());
		newcombination.push_back(count++);
		combinations.push_back(newcombination);
	}

	else {
		for (int i = ii; i < M; i++) {
			combination.push_back(i);
			getUniqueCombinations(combinations, combination, M, m, depth + 1, i + 1, count);
			combination.pop_back();
		}
	}
}

void getUniquePermutations(vector<int*> &permutations, int M, int m) {
	vector<vector<long long int>> keycheck;
	vector<int*> newperm;
	int size = M / m;
	int over_m = M % m;
	int extra_size = over_m > 0 ? 1 : 0;

	for (vector<int*>::iterator it = permutations.begin(); it != permutations.end(); ++it) {

		vector<long long int> keylist;

		int* ptr = *it;

		for (int i = 0; i < size; i++) {
			set<int> keyset;
			for (int j = 0; j < m; j++)
				keyset.insert(ptr[i * m + j]);

			long long int key = integer_keygen(keyset, M);
			keylist.push_back(key);
		}

		if (extra_size) {
			set<int> keyset;
			for (int j = 0; j < over_m; j++)
				keyset.insert(ptr[size * m + j]);

			long long int key = integer_keygen(keyset, M);
			keylist.push_back(key);
		}

		bool flag = false;
		for (vector<vector<long long int>>::iterator it1 = keycheck.begin(); it1 != keycheck.end(); ++it1)
			flag = flag || vector_isequal(keylist, *it1);

		if (!flag) {
			newperm.push_back(ptr);
			keycheck.push_back(keylist);
		}
		else {
			delete[] ptr;
			ptr = 0;
		}
	}

	permutations = newperm;

}

long long int integer_keygen(set<int> digits, int base) {

	long long int key = 0;

	int size = digits.size() - 1;

	for (set<int>::iterator it = digits.begin(); it != digits.end(); ++it) {
		long long int digit = *it;
		if (digit >= base) {
			cerr << "Digit cannot be greater than base" << endl;
			exit(-1);
		}

		long long int mul = 1;
		for (int i = 0; i < size; i++)
			mul *= base;
		key += digit * mul;
		size--;
	}

	return key;
}

bool vector_isequal(vector<long long int> vector1, vector<long long int> vector2) {

	if (vector1.size() != vector2.size())
		return false;

	int size = vector1.size();

	for (int i = 0; i < size; i++)
		if (vector1.at(i) != vector2.at(i))
			return false;

	return true;
}

// This function checks the objective value of a solution for comparison.
void verifyObjectiveValue(Layout *layout, ObjectiveFunction *optimal) {

	int nfcount = layout->nfcount;
	int efcount = layout->barriercount;
	int iocount = layout->iocount;

	float lwidth = layout->width;
	float lheight = layout->height;

	float *efefflow = layout->efef;
	float *efnfflow = layout->efnf;
	float *nfnfflow = layout->nfnf;

	vector<Barrier> barrierlist(layout->barrierlist, layout->barrierlist + efcount);
	vector<IOPoint> iolist(layout->iolist, layout->iolist + iocount);
	vector<Barrier> nflist(layout->nflist, layout->nflist + nfcount);


	vector<int> to_be_placed;
	vector<int> placed;

	for (int i = 0; i < nfcount; i++)
		placed.push_back(i);

	for (int k = 0; k < placed.size(); k++) {
		int nfid = placed.at(k);
		Barrier newb;
		newb.setBarrier(optimal->ptr[2 * nfid], optimal->ptr[2 * nfid + 1], nflist.at(nfid).getWidth(), nflist.at(nfid).getHeight(), nflist.at(nfid).getCongestion());
		barrierlist.push_back(newb);

		IOPoint newio;
		newio.x = optimal->ptr[2 * nfid];
		newio.y = optimal->ptr[2 * nfid + 1];
		newio.isNF = false;
		iolist.push_back(newio);
	}

	vector<Barrier> newnflist;

	float *new_efefflowmatrix = new float[(efcount + placed.size()) * (efcount + placed.size())];
	float *new_efnfflowmatrix = new float[(efcount + placed.size()) * to_be_placed.size()];
	float *new_nfnfflowmatrix = new float[to_be_placed.size() * to_be_placed.size()];

	updateFlowMatrices(new_efefflowmatrix, new_efnfflowmatrix, new_nfnfflowmatrix, efefflow, efnfflow, nfnfflow, efcount, nfcount, placed, to_be_placed);

	ObjectiveFunction lcl_opt2;
	lcl_opt2.ptr = new float[2 * nfcount];
	fill(lcl_opt2.ptr, lcl_opt2.ptr + 2 * nfcount, -1.0);
	lcl_opt2.value = FLT_MAX;
	lcl_opt2.candidatecount = 0;
	lcl_opt2.perm = new int[nfcount];

	evaluateCandidates(&lcl_opt2, optimal->perm, new_efefflowmatrix, new_efnfflowmatrix, new_nfnfflowmatrix, optimal->ptr, barrierlist, newnflist, iolist, lwidth, lheight, false);

	cout << "Test Value: " << lcl_opt2.value << endl;

	delete[] new_efefflowmatrix;
	delete[] new_efnfflowmatrix;
	delete[] new_nfnfflowmatrix;

	delete[] lcl_opt2.ptr;
	delete[] lcl_opt2.perm;

}