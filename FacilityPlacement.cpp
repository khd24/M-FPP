// FacilityPlacement.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>

#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <queue>
#include <float.h>
#include <ctime>
#include "structures_1.h"
#include "structures_2.h"
#include "functions_1.h"
#include "functions_2.h"
#include "functions_3.h"
#include "functions_4.h"
#include "functions_alg.h"
#include "functions_u.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <cctype>
#include <math.h>

using namespace std;

//#define EXTERNAL /* This parameter is used to switch on/off the command line input */

#define MAX_TIME 0;
#define MAX_GAP 0.0f;

int main(int argc, char** argv) {

	// convert now to string form

	float maxGap = MAX_GAP;
	int maxTime = MAX_TIME;
	string extension = ".txt";

	std::string in_path, out_path, layout_path, repair_path;
	int problemsize;

	ofstream myfile, layoutfile;

	int rpr_itn = 15;

#ifndef EXTERNAL

	in_path = "problemset_4_3";
	out_path = "output_4_3";
	layout_path = "layout_4_3";
	int mode = 0;
	int bucketsize = 1;
	bool repair = 0;
	int rpr_bucketsize = 2;


#else
	in_path = argv[1];
	out_path = argv[2];
	layout_path = argv[3];
	int mode = atoi(argv[4]);
	int bucketsize = atoi(argv[5]);
	int repair = atoi(argv[6]);
	int rpr_bucketsize = atoi(argv[7]);


	if (argc != 8) {
		std::cout << "Invalid number of arguments!" << std::endl;
		std::cout << "Usage: arg[1]: input file, arg[2]: log file, arg[3]: layout file, arg[4]: mode, arg[5]: heu bucketsize, arg[6]: repair mode, arg[7]: repair bucketsize" << std::endl;
		exit(-1);
	}
#endif

	if (mode < 0 || mode > 6) {
		std::cerr << "Invalid mode! Should be a number between 0 to 6." << std::endl;
		std::cerr << "Usage: <0>: Optimal, <1>: Heuristic 1, <2>: Heuristic 2, <3>: Heuristic 3, <4>: Heuristic 4, <5>: Heuristic 5, <6>: QSAP LB" << std::endl;
		exit(-1);
	}

	if (repair < 0 || repair > 1) {
		std::cerr << "Invalid repair argument! Should be 0 or 1." << std::endl;
		std::cerr << "Usage: <0>: Repair off, <1>: Repair on" << std::endl;
		exit(-1);
	}

	if (bucketsize > 4 || rpr_bucketsize > 4)
		std::cout << "Large bucket size selected. Solver may take long time." << std::endl;


	////////////////////////// Read problem and generate permutations //////////////////////////

	Layout* problemset = readLayoutFromFile(in_path.append(extension), problemsize);

	int nfcount = problemset[0].nfcount;
	int efcount = problemset[0].barriercount;
	int iocount = problemset[0].iocount;

	float lwidth = problemset[0].width;
	float lheight = problemset[0].height;

	if (nfcount > 10 && mode % 2 == 0)
		std::cout << "Number of permutations will be too large. Solver may take long time." << std::endl;

	std::vector<int*> permutations;

	if (mode % 2 == 0) {
		int *a = new int[nfcount];
		for (int i = 0; i < nfcount; i++)
			a[i] = i;

		permute(permutations, a, 0, nfcount);

		delete[] a;
	}

	///////////////////////////////////////////////////////////////////////////////////////////

	char timestamp[16];
	time_t caltime;
	struct tm* broketime;

	// find current time, convert to broken-down time
	time(&caltime);
	broketime = localtime(&caltime);

	// append timestamp in the format "_yymmdd_hhmmss"
	strftime(timestamp, 16, "_%m%d%y_%H%M%S", broketime);

	if (bucketsize > nfcount)
		bucketsize = nfcount;
	if (bucketsize < 1)
		bucketsize = 1;

	if (rpr_bucketsize > nfcount)
		rpr_bucketsize = nfcount;
	if (rpr_bucketsize < 1)
		rpr_bucketsize = 1;

	std::stringstream sstream1;
	sstream1 << "h" << mode;
	if (mode > 0 && mode < 6) {
		if (mode != 1)
			sstream1 << "_m" << bucketsize;
		if (repair)
			sstream1 << "_rm" << rpr_bucketsize;
	}
	sstream1 << "_" << out_path.c_str();
	out_path = sstream1.str();

#ifdef EXTERNAL

	out_path.append(timestamp);

#endif

	std::stringstream sstream2;
	sstream2 << "h" << mode;
	if (mode > 0 && mode < 6) {
		if (mode != 1)
			sstream2 << "_m" << bucketsize;
		if (repair)
			sstream2 << "_rm" << rpr_bucketsize;
	}
	sstream2 << "_" << layout_path.c_str();
	layout_path = sstream2.str();

#ifdef EXTERNAL

	layout_path.append(timestamp);

#endif

	myfile.open(out_path.append(extension).c_str());
	myfile << "Problem size: " << problemsize << " | Layout size: " << problemset->width << " X " << problemset->height << " | No. of EFs: " << problemset->barriercount << " | No. of NFs: " << problemset->nfcount << " | Max Gap: " << maxGap << " | Max Time: " << maxTime << "\n" << endl;
	if (mode < 6) {
		if (repair && mode > 0)
			myfile << "Pr. No.\tTotal Permutations\tOptimal Permutation\tFeasible Candidates\tOptimal Placement\tObjective Value\tTime (s)\tMode\tImproved Objective Value\tTime(s)" << endl;
		else
			myfile << "Pr. No.\tTotal Permutations\tOptimal Permutation\tNodes Processed\tFeasible Candidates\tOptimal Placement\tObjective Value\tBest Lower Bound\tTime (s)\tMode" << endl;
	}
	else
		myfile << "Pr. No.\tObjective Value\tTime (s)\tMode" << endl;

	myfile.flush();

	if (mode < 6) {
		layoutfile.open(layout_path.append(extension).c_str());
		layoutfile << problemsize << std::endl;
		layoutfile << lwidth << std::endl;
		layoutfile << lheight << std::endl;
		layoutfile << efcount << std::endl;
		layoutfile << nfcount << std::endl;
		layoutfile.close();
	}

	for (int pn = 0; pn < problemsize; pn++)
	{

		time_t now = time(0);

		char* dt = ctime(&now);

		dt[strlen(dt) - 1] = 0;

		Layout* layout = problemset + pn;

		std::cout << "PROB # " << pn + 1 << "\tTime: " << dt << std::endl;

		////////////////////////////////// OPTIMAL ALGORITHM /////////////////////////////////////////////

		float etime_gbl_opt[3];
		ObjectiveFunction gbl_opt;


		gbl_opt.ptr = new float[2 * nfcount];
		std::fill(gbl_opt.ptr, gbl_opt.ptr + 2 * nfcount, -1.0);
		gbl_opt.value = FLT_MAX;
		gbl_opt.candidatecount = 0;
		gbl_opt.feasibleCandidates = 0;
		gbl_opt.perm = new int[nfcount];
		std::fill(gbl_opt.perm, gbl_opt.perm + nfcount, -1);

		/////////////////// GLOBAL OPTIMAL /////////////////////////////////

		int *permutation_aldep = new int[nfcount];

		if (mode == 0)
			//executeComprehensiveOptimal(layout, permutations, etime_gbl_opt, &gbl_opt);
		//			executeCornerOptimalHeuristic_allPerm(layout, permutations, etime_gbl_opt, &gbl_opt);

			executeTrueComprehensiveOptimal(layout, permutations, etime_gbl_opt, &gbl_opt, maxGap, maxTime);

				////////////////////////////////////////////////////////////////

				//////////////////////// GIVEN PERM ////////////////////////////////////

		/*		if (mode == 1) {
					std::cout << "Enter a permutation" << std::endl;
					std::vector<int> perm;
					while (perm.size() != nfcount) {
						int p;
						cin >> p;
						std::vector<int>::iterator it;
						it = find(perm.begin(), perm.end(), p);
						if (it == perm.end() && p > 0 && p <= nfcount)
							perm.push_back(p);
						else
							std::cout << "Enter valid facility number" << std::endl;
					}
					for (int i = 0; i < perm.size(); i++)
						permutation_aldep[i] = perm.at(i) - 1;

					executeCornerOptimalHeuristic_givenPerm(layout, permutation_aldep, etime_gbl_opt, &gbl_opt);
				}*/

				//////////////////////// ALDEP ////////////////////////////////////

		if (mode == 1) {
			generatePermutation_ALDEP(layout, permutation_aldep);
			executeCornerOptimalHeuristic_givenPerm(layout, permutation_aldep, etime_gbl_opt, &gbl_opt);

		}

		////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////

		if (mode == 2)
			executeClass1(layout, permutations, etime_gbl_opt, &gbl_opt, bucketsize);

		////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////

		if (mode == 3) {
			generatePermutation_ALDEP(layout, permutation_aldep);
			std::vector<int*> perm;
			perm.push_back(permutation_aldep);
			executeClass1(layout, perm, etime_gbl_opt, &gbl_opt, bucketsize);
		}

		////////////////////////////////////////////////////////////////

		if (mode == 4)
			executeClass2(layout, permutations, etime_gbl_opt, &gbl_opt, bucketsize);

		////////////////////////////////////////////////////////////////

		if (mode == 5) {
			generatePermutation_ALDEP(layout, permutation_aldep);
			std::vector<int*> perm;
			perm.push_back(permutation_aldep);
			executeClass2(layout, perm, etime_gbl_opt, &gbl_opt, bucketsize);
		}

		////////////////////////////////////////////////////////////////

		//////////////////////// LB CURRENTLY NOT IN USE ////////////////////////////////////

		if (mode == 6)
			lb_finite(layout, permutations, etime_gbl_opt, &gbl_opt, false, false);

		////////////////////////////////////////////////////////////////



		//////////////////////////// PRINT RESULTS ////////////////////////////////////

		if (mode < 6) {
			myfile << pn + 1 << "\t" << ((mode % 2 == 0) ? permutations.size() : 1);

			if (gbl_opt.candidatecount > 0) {

				std::stringstream s;
				if (mode == 0)
					s << "OPT";
				else if (mode == 1)
					s << "GIVEN PERM";
				else if (mode == 2)
					s << "HEU 2 BUCKET " << bucketsize;
				else if (mode == 3)
					s << "HEU 3 BUCKET " << bucketsize;
				else if (mode == 4)
					s << "HEU 4 BUCKET " << bucketsize;
				else
					s << "HEU 5 BUCKET " << bucketsize;

				myfile << "\t[";

				for (int nfid = 0; nfid < nfcount - 1; nfid++)
					myfile << gbl_opt.perm[nfid] + 1 << ", ";

				myfile << gbl_opt.perm[nfcount - 1] + 1 << "]\t" << gbl_opt.candidatecount << "\t" << gbl_opt.feasibleCandidates << "\t";

				for (int nfid = 0; nfid < nfcount; nfid++)
					myfile << "(" << gbl_opt.ptr[2 * nfid] << ", " << gbl_opt.ptr[2 * nfid + 1] << ")";

				myfile << "\t" << gbl_opt.value << "\t" << gbl_opt.bestbound << "\t" << etime_gbl_opt[2];

				if (repair &&  mode > 0 && mode < 6) {
					executeRepair(layout, etime_gbl_opt, &gbl_opt, rpr_bucketsize, rpr_itn);
					s << " W/ REPAIR BUCKET " << rpr_bucketsize << "\t" << gbl_opt.value << "\t" << etime_gbl_opt[2];
				}

				myfile << "\t" << s.str();

			}

			myfile << std::endl;

			myfile.flush();


			//////////////////////////////////////////////////////////////////////////////

			for (int nfid = 0; nfid < layout->nfcount; nfid++) {
				Barrier* bptr = &(layout->nflist[nfid]);
				bptr->setX(gbl_opt.ptr[2 * nfid]);
				bptr->setY(gbl_opt.ptr[2 * nfid + 1]);
			}

			writeLayoutToFile(layout_path, layout->barriercount, layout->nfcount, layout->barrierlist, layout->nflist, layout->iolist);
		}

		else {
			myfile << pn + 1 << "\t" << gbl_opt.value << "\t" << etime_gbl_opt[2] << "\tLB";
			myfile << std::endl;
			myfile.flush();
		}


		//////////////////////////////////////////////////////////////////////////////


		delete[] permutation_aldep;

		delete[] gbl_opt.ptr;
		delete[] gbl_opt.perm;

		delete[] layout->barrierlist;
		delete[] layout->iolist;
		delete[] layout->efef;
		delete[] layout->efnf;
		delete[] layout->nfnf;

	}

	for (int i = 0; i < permutations.size(); i++)
		delete[] permutations.at(i);

	myfile.close();

	//	MPI_Finalize();

	return 0;
}