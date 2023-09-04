/*
 * functions_lb.cpp
 *
 *  Created on: Jun 13, 2020
 *      Author: ketandat, aanandan
 */

#include "functions_lb.h"

void allPairSP(float *shortestpaths, Node *network, int _hgridline_count, int _vgridline_count)
{
	for (int i = 0; i < (_hgridline_count)*(_vgridline_count); i++) {

		float *dijdist = &shortestpaths[i * (_hgridline_count)*(_vgridline_count)];
		fill_n(dijdist, (_hgridline_count)*(_vgridline_count), FLT_MAX);

		executeDijkstra(network, network + i, dijdist, (_hgridline_count), (_vgridline_count));
	}
}

void constructXCosts(float * xcosts, float *efnfflow, float* shortestpaths, vector<IOPoint> iolist, int iocount, int nfcount, int _hgridlinecount, int _vgridlinecount) {

	float *nfefflow = new float[nfcount * iocount];
	for (int i = 0; i < iocount; i++)
		for (int j = 0; j < nfcount; j++)
			nfefflow[j * iocount + i] = efnfflow[i * nfcount + j];

	float *dist = new float[iocount * (_hgridlinecount *_vgridlinecount)];

	for (int i = 0; i < iocount; i++) {
		int io = iolist.at(i).nodeid;
		copy(&shortestpaths[io * (_hgridlinecount *_vgridlinecount)], &shortestpaths[io * (_hgridlinecount *_vgridlinecount)] + (_hgridlinecount *_vgridlinecount), &dist[i * (_hgridlinecount *_vgridlinecount)]);
	}


	for (int m = 0; m < nfcount; m++) {

		for (int n = 0; n < (_hgridlinecount *_vgridlinecount); n++) {

			float sum = 0;

			for (int i = 0; i < iocount; i++) {

				float flow = nfefflow[m * iocount + i];
				float distance = dist[i * (_hgridlinecount *_vgridlinecount) + n];

				if (float_lt(distance, FLT_MAX)) {
					sum += flow * distance;
				}

				else {
					sum = FLT_MAX;
					break;
				}
			}

			xcosts[m * (_hgridlinecount *_vgridlinecount) + n] = sum;
		}
	}

	delete[] nfefflow;
	delete[] dist;
}

void getEFInteraction(float &efinteraction, vector<IOPoint> iolist, float *efefflow, float* shortestpaths, int nodecount) {

	efinteraction = 0;

	int iocount = iolist.size();
	for (int i = 0; i < iocount; i++) {
		int nodeidI = iolist.at(i).nodeid;

		for (int j = i + 1; j < iocount; j++) {
			int nodeidJ = iolist.at(j).nodeid;

			float flow = efefflow[i * iocount + j];
			float dist = shortestpaths[nodeidI * nodecount + nodeidJ];

			efinteraction += flow * dist;
		}
	}
}

//void solveQSAP(float &objval, float *b, float *F, float *D, int M, int N) {
//
//	try
//    {
//        GRBEnv env = GRBEnv();
//		env.set(GRB_IntParam_OutputFlag, 0);
//        GRBModel model = GRBModel(env);
// 
//        // Set variables
//
//		GRBVar *x = new GRBVar[M * N];
//		GRBVar *y = new GRBVar[M * M * N * N];
//
//		for (int i = 0; i < M; i++) {
//			for (int p = 0; p < N; p++) {
//				float cost = b[i * N + p];
//				stringstream s;
//				s << "X_" << i << "_" << p;
//				x[i * N + p] = model.addVar(0.0, 1.0, cost, GRB_BINARY, s.str());
//			}
//		}
//
//		for (int i = 0; i < M; i++) {
//			for (int p = 0; p < N; p++) {
//				for (int j = i + 1; j < M; j++) {
//					for (int q = 0; q < N; q++) {
//						float cost = float_lt(D[p * N + q], FLT_MAX) ? F[i * M + j] * D[p * N + q] : FLT_MAX;
//						stringstream s;
//						s << "Y_" << i << "_" << j << "_" << p << "_" << q;
//						y[i * N * M * N + p * M * N + j * N + q] = model.addVar(0.0, 1.0, cost, GRB_CONTINUOUS, s.str());	
//					}	
//				}
//			}
//		}
//      
//        model.update();
//
//		// X LAP constraints
//
//		for (int i = 0; i < M; i++) {
//			GRBLinExpr lhs = 0;
//			for (int p = 0; p < N; p++)
//				lhs += x[i * N + p];
//
//			stringstream s;
//			s << "XR_" << i;
//			model.addConstr(lhs == 1, s.str());
//		}
//
//		// Y LAP constraints
//
//		for (int i = 0; i < M; i++) {
//			for (int p = 0; p < N; p++) {
//				for (int j = 0; j < M; j++) {
//					if (j != i) {
//						GRBLinExpr lhs = 0;
//						for (int q = 0; q < N; q++){
//							if(i < j) 
//								lhs += y[i * N * M * N + p * M * N + j * N + q];
//							else
//								lhs += y[j * N * M * N + q * M * N + i * N + p];	
//						}
//						
//						stringstream s;
//						s << "YR_" << i << "_" << j << "_" << p;
//						model.addConstr(lhs == x[i * N + p], s.str());
//					}
//				}
//			}
//		}
//
//		model.update();
// 
//        // Run the optimizion
//        model.optimize();
//
//		objval = (float) model.get(GRB_DoubleAttr_ObjVal);
// 
//  //      cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
//
//		delete[] x;
//		delete[] y;
// 
//    }
//    catch (GRBException e) 
//    {
//        cout << "Error code = "
//                  << e.getErrorCode() 
//                  << endl;
//        cout << e.getMessage() << endl;
//    }
//    catch (...)
//    {
//        cout << "Exception during optimization"
//                  << endl;
//    }
//}