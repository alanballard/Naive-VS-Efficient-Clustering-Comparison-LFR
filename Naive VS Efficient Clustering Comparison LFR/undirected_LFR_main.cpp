#include <iostream>
//#include "Network.h"
using namespace std;

void srand_file(void);
//int benchmark(int **arr, int *comm, bool excess, bool defect, int nodenum, double avgdeg, int maxdeg, double expdeg, double comsiz, double mixpar, int rep);
int benchmark(bool excess, bool defect, int nodenum, double avgdeg, int maxdeg, double expdeg, double comsiz, double mixpar, int rep, int numstart);


void main()
{
	srand_file();

/*
To create:
N=250, 500, 750, 1000, 1250, 1500, 1750, 2000
Values from Fig. 9...
avgdeg=10
maxdeg=25
GAMMA expdeg=2 
BETA comsize=1
MYU mixpar=0.1
*/
	/*LFR Graph Generator Parameters*/
	bool excess = false;
	bool defect = false;
	int nodenum = 100;// 1000;
	double min_avgdeg = 10;// 50;// (double)nodenum / 20;
	double max_avgdeg = 10;// 100;// (double)nodenum / 10;
	double step_avgdeg = 1;// 25;//(double)nodenum / 60;
	//double avgdeg = (double)nodenum/20;// nodenum*.02;
	int min_maxdeg = 25;// 200;// nodenum / 5;
	int max_maxdeg = 25;// 00;// nodenum / 5;
	int step_maxdeg = 1;

//	int maxdeg = 100;

	//int maxdeg = nodenum/5;// nodenum*.05;
	double min_expdeg = 2;
	double max_expdeg = 2;
	double step_expdeg = 1;
	//double expdeg = 2;
	double min_comsize = 1;
	double max_comsize = 1;
	double step_comsize = 1;
	//double comsiz = 1;
	double min_mixpar = 0.1;
	double max_mixpar = 0.1;
	double step_mixpar = 0.1;
	//double mixpar = 0.1;
	int rep_nbr = 1;
	int numstart = 0; //Vertices will be numbered started with this value. Should be {0,1}. Community numbers will still start at 1.
		
	int numberofgroups = 0;		//Used to store the true number of groups as determined by LFR
	for (int rep = 1; rep <= rep_nbr; rep++){
		for (double avgdeg = min_avgdeg; avgdeg <= max_avgdeg; avgdeg=avgdeg + step_avgdeg){
			for (int maxdeg = min_maxdeg; maxdeg <= max_maxdeg; maxdeg=maxdeg + step_maxdeg){
				for (double expdeg = min_expdeg; expdeg <= max_expdeg; expdeg=expdeg + step_expdeg){
					for (double comsize = min_comsize; comsize <= max_comsize; comsize=comsize + step_comsize){
						for (double mixpar = min_mixpar; mixpar <= max_mixpar; mixpar=mixpar + step_mixpar){
							//cout << avgdeg << endl << maxdeg << endl << expdeg << endl << comsize << endl << mixpar << endl;
							//cin.get(); cin.get();
							numberofgroups = benchmark(excess, defect, nodenum, avgdeg, maxdeg, expdeg, comsize, mixpar, rep, numstart);
						}
					}
				}
			}
		}
	}
	/*Create LFR graph*/
	//here I am using the TRUE number of assigned groups as number of groups
	//It can be used as an input for testing later
	//numberofgroups = benchmark(Adj, comm, excess, defect, nodenum, avgdeg, maxdeg, expdeg, comsiz, mixpar, rep); //returns true number of groups
	
	//numberofgroups = benchmark(excess, defect, nodenum, avgdeg, maxdeg, expdeg, comsiz, mixpar, rep, numstart); //returns true number of groups
	
}