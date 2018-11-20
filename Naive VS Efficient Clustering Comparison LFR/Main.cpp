#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>	
#include <time.h>		//for calculating run-time
#include <vector>
#include <cctype>		//for yes/no question
#include <iomanip>

void srand_file(void);
int benchmark(bool excess, bool defect, int nodenum, double avgdeg, int maxdeg, double expdeg, double comsiz, double mixpar, int rep, int network_key);


using namespace std;

//FUNCTION PROTOTYPES
vector< vector<pair<int, int>> > read_el(char *filename); //Read in edge list and convert to adjacency list
void print_el(vector< vector<pair<int, int>> > links); //Print adjacency and/or edge list

//Each method will be fed the same starting network and other paramters
//Naive likelihood 
int lik_naive(
	int network_key,
	vector< vector<pair<int, int>> > links,
	int nb_links,
	int N,
	int min_k,
	int max_k,
	int k_int,
	double InitTemp,
	double CR,
	double TL,
	int Max_Success,
	double LimitIT,
	double epsilon,
	int TL_counter,
	int Success_counter,
	double IT
	);
//Iterative/efficient likelihood
int lik_efficient(
	int network_key,
	vector< vector<pair<int, int>> > links,
	int nb_links,
	int N,
	int min_k,
	int max_k,
	int k_int,
	double InitTemp,
	double CR,
	double TL,
	int Max_Success,
	double LimitIT,
	double epsilon,
	int TL_counter,
	int Success_counter,
	double IT
	);
 //Modularity
 int modularity(
	int network_key,
	vector< vector<pair<int, int>> > links,
	int nb_links,
	int N,
	int min_k,
	int max_k,
	int k_int,
	double InitTemp,
	double CR,
	double TL,
	int Max_Success,
	double LimitIT,
	double epsilon,
	int TL_counter,
	int Success_counter,
	double IT
);
 int louvain_convert_2b(vector<vector<pair<int, int/*long double*/> > > source_links, int nb_links, char *filename);
 vector<int> louvain_clustering(char *filename, int network_key, vector< vector<pair<int, int>> > links, int nb_links);

int main()
{
	clock_t start_new = clock();

//STEP: Create LFR graph(s)
	bool NEW = true;
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

/*LFR Graph Generator Parameters. */

	bool excess = false; //These parameters should not be changed.
	bool defect = false; //These parameters should not be changed.

	//Number of nodes in the network
	int nodenum = 100;

	//average degree
	double min_avgdeg = 5;
	double max_avgdeg = 5;//10
	double step_avgdeg = 5;

	//maximum degree
	int min_maxdeg = 20;
	int max_maxdeg = 20;
	int step_maxdeg = 1;

	//degree distribution power law exponent
	double min_expdeg = 2;
	double max_expdeg = 2;//3
	double step_expdeg = 1;

	//community size distribution power law exponent
	double min_comsize = 1;
	double max_comsize = 1;//2
	double step_comsize = 1;

	//mixing parameter
	double min_mixpar = 0.1;
	double max_mixpar = 0.6;//.6
	double step_mixpar = 0.1;

	int rep_nbr = 1;//10

//int minnumber = 0;
//int maxnumber = 0;
//int numberflag = 1;

	int numberofgroups = 0;		//Used to store the true number of groups as determined by LFR

	int network_key = 1;
	for (int rep = 1; rep <= rep_nbr; rep++) {
		for (double avgdeg = min_avgdeg; avgdeg <= max_avgdeg; avgdeg = avgdeg + step_avgdeg) {
			for (int maxdeg = min_maxdeg; maxdeg <= max_maxdeg; maxdeg = maxdeg + step_maxdeg) {
				for (double expdeg = min_expdeg; expdeg <= max_expdeg; expdeg = expdeg + step_expdeg) {
					for (double comsize = min_comsize; comsize <= max_comsize; comsize = comsize + step_comsize) {
						for (double mixpar = min_mixpar; mixpar <= max_mixpar; mixpar = mixpar + step_mixpar) {
							numberofgroups = 0;
							numberofgroups = benchmark(excess, defect, nodenum, avgdeg, maxdeg, expdeg, comsize, mixpar, rep, network_key);
//cout<<"numberofgroups is "<<numberofgroups<<endl;
							//Just FYI. numberofgroups is the true number of clusters as determined by the LFR network generator
							//Can be used for testing, if desired.
//if (numberflag == 1) { minnumber = numberofgroups; maxnumber = numberofgroups; numberflag = 0; }
//else {
//		if (numberofgroups < minnumber) { minnumber = numberofgroups; }
//		if (numberofgroups > maxnumber) { maxnumber = numberofgroups; }
//}

//STEP: Write parameters of current LFR network to file
							cout << "Currently evaluating network " << network_key << " with parameters: " << avgdeg << ", " << maxdeg << ", " << expdeg << ", " << comsize << ", " << mixpar << ", " << rep << endl;
							
							ostringstream PARAMETERS_NAME;
							ofstream NETWORK_PARAMETERS;
							PARAMETERS_NAME << "LFR_PARAMETERS" << ".dat";
							NETWORK_PARAMETERS.open("LFR_PARAMETERS.dat", ios::app);
							NETWORK_PARAMETERS << network_key << "," << avgdeg << ", " << maxdeg << ", " << expdeg << ", " << comsize << ", " << mixpar << ", " << rep;
							NETWORK_PARAMETERS << endl;
							NETWORK_PARAMETERS.close();
//STEP: Read current LFR network
							char filename[256];
							strcpy_s(filename, "undirected_LFR_network.dat");
							ifstream ifile(filename);


	//cout << "Reading in edge list..." << endl;
	//Input network edge list and convert to adjacency list that includes weight information
	//time_t read_start = time(0);//Time until solution
	vector< vector<pair<int, int>> > links;
	links = read_el(filename);

	//Check to see that every vertex in the network has at least 1 edge to another vertex. No singletons allowed.
	int missing = 0;
	for (int i = 0; i < links.size(); i++)
	{
		if (links[i].size() == 0){ missing = 1; }
		//links[i].size()
	}
	if (missing == 1){ 
		cout << "WARNING WARNING WARNING" << endl;
		cout << "*** ADJACENCY LIST IS NOT FULLY POPULATED ***" << endl;
		cout << "THE PROGRAM WILL NOT WORK CORRECTLY" << endl;
		cin.get();
	}

	//Get the number of distinct vertices in the network
	int N = links.size();

	//cout << "Vertex count: " << links.size() << endl;
	//time_t read_end = time(0);
	//double read_time = difftime(read_end, read_start);
	//cout << read_time / 60 << " minutes required to read edge list" << endl;

	//Get the number of DISTINCT edges in the network
	int nb_links = 0;
	for (int i = 0; i < links.size(); i++){
		for (int j = 0; j < links[i].size(); j++)	{
			nb_links = nb_links + 1;
		}
	}
	//This counts the link between i and j, and between j and i. In the undirected network case using an adjacency list, this double counts edges, so:
	nb_links = nb_links / 2;

	//Print network as edge list and adjacency list. Can be commented out when running the real program.
	//print_el(links); //print network as adjacency list, including weights. Can also print out as edge list, but currently commented out

	//Specify the number of clusters in which to cluster the network
	int min_k = numberofgroups;	//minimum number of clusters to consider. Must be >=2
	int max_k = numberofgroups;	//maximum number of clusters to consider. Must be small enough that each cluster can contain atlest 2 vertices
	int k_int =1;   //step size from min_k to max_k. 1=evaluate every cluster size from min_k to max_k. 2=evaluate every other cluster size, and so on.
	//100
		//4-32, 6-22, 5-27, 6-28, 6-24, 6-26, 6-22, 6-25
	//1000
		//75-257, 77-254, 76-252

	//Simulated Annealing Algorithm Parameters
	double InitTemp = 2.5;
	double CR = .9925;
	int TL = 1000;			//Maximum number of reclustering attempts at a given temperature
	int Max_Success = TL;	//Maximum number of successes allowed at a given temperature
							//Note: Setting both these values=1 is the equivalent of not using them at all in the SA algorithm.
	double LimitIT = 1.0e-8;
	double epsilon = 1.0e-6;
	int TL_counter = 0;
	int Success_counter = 0;

	double IT = InitTemp; //IT changes during Simulated Annealing process. InitTemp is used to initialize (or re-initialize) IT at the start of the process when multiple k values are analyzed

//	int seed = 0;	
	srand(time(0));
	//		random_device rd;
	//		mt19937 e2(rd());
	//		uniform_real_distribution<double> dist_uni(0.0, 1.0);

	//cout << "#############################################################" << endl;
	 modularity(
		network_key,
		links,
		nb_links,
		N,
		min_k,
		max_k,
		k_int,
		InitTemp,
		CR,
		TL,
		Max_Success,
		LimitIT,
		epsilon,
		TL_counter,
		Success_counter,
		IT
	);

//cout << "#############################################################" << endl;
	lik_naive(
		network_key,
		links,
		nb_links,
		N,
		min_k,
		max_k,
		k_int,
		InitTemp,
		CR,
		TL,
		Max_Success,
		LimitIT,
		epsilon,
		TL_counter,
		Success_counter,
		IT
	);

//cout << "#############################################################" << endl;
	lik_efficient(
		network_key,
		links,
		nb_links,
		N,
		min_k,
		max_k,
		k_int,
		InitTemp,
		CR,
		TL,
		Max_Success,
		LimitIT,
		epsilon,
		TL_counter,
		Success_counter,
		IT
		);

//	cout << "LOUVAIN METHOD" << endl;
	louvain_convert_2b(links, nb_links, filename);
	char binary_file[256] = "louvain_binary.txt";
	vector<int> louvain_solution = louvain_clustering(binary_file, network_key, links, nb_links);
//	cout << "#############################################################" << endl;


	network_key++;
	//cout << "#############################################################" << endl;


						}
					}
				}
			}
		}
	}//End LFR parameter for loops
//cout << minnumber << "***" << maxnumber << endl; 
//CALCULATE TIME REQUIRED TO RUN ENTIRE PROGRAM
		//cout << fixed;

	cout << setprecision(10);
	clock_t end_new = clock();
	double time = (double)(end_new - start_new) / CLOCKS_PER_SEC;/*(double)(end_new - start_new) * 1000.0 / CLOCKS_PER_SEC  for milliseconds*/

	cout << "Program ended in "<< time <<" seconds. Hit enter to exit.";
	cin.get(); cin.get();

}//END main()