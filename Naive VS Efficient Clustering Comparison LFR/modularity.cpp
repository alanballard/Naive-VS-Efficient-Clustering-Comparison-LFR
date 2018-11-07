#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>	
#include <time.h>		//for calculating run-time
#include <chrono>
#include <string>		//for outputting file name
#include <iomanip>
#include <math.h>       //need for exponent fct.
#include <algorithm>    //need for mininimum fct.
#include <vector>
#include <cctype>		//for yes/no question

#include <random>

using namespace std;

using namespace chrono;

int initial_pop(vector<int>& group_size, vector<int>& cluster_membership, int k, int N/*, mt19937& engine*/); //create initial membership vector
vector<int> mutate_membership(vector<int> current_membership, vector<int> group_size/*, mt19937& engine*/); //mutate membership vector

int modularity(
	//int seed,
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
)
{
	//Values for the "best" solution
	int final_k = 0;
	long double final_BIC = 0;
	long double final_LOGLIK = 0;
	long double final_TESTSTAT = 0;
	long double final_MODULARITY = 0;
	vector<int> final_membership(N, 0);
	long double final_TIME = 0;

	//Repeat Simulated Annealing process for each k in the range of specified cluster numbers
	for (int k = min_k; k <= max_k; k = k + k_int)
	{
//		cout << "MODULARITY METHOD, k=" << k << endl;
//		mt19937 engine(seed);
//		random_device rd;
//		mt19937 e4(rd());
		uniform_real_distribution<double> dist_uni(0.0, 1.0);

		//measure time-to-solution
		clock_t start_new = clock();
		
//Step: Initialize kxk current/proposed OBS and POS matrices, and current/proposed membership vector

		//Note that OBS/POP are only used, at the end, to calculate the likelihood/BIC of the "best" solution identified by maximizing modularity
		//Data types available: http://www.cplusplus.com/doc/tutorial/variables/
		vector< vector <unsigned long long int >> OBS(k, vector<unsigned long long int>(k, 0));
		vector< vector <unsigned long long int >> POS(k, vector<unsigned long long int>(k, 0));
		
		//Initialize and populate cluster assignment and size vectors
		vector<int> group_size(k, 0);		//used in initial_pop to store number of nodes in each group 1:k
		vector<int> cluster_membership;		//holds cluster assignments of N vertices

		vector<int> prop_group_size(k, 0);	//used in store number of nodes in each group 1:k
		vector<int> proposed_membership;	//holds cluster assignments of N vertices

		//vertex to be reassigned and the cluster to which it will be reassigned
		vector<int> reassignment;	//stores results of mutate_membership() call:
		int l_star = 0;				//The vertex to be reassigned
		int clust_l = 0;			//L-star's cluster membership before reassignment
		int clust_j = 0;			//L-star's cluster membership after reassignment

		//Values used to calculate modularity
		long double e_mm = 0; // the total fraction of edges between vertices within the same clusters
		vector<long double> a_m(k, 0); //the squared fraction of edges between vertices of differing clusters, by cluster
		long double a_mm = 0; //the total squared fraction of edges between vertices of differing clusters
		long double mod_0 = 0; //modularity of current membership assignment
		long double mod_1 = 0; //modularity of proposed membership assignment
		long double delta_mod = 0; //difference in modularities

		//Flag to allow creation of initial clustering that will be mutated on later runs.
		int first_run = 1;

		//Reset initial temperature at the start of each k-loop
		IT = InitTemp;
		//Reset count length at the start of each k-loop
		TL_counter = 0;
		//Reset success counter at the start of each k-loop
		Success_counter = 0;

		//START Simulated Annealing Algorithm
		while (IT > LimitIT)
		{
			//Slowly decrease current temperature
			//IT = CR*IT;

			TL_counter++; //Iterate number of attempts at current temperature

			//After TL attempts at reclustering at the current temperature, or Max_Success # of successes at curren temperature, decrease temperature and reset counters
			if (TL_counter >= TL || Success_counter >= Max_Success)
			{	IT = CR*IT; TL_counter = 0; Success_counter = 0;	}

			//Generate initial clustering (1st run only). This is the initial current clustering. On runs >1, the "current clustering" will either be the one generated here, or a proposed clustering that was accepted below.
			while (first_run == 1)
			{
//Step: Make initial membership
				//Populate cluster assignment and size vectors
				initial_pop(group_size, cluster_membership, k, N/*, engine*/);//assigns initial cluster assignments to z, and counts cluster sizes
				first_run = 0; 

//Step: Calculate modularity of current membership 
			//"Finding and evaluating community structure in networks"
			//https://arxiv.org/pdf/cond-mat/0308217.pdf
			//Also: https://github.com/ivanbrugere/matlab-networks-toolbox/blob/master/modularityMetric.m			
			// e_mm=fraction of edges that connect community m to community m
			// a_i=the fraction of edges that connect to vertices in community i to vertices in community i
			// modularity Q =SUM[ (e_mm - a_m^2) ] over all m
			
			//1st step: calculate the fraction of edges between vertices within the same clusters
				//This requires 1 full traversal of adjacency list
				e_mm = 0;	a_mm = 0;	fill(a_m.begin(), a_m.end(), 0);
				for (int i = 0; i < links.size(); i++)
				{
					for (int j = 0; j < links[i].size(); j++)
					{
						//for vertex i
						//cluster_membership[i] is the cluster membership, 1...k
						//cluster_membership[links[i][j].first] is the cluster membership of the vertex connected to i
						//within cluster edges
						if (cluster_membership[i] == cluster_membership[links[i][j].first])
						{
							e_mm = e_mm + 1; //counting how many edges fall between vertices within the same cluster
							a_m[cluster_membership[i] - 1] = a_m[cluster_membership[i] - 1] + 1; //counting how many edges fall between vertices within cluster i
						}
						//between cluster edges
						else 
						{
							//a_m[i-1] contains the number of edges between vertices of cluster i and other clusters (but not between vertices of cluster i)
							a_m[cluster_membership[i] - 1] = a_m[cluster_membership[i] - 1] + 1; //counting how many edges fall between vertices of cluster i and other clusters
						}
					}
				}
				//Note: Since the adjacency list contains 2 entries for every edge and the number of observed edges are counted
				//over the adjacency list, the total edge count needs to be divided by 2.
				//finally, divide by the number of distinct edges to get the fraction of edges that fall between vertices within the same clusters				
				e_mm = e_mm / (2);
				e_mm = e_mm / (nb_links);

			//2nd step: calculate the total squarted fraction of edges that connect vertices of one cluster to another cluster
				for (int i = 0; i < k; i++) 
				{ 
					//a_m[i] is the squared fraction of edges that fall between vertices of cluster i and other clusters
					a_m[i] = a_m[i] / 2; //to account for duplication of edges in adjacency list
					a_m[i] = a_m[i] / nb_links; //divide by total number of distinct edges
					a_m[i] = a_m[i]* a_m[i]; //square, per Modularity formula
				} 
				for (int i = 0; i < k; i++) 
				{ 
					//total squared fraction of between-cluster edges
					a_mm = a_mm + a_m[i]; 
				} 
		mod_0 = e_mm - a_mm; //modularity of current membership assignment
			}//end 1st-run WHILE loop

			 //At this point we have: a cluster membership vector and a group size vector, modularity of current clustering, either generated on the first run, or of the current clustering for subsequent runs.
			//Enter SA Algorithm
//Step: Propose a new clustering by mutating current clustering
			 //Propose a new clustering. 
			 //Select Random Vertex l* for reassignment and select cluster to which it will be moved
			 //This returns a) The vertex to be reassigned and b) the cluster to which it will be moved.
			reassignment = mutate_membership(cluster_membership, group_size/*, engine*/);	//Proposed membership vector due to reassignment
			l_star = reassignment[0];				//The vertex to be reassigned
			clust_l = cluster_membership[l_star];	//l-star's cluster membership before reassignment
			clust_j = reassignment[1];				//l-star's cluster membership after reassignment
			proposed_membership = cluster_membership;
			proposed_membership[l_star] = clust_j;
			prop_group_size = group_size;
			prop_group_size[clust_l - 1]--; //cluster l loses 1 vertex
			prop_group_size[clust_j - 1]++; //cluster j gains 1 vertex

//Step: Calculate modularity of the proposed membership
			//1st step: calculate the fraction of edges between vertices within the same clusters
			e_mm = 0;	a_mm = 0;	fill(a_m.begin(), a_m.end(), 0);
	
			for (int i = 0; i < links.size(); i++)
			{
				for (int j = 0; j < links[i].size(); j++)
				{
					if (proposed_membership[i] == proposed_membership[links[i][j].first])
					{
						e_mm = e_mm + 1;
						a_m[proposed_membership[i] - 1] = a_m[proposed_membership[i] - 1] + 1;
					}
					else 
					{
						a_m[proposed_membership[i] - 1] = a_m[proposed_membership[i] - 1] + 1;
					}
				}
			}
			e_mm = e_mm / (2);
			e_mm = e_mm / (nb_links);

			//2nd step: calculate the total squared fraction of edges that connect vertices of one cluster to another cluster
			for (int i = 0; i < k; i++) 
			{ 
				a_m[i] = a_m[i] / 2;
				a_m[i] = a_m[i] / nb_links;
				a_m[i] = a_m[i] * a_m[i];
			} 
			for (int i = 0; i < k; i++) 
			{ 
				a_mm = a_mm + a_m[i];
			} 
		mod_1 = e_mm - a_mm; //modularity of current membership assignment

//Step: Calculate change in modularity due to reclustering
			delta_mod = mod_1 - mod_0;

			int STAGE = 0; //Will indicate when/if a proposed recustering was accepted or rejected.

			if (delta_mod > epsilon)
			{//ACCEPT PROPOSED CLUSTERING AS CURRENT (best) CLUSTERING
//Step: Compare the two modularities. If accept new mod, then:
			 //Proposed membership vector becomes current membership vector
			 //Proposed group (cluster) size vector becomes current group size vector
			 //Proposed modularity becomes current modularity	
				STAGE = 1;
				cluster_membership = proposed_membership;
				//Update group size vector
				group_size = prop_group_size;
				//Update modularity
				mod_0 = mod_1;
				//Iterate number of success at this temperature
				Success_counter++;

				mod_1 = 0;
				delta_mod = 0;
			}
			else
			{
				long double one = 1;      //used in MIN call to match data types.

				double uni_draw = rand() / double(RAND_MAX);//dist_uni(e4);// dist_uni(engine)/*uniform_rand(rand())*/;
				long double min_val = min(one, exp(delta_mod / IT)); 

				if (uni_draw < min_val) //Accept proposed clustering
				{
//Step:  Compare the two modularities. If accept new mod, then:
					//Proposed membership vector becomes current membership vector
					//Proposed group (cluster) size vector becomes current group size vector
					//Proposed mod becomes current mod. 		
					STAGE = 2;
					cluster_membership = proposed_membership;
					//Update group size vector
					group_size = prop_group_size;
					//Update modularity
					mod_0 = mod_1;

					//Iterate number of success at this temperature
					Success_counter++;

					mod_1 = 0;
					delta_mod = 0;
				}
				else //Reject proposed clustering
				{//Current cluster values are not updated. The current cluster is retained and the process repeats itself.	
					STAGE = 3;

					mod_1 = 0;
					delta_mod = 0;
				}
			}
			//At this point, the proposed clustering has:
			//1. Been Accepted. membership vector/cluster counts have beeen updated to reflect new clustering
			//2. Rejected. memberhsip vector/cluster counts still correspond to the current clustering.

		}//End SA LOOP

//CALCULATE TIME REQUIRE TO REACH SOLUTION
		//cout << fixed;
		cout << setprecision(10);		
		clock_t end_new = clock();
		double time = (double)(end_new - start_new) / CLOCKS_PER_SEC;/*(double)(end_new - start_new) * 1000.0 / CLOCKS_PER_SEC  for milliseconds*/


//In the case of modularity, we need to populate the OBS/POS matrices for the final membership here in order to calculate the likelihood for the modularity solution, so we can calculate BIC and a test statistic.
		for (int i = 0; i < links.size(); i++)
		{
			for (int j = 0; j < links[i].size(); j++)
			{
				//cout << "i=" << i << ", j=" << links[i][j].first << ", ct. of i's nbrs=" << links[i].size() << ", i's clust=" << cluster_membership[i] << ", j's cluster=" << cluster_membership[links[i][j].first] << ", weight=" << links[i][j].second;
				if (cluster_membership[i] == cluster_membership[links[i][j].first])
				{
					OBS[cluster_membership[i] - 1]
						[cluster_membership[i] - 1] =
						OBS[cluster_membership[i] - 1]
						[cluster_membership[i] - 1] + links[i][j].second;

				}
				else
				{
					OBS[cluster_membership[i] - 1]
						[cluster_membership[links[i][j].first] - 1] =
						OBS[cluster_membership[i] - 1]
						[cluster_membership[links[i][j].first] - 1] + links[i][j].second;
					//And to the corresponding lower triangle
					//OBS[cluster_membership[j] - 1][cluster_membership[i] - 1] = OBS[cluster_membership[j] - 1][cluster_membership[i] - 1] + links[i][j].second;
				}
			}//j'th neighbor to i'th vertex
		}//i'th vertex in adjacency list

		 //ADDRESS DOUBLE-COUNTS
		 //Note: Since the adjacency list contains 2 entries for every edge and the number of observed edges are counted
		 //over the adjacency list, the total observed count needs to be divided by 2.
		for (int i = 0; i < k; i++) { OBS[i][i] = OBS[i][i] / 2; }

		//Populate possible edge count vector. This is done using the cluster size information
		for (int i = 0; i < k; i++)
		{
			POS[i][i] = ((unsigned long long int) group_size[i] * (group_size[i] - 1)) / 2;
		}
		//Since we can use j=i+1 here, we don't have double counts.
		for (int i = 0; i < k - 1; i++) {
			for (int j = i + 1; j < k; j++)
			{
				//add number of edges possible between cluster i and cluster j
				POS[i][j] = POS[i][j] + (unsigned long long int) group_size[i] * group_size[j];
				POS[j][i] = POS[j][i] + (unsigned long long int) group_size[j] * group_size[i];
			}
		}


		 //At this point, a solution clustering has been identified.
		long double full_loglik_curr = 0;	//loglikelihood for the solution clustering
		long double LO = 0;					//loglikelihood under null hypothesis: # of clusters=1
		long double BIC = 0;				//Bayesian Information Criterion
		long double TESTSTAT = 0;

		//Get loglikelihood of solution clustering and loglikelihood under null hypothesis: # of clusters=1
		double l_0_w = 0; //loglikelihood of within-clusters
		double l_0_b = 0; //loglikelihood of between-clusters
		//l_0	=log(product of within and between-cluster likelihoods)
		//		=sum(within and between-cluster loglikelihoods)
		//		=(loglikelihood of within-clusters)+(loglikelihood of between-clusters);
		//		=l_0_w+l_0_b

		//Get loglikelihood of within-cluster edges
		for (int i = 0; i < k; i++) {
			if (OBS[i][i] != 0 && (OBS[i][i] != POS[i][i])) {
				l_0_w = l_0_w +
					(OBS[i][i])*log(OBS[i][i] / (long double)POS[i][i]) + (POS[i][i] - OBS[i][i])*log(1 - OBS[i][i] / (long double)POS[i][i]);
			}
		}

		//Get loglikelihood of between-cluster edges
		unsigned long long int TOT_OBS = 0;
		unsigned long long int TOT_POS = 0;
		for (int i = 0; i < k; i++) {
			for (int j = i + 1; j < k; j++) {

				if (OBS[i][j] != 0 ) {
					//Modeling between-cluster vertices with a single parameter. Likelihood calculate below, outside of i/j loop
					TOT_OBS = OBS[i][j];
					TOT_POS = POS[i][j];
					//To model between-cluster vertices individually with kchoose2 paramters.
					//	l_0_b = l_0_b +
					//		(OBS[i][j])*log(OBS[i][j] / (long double)POS[i][j]) + (POS[i][j] - OBS[i][j])*log(1 - OBS[i][j] / (long double)POS[i][j]);
				}
			}
		}
		//Loglikelihood of between cluster edges when modeled with a single paramter
		l_0_b = (TOT_OBS)*log(TOT_OBS / (long double)TOT_POS) + ((long double)TOT_POS - TOT_OBS)*log(1 - (TOT_OBS / (long double)TOT_POS));

		//loglikelihood of solution clustering
		full_loglik_curr = l_0_w + l_0_b;

		//loglikelihood under null assumption of only 1 cluster in network
		long long int LO_OBS = 0;
		long long int LO_POS = 0;
		for (int i = 0; i < k; i++)
		{
			for (int j = i; j < k; j++)
			{
				LO_OBS = LO_OBS + OBS[i][j];
				//cout << "OBS[i][j]="<<OBS[i][j] <<", LO_OBS="<< LO_OBS << endl;
			}
		}
		LO_POS = N*(N - 1) / 2; //N choose 2
		
		LO = (LO_OBS)*log(LO_OBS / (long double)LO_POS) + (LO_POS - LO_OBS)*log(1 - (LO_OBS / (long double)LO_POS));
		
		//Indicates how many, if any, of the clusters have a observed edge weight/count=0.
		int flag = 0;
		for (int i = 0; i < k; i++) {
			for (int j = i; j < k; j++) {
				if (OBS[i][j] == 0) { flag++; }
			}
		}

		//Calculate Bayesian Information Criterion for solution clustering
		//In old code, BIC = -2 * (LO + 0.5*(-bestfit)) + (k + 1)*log(N*(N - 1) / 2);
		//but bestfit = -(-2 * (LO - full_loglik_curr)), so this can be simplified:		
		//number of parameters=(k+1), number of observations=Nchoose2=N(N-1)/2, maximized loglikelihood=full_loglik_curr
//NOTE: In the case of modularity, full_loglik_curr is actually just the likelihood of the modularity solution, which may be different than the solution that maximizes likelihood
//Therefore the BIC, loglikelihood and test statistic for the solution clustering under modularity are not particularly meaningful. Just calculating here to be thorough.
		BIC = -2 * (full_loglik_curr)+(k + 1)*log(N*(N - 1) / (long double)2);
		TESTSTAT = -2 * (LO - full_loglik_curr);

//cout << "at k=" << k << ", modularity= "<<mod_0 <<", loglik=" << full_loglik_curr << ", BIC=" << BIC << ", TESTAT="<<TESTSTAT<< ", time=" <<time<< " sec.s" << endl;

		//Identify the optimum number of clusters, based on minimum Bayesian Information Criterion
		//populate variables on 1st run (min k# of clusters evaluated)
		if (k == min_k) {
			final_k = k;
			final_BIC = BIC;
			final_LOGLIK = full_loglik_curr;
			final_TESTSTAT = TESTSTAT;
			final_MODULARITY = mod_0;
			final_membership = cluster_membership;
			final_TIME = time;
		}

		if (mod_0 > final_MODULARITY) {
			final_k = k;
			final_BIC = BIC;
			final_LOGLIK = full_loglik_curr;
			final_TESTSTAT = TESTSTAT;
			final_MODULARITY = mod_0;
			final_membership = cluster_membership;
			final_TIME = time;
		}

		ostringstream atkstring;
		ofstream atkinfo;
		// If old file exists from a previous run, delete it.
		//remove("completion_mod_all_sols.txt");
		atkstring << "completion_mod_all_sols" << ".txt";
		atkinfo.open("completion_mod_all_sols.txt", ios::app);
		atkinfo << network_key << "," << N << "," << nb_links << "," << k << "," << mod_0 << "," << full_loglik_curr <<","<< BIC << "," << TESTSTAT << "," << time << endl;
		atkinfo.close();


		//cout << "**********************************" << endl;
	}//end for-k loop
	 
	ostringstream atkstring2;
	ofstream atkinfo2;
	// If old file exists from a previous run, delete it.
	//remove("completion_mod_final_sol.txt");
	atkstring2 << "completion_mod_final_sol" << ".txt";
	atkinfo2.open("completion_mod_final_sol.txt", ios::app);
	atkinfo2 << network_key << "," << N << "," << nb_links << "," << final_k << "," << final_MODULARITY << "," << final_LOGLIK << "," << final_BIC << "," << final_TESTSTAT << "," << final_TIME << endl;
	atkinfo2.close();


	 //PRINT OUT MODULARITY SOLUTION TO FILE
	 ostringstream bestzsol;
	 ofstream bestzstorage;
	 // If old file exists from a previous run, delete it.
	 //remove("solution_mod.txt");
	 bestzsol << "solution_mod" << ".txt";
	 bestzstorage.open("solution_mod.txt", ios::app);
	 for (int i = 0; i < N; i++) {
		 bestzstorage << final_membership[i];
		 if (i < N - 1) { bestzstorage << ","; }
	 }
	 bestzstorage << endl;
	 bestzstorage.close();
	 
	return 0;
}
