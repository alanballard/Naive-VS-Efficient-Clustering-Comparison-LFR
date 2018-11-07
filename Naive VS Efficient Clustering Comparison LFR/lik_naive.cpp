#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>	
#include <time.h>		//for calculating run-time
#include <string>		//for outputting file name
#include <iomanip>
#include <math.h>       //need for exponent fct.
#include <algorithm>    //need for mininimum fct.
#include <vector>
#include <cctype>		//for yes/no question
#include <random>

using namespace std;

int initial_pop(vector<int>& group_size, vector<int>& cluster_membership, int k, int N); //create initial membership vector
vector<int> mutate_membership(vector<int> current_membership, vector<int> group_size); //mutate membership vector

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
)
{
	int final_k = 0;
	long double final_BIC = 0;
	long double final_LOGLIK = 0;
	long double final_MODULARITY = 0;
	long double final_TESTSTAT = 0;
	vector<int> final_membership(N, 0);
	long double final_TIME = 0;

	//Repeat Simulated Annealing process for each k in the range of specified cluster numbers
	for (int k = min_k; k <= max_k; k = k + k_int)
	{
		//long double diff = 1.0e-10; //Used to check for absolute differences between 2 numbers

		int STG1 = 0; int STG2 = 0; int STG3 = 0;
		//cout << "NAIVE LIKELIHOOD METHOD, k=" << k << endl;

		clock_t start_new = clock(); //Time until solution

		//Step 1: Initialize kxk current/proposed OBS and POS matrices, and current/proposed membership matrices
				//Data types available: http://www.cplusplus.com/doc/tutorial/variables/
		vector< vector <unsigned long long int >> OBS(k, vector<unsigned long long int>(k, 0));
		vector< vector <unsigned long long int >> POS(k, vector<unsigned long long int>(k, 0));

		vector< vector <unsigned long long int >> PROP_OBS(k, vector<unsigned long long int>(k, 0));
		vector< vector <unsigned long long int >> PROP_POS(k, vector<unsigned long long int>(k, 0));

		//Initialize and populate cluster assignment and size vectors
		vector<int> group_size(k, 0);		//used in initial pop to store number of nodes in each group 1:k
		vector<int> cluster_membership;		//holds cluster assignments of N vertices

		vector<int> prop_group_size(k, 0);	//used in store number of nodes in each group 1:k
		vector<int> proposed_membership;	//holds cluster assignments of N vertices

		//vertex to be reassigned and the cluster to which it will be reassigned
		vector<int> reassignment;	//stores results of mutate_membership() call:
		int l_star;				//The vertex to be reassigned
		int clust_l;			//L-star's cluster membership before reassignment
		int clust_j;			//L-star's cluster membership after reassignment

		long double l_0_w = 0; //loglikelihood of within-clusters
		long double l_0_b = 0; //loglikelihood of between-clusters
		long double param_i_0 = 0; //Check validity of between-cluster parameter i in current membership
		//long double check_param = 0;
		long double check_param = 0; //Used to check validity of individual between-cluster observed/possible ratios in both current and proposed membership
		long double param_bw_0 = 0; //Between-cluster parameter in current membership
		long double l_1_w = 0; //loglikelihood of within-clusters
		long double l_1_b = 0; //loglikelihood of between-clusters
		long double param_i_1 = 0; //Check validity of between-cluster parameter i in proposed membership
		long double param_bw_1 = 0; //Between-cluster parameter in proposed membership

		unsigned long long int TOT_OBS = 0; //Total between-cluster observed edge count for current membership
		unsigned long long int TOT_POS = 0; //Total possible number of between-cluster edges for current membership
		unsigned long long int TOT_OBS_PROP = 0; //Total between-cluster observed edge count for proposed membership
		unsigned long long int TOT_POS_PROP = 0; //Total possible number of between-cluster edges for proposed membership

		//loglikelihoods
		long double l_0; //loglikelihood of current membership assignment
		long double l_1; //loglikelihood of proposed membership assignment
		long double delta_l; //difference in loglikelihoods

		//Values used to calculate modularity
		long double e_mm = 0; // the total fraction of edges between vertices within the same clusters
		vector<long double> a_m(k, 0); //the squared fraction of edges between vertices of differing clusters, by cluster
		long double a_mm = 0; //the total squared fraction of edges between vertices of differing clusters
		long double mod_0 = 0; //modularity of solution membership assignment

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
			TL_counter++; //Iterate number of attempts at current temperature

			//After TL attempts at reclustering at the current temperature, or Max_Success # of successes at curren temperature
			//Decrease temperature
			if (TL_counter >= TL || Success_counter >= Max_Success)
			{
				IT = CR * IT;
				TL_counter = 0;
				Success_counter = 0;
			}

			//Generate initial clustering (1st run only). This is the initial "current" clustering.
			//On runs >1, the "current" clustering will either be the one generated here, or a proposed clustering that was accepted below.
			while (first_run == 1)
			{
//Step: Make initial membership
				//Populate cluster assignment and size vectors
				initial_pop(group_size, cluster_membership, k, N);//assigns initial cluster assignments to cluster_membership, and counts cluster sizes
				first_run = 0;

//Step: Populate OBS and POS matrices based on this membership
				//Populate observed OBS edge weight vector. This requires 1 full traversal of adjacency list

				for (int i = 0; i < links.size(); i++)
				{
					for (int j = 0; j < links[i].size(); j++)
					{
						//for vertex i, cluster_membership[i] is the cluster membership, 1...k
						//cluster_membership[links[i][j].first] is the cluster membership of the vertex connected to i
						//links[i][j].second is the weight of the edge connecting to vertex i (i if unweighted)
						//cout << "<" << cluster_membership[i] << ", " << cluster_membership[links[i][j].first] << ">" <<endl;
						//within cluster observed edges
						if (cluster_membership[i] == cluster_membership[links[i][j].first]) {

							OBS[cluster_membership[i] - 1][cluster_membership[i] - 1]
								= OBS[cluster_membership[i] - 1][cluster_membership[i] - 1] 
									+ links[i][j].second;
						}
						//between cluster observed edges
						else
						{
							OBS[cluster_membership[i] - 1][cluster_membership[links[i][j].first] - 1]
								= OBS[cluster_membership[i] - 1][cluster_membership[links[i][j].first] - 1]
									+ links[i][j].second;
						}
					}
				}
				//Note: Since the adjacency list contains 2 entries for every edge and the number of observed edges are counted
				//over the adjacency list, the total observed count needs to be divided by 2.
				for (int i = 0; i < k ; i++) { OBS[i][i] = OBS[i][i] / 2; }

				//Populate possible within-cluster edge count elements. This is done using the cluster size information
				for (int i = 0; i < k; i++) { POS[i][i] = (unsigned long long int)( group_size[i] * (group_size[i] - 1) / 2); }

				//Populate possible between-cluster edge count elements.
				for (int i = 0; i < k - 1; i++) {
					for (int j = i + 1; j < k; j++)
					{	//add number of edges possible between cluster i and cluster j
						POS[i][j] = POS[i][j] + (unsigned long long int) (group_size[i] * group_size[j]);
						POS[j][i] = POS[j][i] + (unsigned long long int) (group_size[j] * group_size[i]);
					}
				}

/*cout << "OBS matrix" << endl;
				for (int i = 0; i < k; i++)
				{
					for (int j = 0; j < k; j++)
					{
						cout << OBS[i][j] << "...";
					} cout << endl;
				} cout << endl;
				cout << "*******************" << endl;
				cout << "*******************" << endl;
cout << "POS matrix" << endl;
				for (int i = 0; i < k; i++)
				{
					for (int j = 0; j < k; j++)
					{
						cout << POS[i][j] << "...";
					} cout << endl;
				} cout << endl;
				cout << "*******************" << endl;
/*int OBS_CHECK=0;
for (int i = 0; i < k; i++) {
	for (int j = i; j < k; j++) { 
		OBS_CHECK = OBS_CHECK + OBS[i][j]; 
		cout << OBS_CHECK << endl;
	}	
}
if (OBS_CHECK != nb_links) { cout << "OBS1=" << OBS_CHECK << endl; cin.get(); } //Should be 78
int POS_CHECK = 0;
for (int i = 0; i < k; i++) { for (int j = i; j < k; j++) { POS_CHECK = POS_CHECK + POS[i][j];cout << POS_CHECK << endl;
} }
if (POS_CHECK != N * (N - 1) / 2) { cout << "POS1=" << POS_CHECK << endl; cin.get(); } //Should be 34 choose 2 = 561
*/

//Step: Calculate full loglikelihood of this membership (Binomial)
				l_0_w = 0; //loglikelihood of within-clusters
				l_0_b = 0; //loglikelihood of between-clusters
				//Get loglikelihood of within-cluster edges						
				for (int i = 0; i < k; i++)
				{
					param_i_0 = 0;
					param_i_0 = OBS[i][i] / (long double)POS[i][i];// parameter for within-cluster i density
					//if (fabs(param_i_0 - 0) < diff || fabs(param_i_0 - 1) < diff) { l_0_w = l_0_w; }
					if (param_i_0 == 0 || param_i_0 == 1) { l_0_w = l_0_w; }
					else 
					{
						l_0_w = l_0_w + (OBS[i][i])*log(param_i_0) + (POS[i][i] - OBS[i][i])*log(1 - param_i_0);
					}
//cout<<"P["<<i<<"]="<< param_i_0 <<", ";
				}
//cin.get();


				//Get loglikelihood of between-cluster edges
				TOT_OBS = 0;
				TOT_POS = 0;
				for (int i = 0; i < k - 1; i++)
				{
					for (int j = i + 1; j < k; j++)
					{
/*						check_param = 0;
						check_param = OBS[i][j] / (long double)POS[i][j];
							//if (fabs(check_param - 0) < diff || fabs(check_param - 1) < diff)
							if (check_param == 0 || check_param == 1)
							{ 
								TOT_OBS = TOT_OBS; 
								TOT_POS = TOT_POS;
							}
							else*/
							{//Modeling between-cluster vertices with a single parameter. Loglikelihood calculate below, outside of i/j loop
								TOT_OBS = TOT_OBS + OBS[i][j];
								TOT_POS = TOT_POS + POS[i][j];
							}
					}
				}
				param_bw_0 = 0;
				param_bw_0 = TOT_OBS / (long double)TOT_POS;
				//if (fabs(param_bw_0 - 0) < diff || fabs(param_bw_0 - 1) < diff) { l_0_b = 0; }
				if (param_bw_0 == 0 || param_bw_0 == 1) { l_0_b = 0; }
				else {
					//Loglikelihood of between cluster edges when modeled with a single parameter
					l_0_b = (TOT_OBS)*log(param_bw_0) + ((long double)TOT_POS - TOT_OBS)*log(1 - param_bw_0);
				}
//cout<<"P[bw]="<< param_bw_0 <<", "<<endl;

				l_0 = l_0_w + l_0_b;
			}//end while 1st run loop

			//Enter SA Algorithm
//Step: Propose a new clustering by mutating current clustering
			//Select Random Vertex l* for reassignment and select cluster to which it will be moved
			//This returns a) The vertex to be reassigned and b) the cluster to which it will be moved.
			reassignment = mutate_membership(cluster_membership, group_size);	//Proposed membership vector due to reassignment
			l_star = reassignment[0];				//The vertex to be reassigned
			clust_l = cluster_membership[l_star];	//l-star's cluster membership before reassignment
			clust_j = reassignment[1];				//l-star's cluster membership after reassignment

			proposed_membership = cluster_membership;
			proposed_membership[l_star] = clust_j;

			prop_group_size = group_size;
			prop_group_size[clust_l - 1]--; //cluster l loses 1 vertex
			prop_group_size[clust_j - 1]++; //cluster j gains 1 vertex

//Step: Populate OBS and POS matrices based on this membership
				//Populate observed OBS edge weight vector. This requires 1 full traversal of adjacency list
			for (int i = 0; i < links.size(); i++)
			{
				for (int j = 0; j < links[i].size(); j++)
				{
					//for vertex i
					//cluster_membership[i] is the cluster membership, 1...k
					//cluster_membership[links[i][j].first] is the cluster membership of the vertex connected to i
					//links[i][j].second is the weight of the edge connecting to vertex i (i if unweighted)
					//cout << "<" << cluster_membership[i] << ", " << cluster_membership[links[i][j].first] << ">" <<endl;
					//within cluster observed edges
					if (proposed_membership[i] == proposed_membership[links[i][j].first]) {

						PROP_OBS[proposed_membership[i] - 1][proposed_membership[i] - 1]
							= PROP_OBS[proposed_membership[i] - 1][proposed_membership[i] - 1]
							+ links[i][j].second;
					}
					//between cluster observed edges
					else
					{
						PROP_OBS[proposed_membership[i] - 1][proposed_membership[links[i][j].first] - 1]
							= PROP_OBS[proposed_membership[i] - 1][proposed_membership[links[i][j].first] - 1]
							+ links[i][j].second;
					}
				}
			}
			//Note: Since the adjacency list contains 2 entries for every edge and the number of observed edges are counted
			//over the adjacency list, the total observed count needs to be divided by 2.
			for (int i = 0; i < k; i++) { PROP_OBS[i][i] = PROP_OBS[i][i] / 2; }

			//Populate possible within-cluster edge count elements. This is done using the cluster size information
			for (int i = 0; i < k; i++) { PROP_POS[i][i] = (unsigned long long int) (prop_group_size[i] * (prop_group_size[i] - 1)) / 2; }

			//Populate possible between-cluster edge count elements.
			for (int i = 0; i < k - 1; i++) {
				for (int j = i + 1; j < k; j++)
				{	//add number of edges possible between cluster i and cluster j
					PROP_POS[i][j] = PROP_POS[i][j] + (unsigned long long int) (prop_group_size[i] * prop_group_size[j]);
					PROP_POS[j][i] = PROP_POS[j][i] + (unsigned long long int) (prop_group_size[j] * prop_group_size[i]);
				}
			}
/*cout << "OBS matrix" << endl;
				for (int i = 0; i < k; i++)
				{
					for (int j = 0; j < k; j++)
					{
						cout << PROP_OBS[i][j] << "...";
					} cout << endl;
				} cout << endl;
				cout << "*******************" << endl;
				cout << "*******************" << endl;
cout << "POS matrix" << endl;
				for (int i = 0; i < k; i++)
				{
					for (int j = 0; j < k; j++)
					{
						cout << PROP_POS[i][j] << "...";
					} cout << endl;
				} cout << endl;
				cout << "*******************" << endl;
int OBS_CHECK=0;
for (int i = 0; i < k; i++) {
	for (int j = i; j < k; j++) { 
		OBS_CHECK = OBS_CHECK + PROP_OBS[i][j];
		cout << OBS_CHECK << endl;
	}	
}
if (OBS_CHECK != nb_links) { cout << "OBS2=" << OBS_CHECK << endl; cin.get(); } //Should be 78
int POS_CHECK = 0;
for (int i = 0; i < k; i++) { for (int j = i; j < k; j++) { POS_CHECK = POS_CHECK + PROP_POS[i][j];cout << POS_CHECK << endl;
} }
if (POS_CHECK != N * (N - 1) / 2) { cout << "POS2=" << POS_CHECK << endl; cin.get(); } //Should be 34 choose 2 = 561
cin.get();
*/

//Step: Calculate full loglikelihood of proposed membership (Binomial)
				l_1_w = 0; //loglikelihood of within-clusters
				l_1_b = 0; //loglikelihood of between-clusters
				//Get loglikelihood of within-cluster edges
						
				for (int i = 0; i < k; i++)
				{
					param_i_1 = 0;
					param_i_1 = PROP_OBS[i][i] / (long double)PROP_POS[i][i];// parameter for within-cluster i density
					
					//if (fabs(param_i_1 - 0) < diff || fabs(param_i_1 - 1) < diff) { l_1_w = l_1_w; }
					if (param_i_1 == 0 || param_i_1 == 1) { l_1_w = l_1_w; }
					else 
					{
						l_1_w = l_1_w + (PROP_OBS[i][i])*log(param_i_1) + (PROP_POS[i][i] - PROP_OBS[i][i])*log(1 - param_i_1);
					}
//cout<<"P["<<i<<"]="<< param_i_1 <<", ";

				}


				//Get loglikelihood of between-cluster edges
				TOT_OBS_PROP = 0;
				TOT_POS_PROP = 0;
				for (int i = 0; i < k - 1; i++)
				{
					for (int j = i + 1; j < k; j++)
					{
/*						check_param = 0;
						check_param = PROP_OBS[i][j] / (long double)PROP_POS[i][j];
							//if(fabs(check_param-0)<diff|| fabs(check_param - 1) < diff)
							if (check_param == 0 || check_param == 1)
							{ 
								TOT_OBS_PROP = TOT_OBS_PROP;
								TOT_POS_PROP = TOT_POS_PROP;
							}
							else*/
							{//Modeling between-cluster vertices with a single parameter. Loglikelihood calculate below, outside of i/j loop
								TOT_OBS_PROP = TOT_OBS_PROP + PROP_OBS[i][j];
								TOT_POS_PROP = TOT_POS_PROP + PROP_POS[i][j];
							}
					}
				}
				param_bw_1 = 0;
				param_bw_1 = TOT_OBS_PROP / (long double)TOT_POS_PROP;

				
				//if (fabs(param_bw_1 - 0) <diff || fabs(param_bw_1 - 1) < diff) { l_1_b = 0; }
				if (param_bw_1 == 0 || param_bw_1 == 1) { l_1_b = 0; }
				else {
					//Loglikelihood of between cluster edges when modeled with a single parameter
					l_1_b = (TOT_OBS_PROP)*log(param_bw_1) + ((long double)TOT_POS_PROP - TOT_OBS_PROP)*log(1 - param_bw_1);
				}

//cout<<"P[bw]="<< param_bw_1 <<", "<<endl;
//cin.get();
				l_1 = l_1_b + l_1_w;

//Step: Calculate change in log-likelihood (Binomial) due to reclustering
			//log(new/old)=log(new)-log(old)

			//If all densities in either the new or the old likelihood were discarded (and thus set to zero),
			//it is not possible to evaluate a change-in-likelihood
			//set entire change-in-likelihood to zero so a new clustering will be considered
delta_l = 0;
			if (l_0 == 0 || l_1 == 0){ delta_l = 0; }
			else{ delta_l = l_1 - l_0; }
			//if delta_l>0, then l_1>l_0, and the reclustered membership is more "likely"
			//Otherwise, l_1<l_0 and the current membership is more "likely"
			
			int STAGE = 0; //Will indicate when/if a proposed recustering was accepted or rejected.

			if (delta_l > epsilon)
			{//ACCEPT PROPOSED CLUSTERING AS CURRENT (best) CLUSTERING
//Step: Compare the two likelihoods. If accept new likelihood, then:
				//Proposed OBS/POS matrices becomes current OBS/POS matrices
				//Proposed membership vector becomes current membership vector
				//Proposed group (cluster) size vector becomes current group size vector
				//Proposed loglikelihood becomes current loglikelihood.		
				STAGE = 1;
				//current clustering <- proposed clustering
				cluster_membership = proposed_membership;
				//Update OBS vector
				OBS = PROP_OBS;
				//Update POS vector
				POS = PROP_POS;
				//Update group size vector
				group_size = prop_group_size;
				//Update loglikelihood
				l_0 = l_1;
				//Iterate number of successes at this temperature
				Success_counter++;
				//zero out PROP_OBS, PROP_POS, l_1 , delta_l	
				//These will be repopulated after the next clustering is proposed
				for (int i = 0; i < k; i++) {
					for (int j = 0; j < k; j++) {
						PROP_OBS[i][j] = 0;
						PROP_POS[i][j] = 0;
					}
				}
				l_1 = 0;
				//delta_l = 0;
STG1++;
			}
			else
			{
				long double one = 1;      //used in MIN call to match data types.

				double uni_draw = rand() / double(RAND_MAX);// dist_uni(e3);// dist_uni(engine)/*uniform_rand(rand())*/;
				double min_val = min(one, exp(delta_l / IT));

				if (uni_draw < min_val) //Accept proposed clustering
				{
//Step: If accept new loglikelihood, then:
										//Proposed OBS/POS matrices becomes current OBS/POS matrices
										//Proposed membership vector becomes current membership vector
										//Proposed group (cluster) size vector becomes current group size vector
										//Proposed loglikelihood becomes current loglikelihood. 	
					STAGE = 2;
					//current clustering <- proposed clustering
					cluster_membership = proposed_membership;
					//Update OBS vector
					OBS = PROP_OBS;
					//Update POS vector
					POS = PROP_POS;
					//Update group size vector
					group_size = prop_group_size;
					//Update loglikelihood
					l_0 = l_1;

					//Iterate number of success at this temperature
					Success_counter++;
					//zero out PROP_OBS, PROP_POS, l_1 , delta_l	
					//These will be repopulated after the next clustering is proposed					
					for (int i = 0; i < k; i++) {
						for (int j = 0; j < k; j++) {
							PROP_OBS[i][j] = 0;
							PROP_POS[i][j] = 0;
						}
					}
					l_1 = 0;
					//delta_l = 0;
STG2++;

				}
				else //Reject proposed clustering
				{//Nothing happens here. The current cluster is retained and the process repeats itself.					
					STAGE = 3;

					//zero out PROP_OBS, PROP_POS, l_1 , delta_l	
					//These will be repopulated after the next clustering is proposed				
					for (int i = 0; i < k; i++) {
						for (int j = 0; j < k; j++) {
							PROP_OBS[i][j] = 0;
							PROP_POS[i][j] = 0;
						}
					}
					l_1 = 0;
					//delta_l = 0;
STG3++;

				}
			}
			//At this point, the proposed clustering has:
			//1. Been Accepted. membership vector/OBS/POS/cluster counts have beeen updated to reflect new clustering
			//2. Rejected. memberhsip vector/OBS/POS/cluster counts still correspond to the current clustering.

		} //end SA/IT>LimitIT loop
			//At this point, the proposed clustering has:
			//1. Been Accepted. membership vector/OBS/POS/cluster counts have beeen updated to reflect new clustering
			//2. Rejected. memberhsip vector/OBS/POS/cluster counts still correspond to the current clustering.	

/*		cout << "Number of success at each stage" << endl;
		cout << "<" << STG1 << ", " << STG2 << ", " << STG3 << ">" << endl;
		cout << "*******************" << endl;
		cout << "Solution cluster sizes" << endl;
		for (int i = 0; i < group_size.size(); i++) { cout << group_size[i] << "..."; }cout << endl;
		cout << "*******************" << endl;
		cout << "Solution OBS matrix" << endl;
		for (int i = 0; i < k; i++)
		{
			for (int j = 0; j < k; j++)
			{
				cout << OBS[i][j] << "...";
			} cout << endl;
		} cout << endl;
		cout << "*******************" << endl;
		cout << "*******************" << endl;
		cout << "Solution POS matrix" << endl;
		for (int i = 0; i < k; i++)
		{
			for (int j = 0; j < k; j++)
			{
				cout << POS[i][j] << "...";
			} cout << endl;
		} cout << endl;
		cout << "*******************" << endl;
		cout << "Solution cluster membership vector" << endl;
		for (int i = 0; i < cluster_membership.size(); i++) { cout << cluster_membership[i] << ", "; }cout << endl;
		cout << "Solution likelihood" << endl;
		cout << l_0 << endl;*/


		//CALCULATE TIME REQUIRED TO REACH SOLUTION
				//cout << fixed;
		cout << setprecision(10);
		clock_t end_new = clock();
		double time = (double)(end_new - start_new) / CLOCKS_PER_SEC; /*(double)(end_new - start_new) * 1000.0 / CLOCKS_PER_SEC  for milliseconds*/

		//At this point, a solution clustering has been identified.
		long double full_loglik_curr = 0;	//loglikelihood for the full clustering
		long double LO = 0;					//loglikelihood under null hypothesis: # of clusters=1
		long double BIC = 0;				//Bayesian Information Criterion
		long double TESTSTAT = 0;			//test statistics for cluster signifigance test
		//Get loglikelihood of solution clustering and loglikelihood under null hypothesis: # of clusters=1
		vector<long double> returned;


		full_loglik_curr = l_0;

//STEP: Calculate likelihood under null hypothesis of k=1 cluster
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
		LO_POS = N * (N - 1) / 2; //N choose 2
		LO = (LO_OBS)*log(LO_OBS / (long double)LO_POS) + (LO_POS - LO_OBS)*log(1 - (LO_OBS / (long double)LO_POS));

//Step: Calculate modularity of solution membership 
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
		mod_0 = e_mm - a_mm; //modularity of solution membership assignment


//Checking the observed and possible edge counts in the final solution
/*for(int i=0; i<k; i++){cout<<OBS[i][i]<<",,,"; }
long long int temp_OBS_BW = 0;
for (int i = 0; i < k; i++)
{
	for (int j = i+1; j < k; j++)
	{
		temp_OBS_BW = temp_OBS_BW + OBS[i][j];
	}
}
cout<< temp_OBS_BW <<endl;
for (int i = 0; i < k; i++) { cout << POS[i][i] << ",,,"; }
long long int temp_POS_BW = 0;
for (int i = 0; i < k; i++)
{
	for (int j = i + 1; j < k; j++)
	{
		temp_POS_BW = temp_POS_BW + POS[i][j];
	}
}
cout << temp_POS_BW << endl;*/

//Checking the null hypothesis stats
//cout << "LO=" << LO << ", LO_OBS=" << LO_OBS << ", LO_POS=" << LO_POS << endl;


//STEP: Calculate Bayesian Information Criterion for solution clustering
		//In old code, BIC = -2 * (LO + 0.5*(-bestfit)) + (k + 1)*log(N*(N - 1) / 2);
		//but bestfit = -(-2 * (LO - full_loglik_curr)), so this can be simplified:		
		//number of parameters=(k+1), number of observations=Nchoose2=N(N-1)/2, maximized loglikelihood=full_loglik_curr
		BIC = -2 * (full_loglik_curr)+(k + 1)*log(N*(N - 1) / (long double)2);
		TESTSTAT = -2 * (LO - full_loglik_curr);

//cout << "at k=" << k << ", modularity= " << mod_0 << ", loglik=" << full_loglik_curr << ", BIC=" << BIC << ", TESTSTAT=" << TESTSTAT<< ", time=" << time << " sec.s"<< endl;
		//Identify the optimal number of clusters, based on minimum Bayesian Information Criterion
		if (k == min_k) {
			final_k = k;
			final_BIC = BIC;
			final_LOGLIK = full_loglik_curr;
			final_TESTSTAT = TESTSTAT;
			final_membership = cluster_membership;
			final_MODULARITY = mod_0;
			final_TIME = time;
		}
		//if BIC < min BIC thus far, update variables
		if (BIC < final_BIC) {
			final_k = k;
			final_BIC = BIC;
			final_LOGLIK = full_loglik_curr;
			final_TESTSTAT = TESTSTAT;
			final_membership = cluster_membership;
			final_MODULARITY = mod_0;
			final_TIME = time;
		}

//STEP: Output solution clustering statistics at each value of k
		ostringstream atkstring;
		ofstream atkinfo;
		// If old file exists from a previous run, delete it.
		//remove("completion_lik_naive_all_sols.txt");
		atkstring << "completion_lik_naive_all_sols" << ".txt";
		atkinfo.open("completion_lik_naive_all_sols.txt", ios::app);
		atkinfo << network_key << "," << N << "," << nb_links << "," << k << "," << mod_0 << "," << full_loglik_curr << ", " << BIC << "," << TESTSTAT << "," << time << endl;
		atkinfo.close();

		//cout << "**********************************" << endl;
	} //end for-k loop

//STEP: Output best solution across all k
	ostringstream atkstring2;
	ofstream atkinfo2;
	// If old file exists from a previous run, delete it.
	//remove("completion_lik_naive_final_sol.txt");
	atkstring2 << "completion_lik_naive_final_sol" << ".txt";
	atkinfo2.open("completion_lik_naive_final_sol.txt", ios::app);
	atkinfo2 << network_key << "," << N << "," << nb_links << "," << final_k << "," << final_MODULARITY << "," << final_LOGLIK << "," << final_BIC << "," << final_TESTSTAT << "," << final_TIME <<endl;
	atkinfo2.close();

	//STEP: Output best solution clustering
	ostringstream bestzsol;
	ofstream bestzstorage;
	// If old file exists from a previous run, delete it.
	//remove("solution_lik_naive.txt");
	bestzsol << "solution_lik_naive" << ".txt";
	bestzstorage.open("solution_lik_naive.txt", ios::app);
	for (int i = 0; i < N; i++) { bestzstorage << final_membership[i];
	if (i < N - 1) { bestzstorage << ","; }
	}
	bestzstorage << endl;
	bestzstorage.close();


	return 0;
}