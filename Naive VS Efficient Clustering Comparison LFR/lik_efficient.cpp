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
	)
{
	int final_k = 0;
	long double final_BIC = 0;
	long double final_LOGLIK = 0;
	long double final_MODULARITY = 0;
	long double final_TESTSTAT = 0;
	vector<int> final_membership(N,0);
	long double final_TIME = 0;

	//Repeat Simulated Annealing process for each k in the range of specified cluster numbers
	for (int k = min_k; k <= max_k; k=k + k_int)
	{
		//cout << "EFFICIENT LIKELIHOOD METHOD, k=" << k << endl;
	
		//random_device rd;
		//mt19937 e3(rd());
		//uniform_real_distribution<double> dist_uni(0.0, 1.0);
		
		clock_t start_new = clock(); //Time until solution
		
//Step: Initialize k+1 current/proposed OBS and POS vectors, and current/proposed membership vector
		vector<unsigned long long int> OBS(k + 1, 0);//elements 0 thru k-1 are within cluster values, element k is between-cluster value
		vector<unsigned long long int> POS(k + 1, 0);//elements 0 thru k-1 are within cluster values, element k is between-cluster value
		//Data types available: http://www.cplusplus.com/doc/tutorial/variables/

		//Initialize and populate cluster assignment and size vectors
		vector<int> group_size(k, 0);//used in initial pop to store number of nodes in each group 1:k
		vector<int> cluster_membership;	//holds cluster assignments of N vertices

		//vertex to be reassigned and the cluster to which it will be reassigned
		vector<int> reassignment;
		int l_star;		//The vertex to be reassigned
		int clust_l;	//L-star's cluster membership before reassignment
		int clust_j;	//L-star's cluster membership after reassignment

		//Total weight of links connected to L*, by cluster (These are the Z_ values in the paper)
		vector<int> Z(k);

		//Total weights and group counts used to develop parameters for proposed clustering
		int Z_l;
		int Z_j;
		int n_l;
		int n_j;

		//Observed edge weights and possible edge counts used to develop parameters for proposed clustering.
		unsigned long long int OBS_l_1;
		unsigned long long int OBS_j_1;
		unsigned long long int OBS_bw_1;
		unsigned long long int POS_l_1;
		unsigned long long int POS_j_1;
		unsigned long long int POS_bw_1;

		//Parameters for current clustering	
		long double param_l_0;
		long double param_j_0;
		long double param_bw_0;

		//Parameters for proposed clustering
		long double param_l_1;
		long double param_j_1;
		long double param_bw_1;

		//Partials loglikelihoods for current and proposed clustering, and the change in loglikelihood.
		long double l_0;
		long double l_1;
		long double delta_l;

		//Densities in the denominator in the change-in-loglik formula [for current clustering]
		long double l_density_0;
		long double j_density_0;
		long double bw_density_0;

		//Densities in the numerator in the change-in-loglik formula [for proposed clustering]
		long double l_density_1;
		long double j_density_1;
		long double bw_density_1;

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
				IT = CR*IT;
				TL_counter = 0;
				Success_counter = 0;
			}

			//Generate initial clustering (1st run only). This is the initial current clustering.
			//On runs >1, the "current clustering" will either be the one generated here, or a proposed clustering that was accepted below.
			while (first_run == 1)
			{
//Step: Make initial membership
				//Populate cluster assignment and size vectors
				initial_pop(group_size, cluster_membership, k, N/*, engine*/);//assigns initial cluster assignments to z, and counts cluster sizes
				first_run = 0;

//Step: Populate OBS and POS vectors based on this membership
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
						if (cluster_membership[i] == cluster_membership[links[i][j].first]){
							OBS[cluster_membership[i] - 1] = OBS[cluster_membership[i] - 1] + links[i][j].second;
						}
						//between cluster observed edges
						else
						{
							OBS[k] = OBS[k] + links[i][j].second;
						}
					}
				}
				//Note: Since the adjacency list contains 2 entries for every edge and the number of observed edges are counted
				//over the adjacency list, the total observed count needs to be divided by 2.
				for (int i = 0; i < k + 1; i++){ OBS[i] = OBS[i] / 2; }

				//Populate possible edge count vector. This is done using the cluster size information
				for (int i = 0; i < k; i++)
				{
					POS[i] = ((unsigned long long int) group_size[i] * (group_size[i] - 1)) / 2;
				}
				for (int i = 0; i < k - 1; i++){
					for (int j = i + 1; j < k; j++){
						//add number of edges possible between cluster i and cluster j
						POS[k] = POS[k] + (unsigned long long int) group_size[i] * group_size[j];
					}
				}
			}//end 1st-run WHILE loop
			
			//Enter SA Algorithm
//Step: Propose a new clustering by mutating current clustering
			//Select Random Vertex l* for reassignment and select cluster to which it will be moved
			//This returns a) The vertex to be reassigned and b) the cluster to which it will be moved.
			reassignment = mutate_membership(cluster_membership, group_size/*, engine*/);	//Proposed membership vector due to reassignment
			l_star = reassignment[0];				//The vertex to be reassigned
			clust_l = cluster_membership[l_star];	//l-star's cluster membership before reassignment
			clust_j = reassignment[1];				//l-star's cluster membership after reassignment

//Step: Calculate abbreviated loglikelihood of the current membership (Binomial)
			//Get denominator of change-in-loglik equation [l,j,bw for current clustering]
			//Calculate parameters for current clustering.
			param_l_0 = OBS[clust_l - 1] / (long double)POS[clust_l - 1];
			param_j_0 = OBS[clust_j - 1] / (long double)POS[clust_j - 1];
			param_bw_0 = OBS[k] / (long double)POS[k];

			l_density_0 = 0;
			j_density_0 = 0;
			bw_density_0 = 0;
			//Calculate l,j,bw densities for current clustering, the denominator in the change-in-loglik formula
			//if unweighted, all edges have weight=1 and Binomial density is appropriate

				//Note: we are calculating the log of the product of 3 binomial densities each for
				//the numerator and denominator of the change-in-loglik formula.
				//If a parameter estimate is 0 or 1, this results in a density estimate of 0 or INF, which throws an error when the log is taken.
				//Thus, we rewrite the log-product as a sum of logs. If a density has a parameter estimate of 0 or 1, 
				//we force that density's log-value to be zero so that it does not throw an error and allows a modified log-likelihood calculation.

				if (param_l_0 == 0 || param_l_0 == 1){ l_density_0 = 0; }
				else{
					l_density_0 = (OBS[clust_l - 1])*log(param_l_0) + (POS[clust_l - 1] - OBS[clust_l - 1])*log(1 - param_l_0);
				}
				if (param_j_0 == 0 || param_j_0 == 1){ j_density_0 = 0; }
				else{
					j_density_0 = (OBS[clust_j - 1])*log(param_j_0) + (POS[clust_j - 1] - OBS[clust_j - 1])*log(1 - param_j_0);
				}
				if (param_bw_0 == 0 || param_bw_0 == 1){ bw_density_0 = 0; }
				else{
					bw_density_0 = (OBS[k])*log(param_bw_0) + (POS[k] - OBS[k])*log(1 - param_bw_0);
				}
				l_0 = l_density_0 + j_density_0 + bw_density_0;
	
//Step: Calculate abbreviated loglikelihood of the proposed membership (Binomial)	
			//Get numerator of change-in-loglik equation [l,j,bw for proposed clustering]
			//Zero out Z vector. Best practice for runs > 1
			for (int i = 0; i < Z.size(); i++){ Z[i] = 0; }
			//Populate total weights of links connected to l*, by cluster (These ares the Z_ values in the paper)
			//links[l_star][j].first is the vertex number of l*'s neighbor, 0~(N-1)
			//cluster_membership[links[l_star][j].first] is the neighbor's cluster #, 1~k
			//Z[cluster_membership[links[l_star][j].first] - 1] is the spot in Z-vector for holding the neighbor's cluster's weights to l*
			for (int j = 0; j < links[l_star].size(); j++){
				Z[cluster_membership[links[l_star][j].first] - 1] = Z[cluster_membership[links[l_star][j].first] - 1] + links[l_star][j].second;
			}
			//cout << "Sum of Edges (Weights) Connected to Vector " << l_star << ", by group" << endl;
			//for (int i = 0; i < k; i++){ cout << "Z_" << i + 1 << ": " << Z[i] << endl; }

			//Calculate parameters for proposed clustering
			Z_l = Z[clust_l - 1];
			Z_j = Z[clust_j - 1];
			n_l = group_size[clust_l - 1];
			n_j = group_size[clust_j - 1];

			//Equations for updating OBS and POS given a reclustering, taken from paper
			OBS_l_1 = OBS[clust_l - 1] - Z_l;
			OBS_j_1 = OBS[clust_j - 1] + Z_j;
			OBS_bw_1 = OBS[k] + Z_l - Z_j;

			POS_l_1 = ((n_l - 1)*(n_l - 2)) / 2;
			POS_j_1 = ((n_j + 1)*n_j) / 2;
			POS_bw_1 = (POS[k] + n_l - n_j - 1);

			param_l_1 = OBS_l_1 / (long double)POS_l_1;
			param_j_1 = OBS_j_1 / (long double)POS_j_1;
			param_bw_1 = OBS_bw_1 / (long double)POS_bw_1;

			l_density_1 = 0;
			j_density_1 = 0;
			bw_density_1 = 0;

			//Calculate l,j,bw densities for PROPOSED clustering
			//if unweighted, all edges have weight=1 and Binomial density is appropriate
				//Calculate modified new log-lik, the numerator in the change-in-loglik formula
			if (param_l_1 == 0 || param_l_1 == 1) { l_density_1 = 0;  }
				else{
					l_density_1 = (OBS_l_1)*log(param_l_1) + (POS_l_1 - OBS_l_1)*log(1 - param_l_1);
				}
				if (param_j_1 == 0 || param_j_1 == 1){ j_density_1 = 0; }
				else{
					j_density_1 = (OBS_j_1)*log(param_j_1) + (POS_j_1 - OBS_j_1)*log(1 - param_j_1);
				}
				if (param_bw_1 == 0 || param_bw_1 == 1){ bw_density_1 = 0; }
				else{
					bw_density_1 = (OBS_bw_1)*log(param_bw_1) + (POS_bw_1 - OBS_bw_1)*log(1 - param_bw_1);
				}
				l_1 = l_density_1 + j_density_1 + bw_density_1;			

//Step: Calculate change in log-likelihood (Binomial) due to reclustering
			//log(new/old)=log(new)-log(old)

			//If all densities in either the new or the old likelihood were discarded (and thus set to zero),
			//it is not possible to evaluate a change-in-likelihood
			//set entire change-in-likelihood to zero so a new clustering will be considered
			if (l_0 == 0 || l_1 == 0){ delta_l = 0; }
			else{ delta_l = l_1 - l_0; }
			//if delta_l>0, then l_1>l_0, and the reclustered membership is more "likely"
			//Otherwise, l_1<l_0 and the current membership is more "likely"
			
			int STAGE = 0; //Will indicate when/if a proposed recustering was accepted or rejected.

			if (delta_l > epsilon)
			{//ACCEPT PROPOSED CLUSTERING AS CURRENT (best) CLUSTERING
//Step: Compare the two likelihoods. If accept new likelihood, then:
				//Update current membership vector to reflect l_star's new cluster membership
				//Update 3 elements (each) of OBS/POSvectors 
				//Update group (cluster) size vector
				//Proposed loglikelihood becomes current loglikelihood.	
				STAGE = 1;

				//current clustering <- proposed clustering
				cluster_membership[l_star] = clust_j;

				//Update OBS vector
				OBS[clust_l - 1] = OBS_l_1;
				OBS[clust_j - 1] = OBS_j_1;
				OBS[k] = OBS_bw_1;

				//Update POS vector
				POS[clust_l - 1] = POS_l_1;
				POS[clust_j - 1] = POS_j_1;
				POS[k] = POS_bw_1;

				//Update group size vector
				group_size[clust_l - 1]--; //cluster l loses 1 vertex
				group_size[clust_j - 1]++; //cluster j gains 1 vertex

				l_0 = l_1;
				l_1 = 0;

				Success_counter++; //Iterate number of success at this temperature
			}
			else
			{
				long double one = 1;      //used in MIN call to match data types.

				double uni_draw = rand() / double(RAND_MAX);// dist_uni(e3);// dist_uni(engine)/*uniform_rand(rand())*/;
				double min_val = min(one, exp(delta_l / IT));
				
				if (uni_draw < min_val) //Accept proposed clustering
				{
//Step: If accept new loglikelihood, then:
					//Update current membership vector to reflect l_star's new cluster membership
					//Update 3 elements (each) of OBS/POSvectors 
					//Update group (cluster) size vector
					//Proposed loglikelihood becomes current loglikelihood.	
					STAGE = 2;
					//current clustering <- proposed clustering
					cluster_membership[l_star] = clust_j;

					//Update OBS vector
					OBS[clust_l - 1] = OBS_l_1;
					OBS[clust_j - 1] = OBS_j_1;
					OBS[k] = OBS_bw_1;

					//Update POS vector
					POS[clust_l - 1] = POS_l_1;
					POS[clust_j - 1] = POS_j_1;
					POS[k] = POS_bw_1;

					//Update group size vector
					group_size[clust_l - 1]--; //cluster l loses 1 vertex
					group_size[clust_j - 1]++; //cluster j gains 1 vertex

					l_0 = l_1;
					l_1 = 0;

					//Iterate number of success at this temperature
					Success_counter++; 
				}
				else //Reject proposed clustering
				{//Nothing happens here. The current cluster is retained and the process repeats itself.					
					STAGE = 3;
				}
			}
			//At this point, the proposed clustering has:
			//1. Been Accepted. membership vector/OBS/POS/cluster counts have beeen updated to reflect new clustering
			//2. Rejected. memberhsip vector/OBS/POS/cluster counts still correspond to the current clustering.
		}//End SA LOOP

//CALCULATE TIME REQUIRED TO REACH SOLUTION
		//cout << fixed;
		cout << setprecision(10);
		clock_t end_new = clock();
		double time = (double)(end_new - start_new) / CLOCKS_PER_SEC;/*(double)(end_new - start_new) * 1000.0 / CLOCKS_PER_SEC  for milliseconds*/

		//At this point, a solution clustering has been identified.
		long double full_loglik_curr = 0; //loglikelihood for the full clustering
		long double LO = 0;				//loglikelihood under null hypothesis: # of clusters=1
		long double BIC = 0;				//Bayesian Information Criterion
		long double TESTSTAT = 0;

//STEP: USING OBS and POS vectors, calculate full loglikelihood of "best" solutions
		double l_0_w = 0; //loglikelihood of within-clusters
		double l_0_b = 0; //loglikelihood of between-clusters
		//Get loglikelihood of within-cluster edges
		for (int i = 0; i < k; i++) {
			if (OBS[i] != 0 && (OBS[i] != POS[i])) {
				l_0_w = l_0_w +
					(
						(OBS[i])*log(OBS[i] / (long double)POS[i]) + (POS[i] - OBS[i])*log(1 - OBS[i] / (long double)POS[i])
					);
			}
		}

//Checking the observed and possible edge counts in the final solution
//for (int i = 0; i <= k; i++) { cout << OBS[i] << "..."; } cout << endl;
//for (int i = 0; i <= k; i++) { cout << POS[i] << "..."; } cout << endl;

		//Get loglikelihood of between-cluster edges
		unsigned long long int TOT_OBS = 0;
		unsigned long long int TOT_POS = 0;
		TOT_OBS = OBS[k];
		TOT_POS = POS[k];

		//Loglikelihood of between cluster edges when modeled with a single parameter
		l_0_b = (TOT_OBS)*log(TOT_OBS / (long double)TOT_POS) + ((long double)TOT_POS - TOT_OBS)*log(1 - (TOT_OBS / (long double)TOT_POS));
		full_loglik_curr = l_0_w + l_0_b;

//STEP: Calculate likelihood under null hypothesis of k=1 cluster
		long long int LO_OBS = 0;
		long long int LO_POS = 0;
		for (int i = 0; i <= k; i++)
		{
				LO_OBS = LO_OBS + OBS[i];			
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
			a_m[i] = a_m[i] * a_m[i]; //square, per Modularity formula
		}
		for (int i = 0; i < k; i++)
		{
			//total squared fraction of between-cluster edges
			a_mm = a_mm + a_m[i];
		}
		mod_0 = e_mm - a_mm; //modularity of solution membership assignment

//Checking the null hypothesis stats
//cout << "LO=" << LO << ", LO_OBS=" << LO_OBS << ", LO_POS=" << LO_POS << endl;

//STEP: Calculate Bayesian Information Criterion for solution clustering
		//In old code, BIC = -2 * (LO + 0.5*(-bestfit)) + (k + 1)*log(N*(N - 1) / 2);
		//but bestfit = -(-2 * (LO - full_loglik_curr)), so this can be simplified:		
		//number of parameters=(k+1), number of observations=Nchoose2=N(N-1)/2, maximized loglikelihood=full_loglik_curr		
		BIC = -2 * (full_loglik_curr)+(k + 1)*log(N*(N - 1) / (long double)2);
		TESTSTAT = -2 * (LO - full_loglik_curr);

//cout << "at k=" << k << ", modularity= " << mod_0 << ", loglik=" << full_loglik_curr << " and BIC=" << BIC << ", TESTSTAT=" << TESTSTAT << ", time=" << time << " sec.s" << endl;;
		//Identify the optimum number of clusters, based on minimum Bayesian Information Criterion
		if (k == min_k){
			final_k = k;
			final_BIC = BIC;
			final_LOGLIK = full_loglik_curr;
			final_TESTSTAT = TESTSTAT;
			final_membership = cluster_membership;
			final_MODULARITY = mod_0;
			final_TIME = time;
		}
		//if BIC < min BIC thus far, update variables
		if (BIC < final_BIC){
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
		//remove("completion_lik_efficient_all_sols.txt");
		atkstring << "completion_lik_efficient_all_sols" << ".txt";
		atkinfo.open("completion_lik_efficient_all_sols.txt", ios::app);
		atkinfo << network_key << "," << N << "," << nb_links << "," << k << "," << mod_0 << "," << full_loglik_curr << "," << BIC << "," << TESTSTAT << "," << time  << endl;
		atkinfo.close();


		//cout << "**********************************" << endl;
	}//end for-k loop

//STEP: Output best solution across all k
	ostringstream atkstring2;
	ofstream atkinfo2;
	// If old file exists from a previous run, delete it.
	//remove("completion_lik_efficient_final_sol.txt");
	atkstring2 << "completion_lik_efficient_final_sol" << ".txt";
	atkinfo2.open("completion_lik_efficient_final_sol.txt", ios::app);
	atkinfo2 << network_key << "," << N << "," << nb_links << "," << final_k << "," << final_MODULARITY << "," << final_LOGLIK << "," << final_BIC << "," << final_TESTSTAT << "," << final_TIME << endl;
	atkinfo2.close();

//STEP: Output best solution clustering
	ostringstream bestzsol;
	ofstream bestzstorage;
	// If old file exists from a previous run, delete it.
	//remove("solution_lik_efficient.txt");
	bestzsol << "solution_lik_efficient" << ".txt";
	bestzstorage.open("solution_lik_efficient.txt", ios::app);
	for (int i = 0; i<N; i++){ bestzstorage << final_membership[i];
	if (i < N - 1) { bestzstorage << ","; }
	}
	bestzstorage << endl;
	bestzstorage.close();


	return 0;
}

