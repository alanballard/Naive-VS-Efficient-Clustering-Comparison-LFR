// File: louvain.cpp
// -- community detection source file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// And based on the article "A Generalized and Adaptive Method for Community Detection"
// Copyright (C) 2014 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume
//
// This file is part of Louvain algorithm.
// 
// Louvain algorithm is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Louvain algorithm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : July 2014
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "louvain.h"

using namespace std;


Louvain::Louvain(int nbp, long double epsq, Quality* q) {
  qual = q;

  neigh_weight.resize(qual->size,-1);
  neigh_pos.resize(qual->size);
  neigh_last = 0;

  nb_pass = nbp;
  eps_impr = epsq;
}

void
Louvain::init_partition(char * filename) {
  ifstream finput;
  finput.open(filename,fstream::in);

  // read partition
  while (!finput.eof()) {
    int node, comm;
    finput >> node >> comm;
    
    if (finput) {
      int old_comm = qual->n2c[node];
      neigh_comm(node);

      qual->remove(node, old_comm, neigh_weight[old_comm]);
      
      int i=0;
      for (i=0 ; i<neigh_last ; i++) {
	int best_comm = neigh_pos[i];
	long double best_nblinks = neigh_weight[neigh_pos[i]];
	if (best_comm==comm) {
	  qual->insert(node, best_comm, best_nblinks);
	  break;
	}
      }
      if (i==neigh_last)
	qual->insert(node, comm, 0);
    }
  }
  finput.close();
}

void
Louvain::neigh_comm(int node) {
  for (int i=0 ; i<neigh_last ; i++)
    neigh_weight[neigh_pos[i]]=-1;
  
  neigh_last = 0;

  pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(node);
  int deg = (qual->g).nb_neighbors(node);

  neigh_pos[0] = qual->n2c[node];
  neigh_weight[neigh_pos[0]] = 0;
  neigh_last = 1;

  for (int i=0 ; i<deg ; i++) {
    int neigh  = *(p.first+i);
    int neigh_comm = qual->n2c[neigh];
    long double neigh_w = ((qual->g).weights.size()==0)?1.0L:*(p.second+i);
    
    if (neigh!=node) {
      if (neigh_weight[neigh_comm]==-1) {
	neigh_weight[neigh_comm] = 0.0L;
	neigh_pos[neigh_last++] = neigh_comm;
      }
      neigh_weight[neigh_comm] += neigh_w;
    }
  }
}

void
Louvain::partition2graph() {
  vector<int> renumber(qual->size, -1);
  for (int node=0 ; node<qual->size ; node++) {
    renumber[qual->n2c[node]]++;
  }

  int end=0;
  for (int i=0 ; i< qual->size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=end++;

  for (int i=0 ; i< qual->size ; i++) {
    pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(i);

    int deg = (qual->g).nb_neighbors(i);
    for (int j=0 ; j<deg ; j++) {
      int neigh = *(p.first+j);
      cout << renumber[qual->n2c[i]] << " " << renumber[qual->n2c[neigh]] << endl;
    }
  }
}

void
Louvain::display_partition() {
  vector<int> renumber(qual->size, -1);
  for (int node=0 ; node < qual->size ; node++) {
    renumber[qual->n2c[node]]++;
  }

  int end=0;
  for (int i=0 ; i < qual->size ; i++)
    if (renumber[i]!=-1)
      renumber[i] = end++;

  for (int i=0 ; i < qual->size ; i++)
    cout << i << " " << renumber[qual->n2c[i]] << endl;
}

//ALAN
vector<int>
Louvain::store_solution(vector<int> current_membership) {
	
	vector<int> updated_membership = current_membership;

	vector<int> renumber(qual->size, -1);
	for (int node = 0; node < qual->size; node++) {
		renumber[qual->n2c[node]]++;
	}

	int end = 0;
	for (int i = 0; i < qual->size; i++)
		if (renumber[i] != -1)
			renumber[i] = end++;

	for (int i = 0; i < qual->size; i++) {
		//cout << i << " " << renumber[qual->n2c[i]] << endl;
		for (int j = 0; j < current_membership.size(); j++) {
			if (current_membership[j] == i) {
				updated_membership[j] = renumber[qual->n2c[i]];
			}
		}
	}

	return updated_membership;
}
//ALAN
void
Louvain::final_output(int network_key, vector< vector<pair<int, int>> > links, vector<int> current_membership, long double quality, int nb_links, long double time) {
	//for (int i = 0; i < current_membership.size(); i++) { cout << i << ", " << current_membership[i] << endl; }
	//cout << quality << endl;

	int final_k = 0;
	for (int i = 0; i < current_membership.size(); i++) { current_membership[i] = current_membership[i] + 1; }//Adjusting cluster numbers to start from 1, rather than zero, to match other objective fct. output
	for (int i = 0; i < current_membership.size(); i++) { if (current_membership[i] > final_k) { final_k = current_membership[i]; } }
	vector<int> group_size(final_k, 0);		//used in initial_pop to store number of nodes in each group 1:final_k
	for (int i = 0; i < current_membership.size(); i++) { 
		group_size[current_membership[i]-1]++; 
	}

//***********************************************
	//In the case of modularity, we need to populate the OBS/POS matrices for the final membership here in order to calculate the likelihood for the modularity solution, so we can calculate BIC and a test statistic.
	vector< vector <unsigned long long int >> OBS(final_k, vector<unsigned long long int>(final_k, 0));
	vector< vector <unsigned long long int >> POS(final_k, vector<unsigned long long int>(final_k, 0));


	for (int i = 0; i < links.size(); i++)
	{
		for (int j = 0; j < links[i].size(); j++)
		{
			//cout << "i=" << i << ", j=" << links[i][j].first << ", ct. of i's nbrs=" << links[i].size() << ", i's clust=" << current_membership[i] << ", j's cluster=" << current_membership[links[i][j].first] << ", weight=" << links[i][j].second;
			if (current_membership[i] == current_membership[links[i][j].first])
			{
				OBS[current_membership[i] - 1]
					[current_membership[i] - 1] =
					OBS[current_membership[i] - 1]
					[current_membership[i] - 1] + links[i][j].second;

			}
			else
			{
				OBS[current_membership[i] - 1]
					[current_membership[links[i][j].first] - 1] =
					OBS[current_membership[i] - 1]
					[current_membership[links[i][j].first] - 1] + links[i][j].second;
				//And to the corresponding lower triangle
				//OBS[current_membership[j] - 1][current_membership[i] - 1] = OBS[current_membership[j] - 1][current_membership[i] - 1] + links[i][j].second;
			}
		}//j'th neighbor to i'th vertex
	}//i'th vertex in adjacency list

	 //ADDRESS DOUBLE-COUNTS
	 //Note: Since the adjacency list contains 2 entries for every edge and the number of observed edges are counted
	 //over the adjacency list, the total observed count needs to be divided by 2.
	for (int i = 0; i < final_k; i++) { OBS[i][i] = OBS[i][i] / 2; }

	//Populate possible edge count vector. This is done using the cluster size information
	for (int i = 0; i < final_k; i++)
	{
		POS[i][i] = ((unsigned long long int) group_size[i] * (group_size[i] - 1)) / 2;
	}
	//Since we can use j=i+1 here, we don't have double counts.
	for (int i = 0; i < final_k - 1; i++) {
		for (int j = i + 1; j < final_k; j++)
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
	for (int i = 0; i < final_k; i++) {
		if (OBS[i][i] != 0 && (OBS[i][i] != POS[i][i])) {
			l_0_w = l_0_w +
				(OBS[i][i])*log(OBS[i][i] / (long double)POS[i][i]) + (POS[i][i] - OBS[i][i])*log(1 - OBS[i][i] / (long double)POS[i][i]);
		}
	}

	//Get loglikelihood of between-cluster edges
	unsigned long long int TOT_OBS = 0;
	unsigned long long int TOT_POS = 0;
	for (int i = 0; i < final_k; i++) {
		for (int j = i + 1; j < final_k; j++) {

			if (OBS[i][j] != 0) {
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
	for (int i = 0; i < final_k; i++)
	{
		for (int j = i; j < final_k; j++)
		{
			LO_OBS = LO_OBS + OBS[i][j];
			//cout << "OBS[i][j]="<<OBS[i][j] <<", LO_OBS="<< LO_OBS << endl;
		}
	}
	LO_POS = current_membership.size() * (current_membership.size() - 1) / 2; //N choose 2

	LO = (LO_OBS)*log(LO_OBS / (long double)LO_POS) + (LO_POS - LO_OBS)*log(1 - (LO_OBS / (long double)LO_POS));

	//Indicates how many, if any, of the clusters have a observed edge weight/count=0.
	int flag = 0;
	for (int i = 0; i < final_k; i++) {
		for (int j = i; j < final_k; j++) {
			if (OBS[i][j] == 0) { flag++; }
		}
	}

	//Calculate Bayesian Information Criterion for solution clustering
	//In old code, BIC = -2 * (LO + 0.5*(-bestfit)) + (final_k + 1)*log(N*(N - 1) / 2);
	//but bestfit = -(-2 * (LO - full_loglik_curr)), so this can be simplified:		
	//number of parameters=(final_k+1), number of observations=Nchoose2=N(N-1)/2, maximized loglikelihood=full_loglik_curr
//NOTE: In the case of modularity, full_loglik_curr is actually just the likelihood of the modularity solution, which may be different than the solution that maximizes likelihood
//Therefore the BIC, loglikelihood and test statistic for the solution clustering under modularity are not particularly meaningful. Just calculating here to be thorough.
	BIC = -2 * (full_loglik_curr)+(final_k + 1)*log(current_membership.size()*(current_membership.size() - 1) / (long double)2);
	TESTSTAT = -2 * (LO - full_loglik_curr);
//************************************************
//	cout << "final k= " <<final_k <<", modularity= " << quality << ", loglik=" << full_loglik_curr << ", BIC=" << BIC << ", TESTAT=" << TESTSTAT << ", time=" << time << " sec.s" << endl;

	ostringstream atkstring2;
	ofstream atkinfo2;
	// If old file exists from a previous run, delete it.
	//remove("completion_louv_final_sol.txt");
	atkstring2 << "completion_louv_final_sol" << ".txt";
	atkinfo2.open("completion_louv_final_sol.txt", ios::app);

	//atkinfo2 << network_key << "," << N << "," << nb_links << "," << final_k << "," << final_MODULARITY << "," << final_LOGLIK << "," << final_BIC << "," << final_TESTSTAT << "," << final_TIME;
	atkinfo2 << network_key << "," << current_membership.size() << "," << nb_links << "," << final_k << "," << quality << "," << full_loglik_curr << "," << BIC << "," << TESTSTAT << "," << time;
	atkinfo2 << endl;
	atkinfo2.close();

	//PRINT OUT LOUVAIN SOLUTION TO FILE
	ostringstream bestzsol;
	ofstream bestzstorage;
	// If old file exists from a previous run, delete it.
	//remove("solution_louv.txt");
	bestzsol << "solution_louv" << ".txt";
	bestzstorage.open("solution_louv.txt", ios::app);
	for (int i = 0; i < current_membership.size(); i++) { bestzstorage << current_membership[i] << "\t"; }
	bestzstorage << endl;
	bestzstorage.close();
}

Graph
Louvain::partition2graph_binary() {
  // Renumber communities
  vector<int> renumber(qual->size, -1);
  for (int node=0 ; node < qual->size ; node++)
    renumber[qual->n2c[node]]++;

  int last=0;
  for (int i=0 ; i < qual->size ; i++) {
    if (renumber[i]!=-1)
      renumber[i] = last++;
  }
  
  // Compute communities
  vector<vector<int> > comm_nodes(last);
  vector<int> comm_weight(last, 0);
  
  for (int node = 0 ; node < (qual->size) ; node++) {
    comm_nodes[renumber[qual->n2c[node]]].push_back(node);
    comm_weight[renumber[qual->n2c[node]]] += (qual->g).nodes_w[node];
  }

  // Compute weighted graph
  Graph g2;
  int nbc = comm_nodes.size();

  g2.nb_nodes = comm_nodes.size();
  g2.degrees.resize(nbc);
  g2.nodes_w.resize(nbc);
  
  for (int comm=0 ; comm<nbc ; comm++) {
    map<int,long double> m;
    map<int,long double>::iterator it;

    int size_c = comm_nodes[comm].size();

    g2.assign_weight(comm, comm_weight[comm]);

    for (int node=0 ; node<size_c ; node++) {
      pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(comm_nodes[comm][node]);
      int deg = (qual->g).nb_neighbors(comm_nodes[comm][node]);
      for (int i=0 ; i<deg ; i++) {
	int neigh = *(p.first+i);
	int neigh_comm = renumber[qual->n2c[neigh]];
	long double neigh_weight = ((qual->g).weights.size()==0)?1.0L:*(p.second+i);

	it = m.find(neigh_comm);
	if (it==m.end())
	  m.insert(make_pair(neigh_comm, neigh_weight));
	else
	  it->second += neigh_weight;
      }
    }

    g2.degrees[comm] = (comm==0)?m.size():g2.degrees[comm-1]+m.size();
    g2.nb_links += m.size();

    for (it = m.begin() ; it!=m.end() ; it++) {
      g2.total_weight += it->second;
      g2.links.push_back(it->first);
      g2.weights.push_back(it->second);
    }
  }

  return g2;
}

bool
Louvain::one_level() {
  bool improvement=false ;
  int nb_moves;
  int nb_pass_done = 0;
  long double new_qual = qual->quality();
  long double cur_qual = new_qual;

  vector<int> random_order(qual->size);
  for (int i=0 ; i < qual->size ; i++)
    random_order[i]=i;
  for (int i=0 ; i < qual->size-1 ; i++) {
    int rand_pos = rand()%(qual->size-i)+i;
    int tmp = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }

  // repeat while 
  //   there is an improvement of quality
  //   or there is an improvement of quality greater than a given epsilon 
  //   or a predefined number of pass have been done
  do {
    cur_qual = new_qual;
    nb_moves = 0;
    nb_pass_done++;

    // for each node: remove the node from its community and insert it in the best community
    for (int node_tmp = 0 ; node_tmp < qual->size ; node_tmp++) {
      int node = random_order[node_tmp];
      int node_comm = qual->n2c[node];
      long double w_degree = (qual->g).weighted_degree(node);

      // computation of all neighboring communities of current node
      neigh_comm(node);
      // remove node from its current community
      qual->remove(node, node_comm, neigh_weight[node_comm]);

      // compute the nearest community for node
      // default choice for future insertion is the former community
      int best_comm = node_comm;
      long double best_nblinks  = 0.0L;
      long double best_increase = 0.0L;
      for (int i=0 ; i<neigh_last ; i++) {
	long double increase = qual->gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree);
	if (increase>best_increase) {
	  best_comm = neigh_pos[i];
	  best_nblinks = neigh_weight[neigh_pos[i]];
	  best_increase = increase;
	}
      }

      // insert node in the nearest community
      qual->insert(node, best_comm, best_nblinks);
     
      if (best_comm!=node_comm)
	nb_moves++;
    }

    new_qual = qual->quality();
    
    if (nb_moves>0)
      improvement=true;

  } while (nb_moves>0 && new_qual-cur_qual > eps_impr);

  return improvement;
}
