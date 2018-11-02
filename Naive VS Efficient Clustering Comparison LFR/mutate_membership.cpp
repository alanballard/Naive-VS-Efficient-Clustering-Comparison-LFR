//#include <algorithm>
//#include <time.h>
#include <iostream>
#include <vector>
#include <random>
using namespace std;

/*This program mutates the cluster membership of a given vertex subject to the following condition:
Each of k clusters must have a minimum of two members after mutation*/

//current_membership: N-element vector containing current group membership of all vertices in Network
//group_size: K-element vector containing number of vertices belonging to each group
vector<int> mutate_membership(vector<int> current_membership, vector<int> group_size)
{	
//cout << "Entering Mutate_membership" << endl;
	//random_device rd;
	//mt19937 e5(rd());
	//uniform_int_distribution<int> dist_vert(0, current_membership.size() - 1);/*Select vertex # 0-((total_number_of_verts=N)-1)*/
	//uniform_int_distribution<int> dist_clust(1, group_size.size());/*select cluster # 1-(total_number_of_clusters=k)*/

	//#1. select a vertex for mutation
	//#2. select a cluster to which to reassign the vertex
	//#3. verify that all clusters still have at least 2 members post-reassignment
	//no? Return to #1.
	//yes? Bind #1 and #2 into a vector and return to Main.

	int vert_prop;	//Vertex selected for reassignment. Int 0 to (N-1)
	int clust_curr; //Original cluster membership of selected vertex. Int 1 to K	
	int clust_prop;	//Cluster selected for reassignment. Int 1 to K;

	//Randomly select a vertex for reassignment, and a cluster to reassign the vertex to.
	vert_prop = rand() % current_membership.size(); // dist_vert(e5);/*dist_vert(engine) rand() % current_membership.size(); */
		clust_prop = rand() % group_size.size() + 1;  //dist_clust(e5);/*dist_clust(engine) rand() % group_size.size() + 1;	*/
//	cout << "V,C=" << vert_prop << ", " << clust_prop << "...";
	//Vertex should be reassigned to a cluster OTHER THAN the cluster to which it is currently assigned. 
	//If they are the same, choose another cluster to reassign the selected vertex.
	while (current_membership[vert_prop] == clust_prop){	
		clust_prop = rand() % group_size.size() + 1;  //dist_clust(e5);/*dist_clust(engine) rand() % group_size.size() + 1;*/
	}

	clust_curr = current_membership[vert_prop];

	//# of times that the reclustering fails to meet out requirements. Used to quit.
	int failure_count = 0;

	//Assume we start with a valid cluster membership vector.
	//If we remove a vector from one cluster and reassign to another, the only cluster that is possibly in danger of 
	//having an insufficient number of members is the initial cluster that lost a member. 
	//Thus, that is the only check required and can be done by subtracting 1  from current cluster size.
	//If the cluster is too small post-reassignment, select a new vertex and cluster for reassignment.
	while ( (group_size[(clust_curr-1)] - 1) <2) 
	{
		failure_count++;
		vert_prop = rand() % current_membership.size(); //dist_vert(e5); /*dist_vert(engine) rand() % current_membership.size();*/
		clust_prop = rand() % group_size.size() + 1; //dist_clust(e5); /*dist_clust(engine) rand() % group_size.size() + 1;*/

		while (current_membership[vert_prop] == clust_prop){	//Do not reassign to currently-assigned cluster
			clust_prop = rand() % group_size.size() + 1;  //dist_clust(e5); /*dist_clust(engine) rand() % group_size.size() + 1;*/			//Random int 0 to (K-1);
		}

		clust_curr = current_membership[vert_prop];
		//cout << failure_count << endl;
		if (failure_count == 10000){ 
			cout << "CANNOT MEET RECLUSTERING REQUIREMENTS" << endl;
			cout << "SEE MUTATE_MEMBERSHIP.cpp"; cin.get(); cin.get();
		}
	}

	vector<int> reassignment;
	reassignment.push_back(vert_prop); //Vertex to be reassigned
	reassignment.push_back(clust_prop);//cluster to which it will be reassigned
	//cout << "Vector to be reassigned: " << reassignment[0] << endl;
	//cout << "Cluster to which it will be reassigned: " << reassignment[1] << endl;
//cout << "Exiting Mutate_membership" << endl;
	return reassignment;
}