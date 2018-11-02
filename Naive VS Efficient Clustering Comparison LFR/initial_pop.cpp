#include <algorithm>
#include <iostream>
#include <time.h>
#include <vector>

#include <random>;
using namespace std;

//FUNCTION: randomly assigns N nodes to k groups and stores the assignments in array[]
//          group[] contains the number of nodes in each group
int initial_pop(vector<int>& group_size, vector<int>& cluster_membership, int k, int N)
{
//cout << "Entering initial_pop" << endl;
	//random_device rd;
	//minstd_rand e6(rand());
	//uniform_int_distribution<int> dist2(1, k);

	//srand(time(0));
	bool check = true;         //true if any group has less than two node members
	//initially set to true so we can enter the while loop

	//srand(time(0));          //seeds random number generator
	
	//# of times that the reclustering fails to meet out requirements. Used to quit.
	int failure_count = 0;
	while (check == true)
	{
		check = false;
		for (int i = 0; i < N; i++) { cluster_membership.push_back(rand() % k + 1); } //cluster_membership.push_back(dist2(e6));/*dist2(engine) rand() % k + 1;*/ 
		//Checks assignment vector z and counts how many nodes are in each group
		int v = 0;
		while (v<k)
		{
			
			v = v + 1;
			int mycount = count(cluster_membership.begin(), cluster_membership.end(), v);
			group_size[v - 1] = mycount;
			if (mycount <= 1 || mycount == N){ check = true; cluster_membership.clear(); }
		}
		if (check == true)//This means at least one of the cluster had insufficient member count
		{
			failure_count++;
		}
		//cout << failure_count << endl;
		if (failure_count == 10000){
			cout << "CANNOT MEET INTITIAL CLUSTERING REQUIREMENTS" << endl;
			cout << "SEE INITIAL_POP.cpp"; cin.get(); cin.get();
		}
	}
//for (int i = 0; i < cluster_membership.size(); i++) { cout << cluster_membership[i] << "-"; }cout << endl;
	return 0;
}
