#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;
//Reads in edge list from text and converts to adjacency list.
//Input  text shoud be set up as: vertex 1 vertex 2
//Vertex numbering must start at 0 and cannot skip numbers
vector< vector<pair<int, int>> > read_el(char *filename)
{

	char weighted = 'N';
	//cout << "Inputting text graph..." << endl;
	ifstream finput;
	finput.open(filename, fstream::in);
	

	int nb_links = 0;

	//vector<vector<int>> links;
	vector< vector<pair<int, int>> > links;
//**********************
//**********************
//Look at the file Graph.CPP on the thumb drive, C++ Code, Generic Louvain, src folder
//That shows how to convert this to a weighted network
//Need to change input format to test for presence of weights, or ask if a weighted network at the start
//If so, input format should be "src	dest	weight"
//If not, then "src	dest" and the weight should be filled in as one.
	while (!finput.eof()) {
		 int src, dest;		
		 int weight;

		 if (weighted == 'Y') //weighted case
			{
				finput >> src >> dest >> weight;
				if (finput) {
					//make room for each new src/dest combination ...
					if (links.size() <= max(src, dest) + 1) //if the size of the vector links is less than (the largest id# of the nodes+1)...
					{
						links.resize(max(src, dest) + 1); //resize the vector links to size (largest id# of the nodes+1) elements
					}
					
					links[src].push_back(make_pair(dest,weight)); //vec <src,dest>
					/*Note: The 2 lines below are necessary if the input edge list only contains one entry for each link.
					i.e., if 2 and 3 are connected, then the link 2 3 will appear in the edge list, but 3 2 will not (or vice versa).
					If 2 3 and 3 2 are present, then the 2 lines below will double the links in the created adjacency list.
					*/
					//			if (src != dest)
					//				links[dest].push_back(src);	//vec <src,dest>
					nb_links++; //increase number of links by one.
				}
			}
			else//unweighted case
			{
				finput >> src >> dest;
				
				if (finput) {
				//make room for each new src/dest combination ...
					if (links.size() <= max(src, dest) + 1) //if the size of the vector links is less than (the largest id# of the nodes+1)...
					{
						links.resize(max(src, dest) + 1); //resize the vector links to size (largest id# of the nodes+1) elements
					}

				links[src].push_back(make_pair(dest,1)); //vec <src,dest>
				/*Note: The 2 lines below are necessary if the input edge list only contains one entry for each link.
				i.e., if 2 and 3 are connected, then the link 2 3 will appear in the edge list, but 3 2 will not (or vice versa).
				If 2 3 and 3 2 are present, then the 2 lines below will double the links in the created adjacency list.
				*/
				//			if (src != dest)
				//				links[dest].push_back(src);	//vec <src,dest>
				nb_links++; //increase number of links by one.
				}
			}
	}
//	cout << "***The number of distinct vertices is: " << links.size() << "***" << endl;
//	cout << "Distinct edge count: " << nb_links/2 << endl;

	finput.close();

	return links;
	//cout << "Press Enter to Continue." << endl; cin.get();
}
