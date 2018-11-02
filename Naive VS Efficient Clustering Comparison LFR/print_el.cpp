#include <iostream>
#include <iomanip>

#include <vector>
#include <algorithm>

using namespace std;

void print_el(vector< vector<pair< int, int>> > links)
{

	/*
		cout << "As an Edge List..." << endl;
		for (unsigned int i = 0; i < links.size(); i++) {				//links.size()=the number of edges in the list, assuming starting at zero. 
			//cout << "links.size(): " << links.size() << endl;
			for (unsigned int j = 0; j < links[i].size(); j++) {		//links[i].size()=number of connections to/from node #i so for 0-2, and 2-1, links[2].size()=2	
				//cout << "links[i].size(): " << links[i].size() << endl;
				int dest = links[i][j].first;
				int weight = links[i][j].second;

				cout << i << " " << dest << " <" << weight << ">" << endl; //I think this assumes that nodes are numbered starting with zero...
			}
		}
	*/	
		//Each row outputs the vertex,all the neighbors of that vetex and the weight of their connecting edges
		cout << "As an Adjacency List..." << endl;
		for (unsigned int i = 0; i < links.size(); i++) {
			cout << "Vertex " << i << ":";
			for (unsigned int j = 0; j < links[i].size(); j++) {
				cout << " " <<  links[i][j].first << " <"<< links[i][j].second <<">";
			}
			cout << endl;
		}
	
	//cout << "Press Enter to Continue..." << endl; cin.get();
}
