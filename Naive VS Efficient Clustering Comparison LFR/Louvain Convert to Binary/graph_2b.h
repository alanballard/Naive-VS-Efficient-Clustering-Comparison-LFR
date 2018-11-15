// File: graph.h
// -- simple graph handling header file
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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
// see README.txt for more details


#ifndef GRAPH_H_2b
#define GRAPH_H_2b

#include <iostream>
#include <sstream> //ALAN: used to create outfile name in louvain_convert_2b
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;


class Graph_b {
 public:
  vector<vector<pair<int, int/*long double*/> > > links;
  
  Graph_b (vector<vector<pair<int, int/*long double*/> > > source_links, int link_count, char *filename, int type);
  
  void clean_b(int type);
  void renumber_b(int type, char *filename);
  void display_b(int type);
  void display_binary_b(char *filename, char *filename_w, int type);
};

#endif // GRAPH_H_2b
