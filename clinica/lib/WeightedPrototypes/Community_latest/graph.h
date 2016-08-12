// File: graph.h
// -- simple graph handling header file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#ifndef GRAPH_H
#define GRAPH_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;

class Graph {
 public:
  vector<vector<pair<int,float> > > links;
  
  Graph (char *filename, int type);
  
  void clean(int type);
  void renumber(int type);
  void display(int type);
  void display_binary(char *filename, char *filename_w, int type);
};

#endif // GRAPH_H
