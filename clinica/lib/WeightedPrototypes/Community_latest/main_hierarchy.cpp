// File: main_hierarchy.cpp
// -- output community structure handling (number of levels, communities of one level)
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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

int display_level = -1;
char *filename = NULL;

void
usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " input_file [options]" << endl << endl;
  cerr << "input_file: read the community tree from this file." << endl;
  cerr << "-l xx\t display the community structure for the level xx." << endl;
  cerr << "\t outputs the community for each node." << endl;
  cerr << "\t xx must belong to [-1,N] if N is the number of levels." << endl;
  cerr << "-n\t displays the number of levels and the size of each level." << endl;
  cerr << "\t equivalent to -l -1." << endl;
  cerr << "-h\tshow this usage message." << endl;
  exit(0);
}

void
parse_args(int argc, char **argv) {
  if (argc<2)
    usage(argv[0], "Bad arguments number\n");

  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 'l':
	display_level = atoi(argv[i+1]);
	i++;
	break;
      case 'n':
	display_level = -1;
	break;
      default:
	usage(argv[0], "Unknown option\n");
      }
    } else {
      if (filename==NULL)
        filename = argv[i];
      else
        usage(argv[0], "More than one filename\n");
    }
  }
  if (filename==NULL)
    usage(argv[0], "No input file has been provided.\n");
}
int
main(int argc, char **argv) {
  parse_args(argc, argv);

  vector<vector<int> >levels;

  ifstream finput;
  finput.open(filename,fstream::in);

  int l=-1;
  while (!finput.eof()) {
    int node, nodecomm;
    finput >> node >> nodecomm;

    if (finput) {
      if (node==0) {
	l++;
	levels.resize(l+1);
      }
      levels[l].push_back(nodecomm);
    }
  }

  if (display_level==-1) {
    cout << "Number of levels: " << levels.size() << endl;
    for (unsigned int i=0 ; i<levels.size();i++)
      cout << "level " << i << ": " << levels[i].size() << " nodes" << endl;
  } else if (display_level<0 || (unsigned)display_level>=levels.size()) {
    cerr << "Incorrect level\n";
  } else {
    vector<int> n2c(levels[0].size());

    for (unsigned int i=0 ; i<levels[0].size() ; i++)
      n2c[i]=i;
    
    for (l=0 ; l<display_level ; l++)
      for (unsigned int node=0 ; node<levels[0].size() ; node++)
	n2c[node] = levels[l][n2c[node]];
    
    for (unsigned int node=0 ; node<levels[0].size() ; node++)
      cout << node << " " << n2c[node] << endl;
  }
}
