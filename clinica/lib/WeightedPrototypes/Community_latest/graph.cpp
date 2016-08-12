// File: graph.cpp
// -- simple graph handling source file
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

#include "graph.h"

using namespace std;

Graph::Graph(char *filename, int type) {
  ifstream finput;
  finput.open(filename,fstream::in);

  int nb_links=0;

  while (!finput.eof()) {
    unsigned int src, dest;
    double weight=1.;

    if (type==WEIGHTED) {
      finput >> src >> dest >> weight;
    } else {
      finput >> src >> dest;
    }
    
    if (finput) {
      if (links.size()<=max(src,dest)+1) {
        links.resize(max(src,dest)+1);
      }
      
      links[src].push_back(make_pair(dest,weight));
      if (src!=dest)
        links[dest].push_back(make_pair(src,weight));

      nb_links++;
    }
  }

  finput.close();
}

void
Graph::renumber(int type) {
  vector<int> linked(links.size(),-1);
  vector<int> renum(links.size(),-1);
  int nb=0;
  
  for (unsigned int i=0 ; i<links.size() ; i++) {
    for (unsigned int j=0 ; j<links[i].size() ; j++) {
      linked[i]=1;
      linked[links[i][j].first]=1;
    }
  }
  
  for (unsigned int i=0 ; i<links.size() ; i++) {
    if (linked[i]==1)
      renum[i]=nb++;
  }

  for (unsigned int i=0 ; i<links.size() ; i++) {
    if (linked[i]==1) {
      for (unsigned int j=0 ; j<links[i].size() ; j++) {
	links[i][j].first = renum[links[i][j].first];
      }
      links[renum[i]]=links[i];
    }
  }
  links.resize(nb);
}

void
Graph::clean(int type) {
  for (unsigned int i=0 ; i<links.size() ; i++) {
    map<int, float> m;
    map<int, float>::iterator it;

    for (unsigned int j=0 ; j<links[i].size() ; j++) {
      it = m.find(links[i][j].first);
      if (it==m.end())
	m.insert(make_pair(links[i][j].first, links[i][j].second));
      else if (type==WEIGHTED)
      	it->second+=links[i][j].second;
    }
    
    vector<pair<int,float> > v;
    for (it = m.begin() ; it!=m.end() ; it++)
      v.push_back(*it);
    links[i].clear();
    links[i]=v;
  }
}

void
Graph::display(int type) {
  for (unsigned int i=0 ; i<links.size() ; i++) {
    for (unsigned int j=0 ; j<links[i].size() ; j++) {
      int dest   = links[i][j].first;
      float weight = links[i][j].second;
      if (type==WEIGHTED)
	cout << i << " " << dest << " " << weight << endl;
      else
	cout << i << " " << dest << endl;
    }
  }
}

void
Graph::display_binary(char *filename, char *filename_w, int type) {
  ofstream foutput;
  foutput.open(filename, fstream::out | fstream::binary);

  unsigned int s = links.size();

  // outputs number of nodes
  foutput.write((char *)(&s),4);

  // outputs cumulative degree sequence
  long tot=0;
  for (unsigned int i=0 ; i<s ; i++) {
    tot+=(long)links[i].size();
    foutput.write((char *)(&tot),8);
  }

  // outputs links
  for (unsigned int i=0 ; i<s ; i++) {
    for (unsigned int j=0 ; j<links[i].size() ; j++) {
      int dest = links[i][j].first;
      foutput.write((char *)(&dest),4);
    }
  }
  foutput.close();

  // outputs weights in a separate file
  if (type==WEIGHTED) {
    ofstream foutput_w;
    foutput_w.open(filename_w,fstream::out | fstream::binary);
    for (unsigned int i=0 ; i<s ; i++) {
      for (unsigned int j=0 ; j<links[i].size() ; j++) {
	float weight = links[i][j].second;
	foutput_w.write((char *)(&weight),4);
      }
    }
    foutput_w.close();
  }
}

