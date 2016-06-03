# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 11:15:38 2016

@author: jacquemont
"""

def log_size(x):
    
    import math
    
    y = -math.log10(1-x)
    
    return y

# size function_definition

def graph_size(G, size_type):
    
    import networkx as nx
    
    if size_type=='extent':
        size = nx.number_of_edges(G)
    elif size_type=='intensity':
        graph = nx.Graph(G)
        size = sum(list(nx.get_edge_attributes(graph, 'weight').values()))
    elif size_type=='log_intensity':
        graph = nx.Graph(G)
        size = sum(map(log_size, list(nx.get_edge_attributes(graph, 'weight').values())))
    else:
        raise IOError('size_type input should set as "extent", "intensity" or "log_intensity".')

    return size