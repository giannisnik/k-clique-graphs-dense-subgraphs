# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 2016

@author: g.nikolentzos

"""

import networkx as nx

def density(G):
	""" 
	Returns the density of the graph
	density = 2*|E(S)|/(|S|*(|S|-1))
    """
	return 2*G.number_of_edges()/float(G.number_of_nodes()*(G.number_of_nodes()-1))


def degree_density(G):
	""" 
	Returns the degree density of the graph
	degree density = 2*|E(S)|/|S|
    """
	return 2*G.number_of_edges()/float(G.number_of_nodes())

	
def edge_surplus(G, alpha):
	""" 
	Returns the edge surplus of the graph
	edge surplus = |E(S)| - alpha*(|S|*(|S|-1))/2
    """
	return G.number_of_edges() - alpha * ((G.number_of_nodes() * (G.number_of_nodes() - 1))/2.0)
	

def triangle_density(G):
	""" 
	Returns the triangle density of the graph
	triangle density = |T(S)|/(|S|*(|S|-1)*(|S|-2))/6
    """
	t = nx.triangles(G)
	return (sum(t.values())/3.0)/((G.number_of_nodes()*(G.number_of_nodes()-1)*(G.number_of_nodes()-2))/6.0)
