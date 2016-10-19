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
	
	
def triangle_graph_density(G):
    """ 
    Returns the density of the triangle-graph created using the input graph
    """
    sum_min_degs = 0
    for node in G.nodes():
        edges = G.edges(node)
        r = 0
        g = 0
        b = 0
        for edge in edges:
            if G[edge[0]][edge[1]]['color'] == 'r':
                r += 1
            elif G[edge[0]][edge[1]]['color'] == 'g':
                g += 1
            elif G[edge[0]][edge[1]]['color'] == 'b':
                b += 1
            else:
                print 'error'
        sum_min_degs += min(r, g, b)

    return float(sum_min_degs) / G.number_of_nodes()
