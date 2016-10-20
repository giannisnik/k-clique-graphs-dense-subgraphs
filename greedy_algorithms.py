# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 2016

@author: g.nikolentzos

"""

import networkx as nx
import sys
from read import read_gml
from evaluation_metrics import (density,degree_density,triangle_density)

def greedy_degree_density(G):
	""" 
	Returns the subgraph with optimal degree density using 
	Charikar's greedy algorithm
    """
	neighbors=G.neighbors_iter
	degrees=G.degree()
	sum_degrees = sum(degrees.values())
	num_nodes = G.number_of_nodes()
	nodes=sorted(degrees,key=degrees.get)
	bin_boundaries=[0]
	curr_degree=0
	for i,v in enumerate(nodes):
		if degrees[v]>curr_degree:
			bin_boundaries.extend([i]*(degrees[v]-curr_degree))
			curr_degree=degrees[v]
	node_pos = dict((v,pos) for pos,v in enumerate(nodes))
	nbrs=dict((v,set(neighbors(v))) for v in G)
		
	max_degree_density = sum_degrees/float(num_nodes)
	ind = 0 
		
	for v in nodes:
		num_nodes -= 1
		while degrees[v] > 0:
			pos=node_pos[v]
			bin_start=bin_boundaries[degrees[v]]
			node_pos[v]=bin_start
			node_pos[nodes[bin_start]]=pos
			nodes[bin_start],nodes[pos]=nodes[pos],nodes[bin_start]
			bin_boundaries[degrees[v]]+=1
			degrees[v]-=1
		
		for u in nbrs[v]:
			
			nbrs[u].remove(v)
			pos=node_pos[u]
			bin_start=bin_boundaries[degrees[u]]
			node_pos[u]=bin_start
			node_pos[nodes[bin_start]]=pos
			nodes[bin_start],nodes[pos]=nodes[pos],nodes[bin_start]
			bin_boundaries[degrees[u]]+=1
			degrees[u]-=1
			sum_degrees -= 2
		
		
		if num_nodes > 0:
			current_degree_density = sum_degrees/float(num_nodes)
			if current_degree_density > max_degree_density:
				max_degree_density = current_degree_density
				ind = G.number_of_nodes()-num_nodes
				
	optimal_nodes = nodes[ind:]
		
	return G.subgraph(optimal_nodes)


def greedy_quasi_cliques(G, alpha):
	""" 
	Returns the subgraph with optimal edge surplus
    """
	neighbors=G.neighbors_iter
	degrees=G.degree()
	sum_degrees = sum(degrees.values())
	num_nodes = G.number_of_nodes()
	nodes=sorted(degrees,key=degrees.get)
	bin_boundaries=[0]
	curr_degree=0
	for i,v in enumerate(nodes):
		if degrees[v]>curr_degree:
			bin_boundaries.extend([i]*(degrees[v]-curr_degree))
			curr_degree=degrees[v]
	node_pos = dict((v,pos) for pos,v in enumerate(nodes))
	nbrs=dict((v,set(neighbors(v))) for v in G)
		
	max_edge_surplus = sum_degrees/2.0 - alpha * ((num_nodes * (num_nodes - 1))/2.0)
	ind = 0 
		
	for v in nodes:
		num_nodes -= 1
		while degrees[v] > 0:
			pos=node_pos[v]
			bin_start=bin_boundaries[degrees[v]]
			node_pos[v]=bin_start
			node_pos[nodes[bin_start]]=pos
			nodes[bin_start],nodes[pos]=nodes[pos],nodes[bin_start]
			bin_boundaries[degrees[v]]+=1
			degrees[v]-=1
		
		for u in nbrs[v]:
			
			nbrs[u].remove(v)
			pos=node_pos[u]
			bin_start=bin_boundaries[degrees[u]]
			node_pos[u]=bin_start
			node_pos[nodes[bin_start]]=pos
			nodes[bin_start],nodes[pos]=nodes[pos],nodes[bin_start]
			bin_boundaries[degrees[u]]+=1
			degrees[u]-=1
			sum_degrees -= 2
		
		
		if num_nodes > 0:
			current_edge_surplus = sum_degrees/2.0 - alpha * ((num_nodes * (num_nodes - 1))/2.0)
			if current_edge_surplus > max_edge_surplus:
				max_edge_surplus = current_edge_surplus
				ind = G.number_of_nodes()-num_nodes
				
	optimal_nodes = nodes[ind:]
		
	return G.subgraph(optimal_nodes)


def find_triangles(G):
	""" 
	For each node returns the number of triangles the node is part of
	and the pair of other nodes that form each of these triangles
    """
	triangles = {}
	nbrs = {}
	for node in G.nodes():
		triangles[node] = 0
		nbrs[node] = set()
		neighbors = G.neighbors(node)
		for i in range(len(neighbors)):
			for j in range(i+1,len(neighbors)):
				if G.has_edge(neighbors[i],neighbors[j]):
					triangles[node] += 1
					if neighbors[i] < neighbors[j]:
						nbrs[node].add((neighbors[i],neighbors[j]))
					else:
						nbrs[node].add((neighbors[j],neighbors[i]))
					
	return triangles, nbrs
	
	
def greedy_triangle_density(G):
	""" 
	Returns the subgraph with optimal triangle density
    """
	triangles, nbrs = find_triangles(G)
	sum_triangles = sum(triangles.values())
	num_nodes = G.number_of_nodes()
	nodes=sorted(triangles,key=triangles.get)
	bin_boundaries=[0]
	curr_triangle_number=0
	for i,v in enumerate(nodes):
		if triangles[v]>curr_triangle_number:
			bin_boundaries.extend([i]*(triangles[v]-curr_triangle_number))
			curr_triangle_number=triangles[v]
	node_pos = dict((v,pos) for pos,v in enumerate(nodes))
		
	max_triangle_density = float(sum_triangles)/num_nodes
	ind = 0 
		
	for v in nodes:

		num_nodes -= 1
		while triangles[v] > 0:
			pos=node_pos[v]
			bin_start=bin_boundaries[triangles[v]]
			node_pos[v]=bin_start
			node_pos[nodes[bin_start]]=pos
			nodes[bin_start],nodes[pos]=nodes[pos],nodes[bin_start]
			bin_boundaries[triangles[v]]+=1
			triangles[v]-=1
			sum_triangles -= 1
		
		for pair in nbrs[v]:
			
			if v < pair[1]:
				nbrs[pair[0]].remove((v,pair[1]))
			else:
				nbrs[pair[0]].remove((pair[1],v))

			pos=node_pos[pair[0]]
			bin_start=bin_boundaries[triangles[pair[0]]]
			node_pos[pair[0]]=bin_start
			node_pos[nodes[bin_start]]=pos
			nodes[bin_start],nodes[pos]=nodes[pos],nodes[bin_start]
			bin_boundaries[triangles[pair[0]]]+=1
			triangles[pair[0]]-=1
			sum_triangles -= 1

			if v < pair[0]:
				nbrs[pair[1]].remove((v,pair[0]))
			else:
				nbrs[pair[1]].remove((pair[0],v))

			pos=node_pos[pair[1]]
			bin_start=bin_boundaries[triangles[pair[1]]]
			node_pos[pair[1]]=bin_start
			node_pos[nodes[bin_start]]=pos
			nodes[bin_start],nodes[pos]=nodes[pos],nodes[bin_start]
			bin_boundaries[triangles[pair[1]]]+=1
			triangles[pair[1]]-=1
			sum_triangles -= 1		
		
		if num_nodes > 0:
			current_triangle_density = float(sum_triangles)/num_nodes
			if current_triangle_density > max_triangle_density:
				max_triangle_density = current_triangle_density
				ind = G.number_of_nodes()-num_nodes
				
	optimal_nodes = nodes[ind:]
		
	return G.subgraph(optimal_nodes)
	
	
def get_triangles(G):
	""" 
	Lists all the triangles of the graph
	"""
	triangles = {}
	edges = {}
	for edge in G.edges():
		if edge[0] < edge[1]:
			edges[(edge[0], edge[1])] = []
		else:
			edges[(edge[1], edge[0])] = []
	ind = 0
	done = set()
	for n in G:
		done.add(n)
		nbrdone = set()
		nbrs = set(G[n])
		for nbr in nbrs:
			if nbr in done:
				continue
			nbrdone.add(nbr)
			for both in nbrs.intersection(G[nbr]):
				if both in done or both in nbrdone:
					continue
				triangles[ind] = sorted((n, nbr, both))
				if n > nbr:
					edges[(nbr, n)].append(ind)
				else:
					edges[(n, nbr)].append(ind)
				if n > both:
					edges[(both, n)].append(ind)
				else:
					edges[(n, both)].append(ind)
				if both > nbr:
					edges[(nbr, both)].append(ind)
				else:
					edges[(both, nbr)].append(ind)
				ind += 1
	return triangles, edges


def generate_triangle_neighbors(triangles, edges):
	""" 
	For each triangle returns the triangles with which it shares a 
	common edge 
	"""
	neighbors = {}
	for triangle in triangles.keys():
		neighbors[triangle] = {}
		neighbors[triangle][(triangles[triangle][0], triangles[triangle][1])] = len(
			edges[(triangles[triangle][0], triangles[triangle][1])]) - 1
		neighbors[triangle][(triangles[triangle][0], triangles[triangle][2])] = len(
			edges[(triangles[triangle][0], triangles[triangle][2])]) - 1
		neighbors[triangle][(triangles[triangle][1], triangles[triangle][2])] = len(
			edges[(triangles[triangle][1], triangles[triangle][2])]) - 1

	return neighbors
    
    
def greedy_triangle_graph_density(G):
	""" 
	Returns the subgraph created by the subgraph of the triangle-graph
	with optimal triangle-graph density
    """
	triangles, edges = get_triangles(G)
	nbrs = generate_triangle_neighbors(triangles, edges)
	
	num_nodes = len(triangles)
	total_nodes = len(triangles)
	min_degs = {}
	for k1 in nbrs.keys():
		s = []
		for k2 in nbrs[k1].keys():
			s.append(nbrs[k1][k2])
		min_degs[k1] = min(s)
	sum_degs = sum(min_degs.values())
	nodes = sorted(min_degs, key=min_degs.get)
	bin_boundaries = [0]
	curr_deg = 0
	for i, v in enumerate(nodes):
		if min_degs[v] > curr_deg:
			bin_boundaries.extend([i] * (min_degs[v] - curr_deg))
			curr_deg = min_degs[v]
	node_pos = dict((v, pos) for pos, v in enumerate(nodes))

	max_degree_density = sum_degs / float(num_nodes)
	ind = 0

	for v in nodes:
		num_nodes -= 1
		sum_degs -= min_degs[v]

		while min_degs[v] > 0:
			pos = node_pos[v]
			bin_start = bin_boundaries[min_degs[v]]
			node_pos[v] = bin_start
			node_pos[nodes[bin_start]] = pos
			nodes[bin_start], nodes[pos] = nodes[pos], nodes[bin_start]
			bin_boundaries[min_degs[v]] += 1
			min_degs[v] -= 1

		edges[(triangles[v][0], triangles[v][1])].remove(v)
		for u in edges[(triangles[v][0], triangles[v][1])]:
			nbrs[u][(triangles[v][0], triangles[v][1])] -= 1
			if len(edges[(triangles[v][0], triangles[v][1])]) - 1 < min_degs[u]:
				pos = node_pos[u]
				bin_start = bin_boundaries[min_degs[u]]
				node_pos[u] = bin_start
				node_pos[nodes[bin_start]] = pos
				nodes[bin_start], nodes[pos] = nodes[pos], nodes[bin_start]
				bin_boundaries[min_degs[u]] += 1
				min_degs[u] -= 1
				sum_degs -= 1

		edges[(triangles[v][0], triangles[v][2])].remove(v)
		for u in edges[(triangles[v][0], triangles[v][2])]:
			nbrs[u][(triangles[v][0], triangles[v][2])] -= 1
			if len(edges[(triangles[v][0], triangles[v][2])]) - 1 < min_degs[u]:
				pos = node_pos[u]
				bin_start = bin_boundaries[min_degs[u]]
				node_pos[u] = bin_start
				node_pos[nodes[bin_start]] = pos
				nodes[bin_start], nodes[pos] = nodes[pos], nodes[bin_start]
				bin_boundaries[min_degs[u]] += 1
				min_degs[u] -= 1
				sum_degs -= 1

		edges[(triangles[v][1], triangles[v][2])].remove(v)
		for u in edges[(triangles[v][1], triangles[v][2])]:
			nbrs[u][(triangles[v][1], triangles[v][2])] -= 1
			if len(edges[(triangles[v][1], triangles[v][2])]) - 1 < min_degs[u]:
				pos = node_pos[u]
				bin_start = bin_boundaries[min_degs[u]]
				node_pos[u] = bin_start
				node_pos[nodes[bin_start]] = pos
				nodes[bin_start], nodes[pos] = nodes[pos], nodes[bin_start]
				bin_boundaries[min_degs[u]] += 1
				min_degs[u] -= 1
				sum_degs -= 1

		if num_nodes > 0:
			current_degree_density = sum_degs / float(num_nodes)
			if current_degree_density > max_degree_density:
				max_degree_density = current_degree_density
				ind = total_nodes - num_nodes

	optimal_nodes = nodes[ind:]

	nodes = set()

	for triangle in optimal_nodes:
		nodes.add(triangles[triangle][0])
		nodes.add(triangles[triangle][1])
		nodes.add(triangles[triangle][2])
		
	subg = G.subgraph(nodes)
	
	return subg

	
def main():
	""" 
	Main method
	"""
	filename = sys.argv[1]
	if filename.split('.')[1] == 'gml':
		G = read_gml('networks/' + filename)
	else:
		G = nx.read_edgelist('networks/' + filename, delimiter='\t', nodetype=int)

	G = G.to_undirected()

	for node in G.nodes_with_selfloops():
		G.remove_edge(node, node)

	G1 = nx.Graph()
	for edge in G.edges():
		u = edge[0]
		v = edge[1]
		if u == v:
			continue
		if not G1.has_edge(u, v):
			G1.add_edge(u, v, weight=1.0)

	G = G1
	print "Number of nodes:", G.number_of_nodes()
	print "Number of edges:", G.number_of_edges()
	print

	subg = greedy_degree_density(G)
	print "----Greedy Degree Density----"
	print "Degree Density: " + str(degree_density(subg))
	print "Density: " + str(density(subg))
	print "Triangle Density: " + str(triangle_density(subg))
	print "# Nodes: " + str(subg.number_of_nodes())
	print

	subg = greedy_quasi_cliques(G, 0.333)
	print "----Greedy Edge Surplus with alpha=1/3----"
	print "Degree Density: " + str(degree_density(subg))
	print "Density: " + str(density(subg))
	print "Triangle Density: " + str(triangle_density(subg))
	print "# Nodes: " + str(subg.number_of_nodes())
	print

	subg = greedy_triangle_density(G)
	print "----Greedy Triangle Density----"
	print "Degree Density: " + str(degree_density(subg))
	print "Density: " + str(density(subg))
	print "Triangle Density: " + str(triangle_density(subg))
	print "# Nodes: " + str(subg.number_of_nodes())
	print

	subg = greedy_triangle_graph_density(G)
	print "----Greedy Triangle-Graph Density----"
	print "Degree Density: " + str(degree_density(subg))
	print "Density: " + str(density(subg))
	print "Triangle Density: " + str(triangle_density(subg))
	print "# Nodes: " + str(subg.number_of_nodes())

if __name__ == "__main__":
    main()
