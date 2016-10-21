# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 2016

@author: g.nikolentzos

"""

import string
import re
import networkx as nx
from nltk.stem.porter import *
from greedy_algorithms import greedy_triangle_graph_density


def load_file(filename):
	""" 
	Read the tab delimited file containing the labels and the docs.
	
	"""
	paper = ""

	with open(filename) as f:
		for line in f:
			paper += unicode(line, errors='ignore') + " "
    
	return paper


def preprocessing(paper):
	""" 
	Permorm data preprocessing.
	
	"""	   
    # Load list of stopwords
	stopwords = set(['i', 'me', 'my', 'myself', 'we', 'our', 'ours', 'ourselves', 'you', 'your', 'yours',
                 'yourself', 'yourselves', 'he', 'him', 'his', 'himself', 'she', 'her', 'hers',
                 'herself', 'it', 'its', 'itself', 'three', 'they', 'them', 'their', 'theirs', 'themselves',
                 'what', 'which', 'who', 'whom', 'this', 'that', 'these', 'those', 'am', 'is', 'are',
                 'was', 'were', 'be', 'been', 'being', 'have', 'has', 'had', 'having', 'do', 'does',
                 'did', 'doing', 'a', 'an', 'the', 'and', 'but', 'if', 'or', 'because', 'as', 'until',
                 'while', 'of', 'at', 'also', 'by', 'for', 'with', 'about', 'against', 'between', 'into',
                 'through', 'during', 'before', 'after', 'above', 'below', 'to', 'from', 'up', 'down',
                 'in', 'out', 'on', 'off', 'over', 'under', 'again', 'further', 'then', 'once', 'here',
                 'there', 'when', 'where', 'why', 'how', 'all', 'any', 'both', 'each', 'few', 'more',
                 'most', 'other', 'some', 'such', 'no', 'nor', 'not', 'one', 'only', 'own', 'same', 'so',
                 'than', 'too', 'very', 's', 't', 'two', 'can', 'will', 'just', 'don', 'should', 'now'])
  	
	punctuation = re.compile(r'([^A-Za-z])')
    #punctuation = set(string.punctuation)

	
	# Remove punctuation and lowercase
	preprocessed_paper = punctuation.sub(" ", paper)
	#current_doc = ''.join([w for w in doc.lower() if w not in punctuation])

	# Turn to lowercase
	preprocessed_paper = preprocessed_paper.lower()

	# Stopword removal
	preprocessed_paper = [w for w in preprocessed_paper.split() if w not in stopwords]  

	return preprocessed_paper

	
def create_graph_of_words(preprocessed_paper, window_size):
	""" 
	Create graphs of words

	"""
	G = nx.Graph()
	for i in range(len(preprocessed_paper)):
		if preprocessed_paper[i] not in G.nodes():
			G.add_node(preprocessed_paper[i])
		for j in range(i+1, i+window_size):
			if j < len(preprocessed_paper):
				G.add_edge(preprocessed_paper[i], preprocessed_paper[j])

	return G


def main():
	""" 
	Main method
	"""
	
	paper = load_file("paper_raw.txt")
	
	preprocessed_paper = preprocessing(paper)
	
	G = create_graph_of_words(preprocessed_paper, 3)
	
	subg = greedy_triangle_graph_density(G)
	print "----Keywords----"
	for node in subg.nodes():
		print node


if __name__ == "__main__":
    main()
