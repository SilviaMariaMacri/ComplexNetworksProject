#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 09:43:16 2020

@author: caterina

"""
######################## FUNZIONI DA CHIAMARE DA FUORI #######################
# PercolationSubgraphs(Gsub, NameVirusFile,VirusNames)
# PercolationVirus(G,NameVirusFile,VirusNames)


import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

'''
FUNZIONE CHE SVOLGE LA PERCOLATION SU UNA SINGOLA NETWORK
##############################################################################

INPUT:
	 G - rete su cui applicare percolation
	 BCdataframe - dataframe relativo a G con due colonne 'Nodi' e 'BC'
	 DegreeDataframe - dataframe relativo a G con due colonne 'Nodi' e 'Degree'
	 namefile - nome del file del grafico in output
	 virus_name - nome del virus relativo a G

OUTPUT: 
	 grafico percolation di G	

###############################################################################	 
'''
def PercolationGraph(G,BCdataframe,DegreeDataframe,namefile,virus_name):
	network=G.copy(as_view=False)
	nodes=np.array(network.nodes())
	
	#Percolation BC
	
	# Sorting BCselection dataframe according to BC value (from largest to lowest)
	BC_sorted = BCdataframe.sort_values(by='BC', ascending=False)
	# reindex BC_sorted 
	BC_sorted.index = np.arange(0,len(BC_sorted),1)
	
	sizeBC=[]
	
	if nx.is_directed(network)==True:
		largest_cc = max(nx.weakly_connected_components(network), key=len)
	else:
		largest_cc = max(nx.connected_components(network), key=len)
	
	sizeBC.append(len(largest_cc))
	removedBC=0
	
	i=0
	while len(list(network.nodes))>1:
		node_to_remove=BC_sorted.iloc[i][0]
		network.remove_node(node_to_remove)
		removedBC=removedBC+1
		
		if nx.is_directed(network)==True:
			largest_cc = max(nx.weakly_connected_components(network), key=len)
		else:
			largest_cc = max(nx.connected_components(network), key=len)
		
		sizeBC.append(len(largest_cc))
		i=i+1
		#print("BC: Rimossi",removedBC, "nodi su", len(BC_sorted), ", size=",len(largest_cc))
	
	
	#Percolation Degree
	network=G.copy(as_view=False)
		
	degree_sorted=DegreeDataframe.sort_values(by='Degree', ascending=False)
	degree_sorted.index=np.arange(0,len(degree_sorted),1)
	
	sizeDegree=[]
	removedDegree=0
	
	if nx.is_directed(network)==True:
		largest_cc = max(nx.weakly_connected_components(network), key=len)
	else:
		largest_cc = max(nx.connected_components(network), key=len)
	
	sizeDegree.append(len(largest_cc))
	i=0
	
	while len(list(network.nodes))>1:
		
		node_to_remove=degree_sorted.iloc[i][0]
		
		network.remove_node(node_to_remove)
		removedDegree=removedDegree+1
		
		if nx.is_directed(network)==True:
			largest_cc = max(nx.weakly_connected_components(network), key=len)
		else:
			largest_cc = max(nx.connected_components(network), key=len)
	
		sizeDegree.append(len(largest_cc))
#		
		i=i+1
		
		#print("Degree: Rimossi",removedDegree, "nodi","size=",len(largest_cc))
	
	#Percolation Random
	network=G.copy(as_view=False)
	removedRandom=[]
	removedRandom=0
	sizeRandom=[]
	
	if nx.is_directed(network)==True:
		largest_cc = max(nx.weakly_connected_components(network), key=len)
	else:
		largest_cc = max(nx.connected_components(network), key=len)
	
	sizeRandom.append(len(largest_cc))
	
	
	random_index = np.arange(len(nodes)) # create an array of indexes for the nodes
	np.random.shuffle(random_index)	
	
	i=0
	while len(list(network.nodes))>1:
		node_to_remove=nodes[random_index[i]]
		network.remove_node(node_to_remove)
		removedRandom=removedRandom+1
		if nx.is_directed(network)==True:
			largest_cc = max(nx.weakly_connected_components(network), key=len)
		else:
			largest_cc = max(nx.connected_components(network), key=len)
		sizeRandom.append(len(largest_cc))
		i=i+1
		#print("Random:Rimossi",removedRandom, "nodi su", len(nodes), ", size=",len(largest_cc))
	
	
	#	PLOT
	
	Removed_nodes=np.arange(1, max(len(sizeBC),len(sizeDegree),len(sizeRandom)+1), step=1)
	
	x_axis =[]
	for i in range (len(Removed_nodes)):
		x_axis.append(Removed_nodes[i]/max(Removed_nodes)*100)
	
	
#	else:
	labeldegree='Degree'
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.5,3.5))
	ax.plot(x_axis, sizeBC, label='BC', linewidth=0.75)
	ax.plot(x_axis, sizeDegree, label=labeldegree, linewidth=0.75)
	ax.plot(x_axis, sizeRandom, label='Random', linewidth=0.75)
	ax.set_title(virus_name)
	ax.set_xlabel('% of removed nodes',labelpad=2)
	ax.set_ylabel('Size of giant component',labelpad=2)
	ax.legend(ncol=1 ,loc='best', fontsize=12)
	
	plt.grid(True)
	plt.savefig(namefile)

#	plt.show()




'''
FUNZIONE CHE ITERA PercolationGraph SU TUTTI I SUBGRAPHS 
##############################################################################

INPUT:
	 Gsub - array di subnetwork della human PPI (undirected networks)
	 NameVirusFile - array con i nomi dei file del virus
	 VirusNames - array con i nomi dei virus completi
OUTPUT: 
	 Per ogni network dell'array chiama la funzione 'PercolationGraph'

###############################################################################	 
'''
def PercolationSubgraphs(Gsub, NameVirusFile,VirusNames):
		
	for j in range(len(NameVirusFile)):
		degree_dict=nx.degree(Gsub[j],weight='weight')
		degree = pd.DataFrame(data=degree_dict) 
		degree.columns=['Nodes', 'Degree']
		bc=pd.DataFrame.from_dict(nx.betweenness_centrality(Gsub[j],weight='weight'),orient='index',columns=['BC'])
		bc.insert(0, 'Nodi', bc.index, allow_duplicates = False)
		
		degree_array=[]
		nodes=np.array(Gsub[j].nodes())
		for i in range (len(nodes)):
			degree_array.append(Gsub[j].degree(nodes[i],weight='weight'))
	
		degree_dict = {'Nodi': nodes, 'Degree': degree_array}
		degree_df = pd.DataFrame(data=degree_dict)
		
		
		print(VirusNames[j])
		Name=f"{'Percol_SUBNETWORK_'}{VirusNames[j]}{'.png'}"
		PercolationGraph(Gsub[j],bc,degree_df,Name,VirusNames[j])






'''
FUNZIONE CHE ITERA PercolationGraph SU TUTTE LE VIRUS-HUMAN PPI
##############################################################################

INPUT:
	 G- array di virus_human interactome (directed networks)
	 NameVirusFile - array con i nomi dei file del virus
	 VirusNames - array con i nomi dei virus completi
OUTPUT: 
	 Per ogni network dell'array chiama la funzione 'PercolationGraph'

###############################################################################	 
'''
def PercolationVirus(G,NameVirusFile,VirusNames):
		
	for j in range (len(NameVirusFile)):
		
		bc=pd.DataFrame.from_dict(nx.betweenness_centrality(G[j],weight='weight'),orient='index',columns=['BC'])
		bc.insert(0, 'Nodi', bc.index, allow_duplicates = False)
		
		degree_array=[]
		virus_nodes=np.array(G[j].nodes())
		for i in range (len(virus_nodes)):
			degree_array.append(G[j].degree(virus_nodes[i],weight='weight'))
			
		degree_dict = {'Nodi': virus_nodes, 'Degree': degree_array}
		degree_df = pd.DataFrame(data=degree_dict)
		#The node degree is the number of edges adjacent to the node. 
		#The weighted node degree is the sum of the edge weights for edges incident to that node.
		
		print(VirusNames[j])
		Name=f"{'Percol_VIRUS_'}{VirusNames[j]}{'.png'}"
		PercolationGraph(G[j],bc,degree_df,Name,VirusNames[j])
		
