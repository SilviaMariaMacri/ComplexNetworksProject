#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 23:20:32 2020

@author: caterina
"""

import os 
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
#import numpy as np
import seaborn as sns

#directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/ComplexNetworksProject/ComplexNetworksProject'
directory= '/home/caterina/Scrivania/Complex_Networks_Project'
os.chdir(directory) 
import Percolation
#%%
# create datadrame of human PPI links
sapiens = pd.read_csv('9606.protein.links.v11.0.txt', sep=" ", 
                    usecols=['protein1', 'protein2', 'combined_score']) 
#select human PPI links depending on the score
sapiens=sapiens[sapiens['combined_score']>400]
#%%
# create a directed graph of human PPI
GH = Percolation.HumanGraph(sapiens)
#%%
# create a copy of GH to store as a reference     
GH_reference = GH.copy(as_view=False)
#%%
NameVirusFile = 'string_interactions_InfluenzaA.tsv'
FileNamePercolation = 'PercolationInfluenzaA.txt'
PlotName = 'plotInfluenzaA.png'
#%%
# create dataframe of virus-human interactions
virus = pd.read_csv(NameVirusFile, sep="\t",usecols=['node1_external_id','node2_external_id', 'combined_score'])

#%%
# create a directed graph of virus-human interactions 
GV = Percolation.VirusGraph(virus)
#%%
#import BC file as dataframe
BC = pd.read_csv('BetweennessCentrality2.csv', sep=",", skiprows=1)
#%%
# create an array of human proteins hit by the virus
hitnodes_non_selected = Percolation.HitnodesNonSelectedString(GV)


#create dataframe of hit nodes and corresponding degrees
ND = Percolation.NodeDegreeDF(GH,hitnodes_non_selected,NameVirusFile)

#create array of hit nodes
hitnodes = Percolation.Hitnodes(GH,hitnodes_non_selected,NameVirusFile)

#create dataframe of hit nodes and corresponding betweenness 
BC_sorted = Percolation.NodeBetweennessDF(GH,BC,hitnodes)
#%%
#cycle that calcolates all the first neighbours 
Neighbors=pd.DataFrame(BC_sorted['nodes'])

#%%
first_neighbors_node=[]#list of 1st neighbors for each BC_sorted node
first_neighbors_selected=[]
first_neighbors_list=[]#list of list 1st neighbors for all BC_sorted node
for i in range(len(BC_sorted)):
	
	first_neighbors_node=[]#list of 1st neighbors for each BC_sorted node
	first_neighbors_node=list(nx.neighbors(GH,BC_sorted.iloc[i][0]))
	
	first_neighbors_selected=[]
	#cycle that allows to select only some neighbors according to score
	count=0
	for a in range (len(first_neighbors_node)):
		
		
		node=first_neighbors_node[a]
#		if  GH.edges[(BC_sorted.iloc[i][0],node)]['weight'] > 900:
#			first_neighbors_selected.append(node)
#			count=count+1
	print ("For hit node number:  ",i," selected ", count, " nodes" )

	first_neighbors_list.append(first_neighbors_selected)

Neighbors['first neighbors']=first_neighbors_list
	#Neighbors.iloc[i][1]=first_neighbors_node
#%%
#cycle that removes all the hit nodes
RemovedNodes=[]
size_BC=[]
for i in range(len(BC_sorted)):
	node_to_remove = BC_sorted.iloc[i][0]
	GH.remove_node(node_to_remove)
	RemovedNodes.append(node_to_remove)
	s=len(max(nx.strongly_connected_components(GH), key=len))
	size_BC.append(s)
	print("Rimosso nodo direttamente colpito numero: ",i, " /", len(BC_sorted)-1, "size: ",s)
	print("###################################################")
#%%
#cycle that removes all the neighbors
for i in range (len(Neighbors)):
	
	for a in range (len(Neighbors.iloc[i][1])):
		node_to_remove=Neighbors.iloc[i][1][a]
		try:
			GH.remove_node(node_to_remove)
		except  nx.NetworkXError:
			#pass
			continue
		
		RemovedNodes.append(node_to_remove)
		s=len(max(nx.strongly_connected_components(GH), key=len))
		size_BC.append(s)
		print("Rimossi ",a,"/",len(Neighbors.iloc[i][1])-1,"primi vicini del nodo numero: ",i, "/", len(BC_sorted)-1, " size:", s)
	print("###################################################################")
#%%
# create a dataframe with hit nodes and corresponding betweenness centrality
	columns_graph = {'sizeG_BC': size_BC}
	graph_df = pd.DataFrame(data=columns_graph) 	

	#create file .txt with different values of the size of the giant component obtained by percolation
	graph_df.to_csv(FileNamePercolation, sep="\t", index=True, index_label = 'removed_nodes') 

	#import file as dataframe
	percol = pd.read_csv(FileNamePercolation, sep="\t")


	sns.set_style('whitegrid')
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,5))
#	x_axis = percol['removed_nodes']/percol['removed_nodes'].iloc[len(percol)-1]
	x_axis = percol['removed_nodes']
	ax.plot(x_axis, percol['sizeG_BC'], label='betweenness', linewidth=0.75)
	ax.set_xlabel('% of removed nodes')
	ax.set_ylabel('Size of giant component')
	ax.legend(ncol=1 ,loc='best', fontsize=14)
	plt.savefig(PlotName)
	plt.show()
	
	
