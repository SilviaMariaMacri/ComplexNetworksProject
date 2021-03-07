#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 21:25:02 2021

@author: caterina
"""

import pandas as pd
import networkx as nx



def HumanGraph(NameHumanFile, threshold_score): 
	
	GH = nx.Graph() 
	
	sapiens = pd.read_csv(NameHumanFile, sep=" ", usecols=['protein1', 'protein2', 'combined_score']) 
	#select human PPI links depending on the score
	sapiens=sapiens[sapiens['combined_score']>threshold_score]
	for i in range(len(sapiens)):
		GH.add_edge(sapiens.iloc[i][0], sapiens.iloc[i][1], weight=sapiens.iloc[i][2])
	
	return GH
	
	


def VirusGraph(NameVirusFile_i,VirusNames_i):
	virus = pd.read_csv(NameVirusFile_i, sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	graph_virus = nx.DiGraph() 
	
	for j in range(len(virus)):
		graph_virus.add_edge(virus.iloc[j][0], virus.iloc[j][1], weight=virus.iloc[j][2])
	
	graph_virus.name=VirusNames_i
		
	print(nx.info(graph_virus),'\n')
	
	return graph_virus


def SubnetworkGraphs (GH,GV_i, VirusNames_i):

	Gsub_i = GH.subgraph(nx.nodes(GV_i))
	
	Gsub_i.name=VirusNames_i
	print(nx.info(Gsub_i),'\n')
		
		
	return Gsub_i

