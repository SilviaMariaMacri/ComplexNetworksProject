#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 21:25:02 2021

@author: caterina
"""

import pandas as pd
import networkx as nx


def NetworksBuilder (NameHumanFile, NameVirusFile, VirusNames):
	
	GH = nx.Graph() 
	GV = [] #array of virus graphs
	Gsub = [] #array of subnetworks of human PPI (one for each virus)

	#HUMAN PPI
	score_threshold=250
	sapiens = pd.read_csv(NameHumanFile, sep=" ", usecols=['protein1', 'protein2', 'combined_score']) 
	#select human PPI links depending on the score
	sapiens=sapiens[sapiens['combined_score']>score_threshold]
	for i in range(len(sapiens)):
		GH.add_edge(sapiens.iloc[i][0], sapiens.iloc[i][1], weight=sapiens.iloc[i][2])
	
	
	#VIRUS NETWORK
	for i in range(len(NameVirusFile)):	
		virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
		graph_virus = nx.DiGraph() 
	for j in range(len(virus)):
		graph_virus.add_edge(virus.iloc[j][0], virus.iloc[j][1], weight=virus.iloc[j][2])
		print(VirusNames[i])
		print(nx.info(graph_virus))
		GV.append(graph_virus)


	#HUMAN SUBNETWORKS
	for i in range(len(GV)):
		Gsub_i = GH.subgraph(nx.nodes(GV[i]))
		Gsub.append(Gsub_i)


	return GH, GV, Gsub

