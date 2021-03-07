#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 18:26:36 2021

@author: caterina
"""
import os
import entropy_canonical_function
import pandas as pd
import networkx as nx
import entropy_canonical_solofunzione
#%%

NameVirusFile = ['string_interactions_WNV.tsv','string_interactions_varicella.tsv'				 
				 ,'string_interactions_SARSCov.tsv','string_interactions_parechovirus2.tsv',
				 'string_interactions_mumps.tsv','string_interactions_MARV.tsv',
				 'string_interactions_lassa.tsv','string_interactions_InfluenzaA.tsv',
				 'string_interactions_HTLV-1.tsv','string_interactions_HPV1a.tsv',
				 'string_interactions_HIV1_553.tsv','string_interactions_hepatitisB.tsv',
				 'string_interactions_ebola.tsv','string_interactions_dengue2.tsv',
				 'string_interactions_cytomegalo.tsv','Covid19.txt']

VirusNames=['WNV','Varicella zoster virus','SARS-CoV','Human parechovirus 2',
			'Mumps virus','MARV','Lassa virus','Influenza A','HTLV1','HPV type 1a',
			'HIV1','Hepatitis B','Ebola','Dengue type 2','Cytomegalo','SARS-CoV-2']

#%%
##%% virus graphs (directed)
#directory= '/home/caterina/Scrivania/CN'
#os.chdir(directory)
#G = []
#
#for i in range(len(NameVirusFile)):	
#	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
#	
#	graph_virus = nx.DiGraph() 
#	for j in range(len(virus)):
#		graph_virus.add_edge(virus.iloc[j][0], virus.iloc[j][1], weight=virus.iloc[j][2])
#	print(VirusNames[i])
#	print(nx.info(graph_virus))
#	G.append(graph_virus)
##%%
##CREATING THE SUBGRAPHS OF EACH VIRUS FROM GH (undirected)
#subgraphs=[] #array of subgraphs that we can modify
#sub=[]#array of the subgraphs
#for a in range (len(NameVirusFile)):
#	# create dataframe of virus-human interactions
#	virus = pd.read_csv(NameVirusFile[a], sep="\t",usecols=['node1_external_id','node2_external_id', 'combined_score'])
#
#	#Selects only human protein hit by the virus
#	hit=set() #set does not allow repeated elements
#	for i in range(len(virus)):
#		if virus.iloc[i,0].startswith('9606.')!=True & virus.iloc[i,1].startswith('9606.')==True:	
#			hit.add(virus.iloc[i][1])
#
#	#creates a subgraph of human network with only proteins hit by the virus
#	
#	sub.append(GH.subgraph(hit))
#	subgraphs.append(nx.Graph(sub[a]))
#
#	
	
#%%
#ENTROPY ANALYSIS VIRUS NETWORKS	
directory= '/home/caterina/Scrivania/CN'
os.chdir(directory)
S0c_array=[]
Src_array=[]	
S0s_array=[]
Srs_array=[]
	


for i in range(len(NameVirusFile)):
	pp = nx.to_pandas_adjacency(GV[i],weight='weight')
	[S_0, S_r, z, P]=entropy_canonical_function.entropy_canonical_c(pp)
	S0c_array.append(S_0)
	Src_array.append(S_r)
	[S_0, S_r, z, P]=entropy_canonical_solofunzione.entropy_canonical_s(pp)
	S0s_array.append(S_0)
	Srs_array.append(S_r)
entropy_dict={'S0 cate':S0c_array, 'S0 silvia':S0s_array, 'Sr cate':Src_array, 'Sr silvia':Srs_array }
entropy_virus=pd.DataFrame(data=entropy_dict)
entropy_virus.to_csv('Entropy_Virus.tsv', sep="\t", index=False) 
#%%
#ENTROPY ANALYSIS VIRUS NETWORKS MADE UNDIRECTED	
directory= '/home/caterina/Scrivania/CN'
os.chdir(directory)
S0c_array=[]
Src_array=[]	
S0s_array=[]
Srs_array=[]
	


for i in range(len(NameVirusFile)):
	pp = nx.to_pandas_adjacency(GVmod[i],weight='weight')
	[S_0, S_r, z, P]=entropy_canonical_function.entropy_canonical_c(pp)
	S0c_array.append(S_0)
	Src_array.append(S_r)
	[S_0, S_r, z, P]=entropy_canonical_solofunzione.entropy_canonical_s(pp)
	S0s_array.append(S_0)
	Srs_array.append(S_r)
entropy_dict={'S0 cate':S0c_array, 'S0 silvia':S0s_array, 'Sr cate':Src_array, 'Sr silvia':Srs_array }
entropy_virus_und=pd.DataFrame(data=entropy_dict)
entropy_virus_und.to_csv('Entropy_Virus_undirected.tsv', sep="\t", index=False) 


#%%
#ENTROPY ANALYSIS SUBNETWORKS	
S0c_array=[]
Src_array=[]	
S0s_array=[]
Srs_array=[]
	
for i in range(len(NameVirusFile)):
	pp = (nx.to_pandas_adjacency(Gsub[i],weight='weight')/1000)
	[S_0, S_r, z, P]=entropy_canonical_function.entropy_canonical_c(pp)
	S0c_array.append(S_0)
	Src_array.append(S_r)
	[S_0, S_r, z, P]=entropy_canonical_solofunzione.entropy_canonical_s(pp)
	S0s_array.append(S_0)
	Srs_array.append(S_r)
entropy_dict={'S0 cate':S0c_array, 'S0 silvia':S0s_array, 'Sr cate':Src_array, 'Sr silvia':Srs_array }
entropy_sub=pd.DataFrame(data=entropy_dict)
entropy_sub.to_csv('Entropy_Sub.tsv', sep="\t", index=False) 