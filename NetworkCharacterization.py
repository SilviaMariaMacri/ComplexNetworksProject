import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import os

directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject'
os.chdir(directory) 

NameVirusFile = ['string_interactions_WNV.tsv','string_interactions_varicella.tsv','string_interactions_SARSCov.tsv','string_interactions_parechovirus2.tsv','string_interactions_mumps.tsv','string_interactions_MARV.tsv','string_interactions_lassa.tsv','string_interactions_InfluenzaA.tsv','string_interactions_HTLV-1.tsv','string_interactions_HPV1a.tsv','string_interactions_HIV1_553.tsv','string_interactions_hepatitisB.tsv','string_interactions_ebola.tsv','string_interactions_dengue2.tsv','string_interactions_cytomegalo.tsv','Covid19.txt']
VirusNames=['WNV','Varicella','SARSCov','Parechovirus2','Mumps','MARV','Lassa','InfluenzaA','HTLV-1','HPV1a','HIV1_553','HepatitisB','Ebola','Dengue2','Cytomegalo','SARSCov2']

	
###FUNZIONE PER DARE I NOMI AI FILE DEI GRAFICI
#il nome del file diventa: Prefisso + nome virus + formato del file
def FileNames(Prefix, FileType): 
	newnames=[]
	for a in range (len(VirusNames)):
		newnames.append(f"{Prefix}{VirusNames[a]}{FileType}")
	return newnames
	
	



#%%

G = []

for i in range(len(NameVirusFile)):	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	
	graph_virus = nx.DiGraph() 
	for j in range(len(virus)):
		graph_virus.add_edge(virus.iloc[j][0], virus.iloc[j][1], weight=virus.iloc[j][2])
	print(nx.info(graph_virus))
	G.append(graph_virus)
	
	
#%%	
#DEGREE IN - DEGREE OUT

#directory dove si salveranno i grafici		
directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/grafici'
os.chdir(directory) 

##nome dei file
NamesDegree=FileNames('DegreeIn_DegreeOut_','.png')

for i in range(len(G)):	
	degreeIN = pd.DataFrame.from_dict(nx.in_degree_centrality(G[i]),orient='index')
	degreeOUT = pd.DataFrame.from_dict(nx.out_degree_centrality(G[i]),orient='index')
	#BC = pd.DataFrame.from_dict(nx.betweenness_centrality(G[i]),orient='index')
	#closeness = pd.DataFrame.from_dict(nx.closeness_centrality(G[i]),orient='index')
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9,7))
	ax.scatter(degreeIN.iloc[:][0],degreeOUT.iloc[:][0],marker='o', linewidths=0.00001, label='human')
	
	for j in range(len(degreeIN)):
		if degreeIN.index[j].startswith('9606.')!=True:
			ax.scatter(degreeIN.iloc[j][0],degreeOUT.iloc[j][0],marker='o', linewidths=0.00001, color='r', label='virus')

	
	
	ax.set_xlabel('Degree in')
	ax.set_ylabel('Degree out')
	
	ax.legend(['human','virus'])
	
	
	plt.savefig(NamesDegree[i])
	#plt.show()
	
#%%

#BETWEENNESS

#directory dove si salveranno i grafici		
directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/grafici'
os.chdir(directory) 

##nome dei file
NamesBC=FileNames('HistoBC_','.png')

for i in range(len(G)):		
	bc=pd.DataFrame.from_dict(nx.betweenness_centrality(G[i],weight='weight'),orient='index',columns=['BC'])
	bc.hist()
	
	
	plt.savefig(NamesBC[i])

#%%

#CLUSTERING COEFFICIENT

#directory dove si salveranno i grafici	
directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/grafici'
os.chdir(directory) 

##nome dei file
NamesCC=FileNames('HistoClusteringCoeff_','.png')
NamesCC_degreeIN=FileNames('ClustCoeff_DegreeIN_','.png')
NamesCC_degreeOUT=FileNames('ClustCoeff_DegreeOUT_','.png')



for i in range(len(G)):		
	#clustering coefficient
	cc=pd.DataFrame.from_dict(nx.clustering(G[i],weight='weight'),orient='index',columns=['CC'])
	
	degreeIN = pd.DataFrame.from_dict(nx.in_degree_centrality(G[i]),orient='index')
	degreeOUT = pd.DataFrame.from_dict(nx.out_degree_centrality(G[i]),orient='index')
	
	#cc.hist() #HISTOGRAM
	#plt.savefig(NamesCC[i])
	
	
	#CLUSTERING COEFFICIENT VS DEGREE IN
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9,7))
	ax.scatter(degreeIN.iloc[:][0],cc.loc[:,'CC'],marker='o', linewidths=0.00001, label='human')
	
	for j in range(len(degreeIN)):
		if degreeIN.index[j].startswith('9606.')!=True:
			ax.scatter(degreeIN.iloc[j][0],cc.iloc[j][0],marker='o', linewidths=0.00001, color='r', label='virus')

	ax.set_xlabel('Degree in')
	ax.set_ylabel('cc')
	
	ax.legend(['human','virus'])
		
	plt.savefig(NamesCC_degreeIN[i])
	
	
	
	
	
	#CLUSTERING COEFFICIENT VS DEGREE OUT
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9,7))
	ax.scatter(degreeOUT.iloc[:][0],cc.loc[:,'CC'],marker='o', linewidths=0.00001, label='human')
	
	for j in range(len(degreeOUT)):
		if degreeOUT.index[j].startswith('9606.')!=True:
			ax.scatter(degreeOUT.iloc[j][0],cc.iloc[j][0],marker='o', linewidths=0.00001, color='r', label='virus')

	ax.set_xlabel('Degree out')
	ax.set_ylabel('cc')
	
	ax.legend(['human','virus'])
		
	plt.savefig(NamesCC_degreeOUT[i])
	
