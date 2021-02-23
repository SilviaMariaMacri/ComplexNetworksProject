import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

NameVirusFile = ['string_interactions_WNV.tsv','string_interactions_varicella.tsv','string_interactions_SARSCov.tsv','string_interactions_parechovirus2.tsv','string_interactions_mumps.tsv','string_interactions_MARV.tsv','string_interactions_lassa.tsv','string_interactions_InfluenzaA.tsv','string_interactions_HTLV-1.tsv','string_interactions_HPV1a.tsv','string_interactions_HIV1_553.tsv','string_interactions_hepatitisB.tsv','string_interactions_ebola.tsv','string_interactions_dengue2.tsv','string_interactions_cytomegalo.tsv']#,'ProteinNamesString_AS_OTHER_VIRUSES_TabSep.txt']


G = []
VirusProtein = [] 


for i in range(len(NameVirusFile)):	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	
	graph_virus = nx.DiGraph() 
	for j in range(len(virus)):
		graph_virus.add_edge(virus.iloc[j][0], virus.iloc[j][1], weight=virus.iloc[j][2])
	print(nx.info(graph_virus))
	G.append(graph_virus)
	
	
	
#	VirusProteinSingle = []	
   
#	for l in range(len(virus)):
#		for n in range(2):
#			if virus.iloc[l,n].startswith('9606.')==False:
#				VirusProteinSingle.append(virus.iloc[l][n])
#	VirusProtein.append(VirusProteinSingle)
	
	
	
#%%	
for i in range(len(G)):
	
	degreeIN = pd.DataFrame.from_dict(nx.in_degree_centrality(G[i]),orient='index')
	degreeOUT = pd.DataFrame.from_dict(nx.out_degree_centrality(G[i]),orient='index')
	#BC = pd.DataFrame.from_dict(nx.betweenness_centrality(G[i]),orient='index')
	#closeness = pd.DataFrame.from_dict(nx.closeness_centrality(G[i]),orient='index')
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9,7))
	ax.scatter(degreeIN.iloc[:][0],degreeOUT.iloc[:][0],marker='o', linewidths=0.00001)
	
	for j in range(len(degreeIN)):
		if degreeIN.index[j].startswith('9606.')!=True:
			ax.scatter(degreeIN.iloc[j][0],degreeOUT.iloc[j][0],marker='o', linewidths=0.00001, color='r')
	
			
			
#	degreeINvirus = pd.DataFrame.from_dict(nx.in_degree_centrality(G[i]),orient='index')
#	degreeOUTvirus = pd.DataFrame.from_dict(nx.out_degree_centrality(G[i],VirusProtein[i]),orient='index')
	
	
	
	
	#ax.scatter(degreeIN.iloc[:][0],BC.iloc[:][0],marker='o', linewidths=0.00001)
	#ax.scatter(degreeOUT.iloc[:][0],closeness.iloc[:][0],marker='o', linewidths=0.00001)
	
	
	
	ax.set_xlabel('Degree')
	ax.set_ylabel('Betweenness')
	#ax.legend(ncol=1 ,loc='best', fontsize=14)
