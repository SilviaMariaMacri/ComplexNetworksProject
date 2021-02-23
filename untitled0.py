#%% 
##############################################################################
# creo file .txt per covid 19 
##############################################################################





# creo dataframe da network PPI human completa
sapiens = pd.read_csv('9606.protein.links.v11.0.txt', sep=" ", 
                    usecols=['protein1', 'protein2', 'combined_score']) 
sapiens=sapiens[sapiens['combined_score']>250]

# creo network diretta PPI human
GH = nx.DiGraph() 
for i in range(len(sapiens)):
	GH.add_edge(sapiens.iloc[i][0], sapiens.iloc[i][1], weight=sapiens.iloc[i][2])
'''DiGraph\nNumber of nodes: 19322\nNumber of edges: 4505762\nAverage 
   in degree: 233.1934\nAverage out degree: 233.1934' '''	



# creo dataframe da file covid 19
covid = pd.read_csv('ProteinNamesString_AS_OTHER_VIRUSES_TabSep.txt', sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
#len(covid)
#439




	
# creo subnetwork di GH con nodi human di covid
G_covid = GH.subgraph(covid['node2_external_id']).copy()
# dovrebbe avere numero di nodi = colonne del dataframe covid-1 (?)
'''nx.info(G_covid)
   'Name: \nType: DiGraph\nNumber of nodes: 438\nNumber of edges: 5992
   \nAverage in degree:  13.6804\nAverage out degree:  13.6804'''


# aggiungo link da dataframe covid a network G_covid 
for i in range(len(covid)):
	G_covid.add_edge(covid.iloc[i][0], covid.iloc[i][1], weight=covid.iloc[i][2])
'''nx.info(G_covid)
   'Name: \nType: DiGraph\nNumber of nodes: 465\nNumber of edges: 6431
   \nAverage in degree:  13.8301\nAverage out degree:  13.8301'	'''


# creo matrix dei link COVID 19
matrix = nx.to_pandas_edgelist(G_covid,source='node1_external_id', target='node2_external_id')
matrix['weight'] = matrix['weight']*0.001 #modifico il peso dei link
matrix = matrix.rename(columns={'weight': 'combined_score'}) #rinomino ultima colonna
	

# salvo la matrice in file .txt
matrix.to_csv('Covid19.txt', sep="\t", index=False) 


#import file as dataframe SOLO PER VERIFICA
#c = pd.read_csv('Covid19.txt', sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	



                                           ### FINE ###
								   
									     ### DA CONTROLLARE ###

















#%% importo tutte network interazione virus-uomo

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

NameVirusFile = ['string_interactions_WNV.tsv','string_interactions_varicella.tsv','string_interactions_SARSCov.tsv','string_interactions_parechovirus2.tsv','string_interactions_mumps.tsv','string_interactions_MARV.tsv','string_interactions_lassa.tsv','string_interactions_InfluenzaA.tsv','string_interactions_HTLV-1.tsv','string_interactions_HPV1a.tsv','string_interactions_HIV1_553.tsv','string_interactions_hepatitisB.tsv','string_interactions_ebola.tsv','string_interactions_dengue2.tsv','string_interactions_cytomegalo.tsv','Covid19.txt']


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
