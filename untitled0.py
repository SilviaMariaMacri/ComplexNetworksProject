#%% 
##############################################################################
# creo file .txt per covid 19 
##############################################################################


import networkx as nx
import pandas as pd
import numpy as np


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

# create a copy of GH to store as a reference     
GH_reference = GH.copy(as_view=False)

#%%

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


















#%% #######################################################################


####### CREAZIONE DI RETE VIRUS, SUBNETWORK DI GH, UNIONE DELLE DUE



#%% virus       SEZIONE A

G = []

for i in range(len(NameVirusFile)):	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	
	graph_virus = nx.DiGraph() 
	for j in range(len(virus)):
		graph_virus.add_edge(virus.iloc[j][0], virus.iloc[j][1], weight=virus.iloc[j][2])
	print(nx.info(graph_virus))
	G.append(graph_virus)







#%% unisco subgraph ricavate da network PPI human grande con quelle dei virus


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

# create a copy of GH to store as a reference     
GH_reference = GH.copy(as_view=False)





#%%            SEZIONE B 

#human    parto da network dei virus e tolgo nodi corrispondenti a proteine non umane
# POSSO USARE QUESTO CODICE ANCHE DA SOLO SENZA UNIRLO CON QUELLA GRANDE

Gh = G.copy()

for i in range(len(Gh)):
	
	nodes = np.array(nx.nodes(G[i]))
	for j in range(len(nodes)):
		if nodes[j].startswith('9606.') != True:
			Gh[i].remove_node(nodes[j])
			


'''	
#%% human  NON SO PERCHè, è diverso dall'altro


G = []

for i in range(len(NameVirusFile)):	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	
	graph_virus = nx.DiGraph() 
	for j in range(len(virus)):
		if virus.iloc[j,0].startswith('9606.')==True:
			graph_virus.add_edge(virus.iloc[j][0], virus.iloc[j][1], weight=virus.iloc[j][2])
	print(nx.info(graph_virus))
	G.append(graph_virus)
	
	'''	
	
#%%    SEZIONE C


#Gh = rete ricavata da file virus togliendo proteine virus

G_unita = []

for i in range(len(Gh)):
	
	
	nodi = np.array(nx.nodes(Gh[i]))
	GH_ridotta = GH.subgraph(nodi).copy()
	#G_unita_singola = nx.compose(GH_ridotta,Gh[i])
	G_unita.append(G_unita_singola)




 
#%%

##################################### SEZIONE A
'''
Name: 
Type: DiGraph
Number of nodes: 251
Number of edges: 3634
Average in degree:  14.4781
Average out degree:  14.4781
Name: 
Type: DiGraph
Number of nodes: 548
Number of edges: 11718
Average in degree:  21.3832
Average out degree:  21.3832
Name: 
Type: DiGraph
Number of nodes: 234
Number of edges: 1596
Average in degree:   6.8205
Average out degree:   6.8205
Name: 
Type: DiGraph
Number of nodes: 208
Number of edges: 5981
Average in degree:  28.7548
Average out degree:  28.7548
Name: 
Type: DiGraph
Number of nodes: 335
Number of edges: 6432
Average in degree:  19.2000
Average out degree:  19.2000
Name: 
Type: DiGraph
Number of nodes: 164
Number of edges: 2073
Average in degree:  12.6402
Average out degree:  12.6402
Name: 
Type: DiGraph
Number of nodes: 100
Number of edges: 896
Average in degree:   8.9600
Average out degree:   8.9600
Name: 
Type: DiGraph
Number of nodes: 561
Number of edges: 16618
Average in degree:  29.6221
Average out degree:  29.6221
Name: 
Type: DiGraph
Number of nodes: 361
Number of edges: 5260
Average in degree:  14.5706
Average out degree:  14.5706
Name: 
Type: DiGraph
Number of nodes: 235
Number of edges: 2358
Average in degree:  10.0340
Average out degree:  10.0340
Name: 
Type: DiGraph
Number of nodes: 623
Number of edges: 22717
Average in degree:  36.4639
Average out degree:  36.4639
Name: 
Type: DiGraph
Number of nodes: 179
Number of edges: 6898
Average in degree:  38.5363
Average out degree:  38.5363
Name: 
Type: DiGraph
Number of nodes: 246
Number of edges: 3066
Average in degree:  12.4634
Average out degree:  12.4634
Name: 
Type: DiGraph
Number of nodes: 212
Number of edges: 1927
Average in degree:   9.0896
Average out degree:   9.0896
Name: 
Type: DiGraph
Number of nodes: 578
Number of edges: 15052
Average in degree:  26.0415
Average out degree:  26.0415
Name: 
Type: DiGraph
Number of nodes: 465
Number of edges: 6431
Average in degree:  13.8301
Average out degree:  13.8301
'''	


######################################### SEZIONE B


		
'''
Name: 
Type: DiGraph
Number of nodes: 241
Number of edges: 3311
Average in degree:  13.7386
Average out degree:  13.7386
Name: 
Type: DiGraph
Number of nodes: 500
Number of edges: 10837
Average in degree:  21.6740
Average out degree:  21.6740
Name: 
Type: DiGraph
Number of nodes: 215
Number of edges: 1318
Average in degree:   6.1302
Average out degree:   6.1302
Name: 
Type: DiGraph
Number of nodes: 206
Number of edges: 5775
Average in degree:  28.0340
Average out degree:  28.0340
Name: 
Type: DiGraph
Number of nodes: 327
Number of edges: 6043
Average in degree:  18.4801
Average out degree:  18.4801
Name: 
Type: DiGraph
Number of nodes: 159
Number of edges: 1914
Average in degree:  12.0377
Average out degree:  12.0377
Name: 
Type: DiGraph
Number of nodes: 99
Number of edges: 797
Average in degree:   8.0505
Average out degree:   8.0505
Name: 
Type: DiGraph
Number of nodes: 550
Number of edges: 15872
Average in degree:  28.8582
Average out degree:  28.8582
Name: 
Type: DiGraph
Number of nodes: 353
Number of edges: 4879
Average in degree:  13.8215
Average out degree:  13.8215
Name: 
Type: DiGraph
Number of nodes: 230
Number of edges: 2002
Average in degree:   8.7043
Average out degree:   8.7043
Name: 
Type: DiGraph
Number of nodes: 610
Number of edges: 21827
Average in degree:  35.7820
Average out degree:  35.7820
Name: 
Type: DiGraph
Number of nodes: 174
Number of edges: 6648
Average in degree:  38.2069
Average out degree:  38.2069
Name: 
Type: DiGraph
Number of nodes: 239
Number of edges: 2810
Average in degree:  11.7573
Average out degree:  11.7573
Name: 
Type: DiGraph
Number of nodes: 206
Number of edges: 1720
Average in degree:   8.3495
Average out degree:   8.3495
Name: 
Type: DiGraph
Number of nodes: 500
Number of edges: 13822
Average in degree:  27.6440
Average out degree:  27.6440
Name: 
Type: DiGraph
Number of nodes: 438
Number of edges: 5992
Average in degree:  13.6804
Average out degree:  13.6804
'''		

############################################## SEZIONE C




'''
Name: 
Type: DiGraph
Number of nodes: 241
Number of edges: 5401
Average in degree:  22.4108
Average out degree:  22.4108
Name: 
Type: DiGraph
Number of nodes: 500
Number of edges: 18703
Average in degree:  37.4060
Average out degree:  37.4060
Name: 
Type: DiGraph
Number of nodes: 215
Number of edges: 2709
Average in degree:  12.6000
Average out degree:  12.6000
Name: 
Type: DiGraph
Number of nodes: 206
Number of edges: 7480
Average in degree:  36.3107
Average out degree:  36.3107
Name: 
Type: DiGraph
Number of nodes: 327
Number of edges: 8018
Average in degree:  24.5199
Average out degree:  24.5199
Name: 
Type: DiGraph
Number of nodes: 159
Number of edges: 2607
Average in degree:  16.3962
Average out degree:  16.3962
Name: 
Type: DiGraph
Number of nodes: 99
Number of edges: 1103
Average in degree:  11.1414
Average out degree:  11.1414
Name: 
Type: DiGraph
Number of nodes: 550
Number of edges: 30328
Average in degree:  55.1418
Average out degree:  55.1418
Name: 
Type: DiGraph
Number of nodes: 353
Number of edges: 7308
Average in degree:  20.7025
Average out degree:  20.7025
Name: 
Type: DiGraph
Number of nodes: 230
Number of edges: 4226
Average in degree:  18.3739
Average out degree:  18.3739
Name: 
Type: DiGraph
Number of nodes: 610
Number of edges: 42731
Average in degree:  70.0508
Average out degree:  70.0508
Name: 
Type: DiGraph
Number of nodes: 174
Number of edges: 10378
Average in degree:  59.6437
Average out degree:  59.6437
Name: 
Type: DiGraph
Number of nodes: 239
Number of edges: 4148
Average in degree:  17.3556
Average out degree:  17.3556
Name: 
Type: DiGraph
Number of nodes: 206
Number of edges: 2995
Average in degree:  14.5388
Average out degree:  14.5388
Name: 
Type: DiGraph
Number of nodes: 500
Number of edges: 23482
Average in degree:  46.9640
Average out degree:  46.9640
Name: 
Type: DiGraph
Number of nodes: 438
Number of edges: 5992
Average in degree:  13.6804
Average out degree:  13.6804
	
