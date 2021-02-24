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











