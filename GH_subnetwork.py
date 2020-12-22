#%%

import os 
import pandas as pd
import networkx as nx


directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/ComplexNetworksProject/ComplexNetworksProject'
#directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/data'
os.chdir(directory) 


import Percolation
import entropy_canonical_solofunzione




NameVirusFile = ['string_interactions_InfluenzaA.tsv','string_interactions_HIV1_553.tsv','string_interactions_SARSCov.tsv','string_interactions_yellowfever.tsv','ProteinNamesString_AS_OTHER_VIRUSES_TabSep.txt']


graphs = []
for i in range(len(NameVirusFile)):
	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	graph_virus = Percolation.HumanGraph(virus)
	
	
	#tolgo nodi corrispondenti a proteine virus
	for i in range(len(virus)):
		if virus.iloc[i,0].startswith('9606.')!=True:
			try:
				graph_virus.remove_node(virus.iloc[i,0])
			except  nx.NetworkXError:
				continue
			

	graphs.append(graph_virus)	

	print(nx.info(graph_virus))
			

#%%	
GH_subnetwork = nx.compose(graphs[0],graphs[1])

for i in range(1,len(graphs)-1):
	GH_subnetwork = nx.compose(GH_subnetwork,graphs[i+1])
	
	
	
	
GV = []
for i in range(len(NameVirusFile)):
	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	GV_single = Percolation.VirusGraph(virus)
	GV.append(GV_single)
	
	

GH_subnetwork_directed = GH_subnetwork.to_directed()

#%%


Sc = []
Ss = []
for i in range(len(GV)):
	Gnew = nx.compose(GV[i],GH_subnetwork_directed)
	pp = nx.to_pandas_adjacency(Gnew,weight='weight')
	Ss_single = entropy_canonical_solofunzione.entropy_canonical_s(pp)
	Sc_single = entropy_canonical_solofunzione.entropy_canonical_c(pp)
	Ss.append(Ss_single)
	Sc.append(Sc_single)
	
	


graph_df.to_csv(FileNamePercolation, sep="\t", index=True, index_label = 'removed_nodes') 



#%% se vuoi calcolare singolarmente


Ss= entropy_canonical_solofunzione.entropy_canonical_s(pp)
Sc = entropy_canonical_solofunzione.entropy_canonical_c(pp)





