#%%

import os 
import pandas as pd
import networkx as nx


directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/ComplexNetworksProject/ComplexNetworksProject'
#directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/data'
os.chdir(directory) 


import Percolation
import entropy_canonical_solofunzione



NameVirusFile = ['string_interactions_InfluenzaA.tsv','string_interactions_HIV1_553.tsv','string_interactions_SARSCov.tsv','string_interactions_HTLV-1.tsv','string_interactions_WNV.tsv','ProteinNamesString_AS_OTHER_VIRUSES_TabSep.txt']


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
	#GV.append(GV_single)
	GV_undirected = GV_single.to_undirected()	
	GV.append(GV_undirected)
	

#GH_subnetwork_directed = GH_subnetwork.to_directed()


'''nx.info(GH_subnetwork_directed)
Out[36]: 'Name: \nType: DiGraph\nNumber of nodes: 1992\nNumber of edges: 
	85720\nAverage in degree:  43.0321\nAverage out degree:  43.0321'
'''
#%%

FileNameAdjacency = ['adjacency_influenzaA_undir.txt','adjacency_HIV1_553_undir.txt','adjacency_SARSCov_undir.txt','adjacency_HTLV-1_undir.txt','adjacency_WNV_undir.txt','adjacency_covid19_undir.txt']


network = []
#Sc = []
Ss = []
for i in range(len(GV)):
	#Gnew = nx.compose(GV[i],GH_subnetwork_directed)
	Gnew = nx.compose(GV[i],GH_subnetwork)
	network.append(Gnew)
	pp = nx.to_pandas_adjacency(Gnew,weight='weight')
	#pp.to_csv(FileNameAdjacency[i], sep="\t") 

	Ss_single = entropy_canonical_solofunzione.entropy_canonical_s(pp)
	#Sc_single = entropy_canonical_solofunzione.entropy_canonical_c(pp)
	Ss.append(Ss_single)
	#Sc.append(Sc_single)
	
	



#%% se vuoi calcolare singolarmente


Ss= entropy_canonical_solofunzione.entropy_canonical_s(pp)
Sc = entropy_canonical_solofunzione.entropy_canonical_c(pp)







Ss[0][0]
Out[41]: 89.85640815826397
matlab: 89.8247

Ss[0][1]
Out[42]: 66.91093566585765
matlab: nan

Ss[1][0]
Out[43]: 90.02233572842047
matlab: 89.7490

Ss[1][1]
Out[44]: 66.95935023349685
matlab: nan

Ss[2][0]
Out[45]: 89.36326604535934

Ss[2][1]
Out[46]: 66.638721897709

Ss[3][0]
Out[47]: 89.54034947518224

Ss[3][1]
Out[48]: 66.81875117315387

Ss[4][0]
Out[49]: 89.5062106914269

Ss[4][1]
Out[50]: 66.69799271636779

Ss[5][0]
Out[51]: nan

Ss[5][1]
Out[52]: 328.8249763194713


