#altro metodo per selezionare la subnetwork (l'aktro è sbagliato)

import Percolation 




	

hitnodes_tot = []

for i in range(len(NameVirusFile)):
	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	GV = Percolation.HumanGraph(virus) #se la volgiamo indiretta
	hitnodes_non_selected = Percolation.HitnodesNonSelectedString(GV)	
	hitnodes = Percolation.Hitnodes(GH,hitnodes_non_selected,NameVirusFile)
	hitnodes_tot = np.hstack((hitnodes, hitnodes_tot))
	
GH_sub = GH.subgraph(hitnodes_tot)
	

'''
nx.info(GH_sub)
Out[10]: 'Name: \nType: Graph\nNumber of nodes: 2239\nNumber of edges: 
	57183\nAverage degree:  51.0791'''















#%%

import os 
import pandas as pd
import networkx as nx



directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/ComplexNetworksProject/ComplexNetworksProject'
#directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/data'
os.chdir(directory) 


import Percolation
import entropy_canonical_solofunzione



NameVirusFile = ['string_interactions_WNV.tsv','string_interactions_varicella.tsv','string_interactions_SARSCov.tsv','string_interactions_parechovirus2.tsv','string_interactions_mumps.tsv','string_interactions_MARV.tsv','string_interactions_lassa.tsv','string_interactions_InfluenzaA.tsv','string_interactions_HTLV-1.tsv','string_interactions_HPV1a.tsv','string_interactions_HIV1_553.tsv','string_interactions_hepatitisB.tsv','string_interactions_ebola.tsv','string_interactions_dengue2.tsv','string_interactions_cytomegalo.tsv','ProteinNamesString_AS_OTHER_VIRUSES_TabSep.txt']




graphs = []



for i in range(len(NameVirusFile)):	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	graph_virus = Percolation.HumanGraph(virus)
	print(nx.info(graph_virus))


	#tolgo nodi corrispondenti a proteine virus
	for j in range(len(virus)):
		if virus.iloc[j,0].startswith('9606.')!=True:
			try:
				graph_virus.remove_node(virus.iloc[j,0])
			except  nx.NetworkXError:
				continue
			

	graphs.append(graph_virus)	

	print(nx.info(graph_virus))
	
	
			

#%%	
GH_subnetwork = nx.compose(graphs[0],graphs[1])

for i in range(1,len(graphs)-1):
	GH_subnetwork = nx.compose(GH_subnetwork,graphs[i+1])
	
	
'''
nx.info(GH_subnetwork)
Out[4]: 'Name: \nType: Graph\nNumber of nodes: 3178\nNumber of edges: 
	80792\nAverage degree:  50.8446'''
	
#%%	
	
	
GV = []
for i in range(len(NameVirusFile)):
	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	GV_single = Percolation.VirusGraph(virus)
	#GV.append(GV_single)
	GV_undirected = GV_single.to_undirected()	
	GV.append(GV_undirected)
	

#GH_subnetwork_directed = GH_subnetwork.to_directed()



#%%

#FileNameAdjacency = ['adjacency_influenzaA_undir.txt','adjacency_HIV1_553_undir.txt','adjacency_SARSCov_undir.txt','adjacency_HTLV-1_undir.txt','adjacency_WNV_undir.txt','adjacency_covid19_undir.txt']


#network = []
#Sc = []
Ss = []
for i in range(len(GV)):
	#Gnew = nx.compose(GV[i],GH_subnetwork_directed)
	Gnew = nx.compose(GV[i],GH_subnetwork)
	#network.append(Gnew)
	pp = nx.to_pandas_adjacency(Gnew,weight='weight')
	#pp.to_csv(FileNameAdjacency[i], sep="\t") 

	Ss_single = entropy_canonical_solofunzione.entropy_canonical_s(pp)
	#Sc_single = entropy_canonical_solofunzione.entropy_canonical_c(pp)
	Ss.append(Ss_single)
	#Sc.append(Sc_single)
	
	







#%%  calcolo delle network undirected NON modificate da string
#con possibilità di togliere proteine virus


GV = []

#Sc=[]
Ss = []
for i in range(len(NameVirusFile)):
	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	GV_single = Percolation.HumanGraph(virus)
	print(nx.info(GV_single))
	
	#tolgo nodi corrispondenti a proteine virus (DA COMMENTARE DI VOLTA IN VOLTA)
#	for i in range(len(virus)):
#		if virus.iloc[i,0].startswith('9606.')!=True:
#			try:
#				GV_single.remove_node(virus.iloc[i,0])
#			except  nx.NetworkXError:
#				continue
			

	
	
	GV.append(GV_single)
	

	pp = nx.to_pandas_adjacency(GV_single,weight='weight')
	Ss_single = entropy_canonical_solofunzione.entropy_canonical_s(pp)
	#Sc_single = entropy_canonical_solofunzione.entropy_canonical_c(pp)
	#Sc_single = entropy_canonical_solofunzione.entropy_canonical_c(pp)
	Ss.append(Ss_single)
	#Sc.append(Sc_single)

#%%

#salvo entropie in array (l'ordine è quello di NameVirusFile)
entropy_0 = []
entropy_r = []
for i in range(len(Ss)):
	entropy_0_single = Ss[i][0]
	entropy_0.append(entropy_0_single)
	entropy_r_single = Ss[i][1]
	entropy_r.append(entropy_r_single)
	


#%% calcolo differenze fra entropia finale e iniziale


differenze_0 = []
differenze_r = []
for i in range(len(entropy_0)):
	#diff_0 = entropy_0[i] - entropy_0_novirus[i]
	#differenze_0.append(diff_0)
	#diff_r = entropy_r[i] - entropy_r_novirus[i]
	#differenze_r.append(diff_r)
	#normalizzate
	diff_0 = (entropy_0[i] - entropy_0_novirus[i])/entropy_0_novirus[i]
	differenze_0.append(diff_0)
	diff_r = (entropy_r[i] - entropy_r_novirus[i])/entropy_r_novirus[i]
	differenze_r.append(diff_r)	
	
#%% grafico
NomiVirus = ['WNV','Varicella','SARSCov','Parechovirus 2','Mumps','MARV','Lassa','Influenza A','HTLV-1','HPV1a','HIV 1','Hepatitis B','Ebola','Dengue 2','Cytomegalo']


import matplotlib.pyplot as plt


#nomi = NameVirusFile.copy()
#nomi = ['string_interactions_WNV.tsv','string_interactions_varicella.tsv','string_interactions_SARSCov.tsv','string_interactions_parechovirus2.tsv','string_interactions_mumps.tsv','string_interactions_MARV.tsv','string_interactions_lassa.tsv','string_interactions_InfluenzaA.tsv','string_interactions_HTLV-1.tsv','string_interactions_HPV1a.tsv','string_interactions_HIV1_553.tsv','string_interactions_hepatitisB.tsv','string_interactions_ebola.tsv','string_interactions_dengue2.tsv','string_interactions_cytomegalo.tsv']


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,5))
ax.plot(NomiVirus, differenze_0 , marker="o", label='$\Delta S_0$', linewidth=0.75)
ax.plot(NomiVirus, differenze_r , marker="o", label='$\Delta S_r$', linewidth=0.75)
plt.grid(True)
plt.xticks(rotation=90)
ax.set_xlabel('virus')
ax.set_ylabel('$\Delta S$')
ax.legend(ncol=1 ,loc='best', fontsize=14)
#	plt.savefig(nameplot)
#	plt.show()


#np.arange(1,17,1)
	
	
#%%	
'''
CON PROTEINE VIRUS:

silvia
entropy_0
Out[16]: 
[20.34831613164889,
 35.558078694864584,
 11.919341432886155,
 30.630890406300598,
 26.971318683418314,
 16.44567160193866,
 10.259825019635064,
 48.84279151102817,
 19.438780986633336,
 15.642679017265841,
 62.056385822523474,
 33.994662155916956,
 21.252919304532337,
 14.67430745760102,
 37.50077142430058,
 nan]

silvia
entropy_r
Out[17]: 
[28.52441088639447,
 44.494987797779835,
 16.008728584117023,
 37.880562420998295,
 33.40127348845125,
 22.188726724590452,
 14.700272964943661,
 67.08643545357711,
 30.912827695945637,
 21.954090317813527,
 77.53870029521114,
 49.526718617853945,
 27.071810574844843,
 20.952940431800712,
 54.47904052600051,
 nan]






SENZA PROTEINE VIRUS

silvia
entropy_0_novirus
Out[26]: 
[19.04867001201345,
 34.36664601105497,
 10.484222763188487,
 29.658312589017807,
 25.873203525289185,
 15.153303679391614,
 9.393486358593599,
 47.087459877316036,
 18.22005952234895,
 13.949192224790345,
 60.46765316485114,
 32.87620611537482,
 19.995697889402386,
 13.450340373443154,
 37.392094217755236,
 nan]

silvia
entropy_r_novirus
Out[27]: 
[27.628490027544828,
 43.36969550201079,
 14.640003816099885,
 37.183327447992895,
 31.91307161024403,
 21.217508984242773,
 13.758252509335767,
 65.91823576307482,
 29.54462816799193,
 20.36684052541546,
 76.38897872047309,
 48.79608187588362,
 26.002982222000487,
 19.99335062844145,
 54.75521566390266,
 nan]
'''

#%%

'''
network indirette virus (come da file string)


CON PROTEINE VIRUS:
	
['string_interactions_WNV.tsv',
 'string_interactions_varicella.tsv',
 'string_interactions_SARSCov.tsv',
 'string_interactions_parechovirus2.tsv',
 'string_interactions_mumps.tsv',
 'string_interactions_MARV.tsv',
 'string_interactions_lassa.tsv',
 'string_interactions_InfluenzaA.tsv',
 'string_interactions_HTLV-1.tsv',
 'string_interactions_HPV1a.tsv',
 'string_interactions_HIV1_553.tsv',
 'string_interactions_hepatitisB.tsv',
 'string_interactions_ebola.tsv',
 'string_interactions_dengue2.tsv',
 'string_interactions_cytomegalo.tsv',
 'ProteinNamesString_AS_OTHER_VIRUSES_TabSep.txt']


Name: 
Type: Graph
Number of nodes: 251
Number of edges: 3634
Average degree:  28.9562
Name: 
Type: Graph
Number of nodes: 548
Number of edges: 11718
Average degree:  42.7664
Name: 
Type: Graph
Number of nodes: 234
Number of edges: 1596
Average degree:  13.6410
Name: 
Type: Graph
Number of nodes: 208
Number of edges: 5981
Average degree:  57.5096
Name: 
Type: Graph
Number of nodes: 335
Number of edges: 6432
Average degree:  38.4000
Name: 
Type: Graph
Number of nodes: 164
Number of edges: 2073
Average degree:  25.2805
Name: 
Type: Graph
Number of nodes: 100
Number of edges: 896
Average degree:  17.9200
Name: 
Type: Graph
Number of nodes: 561
Number of edges: 16618
Average degree:  59.2442
Name: 
Type: Graph
Number of nodes: 361
Number of edges: 5260
Average degree:  29.1413
Name: 
Type: Graph
Number of nodes: 235
Number of edges: 2358
Average degree:  20.0681
Name: 
Type: Graph
Number of nodes: 623
Number of edges: 22717
Average degree:  72.9278
Name: 
Type: Graph
Number of nodes: 179
Number of edges: 6898
Average degree:  77.0726
Name: 
Type: Graph
Number of nodes: 246
Number of edges: 3066
Average degree:  24.9268
Name: 
Type: Graph
Number of nodes: 212
Number of edges: 1927
Average degree:  18.1792
Name: 
Type: Graph
Number of nodes: 578
Number of edges: 15052
Average degree:  52.0830
Name: 
Type: Graph
Number of nodes: 465
Number of edges: 439
Average degree:   1.8882







SENZA PROTEINE:
	
	Name: 
Type: Graph
Number of nodes: 241
Number of edges: 3311
Average degree:  27.4772
Name: 
Type: Graph
Number of nodes: 502
Number of edges: 10837
Average degree:  43.1753
Name: 
Type: Graph
Number of nodes: 216
Number of edges: 1318
Average degree:  12.2037
Name: 
Type: Graph
Number of nodes: 206
Number of edges: 5775
Average degree:  56.0680
Name: 
Type: Graph
Number of nodes: 327
Number of edges: 6043
Average degree:  36.9602
Name: 
Type: Graph
Number of nodes: 159
Number of edges: 1914
Average degree:  24.0755
Name: 
Type: Graph
Number of nodes: 99
Number of edges: 797
Average degree:  16.1010
Name: 
Type: Graph
Number of nodes: 550
Number of edges: 15872
Average degree:  57.7164
Name: 
Type: Graph
Number of nodes: 353
Number of edges: 4879
Average degree:  27.6431
Name: 
Type: Graph
Number of nodes: 230
Number of edges: 2002
Average degree:  17.4087
Name: 
Type: Graph
Number of nodes: 610
Number of edges: 21827
Average degree:  71.5639
Name: 
Type: Graph
Number of nodes: 174
Number of edges: 6648
Average degree:  76.4138
Name: 
Type: Graph
Number of nodes: 239
Number of edges: 2810
Average degree:  23.5146
Name: 
Type: Graph
Number of nodes: 207
Number of edges: 1720
Average degree:  16.6184
Name: 
Type: Graph
Number of nodes: 508
Number of edges: 13822
Average degree:  54.4173
Name: 
Type: Graph
Number of nodes: 438
Number of edges: 0
Average degree:   0.0000

























CONSIDERANDO NETWORK COMUNE


SENZA VIRUS:
entropy_0_novirus: 54.95060112493653
entropy_r_novirus: 77.69881236343981

CON VIRUS:

silvia 
entropy_0
Out[79]: 
[55.032675796579525,
 55.187430762117906,
 54.98951335275074,
 55.055601950795534,
 55.16496144262161,
 55.01706434488959,
 55.00396709600124,
 55.314943837850926,
 55.17149507774756,
 55.15970673926066,
 55.373443767385936,
 55.03140125606253,
 55.04795659409541,
 55.06762014463314,
 54.86446338293951,
 nan]

entropy_r
Out[80]: 
[77.72125856133682,
 77.82817273013227,
 77.67806836442605,
 77.79185157313218,
 77.89698372805499,
 77.74021177658898,
 77.74517534351597,
 78.02485995760487,
 77.86098457785381,
 77.8254277203943,
 78.1162807119109,
 77.7820224012025,
 77.76776028816181,
 77.76133242117604,
 77.53923212448844,
 440.11679912425575]

















