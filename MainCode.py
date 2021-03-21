import os 

directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/ComplexNetworksProject/ComplexNetworksProject'
#directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/data'
os.chdir(directory) 












#%%  create human and human-virus PPI graphs 


# STRING file of human PPI network
NameHumanFile='9606.protein.links.v11.0.txt.gz'

# STRING files of human-virus PPI networks
NameVirusFile = ['string_interactions_cytomegalo.tsv',
				 'string_interactions_dengue2.tsv',
				 'string_interactions_ebola.tsv',
				 'string_interactions_hepatitisB.tsv',
				 'string_interactions_HIV1_553.tsv',
				 'string_interactions_HPV1a.tsv',
				 'string_interactions_HTLV-1.tsv',
				 'string_interactions_InfluenzaA.tsv',
				 'string_interactions_lassa.tsv',
				 'string_interactions_MARV.tsv',
				 'string_interactions_mumps.tsv',
				 'string_interactions_parechovirus2.tsv',
				 'string_interactions_SARSCov.tsv',
				 'string_interactions_varicella.tsv',
				 'string_interactions_WNV.tsv']#,'Covid19.txt']

# list of considered viruses
VirusNames=['Cytomegalo',
			'Dengue type 2',
			'Ebola',
			'Hepatitis B',
			'HIV1',
			'HPV type 1a',
			'HTLV1',
			'Influenza A',
			'Lassa virus',
			'MARV',
			'Mumps virus',
			'Human parechovirus 2',
			'SARS-CoV',
			'Varicella zoster virus',
			'WNV']#,'SARS-CoV-2']





import NetworksBuilder



#human PPI network
GH=NetworksBuilder.HumanGraph(NameHumanFile,250,'human PPI')

#%%

#list of human-virus PPI
GV=[] 
for i in range (len(NameVirusFile)):
	GV_i = NetworksBuilder.VirusGraph(NameVirusFile[i], VirusNames[i])
	GV.append(GV_i)

#list of human PPI subnetworks related to each virus
Gsub=[] 
for i in range (len(NameVirusFile)):
	Gsub_i = NetworksBuilder.SubGraph(GH,GV[i])
	Gsub.append(Gsub_i)
	




 


#%%  Network characterization



import NetworkCharacterization



# dataframe of K, Knn, BC, CLof human PPI network
#GH_centrality = NetworkCharacterization.NetworkCharacterization(GH)


GV_centrality = []
for i in range(len(GV)):
	# list of three dataframes of human-virus PPI network 
	# (complete graph, subgraph of human nodes and subgraph of virus nodes)
	# K, Kin, Kout, Knn_in, Knn_out, BC, CL
	GV_centrality_i = NetworkCharacterization.NetworkCharacterization(GV[i])
	GV_centrality.append(GV_centrality_i)
	print('GV ',VirusNames[i],' completed')
	
	
Gsub_centrality = []
for i in range(len(Gsub)):
	# dataframe of human PPI subgraph
	# K, Knn, BC, CL
	Gsub_centrality_i = NetworkCharacterization.NetworkCharacterization(Gsub[i])
	Gsub_centrality.append(Gsub_centrality_i)
	print('Gsub ',VirusNames[i],' completed')










# characterization plots of human PPI network

NetworkCharacterization.PlotDegreeHist(GH_centrality,'Human PPI network')
NetworkCharacterization.PlotDegreeNNvsDegree(GH_centrality,'Human PPI network')
NetworkCharacterization.PlotINvsOUT(GH_centrality,'Human PPI network')
n=1000
NetworkCharacterization.PlotBcClvsDegree(GH_centrality,'Human PPI network',n)



#%%

# characterization plots of human-virus PPI networks 

for i in range(len(GV)):
	# histogram of degree
	NetworkCharacterization.PlotDegreeHist(GV_centrality[i],VirusNames[i])
	# knn vs K
	NetworkCharacterization.PlotDegreeNNvsDegree(GV_centrality[i],VirusNames[i])
	# Kin vs Kout
	NetworkCharacterization.PlotINvsOUT(GV_centrality[i],VirusNames[i])
	# BC,Cl vs K
	n=21
	NetworkCharacterization.PlotBcClvsDegree(GV_centrality[i],VirusNames[i],n)

#%%


# characterization plots of human PPI subgraphs

for i in range(len(Gsub)):
	# histogram of degree
	NetworkCharacterization.PlotDegreeHist(Gsub_centrality[i],VirusNames[i])
	# Knn vs K
	NetworkCharacterization.PlotDegreeNNvsDegree(Gsub_centrality[i],VirusNames[i])
	# BC,CL vs K
	n=21
	NetworkCharacterization.PlotBcClvsDegree(Gsub_centrality[i],VirusNames[i],n)







#%%  Percolation





import Percolation


for i in range(len(GV)):
	Percolation.Percolation(GV[i], GV_centrality[i], VirusNames[i])


for i in range(len(Gsub)):
	Percolation.Percolation(Gsub[i], Gsub_centrality[i], VirusNames[i])
	
	

	
	
	
#%%  Entropy


import EntropyCanonical


E = EntropyCanonical.EntropyDifference(VirusNames,Gsub,GV)
	

