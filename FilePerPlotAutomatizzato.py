'''
NOTE:
	- errore: GH non era ricreata dopo BC_percolation CORRETTO
	- nodi di human PPI

- calcolare altre misure di centrality?






HOMO SAPIENS = 9606.
Human immunodeficiency virus type 1 group M subtype B (isolate HXB2) = 11676.
Influenza A virus (strain A/Puerto Rico/8/1934 H1N1) = 11320.


similarities at small scales:
	SARS-Cov
	Influenza A
	HAdV
similarities at larger scales:
	HIV1
	HTLV1

'''


import os 
import pandas as pd



directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/ComplexNetworksProject/ComplexNetworksProject'
#directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/data'
os.chdir(directory) 





import Percolation
import IterativePercolation




# create datadrame of human PPI links
sapiens = pd.read_csv('9606.protein.links.v11.0.txt', sep=" ", 
                    usecols=['protein1', 'protein2', 'combined_score']) 
#select human PPI links depending on the score
sapiens=sapiens[sapiens['combined_score']>400]





# create a directed graph of human PPI
GH = Percolation.HumanGraph(sapiens)

# create a copy of GH to store as a reference     
GH_reference = GH.copy(as_view=False)




#import BC file as dataframe
BC = pd.read_csv('BetweennessCentrality2.csv', sep=",", skiprows=1)
 
	
#%%




NumberOfViruses = 2 
# list of names of virus data
#NameVirusFile = 'ProteinNamesString_AS_OTHER_VIRUSES_TabSep.txt'
NameVirusFile = ['string_interactions_HIV_piccola.tsv','string_interactions_influenzaA_piccola.tsv']
# list of names of the file storing the size of the giant component
#FileNamePercolation = 'PercolationCovid.txt'
FileNamePercolation = ['PercolationHIV.txt','PercolationInfluenzaA.txt']
# list of names of final plots
#PlotName = 'PlotCovid.png'
PlotName = ['plotHIV.png','plotInfluenzaA.png']


#'se si usano più file:'
for i in range(NumberOfViruses):
	IterativePercolation.IterativePercolation(BC,GH_reference,NameVirusFile[i],FileNamePercolation[i],PlotName[i])
#'se si usa un solo file:'
#IterativePercolation.IterativePercolation(BC,GH_reference,NameVirusFile,FileNamePercolation,PlotName)









#%%

'''
UNICA PROTEINA UMANA RIPETUTA DUE VOLTE IN 'Preys'

In: covid[covid['Preys']=='Q9Y5L0']
Out: 
                Bait   Preys      MIST
59       SARS-CoV2 M  Q9Y5L0  0.600814
330  SARS-CoV2 orf7a  Q9Y5L0  0.616847
'''