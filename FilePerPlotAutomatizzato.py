'''
NOTE:
	- errore: GH non era ricreata dopo BC_percolation CORRETTO
	- nodi di human PPI







HOMO SAPIENS = 9606.
Human immunodeficiency virus type 1 group M subtype B (isolate HXB2) = 11676.
Influenza A virus (strain A/Puerto Rico/8/1934 H1N1) = 11320.
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


#%%


#quanti file dati di virus abbiamo?
NumberOfViruses = 2 
# list of names of virus data
NameVirusFile = ['string_interactions_HIV_piccola.tsv','string_interactions_influenzaA_piccola.tsv']
# list of names of the file storing the size of the giant component
FileNamePercolation = ['PercolationHIV.txt','PercolationInfluenzaA.txt']
# list of names of final plots
PlotName = ['plotHIV.png','plotInfluenzaA.png']


for i in range(NumberOfViruses):
	IterativePercolation.PercolationCode(GH_reference,NameVirusFile[i],FileNamePercolation[i],PlotName[i])

