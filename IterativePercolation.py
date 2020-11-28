import os 
import pandas as pd

directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/ComplexNetworksProject/ComplexNetworksProject'
#directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/data'
os.chdir(directory) 





import Percolation




def IterativePercolation(BC,GH_reference,NameVirusFile,FileNamePercolation,nameplot):
	
	
	# recreate the original GH graph
	GH = GH_reference.copy(as_view=False)
	
	
	
	# create dataframe of virus-human interactions
	virus = pd.read_csv(NameVirusFile, sep="\t", 
                   usecols=['node1_external_id','node2_external_id', 'combined_score'])


	#covid
	if NameVirusFile == 'ProteinNamesString_AS_OTHER_VIRUSES_TabSep.txt':
		covid = virus['node2_external_id']
		covid = covid.drop(59)
		hitnodes_non_selected = list(covid)
	
	# other viruses
	
	else:
		# create a directed graph of virus-human interactions 
		GV = Percolation.VirusGraph(virus)
		
		hitnodes_non_selected = Percolation.HitnodesNonSelectedString(GV)

	print('len(hitnodes_non_selected): ',len(hitnodes_non_selected))





	#create dataframe of hit nodes and corresponding degrees
	ND = Percolation.NodeDegreeDF(GH,hitnodes_non_selected,NameVirusFile)

	#create array of hit nodes
	hitnodes = Percolation.Hitnodes(GH,hitnodes_non_selected,NameVirusFile)

	#create dataframe of hit nodes and corresponding betweenness 
	BC_sorted = Percolation.NodeBetweennessDF(GH,BC,hitnodes)







	Percolation.PercolationCode(GH_reference,ND,hitnodes,BC_sorted,FileNamePercolation,nameplot)




	return


