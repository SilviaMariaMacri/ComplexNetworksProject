import os 
import pandas as pd


directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/ComplexNetworksProject/ComplexNetworksProject'
#directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/data'
os.chdir(directory) 





import Percolation




def PercolationCode(GH_reference,nomefilevirus,FileNamePercolation,nameplot):
	
	
	# recreate the original GH graph
	GH = GH_reference.copy(as_view=False)
	
	
	
	# create dataframe of virus-human interactions
	virus = pd.read_csv(nomefilevirus, sep="\t", 
                   usecols=['node1_external_id','node2_external_id', 'combined_score'])


	# create a directed graph of virus-human interactions 
	GV = Percolation.VirusGraph(virus)

	


	
	hitnodes = Percolation.HitnodesNonSelectedString(GV)

	#create dataframe of hit nodes and corresponding degrees
	ND = Percolation.NodeDegreeDF(GH,GV,hitnodes)

	#create array of hit nodes
	hitnodes = Percolation.Hitnodes(GH,GV,hitnodes)

	#import BC file as dataframe
	BC = pd.read_csv('BetweennessCentrality2.csv', sep=",", skiprows=1)
 
	#create dataframe of hit nodes and corresponding betweenness 
	BC_sorted = Percolation.NodeBetweennessDF(GH,GV,BC,hitnodes)






	

	# RANDOM PERCOLATION
	sizeG_random = Percolation.PercolationRandom(GH,ND)

	# recreate the original GH graph
	GH = GH_reference.copy(as_view=False)

	# DEGREE-BASED PERCOLATION
	sizeG_degree = Percolation.PercolationDegree(GH,ND,hitnodes)

	# recreate the original GH graph
	GH = GH_reference.copy(as_view=False)

	#BETWEENNESS CENTRALITY - BASED PERCOLATION
	sizeG_BC = Percolation.PercolationBetweenness(GH,BC_sorted)

	






	# create a dataframe with hit nodes and corresponding betweenness centrality
	columns_graph = {'sizeG_random': sizeG_random, 'sizeG_degree': sizeG_degree, 'sizeG_BC': sizeG_BC}
	graph_df = pd.DataFrame(data=columns_graph) 	


	#create file .txt with different values of the size of the giant component obtained by percolation
	graph_df.to_csv(FileNamePercolation, index=True, index_label = 'removed_nodes') 


	#import file as dataframe
	percol = pd.read_csv(FileNamePercolation, sep=",")






	#plot
	Percolation.GraphPercolation(percol,nameplot)
	
	
	
	
	return 





