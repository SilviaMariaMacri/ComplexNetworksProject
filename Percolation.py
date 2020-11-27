import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns







# create a directed graph of human PPI
def HumanGraph(sapiens):
	
	GH = nx.DiGraph() 
	for i in range(len(sapiens)):
		GH.add_edge(sapiens.iloc[i][0], sapiens.iloc[i][1], weight=sapiens.iloc[i][2])
		
	return GH






# NOT FOR COVID DATA
# create a directed graph pf virus-human interactions
def VirusGraph(virus):
	
	# select only links with virus protein as starting node and human protein as arrival node 
	interactionsVH = []    
	for i in range(len(virus)):
		#if virus.iloc[i,0].startswith(string)==True & virus.iloc[i,1].startswith('9606.')==True:
		if virus.iloc[i,0].startswith('9606.')!=True & virus.iloc[i,1].startswith('9606.')==True:	
			interactionsVH.append(virus.iloc[i])
        

	#create a directed graph of virus-human interactions        
	GV = nx.DiGraph() 
	for i in range(len(interactionsVH)): 
		GV.add_edge(interactionsVH[i][0], interactionsVH[i][1], weight=interactionsVH[i][2])
    
	print(nx.info(GV))	
	   
	return GV





# NOT FOR COVID DATA
def HitnodesNonSelectedString(GV):
	
	# create an array of human proteins hit by the virus
	hitnodes=[]
	for i in nx.nodes(GV):
		if i.startswith('9606.')==True:
			hitnodes.append(i)
			
	print('len(HitnodesNonSelectedString): ', len(hitnodes)	)	
			
	return hitnodes







def NodeDegreeDF(GH,hitnodes,NameVirusFile):
	
	# create an array of human proteins hit by the virus
	#hitnodes=[]
	#for i in nx.nodes(GV):
	#	if i.startswith('9606.')==True:
	#		hitnodes.append(i)
		

	
	# create an array of degree values of human proteins hit by the virus
	# (the indexing is the same than the array hitnodes)
	degree=[]
	for i in range(len(hitnodes)):
		degree.append(nx.degree(GH,hitnodes[i]))    
		
		
	# create dataframe of two columns: nodes hit by the virus and the corresponding degree
	columns_nodes_degree = {'nodes': hitnodes, 'degree': degree}
	ND = pd.DataFrame(data=columns_nodes_degree)        
	#len(ND)
	
	
	# if virus is not covid, we have to cancel those protein that are present in GV and absent in GH
	if NameVirusFile != 'ProteinNamesString_AS_OTHER_VIRUSES_TabSep.txt':
		
		# create an array with the indexes of nodes corresponding to human proteins that 
		# are absent in GH because of the initial score-based selection
		index_absent_node=[]
		for i in range(len(ND)): 
			if type(ND['degree'][i])!=int:
				index_absent_node.append(i)
	
    
    
		# erase rows of ND that correspond to absent nodes in GH
		ND = ND.drop(index_absent_node)        
		# reindex ND 
		ND.index = np.arange(0,len(ND),1)
		#len(ND)
	
	print('len(NodeDegreeDF): ',len(ND))
	
	
	return ND







# selection of hitnodes present in GH
def Hitnodes(GH,hitnodes,NameVirusFile):
	
	ND = NodeDegreeDF(GH,hitnodes,NameVirusFile)

	# recreate a new array of hit nodes that are all present in GH
	hitnodes = np.array(ND['nodes'])
	
	print('len(Hitnodes): ', len(hitnodes))
	
	return hitnodes









def NodeBetweennessDF(GH,BC,hitnodesBC):
	
	#hitnodesBC = Hitnodes(GH,hitnodes)
	
	# select only BC values of hit proteins			
	BCnodes = []
	BCbetweenness = []			
	for i in range(len(hitnodesBC)):
		BCrow = BC[BC['Nodi']==hitnodesBC[i]]
		BCnodes.append(BCrow.iloc[0][0])
		BCbetweenness.append(BCrow.iloc[0][1])
	
	
	# create a dataframe with hit nodes and corresponding betweenness centrality
	columns_nodes_BC = {'nodes': BCnodes, 'BC': BCbetweenness}
	BCselection = pd.DataFrame(data=columns_nodes_BC) 	

	# Sorting BCselection dataframe according to BC value (from largest to lowest)
	BC_sorted = BCselection.sort_values(by='BC', ascending=False)
	# reindex BC_sorted 
	BC_sorted.index = np.arange(0,len(BC_sorted),1)
	
	print('len(BC_sorted): ',len(BC_sorted))
	
	
	return BC_sorted









# RANDOM PERCOLATION 
        
def PercolationRandom(GH,ND):
	
	sizeG_single = len(max(nx.strongly_connected_components(GH)))
	sizeG=[] #create an empty array for the size of the giant component over time
	sizeG.append(sizeG_single)
	
	print('sizeG iniziale (RANDOM): ',sizeG_single)
	
	random_index = np.arange(len(ND)) # create an array of indexes for the nodes
	np.random.shuffle(random_index) # shuffle indexes
        
	for i in range(len(ND)):
		# remove nodes following the random order of random_index
		GH.remove_node(ND.iloc[random_index[i]][0])
		# calculate the size of the giant component
		sizeG_single = len(max(nx.strongly_connected_components(GH)))
		sizeG.append(sizeG_single) 
		
		
	print('sizeG finale (RANDOM): ',sizeG_single)
    
	return sizeG










# DEGREE-BASED PERCOLATION

def PercolationDegree(GH,ND,hitnodes):
	
	sizeG_single = len(max(nx.strongly_connected_components(GH)))
	sizeG=[] #create an empty array for the size of the giant component over time
	sizeG.append(sizeG_single)  

	print('sizeG iniziale (DEGREE): ',sizeG_single)       

	# node of maximum degree 
	node_to_remove = ND[ND['degree']==max(ND['degree'])]

	for j in range(len(hitnodes)):
		#remove the node with maximum degree
		GH.remove_node(node_to_remove.iloc[0][0]) 

	    # calculate the size of the giant component
		sizeG_single = len(max(nx.strongly_connected_components(GH)))
		sizeG.append(sizeG_single)
    
	    # erase from ND the row corresponding to the just erased node of GH
		index_removed_node = ND[ND['nodes']==node_to_remove.iloc[0][0]].index
		ND = ND.drop(index_removed_node)
		ND.index =np.arange(0,len(ND),1) # reindex ND
    
	    # calculate again the degree of every node
		if len(ND)>0:
			for i in range(len(ND)): 
				ND.iloc[i][1] = nx.degree(GH,ND.iloc[i][0]) 
			node_to_remove = ND[ND['degree']==max(ND['degree'])]
			
	print('sizeG finale (DEGREE): ',sizeG_single)
    
	
	return sizeG









#  BETWEENNESS CENTRALITY - BASED PERCOLATION 

def PercolationBetweenness(GH,BC_sorted):
	
	sizeG_single = len(max(nx.strongly_connected_components(GH)))
	sizeG = [] #create an empty array for the size of the giant component over time
	sizeG.append(sizeG_single) 
	
	print('sizeG iniziale (BC): ',sizeG_single)


	for i in range (len(BC_sorted)):
		# remove nodes following the descending order of betweenness
		GH.remove_node(BC_sorted.iloc[i][0])
		# calculate the size of the giant component
		sizeG_single = len(max(nx.strongly_connected_components(GH)))
		sizeG.append(sizeG_single) 
		
		
	print('sizeG finale (BC): ',sizeG_single)
		
		
	return sizeG
    










#   GRAPH

def GraphPercolation(percol,nameplot):
	sns.set_style('whitegrid')
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,5))
	x_axis = percol['removed_nodes']/percol['removed_nodes'].iloc[len(percol)-1]
	ax.plot(x_axis, percol['sizeG_random'], label='random', linewidth=0.75)
	ax.plot(x_axis, percol['sizeG_degree'], label='degree', linewidth=0.75)
	ax.plot(x_axis, percol['sizeG_BC'], label='betweenness', linewidth=0.75)
	ax.set_xlabel('% of removed nodes')
	ax.set_ylabel('Size of giant component')
	ax.legend(ncol=1 ,loc='best', fontsize=14)
	plt.savefig(nameplot)
	plt.show()
	
	return











def PercolationCode(GH_reference,ND,hitnodes,BC_sorted,FileNamePercolation,nameplot):	
	
	
	# recreate the original GH graph
	GH = GH_reference.copy(as_view=False)

	# RANDOM PERCOLATION
	sizeG_random = PercolationRandom(GH,ND)

	# recreate the original GH graph
	GH = GH_reference.copy(as_view=False)

	# DEGREE-BASED PERCOLATION
	sizeG_degree = PercolationDegree(GH,ND,hitnodes)

	# recreate the original GH graph
	GH = GH_reference.copy(as_view=False)

	#BETWEENNESS CENTRALITY - BASED PERCOLATION
	sizeG_BC = PercolationBetweenness(GH,BC_sorted)

	






	# create a dataframe with hit nodes and corresponding betweenness centrality
	columns_graph = {'sizeG_random': sizeG_random, 'sizeG_degree': sizeG_degree, 'sizeG_BC': sizeG_BC}
	graph_df = pd.DataFrame(data=columns_graph) 	


	#create file .txt with different values of the size of the giant component obtained by percolation
	graph_df.to_csv(FileNamePercolation, sep="\t", index=True, index_label = 'removed_nodes') 


	#import file as dataframe
	percol = pd.read_csv(FileNamePercolation, sep="\t")


	 



	#plot
	GraphPercolation(percol,nameplot)
	
	
	
	
	return 


