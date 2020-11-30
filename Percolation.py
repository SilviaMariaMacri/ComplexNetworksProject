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
 
'''       
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

'''


def PercolationRandom(GH,ND):
	
	sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
	sizeG=[] #create an empty array for the size of the giant component over time
	sizeG.append(sizeG_single)
	
	NumberOfComponents = []
	NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
	NumberOfComponents.append(NumberOfComponents_single)
	
#	print('sizeG iniziale (RANDOM): ',sizeG_single)
#	print('numero di componenti iniziale: ',NumberOfComponents_single)
	
	random_index = np.arange(len(ND)) # create an array of indexes for the nodes
	np.random.shuffle(random_index) # shuffle indexes
     
	
#	print('numero di nodi iniziali: ',nx.number_of_nodes(GH)) 
	
	
	NumberRemovedNodes = []	#per far vedere quanti nodi sono rimossi a ogni iterazione
	
	# remove the first node following random_index
	node_to_remove = ND.iloc[random_index[0]][0]
			
	print('iterazione',1)  
#	print('nodo colpito: ',node_to_remove)

			
	#calculate first neighbors of the node before the removal
	first_neighbors = nx.neighbors(GH,node_to_remove)
	first_neighbors = list(first_neighbors)
	
	#select neighbors with links with combined_score > 900
	first_neighbors_selected = list()
	for i in range(len(first_neighbors)):
		if GH.edges[(node_to_remove,first_neighbors[i])]['weight'] > 900:
			first_neighbors_selected.append(first_neighbors[i])
			
			
	print('numero di primi vicini: ',len(first_neighbors_selected))
	
		
	
	
			
	GH.remove_node(node_to_remove)
#	print('numero di nodi rimasti: ',nx.number_of_nodes(GH)) 
	
	
	sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
	sizeG.append(sizeG_single)
	
	NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
	NumberOfComponents.append(NumberOfComponents_single)
	
	
	
	NumberRemovedNodes.append(len(node_to_remove))	# quenti nodi sono stati rimossi	
		

	#AGGIUNTO SOLO PER QUEST'ULYTIM prova:
#	#first_neighbors_selected = first_neighbors.copy()

	print(NumberOfComponents_single)

	#19254 è numero di nodi in GH 
	for i in range(1,19254):
		
		print('iterazione' , i+1)
			
		# nodes to remove are the first neighbors of the removed node and the second directly-hit node
		
		#calcolo i nodi in GH
		present_nodes_in_GH = nx.nodes(GH)
		present_nodes_in_GH = list(present_nodes_in_GH)
			
	

			
		node_to_remove =  first_neighbors_selected.copy()
		#hit node
		if i < len(ND):
			NDi = ND.iloc[random_index[i]][0]
			if NDi in present_nodes_in_GH: #funziona solo se NDi è in GH
				if NDi not in node_to_remove: 
					node_to_remove.append(NDi)
#		print('nodo colpito: ', NDi)
#		print('quanti nodi da rimuovere: ',len(node_to_remove))
		# node_to_remove è lista di nodi da rimuovere prima della prossima iterazione
		
		NumberRemovedNodes.append(len(node_to_remove))
		
		#calculate first neighbors of the nodes before the removal
		
		first_neighbors_selected = list() #lista definitiva che andrò a riempire con tutti i vicini
		
		#S = [950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100,50]
		#S = [900,800,700,600,500,400,300,200,100,0]
		S = np.arange(998,-2,-2)
	   
		for s in S:
			
			for k in range(len(node_to_remove)): 
				#calculate first neighbors for each just removed node
				first_neighbors_single = nx.neighbors(GH,node_to_remove[k]) #vicini a uno dei nodi che verranno rimossi
				first_neighbors_single = list(first_neighbors_single)
			
				#seleziono in base a score
				first_neighbors = list()
				for y in range(len(first_neighbors_single)):
					if GH.edges[(node_to_remove[k],first_neighbors_single[y])]['weight'] > s:
						first_neighbors.append(first_neighbors_single[y])
				
#aggiunto per l'occasione
#				first_neighbors	= first_neighbors_single.copy()
							
				#seleziono in modo che non si ripetano i nodi e che non siano all'interno
				#di 'first_neighbors_selected' e che non siano in 'node_to_remove'
				
				
				#for each element of the first neighbors...
				for j in range(len(first_neighbors)): 
					#append to the list 'first_neighbors_selected' only if the node is 
					#not already there and not in 'node_to_remove'
					if first_neighbors[j] not in first_neighbors_selected:
						if first_neighbors[j] not in node_to_remove:
							first_neighbors_selected.append(first_neighbors[j])
							
			if len(first_neighbors_selected) != 0:
				break
			
			
		print('numero di primi vicini: ',len(first_neighbors_selected))	
		
			
		
		if len(first_neighbors_selected) == 0:
			
			for z in range(len(node_to_remove)):
				GH.remove_node(node_to_remove[z])
			
				sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
				sizeG.append(sizeG_single) 
		
				NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
				NumberOfComponents.append(NumberOfComponents_single)
		
		
				NumberRemovedNodes.append(len(node_to_remove))	# quenti nodi sono stati rimossi	
	
			break

			
		#if len(node_to_remove) != 0:
		for z in range(len(node_to_remove)):
			GH.remove_node(node_to_remove[z])
		
#			print('numero di nodi rimasti: ',nx.number_of_nodes(GH)) 
	
		
			
		
		
		#if n == 0: #se numero di nodi in GH = 0 fermo il ciclo
		#	sizeG_single = 0
		#	sizeG.append(sizeG_single)
		#	break
		
		# calculate the size of the giant component
		
			if sizeG_single == 0:
				sizeG.append(sizeG_single)
				break

		
			sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
			sizeG.append(sizeG_single) 
		
			NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
			NumberOfComponents.append(NumberOfComponents_single)
		
			print(NumberOfComponents_single)
		
		
		
#			print('sizeG: ',sizeG_single)
		


    
	return sizeG,NumberOfComponents,NumberRemovedNodes


'''
SAVE

columns_graph = {'sizeG_random': PRandom}
graph_df = pd.DataFrame(data=columns_graph) 	

graph_df.to_csv('Percolation_InfluenzaA_sizeG.txt', sep="\t", index=True, index_label = 'removed_nodes') 






GHRidotta = GH.copy(as_view=False)
In: nx.info(GH)
Out: 'Name: \nType: DiGraph\nNumber of nodes: 11649\nNumber of edges: 261246
      \nAverage in degree:  22.4265\nAverage out degree:  22.4265'



'''





















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
'''
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


def PercolationBetweenness(GH,BC_sorted):
	
	sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
	sizeG=[] #create an empty array for the size of the giant component over time
	sizeG.append(sizeG_single)
	
	NumberOfComponents = []
	NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
	NumberOfComponents.append(NumberOfComponents_single)
	
	print('sizeG iniziale (RANDOM): ',sizeG_single)
	print('numero di componenti iniziale: ',NumberOfComponents_single)
	
	random_index = np.arange(len(BC_sorted)) # create an array of indexes for the nodes
	np.random.shuffle(random_index) # shuffle indexes
     
	
	print('numero di nodi iniziali: ',nx.number_of_nodes(GH)) 
	
	
	NumberRemovedNodes = []	#per far vedere quanti nodi sono rimossi a ogni iterazione
	
	# remove the first node following random_index
	node_to_remove = BC_sorted.iloc[random_index[0]][0]
			
	print('iterazione',1)  
	print('nodo colpito: ',node_to_remove)

			
	#calculate first neighbors of the node before the removal
	first_neighbors = nx.neighbors(GH,node_to_remove)
	first_neighbors = list(first_neighbors)
	
	#select neighbors with links with combined_score > 900
	first_neighbors_selected = list()
	for i in range(len(first_neighbors)):
		if GH.edges[(node_to_remove,first_neighbors[i])]['weight'] > 950:
			first_neighbors_selected.append(first_neighbors[i])
	
		
		
	
	print('quanti primi vicini ci sono: ',len(first_neighbors_selected))
			
	GH.remove_node(node_to_remove)
	print('numero di nodi rimasti: ',nx.number_of_nodes(GH)) 
	
	
	sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
	sizeG.append(sizeG_single)
	
	NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
	NumberOfComponents.append(NumberOfComponents_single)
	
	
	
	NumberRemovedNodes.append(len(node_to_remove))	# quenti nodi sono stati rimossi	
		
	for i in range(1,19254):
		
		print('iterazione' , i+1)
			
		# nodes to remove are the first neighbors of the removed node and the second directly-hit node
		
		#calcolo i nodi in GH
		present_nodes_in_GH = nx.nodes(GH)
		present_nodes_in_GH = list(present_nodes_in_GH)
			
				
		node_to_remove =  first_neighbors_selected.copy()
		#hit node
		if i < len(BC_sorted):
			NDi = BC_sorted.iloc[random_index[i]][0]
			if NDi in present_nodes_in_GH: #funziona solo se NDi è in GH
				if NDi not in node_to_remove: 
					node_to_remove.append(NDi)
		print('nodo colpito: ', NDi)
		print('quanti nodi da rimuovere: ',len(node_to_remove))
		# node_to_remove è lista di nodi da rimuovere prima della prossima iterazione
		
		NumberRemovedNodes.append(len(node_to_remove))
		
		#calculate first neighbors of the nodes before the removal
		
		first_neighbors_selected = list() #lista definitiva che andrò a riempire con tutti i vicini
		
		#S = [950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100,50]
		S = [800,600,400,200,0]
	   
		for s in S:
			
			for k in range(len(node_to_remove)): 
				#calculate first neighbors for each just removed node
				first_neighbors_single = nx.neighbors(GH,node_to_remove[k]) #vicini a uno dei nodi che verranno rimossi
				first_neighbors_single = list(first_neighbors_single)
			
				#seleziono in base a score
				first_neighbors = list()
				for y in range(len(first_neighbors_single)):
					if GH.edges[(node_to_remove[k],first_neighbors_single[y])]['weight'] > s:
						first_neighbors.append(first_neighbors_single[y])
					
							
				#seleziono in modo che non si ripetano i nodi e che non siano all'interno
				#di 'first_neighbors_selected' e che non siano in 'node_to_remove'
				
				
				#for each element of the first neighbors...
				for j in range(len(first_neighbors)): 
					#append to the list 'first_neighbors_selected' only if the node is 
					#not already there and not in 'node_to_remove'
					if first_neighbors[j] not in first_neighbors_selected:
						if first_neighbors[j] not in node_to_remove:
							first_neighbors_selected.append(first_neighbors[j])
							
			if len(first_neighbors_selected) != 0:
				break

		
		print('quanti primi vicini ci sono: ',len(first_neighbors_selected))	

		if len(first_neighbors_selected) == 0:
			
			for z in range(len(node_to_remove)):
				GH.remove_node(node_to_remove[z])
			
			sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
			sizeG.append(sizeG_single) 
		
			NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
			NumberOfComponents.append(NumberOfComponents_single)
		
		
			NumberRemovedNodes.append(len(node_to_remove))	# quenti nodi sono stati rimossi	
	
			break

			
		#if len(node_to_remove) != 0:
		for z in range(len(node_to_remove)):
			GH.remove_node(node_to_remove[z])
		
		print('numero di nodi rimasti: ',nx.number_of_nodes(GH)) 
	
		
			
		
		
		#if n == 0: #se numero di nodi in GH = 0 fermo il ciclo
		#	sizeG_single = 0
		#	sizeG.append(sizeG_single)
		#	break
		
		# calculate the size of the giant component
		
		if sizeG_single == 0:
			sizeG.append(sizeG_single)
			break

		
		sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
		sizeG.append(sizeG_single) 
		
		NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
		NumberOfComponents.append(NumberOfComponents_single)
		
		
		
		
		
		print('sizeG: ',sizeG_single)
		


    
	return sizeG,NumberOfComponents,NumberRemovedNodes


'''


def PercolationBetweenness(GH,BC_sorted):
	
	sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
	sizeG=[] #create an empty array for the size of the giant component over time
	sizeG.append(sizeG_single)
	
	NumberOfComponents = []
	NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
	NumberOfComponents.append(NumberOfComponents_single)
	
#	print('sizeG iniziale (RANDOM): ',sizeG_single)
#	print('numero di componenti iniziale: ',NumberOfComponents_single)
	
	random_index = np.arange(len(BC_sorted)) # create an array of indexes for the nodes
	np.random.shuffle(random_index) # shuffle indexes
     
	
#	print('numero di nodi iniziali: ',nx.number_of_nodes(GH)) 
	
	
	NumberRemovedNodes = []	#per far vedere quanti nodi sono rimossi a ogni iterazione
	
	# remove the first node following random_index
	node_to_remove = BC_sorted.iloc[random_index[0]][0]
			
	print('iterazione',1)  
#	print('nodo colpito: ',node_to_remove)

			
	#calculate first neighbors of the node before the removal
	first_neighbors = nx.neighbors(GH,node_to_remove)
	first_neighbors = list(first_neighbors)
	
	#select neighbors with links with combined_score > 900
	first_neighbors_selected = list()
	for i in range(len(first_neighbors)):
		if GH.edges[(node_to_remove,first_neighbors[i])]['weight'] > 900:
			first_neighbors_selected.append(first_neighbors[i])
			
			
	print('numero di primi vicini: ',len(first_neighbors_selected))
	
		
	
	
			
	GH.remove_node(node_to_remove)
#	print('numero di nodi rimasti: ',nx.number_of_nodes(GH)) 
	
	
	sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
	sizeG.append(sizeG_single)
	
	NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
	NumberOfComponents.append(NumberOfComponents_single)
	
	
	
	NumberRemovedNodes.append(len(node_to_remove))	# quenti nodi sono stati rimossi	
		

	#AGGIUNTO SOLO PER QUEST'ULYTIM prova:
#	#first_neighbors_selected = first_neighbors.copy()

	print(NumberOfComponents_single)

	#19254 è numero di nodi in GH 
	for i in range(1,19254):
		
		print('iterazione' , i+1)
			
		# nodes to remove are the first neighbors of the removed node and the second directly-hit node
		
		#calcolo i nodi in GH
		present_nodes_in_GH = nx.nodes(GH)
		present_nodes_in_GH = list(present_nodes_in_GH)
			
	

			
		node_to_remove =  first_neighbors_selected.copy()
		#hit node
		if i < len(BC_sorted):
			NDi = BC_sorted.iloc[random_index[i]][0]
			if NDi in present_nodes_in_GH: #funziona solo se NDi è in GH
				if NDi not in node_to_remove: 
					node_to_remove.append(NDi)
#		print('nodo colpito: ', NDi)
#		print('quanti nodi da rimuovere: ',len(node_to_remove))
		# node_to_remove è lista di nodi da rimuovere prima della prossima iterazione
		
		NumberRemovedNodes.append(len(node_to_remove))
		
		#calculate first neighbors of the nodes before the removal
		
		first_neighbors_selected = list() #lista definitiva che andrò a riempire con tutti i vicini
		
		#S = [950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100,50]
		#S = [900,800,700,600,500,400,300,200,100,0]
		S = np.arange(998,-2,-2)
	   
		for s in S:
			
			for k in range(len(node_to_remove)): 
				#calculate first neighbors for each just removed node
				first_neighbors_single = nx.neighbors(GH,node_to_remove[k]) #vicini a uno dei nodi che verranno rimossi
				first_neighbors_single = list(first_neighbors_single)
			
				#seleziono in base a score
				first_neighbors = list()
				for y in range(len(first_neighbors_single)):
					if GH.edges[(node_to_remove[k],first_neighbors_single[y])]['weight'] > s:
						first_neighbors.append(first_neighbors_single[y])
				
#aggiunto per l'occasione
#				first_neighbors	= first_neighbors_single.copy()
							
				#seleziono in modo che non si ripetano i nodi e che non siano all'interno
				#di 'first_neighbors_selected' e che non siano in 'node_to_remove'
				
				
				#for each element of the first neighbors...
				for j in range(len(first_neighbors)): 
					#append to the list 'first_neighbors_selected' only if the node is 
					#not already there and not in 'node_to_remove'
					if first_neighbors[j] not in first_neighbors_selected:
						if first_neighbors[j] not in node_to_remove:
							first_neighbors_selected.append(first_neighbors[j])
							
			if len(first_neighbors_selected) != 0:
				break
			
			
		print('numero di primi vicini: ',len(first_neighbors_selected))	
		
			
		
		if len(first_neighbors_selected) == 0:
			
			for z in range(len(node_to_remove)):
				GH.remove_node(node_to_remove[z])
			
				sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
				sizeG.append(sizeG_single) 
		
				NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
				NumberOfComponents.append(NumberOfComponents_single)
		
		
				NumberRemovedNodes.append(len(node_to_remove))	# quenti nodi sono stati rimossi	
	
			break

			
		#if len(node_to_remove) != 0:
		for z in range(len(node_to_remove)):
			GH.remove_node(node_to_remove[z])
		
#			print('numero di nodi rimasti: ',nx.number_of_nodes(GH)) 
	
		
			
		
		
		#if n == 0: #se numero di nodi in GH = 0 fermo il ciclo
		#	sizeG_single = 0
		#	sizeG.append(sizeG_single)
		#	break
		
		# calculate the size of the giant component
		
			if sizeG_single == 0:
				sizeG.append(sizeG_single)
				break

		
			sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
			sizeG.append(sizeG_single) 
		
			NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
			NumberOfComponents.append(NumberOfComponents_single)
		
			print(NumberOfComponents_single)
		
		
		
#			print('sizeG: ',sizeG_single)
		


    
	return sizeG,NumberOfComponents,NumberRemovedNodes









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


