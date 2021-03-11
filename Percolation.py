import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt









'''
FUNZIONE CHE SVOLGE LA PERCOLATION SU UNA SINGOLA NETWORK
##############################################################################

INPUT:
	 - G: rete su cui applicare percolation
	 - centrality: dataframe with nodes, BC and degree as columns 
	 - network_name: title of the plot

OUTPUT: 
	 - percolation plot: size of the giant component vs % of removed nodes
	 - list of three arrays: each one is the size of the giant component at 
	   each node removal respectively in the case of node removal by BC, degree and randomly
	 

###############################################################################	 
'''
def Percolation(G, centrality, network_name):#, namefile
	
	
	network = G.copy(as_view=False)
	nodes = np.array(network.nodes())
	
	
	
	
	if nx.is_directed(network) == True:
		BCdataframe = centrality[0]
		BCdataframe = BCdataframe[['Nodes','BC']]
		DegreeDataframe = centrality[0]
		DegreeDataframe = DegreeDataframe[['Nodes','K']]
	else:
		BCdataframe = centrality[['Nodes','BC']]
		DegreeDataframe = centrality[['Nodes','K']]
		
	
	
	
	
	
	
	#####   Percolation BC
	
	
	# Sorting BCselection dataframe according to BC value (from highest to lowest)
	BC_sorted = BCdataframe.sort_values(by='BC', ascending=False)
	# reindex BC_sorted 
	BC_sorted.index = np.arange(0,len(BC_sorted),1)
	
	
	sizeBC=[]
	
	if nx.is_directed(network)==True:
		largest_cc = max(nx.weakly_connected_components(network), key=len)
	else:
		largest_cc = max(nx.connected_components(network), key=len)
	
	sizeBC.append(len(largest_cc))
	removedBC=0
	
	i=0
	while len(list(network.nodes))>1:
		node_to_remove=BC_sorted.iloc[i][0]
		network.remove_node(node_to_remove)
		removedBC=removedBC+1
		
		if nx.is_directed(network)==True:
			largest_cc = max(nx.weakly_connected_components(network), key=len)
		else:
			largest_cc = max(nx.connected_components(network), key=len)
		
		sizeBC.append(len(largest_cc))
		i=i+1
		#print("BC: Rimossi",removedBC, "nodi su", len(BC_sorted), ", size=",len(largest_cc))
	
	
	
	
	
	
	
	
	
	
	#####   Percolation Degree
	network=G.copy(as_view=False)
		
	degree_sorted=DegreeDataframe.sort_values(by='K', ascending=False)
	degree_sorted.index=np.arange(0,len(degree_sorted),1)
	
	sizeDegree=[]
	removedDegree=0
	
	if nx.is_directed(network)==True:
		largest_cc = max(nx.weakly_connected_components(network), key=len)
	else:
		largest_cc = max(nx.connected_components(network), key=len)
	
	sizeDegree.append(len(largest_cc))
	i=0
	
	while len(list(network.nodes))>1:
		
		node_to_remove=degree_sorted.iloc[i][0]
		
		network.remove_node(node_to_remove)
		removedDegree=removedDegree+1
		
		if nx.is_directed(network)==True:
			largest_cc = max(nx.weakly_connected_components(network), key=len)
		else:
			largest_cc = max(nx.connected_components(network), key=len)
	
		sizeDegree.append(len(largest_cc))
		
		i=i+1
		
		#print("Degree: Rimossi",removedDegree, "nodi","size=",len(largest_cc))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	#Percolation Random
	
	
	network=G.copy(as_view=False)
	removedRandom=[]
	removedRandom=0
	sizeRandom=[]
	
	if nx.is_directed(network)==True:
		largest_cc = max(nx.weakly_connected_components(network), key=len)
	else:
		largest_cc = max(nx.connected_components(network), key=len)
	
	sizeRandom.append(len(largest_cc))
	
	
	random_index = np.arange(len(nodes)) # create an array of indexes for the nodes
	np.random.shuffle(random_index)	
	
	i=0
	while len(list(network.nodes))>1:
		node_to_remove=nodes[random_index[i]]
		network.remove_node(node_to_remove)
		removedRandom=removedRandom+1
		if nx.is_directed(network)==True:
			largest_cc = max(nx.weakly_connected_components(network), key=len)
		else:
			largest_cc = max(nx.connected_components(network), key=len)
		sizeRandom.append(len(largest_cc))
		i=i+1
		#print("Random:Rimossi",removedRandom, "nodi su", len(nodes), ", size=",len(largest_cc))
	
	
	
	
	
	
	
	
	
	
	
	###	PLOT
	
	Removed_nodes=np.arange(1, max(len(sizeBC),len(sizeDegree),len(sizeRandom)+1), step=1)
	
	x_axis =[]
	for i in range (len(Removed_nodes)):
		x_axis.append(Removed_nodes[i]/max(Removed_nodes)*100)
	
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,5))#(4.5,3.5)
	ax.plot(x_axis, sizeBC, label='BC', linewidth=0.75)
	ax.plot(x_axis, sizeDegree, label='Degree', linewidth=0.75)
	ax.plot(x_axis, sizeRandom, label='Random', linewidth=0.75)
	ax.set_title(network_name)
	ax.set_xlabel('% of removed nodes',labelpad=2)
	ax.set_ylabel('Size of giant component',labelpad=2)
	ax.legend(ncol=1 ,loc='best', fontsize=12)
	
	plt.grid(True)
	#plt.savefig(namefile)

	plt.show()





	return  sizeBC, sizeDegree, sizeRandom
