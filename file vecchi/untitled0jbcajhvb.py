





G = nx.cycle_graph(4, create_using=nx.DiGraph())
nx.add_cycle(G, [10, 11, 12])
nx.add_cycle(G, [4, 5, 12])
nx.add_cycle(G, [3, 10, 1])


nx.draw(G,with_labels=True)


a = nx.neighbors(G, 10)
a=list(a)
b = nx.all_neighbors(G, 10)
b=list(b)

a
b

#%%

G = nx.gnp_random_graph(40, 0.2, seed=1, directed=True)
G = nx.erdos_renyi_graph(40,0.2)
G = nx.to_directed(G)


nx.draw(G,with_labels=True)

columns_nodes_degree = {'nodes': [1,4,7,12], 'degree': [3,2,1,3]}
ND = pd.DataFrame(data=columns_nodes_degree)  
	
	
	
G = nx.DiGraph() 
G.add_edge(1, 2, weight=950)
G.add_edge(1, 3, weight=400)
G.add_edge(2, 4, weight=950)
G.add_edge(2, 5, weight=400)
G.add_edge(3, 6, weight=950)
G.add_edge(3, 7, weight=400)
G.add_edge(4, 8, weight=950)
G.add_edge(4, 9, weight=400)
G.add_edge(5, 10, weight=950)
G.add_edge(5, 11, weight=400)
G.add_edge(6, 12, weight=950)
G.add_edge(6, 13, weight=400)
G.add_edge(7, 14, weight=950)
G.add_edge(7, 15, weight=400)


nx.draw(G,with_labels=True)	
	
	
	
	
	
	
#%%

import pandas as pd


#def PercolationRandom(GH,ND):
	
sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
sizeG=[] #create an empty array for the size of the giant component over time
sizeG.append(sizeG_single)
	
NumberOfComponents = []
NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
NumberOfComponents.append(NumberOfComponents_single)
	
print('sizeG iniziale (RANDOM): ',sizeG_single)
print('numero di componenti iniziale: ',NumberOfComponents_single)
	
random_index = np.arange(len(ND)) # create an array of indexes for the nodes
np.random.shuffle(random_index) # shuffle indexes
     
	
print('numero di nodi iniziali: ',nx.number_of_nodes(GH)) 
	
	
NumberRemovedNodes = []	#per far vedere quanti nodi sono rimossi a ogni iterazione
	
# remove the first node following random_index
node_to_remove = ND.iloc[random_index[0]][0]
			
print('iterazione',1)  
print('nodo colpito: ',node_to_remove)

			
#calculate first neighbors of the node before the removal
first_neighbors = nx.neighbors(GH,node_to_remove)
first_neighbors = list(first_neighbors)
	
#select neighbors with links with combined_score > 900
first_neighbors_selected = list()
for i in range(len(first_neighbors)):
	if GH.edges[(node_to_remove,first_neighbors[i])]['weight'] > 800:
		first_neighbors_selected.append(first_neighbors[i])
	
		
		
	
print('quanti primi vicini ci sono: ',len(first_neighbors_selected))
			
GH.remove_node(node_to_remove)
print('numero di nodi rimasti: ',nx.number_of_nodes(GH)) 
	
	
sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
sizeG.append(sizeG_single)
	
NumberOfComponents_single = [len(c) for c in sorted(nx.strongly_connected_components(GH), key=len, reverse=True)]
NumberOfComponents.append(NumberOfComponents_single)
	
	
	
NumberRemovedNodes.append(len(node_to_remove))	# quenti nodi sono stati rimossi	
		




#%%
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
	print('nodo colpito: ', NDi)
	print('quanti nodi da rimuovere: ',len(node_to_remove))
	# node_to_remove è lista di nodi da rimuovere prima della prossima iterazione
		
	NumberRemovedNodes.append(len(node_to_remove))
		
	#calculate first neighbors of the nodes before the removal
		
	first_neighbors_selected = list() #lista definitiva che andrò a riempire con tutti i vicini
		




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
			continue
		
	print('quanti primi vicini ci sono: ',len(first_neighbors_selected))	



			
	#if len(node_to_remove) != 0:
	for z in range(len(node_to_remove)):
		GH.remove_node(node_to_remove[z])
		
	print('numero di nodi rimasti: ',nx.number_of_nodes(GH)) 
	
		
			
		
		
		#if n == 0: #se numero di nodi in GH = 0 fermo il ciclo
		#	sizeG_single = 0
		#	sizeG.append(sizeG_single)
		#	break
		
		# calculate the size of the giant component
		
	sizeG_single = len(max(nx.strongly_connected_components(GH), key=len))
	sizeG.append(sizeG_single) 
		
	NumberOfComponents_single = len(nx.strongly_connected_components(GH))
	NumberOfComponents.append(NumberOfComponents_single)
	
		
	NumberRemovedNodes.append(len(node_to_remove))	# quenti nodi sono stati rimossi	
	
		
		
	print('sizeG: ',sizeG_single)








sizeG,NumberOfComponents,NumberRemovedNodes


columns_graph = {'sizeG_random': sizeG_random[0], 'nodi tolti: ': sizeG_random[2]}
graph_df = pd.DataFrame(data=columns_graph)	
graph_df.to_csv('random.txt', sep="\t", index=True, )


columns_graph = {'sizeG_BC': sizeG_BC[0], 'nodi tolti: ': sizeG_BC[2]}
graph_df = pd.DataFrame(data=columns_graph)	
graph_df.to_csv('BC.txt', sep="\t", index=True)




	#import file as dataframe
percol_random = pd.read_csv('random.txt', sep="\t")
percol_BC = pd.read_csv('BC.txt', sep="\t")

	 






fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,5))

ax.scatter(asse_x, percol_random['b'], label='random', linewidth=0.75)
#ax.plot(x_axis, percol['sizeG_degree'], label='degree', linewidth=0.75)
ax.scatter(asse_x_BC, percol_BC['b'], label='betweenness', linewidth=0.75)
ax.set_xlabel('% of removed nodes')
ax.set_ylabel('Size of giant component')
ax.legend(ncol=1 ,loc='best', fontsize=14)
#plt.savefig(nameplot)
#plt.show()
	
asse_x = []
#asse_x.append(19254)
sss = percol_random.iloc[0,2]
asse_x.append(sss)
for s in range(1,len(percol_random)):
	sss = sss + percol_random.iloc[i,2]
	asse_x.append(sss) 
	
	




percol_random = pd.read_csv('Percolation_InfluenzaA_Random_sizeG.txt', sep="\t")
percol_BC = pd.read_csv('Percolation_InfluenzaA_BC_sizeG.txt', sep="\t")
	

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,5))
#x_axis_random = percol_random['removed_nodes']#/percol_random['removed_nodes'].iloc[len(percol)-1]
x_axis_BC = percol_BC['removed_nodes']#/percol_BC['removed_nodes'].iloc[len(percol)-1]
#ax.plot(x_axis_random, percol['sizeG_random'], label='random', linewidth=0.75)
#ax.plot(x_axis, percol['sizeG_degree'], label='degree', linewidth=0.75)
ax.plot(x_axis_BC, percol_BC['sizeG_BC'], label='betweenness', linewidth=0.75)
ax.set_xlabel('number of removed nodes')#'% of removed nodes')
ax.set_ylabel('Size of giant component')
ax.legend(ncol=1 ,loc='best', fontsize=14)

	
	
	
	# create a dataframe with hit nodes and corresponding betweenness centrality
columns_graph = {'sizeG_degree': sizeG}
graph_df = pd.DataFrame(data=columns_graph) 	

graph_df.to_csv('Percolation_InfluenzaA_sizeG.txt', sep="\t", index=True, index_label = 'removed_nodes') 
	
	