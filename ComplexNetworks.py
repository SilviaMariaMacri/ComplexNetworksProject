'''
HOMO SAPIENS = 9606.
Human immunodeficiency virus type 1 group M subtype B (isolate HXB2) = 11676.
Influenza A virus (strain A/Puerto Rico/8/1934 H1N1) = 11320.



COMMENTI:
- riga 164: a me funziona, cmq ho lasciato anche il ciclo for, è equivalente 
- intorno a 285: non mi funzionava la selezione della betweenness. 
  credo perchè dopo aver creato il dataframe vuoto non c'era un comando che 
  gli facesse aggiungere una riga ad ogni ciclo. Non sono riuscita a trovarlo 
  quindi ho messo quelle altre righe che dovrebbero funzionare lo stesso
- grafico: ho messo una funzione che fa tutto, l'ho copiata da un'altro programma,
  dovrebbe fare grafici un po' più carini 
  PROBLEMA. grafico della percolation nell'articolo: perchè la size normalizzata va a zero?' 
  






DUBBI E PROBLEMI: 
- tenere solo le librerie che servono
- vedi se si possono leggere direttamente i dati da internet
- matrice fatta come i files ha un nome?
- per non considerare nodi presenti in virus ma non in sapiens abbiamo calcolato il degree. 
  va bene effettuare questo criterio anche prima della random selection (il degree non servirebbe)
  o c'è un altro modo migliore?? 
- ho provato a fare un grafico un po' di fretta nel caso della percolation con degree    
  ed è venuta una retta. secondo me dovremmo provare a includere più nodi, come dicevamo ieri

-se considerassimo anche iteratori dalla seconda shell?
'''
#%%
import os 
import numpy as np
#import scipy.stats as st
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
#from scipy.integrate import odeint
import networkx as nx


#directory = 'C:\Users\silvi\Downloads'
#directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/dati'
directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/data'
os.chdir(directory) 



#%% import data

# Dataframes:  
# each row corresponds to a link
# 3 columns:
#    1. starting node    
#    2. arrival node
#    3. score of the link


# create datadrame of human PPI links
sapiens = pd.read_csv('9606.protein.links.v11.0.txt', sep=" ", 
                    usecols=['protein1', 'protein2', 'combined_score']) 
#select human PPI links depending on the score
sapiens=sapiens[sapiens['combined_score']>400]


# create dataframe of virus-human interactions
virus = pd.read_csv('string_interactions_InfluenzaA.tsv', sep="\t", 
                   usecols=['node1_external_id','node2_external_id', 'combined_score'])





 
#%%

# create a directed graph of human PPI
GH = nx.DiGraph() 
for i in range(len(sapiens)): 
    GH.add_edge(sapiens.iloc[i][0], sapiens.iloc[i][1], weight=sapiens.iloc[i][2])

'''
nx.info(GH)
Out: 'Name: \nType: DiGraph\nNumber of nodes: 19254\nNumber of edges: 
    1985196\nAverage in degree: 103.1056\nAverage out degree: 103.1056'''


# create a copy of GH to store as a reference     
GH_reference = GH.copy(as_view=False)


#%% 

# select only links with virus protein as starting node and human protein as arrival node 
interactionsVH = []    
for i in range(len(virus)):
    if virus.iloc[i,0].startswith('11320.')==True & virus.iloc[i,1].startswith('9606.')==True:
        interactionsVH.append(virus.iloc[i])
        

#create a directed graph of virus-human interactions        
GV = nx.DiGraph() 
for i in range(len(interactionsVH)): 
    GV.add_edge(interactionsVH[i][0], interactionsVH[i][1], weight=interactionsVH[i][2])
        



         
#%% 

# create an array of human proteins hit by the virus
hitnodes=[]
for i in nx.nodes(GV):
    if i.startswith('9606.')==True:
        hitnodes.append(i)
        


        
        


# create an array of degree values of human proteins hit by the virus
# (the indexing is the same than the array hitnodes)
degree=[]
for i in range(len(hitnodes)):
    degree.append(nx.degree(GH,hitnodes[i]))    
        
        

     
# create dataframe of two columns: nodes hit by the virus and the corresponding degree
columns_nodes_degree = {'nodes': hitnodes, 'degree': degree}
ND = pd.DataFrame(data=columns_nodes_degree)        
#len(ND)


# create an array with the indexes of nodes corresponding to human proteins that 
# are absent in GH because of the initial score-based selection
index_absent_node=[]
for i in range(len(ND)): 
    if type(ND['degree'][i])!=int:
        index_absent_node.append(i)
#len(indici)
    
    
# erase rows of ND that correspond to absent nodes in GH
ND = ND.drop(index_absent_node)        
# reindex ND 
ND.index = np.arange(0,len(ND),1)
#len(ND)


# recreate a new array of hit nodes that are all present in GH
#for i in range (len(ND)):
#    hitnodes.append(np.array(ND['nodes'])[i])

# recreate a new array of hit nodes that are all present in GH
hitnodes = np.array(ND['nodes'])


#%% RANDOM PERCOLATION 

        
    

sizeG_random=[] #create an empty array for the size of the giant component over time
# size of the giant component before starting the percolation
size_G_single = len(max(nx.strongly_connected_components(GH)))
sizeG_random.append(size_G_single) 

        
random_index = np.arange(len(ND)) # create an array of indexes for the nodes
np.random.shuffle(random_index) # shuffle indexes
        
for i in range(len(ND)): 
    # remove nodes following the random order of random_index 
    GH.remove_node(ND.iloc[random_index[i]][0])
    #print('removed node (with degree): ') 
    #print(ND.iloc[random_index[i]])
    
    # calculate the size of the giant component
    size_G_single = len(max(nx.strongly_connected_components(GH)))
    sizeG_random.append(size_G_single) 
    #print('size of giant component: ')
    #print(len(max(nx.strongly_connected_components(GH))))
   


#%% DEGREE-BASED PERCOLATION


# recreate the original GH graph
GH = GH_reference.copy(as_view=False)




 
sizeG_degree=[] #create an empty array for the size of the giant component over time
# size of the giant component before starting the percolation
size_G_single = len(max(nx.strongly_connected_components(GH)))
sizeG_degree.append(size_G_single)         



# node of maximum degree 
node_to_remove = ND[ND['degree']==max(ND['degree'])]

for j in range(len(hitnodes)): 
    #remove the node with maximum degree
    GH.remove_node(node_to_remove.iloc[0][0]) 
    #print('removed node (with degree): ') 
    #print(node_to_remove.iloc[0])
    
    # calculate the size of the giant component
    size_G_single = len(max(nx.strongly_connected_components(GH)))
    sizeG_degree.append(size_G_single)
    #print('size of giant component: ')
    #print(len(max(nx.strongly_connected_components(GH))))
    
    # erase from ND the row corresponding to the just erased node of GH
    index_removed_node = ND[ND['nodes']==node_to_remove.iloc[0][0]].index
    ND = ND.drop(index_removed_node)  
    ND.index =np.arange(0,len(ND),1) # reindex ND
    
    # calculate again the degree of every node
    for i in range(len(ND)):
        ND.iloc[i][1] = nx.degree(GH,ND.iloc[i][0]) 
    node_to_remove = ND[ND['degree']==max(ND['degree'])]
    #print('node to remove in the next iteration: ',node_to_remove.iloc[0])







#%%  BETWEENNESS CENTRALITY - BASED PERCOLATION 



'''
#calculate BC
BC = nx.betweenness_centrality(GH, weight ='weights') 

#convert dict to dataframe
BCdf = pd.DataFrame.from_dict(BC,orient='index')

#create file .txt with BC values
BCdf.to_csv('BC.txt', index=True) 
'''






# recreate the original GH graph
GH = GH_reference.copy(as_view=False)





#import BC file as dataframe
BC = pd.read_csv('BetweennessCentrality2.csv', sep=",", skiprows=1)
  


'''
#Selecting only Bc values of hit protein
BCselection=pd.DataFrame()
for i in range (len(BC)):
    for j in range (len(hitnodes)):
        if BC.iloc[i][0]==hitnodes[j]:
            BCselection=BC.iloc[i]'''	
			
# select only BC values of hit proteins			
BCnodes = []
BCbetweenness = []			
for i in range(len(hitnodes)):
	BCrow = BC[BC['Nodi']==hitnodes[i]]
	BCnodes.append(BCrow.iloc[0][0])
	BCbetweenness.append(BCrow.iloc[0][1])
	
	
# create a dataframe with hit nodes and corresponding betweenness centrality
columns_nodes_BC = {'nodes': BCnodes, 'BC': BCbetweenness}
BCselection = pd.DataFrame(data=columns_nodes_BC) 	


			


# Sorting BCselection dataframe according to BC value (from largest to lowest)
BC_sorted = BCselection.sort_values(by='BC', ascending=False)
# reindex BC_sorted 
BC_sorted.index = np.arange(0,len(BC_sorted),1)






sizeG_BC = [] #create an empty array for the size of the giant component over time
# size of the giant component before starting the percolation
size_G_single = len(max(nx.strongly_connected_components(GH)))
sizeG_BC.append(size_G_single) 


for i in range (len(BC_sorted)):
	# remove nodes following the descending order of betweenness
    GH.remove_node(BC_sorted.iloc[i][0])
    # calculate the size of the giant component
    size_G_single = len(max(nx.strongly_connected_components(GH)))
    sizeG_BC.append(size_G_single) 
    




#%%   GRAPH

# create a dataframe with hit nodes and corresponding betweenness centrality
columns_graph = {'sizeG_random': sizeG_random, 'sizeG_degree': sizeG_degree, 'sizeG_BC': sizeG_BC}
graph_df = pd.DataFrame(data=columns_graph) 	


#create file .txt with different values of the size of the giant component obtained by percolation
graph_df.to_csv('percolation_nomevirus.txt', index=True) 


#import BC file as dataframe
percol = pd.read_csv('percolation_nomevirus.txt', sep=",")
#RICORDARE DI RINOMINARE LA PRIMA COLONNA DEL FILE


'''   
removed_nodes=np.arange(len(hitnodes)+1)
plt.plot(removed_nodes,sizeG_random, 'ro')
plt.plot(removed_nodes,sizeG_degree, 'bo')
plt.plot(removed_nodes,sizeG_BC, 'go')
plt.show()
'''


def graph_percolation(percol):
	sns.set_style('whitegrid')
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
	x_axis = percol['removed_nodes']/percol['removed_nodes'].iloc[len(percol)-1]
	ax.plot(x_axis, percol['sizeG_random']/percol['sizeG_random'][0], label='random', linewidth=0.75)
	ax.plot(x_axis, percol['sizeG_degree']/percol['sizeG_degree'][0], label='degree', linewidth=0.75)
	ax.plot(x_axis, percol['sizeG_BC']/percol['sizeG_BC'][0], label='betweenness', linewidth=0.75)
	ax.set_xlabel('% of removed nodes')
	ax.set_ylabel('Normalized size of giant component')
	ax.legend(ncol=1 ,loc='best', fontsize=14)
	plt.show()
	
	return

graph_percolation(percol)

