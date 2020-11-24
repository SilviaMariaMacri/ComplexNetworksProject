'''
HOMO SAPIENS = 9606.
Human immunodeficiency virus type 1 group M subtype B (isolate HXB2) = 11676.
Influenza A virus (strain A/Puerto Rico/8/1934 H1N1) = 11320.
'''


'AVVISO: FILE DATI E FILE PERCOLATION.PY DOVREBBERO STARE SULLA STESSA CARTELLA
'        OPPURE CARICARE PRIMA FILE SAPIENS E POI SPOSTARSI NELLA CARTELLA DI PERCOLATION.PY


import os 
import pandas as pd


#directory = 'C:\Users\silvi\Downloads'
#directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/dati'
directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/data'
os.chdir(directory) 





import Percolation






'''
############################ IMPORT DATA #####################################
'''

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
virus = pd.read_csv('string_interactions_prova.tsv', sep="\t", 
                   usecols=['node1_external_id','node2_external_id', 'combined_score'])








'''
############################### CREATE GRAPHS ################################
'''


# create a directed graph of human PPI
GH = Percolation.HumanGraph(sapiens)

# create a copy of GH to store as a reference     
GH_reference = GH.copy(as_view=False)

# create a directed graph of virus-human interactions 
GV = Percolation.VirusGraph(virus,'11676.')









'''
################### define dataframes and array of hit nodes, degree and BC 
'''

ND = Percolation.NodeDegreeDF(GH,GV)
hitnodes = Percolation.Hitnodes(GH,GV)

#import BC file as dataframe
BC = pd.read_csv('BetweennessCentrality2.csv', sep=",", skiprows=1)
  
BC_sorted = Percolation.NodeBetweennessDF(GH,GV,BC)















'''
############################ PERCOLATION #####################################
'''




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







'''
################################ PLOT #########################################
'''


# create a dataframe with hit nodes and corresponding betweenness centrality
columns_graph = {'sizeG_random': sizeG_random, 'sizeG_degree': sizeG_degree, 'sizeG_BC': sizeG_BC}
graph_df = pd.DataFrame(data=columns_graph) 	


#create file .txt with different values of the size of the giant component obtained by percolation
graph_df.to_csv('percolation_nomevirus.txt', index=True, index_label = 'removed_nodes') 


#import file as dataframe
percol = pd.read_csv('percolation_nomevirus.txt', sep=",")




Percolation.GraphPercolation(percol)





