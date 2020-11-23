#%%
import os 
import numpy as np
#import scipy.stats as st
#import pylab as plt
import pandas as pd
#import seaborn as sns
#from scipy.integrate import odeint
import networkx as nx
import matplotlib.pyplot as plt


directory = '/home/caterina/Scaricati'
os.chdir(directory)

#leggo da files e creo due dataframe: uomo e virus

uomo= pd.read_csv('9606.protein.links.v11.0.txt', sep=" ", 
					usecols=['protein1', 'protein2', 'combined_score']) 
#seleziono link HUMAN in base a score
uomo=uomo[uomo['combined_score']>400]
uomo.index = np.arange(0, len(uomo), 1) 
#%%
virus = pd.read_csv('Human_metapneumovirus_score0400_it500.tsv', sep="\t", 
				   usecols=['node1_external_id','node2_external_id', 'combined_score'])


#%%creo graph HUMAN
GH = nx.DiGraph() 
for i in range(len(uomo)): 
	GH.add_edge(uomo.iloc[i][0], uomo.iloc[i][1], weight=uomo.iloc[i][2])
GHcopy=GH.copy(as_view=False)
#%%seleziono da dataframe virus solo interazioni virus-human
interazioniVH = []	
for i in range(len(virus)):
	if virus.iloc[i,0].startswith('162145.')==True & virus.iloc[i,1].startswith('9606.')==True:
		interazioniVH.append(virus.iloc[i]) #aggiungo la riga selezionata da virus
		
#%%#creo graph virus-human interactions		
GV = nx.DiGraph() 
for i in range(len(interazioniVH)): 
	GV.add_edge(interazioniVH[i][0], interazioniVH[i][1], weight=interazioniVH[i][2])
		
	 	
#%%#creo array con nodi HUMAN colpiti da virus
nodi=[]
for i in nx.nodes(GV):
	if i.startswith('9606.')==True:
		nodi.append(i)
	
#%%	
##############CRITERIO PER RIORDINO NODI = DEGREE	
	 
#creo array di degree dei nodi colpiti da virus 
grado=[]
for i in range(len(nodi)):
	grado.append( nx.degree(GH,nodi[i])) #degree è calcolato da rete umana
	
	
#creo dataframe di due colonne nodi-degree
nodidegree = {'nodi': nodi, 'degree': grado} #questo è un dizionario di liste
ND = pd.DataFrame(data=nodidegree)		


#pongo uguale a ZERO il degree dei nodi presenti in virus ma non in human
for i in range(len(ND)):
	if type(ND.iloc[i][1])!=int:
		ND.iloc[i][1]=0

#%%
##############CRITERIO PER RIORDINO NODI = BETWEENNESS CENTRALITY	      
node_list=list() #serve solo per togliere gli zeri da quell'altra
for i in range (len (ND)):
    if ND.iloc[i][1]!=0:
        node_list.append(ND.iloc[i][0])

BC=nx.betweenness_centrality(GH, weight='weight')

#%%TRASFORMO BC DA DICT A DATAFRAME
BCdat=pd.DataFrame.from_dict(BC,orient='index',columns=['BC'])
#salvo in file csv
BCdat.to_csv('BetweennessCentrality1.csv', index=True)
