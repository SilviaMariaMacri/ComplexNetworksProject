import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
import seaborn as sns



'''
 INPUT:
 G = directed or undirected network

 OUTPUT:
 if G undirected: 
 	df = dataframe with degree, degree of nearest neighbors, betweenness 
        centrality and closeness centrality as columns and each row corresponding 
        to a node of the network	
 if G directed:
	df = array of three dataframes as entries: 
	df[0] = dataframe with In degree, Out degree, In degree of nearest neighbors, 
            Out degree of nearest neighbors, betweenness centrality and closeness 
            centrality as columns and each row corresponding to a node 
    df[1] = dataframe with same columns of df[0] and row corresponding only to 
            human nodes
    df[2] = dataframe with same columns of df[0] and row corresponding only to 
            viral nodes
'''


def NetworkCharacterization(G):
	
	if nx.is_directed(G) == False:
		
		degree = pd.DataFrame.from_dict(G.degree(weight='weight'))
		degree.columns=['Nodes','K']
		
		#degNN = pd.DataFrame.from_dict(nx.average_neighbor_degree(G, weight='weight'),orient='index',columns=['K_nn'])
		#degNN = degNN.set_index(np.arange(0,len(degNN),1))
		
		#BC = pd.DataFrame.from_dict(nx.betweenness_centrality(G,weight='weight'),orient='index',columns=['BC']) #è già normalizzata
		#BC = BC.set_index(np.arange(0,len(BC),1))
		
		#CL = pd.DataFrame.from_dict(nx.closeness_centrality(G),orient='index',columns=['CL']) #è già normalizzata
		#CL = CL.set_index(np.arange(0,len(CL),1))
		
	
		#dataframe
		#df = pd.concat([degree,degNN,BC,CL], axis=1)
		df = degree.copy()
		
		
	else:
		
		degreeIN = pd.DataFrame.from_dict(G.in_degree(weight='weight'))
		degreeIN.columns=['Nodes','Kin']
		
		degreeOUT = pd.DataFrame.from_dict(G.out_degree(G,weight='weight'))
		degreeOUT.columns=['Nodes','Kout']
		degreeOUT = degreeOUT['Kout']
		
		degNN_IN = pd.DataFrame.from_dict(nx.average_neighbor_degree(G,source='in', target='in',weight='weight'),orient='index',columns=['Kin_nn'])
		degNN_IN = degNN_IN.set_index(np.arange(0,len(degNN_IN),1))
		
		degNN_OUT = pd.DataFrame.from_dict(nx.average_neighbor_degree(G,source='out', target='out',weight='weight'),orient='index',columns=['Kout_nn'])
		degNN_OUT = degNN_OUT.set_index(np.arange(0,len(degNN_OUT),1))
	
		BC = pd.DataFrame.from_dict(nx.betweenness_centrality(G,weight='weight'),orient='index',columns=['BC'])
		BC = BC.set_index(np.arange(0,len(BC),1))
		
		CL = pd.DataFrame.from_dict(nx.closeness_centrality(G),orient='index',columns=['CL'])
		CL = CL.set_index(np.arange(0,len(CL),1))
		
		
		
		df_complete = pd.concat([degreeIN,degreeOUT,degNN_IN,degNN_OUT,BC,CL], axis=1)


		#split dataframe in human and viral part
		virus_index = []
		for j in range(len(df_complete)):
			if df_complete['Nodes'].iloc[j].startswith('9606.')==True:
				virus_index.append(df_complete.index[j])	


		
		human_index = []
		for j in range(len(df_complete)):
			if df_complete['Nodes'].iloc[j].startswith('9606.')!=True:
				human_index.append(df_complete.index[j])


		df_virus = df_complete.drop(virus_index)	
		df_human = df_complete.drop(human_index)
	
	
	
		#array of dataframes
		df = [df_complete,df_human,df_virus]
		
	return df


#%%

#%%	
NamesBC = FileNames('BC_','.png')
NamesCL = FileNames('CL_','.png')
NamesINvsOUT = FileNames('DegreeINvsOUT_','.png')
NamesDegIN = FileNames('DegreeHistIN_','.png')
NamesDegOUT = FileNames('DegreeHistOUT_','.png')
NamesKNN = FileNames('SUB_degNN_deg_','.png')	
	
for i in range(len(G)):
	#PlotBcClvsDegree(G[i],VirusNames[i],21,NamesBC[i],NamesCL[i])
	#PlotINvsOUT(G[i],VirusNames[i],NamesINvsOUT[i])
	PlotDegreeHist(Gsub[i],VirusNames[i],NamesDegIN[i],NamesDegOUT[i])
	#PlotDegreeNNvsDegree(Gsub[i],VirusNames[i],NamesKNN[i])
	


#%%

Gsub = []
for i in range(len(G)):
	Gsub_i = GH.subgraph(nx.nodes(G[i]))
	Gsub.append(Gsub_i)


	
	#%%






def PlotDegreeHist(G,title,nomiIN,nomiOUT):
	
	
	centrality = NetworkCharacterization(G)
	

	if nx.is_directed(G) == False:
		
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax.set_title(title)
		ax.hist(centrality['K']/1000,bins=30)
		ax.set_xlabel('Degree')
		ax.set_ylabel('# Nodes')
		
		
		
		plt.show()
		
	else:
		
		human = centrality[1]
		virus = centrality[2]
		
		fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax1.set_title(title)		
		ax1.hist([human['Kin'],virus['Kin']],bins=30, histtype='barstacked')
		ax1.set_xlabel('Degree IN')
		ax1.set_ylabel('# Nodes')
		ax1.legend(['Homo sapiens','Virus'])
		
		plt.savefig(nomiIN)
		
		fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax2.set_title(title)	
		ax2.hist([human['Kout'],virus['Kout']],bins=30, histtype='barstacked')
		ax2.set_xlabel('Degree OUT')
		ax2.set_ylabel('# Nodes')
		ax2.legend(['Homo sapiens','Virus'])
		
		plt.savefig(nomiOUT)
		
		plt.show()
		
		
	return
	
	
	






def PlotDegreeNNvsDegree(G,title,nomi):

	
	centrality = NetworkCharacterization(G)
	

	if nx.is_directed(G) == False:
		
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax.set_title(title)
		ax.scatter(centrality['K'],centrality['K_nn'],s=20, alpha=0.4, edgecolors='b')
		ax.set_xlabel('K')
		ax.set_ylabel('$K_{NN}$')
		
		plt.show()
	


	
	
		
	else:
		
		human = centrality[1]
		virus = centrality[2]
		
		fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax1.set_title(title)
		ax1.scatter(human['Kin'],human['Kin_nn'],s=30, alpha=0.5, edgecolors='b')
		ax1.scatter(virus['Kin'],virus['Kin_nn'],s=30, alpha=0.5, edgecolors='r')
	
		ax2.scatter(human['Kout'],human['Kout_nn'],s=30, alpha=0.5, edgecolors='b')
		ax2.scatter(virus['Kout'],virus['Kout_nn'],s=30, alpha=0.5, edgecolors='r')
		ax1.legend(title='$K_{NN,IN}$ vs $K_{IN}$')
		ax2.legend(title='$K_{NN,OUT}$ vs $K_{OUT}$')

		plt.savefig(nomi)
		
		plt.show()
		
		
	return
	
	
	
	
	

def PlotINvsOUT(G,title,nomipersalvataggio):
	
	centrality = NetworkCharacterization(G)
	

	if nx.is_directed(G) == False:
		
		print('Error')
		
	
	else:
		
		human = centrality[1]
		virus = centrality[2]
		
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax.set_title(title)
		ax.scatter(human['Kin'],human['Kout'],s=30, alpha=0.5, edgecolors='b')
		ax.scatter(virus['Kin'],virus['Kout'],s=30, alpha=0.5, edgecolors='r')
	
		ax.set_xlabel('Degree IN')
		ax.set_ylabel('Degree OUT')

		plt.savefig(nomipersalvataggio)
		
		plt.show()
		
		
	return
	
	
	







def Mean(array):
	
	mu = sum(array)/len(array)
	
	err = []
	for i in range(len(array)):
		err_i = (array[i]-mu)**2
		err.append(err_i)
	sigma = np.sqrt(sum(err)/len(array))
	
	return mu,sigma

	





def PlotBcClvsDegree(G,title,n,nomipersalvataggioBC,nomipersalvataggioCL): #n = number of the values to average
	

	if nx.is_directed(G) == False:
		
		centrality = NetworkCharacterization(G)
		centrality = centrality.sort_values('K')
		
		
		
		
		deg = []
		bc = []
		c = []

		deg_err = []
		bc_err = []
		c_err = []


		number_intervals = int(len(centrality)/n)
	
		for i in range(number_intervals+1):
			
			
			mu = Mean(centrality['K'].iloc[i*n:(i+1)*n].to_numpy())
			deg_i = mu[0]
			deg_err_i = mu[1]
		

			mu = Mean(centrality['BC'].iloc[i*n:(i+1)*n].to_numpy())
			bc_i = mu[0]
			bc_err_i = mu[1]
		

			mu = Mean(centrality['CL'].iloc[i*n:(i+1)*n].to_numpy())
			c_i = mu[0]
			c_err_i = mu[1]
		



			deg.append(deg_i)
			bc.append(bc_i)
			c.append(c_i)


			deg_err.append(deg_err_i)
			bc_err.append(bc_err_i)
			c_err.append(c_err_i)

		
		
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax.scatter(deg_in,bc,s=30, alpha=0.5, edgecolors='r',label='BC',facecolors='r')
		#ax.scatter(centrality['Kin'],centrality['BC'])
		ax.errorbar(deg_in,bc, yerr=bc_err, fmt="|",color='r')
		ax.set_title(title)
		ax.set_xlabel('Degree')
		ax.set_ylabel('Averaged Betweenness Centrality')
		
		plt.savefig(nomipersalvataggioBC)	
		
	
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax.scatter(deg_in,c,s=30, alpha=0.5, edgecolors='g',label='CL',facecolors='g')
		ax.errorbar(deg_in,c, yerr=c_err, fmt="|",color='g')
		#ax.scatter(centrality['Kin'],centrality['CL'])
		ax.set_title(title)
		ax.set_xlabel('Degree')
		ax.set_ylabel('Averaged Closeness Centrality')
		
		plt.savefig(nomipersalvataggioCL)	
		
		plt.show()
		
		
	else:
		
		

		
		centrality1, centrality2, centrality3 = NetworkCharacterization(G)
		centrality = centrality1.sort_values('Kin')
		
		
		deg_in = []
		bc = []
		c = []

		deg_in_err = []
		bc_err = []
		c_err = []


		number_intervals = int(len(centrality)/n)
	
		for i in range(number_intervals+1):
			
			
			mu = Mean(centrality['Kin'].iloc[i*n:(i+1)*n].to_numpy())
			deg_in_i = mu[0]
			deg_in_err_i = mu[1]
		

			mu = Mean(centrality['BC'].iloc[i*n:(i+1)*n].to_numpy())
			bc_i = mu[0]
			bc_err_i = mu[1]
		

			mu = Mean(centrality['CL'].iloc[i*n:(i+1)*n].to_numpy())
			c_i = mu[0]
			c_err_i = mu[1]
		



			deg_in.append(deg_in_i)
			bc.append(bc_i)
			c.append(c_i)


			deg_in_err.append(deg_in_err_i)
			bc_err.append(bc_err_i)
			c_err.append(c_err_i)

		
		
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax.scatter(deg_in,bc,s=30, alpha=0.5, edgecolors='r',label='BC',facecolors='r')
		#ax.scatter(centrality['Kin'],centrality['BC'])
		ax.errorbar(deg_in,bc, yerr=bc_err, fmt="|",color='r')
		ax.set_title(title)
		ax.set_xlabel('Degree IN')
		ax.set_ylabel('Averaged Betweenness Centrality')
		
		plt.savefig(nomipersalvataggioBC)	
		
	
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
		sns.set_style('whitegrid')
		ax.scatter(deg_in,c,s=30, alpha=0.5, edgecolors='g',label='CL',facecolors='g')
		ax.errorbar(deg_in,c, yerr=c_err, fmt="|",color='g')
		#ax.scatter(centrality['Kin'],centrality['CL'])
		ax.set_title(title)
		ax.set_xlabel('Degree IN')
		ax.set_ylabel('Averaged Closeness Centrality')
		
		plt.savefig(nomipersalvataggioCL)	
		
		plt.show()
	


#%%



import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
import numpy as np


directory = 'C:/Users/silvi/Desktop/Fisica/ComplexNetworks/progetto/ComplexNetworksProject/ComplexNetworksProject'
#directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject'
os.chdir(directory) 

#%%
NameVirusFile = ['string_interactions_WNV.tsv','string_interactions_varicella.tsv',
				 'string_interactions_SARSCov.tsv','string_interactions_parechovirus2.tsv',
				 'string_interactions_mumps.tsv','string_interactions_MARV.tsv',
				 'string_interactions_lassa.tsv','string_interactions_InfluenzaA.tsv',
				 'string_interactions_HTLV-1.tsv','string_interactions_HPV1a.tsv',
				 'string_interactions_HIV1_553.tsv','string_interactions_hepatitisB.tsv',
				 'string_interactions_ebola.tsv','string_interactions_dengue2.tsv',
				 'string_interactions_cytomegalo.tsv','Covid19.txt']
VirusNames=['WNV','Varicella zoster virus','SARS-CoV','Human parechovirus 2',
			'Mumps virus','MARV','Lassa virus','Influenza A','HTLV1','HPV type 1a',
			'HIV1','Hepatitis B','Ebola','Dengue type 2','Cytomegalo','SARS-CoV-2']





#%%


	
###FUNZIONE PER DARE I NOMI AI FILE DEI GRAFICI
#il nome del file diventa: Prefisso + nome virus + formato del file
def FileNames(Prefix, FileType): 
	newnames=[]
	for a in range (len(VirusNames)):
		newnames.append(f"{Prefix}{VirusNames[a]}{FileType}")
	return newnames
	
	





#%%  CREO RETI VIRUS

G = []

for i in range(len(NameVirusFile)):	
	virus = pd.read_csv(NameVirusFile[i], sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	
	graph_virus = nx.DiGraph() 
	for j in range(len(virus)):
		graph_virus.add_edge(virus.iloc[j][0], virus.iloc[j][1], weight=virus.iloc[j][2])
	print(nx.info(graph_virus))
	G.append(graph_virus)







	






	
	
	

#%%

#%% algoritmo di clustering
from networkx.algorithms.community.centrality import girvan_newman	
#from networkx.algorithms.community import asyn_lpa_communities

#%%
alg = []
for i in range(len(G)):
	print(i+1)
	alg.append(list(girvan_newman(G[i])))




#%%




for i in range(len(alg)):
	print('len alg[',i,']: ', len(alg[i]))
	for j in range(len(alg[i])):
		print('    len alg[',i,'][',j,']: ', len(alg[i][j]))
		for k in range(len(alg[i][j])):
			print('        len alg[',i,'][',j,'][',k,']: ', len(alg[i][j][k]))
		
 
	
	
#%%	
	
grafico bc, degree, closeness in funz di nodi
metti a posto programma
scrivi qualcosa
per size netowrk vuoi fare istogramma?
