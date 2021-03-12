import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns








'''
##### INPUT:
 G = directed or undirected network

##### OUTPUT:
 if G undirected: 
 	df = dataframe with degree, degree of nearest neighbors, betweenness 
        centrality and closeness centrality as columns and each row corresponding 
        to a node of the network	
 if G directed:
	df = array of three dataframes as entries: 
	df[0] = dataframe with degree, In degree, Out degree, average degree of nearest neighbors, 
	        (average In degree of nearest neighbors, average Out degree of nearest 
		    neighbors), betweenness centrality and closeness 
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
		
		degNN = pd.DataFrame.from_dict(nx.average_neighbor_degree(G, weight='weight'),orient='index',columns=['K_nn'])
		degNN = degNN.set_index(np.arange(0,len(degNN),1))
		
		BC = pd.DataFrame.from_dict(nx.betweenness_centrality(G,weight='weight'),orient='index',columns=['BC']) #è già normalizzata
		BC = BC.set_index(np.arange(0,len(BC),1))
		
		CL = pd.DataFrame.from_dict(nx.closeness_centrality(G),orient='index',columns=['CL']) #è già normalizzata
		CL = CL.set_index(np.arange(0,len(CL),1))
		
	
		#dataframe
		df = pd.concat([degree,degNN,BC,CL], axis=1)

		
		
	else:
		
		degree = pd.DataFrame.from_dict(G.degree(weight='weight')) 
		degree.columns=['Nodes','K']
		
		degreeIN = pd.DataFrame.from_dict(G.in_degree(weight='weight'))
		degreeIN.columns=['Nodes','Kin']
		degreeIN = degreeIN['Kin']
		
		degreeOUT = pd.DataFrame.from_dict(G.out_degree(G,weight='weight'))
		degreeOUT.columns=['Nodes','Kout']
		degreeOUT = degreeOUT['Kout']
		
		degNN = pd.DataFrame.from_dict(nx.average_neighbor_degree(G, weight='weight'),orient='index',columns=['K_nn'])
		degNN = degNN.set_index(np.arange(0,len(degNN),1))
		
		#degNN_IN = pd.DataFrame.from_dict(nx.average_neighbor_degree(G,source='in', target='in',weight='weight'),orient='index',columns=['Kin_nn'])
		#degNN_IN = degNN_IN.set_index(np.arange(0,len(degNN_IN),1))
		
		#degNN_OUT = pd.DataFrame.from_dict(nx.average_neighbor_degree(G,source='out', target='out',weight='weight'),orient='index',columns=['Kout_nn'])
		#degNN_OUT = degNN_OUT.set_index(np.arange(0,len(degNN_OUT),1))
	
		BC = pd.DataFrame.from_dict(nx.betweenness_centrality(G,weight='weight'),orient='index',columns=['BC'])
		BC = BC.set_index(np.arange(0,len(BC),1))
		
		CL = pd.DataFrame.from_dict(nx.closeness_centrality(G),orient='index',columns=['CL'])
		CL = CL.set_index(np.arange(0,len(CL),1))
		
		
		
		df_complete = pd.concat([degree,degreeIN,degreeOUT,degNN,#degNN_IN,degNN_OUT,
						                                       BC,CL], axis=1)


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





'''
INPUT:
	centrality = NetworkCharacterization(G)
	title of the plot
OUTPUT:
	plot of histogram degree (distinction in two plots corresponding to IN 
						        and OUT degree in case of direct network)
'''



 
def PlotDegreeHist(centrality,title):#,nomiIN,nomiOUT):
	
	
	# undirected network
	if type(centrality) == pd.core.frame.DataFrame:
		
		
		sns.set_style('whitegrid')
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax.set_title(title)
		ax.hist(centrality['K']/1000,bins=30)#5000
		ax.set_xlabel('Degree')
		ax.set_ylabel('# Nodes')
	
		
		plt.show()
		
		
	# directed network	
	else:
		
		human = centrality[1]
		virus = centrality[2]
		
		sns.set_style('whitegrid')
		fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax1.set_title(title)		
		ax1.hist([human['Kin'],virus['Kin']],bins=30, histtype='barstacked')
		ax1.set_xlabel('Degree IN')
		ax1.set_ylabel('# Nodes')
		ax1.legend(['Homo sapiens','Virus'])
		
		#plt.savefig(nomiIN)
		
		sns.set_style('whitegrid')
		fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax2.set_title(title)	
		ax2.hist([human['Kout'],virus['Kout']],bins=30, histtype='barstacked')
		ax2.set_xlabel('Degree OUT')
		ax2.set_ylabel('# Nodes')
		ax2.legend(['Homo sapiens','Virus'])
		
		#plt.savefig(nomiOUT)
		
		plt.show()
		
		
	return
	
	
	







'''
INPUT:
	centrality = NetworkCharacterization(G)
	title of the plot
OUTPUT:
	plot of average neighbor degree vs degree (distinction in two plots 
	corresponding to IN and OUT degree in case of direct network)  
'''



def PlotDegreeNNvsDegree(centrality,title):#,nomi):


	
	# undirected network
	if type(centrality) == pd.core.frame.DataFrame:
		
		sns.set_style('whitegrid')
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax.set_title(title)
		ax.scatter(centrality['K']/1000,centrality['K_nn'],s=20, alpha=0.4, edgecolors='b')
		ax.set_xlabel('K')
		ax.set_ylabel('$K_{NN}$')
		
		plt.show()
	


	
	
	# directed network	
	else:
		
		human = centrality[1]
		virus = centrality[2]
		
		
		
		
		sns.set_style('whitegrid')
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax.set_title(title)
		ax.scatter(human['K'],human['K_nn'],s=30, alpha=0.5, edgecolors='b', label='Homo Sapiens')
		ax.scatter(virus['K'],virus['K_nn'],s=30, alpha=0.5, edgecolors='r', label='Virus')
		ax.set_xlabel('K')
		ax.set_ylabel('$K_{NN}$')
		ax.legend()
		
		plt.show()
		
		
		
		#sns.set_style('whitegrid')
		#fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		#ax.set_title(title)
		#ax.scatter(human['Kin'],human['Kin_nn'],s=30, alpha=0.5, edgecolors='b', label='Homo Sapiens')
		#ax.scatter(virus['Kin'],virus['Kin_nn'],s=30, alpha=0.5, edgecolors='r', label='Virus')
		#ax.set_xlabel('$K_{IN}$')
		#ax.set_ylabel('$K_{NN,IN}$')
		#ax.legend()
		
		#plt.show()
	
		#sns.set_style('whitegrid')
		#fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		#ax.set_title(title)
		#ax.scatter(human['Kout'],human['Kout_nn'],s=30, alpha=0.5, edgecolors='b', label='Homo Sapiens')
		#ax.scatter(virus['Kout'],virus['Kout_nn'],s=30, alpha=0.5, edgecolors='r', label='Virus')
		#ax.set_xlabel('$K_{OUT}$')
		#ax.set_ylabel('$K_{NN,OUT}$')
		#ax.legend()

		#plt.savefig(nomi)
		
		#plt.show()
		
		
	return
	
	
	
	
	
	
	
	
	
'''
INPUT:
	centrality = NetworkCharacterization(G)
	title of the plot
OUTPUT:
	plot of IN degree vs OUT degree (Error if G undirected)
'''

	
	

def PlotINvsOUT(centrality,title):#,nomipersalvataggio):
	
	
	# undirected network
	if type(centrality) == pd.core.frame.DataFrame:
		
		print('Error')
		
	
	# directed network
	else:
		
		human = centrality[1]
		virus = centrality[2]
		
		sns.set_style('whitegrid')
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax.set_title(title)
		ax.scatter(human['Kin'],human['Kout'],s=30, alpha=0.5, edgecolors='b', label='Homo Sapiens')
		ax.scatter(virus['Kin'],virus['Kout'],s=30, alpha=0.5, edgecolors='r', label='Virus')
	
		ax.set_xlabel('Degree IN')
		ax.set_ylabel('Degree OUT')
		ax.legend()

		#plt.savefig(nomipersalvataggio)
		
		plt.show()
		
		
	return
	
	
	









	
'''
INPUT:
	array of numbers
OUTPUT:
	mu = average 
	sigma = root mean square deviation
'''




def Mean(array):
	
	mu = sum(array)/len(array)
	
	err = []
	for i in range(len(array)):
		err_i = (array[i]-mu)**2
		err.append(err_i)
	sigma = np.sqrt(sum(err)/len(array))
	
	return mu,sigma

	













	
'''
INPUT:
	centrality = NetworkCharacterization(G)
	title = title of the plot
	n = number of nodes used to calculate the average of the correpsonding centrality measures
OUTPUT:
	plot of BC + averaged BC vs degree
	plot of CL + averaged CL vs degree
'''



def PlotBcClvsDegree(centrality,title,n):#,nomipersalvataggioBC,nomipersalvataggioCL): #n = number of the values to average
	

	# undirected network
	if type(centrality) == pd.core.frame.DataFrame:
		
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
			deg_i = mu[0]/1000
			deg_err_i = mu[1]/1000
		

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

		
		
		sns.set_style('whitegrid')
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax.scatter(deg,bc,s=30 ,edgecolors='r',label='averaged BC',facecolors='r')
		ax.scatter(centrality['K']/1000,centrality['BC'],s=30, alpha=0.3, edgecolors='b',label='BC')
		ax.errorbar(deg,bc, yerr=bc_err, fmt="|",color='r')
		ax.set_title(title)
		ax.set_xlabel('Degree')
		ax.set_ylabel('Betweenness Centrality')
		ax.legend()
		#plt.savefig(nomipersalvataggioBC)	
		
	
		sns.set_style('whitegrid')
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax.scatter(deg,c,s=30, edgecolors='g',label='averaged CL',facecolors='g')
		ax.errorbar(deg,c, yerr=c_err, fmt="|",color='g')
		ax.scatter(centrality['K']/1000,centrality['CL'],s=30, alpha=0.3, edgecolors='b',label='CL')
		ax.set_title(title)
		ax.set_xlabel('Degree')
		ax.set_ylabel('Closeness Centrality')
		ax.legend()
		#plt.savefig(nomipersalvataggioCL)	
		
		plt.show()
		
	
		
	# directed network
	else:
		
		

		
		#centrality1, centrality2, centrality3 = NetworkCharacterization(G)
		centrality = centrality[0].sort_values('Kin')
		
		
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

		
		
		sns.set_style('whitegrid')
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax.scatter(deg_in,bc,s=30 ,edgecolors='r',label='averaged BC',facecolors='r')
		ax.scatter(centrality['Kin'],centrality['BC'],s=30, alpha=0.3, edgecolors='b',label='BC')
		ax.errorbar(deg_in,bc, yerr=bc_err, fmt="|",color='r')
		ax.set_title(title)
		ax.set_xlabel('Degree IN')
		ax.set_ylabel('Betweenness Centrality')
		ax.legend()
		
		#plt.savefig(nomipersalvataggioBC)	
		
	
		sns.set_style('whitegrid')
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
		ax.scatter(deg_in,c,s=30, edgecolors='g',label='averaged CL',facecolors='g')
		ax.errorbar(deg_in,c, yerr=c_err, fmt="|",color='g')
		ax.scatter(centrality['Kin'],centrality['CL'],s=30, alpha=0.3, edgecolors='b',label='CL')
		ax.set_title(title)
		ax.set_xlabel('Degree IN')
		ax.set_ylabel('Closeness Centrality')
		ax.legend()
		
		#plt.savefig(nomipersalvataggioCL)	
		
		plt.show()
		
		
	
	
	
	return
	
