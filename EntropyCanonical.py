import numpy as np
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt





'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT:
PP - (weighted) undirected adjacency matrix

OUTPUT:
S_0 - canonical entropy of the undirected networks with the same degree distribution
S_r - canonical entropy with fixed GLOBAL link density
z - lagrange multipliers of degree sequence constraints (#node vector)
P - link probability matrix Pij

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''




# WEIGHTS MUST BE LOWER THAN 1



def EntropyCanonicalUndirected(PP):
	
	precision=10**(-3)
	loops = 5000
	S_0 = 0
	
	pp=(PP.values).copy()
	n = len(pp) 
	connectivity = np.reshape(pp.sum(axis=0),(n,1))#row vector, sum of columns terms of pp
	
	
	avg_conn=np.sum(connectivity)/n 
	
	z=connectivity/(np.sqrt(n*avg_conn)) #column vector
	
	oldz=np.zeros((n,1))
	
	
	
	counts = 0
	
	for k in range(loops):
		
		
		U=np.ones((n,1))@z.T
		D=np.ones((n,n))+z@z.T
		
		X=( U/D - np.diag(np.diag(U/D)) ).T
		#z=connectivity ./(sum( ( U./D - diag  ( diag (U./D)) )' ) )'
		z = connectivity /np.reshape(X.sum(axis=0),(n,1))
		
		z=np.maximum(z,10**(-15))
		
		#if     max(abs((z>0).*(1-z./(oldz+(oldz==0)))))<precision
		if  (np.max(abs((z>0)*(1-z/(oldz+(oldz==0)))))<precision):
			break
		
		oldz=z.copy()
		
		counts = counts+1
		
	print('iterations: ',counts)
		
		
	'''Compute link probability'''
	P = (z@z.T) / (np.ones((n,n)) + z@z.T) # nxn matrix
	P=P-np.diag(np.diag(P)) #null diagonal terms
	
	'''Compute Shannon Entropy per node'''
#	Ptot=P.*     log(P+(P==0))+   (ones(n,n)-P)       .*log(   ones(n,n)-P+((   ones(n,n)-P)==0));
	Ptot=P* np.log(P+(P==0))+(np.ones((n,n))-P)  * np.log(np.ones((n,n))-P+((np.ones((n,n))-P)==0))
#	S_0=(-sum(sum(   triu(Ptot,1))))/n
	S_0=(-((np.triu(Ptot,1)).sum(axis=0)).sum(axis=0))/n
	
	
	'''Random Shannon entropy'''
	L=n*avg_conn/2 
	p=2*L/(n*(n-1))  
	S_r=-(L*np.log(p)+(n*(n-1)/2-L)*np.log(abs(1-p)))/n



	print('S_0: ',S_0)
	print('S_r: ',S_r)


	return S_0, S_r, z, P
	









'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT:
PP - (weighted) directed adjacency matrix

OUTPUT:
S_0 - canonical entropy of the undirected networks with the same degree distribution
S_r - canonical entropy with fixed GLOBAL link density
z - lagrange multipliers of degree sequence constraints (#node vector)
P - link probability matrix Pij

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''






def EntropyCanonicalDirected(PP):
	
	precision=10**(-3)
	loops = 5000
	S_0 = 0
	
	pp=(PP.values).copy()
	n = len(pp) 
	connectivity_in = np.reshape(pp.sum(axis=0),(n,1))
	connectivity_out = np.reshape(pp.T.sum(axis=0),(n,1))
	
	
	avg_conn_in=np.sum(connectivity_in)/n 
	avg_conn_out=np.sum(connectivity_out)/n 
	
	
	z_in=connectivity_in/(np.sqrt(n*avg_conn_in)) 
	z_out=connectivity_out/(np.sqrt(n*avg_conn_out))
	
	oldz_in=np.zeros((n,1))
	oldz_out=np.zeros((n,1))
	
	
	
	counts_in = 0
	
	for k in range(loops):
		
		counts_in = counts_in+1
		
		
		U=np.ones((n,1))@z_in.T
		D=np.ones((n,n))+z_in@z_in.T
		
		X=( U/D - np.diag(np.diag(U/D)) ).T
		z = connectivity_in /np.reshape(X.sum(axis=0),(n,1))
		
		z_in=np.maximum(z_in,10**(-15))
		
		
		if  (np.max(abs((z_in>0)*(1-z_in/(oldz_in+(oldz_in==0)))))<precision):
			break
		
		oldz_in=z_in.copy()



	counts_out = 0

	for k in range(loops):
		
		
		
		counts_out = counts_out+1
		
		U=np.ones((n,1))@z_out.T
		D=np.ones((n,n))+z_out@z_out.T
		
		X=( U/D - np.diag(np.diag(U/D)) ).T
		z_out = connectivity_out /np.reshape(X.sum(axis=0),(n,1))
		
		z_out=np.maximum(z_out,10**(-15))
		
		if  (np.max(abs((z_out>0)*(1-z_out/(oldz_out+(oldz_out==0)))))<precision):
			break
		
		oldz_out=z_out.copy()
		

		
		
	'''Compute link probability'''
	P_in = (z_in@z_in.T) / (np.ones((n,n)) + z_in@z_in.T) 
	P_in=P_in-np.diag(np.diag(P_in)) 
	
	P_out = (z_out@z_out.T) / (np.ones((n,n)) + z_out@z_out.T) 
	P_out=P_out-np.diag(np.diag(P_out)) 
	
	P=np.zeros((len(pp),len(pp)))
	for i in range(len(pp)):
		for j in range(i,len(pp)):
			P[i][j]=P_out[i][j]
	for i in range(len(pp)):
		for j in range(i):
			P[i][j]=P_in[i][j]
		
	
	'''Compute Shannon Entropy per node'''
	Ptot=P* np.log(P+(P==0))+(np.ones((n,n))-P)  * np.log(np.ones((n,n))-P+((np.ones((n,n))-P)==0))
	S_0=(-(Ptot.sum(axis=0)).sum(axis=0))/n
	
	
	'''Random Shannon entropy'''
	L=n*(avg_conn_in+avg_conn_out)/2 
	p=L/(n*(n-1)) 
	S_r=-(L*np.log(p)+(n*(n-1)-L)*np.log(abs(1-p)))/n
	
	


	print('S_0: ',S_0)
	print('S_r: ',S_r)
	

	return S_0, S_r, z, P














'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT:
VirusNames - list of virus considered 
GV - list of the corresponding human-virus networks
Gsub - list of the correpsonding human subgraphs 

OUTPUT:
entropy - dataframe with S_0 and S_r values of entropy for each network
plot of entropy differences for each virus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'''





def EntropyDifference(VirusNames,Gsub,GV):

	
		
		
	# entropy of human subgraphs
	S0_sub = []
	Sr_sub = []	
	#Gsub_number_nodes = []

	for i in range(len(Gsub)):
		
		print('GV ',VirusNames[i])
		PP = (nx.to_pandas_adjacency(Gsub[i],weight='weight')/1000)
		[S_0, S_r, z, P] = EntropyCanonicalUndirected(PP)
		S0_sub.append(S_0)
		Sr_sub.append(S_r)
		#Gsub_number_nodes.append(Gsub[i].number_of_nodes())
	
	# entropy of virus networks
	S0_virus = []
	Sr_virus = []	
	#GV_number_nodes = []

	for i in range(len(GV)):
		
		print('Gsub ',VirusNames[i])
		PP = nx.to_pandas_adjacency(GV[i],weight='weight')
		[S_0, S_r, z, P] = EntropyCanonicalDirected(PP)
		S0_virus.append(S_0)
		Sr_virus.append(S_r)
		#GV_number_nodes.append(GV[i].number_of_nodes())
	
	
	
	
	
	# create dataframe with entropy values for each virus
	
	entropy_dict = {'Virus':VirusNames, 'S0 virus':S0_virus, 'Sr virus':Sr_virus, 'S0 sub':S0_sub, 'Sr sub':Sr_sub}
	entropy = pd.DataFrame(data=entropy_dict)
	#entropy.to_csv('Entropy.tsv', sep="\t", index=False) 
	
	
	deltaS0 = entropy['S0 virus']-entropy['S0 sub']
	deltaSr = entropy['Sr virus']-entropy['Sr sub']
	

	# plot of entropy differences
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
	ax.plot(VirusNames, deltaS0 , marker="o", label='$\Delta S_0$', linewidth=0.75)
	ax.plot(VirusNames, deltaSr , marker="o", label='$\Delta S_r$', linewidth=0.75)
	plt.grid(True)
	plt.xticks(rotation=90)
	ax.set_ylabel('$\Delta S$')
	ax.legend(ncol=1 ,loc='best', fontsize=14)
	plt.show()
	
	
	#plot of entropy values GV
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
	ax.plot(VirusNames, entropy['S0 virus'], marker="o", label='$S_0$', linewidth=0.75)
	ax.plot(VirusNames, entropy['Sr virus'], marker="o", label='$S_r$', linewidth=0.75)
	plt.grid(True)
	plt.xticks(rotation=90)
	ax.set_title('Human-Virus PPI graphs')
	ax.set_ylabel('Entropy')
	ax.legend(ncol=1 ,loc='best', fontsize=14)
	plt.show()
	
	
	#plot of entropy values Gsub
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
	ax.plot(VirusNames, entropy['S0 sub'], marker="o", label='$S_0$', linewidth=0.75)
	ax.plot(VirusNames, entropy['Sr sub'], marker="o", label='$S_r$', linewidth=0.75)
	plt.grid(True)
	plt.xticks(rotation=90)
	ax.set_title('Human PPI subgraphs')
	ax.set_ylabel('Entropy')
	ax.legend(ncol=1 ,loc='best', fontsize=14)
	plt.show()





	return entropy

	

	

