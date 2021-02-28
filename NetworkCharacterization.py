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





#%% CREO ARRAY DI DATAFRAMES CON INFO CENTRALITY


#array di dataframe misure di centrality
centrality = []
#array di dataframe solo indici=proteine virus
virus = []
#array di dataframe solo indici=proteine umane
human = []

for i in range(len(G)):
	
	degreeIN = pd.DataFrame.from_dict(nx.in_degree_centrality(G[i]),orient='index',columns=['IN'])
	degreeOUT = pd.DataFrame.from_dict(nx.out_degree_centrality(G[i]),orient='index', columns=['OUT'])
	BC=pd.DataFrame.from_dict(nx.betweenness_centrality(G[i],weight='weight'),orient='index',columns=['BC'])
	CC=pd.DataFrame.from_dict(nx.clustering(G[i],weight='weight'),orient='index',columns=['CC'])
	
	clos = pd.DataFrame.from_dict(nx.closeness_centrality(G[i]),orient='index',columns=['clos'])
	neigh_deg_in = pd.DataFrame.from_dict(nx.average_neighbor_degree(G[i],source='in', target='in',weight='weight'),orient='index',columns=['avg IN'])
	neigh_deg_out = pd.DataFrame.from_dict(nx.average_neighbor_degree(G[i],source='out', target='out',weight='weight'),orient='index',columns=['avg OUT'])
	
	
	df = pd.concat([degreeIN,degreeOUT,BC,CC,clos,neigh_deg_in,neigh_deg_out], axis=1)
	
	centrality.append(df)
	
	
	#separa dataframe
	virus_index = []
	for j in range(len(df)):
		if df.index[j].startswith('9606.')==True:
			virus_index.append(df.index[j])
		
	human_index = []
	for j in range(len(df)):
		if df.index[j].startswith('9606.')!=True:
			human_index.append(df.index[j])


	df_virus = df.drop(virus_index)	
	df_human = df.drop(human_index)
	
	virus.append(df_virus)
	human.append(df_human)




#%% algoritmo di clustering
from networkx.algorithms.community.centrality import girvan_newman	
from networkx.algorithms.community import asyn_lpa_communities

#%%
alg = []
for i in range(len(G)):
	print(i+1)
	alg.append(asyn_lpa_communities(G[i]))
	
	
#%%	
	
grafico bc, degree, closeness in funz di nodi
metti a posto programma
scrivi qualcosa
per size netowrk vuoi fare istogramma?
	
#%%

#BETWENNESS HISTOGRAM

NamesBC=FileNames('HistoBC_','.png')


for i in range(len(G)):
	bc=pd.DataFrame.from_dict(nx.betweenness_centrality(G[i],weight='weight'),orient='index',columns=['BC'])
	
	BCmax=max(bc['BC'])
	
	BCvirus=[]#normalized
	BChuman=[]#normalized
	
	for j in range (len(bc)):
		if bc.index[j].startswith('9606.')==True:
			BChuman.append(bc.iloc[j][0]/BCmax)
		else:
			BCvirus.append(bc.iloc[j][0]/BCmax)
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
	ax.hist([BChuman,BCvirus], bins=30, log =False,color=['dodgerblue','orangered'], histtype='barstacked') #,stacked=False) #HISTOGRAM
	
	#ax.set_title(f"BC {VirusNames[i]}")
	ax.set_title(VirusNames[i])
	ax.set_xlabel('BC/BCmax',labelpad=2)
	#ax.set_ylabel('Entries')
	ax.legend(['human','virus'])
	plt.grid(True)
	#plt.show()
	
	plt.savefig(NamesBC[i])


#%% CLOSENESS HISTOGRAM
NamesC=FileNames('HistoCloseness_','.png')

for i in range(len(G)):		
	c=pd.DataFrame.from_dict(nx.closeness_centrality(G[i]),orient='index',columns=['C'])
	Cmax=max(c['C'])
	
	Cvirus=[]#normalized
	Chuman=[]#normalized
	for j in range (len(c)):
		if c.index[j].startswith('9606.')==True:
			Chuman.append(c.iloc[j][0]/Cmax)
		else:
			Cvirus.append(c.iloc[j][0]/Cmax)

	
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
	ax.hist([Chuman,Cvirus], bins=30, log =False,color=['dodgerblue','orangered'], histtype='barstacked') #HISTOGRAM
	
	ax.set_title(VirusNames[i])
	ax.set_xlabel('C/Cmax',labelpad=2)
	#ax.set_ylabel('Entries')
	ax.legend(["human","virus"])
	plt.grid(True)
	plt.savefig(NamesC[i])
#%% ISTOGRAMMA CLUSTERING COEFFICIENT
	
	
#directory dove si salveranno i grafici	
directory= '/home/caterina/Documenti/GitHub/ComplexNetworksProject/grafici'
os.chdir(directory) 

##nome dei file
NamesCC=FileNames('HistoClusteringCoeff_','.png')

for i in range(len(G)):		
	#clustering coefficient
	cc=pd.DataFrame.from_dict(nx.clustering(G[i],weight='weight'),orient='index',columns=['CC'])
	CCmax=max(cc['CC'])
	
	CCvirus=[]#normalized
	CChuman=[]#normalized
	for j in range (len(cc)):
		if cc.index[j].startswith('9606.')==True:
			CChuman.append(cc.iloc[j][0]/CCmax)
		else:
			CCvirus.append(cc.iloc[j][0]/CCmax)

	
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
	ax.hist([CChuman,CCvirus], bins=30) #HISTOGRAM
	
	ax.set_title(VirusNames[i])
	ax.set_xlabel('CC/CCmax')
	ax.set_ylabel('Entries')
	ax.legend(["human","virus"])
	plt.savefig(NamesCC[i])
	
	
	
	




	
	
#%%%	DEGREE IN vs DEGREE OUT



NamesDegree = FileNames('DegreeINvsOUT_','.png')
	
for i in range(len(G)):
	
	
	#df = centrality[i]#.sort_values('IN')	
	h = human[i]
	v = virus[i]
	
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
	#sns.set_style('whitegrid')
	ax.set_title(VirusNames[i])
	ax.scatter(h['IN'],h['OUT'],s=30, alpha=0.5, edgecolors='b')
	ax.scatter(v['IN'],v['OUT'],s=30, alpha=0.5, edgecolors='r')
	
	ax.set_xlabel('Degree IN')
	ax.set_ylabel('Degree OUT')


	#plt.savefig(NamesDegree[i])	
	
	
	
#%%	  HISTOGRAMMA DEGREE 

NamesDegreeHist = FileNames('DegreeHistLOG_','.png')
	
for i in range(len(G)):
	
	
	df = centrality[i]#.sort_values('IN')	
	#h = human[i]
	#v = virus[i]
	
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
	sns.set_style('whitegrid')
	ax.set_title(VirusNames[i])
	
	#ax.hist([h['IN'],h['OUT'],v['IN'],v['OUT']],bins=10, color=['blue','dodgerblue','orangered','orange'], log =True)
	ax.hist([df['IN'],df['OUT']],bins=30)
	
	#ax.legend(['Human IN', 'Human OUT','Virus IN', 'Virus OUT'])
	ax.legend(['IN','OUT'])
	
	ax.set_xlabel('Degree')
	ax.set_ylabel('Entries')


	#plt.savefig(NamesDegreeHist[i])	
	




#%%       degree neighbours vs degree



NamesDegree = FileNames('degNN_deg_','.png')
	
for i in range(len(G)):
	
	
	#df = centrality[i]#.sort_values('IN')	
	h = human[i]
	v = virus[i]
	
	
	fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, figsize=(4,3))
	sns.set_style('whitegrid')
	ax1.set_title(VirusNames[i])
	ax1.scatter(h['IN'],h['avg IN'],s=30, alpha=0.5, edgecolors='b')
	ax1.scatter(v['IN'],v['avg IN'],s=30, alpha=0.5, edgecolors='r')
	
	ax2.scatter(h['OUT'],h['avg OUT'],s=30, alpha=0.5, edgecolors='b')
	ax2.scatter(v['OUT'],v['avg OUT'],s=30, alpha=0.5, edgecolors='r')
	#ax.hist([],bins=30)
	ax1.legend(title='$K_{IN}$ vs $K_{NN,IN}$')
	ax2.legend(title='$K_{OUT}$ vs $K_{NN,OUT}$')
	
	#ax1.set_xlabel('$K_{IN}$')
	#ax1.set_ylabel('$K_{NN,IN}$')
	
	#ax2.set_xlabel('$K_{OUT}$')
	#ax2.set_ylabel('$K_{NN,OUT}$')


	plt.savefig(NamesDegree[i])	
	




#%%  centrality mediate


def media(array):
	
	mu = sum(array)/len(array)
	
	scarti = []
	for i in range(len(array)):
		scartoi = (array[i]-mu)**2
		scarti.append(scartoi)
	err = sum(scarti)/len(scarti)
	
	return mu,err

#%%

for j in range(len(G)):  #  problema= 6 lassa virus
	#j=6
	
	df = centrality[j].sort_values('IN')
	
	x = []
	deg_in = []
	deg_in_prova = []
	bc = []
	c = []
	cc = []
	
	x_s = []
	deg_in_s = []
	bc_s = []
	c_s = []
	cc_s = []
	
	assex = np.arange(1,len(df)+1,1)
	#num_divisioni = 10
	l = 21#int(len(assex)/num_divisioni)
	num_divisioni = int(len(assex)/l)
	
	for i in range(num_divisioni+1):
		
		#media e standard deviation
		mu_x = assex[i*l:(i+1)*l].mean()
		s_x = assex[i*l:(i+1)*l].std()
		
		
		mu = media(df['IN'].iloc[i*l:(i+1)*l].to_numpy())
		mu_deg_in = mu[0]
		s_deg_in = mu[1]
		
		#mu_deg_in_prova = sum(df['IN'].iloc[i*l:(i+1)*l])/l

		
		bc_norm = df['BC']/max(df['BC'])
		mu = media(bc_norm.iloc[i*l:(i+1)*l].to_numpy())
		mu_bc = mu[0]
		s_bc = mu[1]
		
		
		c_norm = df['clos']/max(df['clos'])
		mu_c = c_norm.iloc[i*l:(i+1)*l].to_numpy().mean()
		s_c = c_norm.iloc[i*l:(i+1)*l].to_numpy().std()
		
		cc_norm = df['CC']/max(df['CC'])
		mu_cc = cc_norm.iloc[i*l:(i+1)*l].to_numpy().mean()
		s_cc = cc_norm.iloc[i*l:(i+1)*l].to_numpy().std()
		
		

		x.append(mu_x)
		deg_in.append(mu_deg_in)
		deg_in_prova.append(mu_deg_in_prova)
		#print(deg_in)
		#print(deg_in_prova)
		bc.append(mu_bc)
		c.append(mu_c)
		cc.append(mu_cc)
		
		x_s.append(s_x)
		deg_in_s.append(s_deg_in)
		bc_s.append(s_bc)
		c_s.append(s_c)
		cc_s.append(s_cc)
		
		
		
		
		
	
		
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
	sns.set_style('whitegrid')

	
	#ax.plot(deg_in,cc)


	#ax.plot(x,deg_in,color='b',label='IN')
	#ax.plot(x,bc,color='r',label='BC')
	#ax.plot(x,c,color='g',label='clos')
	#ax.plot(x,cc,label='CC')
	#ax.scatter(x,deg_in,color='b',s=30, alpha=0.5, edgecolors='b')
	#ax.scatter(x,bc,color='r',s=30, alpha=0.5, edgecolors='r')
	#ax.scatter(x,c,color='r',s=30, alpha=0.5, edgecolors='g')
	
	ax.plot(x,deg_in,color='b',label='IN')
	ax.plot(x,bc,color='r',label='BC')
	ax.plot(x,c,color='g',label='clos')
	#ax.scatter(x,deg_in,color='b',s=30, alpha=0.5, edgecolors='b')
	#ax.scatter(x,bc,color='r',s=30, alpha=0.5, edgecolors='r')
	#ax.scatter(x,c,color='r',s=30, alpha=0.5, edgecolors='g')
	ax.errorbar(x,deg_in, xerr=deg_in_s, fmt="|")
	ax.errorbar(x,bc, xerr=bc_s, fmt="|")
	ax.errorbar(x,c, xerr=c_s, fmt="|")
	
	
	ax.set_title(VirusNames[j])
	ax.set_xlabel('Nodes')
	ax.set_ylabel('Averaged centrality measures')
	ax.legend(ncol=1 ,loc='best', fontsize=10)
	
	
	
	
	
	
#%%	
	
	
	#asse x
	assex = np.arange(1,len(df)+1,1)
	x=[]			
	l_x = len(assex)/5
	l_bc = (max(df['BC'])-min(df['BC']))/(max(df['BC'])*5)

	bc=[]
	
	for i in range(5):
		x_single=[]
		bc_single=[]
		for j in range(len(df)):
			if assex[j] < (i+1)*l_x:
				if assex[j] > i*l_x:
					x_single.append(assex[j])
		x.append(sum(x_single)/len(x_single))
		for j in range(len(df)):
			if df['BC'][j] < (i+1)*l_bc:
				if df['BC'][j] > i*l_bc:
					bc_single.append(df['BC'][j])
		bc.append(sum(bc_single)/len(bc_single))
	
	
df.to_numpy().mean()
std()
describe()

		
						
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
	sns.set_style('whitegrid')
	
	ax.plot(x,bc)
	
		
	#ax.scatter(assex,df['BC']/max(df['BC']),s=30, alpha=0.5, edgecolors='b')
	#ax.scatter(assex,df['clos']/max(df['clos']),s=30, alpha=0.5, edgecolors='r')
	#ax.scatter(assex,df['IN'],s=30, alpha=0.5, edgecolors='b')
	#ax.scatter(assex,df['OUT'],s=30, alpha=0.5, edgecolors='g')
	


	

#%%
#	cmap=plt.cm.BuGn_r),,s=30, alpha=0.5, edgecolors='b'facecolors='none', edgecolors='b')#
	#ax.hist2d(df['IN'],df['OUT'],bins=50,cmin = 1, cmap=plt.cm.jet)#plt.cm.Reds)
	
	#ax.hexbin(df['IN'],df['OUT'], gridsize=50, cmap=plt.cm.BuGn_r)

	

	#assex = np.arange(1,len(df)+1,1)
	#ax.plot(assex,df['IN'],label='degree IN')
	#ax.plot(assex,df['OUT'],label='degree OUT')
	#ax.scatter(assex,df['BC']/max(df['BC']),label='betweenness',linewidths=0.001,color='r')
	#ax.plot(assex,df['CC'],label='clustering coeff')
	#ax.scatter(assex,df['clos'],label='closeness')
	
	
	#ax.scatter(df['IN'],df['BC']/max(df['BC']),linewidths=0.001)
	
	#ax.legend(ncol=1 ,loc='best', fontsize=12)


	

#%% divido dataframe degree IN e OUT + istogramma degree IN e OUT
#    + istogramma degree in e out
	



for i in range(len(G)):
	
	df = centrality[i].sort_values('IN')
	
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
	
	#ax.scatter(df['IN'],df['BC'],s=30, alpha=0.5, edgecolors='b')
	#ax.scatter(df['IN'],df['clos'],s=30, alpha=0.5, edgecolors='r')
	
	#ax.scatter(df['OUT'],df['BC'],s=30, alpha=0.5, edgecolors='b')
	#ax.scatter(df['OUT'],df['clos'],s=30, alpha=0.5, edgecolors='r')
	
	assex = np.arange(1,len(df)+1,1)
	ax.scatter(assex,df['BC']/max(df['BC']),s=30, alpha=0.5, edgecolors='b')
	ax.scatter(assex,df['clos']/max(df['clos']),s=30, alpha=0.5, edgecolors='r')
	#ax.scatter(assex,df['IN'],s=30, alpha=0.5, edgecolors='b')
	ax.scatter(assex,df['OUT'],s=30, alpha=0.5, edgecolors='g')
	
	


#istogramma ============= rwidth=0.5

