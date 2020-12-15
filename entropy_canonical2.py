#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 11:32:44 2020

@author: caterina
"""
import numpy as np
import networkx as nx

s = [1, 1, 1, 2, 2, 3];
t = [2, 3, 4, 5, 6, 7];
p = [100, 500, 900, 300, 200, 400]

G = nx.Graph()
#G.add_edge(1,2,weight=100)
#G.add_edge(1,3,weight=500)
#G.add_edge(1,4,weight=900)
#G.add_edge(2,5,weight=300)
#G.add_edge(2,6,weight=200)
#G.add_edge(3,7,weight=400)

G.add_edge(1,2)
G.add_edge(1,3)
G.add_edge(1,4)
G.add_edge(2,5)
G.add_edge(2,6)
G.add_edge(3,7)

#pp = nx.to_pandas_adjacency(G,weight='weight')
pp = nx.to_pandas_adjacency(G)

#pp=nx.to_numpy_matrix(G)


#%%
def entropy_canonical(pp):
	
	precision=10**(-3)
	loops = 5000
	S_0 = 0
		
	n = len(pp) #int = numero di nodi di GH, ovvero numero di righe (e colonne) di pp
	ss=[]
	for i in range(len(pp)):
		ss.append(sum(pp.iloc[i]))
	
	connectivity = ss
	l=len(connectivity)
	avg_conn=sum(connectivity)/n
	z=np.reshape(connectivity/(np.sqrt(n*avg_conn)),(l,1))
	z1=np.reshape(z,(1,l))
	oldz=np.zeros((n,1))
	
	
	limit_single=[]
	for k in range(loops):
		
		U=np.ones((n,1))*z1
		D=np.ones((n,n))+z*z1
		
		UD_nulldiag = U / D
		for i in range(len(UD_nulldiag)):
			UD_nulldiag[i][i] = 0
		z = np.reshape((connectivity / sum(UD_nulldiag.T)), (l,1))
		for i in range(len(z)):
			if z[i][0] < 10**(-15):
				z[i][0] = 10**(-15)
		z1=np.reshape(z,(1,l))
		
		
		for i in range (len(z)):
			if z[i][0]>0:
				if oldz[i][0]==0:
					x=1*(1-z[i][0])/(oldz[i][0]+1)
				#print(z[i][0])
					limit_single.append(abs(x))
				else:
					x=1*(1-z[i][0])/(oldz[i][0])
					limit_single.append(abs(x))
			else: 
				limit_single.append(0)
		
		
		if max(limit_single)<precision:
			break
		oldz=z
		
	
	#Compute link probability
	P = (z*z1) /( (np.ones((n,n)) + z*z1))
	
	for i in range(len(P)):
		P[i][i] = 0
	
	
	#Compute shannon Entropy per node
	import math 	
	
	Primo_termine=0
	Secondo_termine=0
	for i in range(len(P)):
		for j in range (len(P)):
			if i<j:
				Primo_termine=Primo_termine+(P[i][j]*math.log(P[i][j]))
				Secondo_termine=Secondo_termine+((1-P[i][j])*math.log(1-P[i][j]))
	
	S_0=(-Primo_termine-Secondo_termine)/n
	
	
	#Random Shannon entropy
	L=n*avg_conn/2
	p=2*L/(n*(n-1))
	S_r=-(L*math.log(p)+(n*(n-1)/2-L)*math.log(abs(1-p)))/n

	print("S_0")
	print(S_0)
	print("S_r")
	print(S_r)
	print("z")
	print(z)
	print("P")
	print(P)
	#return S_0, S_r, z, P
	return 	

