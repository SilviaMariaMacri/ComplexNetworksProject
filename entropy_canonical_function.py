'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% pp - (weighted) undirected adjacency matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
% S_0 - canonical entropy of the undirected networks with the same degree distribution
% S_r - canonical entropy with fixed GLOBAL link density
% z - lagrange multipliers of degree sequence constraints (#node vector)
% P - link probability matrix Pij
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'''
# WEIGHTS MUST BE LOWER THAN 1

import numpy as np

def entropy_canonical_c(PP):
	
	precision=10**(-3)
	loops = 5000
	S_0 = 0
	
	pp=(PP.values).copy()
	n = len(pp) 
	connectivity = np.reshape(pp.sum(axis=0),(n,1))#row vector, sum of columns terms of pp
	
	avg_conn=np.sum(connectivity)/n
	
	z=connectivity/(np.sqrt(n*avg_conn)) #column vector
	
	oldz=np.zeros((n,1))
	
	
	for k in range(loops):
		
		print('iterazione: ', k+1)
		
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
		
		#print (z)
		
		
	'''Compute link probability'''
	P = (z@z.T) / (np.ones((n,n)) + z@z.T) # nxn matrix
	P=P-np.diag(np.diag(P)) #null diagonal terms
	
	'''Compute Shannon Entropy per node'''
#	Ptot=P.*     log(P+(P==0))+   (ones(n,n)-P)       .*log(   ones(n,n)-P+((   ones(n,n)-P)==0));
	Ptot=P* np.log(P+(P==0))+(np.ones((n,n))-P)  * np.log(np.ones((n,n))-P+((np.ones((n,n))-P)==0))
#	S_0=(-sum(sum(   triu(Ptot,1))))/n
	S_0=(-((np.triu(Ptot,1)).sum(axis=0)).sum(axis=0))/n
	
	
	'''Random Shannon entropy'''
	L=n*avg_conn/2 #scalar
	p=2*L/(n*(n-1))
	S_r=-(L*np.log(p)+(n*(n-1)/2-L)*np.log(abs(1-p)))/n

	return S_0, S_r, z, P
	

