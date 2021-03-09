
import numpy as np




def entropy_canonical_s(pp):
	
	precision = 10**(-3)
	loops = 5000
	S_0 = 0
	
	
	
	n = len(pp) #int = numero di nodi di GH, ovvero numero di righe (e colonne) di pp
	connectivity = pp.sum(axis=1) #dataframe = somma dei pesi dei link per ogni nodo
	connectivity = np.array(connectivity) #convert to array
	l=len(connectivity)
	
	avg_conn = sum(connectivity)/n #float = valore medio del peso per nodo
	z = connectivity/(np.sqrt(n*avg_conn)) #indica il peso dei link legati a ogni nodo
	oldz = np.zeros(n) # è un array!!!!
	z1 = np.reshape(z,(l,1)) #=z'

	#count = 0
	for k in range(loops):
		
		print('iterazione: ',k+1)
		
		  
		U = np.ones(n)*z1 #valori uguali sono nella stessa riga (diverso da matlab)
		D = np.ones((n,n)) + z*z1

		UD_nulldiag = U.T / D
		for i in range(len(UD_nulldiag)):
			UD_nulldiag[i][i] = 0
		z = connectivity / sum(UD_nulldiag.T)
	
		for i in range(len(z)):
			if z[i] < 10**(-15):
				z[i] = 10**(-15)
		
		
		#prova, in realtà funzionerebbe anche come scritto in matlab
		limit_single = np.empty(n)
		for i in range(len(z)):
			if z[i] > 0:
				if oldz[i] == 0:
					limit_single[i] = abs(1-z[i])
				else:
					limit_single[i] = abs(1-z[i]/oldz[i])
			else:
				limit_single[i] = 0
				
		limit = max(limit_single)
			
		#limit = max(abs((z>0)*(1-z/(oldz+(oldz==0)))))
		if limit < precision:
			break
				
		oldz = z
		z1 = np.reshape(z,(l,1))
		print(z)
		#count = count+1
	
	
	#Compute link probability
	#print(z)
	P = (z*z1) / (np.ones((n,n)) + z*z1)
	for i in range(len(P)):
		P[i][i] = 0
	#print(P)
	
	
	
	#Compute Shannon entropy per node
	Ptot=P*np.log(P+(P==0))+(np.ones((n,n))-P)*np.log(np.ones((n,n))-P +((np.ones((n,n))-P)==0))
	element = 0
	for i in range(len(Ptot)):
		for j in range (len(Ptot)):
			if i<j:
				element = element + Ptot[i][j]
	S_0 = -element/n		
	
	
	
	
	#Compute Random Shannon entropy per node - teorico
	L=n*avg_conn/2
	p=2*L/(n*(n-1))
	S_r=-(L*np.log(p)+(n*(n-1)/2 - L)*np.log(1-p))/n
	
		
		
	return S_0, S_r, z, P







def entropy_canonical_c(pp):
	
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
	 #perchè qui?
	for k in range(loops):
		
		limit_single=[] 
		print('iterazione: ', k+1)
		
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
		
		print(z)
	
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

