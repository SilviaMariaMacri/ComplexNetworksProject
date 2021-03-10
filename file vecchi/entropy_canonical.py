'''##################################
cd C:\Users\silvi\Desktop\Fisica\ComplexNetworks\05_NetworkProperties\PropertiesLab2

U = pd.read_csv('US_largest500_airportnetwork.txt', sep="  ",names = ['partenza','arrivo','vuota','score'], header=None)
colonne = ['partenza','arrivo','score']
U = U[colonne]
#len(U) = 2980
#nx.info(G)
#'Name: \nType: Graph\nNumber of nodes: 500\nNumber of edges: 2980\nAverage degree:  11.9200'
'''
##################################
### NETWORK DI PROVA###
#######################

#MATLAB
s = [1 1 1 2 2 3];
t = [2 3 4 5 6 7];
p = [100. 500. 900. 300. 200. 400.]
G = graph(s,t,p)

pp = adjacency(G,'weighted')


#PYTHON

s = [1, 1, 1, 2, 2, 3];
t = [2, 3, 4, 5, 6, 7];
p = [100, 500, 900, 300, 200, 400]

G = nx.Graph()
G.add_edge(1,2,weight=100)
G.add_edge(1,3,weight=500)
G.add_edge(1,4,weight=900)
G.add_edge(2,5,weight=300)
G.add_edge(2,6,weight=200)
G.add_edge(3,7,weight=400)

pp = nx.to_pandas_adjacency(G,weight='weight')



# nx.draw(G,with_labels=True)


	
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

#dovrebbe andare bene:
	# pp = nx.to_pandas_adjacency(GH,weight='weight')
	# pp è un dataframe


import numpy as np


'''
1. perchè z definita tramite rapporto con radice quadrata?
   (tipo connettività normalizzata?)
'''


def entropy_canonical(pp):
	
	precision = 10**(-3)
	loops = 5000
	S_0 = 0
	
	
	
	n = len(pp) #int = numero di nodi di GH, ovvero numero di righe (e colonne) di pp
	connectivity = pp.sum(axis=1) #dataframe = somma dei pesi dei link per ogni nodo
	connectivity = np.array(connectivity) #convert to array
	
	avg_conn = sum(connectivity)/n #float = valore medio del peso per nodo
	z = connectivity/(np.sqrt(n*avg_conn)) #indica il peso dei link legati a ogni nodo
	oldz = np.zeros(n) # è un array!!!!
	z1 = np.reshape(z,(7,1)) #=z'


	count = 0
	for k in range(loops):
		
		  
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
		z1 = np.reshape(z,(7,1))
		
		count = count+1
	
	
	#Compute link probability
	#print(z)
	P = (z*z1) / (np.ones((n,n)) + z*z1)
	for i in range(len(P)):
		P[i][i] = 0
	#print(P)
	
	
	#Compute Shannon entropy per node
	

	
		
		
	return S_0, S_r, z, P






###############################

######################## MATLAB
#CODICE COPIATO
#format long 




precision=10^(-3);
loops=5000;
S_0=0;

n=max(size(pp));
connectivity=sum(pp)';

avg_conn=sum(connectivity)/n;
z=connectivity/(sqrt(n*avg_conn));
oldz=zeros(n,1);

count = 0
for kk=1:70
    
    U=ones(n,1)*z';
    D=ones(n,n) + z*z';
    z=connectivity ./ (sum( ( U./D - diag(diag(U./D)) )' ) )';
    z=max(z,10^(-15));
    
	
    if max(abs((z>0).*(1-z./(oldz+(oldz==0)))))<precision
        break
    end
    oldz=z;
	count = count+1;
end





%z
%Compute link probability
P=(z*z')./(ones(n,n)+z*z');
P=P-diag(diag(P));
%P

%Compute Shannon entropy per node
Ptot=P.*log(P+(P==0))+(ones(n,n)-P).*log(ones(n,n)-P +((ones(n,n)-P)==0));
S_0=(-sum(sum(triu(Ptot,1))))/n;

%Compute Random Shannon entropy per node - teorico
L=n*avg_conn/2;
p=2*L/(n*(n-1));
S_r=-(L*log(p)+(n*(n-1)/2 - L)*log(1-p))/n;

return
