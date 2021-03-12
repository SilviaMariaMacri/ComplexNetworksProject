# Structure of the project

### Data files downloaded from STRING database:
- '9606.protein.links.v11.0.txt.gz', representing Human PPI network. 
(https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Homo+sapiens)
- 15 files .tsv, representing Human-Virus PPI networks of the following viruses: Cytomegalo, Dengue type 2, Ebola, Hepatitis B, HIV1, HPV type 1a, HTLV1, Influenza A, Lassa virus, MARV, Mumps virus, Human parechovirus 2, SARS-CoV, Varicella zoster virus, WNV. (http://viruses.string-db.org)



### Python codes:
- NetworkBuilder.py
- NetworkCharacterization.py
- Percolation.py
- EntropyCanonical.py
- MainCode.py


# Description of the codes

#### NetworkBuilder.py
Creation of the graph from the STRING file. It consists of three functions:
- *HumanGraph* takes as input the STRING file, the lower limit on the link score we want to cosider and the name of the network that we want to create; the output is an undirected graph. It is used to create the Human PPI network.
- *VirusGraph* takes as input the STRING file and the name of the network that we want to create; the output is a directed graph. It is used to create a Human-Virus PPI network. 
- *SubGraph* takes as input two graphs and returns a subgraph of the first input network consisting of the nodes of the second input network. It is used to create a Human PPI subgraph only with the nodes linked to viral nodes of a Human-Virus PPI network.


#### NetworkCharacterization.py
Analysis of a network through measures of degree K, average neighbor degree Knn, betweenness centrality BC and closeness centrality CL. It consists of six functions:
- *NetworkCharacterization* takes as input a graph; the output is different in case of undirected and directed graph. In case of undirected graph (in our project representing a Human-Virus PPI network), it returns an array of three dataframes: the first one has K, vKin, Kout, Knn,in, Knn,out, BC and CL as column and each row corresponding to a node; the second one is a selection of the rows of the first dataframe corresponding to human proteins; the third one is a selection of the rows of the first dataframe corresponding to virus protein. In case of directed graph (in our project representing a Human PPI subgraph), it returns a dataframe with K, Knn, BC anc CL as columns and each row corresponding to a node.
- *PlotDegreeHist* returns the plot of histogram degree; it takes as input the dataframe defined by the first function and the title of the plot. If the dataframe refers to a directed network, two plots of in degree histogram and out degree histogram are distinguished. 
- *PlotDegreeNNvsDegree* returns a scatter plot of Knn vs K; it takes as input the dataframe defined by the first function and the title of the plot. If the dataframe refers to a directed network, two plots of Knn,in vs Kin and Knn,out vs Kout are distinguished.
- *PlotINvsOUT* returns a scatter plot of out degree vs in degree; it takes as input the dataframe defined by the first function and the title of the plot. It raises Error if the dataframe refers to an undirected network.
- *PlotBcClvsDegree* returns two scatter plots, the first one of BC vs K and the second one of CL vs K. It takes as input the dataframe defined by the first function, the title of the plot and an integer number n. The values of K, BC and CL are also averaged over n number of nodes (ordered by degree). If the dataframe refers to an undirected network, Kin is considered. 
- *Mean* takes as input an array of numbers and returns the mean with the root mean square deviation. It is used to calculate the average values of K, BC and CL in the previous function.


#### Percolation.py
It consists of the function *Percolation*. It takes as input a graph, the corresponding dataframe obtained by the function *NetworkCharacterization* of the previous file, and the name of the network and performs a percolation analysis on the given network. The size of the giant component is calculated after each removal of a single node. The percolation analysis is performed three times, each one following a different node removal order; in particular the nodes are ordered by descending betweenness centrality, descending degree and randomly. The function returns a plot of the size of the giant component vs the percentage of removed nodes and a three arrays (one for each removal method) containing the size of the giant component after each removal. 


#### EntropyCanonical.py
It calculates the Shannon entropy values of a series of networks and compares them. It consists of three functions:
- *EntropyCanonicalUndirected* takes as input the adjacency matrix of a (weighted) undirected graph. It gives as output the canonical Shannon entropy S0, calculated by fixing the degree distribution, the canonical Shannon entropy Sr, calculated by fixing the total link probability, the lagrange multipliers of degree sequence constraints and the link probability matrix. The code is the translation in Python language of the MATLAB code *entropy_canonical.m* saw at lesson.
- *EntropyCanonicalDirected* is the adaptation of the previous function in the case of a (weighted) directed network. In this case two sets of lagrange multipliers are calculated, one corresponding to the in connectivity and the other corresponding to the out connectivity; the link probability matrix, used to calculate S0, is obtained by the superposition of the upper triangular link probability matrix corresponding to the out connectivity and the lower triangular one corresponding to the in connectivity. Sr in calculated by fixing the average connectivity as the mean of the in connectivity and the out connectivity.
- *EntropyDifference* calculates Shannon entropy values of different networks, by using the two previous functions, and compares them. It takes as input two lists of networks that we want to compare (in our project they are respectively the list of Human-Virus PPI networks and the list of the corresponding Human PPI subnetworks) and the corresponding list of names of each pair of networks. It returns an array with S0 and Sr values for each group of network as columns and each row corresponding to a particular network pair. Moreover, it returns three plots: the first one shows the difference of the entropy values between each pair of network, the second and the third ones show the absolute entropy values for each graph of the two groups of networks.



#### MainCode.py
It imports all the functions of the previously described Python files, in order to obtain a complete analysis of all the networks. It defines the list of STRING files that we want to use and the corresponding list of virus names and creates all the networks through iterations of NetworkBuikder.py functions; Later, it iterates the functions of NetworkCharacterization.py, Percolation.py and EntropyCanonical.py over all the graphs.
