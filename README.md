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
- *Mean* takes as input an array of numbers and returns the mean with error. It is used to calculate the average values of K, BC and CL in the previous function.+ 
#### Percolation.py
#### EntropyCanonical.py
#### MainCode.py


