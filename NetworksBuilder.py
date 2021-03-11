import networkx as nx
import pandas as pd





  
'Create Human PPI UNDIRECTED GRAPH'  

# Input:  STRING file 
#         minimum link score (= integer between 0 and 1000) 
  
def HumanGraph(NameHumanFile, threshold_score, graph_name): 
 	
	
	GH = nx.Graph() 
	
	sapiens = pd.read_csv(NameHumanFile, sep=" ", usecols=['protein1', 'protein2', 'combined_score']) 
	sapiens=sapiens[sapiens['combined_score']>threshold_score]
	
	for i in range(len(sapiens)):
		GH.add_edge(sapiens.iloc[i][0], sapiens.iloc[i][1], weight=sapiens.iloc[i][2])
	
	GH.name = graph_name
	
	print(nx.info(GH))
	
	return GH
	
	
	
	
	
'Create Virus DIRECTED GRAPH'

# Input:   STRING file 

def VirusGraph(NameVirusFile, graph_name):
	
	
	GV = nx.DiGraph() 
	
	virus = pd.read_csv(NameVirusFile, sep="\t", usecols=['node1_external_id','node2_external_id', 'combined_score'])
	
	for i in range(len(virus)):
		GV.add_edge(virus.iloc[i][0], virus.iloc[i][1], weight=virus.iloc[i][2])
	
	GV.name = graph_name
	
	print(nx.info(GV))
	
	return GV





'Create SUBGRAPH of a network (GH) with nodes of an other network (GV)'
 
# Input:   GH = directed/undirected graph
#          GV = directed/undirected graph


 

def SubGraph(GH,GV):
 
	Gsub = GH.subgraph(nx.nodes(GV))
	Gsub.name = GV.name
	
	print(nx.info(Gsub))
	
	return Gsub
