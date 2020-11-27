'leggere dati covid'

import numpy as np
import pandas as pd
import networkx as nx



# import covid data as dataframe
covid = pd.read_excel('2020-03-18_Krogan_SARSCoV2_27baits_LowThreshold.xlsx', sep=",", 
                    usecols=['Bait', 'Preys', 'MIST']) 

# create network of covid-human protein interactions
GC = nx.DiGraph() 
for i in range(len(covid)):
	GC.add_edge(covid.iloc[i][0], covid.iloc[i][1], weight=covid.iloc[i][2])


# control if nodes are repeated in the dataframe	
#covid.Bait.value_counts()
#covid.Preys.value_counts()



'''
UNICA PROTEINA UMANA RIPETUTA DUE VOLTE IN 'Preys'

In: covid[covid['Preys']=='Q9Y5L0']
Out: 
                Bait   Preys      MIST
59       SARS-CoV2 M  Q9Y5L0  0.600814
330  SARS-CoV2 orf7a  Q9Y5L0  0.616847
'''





#### NON SERVE PIù:

# create dataframe and/or array of hit human proteins	
#dataframe	
hitnodes_covid = covid['Preys'] #it is still a dataframe
hitnodes_covid = hitnodes_covid.drop(59) #remove the only repeated protein
#array
hitnodes_covid_array = np.array(hitnodes_covid) #turn it into array





# import as dataframe the file with protein names
name_protein = pd.read_csv('ProteinNamesString.txt', sep=",") 



# ALTRO PROGRAMMA E FILE DA CONTROLLARE !!!!!!
