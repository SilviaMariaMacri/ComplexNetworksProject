#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:56:03 2020

@author: caterina
"""

'''
Modifica del codice preso da string per convertire i nomi delle proteine del covid

'''

#!/usr/bin/env python3

##########################################################
## For a given list of proteins the script resolves them
## (if possible) to the best matching STRING identifier
## and prints out the mapping on screen in the TSV format
##
## Requires requests module:
## type "python -m pip install requests" in command line
## (win) or terminal (mac/linux) to install the module
###########################################################

import requests ## python -m pip install requests
import pandas as pd
import os




#Legge file.csv e lo mette in un dataframe
directory= '/home/caterina/Scaricati'
os.chdir(directory)

string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "get_string_ids"

##
## Set parameters
##

 

covid=pd.read_excel('2020-03-18_Krogan_SARSCoV2_27baits.xlsx',usecols=['Preys'])

protein_covid=[]


#mette la colonna delle proteine del file del covid in una lista

for i in range(len(covid)):
    protein_covid.append(covid.iloc[i][0]) 

params = {

  #"identifiers" : "\r".join(["p53", "BRCA1", "cdk2", "Q99835"]), # your protein list
  "identifiers" : "\r".join(protein_covid), # passo la nostra lista
  "species" : 9606, # species NCBI identifier 
    "limit" : 1, # only one (best) identifier per input protein
    "echo_query" : 1, # see your input identifiers in the output
    "caller_identity" : "www.awesome_app.org" # your app name

}

##
## Construct URL
##


request_url = "/".join([string_api_url, output_format, method])

##
## Call STRING
##

results = requests.post(request_url, data=params)

##
## Read and parse the results
##

#array che conterrà tutte le proteine nella notazione 9606....
covid_hitproteins=[]

for line in results.text.strip().split("\n"):
    l = line.split("\t")
    input_identifier, string_identifier = l[0], l[2]
    print("Input:", input_identifier, "STRING:", string_identifier, sep="\t")
   
    
    covid_hitproteins.append(string_identifier)