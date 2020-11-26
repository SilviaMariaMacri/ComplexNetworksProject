'COPIATO DA STRING'



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

string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "get_string_ids"

##
## Set parameters
##

params = {

    "identifiers" : "\r".join(hitnodes_covid_array),#["p53", "BRCA1", "cdk2", "Q99835"]), # your protein list
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


'FIN QUI COPIATO, LE RIGHE SEGUENTI LEGGERMENTE MODIFICATE'


##
## Read and parse the results
##
name_identifier = []
name_string = []
for line in results.text.strip().split("\n"):
	l = line.split("\t")
	input_identifier, string_identifier = l[0], l[2]
	name_identifier.append(input_identifier)
	name_string.append(string_identifier)
	print("Input:", input_identifier, "STRING:", string_identifier, sep="\t")
	
	
	





# save everything on a dataframe
columns_names = {'identifier': name_identifier, 'string_name': name_string}
names = pd.DataFrame(data=columns_names)

	
	
	
# save dataframe in a file .txt	
names.to_csv('ProteinNamesString.txt', index=False)
	
	
	
	
	
	
	
	
	
	