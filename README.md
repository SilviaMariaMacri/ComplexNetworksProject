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
- *HumanGraph* takes as input the STRING file, the lower limit on the link score we want to cosider and the name of the network that we want to create; the output is an undirected graph.
- VirusGraph 

#### NetworkCharacterization.py
#### Percolation.py
#### EntropyCanonical.py
#### MainCode.py


