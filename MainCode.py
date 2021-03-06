#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 21:33:14 2021

@author: caterina
"""
import NetworksBuilder

NameHumanFile='9606.protein.links.v11.0.txt'

NameVirusFile = ['string_interactions_WNV.tsv','string_interactions_varicella.tsv'				 
				 ,'string_interactions_SARSCov.tsv','string_interactions_parechovirus2.tsv',
				 'string_interactions_mumps.tsv','string_interactions_MARV.tsv',
				 'string_interactions_lassa.tsv','string_interactions_InfluenzaA.tsv',
				 'string_interactions_HTLV-1.tsv','string_interactions_HPV1a.tsv',
				 'string_interactions_HIV1_553.tsv','string_interactions_hepatitisB.tsv',
				 'string_interactions_ebola.tsv','string_interactions_dengue2.tsv',
				 'string_interactions_cytomegalo.tsv','Covid19.txt']

VirusNames=['WNV','Varicella zoster virus','SARS-CoV','Human parechovirus 2',
			'Mumps virus','MARV','Lassa virus','Influenza A','HTLV1','HPV type 1a',
			'HIV1','Hepatitis B','Ebola','Dengue type 2','Cytomegalo','SARS-CoV-2']

GH=NetworksBuilder.HumanGraph(NameHumanFile,250) #human PPI network

GV=[] #list of human-virus PPI
for i in range (len(NameVirusFile)):
	GV.append(NetworksBuilder.VirusGraph(NameVirusFile[i],VirusNames[i]))


Gsub=[] #list of human PPI subnetworks related to each virus
for i in range (len(NameVirusFile)):
	Gsub.append(NetworksBuilder.SubnetworkGraphs (GH,GV[i], VirusNames[i]))