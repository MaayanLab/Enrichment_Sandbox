import requests
import csv
import os
import h5py
import pickle
import numpy as np
import pandas as pd
from joblib import Parallel, delayed 
from setup import open_csv, convert_gmt

def edgelist_to_gmt(file):
	#check if this has already been done. 
	output_fname = file.partition('.csv')[0]+ '_transformed.csv'
	if os.path.isfile(output_fname): 
		print(output_fname, 'already created.')
		return

	#otherwise, proceed.
	print(file)

	#read file.
	f = pd.read_csv(file, sep=',', encoding='latin1')

	#create a new dataframe using the drug IDs as the index. 
	#I chose to use drug IDs instead of drug names so that I don't have to worry
		#about parsing strings.
	drugs = set(f['DrugID_ChEMBL'])
	new_df = pd.DataFrame(columns=drugs, dtype=bool)

	#we are essentially going to one-hot encode the targets of each drug by
		#iterating over the rows of the edgelist. 
	for i in range(f.shape[0]):
		#(for checking on progress)
		if i % 100 == 0: print(i)
		new_df.at[f.at[i,'TargetGeneSymbol_Entrez'], 
			f.at[i,'DrugID_ChEMBL']] = True

	#save the results.
	new_df.to_csv(file.partition('.csv')[0]+ '_transformed.csv', sep='\t')
	return

libs = ('1_DrugBank_EdgeList_10-05-17.csv', 
	'2_TargetCentral_EdgeList_10-05-17.csv',
	'3_EdgeLists_Union_10-05-17.csv', 
	'4_EdgeLists_Intersection_10-05-17.csv')

os.chdir('libs')
for lib in libs:
	edgelist_to_gmt(lib)
os.chdir('..')