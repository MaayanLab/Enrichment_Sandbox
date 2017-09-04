import csv
import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import enrichment_methods as m
from setup import convert_gmt
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, \
	RandomTreesEmbedding, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC
from xgboost import XGBClassifier
from random import uniform as rand
from setup import open_csv
import scipy.stats as stats
from get_classifiers import get_classifiers
import h5py

def clean(tf):
	'''Extracts the transcription factor name from the name of a gmt experiment.'''
	return str(tf).partition('_')[0].partition(' ')[0].upper()

def get_overlaps(label_tfs, feature_tfs):
	'''Returns a list of label library experiments whose transcription factors are also in any feature library experiment.
	label_tfs: list-like object
	feature_tfs: list-like object
	'''
	cleaned_overlaps = {clean(x) for x in feature_tfs} & {clean(x) for x in label_tfs}
	l_in_overlaps = [x for x in label_tfs if clean(x) in cleaned_overlaps]
	return l_in_overlaps

def clean_CREEDS_Drugs(drug):
	'''Like clean(), but for drugs.'''
	return str(drug).partition(' GSE')[0].partition(' ')[0].partition(' ')[0].upper()

def get_drug_overlaps(pair):
	'''Returns a list of label library experiments whose corresponding drugs are also in any feature library experiment.'''
	if pair['l'].index.name == 'CREEDS_Drugs': creeds, drugbank = pair['l'], pair['f']
	else: drugbank, creeds = pair['l'], pair['f']
	print(creeds.shape, drugbank.shape)

	label_drugs = pair['l'].columns.values

	cleaned_overlaps = {str(x).upper() for x in drugbank.columns.values} & {clean_CREEDS_Drugs(x) for x in creeds.columns.values}
	if pair['l'].index.name == 'CREEDS_Drugs': l_in_overlaps = [x for x in label_drugs if clean_CREEDS_Drugs(x) in cleaned_overlaps]
	else: l_in_overlaps = [x for x in label_drugs if x.upper() in cleaned_overlaps]
	return l_in_overlaps

def get_methods_and_params(l,f, l_name, f_name):
	'''
	Returns a dataframe with methods and their parameters.
	l : pandas.DataFrame
		the "label" gmt library, from which tfs are being used as input gene sets
	f: pandas.DataFrame
		the "feature" gmt library, whose tfs are being ranked by the enrichment method
	l_name : str
		the name of l, for example "ChEA_2016"
	f_name: str
		the name of f, for example "CREEDS"
	'''

	#(Define any other variables, as necessary, here.)
	train_group = f
	features = f.columns.values

	#Create a classifier - only necessary for ML_fisher_features.
	#Otherwise, comment out. 
	#classifier = get_classifiers(l_name, f_name)

	#=====================================================================================================================================
	#This is where you choose which enrichment methods to run, and which paramaters to use!
	#Create a dataframe with index ['method', 'params'], and columns as the names of each enrichment method/param combo.
	#For example, a column with name 'Foo3' might have 'method' m.Foo() and 'params' bootstrap = False.
	#See enrichment_methods.py for the available methods and necessary params.
	#You must specify ALL the params EXCEPT for l_tf_genes, which is initialized and called later, in enrichment_wrapper().
	#=====================================================================================================================================
	df = pd.DataFrame(index=['func', 'params'])
	df['Pairwise_Gini'] = [m.pairwise_impurity, [f, 'Gini']] 

	return df

def enrichment_wrapper(pair):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : dict
		key 'l' contains the label library df, and key 'f' contains the feature library df. 
	'''
	l_name, f_name = pair['l'].index.name, pair['f'].index.name
	print('Beginning enrichment analysis from', l_name, 'to', f_name)

	#Get the label library experiments whose transcription factors (or drugs) are also in any feature library experiment.
	if 'Drugs' in l_name: overlaps = get_drug_overlaps(pair)
	else: overlaps = get_overlaps(pair['l'].columns.values, pair['f'].columns.values)

	#Get the methods and parameters with which to perform enrichment. 
	methods_and_params = get_methods_and_params(pair['l'], pair['f'], l_name, f_name)
	output_heading = 'from_' + l_name + '_to_' + f_name

	#Iterate over each method.
	for column in methods_and_params:
		mp = methods_and_params[column]

		#Some methods actually return multiple results. These will need multiple output files.
		if mp.name == 'ZAndCombined': output_fnames = (output_heading + '_Z.csv', output_heading + '_Combined.csv')
		elif mp.name == 'Pairwise_Gini': output_fnames = (output_heading + '_Pair_Gini_1.csv', 
			output_heading + '_Pair_Gini_2.csv',output_heading + '_Pair_Gini_3.csv',output_heading + '_Pair_Gini_4.csv')
		else: output_fnames = (output_heading + '_' + mp.name + '.csv',)

		#Check if the file has already been created.
		if os.path.isfile(output_fnames[0]): print('ranking file already created for', mp.name)

		#If not, create it. 
		else:
			#Use dataframes to store results after each tf iteration.
			dfs = [pd.DataFrame() for n in range(len(output_fnames))]
			#Iterate over each tf in the overlaps.
			for l_tf in overlaps:
				print(mp.name, l_tf) #for diagnostics.
				l_tf_genes = pair['l'].index[pair['l'][l_tf]]
				result = mp['func'](l_tf_genes, *mp['params'])
				for x in range(len(dfs)): 
					df = dfs[x]
					#Store result as a column in the df.
					if len(dfs) == 1: df[l_tf] = result
					else: df[l_tf] = result[x]
			for x in range(len(dfs)): 
				#Save the dfs.
				df = dfs[x]
				df.index = pair['f'].columns
				df.to_csv(output_fnames[x], sep='\t')
	return

if __name__ == '__main__':
	all_libs = ['ENCODE_TF_ChIP-seq_2015', 'ChEA_2016', 'CREEDS']
	#all_libs = ['DrugBank','CREEDS_Drugs']

	#Get dataframes of each gmt library in all_libs
	os.chdir('libs')
	all_dfs = {x:convert_gmt('df', x) for x in all_libs}
	os.chdir('..')
	if not os.path.isdir('results'): os.makedirs('results')
	os.chdir('results')

	#Iterate over each gmt pair.
	lib_df_pairs = [{'l':all_dfs[a], 'f':all_dfs[b]} for a in all_libs for b in all_libs if a != b]
	Parallel(n_jobs=1, verbose=0)(delayed(enrichment_wrapper)(pair)for pair in lib_df_pairs)