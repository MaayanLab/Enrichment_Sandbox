import csv
import os
import pickle
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

def clean_wrapper(i, lib_name):
	'''More general version of `clean()` which works with drug libraries.'''
	if lib_name == 'CREEDS_Drugs': 
		#The library is CREEDS_Drugs.
		return str(i).partition(' GSE')[0].partition(' ')[0].partition(' ')[0].lower() 
	elif ('10-05-17' in lib_name) or ('Drug' in lib_name):
		#The library is some other drug library.
		return str(i).lower() 
	else:
		#The library is a transcription factor library.
		return clean(i) 

def get_overlaps(l_lib_name, l_tfs, f_lib_name, f_tfs):
	'''Return the transcription factors in the label library which have matches in the feature library.'''
	l = {clean_wrapper(i, l_lib_name) for i in l_tfs}
	f = {clean_wrapper(i, f_lib_name) for i in f_tfs}
	overlaps = l & f
	overlaps_in_l = [i for i in l_tfs if clean_wrapper(i, l_lib_name) in overlaps]
	return overlaps_in_l

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

	#======================================================================================================================
	#This is where you choose which enrichment methods to run, and which paramaters to use!
	#Create a dataframe with index ['method', 'params'], and columns as the names of each enrichment method/param combo.
	#For example, a column with name 'Foo3' might have 'method' m.Foo() and 'params' bootstrap = False.
	#See enrichment_methods.py for the available methods and necessary params.
	#You must specify ALL the params EXCEPT for l_tf_genes, which is initialized and called later, in enrichment_wrapper().
	#======================================================================================================================
	df = pd.DataFrame(index=['func', 'params'])
	df['Fisher'] = [m.Fisher, [f]] 
	df['RandomForest'] = [m.ML_wrapper, [RandomForestClassifier, train_group, features, 101317]]
	return df
	#======================================================================================================================

def enrichment_wrapper(pair):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : dict
		key 'l' contains the label library df, and key 'f' contains the feature library df. 
	'''
	l_name, f_name = pair['l'].index.name, pair['f'].index.name
	print('Beginning enrichment analysis from', l_name, 'to', f_name)

	#Get the label library experiments whose transcription factors (or drugs) are also in any feature library experiment.
	overlaps = get_overlaps(l_name, pair['l'].columns.values, f_name, pair['f'].columns.values)
	print(str(len(overlaps)), 'overlaps')

	#Get the methods and parameters with which to perform enrichment. 
	methods_and_params = get_methods_and_params(pair['l'], pair['f'], l_name, f_name)
	output_heading = 'from_' + l_name + '_to_' + f_name

	#Iterate over each method.
	for column in methods_and_params:
		mp = methods_and_params[column]

		#Some methods actually return multiple results. These will need multiple output files.
		if mp.name == 'ZAndCombined': output_fnames = (output_heading + '_Z.csv', output_heading + '_Combined.csv')
		elif mp.name == 'Pairwise_Gini': output_fnames = (output_heading + '_Pair_Gini_ltf100_1.csv', 
			output_heading + '_Pair_Gini_ltf100_5.csv',output_heading + '_Pair_Gini_ltf100_10.csv',output_heading + '_Pair_Gini_ltf100_25.csv')
		elif mp.name == 'Pairwise_Gini_weighted': output_fnames = (output_heading + '_Pair_Gini_ltf100_w_1.csv', 
			output_heading + '_Pair_Gini_ltf100_w_5.csv',output_heading + '_Pair_Gini_ltf100_w_10.csv',output_heading + '_Pair_Gini_ltf100_w_25.csv')
		elif mp.name == 'Pairwise_Gini_nocset4': output_fnames = (output_heading + '_Pair_Gini_ltf100_n_1.csv', 
			output_heading + '_Pair_Gini_ltf100_n_5.csv',output_heading + '_Pair_Gini_ltf100_n_10.csv',output_heading + '_Pair_Gini_ltf100_n_25.csv')
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
	tf_libs = ('ENCODE_TF_ChIP-seq_2015_abridged', 'ChEA_2016_abridged', 'CREEDS_abridged')
	drug_libs = ('1_DrugBank_EdgeList_10-05-17', 
		'2_TargetCentral_EdgeList_10-05-17',
		'3_EdgeLists_Union_10-05-17', 
		'4_EdgeLists_Intersection_10-05-17',
		'DrugBank',
		'CREEDS_Drugs')

	#========================================================
	#Choose which libraries with which to perform enrichment.
	#========================================================
	libs = drug_libs
	#========================================================

	#Get dataframes of each gmt library in libs
	os.chdir('libs')
	all_dfs = {x:convert_gmt('df', x) for x in libs}
	os.chdir('..')
	if not os.path.isdir('results'): os.makedirs('results')
	os.chdir('results')

	#Iterate over each gmt pair.
	lib_df_pairs = [{'l':all_dfs[a], 'f':all_dfs[b]} for a in libs for b in libs if a != b]
	Parallel(n_jobs=4, verbose=0)(delayed(enrichment_wrapper)(pair) for pair in lib_df_pairs)