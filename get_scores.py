import csv
import os
import pickle
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import enrichment_functions as m
from setup import convert_gmt
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, \
	RandomTreesEmbedding, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC
from xgboost import XGBClassifier
from random import uniform as rand
from setup import open_csv
import scipy.stats as stats
#from get_classifiers import get_classifiers #REMOVED. the script was put in the old or unused scripts folder.
import h5py

def clean_tf(annot):
	'''Extracts the transcription factor name from an annotation.'''
	return str(annot).partition('_')[0].partition(' ')[0].upper()

def clean(annot, libname):
	'''
	More general version of `clean_tf()` which checks if the lib has tf or drug annotations based on its name,
	and then extracts either the tf name or the drug name from the inputted annotation.
	'''
	if libname == 'CREEDS_Drugs': 
		return str(annot).partition(' GSE')[0].partition(' ')[0].partition(' ')[0].lower() 
	elif 'DrugMatrix' in libname:
		for parser in ('mg/kg', '_uM_', 'ng/ml'):
			if parser in annot:
				return str(annot).partition(parser)[0].rpartition('-')[0].lower()
		raise ValueError(annot + ' could not be parsed.')
	elif ('10-05-17' in libname) or ('Drug' in libname):
		#The library is some other drug library.
		return str(annot).lower() 
	else:
		#The library is a transcription factor library.
		return clean_tf(annot)

def get_overlapping_ilib_annots(ilib_name, ilib_annots, slib_name, slib_annots):
	'''Return the annotations in the label library which have matches in the feature library.'''
	cleaned_ilib_annots = {clean(annot, ilib_name) for annot in ilib_annots}
	cleaned_slib_annots = {clean(annot, slib_name) for annot in slib_annots}
	cleaned_overlaps = cleaned_ilib_annots & cleaned_slib_annots
	overlapping_ilib_annots = [annot for annot in ilib_annots if clean(annot, ilib_name) in cleaned_overlaps]
	return overlapping_ilib_annots

def get_enrichment_algorithms(ilib_gvm, slib_gvm, ilib_name, slib_name):
	'''
	Returns a dataframe with methods and their parameters.
	ilib_gvm : pandas.DataFrame
		the input library's gene vector matrix, from which column vectors are being used as input gene sets
	slib_gvm: pandas.DataFrame
		the search library's gvm, whose annotations are being ranked by the enrichment algorithm based on their gene sets
	ilib_name : str
		the name of ilib_gvm, for example "ChEA_2016"
	slib_name: str
		the name of slib_gvm, for example "CREEDS"
	'''
	#=================================================
	#Define any other variables, as necessary, here.
	#=================================================
	train_group = slib_gvm
	features = slib_gvm.columns.values
	#REMOVED-DOES NOT WORK. see get_classifiers.py in the old or unused scripts folder.
	##Create a classifier - only necessary for ML_fisher_features (REMOVED).
	##Otherwise, comment out. 
	#classifier = get_classifiers(ilib_name, slib_name)

	#======================================================================================================================
	#This is where you specify which enrichment methods to run.
	#Specify each method as a column in a dataframe, where the first row is the enrichment function,
	#	the second row is a tuple of parameters, and the column name is the chosen name of the enrichment method.
	#For example, for a method using the RandomForestClassifier with max_depth=3, you might say:
	#	df['RF_md3'] = (m.ML_wrapper, (RandomForestClassifier, train_group, features, 10304, max_depth=3))
	#See enrichment_functions.py for the available methods and necessary params.
	#IMPORTANT: DO NOT specify the input_geneset parameter:
	#	this is created and called later, in perform_enrichment().
	#======================================================================================================================
	enrichment_algorithms = pd.DataFrame(index=['func', 'params'])
	enrichment_algorithms['Fisher'] = [m.Fisher, [slib_gvm]] 
	enrichment_algorithms['RandomForest'] = [m.ML_wrapper, [RandomForestClassifier, train_group, features, 101317]]
	return enrichment_algorithms
	#======================================================================================================================

def perform_enrichment(pair):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : dict
		key 'i' gives the input library gvm, and key 's' gives the search library gvm. 
	'''
	ilib_name, slib_name = pair['i'].index.name, pair['s'].index.name
	print('Beginning enrichment analysis inputting', ilib_name, 'into', slib_name)

	#Get the input library annotations whose corresponding tf/drug also corresponds to 
	#	at least one search library annotation.
	overlapping_ilib_annots = get_overlapping_ilib_annots(ilib_name, pair['i'].columns.values, 
		slib_name, pair['s'].columns.values)
	print(str(len(overlapping_ilib_annots)), 'overlaps')

	#Get the algorithms with which to perform enrichment. 
	enrichment_algorithms = get_enrichment_algorithms(pair['i'], pair['s'], ilib_name, slib_name)
	prefix = 'input_' + ilib_name + '_into_' + slib_name

	#Iterate over each algorithm (i.e. each column).
	for algorithm_name in enrichment_algorithms:
		algorithm = enrichment_algorithms[algorithm_name]

		#Some methods actually return multiple results. These will need multiple score files.
		if algorithm_name == 'ZAndCombined': output_fnames = (prefix + '_Z.csv', prefix + '_Combined.csv')
		#elif algorithm_name == 'Pairwise_Gini': output_fnames = (prefix + '_Pair_Gini_ltf100_1.csv', 
		#	prefix + '_Pair_Gini_ltf100_5.csv', prefix + '_Pair_Gini_ltf100_10.csv', prefix + '_Pair_Gini_ltf100_25.csv')
		#========================================================================================
		#If using a algorithm which returns multiple results, add the appropriate elif statement here.
		#========================================================================================
		else: output_fnames = (prefix + '_' + algorithm_name + '.csv',)

		#Check if the file has already been created.
		if os.path.isfile(output_fnames[0]): print('score file already created for', algorithm_name)

		#If not, create it. 
		else:
			#Results will be stored after each tf iteration.
			score_dfs = [pd.DataFrame() for n in range(len(output_fnames))]
			#Iterate over each overlapping ilib annotation.
			for annot in overlapping_ilib_annots:
				print(algorithm_name, annot) #for diagnostics.
				input_geneset = pair['i'].index[pair['i'][annot]]
				#Get scores for all the slib annotations.
				result = algorithm['func'](input_geneset, *algorithm['params'])
				for x in range(len(score_dfs)): 
					df = score_dfs[x]
					#Store this result as a column in the score df.
					if len(score_dfs) == 1: df[annot] = result
					else: df[annot] = result[x]
			for x in range(len(score_dfs)): 
				#Save the score_dfs as csv files.
				df = score_dfs[x]
				df.index = pair['s'].columns
				df.to_csv(output_fnames[x], sep='\t')
	return

if __name__ == '__main__':
	tf_libs = ('ENCODE_TF_ChIP-seq_2015.txt', 'ChEA_2016.txt', 'CREEDS_tfs.txt')
	drug_libs = ('1_DrugBank_EdgeList_10-05-17.txt', 
		'2_TargetCentral_EdgeList_10-05-17.txt',
		'3_EdgeLists_Union_10-05-17.txt', 
		'4_EdgeLists_Intersection_10-05-17.txt',
		'DrugBank',
		'CREEDS_Drugs',
		'DrugMatrix_Union',
		)

	#========================================================
	#Choose which libraries with which to perform enrichment.
	#========================================================
	libs = tf_libs
	#========================================================

	#Get dataframes of each gmt library in libs
	os.chdir('libs')
	all_dfs = {lib:convert_gmt(lib, 'gvm') for lib in libs}
	os.chdir('..')
	if not os.path.isdir('results'): os.makedirs('results')
	os.chdir('results')

	#Iterate over each gmt pair.
	lib_df_pairs = [{'i':all_dfs[a], 's':all_dfs[b]} for a in libs for b in libs if a != b]
	Parallel(n_jobs=3, verbose=0)(delayed(perform_enrichment)(pair) for pair in lib_df_pairs)