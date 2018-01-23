import csv
import os
import pickle
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import enrichment_functions as m
from setup import open_gvm
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, \
	RandomTreesEmbedding, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC
from xgboost import XGBClassifier
from random import uniform as rand
import scipy.stats as stats
#from get_classifiers import get_classifiers #REMOVED. the script was put in the old or unused scripts folder.
import h5py

def list_of_drug_libs():
	'''Generates and returns a tuple of drug library names.'''
	expanded_drug_libs = (
		'repurposing_drugs_20170327',
		'interactions', 
		'1_DrugBank_Edgelist_10-05-17', 
		'2_TargetCentral_Edgelist_10-05-17',
		'3_Edgelists_Union_10-05-17', 
		'4_EdgeLists_Intersection_10-05-17')
	ppi_libs = (
		'hu.MAP',
		'BioGRID',
		'ARCHS4')

	non_expanded_drug_libs = (
		'CREEDS_Drugs',
		'DrugMatrix_Union')

	drug_libs = tuple(
		edl + '_expanded_with_' + ppi for edl in expanded_drug_libs for ppi in ppi_libs) + non_expanded_drug_libs

	return drug_libs

def lib_type(lib_name):
	'''Uses a library's name to determine its annotation type, i.e. transcription factor or drug.'''
	tf_libs = (
		'CREEDS_TFs', 
		'ChEA_2016',
		'ENCODE_TF_ChIP-seq_2015',
		'Single_Gene_Perturbations_from_GEO_down',
		'Single_Gene_Perturbations_from_GEO_up')

	drug_libs = list_of_drug_libs()

	if lib_name in tf_libs: return 'tf'
	elif lib_name in drug_libs: return'drug'
	else: raise ValueError(lib_name, 'is not a suppored gene set library.') 

def clean(annot, lib_type):
	'''
	Extracts either the transcription factor name or the drug name from the annotation.
	annot : str
		The annotation to be cleanend.
	lib_type : str
		Either 'tf' or 'drug', depending on the library from which the annotation was taken.
	'''
	annot = str(annot)
	if lib_type == 'tf': return annot.partition('_')[0].partition(' ')[0].upper()
	elif lib_type == 'drug': return annot.partition('|||')[0]
	else: raise ValueError(lib_type, 'is not a valid library type.')

def get_common_ilib_annots(ilib_name, ilib_annots, slib_name, slib_annots):
	'''Return the annotations in the input library which have matches in the search library.'''
	cleaned_ilib_annots = {clean(annot, lib_type(ilib_name)) for annot in ilib_annots}
	cleaned_slib_annots = {clean(annot, lib_type(slib_name)) for annot in slib_annots}
	cleaned_overlaps = cleaned_ilib_annots & cleaned_slib_annots
	common_ilib_annots = [annot for annot in ilib_annots if clean(annot, lib_type(ilib_name)) in cleaned_overlaps]
	return common_ilib_annots

def get_enrichment_algorithms(ilib_gvm, ilib_name, slib_gvm, slib_name):
	'''
	Returns a dataframe with methods and their parameters.
	ilib_gvm : pandas.DataFrame
		the input library's gene vector matrix, from which column vectors are being used as input gene sets
	slib_gvm: pandas.DataFrame
		the search library's gvm, whose annotations are being ranked by the enrichment algorithm based on their gene sets
	ilib_name : str
		the name of ilib_gvm, for example "ChEA_2016"
	slib_name : str
		the name of slib_gvm, for example "CREEDS_TFs"
	algorithm_name_only : bool
		If True, returns only the names of the algorithms. 
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
	#==================================================

	#======================================================================================================================
	#This is where you specify which enrichment methods to run.
	#Specify each method as a column in a dataframe, where the first row is the enrichment function,
	#	the second row is a tuple of parameters, and the column name is the chosen name of the enrichment method.
	#For example, for a method using the RandomForestClassifier with max_depth=3, you might say:
	#	df['RF_md3'] = (m.ML_wrapper, (RandomForestClassifier, train_group, features, 10304, max_depth=3))
	#See enrichment_functions.py for the available methods and necessary params.
	#IMPORTANT: DO NOT specify the input_geneset parameter:
	#	this is created and called later, in enrichment().
	#======================================================================================================================
	enrichment_algorithms = pd.DataFrame(index=['func', 'params'])
	enrichment_algorithms['Fisher'] = [m.Fisher, [slib_gvm]] 
	#enrichment_algorithms['RandomForest'] = [m.ML_wrapper, [RandomForestClassifier, train_group, features, 101317]]
	#======================================================================================================================
	return enrichment_algorithms

def enrichment(pair):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : tuple of str
		(input library fname, search library fname)
	'''
	ilib_fname, slib_fname = pair
	ilib_name, slib_name = (fname.rpartition('\\')[2].partition('_gvm')[0] for fname in pair)
	prefix = 'input_' + ilib_name + '_into_' + slib_name

	script_dir = os.path.dirname(os.path.abspath(__file__))
	results_dir = os.path.join(script_dir, 'results\\')

	#==============================================================================================================
	#OPTIONAL: To speedily allow this script to continue working from where it left off when it was previously run,
	#	specify ALL the enrichment algorithm names again here.
	#==============================================================================================================
	algorithm_names = ('Fisher',) #e.g. `('Fisher', 'RandomForest')` or `None`. 
	#==============================================================================================================

	#Exit if enrichment between this library pair has already been done.
	#See above chunk: `algorithm_names` must be specified.
	if algorithm_names is not None:
		output_fnames = tuple(results_dir + prefix + '_' + name + '.csv' for name in algorithm_names)
		if all((os.path.isfile(fname) for fname in output_fnames)):
			print('Already done enrichment for', prefix)
			return

	#Otherwise, begin by loading the gvms.
	print('Beginning enrichment analysis inputting', ilib_name, 'into', slib_name)
	ilib_gvm, slib_gvm = (open_gvm(fname) for fname in pair)

	#Get the algorithms with which to perform enrichment. 
	enrichment_algorithms = get_enrichment_algorithms(ilib_gvm, ilib_name, slib_gvm, slib_name)

	#Get the input library annotations whose corresponding tf/drug also corresponds to 
	#	at least one search library annotation.
	common_ilib_annots = get_common_ilib_annots(
		ilib_name, ilib_gvm.columns.values, 
		slib_name, slib_gvm.columns.values)
	print(str(len(common_ilib_annots)), 'overlaps')

	#Iterate over each algorithm (i.e. each column).
	for algorithm_name in enrichment_algorithms:
		algorithm = enrichment_algorithms[algorithm_name]

		#Some methods actually return multiple results. These will need multiple score files.
		if algorithm_name == 'ZAndCombined': output_fnames = (
			results_dir + prefix + '_Z.csv', 
			results_dir + prefix + '_Combined.csv')
		#========================================================================================
		#If using a algorithm which returns multiple results, add the appropriate elif statement here.
		#========================================================================================
		#elif algorithm_name == 'Pairwise_Gini': output_fnames = (prefix + '_Pair_Gini_ltf100_1.csv', 
		#	prefix + '_Pair_Gini_ltf100_5.csv', prefix + '_Pair_Gini_ltf100_10.csv', prefix + '_Pair_Gini_ltf100_25.csv')
		#========================================================================================
		else: output_fnames = (results_dir + prefix + '_' + algorithm_name + '.csv',)

		#Exit if this algorithm has already been run for this library pair.
		if os.path.isfile(output_fnames[0]): 
			print('score file already created for', algorithm_name)
			return

		#Otherwise, run enrichment.
		#Results will be stored after iteration over the common ilib annotations.
		score_dfs = [pd.DataFrame() for n in range(len(output_fnames))]
		for annot in common_ilib_annots:
			print(algorithm_name, annot) #for checking progress.
			input_geneset = ilib_gvm.index[ilib_gvm[annot]]
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
			df.index = slib_gvm.columns
			df.to_csv(output_fnames[x], sep='\t')
	return

if __name__ == '__main__':
	tf_libs = ('ENCODE_TF_ChIP-seq_2015', 'ChEA_2016', 'CREEDS_TFs')

	drug_libs = list_of_drug_libs()

	#======================================================
	#Choose the libraries with which to perform enrichment.
	#======================================================
	libs = tf_libs
	#======================================================

	#Perform enrichment over all non-identical pairs of the chosen libraries.
	script_dir = os.path.dirname(os.path.abspath(__file__))
	libs_dir = os.path.join(script_dir, 'libs\\')
	libs_fnames = tuple(libs_dir+lib+'_gvm.csv' for lib in libs)
	lib_pairs = tuple((a,b) for a in libs_fnames for b in libs_fnames if a != b)
	Parallel(n_jobs=1, verbose=0)(delayed(enrichment)(pair) for pair in lib_pairs)