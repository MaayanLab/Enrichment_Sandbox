import os
import pandas as pd

#========================================================================================================
'''
The purpose of this script is to adjust the score or ranking files for a certain method.
'''
#========================================================================================================

def rename(str_to_change, new_str):
	'''
	Replaces str_to_change to new_str for all score files. 
	You should call remove_old_methods_from_ranking_files() after.
	'''
	for x in os.listdir(os.getcwd()):
		if (str_to_change) in x: os.rename(x, x.replace(str_to_change, new_str))

def invert(method):
	'''
	Takes the negation of all scores in score files for the given method.
	This may be useful if you accidentally write a method where highest scores are the best. 
	(evaluate_scores.py will always reward low scores.)
	'''
	for x in os.listdir(os.getcwd()):
		if ('_' + method + '.csv') in x:
			f = pd.read_csv(x, sep='\t', index_col=0)
			for col in f:
				print(col)
				f[col] = [-i for i in f[col]]
			f.to_csv(x, sep='\t')

def remove_old_methods_from_ranking_files():
	'''For all ranking files, removes cols corresponding to methods which do not have score files.'''
	for x in os.listdir(os.getcwd()):
		if 'rankings_from_' in x:
			f = pd.read_csv(x, sep='\t', index_col=0)
			for col in f:
				if not any(col.partition(',')[0] in scores for scores in os.listdir(os.getcwd())):
					f.drop(col, axis=1, inplace=True)
			f.to_csv(x, sep='\t')

def remove_this_method_from_ranking_files(method):
	'''For all ranking files, removes cols corresponding to the specified method.'''
	for x in os.listdir(os.getcwd()):
		if 'rankings_from_' in x:
			f = pd.read_csv(x, sep='\t', index_col=0)
			for col in f:
				if col.partition(',')[0] == method:
					f.drop(col, axis=1, inplace=True)
			f.to_csv(x, sep='\t')

def comma_sep_to_tab_sep(file):
	x = pd.read_csv(file, index_col=0, sep=',', low_memory=False, keep_default_na=False, encoding='Latin-1')
	x.to_csv(file, sep='\t')
	
# os.chdir('results')
# rename('RandomForest.csv','BadRForest.csv')
# rename('RandomForest_rep.csv','RandomForest.csv')
# #invert('InfoGainEntropy')
# #remove_this_method_from_ranking_files('ExtraTrees')
# #rename('Classifier','')
# remove_old_methods_from_ranking_files()

os.chdir('libs')
for x in os.listdir(os.getcwd()):
	if '_fisher_classifier_data' in x:
		print(x)
		comma_sep_to_tab_sep(x)