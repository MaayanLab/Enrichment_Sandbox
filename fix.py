import os
import pandas as pd

#========================================================================================================
'''
The purpose of this script is to adjust the score files for a certain method.
'''
#========================================================================================================

def rename(str_to_change, new_str):
	for x in os.listdir(os.getcwd()):
		if (str_to_change) in x: os.rename(x, x.replace(str_to_change, new_str))


def invert(method):
	'''
	This may be useful if you write a method where highest scores are the best. 
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
	for x in os.listdir(os.getcwd()):
		if 'rankings_from_' in x:
			f = pd.read_csv(x, sep='\t', index_col=0)
			for col in f:
				if not any(col.partition(',')[0] in scores for scores in os.listdir(os.getcwd())):
					f.drop(col, axis=1, inplace=True)
			f.to_csv(x, sep='\t')

def remove_this_method_from_ranking_files(method):
	for x in os.listdir(os.getcwd()):
		if 'rankings_from_' in x:
			f = pd.read_csv(x, sep='\t', index_col=0)
			for col in f:
				if col.partition(',')[0] == method:
					f.drop(col, axis=1, inplace=True)
			f.to_csv(x, sep='\t')

os.chdir('results')

#invert('InfoGainEntropy')
remove_this_method_from_ranking_files('ExtraTrees')
#rename('Classifier','')
#remove_old_methods_from_ranking_files()