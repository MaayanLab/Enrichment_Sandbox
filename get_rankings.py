import csv
import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import gsea_methods as m
from setup import convert_gmt

def clean(tf):
	'''Input: a list or set of transcription factor names. Output: a set of the 'cleaned', i.e. simplified tf names.'''
	return str(tf).partition('_')[0].partition(' ')[0].upper()

def get_overlaps(label_tfs, feature_tfs):
	'''Input: list with feature_lib and label_lib names; list or set of feature library tf names; list or set of label library tf names.
	Output: list of tfs in label library if the clean version of that tf is also found in the feature library'''
	cleaned_overlaps = {clean(x) for x in feature_tfs} & {clean(x) for x in label_tfs}
	l_in_overlaps = [x for x in label_tfs if clean(x) in cleaned_overlaps]
	return l_in_overlaps

def get_methods(l,f, l_name, f_name):
	'''This is where you choose which gsea methods to run, and which paramaters to use!
	Return a dataframe indexed by ['method', 'params'], and columns as the name of each gsea method.
	For example, a column with name 'Forest3' might have 'method' Forest and 'params' bootstrap = False.
	See gsea_methods for the available methods and necessary params.
	You must specify ALL the params EXCEPT for l_tf_genes, which is defined and called later, in enrichment_wrapper.
	l : pd.DataFrame
		the "label" gmt library, from which tfs are being used as input
	f: pd.DataFrame
		the "feature" gmt library, whose tfs are being ranked by the gsea method
	l_name : str
		the name of l, for example "ChEA_2016"
	f_name: str
		the name of f, for example "CREEDS"
	'''
	df = pd.DataFrame(index=['method', 'params'])
	df['Control'] = [m.Control, ([f.columns.values])]
	df['Fisher'] = [m.Fisher, ([f])]
	df['FisherAdjusted'] = [m.FisherAdjusted, (f, l_name, f_name)]
	df['ZAndCombined'] = [m.ZAndCombined, (f_name, l_name, f.columns.values)]

	train_group, features,  = f, f.columns.values, 
	random_state = 70317 #int(time.strftime('%m_%d_%y'))

	Tree = {m.Forest, m.ForestDrop}
	for x in Tree:
		n = 1
		for mf in ['auto', 'log2']:
			for bs in [True, False]:
				for cw in [None, 'balanced']:
					if x == m.ForestDrop: maxdepths = [7,5]
					else: maxdepths = [None, 7]
					for md in maxdepths:
						xname = x.__name__ + str(n)
						df[xname] = [x, (train_group, features, random_state, mf, bs, cw, md)]
						n += 1

	ML = {m.ExtraTreesClassifier, m.RandomTreesEmbedding, m.AdaBoostClassifier, 
		m.GradientBoostingClassifier}
	for x in ML:
		xname = x.__name__.partition('Classifier')[0]
		df[xname] = [m.ML_wrapper, (x, train_group, features, random_state)]

	return df

def enrichment_wrapper(pair):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : dict
		key 'l' contains the label library df, and key 'f' contains the feature library df. 
	'''
	output_dir = 'from_' + pair['l'].name + '_to_' + pair['f'].name
	print('Beginning enrichment analysis. Storing results in ' + output_dir)

	#Get the intersection of tfs in the label library and tfs in the feature library
	overlaps = get_overlaps(pair['l'].columns.values, pair['f'].columns.values)

	#Get the methods and parameters with which to perform GSEA. 
	method_params = get_methods(pair['l'], pair['f'], pair['l'].name, pair['f'].name)
	#ARCHS4

	#Iterate over each method
	for m_name in method_params:
		#Special case for ZAndCombined
		if m_name == 'ZAndCombined': 
			output_fnames = (output_dir + '_Z.csv', output_dir + '_Combined.csv')
			#dfs will store rankings from each iteration over the overlaps
			dfs = (pd.DataFrame(), pd.DataFrame())
		else: 
			output_fnames = (output_dir + '_' + m_name + '.csv',)
			#dfs will store rankings from each iteration over the overlaps
			dfs = (pd.DataFrame(),)
		#Check if the file has already been created
		if os.path.isfile(output_fnames[0]): print('ranking file already created for', m_name)
		else:
			#Iterate over each tf in the overlaps
			for l_tf in list(overlaps):
				#print(m_name, l_tf) #for diagnostics
				l_tf_genes = set(l[l_tf].dropna().index)
				result = mp['method'](l_tf_genes, *mp['params'])
				for x in range(len(dfs)): dfs[x][l_tf] = result[x]
			#Send the results to the csv
			for x in range(len(dfs)): dfs[x].to_csv(output_fnames[x], sep='\t')
			print('done ' + m_name)
	return

if __name__ == '__main__':
	all_libs = ['CREEDS', 'ChEA_2016', 'ENCODE_TF_ChIP-seq_2015', 'ENCODE_2017']

	os.chdir('libs')
	#Get dataframes of each gmt library in all_libs
	all_dfs = {x:convert_gmt('df', x) for x in all_libs}
	os.chdir('..')
	os.chdir('results')

	#Iterate over each gmt pair.
	lib_df_pairs = [{'l':all_dfs[a], 'f':all_dfs[b]} for a in all_dfs for b in all_dfs if a != b]
	Parallel(n_jobs=6, verbose=0)(delayed(enrichment_wrapper)(pair) for pair in lib_df_pairs)