import csv
import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import gsea_methods as m
from setup import convert_gmt
import h5py

def clean(tf):
	'''Input: a list or set of transcription factor names. Output: a set of the 'cleaned', i.e. simplified tf names.'''
	return tf.partition('_')[0].partition(' ')[0].upper()

def get_overlaps(label_tfs, feature_tfs):
	'''Input: list with feature_lib and label_lib names; list or set of feature library tf names; list or set of label library tf names.
	Output: list of tfs in label library if the clean version of that tf is also found in the feature library'''
	cleaned_overlaps = {clean(x) for x in feature_tfs} & {clean(x) for x in label_tfs}
	l_in_overlaps = [x for x in label_tfs if clean(x) in cleaned_overlaps]
	return l_in_overlaps

def get_methods(l,f, l_name, f_name):
	'''This is where you choose which gsea methods to run, and which paramaters to use!
	Return a dataframe indexed by ['func', 'params'], and columns as the name of each gsea method.
	For example, a column with name 'Forest3' might have 'func' Forest and 'params' bootstrap = False.
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
	os.chdir('..')
	os.chdir('libs')
	#Get the correlation data
	ARCHS4 = h5py.File(l_name + '_ARCHS4_corr.h5', 'r+')
	os.chdir('..')
	os.chdir('results')
	h_genes = ARCHS4['human']['meta']['genes']
	m_genes = ARCHS4['mouse']['meta']['genes']
	ARCHS4_genes_dict = {'human': pd.Series(np.arange(len(h_genes)), index=h_genes[...]),
		'mouse': pd.Series(np.arange(len(m_genes)), index=m_genes[...])}
	ARCHS4.close()

	df = pd.DataFrame(index=['func', 'params'])
	#df['Control'] = [m.Control, ([f.columns.values])]
	#df['Fisher'] = [m.Fisher, ([f])]
	df['FAV'] = [m.FisherAdjusted, (f, l_name, f_name, ARCHS4_genes_dict)]
	#df['RandomForest'] = [m.Forest, (f, f.columns.values, 73017, None, True, None, None)]

	return df

def enrichment_wrapper(pair, l_name, f_name):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : dict
		key 'l' contains the label library df, and key 'f' contains the feature library df. 
	'''
	output_heading = 'from_' + l_name + '_to_' + f_name
	print('Beginning enrichment analysis. Output file names will begin with', output_heading)

	#Get the intersection of tfs in the label library and tfs in the feature library.
	overlaps = get_overlaps(pair['l'].columns.values, pair['f'].columns.values)

	#Get the methods and parameters with which to perform GSEA. 
	method_params = get_methods(pair['l'], pair['f'], l_name, f_name)
	#ARCHS4

	#Iterate over each method.
	for m_name in method_params:
		#Special case for ZAndCombined.
		if m_name == 'ZAndCombined': 
			output_fnames = (output_heading + '_Z.csv', output_heading + '_Combined.csv')
		elif m_name == 'FAV':
			output_fnames = [output_heading + '_FisherAdjusted' + str(x) + '.csv' for x in range(2,11)]
		else: 
			output_fnames = (output_heading + '_' + m_name + '.csv',)
		#Check if the file has already been created
		if os.path.isfile(output_fnames[0]): print('ranking file already created for', m_name)
		else:
			#Use dataframes to store rankings from each tf iteration.
			dfs = {x:pd.DataFrame() for x in range(2,11)}
			#Iterate over each tf in the overlaps.
			for l_tf in list(overlaps):
				print(m_name, l_tf) #for diagnostics
				l_tf_genes = set(pair['l'].index[pair['l'][l_tf]])
				result = method_params.at['func', m_name](l_tf_genes, *method_params.at['params', m_name])
				if len(dfs) == 1: dfs[0][l_tf] = result
				else: 
					for x in dfs: 
						dfs[x][l_tf] = result[x-2]
						#Send the results to the csv.
			for x in dfs: dfs[x].to_csv(output_fnames[x-2], sep='\t')
	return

if __name__ == '__main__':
	all_libs = ['CREEDS', 'ENCODE_TF_ChIP-seq_2015', 'ChEA_2016']

	os.chdir('libs')
	#Get dataframes of each gmt library in all_libs
	all_dfs = {x:convert_gmt('df', x) for x in all_libs}
	os.chdir('..')
	if not os.path.isdir('results'): os.makedirs('results')
	os.chdir('results')

	#Iterate over each gmt pair.
	#lib_df_pairs = [{'l':all_dfs[a], 'f':all_dfs[a]} for a in all_dfs if a != 'CREEDS']
	lib_df_pairs = [{'l':all_dfs[a], 'f':all_dfs[b]} for a in all_dfs for b in all_dfs if a.partition('_')[0] != b.partition('_')[0]]
	Parallel(n_jobs=6, verbose=0)(delayed(enrichment_wrapper)(pair, pair['l'].index.name, pair['f'].index.name) for pair in lib_df_pairs)