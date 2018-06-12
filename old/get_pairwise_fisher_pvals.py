import csv
import os
import numpy as np
import pandas as pd
from setup import convert_gmt
from setup import open_csv
import scipy.stats as stats
from get_classifiers import get_classifiers
from joblib import Parallel, delayed
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

def fisher_pairwise_pvals(l_tf, l_tf_genes, f_lib_genes, f_tf_pairs, f_tf_dict, output_fname):

	l = f_lib_genes & l_tf_genes
	#let's call the genes in `l_tf_genes` but not in `f_lib` the residuals of l, `l_res`.
	#we only need its length. 
	l_res = len(l_tf_genes - l)

	if len(l) == 0:
		print('no overlap')
		pair_results = pd.DataFrame(0, index=pd.MultiIndex.from_tuples(f_tf_pairs, names=['tfi','tfj']), 
			columns=['intersection', 'union', 'i minus j', 'j minus i'])
	else:
		#Store the score for each pair in ```pair_results```.
		pair_results = pd.DataFrame(index=pd.MultiIndex.from_tuples(f_tf_pairs, names=['tfi','tfj']), 
			columns=['intersection', 'union', 'i minus j', 'j minus i'])

		notl = f_lib_genes - l
		len_f_df_genes = len(f_lib_genes)
		len_l = len(l)
		len_notl = len(notl)

		#For each feature library tf pair, calculate the info gain resulting from the split.
		#(The split will produce four subsets whose union is the original sample, l_tf_genes).
		for tf_i in set(pair_results.index.get_level_values('tfi')): #iterate over all unique tfs
			#print(tf_i)
			i = f_tf_dict[tf_i][0]
			noti = f_tf_dict[tf_i][1]

			l_and_i = l & i
			notl_and_i = i - l_and_i

			l_and_noti = l - l_and_i
			notl_and_noti = len(noti) - len(l_and_noti)

			### if zeroes, just do noti
			for tf_j in pair_results.loc[tf_i,:].index:
				j = f_tf_dict[tf_j][0]
				notj = f_tf_dict[tf_j][1]

				l_and_i_and_j = len(l_and_i & j)
				l_and_i_and_notj = len(l_and_i) - l_and_i_and_j

				l_and_noti_and_j = len(l_and_noti & j)
				l_and_noti_and_notj = len(l_and_noti) - l_and_noti_and_j

				notl_and_i_and_j = len(notl_and_i & j)
				notl_and_i_and_notj = len(notl_and_i) - notl_and_i_and_j

				notl_and_noti_and_notj = len_f_df_genes - l_and_i_and_notj - l_and_noti_and_notj - notl_and_i_and_notj - len(j)
				notl_and_noti_and_j = notl_and_noti - notl_and_noti_and_notj

				#for debugging
				#print(l_tf, tf_i, tf_j)
				#print(len(l), len(i), len(j), len(i&j))
				#print(l_and_i_and_j, notl_and_i_and_j, len_l - l_and_i_and_j + l_res)
				#print(len_l - l_and_noti_and_notj + l_res, len_notl - notl_and_noti_and_notj, l_and_noti_and_notj)
				#print(l_and_i_and_notj, notl_and_i_and_notj, len_l - l_and_i_and_notj + l_res)
				#print(l_and_noti_and_j, notl_and_noti_and_j, len_l - l_and_noti_and_j + l_res)

				F_i_and_j = stats.fisher_exact([[l_and_i_and_j,notl_and_i_and_j],
					[len_l - l_and_i_and_j + l_res, 20000]], alternative='greater')[1]
				F_i_union_j = stats.fisher_exact([[len_l - l_and_noti_and_notj,
					len_notl - notl_and_noti_and_notj],[l_and_noti_and_notj + l_res,20000]], alternative='greater')[1]
				F_i_minus_j = stats.fisher_exact([[l_and_i_and_notj,notl_and_i_and_notj],[
					len_l - l_and_i_and_notj + l_res, 20000]], alternative='greater')[1] #d
				F_j_minus_i = stats.fisher_exact([[l_and_noti_and_j, notl_and_noti_and_j],[
					len_l - l_and_noti_and_j + l_res, 20000]], alternative='greater')[1] #d

				pair_results.at[(tf_i, tf_j),'intersection'] = F_i_and_j
				pair_results.at[(tf_i, tf_j),'union'] = F_i_union_j
				pair_results.at[(tf_i, tf_j),'i minus j'] = F_i_minus_j
				pair_results.at[(tf_i, tf_j),'j minus i'] = F_j_minus_i

	pair_results.to_csv(output_fname, sep='\t')

def enrichment_wrapper(gmt_pair):
	'''This function is called for each lib pair, and iterates over each method and each tf. 
	pair : dict
		key 'l' contains the label library df, and key 'f' contains the feature library df. 
	'''
	l, f = gmt_pair['l'], gmt_pair['f']
	l_name, f_name = l.index.name, f.index.name
	print('Beginning pairwise p value table creation from', l_name, 'to', f_name)

	#Get the label library experiments whose transcription factors (or drugs) are also in any feature library experiment.
	if 'Drugs' in l_name: overlaps = get_drug_overlaps(gmt_pair)
	else: overlaps = get_overlaps(l.columns.values, f.columns.values)
	print(len(overlaps), 'overlaps found')

	f_lib_genes = set(f.index)
	f_tfs = list(f.columns)
	f_tf_pairs = [(i,j) for i in f_tfs for j in f_tfs if str(i) > str(j)]

	f_tf_dict = {tf:[set(f.index[f[tf]]), set(f.index[~f[tf]])] for tf in f}

	for l_tf in overlaps:
		l_tf_genes = set(l.index[l[l_tf]])
		l_tf_output_fname = 'from_' + l_name + '_to_' + f_name + '_' + l_tf.replace('/', 'fslash').replace('\\', 'bslash') + '_pair_fisher_pvals.csv'
		print(l_tf_output_fname)
		if os.path.isfile(l_tf_output_fname): print('l_tf pairwise fisher p value file already created for', l_tf)
		else:
			fisher_pairwise_pvals(l_tf, l_tf_genes, f_lib_genes, f_tf_pairs, f_tf_dict, l_tf_output_fname)
	return

if __name__ == '__main__':
	all_libs = ['ENCODE_TF_ChIP-seq_2015_abridged', 'ChEA_2016_abridged', 'CREEDS_abridged']

	#Get dataframes of each gmt library in all_libs
	os.chdir('libs')
	all_dfs = {x:convert_gmt('df', x) for x in all_libs}
	os.chdir('..')
	if not os.path.isdir('fisher_pairwise_results'): os.makedirs('fisher_pairwise_results')
	os.chdir('fisher_pairwise_results')

	#Iterate over each gmt pair.
	lib_df_pairs = [{'l':all_dfs[a], 'f':all_dfs[b]} for a in all_libs for b in all_libs if a != b]
	Parallel(n_jobs=1, verbose=0)(delayed(enrichment_wrapper)(pair)for pair in lib_df_pairs)