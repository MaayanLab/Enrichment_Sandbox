import csv
import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import h5py
import math
from joblib import Parallel, delayed
from setup import convert_gmt
import open_csv
import matplotlib as plt

#====================================================================================================================================
'''
The purpose of this script is to collect and visualize the variables which will be used for the modified Fisher's exact test using
	correlation data from ARCHS4. 

p : the Fisher's exact test p value with upper tail
r : the average Pearson's correlation coefficient among non-identical pairs in cell a.
o_frac : the proportion of genes in cell a which are in ARCHS4.
r2 : the average Pearson's correlation coefficient among non-identical pairs between the input gene set and the gene set of interest.
	(in other words, non-identical pairs between (a|b) and (a|c).)
o_frac2 : the proportion of genes in cells a, b and c which are in ARCHS4.
r3 : the average Pearson's correlation coefficient among pairs between cells b and c.
o_frac3 : the proportion of genes in cells b and c which are in ARCHS4.
'''
#====================================================================================================================================

def clean(tf):
	'''Extracts the transcription factor name from the name of a gmt experiment.'''
	return str(tf).partition('_')[0].partition(' ')[0].upper()

def get_overlaps(label_tfs, feature_tfs):
	'''Returns a list of label library experiments whose transcription factors are also in any feature library experiment.
	label_tfs: list-like 
	feature_tfs: list-like 
	'''
	cleaned_overlaps = {clean(x) for x in feature_tfs} & {clean(x) for x in label_tfs}
	l_in_overlaps = [x for x in label_tfs if clean(x) in cleaned_overlaps]
	return l_in_overlaps

def get_fav_vars(l_tf, l_tf_genes, f_matrix, l_lib, f_lib, ARCHS4_genes_dict):
	'''For each tf pair, saves the variables for the adjusted Fisher's test to a csv.'''
	cwd = os.getcwd()
	os.chdir('..')
	os.chdir('libs')
	ARCHS4 = h5py.File(l_lib + '_ARCHS4_corr.h5', 'r+')
	os.chdir('..')
	os.chdir(cwd)

	#Get variables which will be the same for each iteration, or only depend on the organism. 
	l_tf_genes = set(l_tf_genes)
	c_overlap_dict = {'human':{str(x).encode('utf-8') for x in l_tf_genes} & set(ARCHS4_genes_dict['human'].index),
		'mouse':{str(x).encode('utf-8') for x in l_tf_genes} & set(ARCHS4_genes_dict['mouse'].index)}
	if f_lib == 'ENCODE_TF_ChIP-seq_2015': organism_dict = {'hg19': 'human', 'mm9': 'mouse'}

	#For each tf, we will store the information collected in 'info'.
	info = pd.DataFrame(index=['p', 'o_frac', 'r', 'o_frac2', 'r2', 'o_frac3', 'r3'], columns = f_matrix.columns)
	#Default value for overlap fraction is zero. 
	info.loc['o_frac',:] = 0
	info.loc['o_frac2',:] = 0
	info.loc['o_frac3',:] = 0

	#Iterate over each f library tf. 
	for tf in list(f_matrix.columns):
		#Get the regular p val, just as we do in Fisher().
		f_tf_genes = set(f_matrix.index[f_matrix[tf]])
		#'a_genes' are the genes in both the feature library tf and the label library tf.
		#In other words, 'a_genes' is the intersection cell of the 2x2 contingency table. 
		a_genes = f_tf_genes & l_tf_genes
		a = len(a_genes)
		b = len(f_tf_genes) - a
		c = len(l_tf_genes) - a
		d = 20000 - a - b - c
		info.at['p', tf] = stats.fisher_exact([[a,b],[c,d]], alternative='greater')[1]

		#Determine which organism this tf data came from. Doing this depends on the gmt file. 
		if f_lib == 'CREEDS': organism = tf.partition(' GSE')[0].rpartition(' ')[2]
		elif f_lib == 'ChEA_2016': 
			organism = tf.rpartition('_')[2].lower()
			if organism in ['ovary', 'hela', 'neurons', 'gbm']: organism = 'human'
		elif f_lib == 'ENCODE_TF_ChIP-seq_2015': organism = organism_dict[tf.rpartition('_')[2].lower()]
		else: print('invalid lib name!')
		if organism == 'rat': organism = 'mouse'
		if organism in ['human', 'mouse']:

			#Use 'ARCHS4_genes_dict' to get the appropriate list of ARCHS4 genes.
			ARCHS4_genes = ARCHS4_genes_dict[organism]

			#b_overlap will be the genes in f_tf_genes which are also in ARCHS4: (a|b) & ARCHS4_genes.
			#c_overlap will be the genes in l_tf_genes which are also in ARCHS4: (a|c) & ARCHS4_genes.
			#a_overlap will be the genes in f_tf_genes, l_tf_genes and ARCHS4: a & ARCHS4_genes.
			b_overlap = {str(x).encode('utf-8') for x in f_tf_genes} & set(ARCHS4_genes.index)
			c_overlap = c_overlap_dict[organism]
			a_overlap = b_overlap & c_overlap

			#only_b: b & ARCHS4_genes.
			#only_c: c & ARCHS4_genes.
			only_b = b_overlap - a_overlap
			only_c = c_overlap - a_overlap

			a_l_o = len(a_overlap)
			b_l_o = len(b_overlap)
			c_l_o = len(c_overlap)
			only_b_l_o = len(only_b)
			only_c_l_o = len(only_c)
			
			if b_l_o > 0 or c_l_o > 0: 
				info.at['o_frac',tf] = (a_l_o + only_b_l_o + only_c_l_o) / (a+b+c)
				if b_l_o > 0 and c_l_o > 0:
					#Create series to use to index ARCHS4
					b_list = list(b_overlap)
					b_series = pd.Series(ARCHS4_genes[b_list], index=b_list)
					b_series.sort_values(inplace=True)
					c_list = list(c_overlap)
					c_series = pd.Series(ARCHS4_genes[c_list], index=c_list)
					c_series.sort_values(inplace=True)

					r_vals_abc = pd.DataFrame(ARCHS4[organism]['data']['correlation'][b_series.values,:][:,c_series.values], index=b_series.index, columns=c_series.index)
					#Get the average r value for pairs, not including matches.
					info.at['r',tf] = (r_vals_abc.values.sum() - a_l_o)/(b_l_o*c_l_o - a_l_o)

					if a_l_o > 1:
						info.at['o_frac2',tf] = a_l_o / (a)
						r_vals_a = r_vals_abc.loc[a_overlap, a_overlap]
						#Get the average r value for non-diagonal i.e. pairwise entries. 
						#(Each pair is duplicated across the diagonal, but this does not affect the result.)
						info.at['r2',tf] = (r_vals_a.values.sum() - a_l_o)/(a_l_o*a_l_o - a_l_o)

					if only_c_l_o > 0 and only_b_l_o > 0:
						info.at['o_frac3',tf] = (only_b_l_o + only_c_l_o) / (b+c)
						r_vals_bc_only = r_vals_abc.loc[only_b, only_c]
						#Get the average r value for all pairs. Note there will be no matches between b and c.
						info.at['r3',tf] = (r_vals_bc_only.values.sum())/(only_b_l_o*only_c_l_o)
		#If the organism cannot be identified, ARCHS4 cannot be used.
		else: print('weird organism:', organism, tf)
	ARCHS4.close()
	info.to_csv('from_' + l_lib + '_to_' + f_lib + '_' + l_tf.replace('/','slash') + '_fisher_adjusted_vars.csv', sep='\t')

def fav_wrapper(pair):
	'''This function is called for each lib pair, and gets the fisher adjusted variables.
	pair : dict
		key 'l' contains the label library df, and key 'f' contains the feature library df. 
	'''
	l_name, f_name = pair['l'].index.name, pair['f'].index.name
	print('Getting fisher adjusted variables from', l_name, 'to', f_name)

	os.chdir('..')
	os.chdir('libs')
	ARCHS4 = h5py.File(l_name + '_ARCHS4_corr.h5', 'r+')
	os.chdir('..')
	os.chdir('fisher_ARCHS4_values')
	h_genes = ARCHS4['human']['meta']['genes']
	m_genes = ARCHS4['mouse']['meta']['genes']
	ARCHS4_genes_dict = {'human': pd.Series(range(len(h_genes)), index=h_genes[...]),
		'mouse': pd.Series(range(len(m_genes)), index=m_genes[...])}
	ARCHS4.close()

	#Get the label library experiments whose transcription factors are also in any feature library experiment.
	overlaps = get_overlaps(pair['l'].columns.values, pair['f'].columns.values)

	#Iterate over each tf in the overlaps.
	for l_tf in overlaps:
		print('fisher_ARCHS4_values', l_tf) #for diagnostics
		l_tf_genes = pair['l'].index[pair['l'][l_tf]]
		get_fav_vars(l_tf, l_tf_genes, pair['f'], l_name, f_name, ARCHS4_genes_dict)
	return

def aggregate_fav_vars(lib_pairs):
	'''For each library pair, this function looks at each Fisher adjusted variable file, and 
		appends each of its columns to either the hit file or the miss file.'''
	for pair in lib_pairs:
		prefix = 'from_' + pair['l'] + '_to_' + pair['f'] + '_'
		hit_output_fname = prefix + 'ARCHS4_vars_hits.csv'
		hit_df = pd.DataFrame(index=['p','o_frac','r','o_frac2','r2','o_frac3','r3'])
		miss_output_fname = prefix + 'ARCHS4_vars_misses.csv'
		miss_df = pd.DataFrame(index=['p','o_frac','r','o_frac2','r2','o_frac3','r3'])

		for file in os.listdir(os.getcwd()):
			if file.startswith(prefix):
				df = open_csv(file)
				file_tf = clean(file.partition(prefix)[2].partition('_fisher_adjusted_vars.csv')[0])
				for col in df:
					if clean(col) == file_tf:
						hit_df[file_tf + ', to ' + col] = df[col]
						print('match:	', col, file_tf)
					else:
						miss_df[file_tf + ', to ' + col] = df[col]
		hit_df.to_csv(hit_output_fname)
		miss_df.to_csv(miss_output_fname)
	return

def compare_vars(lib_pairs, all_libs):
	'''This function plots comparisons of variables between the hit file and miss file.'''
	for x in [['p'],['o_frac','o_frac2','o_frac3'],['r','r2','r3']]:
		f, axarr = plt.subplots(nrows=len(all_libs),ncols=len(all_libs), figsize=(8,7))
		font = {'size':11}
		plt.rc('font', **font)
		methods = pd.Series()

		#Create the grid by iterating over all_libs.
		for i in range(len(all_libs)):
			for j in range(len(all_libs)):
				rlib = all_libs[i]
				clib = all_libs[j]
				if rlib != clib:
					subplot = axarr[i,j]
					prefix = 'from_' + rlib + '_to_' + clib
					hits = open_csv(prefix + '_ARCHS4_vars_hits.csv')
					misses = open_csv(prefix + '_ARCHS4_vars_misses.csv')
					it = 0
					for var in x:
						subplot.hist(hits.loc[var,:], bins=50, alpha=.6, color='C'+it, label=var + ' hits', hatch='//')
						subplot.hist(misses.loc[var,:], bins=50, alpha=.4, color='C'+it, label=var + ' misses', hatch='\\')
						it += 1
					#Only show y-axis on left-most plots.
					if j != 0: subplot.yaxis.set_visible(False)
					#Do not show ticks -- although axis='both', this only seems to affect the x-axis.
					subplot.tick_params(axis='both', which='both', bottom='off', top='off',labelbottom='off')
					subplot.axes.get_yaxis().set_ticks([])

		#Label the rows and columns of the figure.
		lib_titles = [x.replace('Single_Gene_Perturbations_from_GEO_up', 'CREEDS_up/dn_sep') for x in all_libs]
		for ax, col in zip(axarr[0], lib_titles): ax.set_title(col)
		for ax, row in zip(axarr[:,0], lib_titles): ax.set_ylabel(row, size='large')
		#Leave some space between the subplots.
		f.subplots_adjust(hspace=.15, wspace=.1)
		#Create a legend in the last cell (should be empty, as it is a diagonal).
		plt.legend([x for sublist in methods.values for x in sublist], methods.index)
		plt.suptitle('Histograms', fontsize=15)
		plt.show()
	return

if __name__ == '__main__':
	all_libs = ['CREEDS', 'ENCODE_TF_ChIP-seq_2015', 'ChEA_2016']
	lib_pairs = [{'l':a, 'f':b} for a in all_libs for b in all_libs if a != b]
	lib_df_pairs = [{'l':all_dfs[a], 'f':all_dfs[b]} for a in all_libs for b in all_libs if a != b]

	#Get dataframes of each gmt library in all_libs.
	os.chdir('libs')
	all_dfs = {x:convert_gmt('df', x) for x in all_libs}
	os.chdir('..')
	if not os.path.isdir('fisher_ARCHS4_variables'): os.makedirs('fisher_ARCHS4_variables')
	os.chdir('fisher_ARCHS4_variables')

	#Get files for the fisher adjusted variables for each pair for each overlapping label library tf. 
	Parallel(n_jobs=6, verbose=0)(delayed(fav_wrapper)(pair)for pair in lib_df_pairs)

	#For each pair, aggregate the files into a file for hits and a file for misses.
	aggregate_fav_vars(lib_pairs)

	#Analyze variables for the hits vs the misses.
	compare_vars(lib_pairs, all_libs)

