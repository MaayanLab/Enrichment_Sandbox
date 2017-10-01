import os
import csv
import pandas as pd
import numpy as np
from get_scores import clean, clean_CREEDS_Drugs

def single_results(prefix):
	os.chdir('results')
	df = pd.read_csv(prefix + '_Fisher.csv', sep='\t')
	df.index = [clean(tf) for tf in df.index]
	for col in df:
		tf = clean(col)
		hits += list(df[col,tf])
		misses += list(df.drop(tf, axis=0)[col]) 
	os.chdir('..')
	return hits, misses

def minus_results(prefix, type):
	hits, misses = [], []
	os.chdir('fisher_pairwise_results')
	for file in os.listdir(os.getcwd()):
		if file.startswith(prefix):
			df = pd.open_csv(file, sep='\t', index_col=False)
			df.index = [clean(x) for x in df.index]
			tf = clean(file.partition(prefix)[2].partition('_pair_fisher_pvals.csv')[0])
			if type == 'minuend':
				hits += list(df.loc[df['i'] == tf]['i minus j']) + list(df.loc[df['j'] == tf]['j minus i']) 
				misses += list(df.loc[df['i'] != tf]['i minus j']) + list(df.loc[df['j'] != tf]['j minus i']) 
			elif type == 'subtrahend':
				hits += list(df.loc[df['i'] == tf]['j minus i']) + list(df.loc[df['j'] == tf]['i minus j']) 
				misses += list(df.loc[df['i'] != tf]['j minus i']) + list(df.loc[df['j'] != tf]['i minus j']) 
			else: raise ValueError('invalid type for minus results')
	os.chdir('..')
	return hits, misses

def regular_results(prefix, result_type):
	hits, misses = [], []
	os.chdir('fisher_pairwise_results')
	for file in os.listdir(os.getcwd()):
		if file.startswith(prefix):
			df = pd.open_csv(file, sep='\t', index_col=False)
			df.index = [clean(x) for x in df.index]
			tf = clean(file.partition(prefix)[2].partition('_pair_fisher_pvals.csv')[0])
			hits += list(df.loc[df['i'] == tf][result_type]) + list(df.loc[df['j'] == tf][result_type]) 
			misses += list(df.loc[df['i'] != tf][result_type]) + list(df.loc[df['j'] != tf][result_type]) 
	os.chdir('..')
	return hits, misses

def subplots(lib_pairs, all_libs, result_type):
	'''
	Like pairwise_plots, but shows all plots as different subplots in the same figure, organized in a grid.
	lib_pairs : dict
		keys are names of the gmt libraries; values are their dfs. 
	all_libs : list-like
		the gmt libraries in lib_pairs
	'''
	f, axarr = plt.subplots(nrows=len(all_libs),ncols=len(all_libs), figsize=(15,15))
	font = {'size':20}
	plt.rc('font', **font)

	#Collect all the methods found so that we can create a legend at the end.
	#IMPORTANT: this only works if each lib_pair has the EXACT same plots, e.g. Fisher and Control.
	methods = pd.Series()

	#Create the grid by iterating over all_libs.
	for i in range(len(all_libs)):
		for j in range(len(all_libs)):
			f_lib = all_libs[i]
			l_lib = all_libs[j]
			subplot = axarr[i,j]
			#Check if you want to plot this pair (e.g. you dont want to if f_lib == l_lib).
			if {'f':f_lib, 'l':l_lib} in lib_pairs:
				hits, misses = [], []
				prefix = 'from_' + f_lib + '_to_' + l_lib

				if result_type == 'single': hits, misses = single_results(prefix)
				elif result_type == 'minuend': hits, misses = minus_results(prefix, 'minuend')
				elif result_type == 'subtrahend': hits, misses = minus_results(prefix, 'subtrahend')
				else: hits, misses = regular_results(prefix, result_type)

				####left off here
				if name in color_dict: 
					methods[name] = subplot.plot(x_vals, y_vals, label=name, color=color_dict[name], linewidth=linewidth)
				else:
					methods[name] = subplot.plot(x_vals, y_vals, label=name, linewidth=linewidth)
					#If you want to view legends for each subplot (e.g. to see the AUC), you will need to un-comment this line.
					#subplot.legend(fontsize=12)


			#Uncomment below to scale all subplots equally (to compare relative sizes between subplots).
			subplot.set_ylim([-.1,.5])
			if top_10: subplot.set_xlim([0,.10])
			#Only show y-axis on left-most subplots.
			if j != 0: subplot.yaxis.set_visible(False)
			#Remove box around the bottom-right subplot (legend will go here.)
			if j == 2 and i == 2: plt.axis('off')
			#Hide ticks -- although axis='both', this only seems to affect the x-axis.
			subplot.tick_params(axis='x', which='both', bottom='off', top='off',labelbottom='off')
			#Hide ticks on the y axis.
			#subplot.axes.get_yaxis().set_ticks([])

	#Label the rows and columns of the figure.
	print(all_libs)
	lib_titles = [x.replace('Single_Gene_Perturbations_from_GEO_up', 'CREEDS_sep').replace('ENCODE_TF_ChIP-seq_2015', 'ENCODE').replace('_2016','') for x in all_libs]
	for ax, col in zip(axarr[0], lib_titles): ax.set_title(col)
	for ax, row in zip(axarr[:,0], lib_titles): ax.set_ylabel(row, size='large')
	#Leave some space between the subplots.
	f.subplots_adjust(hspace=.15, wspace=.1)
	#Create a legend in the last cell (should be empty, because it is a diagonal).
	plt.legend([x for sublist in methods.values for x in sublist], methods.index, loc='upper left', ncol=1)
	#Title the plot.
	if top_10: plt.suptitle('Bridge plots from [row] to [column], top 10 percentile of ranks', fontsize=28)
	else: plt.suptitle('Bridge plots from [row] to [column]', fontsize=32)
	plt.show()
	return 


def minus_results()

all_libs = ['ENCODE_TF_ChIP-seq_2015_abridged', 'ChEA_2016_abridged', 'CREEDS_abridged']
lib_pairs = [{'f':i,'l':j} for i in libs for j in libs if i != j]
for result_type in ('single','intersection','union','minuend', 'subtrahend'):
	subplots(lib_pairs, all_libs, result_type)