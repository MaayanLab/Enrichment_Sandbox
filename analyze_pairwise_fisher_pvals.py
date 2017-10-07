import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from get_scores import clean, clean_CREEDS_Drugs

color_dict = {'single':'C0', 'intersection':'C1', 'union':'C2', 'minuend':'C3', 'subtrahend':'C4'}

def to_list(series_or_number):
	if type(series_or_number) == type(pd.Series()): return list(series_or_number)
	else: return [series_or_number]

def single_results(prefix):
	os.chdir('results')
	df = pd.read_csv(prefix + '_Fisher.csv', sep='\t', index_col=0)
	df.index = [clean(tf) for tf in df.index]
	hits, misses = [], []
	for col in df:
		tf = clean(col)
		if tf in df.index: 
			#to_hits = df.loc[tf,col]
			#if type(to_hits) != type(dummy_series): to_hits = [to_hits]
			#else: to_hits = list(to_hits)
			hits += to_list(df.loc[tf,col])
			misses += list(df.drop(tf, axis=0)[col])
		else: misses += list(df[col])
	os.chdir('..')
	return hits, misses

def minus_results(prefix, type):
	hits, misses = [], []
	os.chdir('fisher_pairwise_results')
	for file in os.listdir(os.getcwd()):
		if file.startswith(prefix):
			df = pd.read_csv(file, sep='\t', index_col=False)
			df['tfi'] = df['tfi'].apply(clean)
			df['tfj'] = df['tfj'].apply(clean)
			tf = clean(file.partition(prefix + '_')[2].partition('_pair_fisher_pvals.csv')[0])
			if type == 'minuend':
				hits += to_list(df.loc[df['tfi'] == tf]['i minus j']) + to_list(df.loc[df['tfj'] == tf]['j minus i']) 
				misses += list(df.loc[df['tfi'] != tf]['i minus j']) + list(df.loc[df['tfj'] != tf]['j minus i']) 
			elif type == 'subtrahend':
				hits += to_list(df.loc[df['tfi'] == tf]['j minus i']) + to_list(df.loc[df['tfj'] == tf]['i minus j']) 
				misses += list(df.loc[df['tfi'] != tf]['j minus i']) + list(df.loc[df['tfj'] != tf]['i minus j']) 
			else: raise ValueError('invalid type for minus results')
	os.chdir('..')
	return hits, misses

def regular_results(prefix, result_type):
	hits, misses = [], []
	os.chdir('fisher_pairwise_results')
	for file in os.listdir(os.getcwd()):
		if file.startswith(prefix):
			df = pd.read_csv(file, sep='\t', index_col=False)
			df['tfi'] = df['tfi'].apply(clean)
			df['tfj'] = df['tfj'].apply(clean)
			tf = clean(file.partition(prefix + '_')[2].partition('_pair_fisher_pvals.csv')[0])
			hits += to_list(df.loc[df['tfi'] == tf][result_type]) + to_list(df.loc[df['tfj'] == tf][result_type]) 
			misses += list(df.loc[df['tfi'] != tf][result_type]) + list(df.loc[df['tfj'] != tf][result_type]) 
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

	#Create the grid by iterating over all_libs.
	for i in range(len(all_libs)):
		for j in range(len(all_libs)):
			f_lib = all_libs[i]
			l_lib = all_libs[j]
			subplot = axarr[i,j]
			#Check if you want to plot this pair (e.g. you dont want to if f_lib == l_lib).
			if {'f':f_lib, 'l':l_lib} in lib_pairs:
				prefix = 'from_' + f_lib + '_to_' + l_lib

				if result_type == 'single': hits, misses = single_results(prefix)
				elif result_type == 'minuend': hits, misses = minus_results(prefix, 'minuend')
				elif result_type == 'subtrahend': hits, misses = minus_results(prefix, 'subtrahend')
				else: hits, misses = regular_results(prefix, result_type)

				hit_weights = np.ones_like(hits)/float(len(hits))
				miss_weights = np.ones_like(misses)/float(len(misses))
				subplot.hist(hits, bins=50, range=(0,1), alpha=0.5, color='C0', label=result_type, weights=hit_weights)
				subplot.hist(misses, bins=50, range=(0,1), alpha=0.5, color='C1', label=result_type, weights=miss_weights)

			#Only show y-axis on left-most subplots.
			if j != 0: subplot.yaxis.set_visible(False)
			#Remove box around the bottom-right subplot (legend will go here.)
			if j == i: 
				if j == 0: 
					for side in ['top', 'right', 'bottom', 'left']: subplot.spines[side].set_visible(False)
					subplot.get_xaxis().set_ticks([])
				else: subplot.axis('off')
			#if j == 2 & i == 2: plt.axis('off')
			#if j == 1 & i == 1: plt.axis('off')
			subplot.axes.set_ylim(0,.8)
			#Hide ticks -- although axis='both', this only seems to affect the x-axis.
			#subplot.tick_params(axis='x', which='both', bottom='off', top='off',labelbottom='off')
			#Hide ticks on the y axis.
			subplot.axes.get_yaxis().set_ticks([0, .4, .8])

	#Label the rows and columns of the figure.
	lib_titles = [x.replace('Single_Gene_Perturbations_from_GEO_up', 'CREEDS_sep').replace('ENCODE_TF_ChIP-seq_2015', 'ENCODE').replace('_2016','') for x in all_libs]
	for ax, col in zip(axarr[0], lib_titles): ax.set_title(col)
	for ax, row in zip(axarr[:,0], lib_titles): ax.set_ylabel(row, size='large')
	#Leave some space between the subplots.
	f.subplots_adjust(hspace=.3, wspace=.1)
	#Create a legend in the last cell (should be empty, because it is a diagonal).
	handles = [Rectangle((0,0),1,1,color=c,ec="k", alpha=0.5) for c in ['C0','C1']]
	labels= ['hits','misses']
	plt.legend(handles, labels)
	plt.legend()#[x for sublist in histograms.values for x in sublist], histograms.index, loc='upper left', ncol=1)
	#Title the plot.
	plt.suptitle('p-value histograms for ' + result_type +' method, from [row] to [column]', fontsize=32)
	plt.show()
	return 

all_libs = ['ENCODE_TF_ChIP-seq_2015_abridged', 'ChEA_2016_abridged', 'CREEDS_abridged']
lib_pairs = [{'f':i,'l':j} for i in all_libs for j in all_libs if i != j]
for result_type in ('single','intersection','union','minuend', 'subtrahend'):
	subplots(lib_pairs, all_libs, result_type)