import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from get_scores import clean
from math import sqrt
from sklearn.metrics import auc
from joblib import Parallel, delayed
from collections import Counter
from setup import open_csv

def plot_curve(df, col, prefix):
	'''
	Plots a single curve.
	(Helper function to pairwise_plots and combined_plot.)
	df : pandas.DataFrame
		Columns looks like this: ['Fisher,x', 'Fisher,y', 'Foo3,x', 'Foo3,y']
	col : tuple
		Contains the method name and axis for the column to plot, e.g. ('Fisher','x')
	prefix : str
		Prefix of the file for the column being plotted. Contains the name of the gmt libs. 
	'''
	name = col[0] 
	#Scale the x_vals here
	x_vals = [a/len(df[name + ',x']) for a in df[name + ',x']]
	y_vals = df[name + ',y']
	plt.plot(x_vals, y_vals, label=prefix + ' ' + name 
		+ '    ' + 'AUC: ' + str(np.round(auc(x_vals, y_vals), 4)))

def pairwise_plots(pair):
	'''Plots a graph for the specified library pair.
	pair : dict
		Key 'l' contains the label library name, and key 'f' contains the feature library name. 
	'''

	def get_ranks(file):
		'''Aggregates the ranks of feature lib experiments whose tfs match that which was used to rank it.'''
		ranks_collection = []
		#if 'Drug' in file: scores = pd.read_csv(file, sep='\t', index_col=0, encoding='cp1252')
		scores = open_csv(file)
		#Recall that the index of scores will be the feature library experiments, and the column will be label library tfs. 
		for column in scores: 
			#Sort the feature library experiments by their scores.
			ordered_tfs = scores[column].sort_values().index
			#Get the indices of the feature library experiments which match the label library tf.
			these_ranks = [ordered_tfs.get_loc(x) for x in ordered_tfs if clean(x) == clean(column)]
			ranks_collection += these_ranks
		#Return scores.shape too, which will be needed to make the graph.
		#(scores.shape should be identical between the different methods)
		return ranks_collection, scores.shape

	def get_ranks_up_dn(up_file, dn_file):
		'''Aggregates the ranks of feature lib experiments whose tfs match that which was used to rank it.'''
		ranks_collection = []
		up = open_csv(up_file)
		dn = open_csv(dn_file)
		for column in up:
			#Get the BEST score for each feature library experiment.
			combined = pd.Series([min(up[column][x], dn[column][x]) for x in range(up.shape[0])], index=up.index) 
			#Sort the feature library experiments by their scores.
			ordered_tfs = combined.sort_values().index
			#Get the indices of the feature library experiments which match the label library tf.
			these_ranks = [ordered_tfs.get_loc(x) for x in ordered_tfs if clean(x) == clean(column)]
			ranks_collection += these_ranks
		#Return scores.shape too, which will be needed to make the graph.
		#(scores.shape should be identical between the different methods)
		return ranks_collection, up.shape

	def get_bridge_coords(hit_ranks, ranks_range, n_rankings):
		'''Get the coordinates of the curve which will be plotted.
		hit_ranks : list
			Aggregated ranks of the feature lib experiments whose tfs match that which was used to rank it.
		ranks_range : int
			The number of feature library experiments, i.e. the range of possible ranks.
		n_rankings : int
			The number of label library tfs, i.e. the number of rankings that were created. 
		'''
		down_const = 1/ (ranks_range - 1)
		vert_scale = len(hit_ranks)

		hits = Counter(hit_ranks)
		coords = pd.Series(0.0, index=range(ranks_range))
		coords[0] = hits[0] / vert_scale
		for x in range(1, ranks_range):
			coords[x] = coords[x - 1] + hits[x] / vert_scale - down_const
		return coords.index.values, coords.values, hits, coords

	prefix = 'from_' + pair['l'] + '_to_' + pair['f']
	print(prefix)
	rank_fname = 'rankings_' + prefix + '.csv'
	
	#Load saved plot coordinates, if any.
	if os.path.isfile(rank_fname): agg_c = pd.read_csv(rank_fname, sep='\t', index_col=0)
	else: agg_c = pd.DataFrame()

	#Get the rankings from each enrichment method.
	agg_r = pd.DataFrame()
	for file in os.listdir(os.getcwd()):
		if file.startswith(prefix):
			print('found', file)
			#Get the enrichment method name.
			m_name = str(file).partition(prefix + '_')[2].partition('.csv')[0]
			#Skip if the results are already saved .
			if str(m_name + ',x') in agg_c.columns.values: continue
			#Get and store the ranking results.
			if '_down' in prefix: continue
			elif '_up' in prefix: r, r_shape = get_ranks_up_dn(file, file.replace('up','down'))
			else: r, r_shape = get_ranks(file)
			print(m_name, '   ', len(r), 'hits', '   ', r_shape)
			agg_r[m_name] = r

	#If there were any new ranking files, we need to get and save their plot coordinates for later.
	if not agg_r.empty:
		ranks_range, n_rankings = r_shape[0], r_shape[1]
		for column in agg_r:
			#Get and store the plot coordinates.
			x, y, hits, coords = get_bridge_coords(agg_r[column].values, ranks_range, n_rankings)
			agg_c[column + ',x']=x
			agg_c[column + ',y']=y
	#Save the plot coordinate file.
	agg_c.to_csv(rank_fname, sep='\t')

	#Plot the results for all enrichment methods, if any.
	if not agg_c.empty:
		#BRIDGE PLOT 
		plt.figure(3, figsize=(5.7,5))
		font = {'size': 11}
		plt.rc('font', **font)
		#Plot each enrichment method.
		for column in agg_c:
			col = (column.partition(',')[0], column.partition(',')[2])
			if col[1] == 'x':
				plot_curve(agg_c, col, '')
		plt.title(pair['l'].replace('_up', '_up/dn') + ' to ' + pair['f'].replace('_up', '_up/dn') + ' Bridge Plot')
		plt.xlabel('Rank')
		#Un-comment the line below to view only the first few ranks.
		#plt.gca().set_xlim([0,10])
		plt.legend(prop={'size':9}, frameon=False)
		plt.show()

	return

def combined_plot(lib_df_pairs):
	'''
	Plots a single graph for all results, across all libraries and methods.
	(You can filter out which plots to view.)
	lib_df_pairs : dict
		keys are names of the gmt libraries; values are their dfs. 
	'''
	plt.figure(3, figsize=(12,12))
	font = {'size': 12}
	plt.rc('font', **font)
	for pair in lib_df_pairs:
		prefix = 'from_' + pair['l'] + '_to_' + pair['f']
		rank_fname = 'rankings_' + prefix + '.csv'
		if os.path.isfile(rank_fname): 
			agg_c = open_csv(rank_fname)
			for column in agg_c:
				col = (column.partition(',')[0], column.partition(',')[2])
				#Use the 'if' statement below to filter out which results you want to view.
				if 'Fisher' in col[0] and col[1] == 'x':
					plot_curve(agg_c, col, prefix)
	plt.title(pair['l'].replace('_up', '_up/dn') + ' to ' + pair['f'].replace('_up', '_up/dn') + ' Bridge Plot')
	plt.xlabel('Rank')
	plt.gca().set_xlim([0,.5])
	plt.legend(prop={'size':10}, frameon=False)
	plt.show()
	return

if __name__ == '__main__':
	creeds_libs = ['Single_Gene_Perturbations_from_GEO_up', 'Single_Gene_Perturbations_from_GEO_down']
	other_libs = ['ENCODE_TF_ChIP-seq_2015', 'ChEA_2016']
	all_libs = creeds_libs + other_libs
	lib_pairs = [{'l':a, 'f':b} for a in all_libs for b in all_libs if a != b]
	os.chdir('results')
	x = Parallel(n_jobs=2, verbose=0)(delayed(pairwise_plots)(pair) for pair in lib_pairs)
	combined_plot(lib_pairs)