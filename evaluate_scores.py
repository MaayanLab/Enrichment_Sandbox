import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from get_scores import clean
from math import sqrt, ceil
from sklearn.metrics import auc
from joblib import Parallel, delayed
from collections import Counter
from setup import open_csv

def plot_curve(df, col, prefix):
	'''
	Plots a single bridge plot curve.
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
	'''Creates a bridge plot for enrichment between the specified library pair.
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
			combined = pd.Series([min(*l) for l in zip(up[column], dn[column])], index=up.index) 
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

	#If any new ranking files were found, we need to get and save their plot coordinates for later.
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
		plt.figure(1, figsize=(10,10))
		font = {'size': 11}
		plt.rc('font', **font)
		#Plot each enrichment method.
		for column in agg_c:
			col = (column.partition(',')[0], column.partition(',')[2])
			#Filter for only certain enrichment methods here using the if statement.
			if col[1] == 'x' and (col[0] in ['Fisher','CombinedFF', 'RandomForest', 'CombinedFF2']):
				plot_curve(agg_c, col, '')
		plt.title(pair['l'].replace('_up', '_up/dn') + ' to ' + pair['f'].replace('_up', '_up/dn') + ' Bridge Plot')
		plt.xlabel('Rank')
		#Uncomment the line below to view only the first few ranks.
		#plt.gca().set_xlim([0,.10])
		plt.legend(prop={'size':9}, frameon=False)
		plt.show()

	return

def combined_plot(lib_df_pairs):
	'''
	Plots a single graph for all results, across all libraries and methods.
	lib_df_pairs : dict
		keys are names of the gmt libraries; values are their dfs. 
	'''
	plt.figure(2, figsize=(12,12))
	font = {'size': 12}
	plt.rc('font', **font)
	for pair in lib_df_pairs:
		prefix = 'from_' + pair['l'] + '_to_' + pair['f']
		rank_fname = 'rankings_' + prefix + '.csv'
		if os.path.isfile(rank_fname): 
			agg_c = open_csv(rank_fname)
			for column in agg_c:
				col = (column.partition(',')[0], column.partition(',')[2])
				#Filter for only certain enrichment methods here using the if statement.
				if col[1] == 'x' and (col[0] in ['Fisher','RandomForest']):
					plot_curve(agg_c, col, prefix)
	plt.title(pair['l'].replace('_up', '_up/dn') + ' to ' + pair['f'].replace('_up', '_up/dn') + ' Bridge Plot')
	plt.xlabel('Rank')
	#Uncomment the line below to view only the first few ranks.
	#plt.gca().set_xlim([0,.10])
	plt.legend(prop={'size':10}, frameon=False)
	plt.show()
	return

def subplots(lib_pairs, all_libs):
	'''
	Like pairwise_plots, but shows all plots as different subplots in the same figure, organized in a grid.
	lib_pairs : dict
		keys are names of the gmt libraries; values are their dfs. 
	all_libs : list-like
		the gmt libraries in lib_pairs
	'''
	f, axarr = plt.subplots(nrows=len(all_libs),ncols=len(all_libs), figsize=(8,7))
	font = {'size':11}
	plt.rc('font', **font)

	#IMPORTANT: this only works if each lib_pair has the EXACT same plots, e.g. Fisher and Control.
	methods = pd.Series()

	#Create the grid by iterating over all_libs.
	for i in range(len(all_libs)):
		for j in range(len(all_libs)):
			rlib = all_libs[i]
			clib = all_libs[j]
			subplot = axarr[i,j]
			#Check if you want to plot this pair (e.g. you dont if rlib and clib are identical).
			if {'l':rlib, 'f':clib} in lib_pairs:
				prefix = 'from_' + rlib + '_to_' + clib
				rank_fname = 'rankings_' + prefix + '.csv'
				if os.path.isfile(rank_fname): 
					agg_c = open_csv(rank_fname)
					for column in agg_c:
						col = (column.partition(',')[0], column.partition(',')[2])
						#Use the 'if' statement below to filter out which results you want to view.
						if col[1] == 'x' and (col[0] in ['Fisher','RandomForest', 'Control', 'CombinedFF']):
							name = col[0] 
							x_vals = [a/len(agg_c[name + ',x']) for a in agg_c[name + ',x']]
							y_vals = agg_c[name + ',y']
							methods[name] = subplot.plot(x_vals, y_vals, label= name)
			subplot.set_ylim([-.1,.4])
			#uncomment below to see just top 10%.
			#subplot.set_xlim([0,.1])
			#Only show y-axis on left-most plots.
			if j != 0: subplot.yaxis.set_visible(False)
			#Do not show ticks -- although axis='both', this only seems to affect the x-axis.
			subplot.tick_params(axis='both', which='both', bottom='off', top='off',labelbottom='off')

	#Label the rows and columns of the figure
	lib_titles = [x.replace('Single_Gene_Perturbations_from_GEO_up', 'CREEDS_up/dn_sep') for x in all_libs]
	for ax, col in zip(axarr[0], lib_titles): ax.set_title(col)
	for ax, row in zip(axarr[:,0], lib_titles): ax.set_ylabel(row, size='large')
	#Leave some space between the subplots
	f.subplots_adjust(hspace=.15, wspace=.1)
	#Create a legend in the last cell (should be empty, as it is a diagonal).
	plt.legend([x for sublist in methods.values for x in sublist], methods.index)
	plt.suptitle('Bridge Plots From "Row" to "Column", Top 10 Percentile', fontsize=15)
	plt.show()
	return methods

def hexbin_method_comparison(libs, m1, m2):
	'''
	Creates a grid of heatmaps comparing how two different methods rank the "hits"/"matches".
	libs : list-like
		the gmt libraries to create heatmaps for
	m1 : str
		the first method, to correspond with the x-axis
	m2 : str
		the second method, to correspond with the y-axis
	'''

	def get_coords(f1, f2):
		coords_collection = []
		s1, s2 = open_csv(f1), open_csv(f2)
		s1len = s1.shape[0]
		s2len = s2.shape[0]
		for column in s1: 
			ordered_s1 = s1[column].sort_values().index
			ordered_s2 = s2[column].sort_values().index
			these_coords = [(ordered_s1.get_loc(x) / s1len, ordered_s2.get_loc(x) / s2len) for 
				x in ordered_s1 if clean(x) == clean(column)]
			coords_collection += these_coords
		x,y = zip(*coords_collection)
		return x,y

	f, axarr = plt.subplots(nrows=len(libs),ncols=len(libs), figsize=(9,9))
	font = {'size':11}
	plt.rc('font', **font)

	#Create the grid by iterating over libs.
	for i in range(len(libs)):
		for j in range(len(libs)):
			rlib = libs[i]
			clib = libs[j]
			subplot = axarr[i,j]
			if {'l':rlib, 'f':clib} in lib_pairs:
				prefix = 'from_' + rlib + '_to_' + clib
				fname1 = prefix + '_' + m1 + '.csv'
				fname2 = prefix + '_' + m2 + '.csv'
				if fname1 in os.listdir(os.getcwd()) and fname2 in os.listdir(os.getcwd()):
					x,y = get_coords(fname1, fname2)
					subplot.hexbin(x, y, gridsize=30, cmap='binary', bins=[pow(1.3,x) for x in range(1,22)])
					#subplot.plot(x,y, 'o') #if you want a regular dot plot
			if j != 0: subplot.yaxis.set_visible(False)
			subplot.tick_params(axis='x', which='both', bottom='off', top='off',labelbottom='off')
			subplot.axes.get_yaxis().set_ticks([])

	lib_titles = [x.replace('Single_Gene_Perturbations_from_GEO_up', 'CREEDS_up/dn_sep') for x in libs]
	for ax, col in zip(axarr[0], lib_titles): ax.set_title(col)
	for ax, row in zip(axarr[:,0], lib_titles): ax.set_ylabel(row, size='large')
	f.subplots_adjust(hspace=.03, wspace=.03)
	plt.suptitle('HexBin Plots, ' + m1 + ' (x) to ' + m2 + ' (y)', fontsize=15)
	plt.show()
	return 

if __name__ == '__main__':
	libs = ['CREEDS', 'ENCODE_TF_ChIP-seq_2015', 'ChEA_2016']
	lib_pairs = [{'l':a, 'f':b} for a in libs for b in libs if a != b]
	CREEDS_sep_pairs = [{'l':a, 'f':'Single_Gene_Perturbations_from_GEO_up'} for a in libs if a != 'CREEDS'] + [
	{'l':'Single_Gene_Perturbations_from_GEO_up', 'f':b} for b in libs if b != 'CREEDS']
	all_pairs = lib_pairs + CREEDS_sep_pairs
	all_libs = ['Single_Gene_Perturbations_from_GEO_up'] + libs
	os.chdir('results')
	#Parallel(n_jobs=1, verbose=0)(delayed(pairwise_plots)(pair) for pair in all_pairs)
	#combined_plot(all_pairs)
	#subplots(lib_pairs, libs)
	hexbin_method_comparison(libs, 'Fisher', 'RandomForest')