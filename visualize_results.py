import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from perform_enrichment import clean_tf, clean
from sklearn.metrics import auc
from joblib import Parallel, delayed
from collections import Counter
from setup import open_csv

#This is used to assign each method a color for the bridge plots.
color_dict = {
	'Fisher': 'black',
	'Fisher_20': 'black',
	'Fisher_50': 'black',
	'Fisher_ltf100': 'black',
	'Control': 'gray',
	'Binomial_proportion': 'C0',
	'RandomForest': 'C2',
	'BadRForest': 'C0',
	'ML_Fisher_features': 'C0',
	'Combined': 'C8',
	'CombinedFF':'C1',
	'CombinedFF2':'C0',
	'CombinedFF3':'C2',
	'CombinedFF4':'C3',
	'Z': 'C9',
	'Entropy': 'C3',
	'Gini': 'C4',
	'Forest_ne_5': '#6bd66b',
	'Forest_ne_25': '#289128',
	'Forest_ne_50': '#175317',
	'Forest_md_2': '#c8f0c8',
	'Forest_md_4': '#a9e7a9',
	'Forest_md_6': '#8adf8a',
	'Forest_md_8': '#6bd66b',
	'Forest_md_10': '#3dca3d',
	'ForestDrop': '#662ca0',
	'ForestDrop5': '#a02c66',
	'ExtraTrees': '#b8e7a1',
	'ExtraTreesDrop': '#72d043',
	'ExtraTreesDrop5': '#2f5e18',
	'ForestFisherCutoff.05': '#decaa6',
	'ForestFisherCutoff.10': '#c9a86d',
	'ForestFisherCutoff.25': '#9a7839',
	'ForestFisherCutoff.5': '#614b24',
	'ForestFisherCutoffV2.05': '#decaa6',
	'ForestFisherCutoffV2.10': '#c9a86d',
	'ForestFisherCutoffV2.25': '#9a7839',
	'ForestFisherCutoffV2.5': '#614b24',
	'FisherForestCutoff.05': '#e1c6cc',
	'FisherForestCutoff.10': '#d0a5ad',
	'FisherForestCutoff.25': '#be838f',
	'FisherForestCutoff.5': '#ad6270',
	'AdaBoost': 'C0',
	'GradientBoosting': 'C1',
	'RandomTreesEmbedding': 'C3',
	'XGBoost': 'C5',
	'LinearSVC':'C6',
}

def shorten_libnames(str_with_libnames):
	str_with_libnames = str_with_libnames.replace(
		'Single_Gene_Perturbations_from_GEO_up', 'CREEDS_sep').replace(
		'ENCODE_TF_ChIP-seq_2015', 'ENCODE').replace(
		'_2016','').replace(
		'_10-05-17', '')
	return str_with_libnames

def plot_curve(df, alg_info, prefix):
	'''
	This helper function plots a single bridge plot curve.
	df : pandas.DataFrame
		Columns looks like this: ['Fisher,x', 'Fisher,y', 'Foo3,x', 'Foo3,y']
	alg_info : tuple
		Contains the algorithm name and axis for the algorithm result to plot, e.g. ('Fisher','x')
	prefix : str
		Prefix of the file for the algorithm result being plotted, e.g. "input_ChEA_2016_into_CREEDS_tfs". 
	'''
	print('plotting')
	prefix = shorten_libnames(prefix)
	algorithm = alg_info.partition(',')[0] 
	#Scale the x_vals here.
	x_vals = [a/len(df[algorithm + ',x']) for a in df[algorithm + ',x']]
	y_vals = df[algorithm + ',y']

	#==========
	#May insert control statements here to change color, linestyle, etc., as long as plt.plot() is also modified accordingly. 
	#==========
	linewidth=2
	#==========

	if algorithm in color_dict:
		plt.plot(x_vals, y_vals, label=prefix + algorithm + '    ' + 'AUC: ' 
			+ str(np.round(auc(x_vals, y_vals), 4)), color=color_dict[algorithm], linewidth=linewidth)
	else: 
		print('plotting', algorithm)
		plt.plot(x_vals, y_vals, label=prefix + algorithm + '    ' + 'AUC: ' 
			+ str(np.round(auc(x_vals, y_vals), 4)), linewidth=linewidth)

def pairwise_plots(pair):
	'''Creates a bridge plot for enrichment results between the specified library pair.
	pair : dict
		Key 'i' contains the input library name, and key 's' contains the search library name. 
	'''
	def get_ranks(file, dn_file):
		'''
		Collects the "hit" ranks:
			the ranks where the search library annotation matches 
			the input library annotation which was used to get the enrichment score.
		Normally, `dn_file == None`. 
			However, this function can also be used to take the best rank between two CREEDS score files.
			In this case, `file` should use the up-reguated genes, and `dn_file` should use the down-reguated genes.
		'''

		i_lib = file.partition('_into_')[0].partition('input_')[2]
		s_lib = file.rpartition('_')[0].partition('_into_')[2]

		hit_ranks_collection = []
		scores = open_csv(file)
		if dn_file is not None: dn_scores = open_csv(dn_file)
		#Recall that the columns of `scores` will be input library annotations,
		#	the index of `scores` will be the search library annotations, and 
		#	the cell values will be the corresponding enrichment scores.

		for input_annot in scores:
			#Rank the search library annotations from best to worst score, as given by the 
			#	enrichment algorithm when `input_annot` is the input. 
			if dn_file is not None:
				#Get the BEST score for each feature library experiment.
				best_scores = pd.Series([min(*l) for l in zip(scores[input_annot], dn_scores[input_annot])], index=scores.index) 
				ordered_annots = best_scores.sort_values().index
			else:
				ordered_annots = scores[input_annot].sort_values().index

			#Collect the rank values for search library annotations whose corresponding tf/drug matches
			#	that of the input library annotation. (These are the "hit" ranks)
			col_clean = clean(input_annot, i_lib)
			hit_ranks = [ordered_annots.get_loc(x) for x in ordered_annots if clean(x, s_lib) == col_clean]
			hit_ranks_collection += hit_ranks

		#Return scores.shape too, which will be needed to make the graph.
		#(scores.shape should be identical between the different methods)
		return hit_ranks_collection, scores.shape

	def get_bridge_coords(hit_ranks, ranks_range, n_rankings):
		'''From the "hit" ranks, get the coordinates of the bridge plot curve.
		hit_ranks : list
			Aggregated ranks of the search lib annotation 
			whose corresponding tf/drug matches that of the input lib annotation.
		ranks_range : int
			The number of search library annotations, i.e. the range of possible ranks.
		n_rankings : int
			The number of input library annotations whose corresponding tf/drug also
			corresponds to at least one search library annotation, 
			i.e. the number of rankings that were created. 
		'''
		down_const = 1/(ranks_range - 1)
		vert_scale = len(hit_ranks)

		hits = Counter(hit_ranks)
		coords = pd.Series(0.0, index=range(ranks_range))
		coords[0] = hits[0] / vert_scale
		for x in range(1, ranks_range):
			coords[x] = coords[x - 1] + hits[x] / vert_scale - down_const
		return coords.index.values, coords.values, hits, coords

	prefix = 'input_' + pair['l'] + '_into_' + pair['f']
	print(prefix)
	rank_fname = 'rankings_' + prefix + '.csv'
	
	#all_coords will the contain all the different methods' bridge plot coordinates for this library pair.

	#Load the saved bridge plot coordinates, if any.
	if os.path.isfile(rank_fname): all_coords = pd.read_csv(rank_fname, sep='\t', index_col=0)
	else: all_coords = pd.DataFrame()
	
	#new_ranks will contain all the NEW algorithms' hit ranks for this pair.
	new_ranks = pd.DataFrame()
	#Let's begin by looking for new score files.
	for file in os.listdir(os.getcwd()):
		if file.startswith(prefix):
			#print('found', file)
			#Get the enrichment algoritm name.
			algorithm_name = str(file).partition(prefix + '_')[2].partition('.csv')[0]
			#If the algoritm is new, get and store its hit ranks.
			if str(algorithm_name + ',x') in all_coords.columns.values: continue
			if '_down' in prefix: continue
			elif '_up' in prefix: hit_ranks, r_shape = get_ranks(file, file.replace('up','down'))
			else: hit_ranks, r_shape = get_ranks(file, None)
			print(algorithm_name, '   ', len(hit_ranks), 'hits', '   ', r_shape)
			new_ranks[algorithm_name] = hit_ranks
			ranks_range, n_rankings = r_shape[0], r_shape[1]

	#If any new ranking files were found, we need to get and save their plot coordinates to all_coords for later.
	if not new_ranks.empty:
		for algorithm in new_ranks:
			#Get and store the plot coordinates.
			x, y, hits, coords = get_bridge_coords(new_ranks[algorithm].values, ranks_range, n_rankings)
			all_coords[algorithm + ',x']=x
			all_coords[algorithm + ',y']=y
		all_coords.index=range(len(x))
	#Save the plot coordinate file.
	all_coords.to_csv(rank_fname, sep='\t')

	#Plot the results for all enrichment methods, if any.
	if not all_coords.empty:
		plt.figure(1, figsize=(10,10))
		font = {'size': 12}
		plt.rc('font', **font)
		#Plot each enrichment method.
		for alg_info in all_coords:
			algorithm_name, axis = (alg_info.partition(',')[0], alg_info.partition(',')[2])
			#Filter for only certain enrichment methods here using the below if statement.
			if axis == 'x':
				plot_curve(all_coords, alg_info, '')
		plt.title(pair['l'].replace('_up', '_up/dn') + ' to ' + pair['f'].replace('_up', '_up/dn') + ' Bridge Plot')
		plt.xlabel('Rank')
		#Uncomment the line below to view only the first few ranks.
		#plt.gca().set_xlim([0,.10])
		plt.legend(prop={'size':12}, frameon=False)
		plt.show()
	return

def combined_plot(lib_df_pairs):
	'''
	Plots a single graph for all results, across all libraries and methods.
	lib_df_pairs : dict
		keys are names of the gmt libraries; values are their dfs. 
	'''
	plt.figure(2, figsize=(10,10))
	font = {'size': 12}
	plt.rc('font', **font)
	for pair in lib_df_pairs:
		prefix = 'input_' + pair['l'] + '_into_' + pair['f']
		rank_fname = 'rankings_' + prefix + '.csv'
		if os.path.isfile(rank_fname): 
			all_coords = open_csv(rank_fname)
			for alg_info in all_coords:
				algorithm_name, axis = alg_info.partition(',')[0], alg_info.partition(',')[2]
				#===========================================================================================
				#Filter for only certain enrichment methods here using an if statement.
				#===========================================================================================
				if axis == 'x' and (algorithm_name in ['RandomForest']): plot_curve(all_coords, alg_info, prefix + ' ')
				#===========================================================================================

	plt.title('Combined Bridge Plot')
	plt.xlabel('Rank')
	#Uncomment the line below to view only the first few ranks.
	#plt.gca().set_xlim([0,.10])
	plt.legend(prop={'size':10}, frameon=False)#, bbox_to_anchor=(1.05, 1), loc=2)
	plt.show()
	return

def subplots(lib_pairs, all_libs, top_10):
	'''
	Like pairwise_plots, but shows all plots as different subplots in the same figure, organized in a grid.
	lib_pairs : dict
		keys are names of the gmt libraries; values are their dfs. 
	all_libs : list-like
		the gmt libraries in lib_pairs
	'''
	f, axarr = plt.subplots(nrows=len(all_libs),ncols=len(all_libs), figsize=(35,35))
	font = {'size':20}
	plt.rc('font', **font)

	#Collect all the algorithms found so that we can create a legend at the end.
	#IMPORTANT: the legend only works if each lib_pair has the EXACT same algorithms.
	algorithms = pd.Series()

	#Create the grid by iterating over all_libs.
	for i in range(len(all_libs)):
		for j in range(len(all_libs)):
			rlib = all_libs[i]
			clib = all_libs[j]
			subplot = axarr[i,j]
			#Check if you want to plot this pair (e.g. you dont want to if rlib == clib).
			if {'l':rlib, 'f':clib} in lib_pairs:
				prefix = 'input_' + rlib + '_into_' + clib
				rank_fname = 'rankings_' + prefix + '.csv'
				if os.path.isfile(rank_fname): 
					all_coords = open_csv(rank_fname)
					for alg_info in all_coords:
						algorithm_name, axis = alg_info.partition(',')[0], alg_info.partition(',')[2]
						#===========================================================================================
						#Filter for only certain enrichment algorithms here using an if statement.
						#===========================================================================================
						if axis == 'x': 
						#===========================================================================================
							x_vals = [a/len(all_coords[algorithm_name + ',x']) for a in all_coords[algorithm_name + ',x']]
							y_vals = all_coords[algorithm_name + ',y']

							linewidth = 2.5
							if algorithm_name in color_dict: 
								algorithms[algorithm_name] = subplot.plot(x_vals, y_vals, label=algorithm_name + ' ' + str(np.round(auc(x_vals, y_vals), 4)), color=color_dict[algorithm_name], linewidth=linewidth)
							else:
								algorithms[algorithm_name] = subplot.plot(x_vals, y_vals, label=algorithm_name + ' ' + str(np.round(auc(x_vals, y_vals), 4)), linewidth=linewidth)
							#If you want to view legends for each subplot (e.g. to see the AUC), you will need to un-comment this line.
							subplot.legend(fontsize=12, loc='upper right')
			#Uncomment below to scale all subplots equally (to compare relative sizes between subplots).
			subplot.set_ylim([-.2,1])
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
	lib_titles = [shorten_libnames(x) for x in all_libs]
	for ax, col in zip(axarr[0], lib_titles): ax.set_title(col)
	for ax, row in zip(axarr[:,0], lib_titles): ax.set_ylabel(row, size='large')
	#Leave some space between the subplots.
	f.subplots_adjust(hspace=.15, wspace=.1)
	#Create a legend in the last cell (should be empty, because it is a diagonal).
	plt.legend([x for sublist in algorithms.values for x in sublist], algorithms.index, loc='upper left', ncol=1)
	#Title the plot.
	if top_10: plt.suptitle('Bridge plots inputting [row] into [column], top 10 percentile of ranks', fontsize=28)
	else: plt.suptitle('Bridge plots inputting [row] into [column]', fontsize=32)
	plt.show()
	return 

def hexbin_method_comparison(libs, m1, m2):
	'''
	Creates a grid of heatmaps comparing how two different methods rank the "hits", or "matches".
	libs : list-like
		the gmt libraries to create heatmaps for
	m1 : str
		the first method, to correspond with the x-axis
	m2 : str
		the second method, to correspond with the y-axis
	'''
	def get_coords(f1, f2):
		'''Get the coordinates of the hits.'''
		#Store coordinates here.
		coords_collection = []
		#Open the two score files.
		s1, s2 = open_csv(f1), open_csv(f2)
		#Get the maximum rank possible (to scale our plot).
		#This should be the same for s1 and s2, but just get s2len for completeness.
		s1len = s1.shape[0]
		s2len = s2.shape[0]
		#Iterate over each label library tf. 
		for column in s1: 
			ordered_s1 = s1[column].sort_values().index
			ordered_s2 = s2[column].sort_values().index
			#Get the coordinates for each hit. 
			###You can change == to != in order to view the hexbin for all the non-hits i.e. misses.###
			these_coords = [(ordered_s1.get_loc(x) / s1len, ordered_s2.get_loc(x) / s2len) for 
				x in ordered_s1 if clean(x) != clean(column)]
			#Store them. 
			coords_collection += these_coords
		x,y = zip(*coords_collection)
		#Return the list of x values, and the list of y values. 
		return x,y

	f, axarr = plt.subplots(nrows=len(libs),ncols=len(libs), figsize=(10,10))
	font = {'size':8}
	plt.rc('font', **font)

	#Create the grid by iterating over libs.
	for i in range(len(libs)):
		for j in range(len(libs)):
			rlib = libs[i]
			clib = libs[j]
			subplot = axarr[i,j]
			if {'l':rlib, 'f':clib} in lib_pairs:
				prefix = 'input_' + rlib + '_into_' + clib
				fname1 = prefix + '_' + m1 + '.csv'
				fname2 = prefix + '_' + m2 + '.csv'
				if fname1 in os.listdir(os.getcwd()) and fname2 in os.listdir(os.getcwd()):
					x,y = get_coords(fname1, fname2)
					subplot.hexbin(x, y, gridsize=30, cmap='binary', bins=[pow(1.3,x) for x in range(1,22)]) #Hexbin plot.
					#subplot.plot(x,y, 'o') #If you want a regular dot plot.
			#Only show y-axis on left-most subplots.
			if j != 0: subplot.yaxis.set_visible(False)
			#Remove box around the bottom-right subplot (legend will go here.)
			subplot.tick_params(axis='x', which='both', bottom='off', top='off',labelbottom='off')
			#Hide ticks on the y axis.
			subplot.axes.get_yaxis().set_ticks([])

	lib_titles = [shorten_libnames(x) for x in all_libs]
	for ax, col in zip(axarr[0], lib_titles): ax.set_title(col)
	for ax, row in zip(axarr[:,0], lib_titles): ax.set_ylabel(row, size='large')
	#Leave some space between the subplots.
	f.subplots_adjust(hspace=.03, wspace=.03)
	#Title the plot.
	plt.suptitle('HexBin Plots, ' + m1 + ' (x) to ' + m2 + ' (y)', fontsize=12)
	plt.show()
	return 

if __name__ == '__main__':
	tf_libs = ['ENCODE_TF_ChIP-seq_2015', 'ChEA_2016', 'CREEDS']
	drug_libs = ['1_DrugBank_EdgeList_10-05-17', 
		'2_TargetCentral_EdgeList_10-05-17',
		'3_EdgeLists_Union_10-05-17', 
		'4_EdgeLists_Intersection_10-05-17',
		'DrugBank',
		'CREEDS_Drugs',
		'DrugMatrix']

	#=============================================================
	#Choose which libraries with which to view enrichment results.
	#=============================================================
	libs = tf_libs
	#=============================================================

	lib_pairs = [{'l':a, 'f':b} for a in libs for b in libs if a != b]

	##Handle the separate results for CREEDS, if applicable. 
	#CREEDS_sep_pairs = [{'l':a, 'f':'Single_Gene_Perturbations_from_GEO_up'} for a in libs if a != 'CREEDS'] + [
	#{'l':'Single_Gene_Perturbations_from_GEO_up', 'f':b} for b in libs if b != 'CREEDS']
	#all_pairs = lib_pairs + CREEDS_sep_pairs
	#all_libs = ['Single_Gene_Perturbations_from_GEO_up'] + libs

	os.chdir('results')

	#=============================================================
	#Choose how to visualize the results.
	#=============================================================
	Parallel(n_jobs=1, verbose=0)(delayed(pairwise_plots)(pair) for pair in lib_pairs)
	#combined_plot(all_pairs)
	subplots(lib_pairs, libs, top_10=False)
	subplots(lib_pairs, libs, top_10=True)
	#hexbin_method_comparison(libs, 'Pair_Gini_ltf100_25', 'Pair_Gini_ltf100_w_25')
	#=============================================================

