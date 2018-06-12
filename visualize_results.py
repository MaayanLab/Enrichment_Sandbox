import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from perform_enrichment import clean
from sklearn.metrics import auc
from joblib import Parallel, delayed
from collections import Counter

GET_DTARGET_LIBNAME = {
	'DrugBank_5OrMoreTargets':'DrugBank',
	'TargetCentral_5OrMoreTargets':'TargetCentral', 
	'RepurposeHub_5OrMoreTargets':'RepurposeHub', 
	'DGIdb_5OrMoreTargets':'DGIdb', 
	'DrugCentral_5OrMoreTargets':'DrugCentral',
	'STITCH_500Cutoff_5OrMoreTargets_noENSP':'STITCH_500Cutoff',
	'STITCH_600Cutoff_5OrMoreTargets_noENSP':'STITCH_600Cutoff',
	'STITCH_700Cutoff_5OrMoreTargets_noENSP':'STITCH_700Cutoff',
	'STITCH_800Cutoff_5OrMoreTargets_noENSP':'STITCH_800Cutoff',
	'DTCommons_5OrMoreTargets':'DTCommons',
	}

def open_csv(fname):
	return pd.read_csv(fname, index_col=0, sep='\t', low_memory=False, encoding='Latin-1')

def plot_curve(df, alg_info):
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
	algorithm = alg_info.partition(',')[0] 
	#Scale the x_vals here.
	x_vals = [a/len(df[algorithm + ',x']) for a in df[algorithm + ',x']]
	y_vals = df[algorithm + ',y']

	#==========
	#May insert control statements here to change color, linestyle, etc., as long as plt.plot() is also modified accordingly. 
	#==========
	linewidth=2
	#==========

	print('plotting', algorithm)
	plt.plot(x_vals, y_vals, label=str(np.round(auc(x_vals, y_vals), 4)), linewidth=linewidth)

def pairwise_plots(pair):
	'''Creates a bridge plot for enrichment results between the specified library pair.
	pair : dict
		Key 'i' contains the input library name, and key 's' contains the search library name. 
	'''
	def get_ranks(file, dn_file):
		'''
		Collects the "match" ranks:
			the ranks where the search library annotation matches 
			the input library annotation which was used to get the enrichment score.
		Normally, `dn_file == None`. 
			However, this function can also be used to take the best rank between two CREEDS score files.
			In this case, `file` should use the up-reguated genes, and `dn_file` should use the down-reguated genes.
		'''

		i_lib = file.partition('_into_')[0].partition('input_')[2]
		s_lib = file.rpartition('_')[0].partition('_into_')[2]

		match_ranks_collection = []
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
			#	that of the input library annotation. (These are the "match" ranks)
			col_clean = clean(input_annot)
			match_ranks = [ordered_annots.get_loc(x) for x in ordered_annots if clean(x) == col_clean]
			match_ranks_collection += match_ranks

		#Return scores.shape too, which will be needed to make the graph.
		#(scores.shape should be identical between the different methods)
		return match_ranks_collection, scores.shape

	def get_bridge_coords(match_ranks, ranks_range, n_rankings):
		'''From the "match" ranks, get the coordinates of the bridge plot curve.
		match_ranks : list
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
		vert_scale = len(match_ranks)

		matches = Counter(match_ranks)
		coords = pd.Series(0.0, index=range(ranks_range))
		coords[0] = matches[0] / vert_scale
		for x in range(1, ranks_range):
			coords[x] = coords[x - 1] + matches[x] / vert_scale - down_const
		return coords.index.values, coords.values, matches, coords

	ilib_fname, slib_fname = pair
	prefix = 'input_' + ilib_fname + '_into_' + slib_fname
	rank_fname = 'rankings_' + prefix + '.csv'
	print(rank_fname)
	
	#all_coords will the contain all the different methods' bridge plot coordinates for this library pair.

	#Load the saved bridge plot coordinates, if any.
	if os.path.isfile(rank_fname): all_coords = pd.read_csv(rank_fname, sep='\t', index_col=0)
	else: all_coords = pd.DataFrame()
	
	#new_ranks will contain all the NEW algorithms' match ranks for this pair.
	new_ranks = pd.DataFrame()
	#Let's begin by looking for new score files.
	for file in os.listdir(os.getcwd()):
		if file.startswith(prefix):
			#print('found', file)
			#Get the enrichment algoritm name.
			algorithm_name = str(file).partition(prefix + '_')[2].partition('.csv')[0]
			#If the algoritm is new, get and store its match ranks.
			if str(algorithm_name + ',x') in all_coords.columns.values: continue
			if '_down' in prefix: continue
			elif '_up' in prefix: match_ranks, r_shape = get_ranks(file, file.replace('up','down'))
			else: match_ranks, r_shape = get_ranks(file, None)
			print(algorithm_name, '   ', len(match_ranks), 'matches', '   ', r_shape)
			new_ranks[algorithm_name] = match_ranks
			ranks_range, n_rankings = r_shape[0], r_shape[1]

	#If any new ranking files were found, we need to get and save their plot coordinates to all_coords for later.
	if not new_ranks.empty:
		for algorithm in new_ranks:
			#Get and store the plot coordinates.
			x, y, matches, coords = get_bridge_coords(new_ranks[algorithm].values, ranks_range, n_rankings)
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
				plot_curve(all_coords, alg_info)
		plt.title(ilib_fname.replace('_up', '_up/dn') + ' to ' + slib_fname.replace('_up', '_up/dn') + ' Bridge Plot')
		plt.xlabel('Rank')
		#Uncomment the line below to view only the first few ranks.
		#plt.gca().set_xlim([0,.10])
		plt.legend(prop={'size':12}, frameon=False)
		plt.show()
	return

#Example of a function used to create a custom plot.
def druglib_target_comparison(top_pct=None):
	dtarget_libs = (
		'DrugBank_5OrMoreTargets',
		'TargetCentral_5OrMoreTargets', 
		'RepurposeHub_5OrMoreTargets', 
		'DGIdb_5OrMoreTargets', 
		'DrugCentral_5OrMoreTargets',
		'STITCH_800Cutoff_5OrMoreTargets_noENSP',
		'STITCH_700Cutoff_5OrMoreTargets_noENSP',
		'STITCH_600Cutoff_5OrMoreTargets_noENSP',
		'STITCH_500Cutoff_5OrMoreTargets_noENSP',
		'DTCommons_5OrMoreTargets'
		)
	ppi_libs = (
		'huMAP',
		'BioGRID',
		'ARCHS4_human')
	dperturb_libs = ('CREEDS', 'CREEDS', 'LINCS', 'LINCS')
	dperturb_lib_titles = ('CREEDS as search', 'CREEDS as input', 'LINCS as search', 'LINCS as input')

	f, axarr = plt.subplots(nrows=len(ppi_libs),ncols=len(dperturb_libs), figsize=(30,20))
	font = {'size':25}
	plt.rc('font', **font)

	#Collect all the algorithms found so that we can create a legend at the end.
	#IMPORTANT: the legend only works if each lib_pair has the EXACT same algorithms.
	algorithms = pd.Series()

	#Create the grid by iterating over all_libs.
	for i in range(len(ppi_libs)):
		for j in range(len(dperturb_libs)):
			ppi_lib = ppi_libs[i]
			dperturb_lib= dperturb_libs[j]

			if 'as search' in dperturb_lib_titles[j]: enrichment_direction = 'fwd'
			else: enrichment_direction = 'rev'

			subplot = axarr[i,j]

			for dtarget_lib in dtarget_libs:
				if enrichment_direction == 'fwd':
					prefix = 'input_' + dtarget_lib + '_expanded_with_' + ppi_lib + '_into_' + dperturb_lib
				else:
					prefix = 'input_' + dperturb_lib + '_into_' + dtarget_lib + '_expanded_with_' + ppi_lib
				rank_fname = 'rankings_' + prefix + '.csv'
				print(rank_fname)
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
							algorithms[GET_DTARGET_LIBNAME[dtarget_lib]] = subplot.plot(x_vals, y_vals, label= str(np.round(auc(x_vals, y_vals), 3)), linewidth=linewidth)
							#If you want to view legends for each subplot (e.g. to see the AUC), you will need to un-comment this line.
							subplot.legend(fontsize=20, loc='upper right')
							#Uncomment below to scale all subplots equally (to compare relative sizes between subplots).
							subplot.set_ylim([-.2,1])
			if top_pct is not None: subplot.set_xlim([0,top_pct])
			#Hide x ticks.
			subplot.tick_params(axis='x', which='both', bottom='off', top='off',labelbottom='off')
			#Hide ticks on the y axis.
			if j != 0: subplot.axes.get_yaxis().set_ticks([])

	#Label the rows and columns of the figure.
	for ax, row in zip(axarr[:,0], [x for x in ppi_libs]): ax.set_ylabel(row, size='large')
	for ax, col in zip(axarr[0], dperturb_lib_titles): ax.set_title(col)
	#Leave some space between the subplots.
	f.subplots_adjust(hspace=.13, wspace=.1, left=.3)
	#Create a legend in the last cell (should be empty, because it is a diagonal).
	plt.figlegend([x for sublist in algorithms.values for x in sublist], algorithms.index, loc=10, bbox_to_anchor=(.1,.5))
		#bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	#Title the plot.
	if top_pct is not None: plt.suptitle('Drug library bridge plots, top ' + 
		str(top_pct).partition('.')[2] + ' percentile of ranks', fontsize=35)
	else: plt.suptitle('Drug library bridge plots', fontsize=35)
	#Save and display results.
	if top_pct is None: plt.savefig('druglib_bridgeplot.png',bbox_inches='tight')
	else: plt.savefig('druglib_bridgeplot_' + str(top_pct).partition('.')[2] + 'th_pctile.png', bbox_inches='tight')
	plt.show()
	return 

def combined_plot(file_patterns):
	'''
	Plots a single graph for all results, across all libraries and methods.
	file_patterns : tuple
		If all in file name, plots the results.
	'''
	algorithms = pd.Series()

	plt.figure(2, figsize=(10,10))
	font = {'size': 25}
	plt.rc('font', **font)

	for dtarget_lib in (
		'DrugBank_5OrMoreTargets',
		'TargetCentral_5OrMoreTargets', 
		'RepurposeHub_5OrMoreTargets', 
		'DGIdb_5OrMoreTargets', 
		'DrugCentral_5OrMoreTargets',
		'STITCH_500Cutoff_5OrMoreTargets_noENSP',
		'STITCH_600Cutoff_5OrMoreTargets_noENSP',
		'STITCH_700Cutoff_5OrMoreTargets_noENSP',
		'STITCH_800Cutoff_5OrMoreTargets_noENSP',
		'DTCommons_5OrMoreTargets'
	):

		fname = 'rankings_input_LINCS_L1000_Chem_Pert_into_' + dtarget_lib + '_expanded_with_BioGRID.csv'
		prefix = fname.partition('rankings')[2].partition('.csv')[0]
		all_coords = open_csv(fname)
		for alg_info in all_coords:
			algorithm_name, axis = alg_info.partition(',')[0], alg_info.partition(',')[2]
			if axis != 'x': continue
			dtarget_lib = prefix.partition('_into_')[2].partition('_expanded_with_')[0]
			x_vals = [a/len(all_coords[algorithm_name + ',x']) for a in all_coords[algorithm_name + ',x']]
			y_vals = all_coords[algorithm_name + ',y']
			#===========================================================================================
			#Filter for only certain enrichment methods here using an if statement.
			#===========================================================================================
			algorithms[GET_DTARGET_LIBNAME[dtarget_lib]] = plt.plot(
                    x_vals, y_vals, label= GET_DTARGET_LIBNAME[dtarget_lib] + ', AUC: ' + str(
                            np.round(auc(x_vals, y_vals), 3)), linewidth=3.5)
			#===========================================================================================
	plt.legend(fontsize=18, loc=2, bbox_to_anchor=(.3, .32))
	plt.title('Bridge Plots, Top 50 Percentile')
	plt.xlabel('Rank')
	#Uncomment the line below to view only the first few ranks.
	plt.gca().set_xlim([0,.50])
	#plt.legend(prop={'size':10}, frameon=False)#, bbox_to_anchor=(1.05, 1), loc=2)
	plt.savefig('bridgeplot_top50.eps', bbox_inches='tight', format='eps', dpi=1000)
	plt.show()
	return

if __name__ == '__main__':

	#========================================================================================
	#Edit this section to identify what pairs to perform enrichment analysis with.
	#========================================================================================
	gvm_fnames = ['gvms//' + fname for fname in os.listdir('gvms') if 'gvm.' in fname]
	gvm_fnames = gvm_fnames + ['gvms//expanded//' + fname for fname in os.listdir('gvms//expanded') 
		if 'gvm.' in fname]

	target_libs = [fname for fname in gvm_fnames if ('CREEDS' not in fname and 'LINCS' not in fname)]
	perturb_libs = [fname for fname in gvm_fnames if ('CREEDS' in fname or 'LINCS' in fname)]

	#We will iterate over each (target, perturb) and (perturb, target) pair.
	fwd_pairs = [(target,perturb) for target in target_libs for perturb in perturb_libs]
	bck_pairs = [(perturb,target) for target in target_libs for perturb in perturb_libs]
	lib_pairs = fwd_pairs + bck_pairs
	#========================================================================================

	os.chdir('results')

	for pair in lib_pairs: pairwise_plots(pair)
	#Parallel(n_jobs=2, verbose=0)(delayed(pairwise_plots)(pair) for pair in lib_pairs)

	#druglib_target_comparison(top_pct=None)
	#druglib_target_comparison(top_pct=.15)
	#combined_plot(['input_LINCS','_expanded_with_BioGRID'])