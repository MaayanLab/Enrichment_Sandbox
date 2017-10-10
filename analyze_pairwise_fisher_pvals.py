import os
import csv
import pandas as pd
import numpy as np
import pickle
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

def minus_results(prefix, minus_type):
	hits, misses = [], []
	os.chdir('fisher_pairwise_results')
	for file in os.listdir(os.getcwd()):
		if file.startswith(prefix):
			df = pd.read_csv(file, sep='\t', index_col=False)
			df['tfi'] = df['tfi'].apply(clean)
			df['tfj'] = df['tfj'].apply(clean)
			tf = clean(file.partition(prefix + '_')[2].partition('_pair_fisher_pvals.csv')[0])
			if minus_type == 'minuend':
				hits += to_list(df.loc[df['tfi'] == tf]['i minus j']) + to_list(df.loc[df['tfj'] == tf]['j minus i']) 
				misses += list(df.loc[df['tfi'] != tf]['i minus j']) + list(df.loc[df['tfj'] != tf]['j minus i']) 

				###temporary
				#hit_dfsubset1 = df.loc[(df['tfi'] == tf) & (df['i minus j'] < 1e-100)][['tfi', 'tfj', 'i minus j']]
				#if not hit_dfsubset1.empty: 
					#print(file)
					#print(hit_dfsubset1)
				#hit_dfsubset2 = df.loc[(df['tfj'] == tf) & (df['j minus i'] < 1e-100)][['tfi', 'tfj', 'j minus i']]
				#if not hit_dfsubset2.empty:
					#print(file)
					#print(hit_dfsubset2)
				###temporary

			elif minus_type == 'subtrahend':
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

def subplots(lib_pairs, all_libs, result_type, plot_type='none'):
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

				if plot_type=='bridge':
					fname = 'coords_' + prefix + '_' + result_type + '.csv'
					if fname in os.listdir('fisher_pairwise_results'):
						os.chdir('fisher_pairwise_results')
						final_coords = pd.read_csv(fname, sep='\t', index_col=0, header=None)
						os.chdir('..')
					else:
						print('ERROR')
						#continue #temp: make sure below code doesen't run, because I know files are already created.
						series_hits = pd.Series(hits, index=['hit'] * len(hits))
						series_misses = pd.Series(misses, index=['miss'] * len(misses))
						hits_and_misses = pd.DataFrame(pd.concat([series_hits, series_misses], axis=0), columns=['x coords'])
						hits_and_misses.sort_values(by = 'x coords', inplace=True)

						final_coords = pd.Series(index=sorted(set(hits_and_misses['x coords'])))
						down = 1/(len(misses))
						up = 1/(len(hits))
						hits_and_misses['y coords'] = 0.0
						if hits_and_misses.index[0] == 'hit': hits_and_misses['y coords'][0] = up
						else: hits_and_misses['y coords'][0] = -down
						final_coords[0] = hits_and_misses['y coords'][0]
						for k in range(1, hits_and_misses.shape[0]):
							if hits_and_misses.index[k] == 'hit': 
								hits_and_misses.iat[k, 1] = hits_and_misses.iat[k-1,1] + up
							else: hits_and_misses.iat[k, 1] = hits_and_misses.iat[k-1,1] - down
							final_coords[hits_and_misses.iat[k,0]] = hits_and_misses.iat[k,1]

						os.chdir('fisher_pairwise_results')
						print(fname)
						final_coords.to_csv(fname,sep='\t')
						os.chdir('..')

				if plot_type=='bridge':
					subplot.plot(final_coords.index, final_coords.values)
				else:
					if plot_type =='log':
						hits = [max(x, 1e-200) for x in hits]
						misses = [max(x, 1e-200) for x in misses]
						xrange_value=(1e-200,1)
						bins = np.logspace(np.log10(1e-200),np.log10(1),200, base=10)
						hit_weights = None
						miss_weights = None
					elif plot_type =='lessthan.01':
						hits = [x for x in hits if x < .1]
						misses = [x for x in misses if x < .1]
						xrange_value = (0,.01)
						bins=50
						hit_weights = np.ones_like(hits)/float(len(hits))
						miss_weights = np.ones_like(misses)/float(len(misses))
					else: 
						xrange_value = (0,1)
						bins=50
						hit_weights = np.ones_like(hits)/float(len(hits))
						miss_weights = np.ones_like(misses)/float(len(misses))

					subplot.hist(hits, bins=bins, range=xrange_value, alpha=0.5, color='C0', label=result_type,weights=hit_weights)
					subplot.hist(misses, bins=bins, range=xrange_value, alpha=0.5, color='C1', label=result_type,weights=miss_weights)

			#Only show y-axis on left-most subplots.
			if plot_type=='log':
				subplot.axes.set_ylim(1e0,1e6)
				subplot.set_yscale('log')
				subplot.axes.get_yaxis().set_ticks([1e0,1e3,1e6])
				subplot.set_xscale('log')
				subplot.get_xaxis().set_ticks([1e-200, 1e-100, 1])
			elif plot_type=='lessthan.01':
				subplot.axes.set_ylim(0,.8)
				subplot.axes.get_yaxis().set_ticks([.4,.8])
				subplot.get_xaxis().set_ticks([.005, .01])
			else:
				subplot.axes.set_ylim(-.2,.5)
				subplot.axes.get_yaxis().set_ticks([0,.2,.4])
				subplot.get_xaxis().set_ticks([.05, 1])	

			if j != 0: subplot.yaxis.set_visible(False)
			#Remove box around the bottom-right subplot (legend will go here.)
			if j == i: 
				if j == 0: 
					for side in ['top', 'right', 'bottom', 'left']: subplot.spines[side].set_visible(False)
					subplot.get_xaxis().set_ticks([])
					subplot.xaxis.set_visible(False)
				else: subplot.axis('off')

			# if j == 2 & i == 2: plt.axis('off')
			# if j == 1 & i == 1: plt.axis('off')
			# subplot.axes.set_ylim(0,.5)
			# Hide ticks -- although axis='both', this only seems to affect the x-axis.
			# subplot.tick_params(axis='x', which='both', bottom='off', top='off',labelbottom='off')
			# Hide ticks on the y axis.
			# subplot.axes.get_yaxis().set_ticks([0, .25, .5])

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
	#os.chdir('..')
	return 

def add_single_pvals(lib_pairs, all_libs):
	for pair in lib_pairs:
		f_lib = pair['f']
		l_lib = pair['l']
		prefix = 'from_' + f_lib + '_to_' + l_lib
		os.chdir('results')
		single = pd.read_csv(prefix + '_Fisher.csv', sep='\t', index_col=0)
		#print(single.index)
		#print(single.columns)
		os.chdir('..')
		os.chdir('fisher_pairwise_results')
		for file in os.listdir(os.getcwd()):
			if file.startswith(prefix):
				tf = file.partition(prefix + '_')[2].partition('_pair_fisher_pvals.csv')[0]
				tf = tf.replace('fslash', '/')
				df = pd.read_csv(file, sep='\t', index_col=False)
				if 'i' in df.columns: 
					print(file, 'already done')
					continue
				print(file)
				for row in df.index:
					df.at[row,'i'] = single.at[df.at[row,'tfi'], tf]
					df.at[row,'j'] = single.at[df.at[row,'tfj'], tf]
				df.to_csv(file, sep='\t', index=False)
		os.chdir('..')

def combined_bridge_subplot(lib_pairs, all_libs):
	f, axarr = plt.subplots(nrows=len(all_libs),ncols=len(all_libs), figsize=(15,15))
	font = {'size':20}
	plt.rc('font', **font)
	methods=pd.Series()
	#Create the grid by iterating over all_libs.
	for i in range(len(all_libs)):
		for j in range(len(all_libs)):
			f_lib = all_libs[i]
			l_lib = all_libs[j]
			subplot = axarr[i,j]
			#Check if you want to plot this pair (e.g. you dont want to if f_lib == l_lib).
			if {'f':f_lib, 'l':l_lib} in lib_pairs:
				prefix = 'from_' + f_lib + '_to_' + l_lib
				print(prefix)
				for result_type in ('single', 'union', 'intersection', 'subtrahend', 'minuend'):
					fname = 'coords_' + prefix + '_' + result_type + '.csv'
					if fname in os.listdir('fisher_pairwise_results'):
						os.chdir('fisher_pairwise_results')
						final_coords = pd.read_csv(fname, sep='\t', index_col=0, header=None)
						methods[result_type] = subplot.plot(final_coords.index, final_coords.values, label=result_type)
						os.chdir('..')
					else: print('ERROR', result_type)

			subplot.axes.set_ylim(-.2,.4)
			subplot.axes.get_yaxis().set_ticks([-.2,0,.2,.4])
			subplot.axes.set_xlim(0,.1)
			subplot.get_xaxis().set_ticks([0, .05, .1])

			if j != 0: subplot.yaxis.set_visible(False)
			#Remove box around the bottom-right subplot (legend will go here.)
			if j == i: 
				if j == 0: 
					for side in ['top', 'right', 'bottom', 'left']: subplot.spines[side].set_visible(False)
					subplot.get_xaxis().set_ticks([])
					subplot.xaxis.set_visible(False)
				else: subplot.axis('off')

	#Label the rows and columns of the figure.
	lib_titles = [x.replace('Single_Gene_Perturbations_from_GEO_up', 'CREEDS_sep').replace('ENCODE_TF_ChIP-seq_2015', 'ENCODE').replace('_2016','') for x in all_libs]
	for ax, col in zip(axarr[0], lib_titles): ax.set_title(col)
	for ax, row in zip(axarr[:,0], lib_titles): ax.set_ylabel(row, size='large')
	#Leave some space between the subplots.
	f.subplots_adjust(hspace=.3, wspace=.2)
	#Create a legend in the last cell (should be empty, because it is a diagonal).
	plt.legend([x for sublist in methods.values for x in sublist], methods.index, loc='upper left', ncol=1)
	#Title the plot.
	plt.suptitle('p-value <.10 bridge plots, from [row] to [column]', fontsize=32)
	plt.show()
	return methods
	#os.chdir('..')


all_libs = ['ENCODE_TF_ChIP-seq_2015_abridged', 'ChEA_2016_abridged', 'CREEDS_abridged']
lib_pairs = [{'f':i,'l':j} for i in all_libs for j in all_libs if i != j]
#for result_type in ('single','intersection','union','minuend', 'subtrahend'):
    #subplots(lib_pairs, all_libs, result_type, plot_type='bridge')
m = combined_bridge_subplot(lib_pairs, all_libs)
#add_single_pvals(lib_pairs, all_libs)