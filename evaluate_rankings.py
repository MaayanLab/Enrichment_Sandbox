import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from get_rankings import clean
from math import sqrt
from sklearn.metrics import auc
from joblib import Parallel, delayed

def pairwise_plots(pair):

	def get_ranks(l_name, f_name):
		ranks = pd.read_csv(file, sep='\t', index_col=0)
		return [x for x in ranks.index for l_tf in ranks.columns if clean(l_tf) == clean(ranks.at[x, l_tf])], ranks.shape

	def get_bridge_coords(hit_ranks, ranks_range, n_rankings):
		N = ranks_range * n_rankings
		G = len(hit_ranks)

		hits = dict.fromkeys(range(0, ranks_range), 0.000)

		for x in hit_ranks:
			hits[x] += 1.000
		coords = pd.Series(0.0, index=np.arange(ranks_range))
		coords[0] = - (sqrt(G/N) * n_rankings) / (sqrt(N/G) * G) + (sqrt(N/G) * hits[0]) / (sqrt(N/G) * G)
		for x in range(1, ranks_range):
			coords[x] = coords[x - 1] - (sqrt(G/N) * n_rankings) / (sqrt(N/G) * G) + (sqrt(N/G) * hits[x]) / (sqrt(N/G) * G)
		return coords.index.values, coords.values, hits, coords

	prefix = 'from_' + pair['l'] + '_to_' + pair['f']
	rank_fname = 'rankings_' + prefix + '.csv'
	
	if os.path.isfile(rank_fname): agg_c = pd.read_csv(rank_fname, sep='\t', index_col=0)
	else: agg_c = pd.DataFrame()

	agg_r = pd.DataFrame()
	for file in os.listdir(os.getcwd()):
		if file.startswith(prefix):
			m_name = str(file).partition(prefix + '_')[2].partition('.csv')[0]
			if str(m_name + ',x') in agg_c.columns.values: continue
			r, r_shape = get_ranks(pair['l'], pair['f'])
			print(m_name, r_shape)
			agg_r[m_name] = r

	if not agg_r.empty:
		ranks_range, n_rankings = r_shape[0], r_shape[1]
		for column in agg_r:
			x, y, hits, coords = get_bridge_coords(agg_r[column].values, ranks_range, n_rankings)
			agg_c[column + ',x']=x
			agg_c[column + ',y']=y
	agg_c.to_csv(rank_fname, sep='\t')

	#bridge
	plt.figure(3, figsize=(12,12))
	font = {'size': 12}
	plt.rc('font', **font)
	for column in agg_c:
		col = (column.partition(',')[0], column.partition(',')[2])
		if col[1] == 'x':
			name = col[0] #+ '    ' + 'AUC: ' + str(int(auc(x,y)))
			plt.plot(agg_c[name + ',x'], agg_c[name + ',y'], label=name)
	plt.title(pair['l'] + ' to ' + pair['f'] + ' Bridge Plot')
	plt.xlabel('Rank')
	#plt.gca().set_xlim([0,50])
	plt.legend(bbox_to_anchor=(1, 1), prop={'size':10}, frameon=False)
	plt.show()

	# #histogram DEPRECATED
	# plt.figure(1)
	# for column in agg_rs:
	# 	x = agg_r[column]
	# 	plt.hist(x, alpha=.5, bins=20, label=column)
	# plt.title('histogram')
	# plt.ylabel('rank')
	# plt.legend()
	# plt.show()
	#agg_rs.plot.hist(alpha=.5, bins=20)

	# #ecdf DEPRECATED
	# plt.figure(2)
	# for column in agg_r:
	# 	x = agg_r[column]
	# 	plt.hist(x, bins=20, normed=1, histtype='step', cumulative=True, label=column)
	# plt.title('ecdf')
	# plt.ylabel('rank')
	# plt.legend(loc=4)
	# plt.show()

	return

def combined_plot(lib_df_pairs):

	plt.figure(3, figsize=(12,12))
	font = {'size': 12}
	plt.rc('font', **font)
	for pair in lib_df_pairs:
		prefix = 'from_' + pair['l'] + '_to_' + pair['f']
		rank_fname = 'rankings_' + prefix + '.csv'
		if os.path.isfile(rank_fname): 
			agg_c = pd.read_csv(rank_fname, sep='\t', index_col=0)
			for column in agg_c:
				col = (column.partition(',')[0], column.partition(',')[2])
				if col[0] != 'Control' and col[1] == 'x':
					name = col[0] #+ '    ' + 'AUC: ' + str(int(auc(x,y)))
					x_vals = [a / len(agg_c[name + ',x']) for a in agg_c[name + ',x']]
					plt.plot(x_vals, agg_c[name + ',y'], label=name + ' ' + prefix)
	plt.title('Aggregated Bridge Plot')
	plt.xlabel('Rank')
	#plt.gca().set_xlim([0,.05])
	plt.legend(bbox_to_anchor=(1, 1), prop={'size':10}, frameon=False)
	plt.show()
	return

if __name__ == '__main__':
	all_libs = ['CREEDS', 'ENCODE_TF_ChIP-seq_2015', 'ChEA_2016']
	lib_pairs = [{'l':a, 'f':b} for a in all_libs for b in all_libs if a != b]
	os.chdir('results')
	Parallel(n_jobs=1, verbose=0)(delayed(pairwise_plots)(pair) for pair in lib_pairs)
	combined_plot(lib_pairs)