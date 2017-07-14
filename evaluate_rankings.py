import csv
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import mlab
from get_rankings import gmt_to_df, clean, get_overlaps
from math import sqrt
from sklearn.metrics import auc

def get_ranks(file):
	r = []
	ranks = pd.DataFrame.from_csv(file, sep='\t')
	print(file, ranks.shape)
	for column in ranks:
		for x in ranks.index:
			if clean(column) == clean(ranks.at[x, column]):
				#print(column, ranks.at[x,column])
				r.append(x)
	return r, ranks.shape[0], ranks.shape[1]

def get_bridge_coords(ranks, n_f_tfs, n_overlaps):
	N = n_f_tfs * n_overlaps
	G = len(ranks)

	hits = dict.fromkeys(range(0, n_f_tfs), 0.000)

	for x in ranks:
		hits[x] += 1.000
	coords = pd.Series(0.0, index=np.arange(n_f_tfs))
	coords[0] = - (sqrt(G/N) * n_overlaps) + (sqrt(N/G) * hits[0])
	for x in range(1, n_f_tfs):
		coords[x] = coords[x - 1] - (sqrt(G/N) * n_overlaps) + (sqrt(N/G) * hits[x])
	return coords.index.values, coords.values, hits, coords

	# #histogram
	# plt.figure(1)
	# for column in dfrs:
	# 	x = dfr[column]
	# 	plt.hist(x, alpha=.5, bins=20, label=column)
	# plt.title('histogram')
	# plt.ylabel('rank')
	# plt.legend()
	# plt.show()

	# #ecdf
	# plt.figure(2)
	# for column in dfr:
	# 	x = dfr[column]
	# 	plt.hist(x, bins=20, normed=1, histtype='step', cumulative=True, label=column)
	# plt.title('ecdf')
	# plt.ylabel('rank')
	# plt.legend(loc=4)
	# plt.show()

	return coords

plt.figure(1, figsize=(10,10))
for file in os.listdir(os.getcwd()):
	if file.endswith('rankings.csv'):
		ranks, n_f_tfs, n_overlaps = get_ranks(file)
		x, y, hits, coords = get_bridge_coords(ranks, n_f_tfs, n_overlaps)
		name = file.partition('_rankings.csv')[0] + '   AUC: ' + str(auc(x,y)).expandtabs()
		plt.plot(x,y, label=name)
plt.title(mydir + ' bridge plot')
plt.xlabel('rank')
#plt.gca().set_xlim([0,50])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()