import csv
import os
import numpy as np
import pandas as pd
from setup import open_csv

os.chdir('libs')
for gmt in ('ChEA_2016', 'ENCODE_TF_ChIP-seq_2015', 'CREEDS'):
	print(gmt)
	f = open_csv(gmt + '_transformed.csv')
	f = f.iloc[:,:200]
	f.fillna(0, inplace=True)
	f = f[(f.T != 0).any()]
	f = f.replace(0, np.nan)
	f.index.name = gmt + '_abridged'
	f.to_csv(gmt + '_abridged_transformed.csv',sep='\t')
