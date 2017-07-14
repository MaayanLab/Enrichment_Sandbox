import csv
import os
import pickle
import numpy as np
import pandas as pd

#Output: A single composite csv file

def gmt_to_dict(lib_name):
	'''Input: gmt file from Enrichr. Output: dict with tfs as keys and genes as values'''

	print('converting gmt file to dict')
	with open(lib_name) as f:
		reader = csv.reader(f)
		row_count = sum(1 for row in reader)
	for row in np.arange(row_count):
		#read one row at a time, since rows are of uneven length
		raw_data = (pd.read_csv(lib_name, skiprows=row, nrows=1, sep='\t', header=None, index_col=None).fillna('').values).tolist()[0]
		#remove null values from the row data
		raw_data = [str(w).replace(',1.0', '') for w in raw_data]
		data = [x for x in raw_data if x != '']
		print(str(row) + ' out of ' + str(row_count) + ' in ' + lib_name)
		#use the tf name as the key, and all its correspondnig genes as the respective value
		result[data[0]] = set(data[1:])
		#remove null values again, just in case
		result[data[0]].discard('')
	return result

#create a list of names for each library comprising the CREEDS dataset
libs = ['Single_Gene_Perturbations_from_GEO_down.txt','Single_Gene_Perturbations_from_GEO_up.txt']

#convert each gmt library to a dict
dicts = [gmt_to_dict(x) for x in libs]

#combine the dicts into a single dict
#(this works because both libraries have the same transcription factors)
combined = {}
for k in dicts[0]:
	combined[k] = list(set(dicts[0][k]) | set(dicts[1][k]))

#convert the dict to a dataframe, one key at a time 
#this appears to be the fastest method I have found
df = pd.DataFrame(False, index = [''], columns = [''], dtype=bool)
for k in combined:
	print(k)
	s = pd.DataFrame(True, index = combined[k], columns = [k], dtype=bool)
	df = pd.concat([df,s], axis=1)

#once again, ensure no null values in data
#then, save as csv
df = df[pd.notnull(df.index)].fillna(False)
df = df.loc[pd.notnull(df.index)]
df.drop('', inplace=True)
df.drop('', axis=1, inplace=True)
df.to_csv('CREEDS_transformed.csv', sep='\t')