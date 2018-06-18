import pandas as pd
import numpy as np
import time
import sys
import os
import csv
import timeit
from convert_files_scripts import *

ChEA = get_genesetlist('original_tf-gene_libs/ChEA_2016.txt', item_type='gmt_fname')
convert_genesetlist(ChEA, to='gvm', output_fname='gvms/ChEA.csv')

ENCODE = get_genesetlist('original_tf-gene_libs/ENCODE_TF_ChIP-seq_2015.txt', item_type='gmt_fname')
convert_genesetlist(ENCODE, to='gvm', output_fname='gvms/ENCODE.csv')

CREEDS_up = get_genesetlist('original_tf-gene_libs/Single_Gene_Perturbations_from_GEO_up.txt', item_type='gmt_fname')
CREEDS_dn = get_genesetlist('original_tf-gene_libs/Single_Gene_Perturbations_from_GEO_down.txt', item_type='gmt_fname')
CREEDS_union = combine_genesetlists(CREEDS_up, CREEDS_dn, merge_type = 'union')
convert_genesetlist(CREEDS_TFs_union, to='gvm', output_fname='gvms/CREEDS_TFs_union.csv')