This is still very much a work in progress! 

Briefly, to use:
*run transform_CREEDS in /libs to create CREEDS_transformed.csv
*download mouse_matrix.h5 and human_matrix.h5 from the ARCHS4 website and put in /libs. 
*run ARCHS4_reader in /libs to create .h5 files for the gene correlation matrices
*now in the main folder, run get_rankings.py
*[does not work, but almost there] run evaluate_rankings.py to get the results