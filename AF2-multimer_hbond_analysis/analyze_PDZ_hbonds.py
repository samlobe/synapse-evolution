#%%
from glob import glob
from get_charges_and_hydrogens import get_charges_and_hydrogens
from get_hbonds import get_hbonds
import os
import pandas as pd

# Define the pattern
working_dir = "pdbs"
pattern = os.path.join(working_dir, '*unrelaxed_rank_001*seed_000.pdb')
# Find all files that match the pattern
files = glob(pattern)

dfs = []
for file in files:
    # get the charges and hydrogens
    get_charges_and_hydrogens(file)
    # get the hydrogen bonds
    dfs.append(get_hbonds(file))

#%%
# combine the dataframes
hbond_df = pd.concat(dfs, axis=1) 
# remove the part before the / in the column names
hbond_df.columns = hbond_df.columns.str.split('/').str[-1]
hbond_df.to_csv('all_PDZ_hbonds.csv')


#%%
# file = 'HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb'
# get_charges_and_hydrogens(file)
# # get the hydrogen bonds
# dfs.append(get_hbonds(file))