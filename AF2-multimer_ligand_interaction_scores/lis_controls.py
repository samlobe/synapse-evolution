#%%

import os
import json
import numpy as np
import pandas as pd
from MDAnalysis import Universe

def match_control_files(directory_path):
    json_files = {}
    pdb_files = {}
    
    for file in os.listdir(directory_path):
        if file.endswith(".json"):
            prefix = file.split("_predicted")[0]
            json_files[prefix] = file
        elif file.endswith(".pdb"):
            prefix = file.split("_unrelaxed")[0]
            pdb_files[prefix] = file

    matched_files = [(json_files[key], pdb_files[key]) for key in json_files if key in pdb_files]
    
    print(f"Matched Files: {matched_files}")
    
    return matched_files

def get_LIS(pae, interchain_mask):
    pae_low_inter = pae < 12
    pae_low_inter[~interchain_mask] = False
    pae_subset = pae[pae_low_inter]
    
    if len(pae_subset) > 0:
        pae_subset_rescaled = pae_subset / 12 * -1 + 1
        return np.mean(pae_subset_rescaled)
    return 0

def process_control_files():
    script_directory = os.path.dirname(os.path.abspath(__file__))  
    directory_path = os.path.join(script_directory, "control_files")  
    
    matched_files = match_control_files(directory_path)
    lis_scores = {}
    
    if not matched_files:
        print("No matched files found. Check your filenames and directory.")
        return
    
    for json_file, pdb_file in matched_files:
        try:
            # Load the PDB file
            u = Universe(os.path.join(directory_path, pdb_file))
            total_residues = len(u.select_atoms('all').residues)
            n_residues_per_chain_B = 10  # Chain B has 10 amino acids
            n_residues_per_chain_A = total_residues - n_residues_per_chain_B  # Chain A length
            
            # Load the JSON file
            with open(os.path.join(directory_path, json_file), 'r') as file:
                data = json.load(file)
            pae_matrix = np.array(data['predicted_aligned_error'])
            
            # Create interchain mask
            interchain_mask = np.ones((total_residues, total_residues), dtype=bool)
            interchain_mask[:n_residues_per_chain_A, :n_residues_per_chain_A] = False
            interchain_mask[-n_residues_per_chain_B:, -n_residues_per_chain_B:] = False
            
            # Calculate LIS score
            LIS = get_LIS(pae_matrix, interchain_mask)
            lis_scores[pdb_file] = LIS
        except Exception as e:
            print(f"Error processing {json_file} and {pdb_file}: {e}")
    
    if not lis_scores:
        print("No LIS scores calculated. Check your file formats and data.")
        return
    
    # Save results to CSV in the script directory
    output_path = os.path.join(script_directory, 'control_lis.csv') 
    df = pd.DataFrame(list(lis_scores.items()), columns=['PDB File', 'LIS Score'])
    df.to_csv(output_path, index=False)
    
    print(f"LIS scores saved to {output_path}")
    return lis_scores

process_control_files()

# %%
