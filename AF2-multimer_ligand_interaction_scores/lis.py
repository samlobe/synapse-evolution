#%%
import os
import json
import numpy as np
from MDAnalysis import Universe

# Specify the directory path
directory_path = './'

# List of updated hardcoded file names
hardcoded_files = [
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_010.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb"),
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_011.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb"),
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_012.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_012.pdb"),
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_015.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb"),
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_025.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_025.pdb"),
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_032.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb"),
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_048.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb"),
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_058.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_058.pdb"),
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_071.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb"),
    ("AQPDZ3_AQCRIPT_predicted_aligned_error_v1_seed_088.json", "AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_010.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_010.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_011.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_011.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_012.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_012.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_015.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_015.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_025.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_025.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_032.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_032.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_048.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_048.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_058.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_058.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_071.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_071.pdb"),
    ("CO1PDZ3_COCRIPT_predicted_aligned_error_v1_seed_088.json", "CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_088.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_010.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_011.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_012.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_015.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_025.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_025.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_032.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_048.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_058.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_071.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb"),
    ("CO2PDZ3_COCRIPT_predicted_aligned_error_v1_seed_088.json", "CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_010.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_011.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_011.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_012.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_015.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_015.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_025.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_025.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_032.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_048.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_048.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_058.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_071.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_071.pdb"),
    ("HSPDZ3_CRIPT_predicted_aligned_error_v1_seed_088.json", "HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb")
]

lis_scores = {}


# Extract LIS scores
for json_file, pdb_file in hardcoded_files:
    # Load the PDB file
    u = Universe(os.path.join(directory_path, pdb_file))

    # Calculate the length of Chain A
    total_residues = len(u.select_atoms('all').residues)
    n_residues_per_chain_B = 10  # Chain B has 10 amino acids
    n_residues_per_chain_A = total_residues - n_residues_per_chain_B  # Calculate length of Chain A

    # Load the JSON file
    with open(os.path.join(directory_path, json_file), 'r') as file:
        data = json.load(file)

    # Extract the pae matrix
    pae_matrix = np.array(data['predicted_aligned_error'])

    # Load interchain mask
    n_chains = 2  # Two chains: Chain A and Chain B
    interchain_mask = np.ones((n_residues_per_chain_A + n_residues_per_chain_B, n_residues_per_chain_A + n_residues_per_chain_B), dtype=bool)
    interchain_mask[:n_residues_per_chain_A, :n_residues_per_chain_A] = False
    interchain_mask[-n_residues_per_chain_B: , -n_residues_per_chain_B:] = False

    def get_LIS(pae):
        #print(f"QC: {pae}")
        #print(f"QC2: {pae[~interchain_mask]}")
        #print(f"pae shape: {pae.shape}")
        pae_low_inter = pae < 12
        #print(f"pae_low_inter before masking: {pae_low_inter}")
        pae_low_inter[~interchain_mask] = False
        #print(f"pae_low_inter after masking: {pae_low_inter}")
        pae_subset = pae[pae_low_inter]
        #print(f"pae_subset: {pae_subset}")  # Print pae_subset values
        if len(pae_subset) > 0:
            pae_subset_rescaled = pae_subset / 12 * -1 + 1
            LIS = np.mean(pae_subset_rescaled)
        else:
            LIS = 0
        return LIS


        # Calculate the LIS score
    LIS = get_LIS(pae_matrix) 

        # Add the LIS score to the dictionary
    lis_scores[pdb_file] = LIS

# Print or save the LIS scores as needed
for pdb_file, LIS in lis_scores.items():
    print(f"{pdb_file}: {LIS}")
    




# %%
    
#PLOTS HISTOGRAMS 
    
import re
import numpy as np
import matplotlib.pyplot as plt

# Input data
data = """
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.661041504539559
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6582037325038881
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_012.pdb: 0.6684457671957672
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.6613681204569055
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_025.pdb: 0.6689362398533647
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6624713765287537
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.667408631415241
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_058.pdb: 0.6610323965651834
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.669603111814346
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6593241316744427
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_010.pdb: 0.6248063623789765
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_011.pdb: 0.6248345791805093
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_012.pdb: 0.6221079246328622
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_015.pdb: 0.6226114341085272
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_025.pdb: 0.6264850869925435
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_032.pdb: 0.6230463117027176
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_048.pdb: 0.623171373200443
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_058.pdb: 0.6251224066390042
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_071.pdb: 0.6318927488282328
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_088.pdb: 0.6247199170124481
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.684631653511435
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6762743309002434
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.687462489862125
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.686888551165147
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_025.pdb: 0.6851326447954637
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6824080845013204
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.6789355231143552
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6734563929738562
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.6926524328689683
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6836160081053697
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.6841392390289449
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_011.pdb: 0.6753383668903803
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.6830945725599814
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_015.pdb: 0.6748519239083441
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_025.pdb: 0.6739291967607792
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6850753295668549
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_048.pdb: 0.6861954854262546
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6791499218575575
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_071.pdb: 0.6771690616155019
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6734230856623494

"""  # Replace ... with your entire text.

# Parse the data
categories = {'HS': [], 'AQ': [], 'CO1': [], 'CO2': []}
for line in data.strip().split('\n'):
    match = re.match(r"(HS|AQ|CO1|CO2).*?: ([\d.]+)", line)
    if match:
        category, value = match.groups()
        categories[category].append(float(value))

# Compute statistics
stats = {}
for category, values in categories.items():
    avg = np.mean(values) if values else 0
    n = len(values)
    std_dev = np.std(values) if values else 0
    stats[category] = {'avg': avg, 'n': n, 'std_dev': std_dev}

# Prepare data for plotting
x_labels = list(stats.keys())
y_values = [stats[cat]['avg'] for cat in x_labels]
std_errors = [stats[cat]['std_dev'] for cat in x_labels]

# Plot histogram with error bars
plt.figure(figsize=(8, 5))
plt.bar(x_labels, y_values, yerr=std_errors, capsize=5, color=['skyblue', 'lightgreen', 'salmon', 'orange'])
plt.xlabel('Species PDZ')
plt.ylabel('Average LIS Score For 10 Random Seeds')
plt.title('Average LIS Scores with Standard Deviation')
plt.show()

# Print results
for cat, stat in stats.items():
    print(f"{cat}: Average={stat['avg']:.4f}, n={stat['n']}, Std Dev={stat['std_dev']:.4f}")


# %%

# AESTHETIC 
    
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Input data
data = """
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.661041504539559
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6582037325038881
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_012.pdb: 0.6684457671957672
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.6613681204569055
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_025.pdb: 0.6689362398533647
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6624713765287537
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.667408631415241
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_058.pdb: 0.6610323965651834
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.669603111814346
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6593241316744427
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_010.pdb: 0.6248063623789765
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_011.pdb: 0.6248345791805093
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_012.pdb: 0.6221079246328622
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_015.pdb: 0.6226114341085272
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_025.pdb: 0.6264850869925435
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_032.pdb: 0.6230463117027176
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_048.pdb: 0.623171373200443
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_058.pdb: 0.6251224066390042
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_071.pdb: 0.6318927488282328
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_088.pdb: 0.6247199170124481
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.684631653511435
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6762743309002434
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.687462489862125
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.686888551165147
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_025.pdb: 0.6851326447954637
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6824080845013204
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.6789355231143552
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6734563929738562
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.6926524328689683
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6836160081053697
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.6841392390289449
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_011.pdb: 0.6753383668903803
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.6830945725599814
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_015.pdb: 0.6748519239083441
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_025.pdb: 0.6739291967607792
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6850753295668549
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_048.pdb: 0.6861954854262546
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6791499218575575
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_071.pdb: 0.6771690616155019
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6734230856623494

"""

# Parse the data
categories = {'CO1': [], 'CO2': [], 'AQ': [],'HS': []}
for line in data.strip().split('\n'):
    match = re.match(r"(CO1|CO2|AQ|HS).*?: ([\d.]+)", line)
    if match:
        category, value = match.groups()
        categories[category].append(float(value))

# Compute statistics
stats = {cat: {'avg': np.mean(vals), 'std_dev': np.std(vals), 'n': len(vals)} for cat, vals in categories.items() if vals}

# Prepare data for plotting
x_labels = list(stats.keys())
y_values = [stats[cat]['avg'] for cat in x_labels]
std_errors = [stats[cat]['std_dev'] for cat in x_labels]

# Aesthetic adjustments
sns.set(style="white")
plt.figure(figsize=(8, 5))
colors = sns.color_palette("Blues", len(x_labels))  # Professional colors

bars = plt.bar(x_labels, y_values, yerr=std_errors, capsize=5, color=colors, edgecolor='black', alpha=0.85)

# Add precise values on top of bars
for bar, value in zip(bars, y_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f'{value:.3f}', ha='center', fontsize=10, fontweight='bold')

plt.xlabel('Species PDZ', fontsize=12)
plt.ylabel('Average LIS Score For 10 Random Seeds', fontsize=12)
plt.title('Comparison of PDZ Binding Scores', fontsize=14, fontweight='bold')
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.show()

# Print results
for cat, stat in stats.items():
    print(f"{cat}: Average={stat['avg']:.4f}, n={stat['n']}, Std Dev={stat['std_dev']:.4f}")




# %%


import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Input data
data = """
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.661041504539559
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6582037325038881
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_012.pdb: 0.6684457671957672
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.6613681204569055
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_025.pdb: 0.6689362398533647
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6624713765287537
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.667408631415241
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_058.pdb: 0.6610323965651834
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.669603111814346
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6593241316744427
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_010.pdb: 0.6248063623789765
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_011.pdb: 0.6248345791805093
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_012.pdb: 0.6221079246328622
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_015.pdb: 0.6226114341085272
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_025.pdb: 0.6264850869925435
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_032.pdb: 0.6230463117027176
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_048.pdb: 0.623171373200443
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_058.pdb: 0.6251224066390042
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_071.pdb: 0.6318927488282328
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_088.pdb: 0.6247199170124481
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.684631653511435
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6762743309002434
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.687462489862125
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.686888551165147
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_025.pdb: 0.6851326447954637
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6824080845013204
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.6789355231143552
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6734563929738562
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.6926524328689683
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6836160081053697
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.6841392390289449
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_011.pdb: 0.6753383668903803
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.6830945725599814
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_015.pdb: 0.6748519239083441
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_025.pdb: 0.6739291967607792
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6850753295668549
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_048.pdb: 0.6861954854262546
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6791499218575575
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_071.pdb: 0.6771690616155019
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6734230856623494

"""

# Parse the data
categories = {'CO1': [], 'CO2': [], 'AQ': [], 'HS': []}
for line in data.strip().split('\n'):
    match = re.match(r"(CO1|CO2|AQ|HS).*?: ([\d.]+)", line)
    if match:
        category, value = match.groups()
        categories[category].append(float(value))

# Positive and Negative Controls
neg_controls = [0.36, 0.38, 0.36, 0.67, 0.48, 0.41]
pos_controls = [0.65, 0.66, 0.61, 0.68, 0.75]

# Compute statistics
stats = {cat: {'avg': np.mean(vals), 'std_dev': np.std(vals), 'n': len(vals)} for cat, vals in categories.items() if vals}

# Prepare data for plotting
x_labels = list(stats.keys())
y_values = [stats[cat]['avg'] for cat in x_labels]
std_errors = [stats[cat]['std_dev'] for cat in x_labels]

# Aesthetic adjustments
sns.set(style="white")
plt.figure(figsize=(8, 5))
colors = sns.color_palette("Blues", len(x_labels))  # Professional colors

# Plot bars
bars = plt.bar(x_labels, y_values, yerr=std_errors, capsize=5, color=colors, edgecolor='black', alpha=0.85)

# Add precise values on top of bars
for bar, value in zip(bars, y_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f'{value:.3f}', ha='center', fontsize=10, fontweight='bold')

# Scatter plot columns for controls
plt.scatter([-2] * len(neg_controls), neg_controls, color='red', label='Negative Controls', alpha=0.7, edgecolor='black')
plt.scatter([-1] * len(pos_controls), pos_controls, color='blue', label='Positive Controls', alpha=0.7, edgecolor='black')

# Adjust x-axis labels
plt.xticks(ticks=range(-2, len(x_labels)), labels=['Negative Control', 'Positive Control'] + x_labels, fontsize=10)
plt.xlabel('Species PDZ', fontsize=12)
plt.ylabel('Average LIS Score For 10 Random Seeds', fontsize=12)
plt.title('Comparison of PDZ Binding Scores', fontsize=14, fontweight='bold')
#plt.legend()
plt.tight_layout()
plt.show()
# %%

# updated HISTOGRAMS

import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Input data
data = """
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.661041504539559
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6582037325038881
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_012.pdb: 0.6684457671957672
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.6613681204569055
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_025.pdb: 0.6689362398533647
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6624713765287537
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.667408631415241
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_058.pdb: 0.6610323965651834
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.669603111814346
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6593241316744427
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_010.pdb: 0.6248063623789765
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_011.pdb: 0.6248345791805093
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_012.pdb: 0.6221079246328622
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_015.pdb: 0.6226114341085272
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_025.pdb: 0.6264850869925435
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_032.pdb: 0.6230463117027176
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_048.pdb: 0.623171373200443
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_058.pdb: 0.6251224066390042
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_071.pdb: 0.6318927488282328
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_088.pdb: 0.6247199170124481
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.684631653511435
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6762743309002434
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.687462489862125
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.686888551165147
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_025.pdb: 0.6851326447954637
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6824080845013204
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.6789355231143552
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6734563929738562
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.6926524328689683
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6836160081053697
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.6841392390289449
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_011.pdb: 0.6753383668903803
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.6830945725599814
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_015.pdb: 0.6748519239083441
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_025.pdb: 0.6739291967607792
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6850753295668549
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_048.pdb: 0.6861954854262546
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6791499218575575
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_071.pdb: 0.6771690616155019
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6734230856623494

"""


# Parse the data
categories = {'CO1': [], 'CO2': [], 'AQ': [], 'HS': []}
for line in data.strip().split('\n'):
    match = re.match(r"(CO1|CO2|AQ|HS).*?: ([\d.]+)", line)
    if match:
        category, value = match.groups()
        categories[category].append(float(value))

# Positive and Negative Controls
neg_controls = [0.36, 0.38, 0.36, 0.67, 0.48, 0.41]
pos_controls = [0.65, 0.66, 0.61, 0.68, 0.75]

# Compute statistics
stats = {cat: {'avg': np.mean(vals), 'std_dev': np.std(vals)} for cat, vals in categories.items() if vals}
neg_avg, neg_std = np.mean(neg_controls), np.std(neg_controls)
pos_avg, pos_std = np.mean(pos_controls), np.std(pos_controls)

# Prepare data for plotting
x_labels = ['Negative Control', 'Positive Control'] + list(stats.keys())
y_values = [neg_avg, pos_avg] + [stats[cat]['avg'] for cat in stats]
std_errors = [neg_std, pos_std] + [stats[cat]['std_dev'] for cat in stats]

# Aesthetic adjustments
sns.set(style="white")
plt.figure(figsize=(8, 5))
bar_color = 'royalblue'  # Single shade of blue for all bars

# Plot bars
bars = plt.bar(x_labels, y_values, yerr=std_errors, capsize=5, color=bar_color, edgecolor='black', alpha=0.85)

# Add precise values on top of bars
for bar, value in zip(bars, y_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f'{value:.3f}', 
             ha='center', fontsize=10, fontweight='bold')

# Adjust axes
plt.xticks(fontsize=10)
plt.xlabel('Species PDZ', fontsize=12)
plt.ylabel('Average LIS Score For 10 Random Seeds', fontsize=12)
plt.title('Comparison of PDZ Binding Scores', fontsize=14, fontweight='bold')
plt.ylim(0.3, max(y_values) + 0.05)  # Set y-axis to start at 0.3

plt.tight_layout()
plt.show()

# %%

# Used in manuscript HISTOGRAMS

import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Input data
data = """
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.661041504539559
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6582037325038881
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_012.pdb: 0.6684457671957672
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.6613681204569055
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_025.pdb: 0.6689362398533647
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6624713765287537
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.667408631415241
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_058.pdb: 0.6610323965651834
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.669603111814346
AQPDZ3_AQCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6593241316744427
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_010.pdb: 0.6248063623789765
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_011.pdb: 0.6248345791805093
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_012.pdb: 0.6221079246328622
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_015.pdb: 0.6226114341085272
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_025.pdb: 0.6264850869925435
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_032.pdb: 0.6230463117027176
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_048.pdb: 0.623171373200443
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_058.pdb: 0.6251224066390042
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_071.pdb: 0.6318927488282328
CO1PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_088.pdb: 0.6247199170124481
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.684631653511435
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_011.pdb: 0.6762743309002434
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.687462489862125
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_015.pdb: 0.686888551165147
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_025.pdb: 0.6851326447954637
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6824080845013204
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_048.pdb: 0.6789355231143552
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6734563929738562
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_071.pdb: 0.6926524328689683
CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6836160081053697
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_010.pdb: 0.6841392390289449
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_011.pdb: 0.6753383668903803
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_012.pdb: 0.6830945725599814
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_015.pdb: 0.6748519239083441
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_025.pdb: 0.6739291967607792
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_032.pdb: 0.6850753295668549
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_048.pdb: 0.6861954854262546
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_058.pdb: 0.6791499218575575
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_071.pdb: 0.6771690616155019
HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_088.pdb: 0.6734230856623494

"""

# Parse the data
categories = {'CO1': [], 'CO2': [], 'AQ': [], 'HS': []}
for line in data.strip().split('\n'):
    match = re.match(r"(CO1|CO2|AQ|HS).*?: ([\d.]+)", line)
    if match:
        category, value = match.groups()
        categories[category].append(float(value))

# Positive and Negative Controls
neg_controls = [0.36, 0.38, 0.36, 0.67, 0.48, 0.41]
pos_controls = [0.65, 0.66, 0.61, 0.68, 0.75]

# Compute statistics
stats = {cat: {'avg': np.mean(vals), 'std_dev': np.std(vals), 'n': len(vals)} for cat, vals in categories.items() if vals}

# Prepare data for plotting
x_labels = list(stats.keys())
y_values = [stats[cat]['avg'] for cat in x_labels]
std_errors = [stats[cat]['std_dev'] for cat in x_labels]

# Aesthetic adjustments
sns.set(style="white")
plt.figure(figsize=(8, 5))
color = "royalblue"

# X-axis positions
x_positions = np.arange(len(x_labels))

# Plot bars
bars = plt.bar(x_positions, y_values, yerr=std_errors, capsize=5, color=color, edgecolor='black', alpha=0.85)

# Add precise values on top of bars
for bar, value in zip(bars, y_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f'{value:.3f}', ha='center', fontsize=10, fontweight='bold')

# Compute control statistics
neg_avg, neg_std = np.mean(neg_controls), np.std(neg_controls)
pos_avg, pos_std = np.mean(pos_controls), np.std(pos_controls)

# Add bars for controls
neg_x, pos_x = -2, -1
plt.bar([neg_x, pos_x], [neg_avg, pos_avg], yerr=[neg_std, pos_std], capsize=5, color=color, edgecolor='black', alpha=0.85)

# Overlay control data points centered on bars
plt.scatter([neg_x] * len(neg_controls), neg_controls, color='red', edgecolor='black', alpha=0.7, label='Negative Controls')
plt.scatter([pos_x] * len(pos_controls), pos_controls, color='blue', edgecolor='black', alpha=0.7, label='Positive Controls')

# Adjust x-axis labels
plt.xticks(ticks=np.concatenate(([neg_x, pos_x], x_positions)), labels=['Negative Controls', 'Positive Controls'] + x_labels, fontsize=10)
plt.xlabel('Species PDZ', fontsize=12)
plt.ylabel('AF2 Multimer LIS Score', fontsize=12)
plt.title('Comparison of PDZ Binding Scores', fontsize=14, fontweight='bold')
plt.ylim(0.3, None)
plt.tight_layout()
plt.show()

# %%
