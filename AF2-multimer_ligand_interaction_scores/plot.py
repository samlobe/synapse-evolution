#%%

import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

lis_data = pd.read_csv("lis_scores.csv")
controls_data = pd.read_csv("control_lis.csv")

pos_controls = controls_data[controls_data['PDB File'].str.startswith((
    "PDZ3_CITRON", "PDZ3_KALIRIN7", "PDZ3_SYNGAP", "PDZ3_VANGL2", "PDZ3_HEXAPEPTIDE"
))]['LIS Score'].values

neg_controls = controls_data[controls_data['PDB File'].str.startswith((
    "HSPDZ_TNF", "HSPDZ_ALANINE", "HSPDZ_ENDOZEPINE", "HSPDZ_GLYCINE", "HSPDZ_MOTIF", "HSPDZ_NMDA"
))]['LIS Score'].values

categories = {'CO1': [], 'CO2': [], 'AQ': [], 'HS': []}
for _, row in lis_data.iterrows():
    pdb_file = row['PDB File']  
    score = row['LIS Score']  
    
    if pdb_file.startswith('CO1'):
        categories['CO1'].append(score)
    elif pdb_file.startswith('CO2'):
        categories['CO2'].append(score)
    elif pdb_file.startswith('AQ'):
        categories['AQ'].append(score)
    elif pdb_file.startswith('HS'):
        categories['HS'].append(score)

stats = {cat: {'avg': np.mean(vals), 'std_dev': np.std(vals), 'n': len(vals)} for cat, vals in categories.items() if vals}

x_labels = list(stats.keys())
y_values = [stats[cat]['avg'] for cat in x_labels]
std_errors = [stats[cat]['std_dev'] for cat in x_labels]

sns.set(style="white")
plt.figure(figsize=(8, 5))
color = "royalblue"

x_positions = np.arange(len(x_labels))

bars = plt.bar(x_positions, y_values, yerr=std_errors, capsize=5, color=color, edgecolor='black', alpha=0.85)

for bar, value in zip(bars, y_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f'{value:.3f}', ha='center', fontsize=10, fontweight='bold')

neg_avg, neg_std = np.mean(neg_controls), np.std(neg_controls)
pos_avg, pos_std = np.mean(pos_controls), np.std(pos_controls)

neg_x, pos_x = -2, -1
plt.bar([neg_x, pos_x], [neg_avg, pos_avg], yerr=[neg_std, pos_std], capsize=5, color=color, edgecolor='black', alpha=0.85)

plt.scatter([neg_x] * len(neg_controls), neg_controls, color='red', edgecolor='black', alpha=0.7, label='Negative Controls')
plt.scatter([pos_x] * len(pos_controls), pos_controls, color='blue', edgecolor='black', alpha=0.7, label='Positive Controls')

plt.xticks(ticks=np.concatenate(([neg_x, pos_x], x_positions)), labels=['Negative Controls', 'Positive Controls'] + x_labels, fontsize=10)
plt.xlabel('Species PDZ', fontsize=12)
plt.ylabel('AF2 Multimer LIS Score', fontsize=12)
plt.title('Comparison of PDZ Binding Scores', fontsize=14, fontweight='bold')
plt.ylim(0.3, None)
plt.tight_layout()
plt.show()

# %%
