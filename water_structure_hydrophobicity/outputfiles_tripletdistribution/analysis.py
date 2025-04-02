#%%

import pandas as pd
import numpy as np
import os

def calculate_tetrahedral_fraction(file_path, angle_range=np.arange(102.5, 120, 5)):
    """
    Calculates the tetrahedral fraction using the trapezoidal rule
    for angles between 100° and 120° from a given CSV file.
    """
    # Read CSV
    df = pd.read_csv(file_path)
    
    # Extract relevant columns
    valid_columns = [str(angle) for angle in angle_range if str(angle) in df.columns]
    if not valid_columns:
        raise ValueError(f"No valid angle columns found in {file_path}")
    
    # Convert column names to numeric angles
    angles = np.array([float(col) for col in valid_columns])
    values = df[valid_columns].sum().values  # Sum over all rows for each angle bin
    
    # Compute trapezoidal integration
    tetrahedral_fraction = np.trapz(values, angles)
    return tetrahedral_fraction

if __name__ == "__main__":
    # Directory containing CSV files
    data_dir = "/Users/riyanilkant/desktop/hydrophobicity_hbond"
    
    # List of CSV files to process
    file_names = [
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed52.csv"
    ]
    
    # Compute tetrahedral fractions for each file
    results = {}
    for file_name in file_names:
        file_path = os.path.join(data_dir, file_name)
        if os.path.exists(file_path):
            results[file_name] = calculate_tetrahedral_fraction(file_path)
        else:
            print(f"Warning: File {file_name} not found.")
    
    # Print results
    for filename, fraction in results.items():
        print(f"{filename}: {fraction:.4f}")

# %%

import pandas as pd
import numpy as np
import os
import re

def calculate_tetrahedral_fraction(file_path, angle_range=np.arange(102.5, 120, 5)):
    """
    Calculates the tetrahedral fraction using the trapezoidal rule
    for angles between 100° and 120° from a given CSV file.
    """
    # Read CSV
    df = pd.read_csv(file_path)
    
    # Extract relevant columns
    valid_columns = [str(angle) for angle in angle_range if str(angle) in df.columns]
    if not valid_columns:
        raise ValueError(f"No valid angle columns found in {file_path}")
    
    # Convert column names to numeric angles
    angles = np.array([float(col) for col in valid_columns])
    values = df[valid_columns].sum().values  # Sum over all rows for each angle bin
    
    # Compute trapezoidal integration
    tetrahedral_fraction = np.trapz(values, angles)
    return tetrahedral_fraction

if __name__ == "__main__":
    # Directory containing CSV files
    data_dir = "/Users/riyanilkant/desktop/hydrophobicity_hbond"
    
    # List of CSV files to process
    file_names = [
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed52.csv"
    ]
    
    # Compute tetrahedral fractions for each file
    results = {}
    for file_name in file_names:
        file_path = os.path.join(data_dir, file_name)
        if os.path.exists(file_path):
            results[file_name] = calculate_tetrahedral_fraction(file_path)
        else:
            print(f"Warning: File {file_name} not found.")
    
    # Compute averages and standard deviations for each residue
    residue_stats = {}
    for file_name, fraction in results.items():
        base_name = re.sub(r'_seed\d+\.csv$', '', file_name)  # Remove seed number to get residue identifier
        if base_name not in residue_stats:
            residue_stats[base_name] = []
        residue_stats[base_name].append(fraction)
    
    # Calculate statistics
    for residue, values in residue_stats.items():
        avg = np.mean(values)
        std = np.std(values)
        print(f"{residue}: Average = {avg:.4f}, Std Dev = {std:.4f}")

# %%

import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt

def calculate_tetrahedral_fraction(file_path, angle_range=np.arange(102.5, 120, 5)):
    """
    Calculates the tetrahedral fraction using the trapezoidal rule
    for angles between 100° and 120° from a given CSV file.
    """
    df = pd.read_csv(file_path)
    valid_columns = [str(angle) for angle in angle_range if str(angle) in df.columns]
    if not valid_columns:
        raise ValueError(f"No valid angle columns found in {file_path}")
    angles = np.array([float(col) for col in valid_columns])
    values = df[valid_columns].sum().values
    return np.trapz(values, angles)

if __name__ == "__main__":
    data_dir = "/Users/riyanilkant/desktop/hydrophobicity_hbond"
    file_names = [
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed52.csv"
    ]
    
    results = {}
    for file_name in file_names:
        file_path = os.path.join(data_dir, file_name)
        if os.path.exists(file_path):
            results[file_name] = calculate_tetrahedral_fraction(file_path)
        else:
            print(f"Warning: File {file_name} not found.")
    
    residue_stats = {}
    for file_name, fraction in results.items():
        base_name = re.sub(r'_and_name_H_triplet_data_seed\d+\.csv$', '', file_name)
        if base_name not in residue_stats:
            residue_stats[base_name] = []
        residue_stats[base_name].append(fraction)
    
    averages = {}
    std_devs = {}
    for residue, values in residue_stats.items():
        averages[residue] = np.mean(values)
        std_devs[residue] = np.std(values)
    
    clusters = [["AQPDZ3_resid_12", "CO2PDZ_resid_23", "CO1PDZ_resid_14", "HSPDZ3_resid_23"],
                ["AQPDZ3_resid_13", "CO2PDZ_resid_24", "CO1PDZ_resid_15", "HSPDZ3_resid_24"],
                ["AQPDZ3_resid_14", "CO2PDZ_resid_25", "CO1PDZ_resid_16", "HSPDZ3_resid_25"]]
    
    fig, ax = plt.subplots()
    width = 0.2
    x = np.arange(len(clusters))
    colors = ['#1f77b4', '#6baed6', '#9ecae1', '#c6dbef']
    
    for j, (category, color) in enumerate(zip(["AQPDZ3", "CO2PDZ", "CO1PDZ", "HSPDZ3"], colors)):
        means = [averages[cluster[j]] for cluster in clusters]
        stds = [std_devs[cluster[j]] for cluster in clusters]
        ax.bar(x + j * width, means, width, yerr=stds, label=category, color=color, capsize=5)
    
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels(["Residue L/F", "Residue G/G", "Resiude F/A"])
    ax.set_ylabel('Tetrahedral Fraction')
    ax.set_ylim(0.16, max(averages.values()) + 0.01)
    ax.legend()
    plt.show()


# %%

import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt

def calculate_tetrahedral_fraction(file_path, angle_range=np.arange(102.5, 120, 5)):
    """
    Calculates the tetrahedral fraction using the trapezoidal rule
    for angles between 100° and 120° from a given CSV file.
    """
    df = pd.read_csv(file_path)
    valid_columns = [str(angle) for angle in angle_range if str(angle) in df.columns]
    if not valid_columns:
        raise ValueError(f"No valid angle columns found in {file_path}")
    angles = np.array([float(col) for col in valid_columns])
    values = df[valid_columns].sum().values
    return np.trapz(values, angles)

if __name__ == "__main__":
    data_dir = "/Users/riyanilkant/desktop/hydrophobicity_hbond"
    file_names = [
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed52.csv"
    ]
    
    results = {}
    for file_name in file_names:
        file_path = os.path.join(data_dir, file_name)
        if os.path.exists(file_path):
            results[file_name] = calculate_tetrahedral_fraction(file_path)
        else:
            print(f"Warning: File {file_name} not found.")
    
    residue_stats = {}
    for file_name, fraction in results.items():
        base_name = re.sub(r'_and_name_H_triplet_data_seed\d+\.csv$', '', file_name)
        if base_name not in residue_stats:
            residue_stats[base_name] = []
        residue_stats[base_name].append(fraction)
    
    averages = {}
    std_devs = {}
    for residue, values in residue_stats.items():
        averages[residue] = np.mean(values)
        std_devs[residue] = np.std(values)
    
    clusters = [["CO1PDZ_resid_14", "CO2PDZ_resid_23", "AQPDZ3_resid_12", "HSPDZ3_resid_23"],
                ["CO1PDZ_resid_15","CO2PDZ_resid_24", "AQPDZ3_resid_13", "HSPDZ3_resid_24"],
                ["CO1PDZ_resid_16", "CO2PDZ_resid_25", "AQPDZ3_resid_14", "HSPDZ3_resid_25"]]
    
    fig, ax = plt.subplots()
    width = 0.2
    x = np.arange(len(clusters))
    colors = ['#1f77b4', '#6baed6', '#9ecae1', '#c6dbef']
    
    for j, (category, color) in enumerate(zip(["CO1PDZ", "CO2PDZ", "AQPDZ3","HSPDZ3"], colors)):
        means = [averages[cluster[j]] for cluster in clusters]
        stds = [std_devs[cluster[j]] for cluster in clusters]
        ax.bar(x + j * width, means, width, yerr=stds, label=category, color=color, capsize=5)
    
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels(["Residue L/F", "Residue G/G", "Resiude F/A"])
    ax.set_ylabel('Average Tetrahedral Water Fraction')
    ax.set_title("Comparison of Tetrahedral Water Fraction for PDZ proteins")
    ax.set_ylim(0.16, max(averages.values()) + 0.01)
    
    legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4, fontsize='small')
    plt.show()

# %%
    
# IN MANUSCRIPT 
    
import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt

def calculate_tetrahedral_fraction(file_path, angle_range=np.arange(102.5, 120, 5)):
    """
    Calculates the tetrahedral fraction using the trapezoidal rule
    for angles between 100° and 120° from a given CSV file.
    """
    df = pd.read_csv(file_path)
    valid_columns = [str(angle) for angle in angle_range if str(angle) in df.columns]
    if not valid_columns:
        raise ValueError(f"No valid angle columns found in {file_path}")
    angles = np.array([float(col) for col in valid_columns])
    values = df[valid_columns].sum().values
    return np.trapz(values, angles)

if __name__ == "__main__":
    data_dir = "/Users/riyanilkant/desktop/hydrophobicity_hbond"
    file_names = [
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_12_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_23_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_13_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_24_and_name_H_triplet_data_seed52.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed50.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed50.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed51.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed51.csv",
        "AQPDZ3_resid_14_and_name_H_triplet_data_seed52.csv", "CO2PDZ_resid_25_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_14_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_23_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_15_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_24_and_name_H_triplet_data_seed52.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed50.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed50.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed51.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed51.csv",
        "CO1PDZ_resid_16_and_name_H_triplet_data_seed52.csv", "HSPDZ3_resid_25_and_name_H_triplet_data_seed52.csv"
    ]
    
    results = {}
    for file_name in file_names:
        file_path = os.path.join(data_dir, file_name)
        if os.path.exists(file_path):
            results[file_name] = calculate_tetrahedral_fraction(file_path)
        else:
            print(f"Warning: File {file_name} not found.")
    
    residue_stats = {}
    for file_name, fraction in results.items():
        base_name = re.sub(r'_and_name_H_triplet_data_seed\d+\.csv$', '', file_name)
        if base_name not in residue_stats:
            residue_stats[base_name] = []
        residue_stats[base_name].append(fraction)
    
    averages = {}
    std_devs = {}
    for residue, values in residue_stats.items():
        averages[residue] = np.mean(values)
        std_devs[residue] = np.std(values)
    
    clusters = [["CO1PDZ_resid_14", "CO2PDZ_resid_23", "AQPDZ3_resid_12", "HSPDZ3_resid_23"],
                ["CO1PDZ_resid_15","CO2PDZ_resid_24", "AQPDZ3_resid_13", "HSPDZ3_resid_24"],
                ["CO1PDZ_resid_16", "CO2PDZ_resid_25", "AQPDZ3_resid_14", "HSPDZ3_resid_25"]]
    
    fig, ax = plt.subplots()
    width = 0.2
    x = np.arange(len(clusters))
    colors = ['#6a0dad', '#8b5cf6', '#b794f4', '#d8b4fe']
    
    for j, (category, color) in enumerate(zip(["CO1PDZ", "CO2PDZ", "AQPDZ3","HSPDZ3"], colors)):
        means = [averages[cluster[j]] for cluster in clusters]
        stds = [std_devs[cluster[j]] for cluster in clusters]
        ax.bar(x + j * width, means, width, yerr=stds, label=category, color=color, capsize=5)
    
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels(["Residue L/F", "Residue G/G", "Resiude F/A"])
    ax.set_ylabel('Tetrahedral Water Fraction')
    #ax.set_title("Comparison of Tetrahedral Water Fraction for PDZ proteins")
    ax.set_ylim(0.16, max(averages.values()) + 0.01)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    
    legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4, fontsize='small')
    plt.show()

# %%
