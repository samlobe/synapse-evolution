#%%
#PDZ ANOVA for each PC

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats

# Load your dataset (adjust file path as necessary)
data = pd.read_csv('PDZs_tree2_embeddings.csv')  # Ensure the file name is a string

# Step 1: Define which columns are embeddings
embedding_columns = [col for col in data.columns if 'embedding_' in col]  # Automatically detects embedding columns

# Step 2: Manually create a 'Group' column based on the 'Unnamed: 0' column or your specific identifier
data['Group'] = data['Unnamed: 0'].apply(lambda x: 'PDZ_1' if 'PDZ_1' in x else 'PDZ_2' if 'PDZ_2' in x else 'PDZ_3')

# Step 3: Standardize the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(data[embedding_columns])

# Step 4: Perform PCA to reduce dimensionality
pca = PCA(n_components=10)  # Adjust number of components based on explained variance
principal_components = pca.fit_transform(scaled_data)

# Step 5: Create a DataFrame for PCA results and add group labels
pca_df = pd.DataFrame(principal_components, columns=[f'PC{i}' for i in range(1, 11)])
pca_df['Group'] = data['Group']  # The newly created 'Group' column

# Step 6: Run ANOVA on the principal components
for i in range(1, 11):  # Test the first 10 PCs
    pc_name = f'PC{i}'
    anova_result = stats.f_oneway(
        pca_df[pca_df['Group'] == 'PDZ_1'][pc_name],
        pca_df[pca_df['Group'] == 'PDZ_2'][pc_name],
        pca_df[pca_df['Group'] == 'PDZ_3'][pc_name]
    )
    print(f"ANOVA result for {pc_name}: {anova_result}")

# %%
#Plot PDZ data for PC2 vs PC7
    
import matplotlib.pyplot as plt
import seaborn as sns

# Create a scatter plot for PC2 vs PC7
plt.figure(figsize=(10, 6))
sns.scatterplot(x=pca_df['PC2'], y=pca_df['PC7'], hue=pca_df['Group'], palette='Set1', s=100)

# Add labels and title
plt.title('Scatter plot of PC2 vs PC7 by PDZ Group')
plt.xlabel('PC2')
plt.ylabel('PC7')
plt.legend(title='Group')
plt.grid(True)
plt.show()

#%%

#Updated ANOVA and post hoc tukey HSD test 

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import numpy as np

# Load your dataset (adjust file path as necessary)
data = pd.read_csv('PDZs_tree2_embeddings.csv')

# Define which columns are embeddings
embedding_columns = [col for col in data.columns if 'embedding_' in col]

# Create a 'Group' column based on the 'Unnamed: 0' column
data['Group'] = data['Unnamed: 0'].apply(lambda x: 'PDZ_1' if 'PDZ_1' in x else 'PDZ_2' if 'PDZ_2' in x else 'PDZ_3')

# Standardize the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(data[embedding_columns])

# Perform PCA to reduce dimensionality
pca = PCA(n_components=10)  # Adjust number of components based on explained variance
principal_components = pca.fit_transform(scaled_data)

# Create a DataFrame for PCA results and add group labels
pca_df = pd.DataFrame(principal_components, columns=[f'PC{i}' for i in range(1, 11)])
pca_df['Group'] = data['Group']

# Run ANOVA on the principal components and perform Tukey's HSD post-hoc test
for i in range(1, 11):  # Test the first 10 PCs
    pc_name = f'PC{i}'
    
    # Run ANOVA
    anova_result = stats.f_oneway(
        pca_df[pca_df['Group'] == 'PDZ_1'][pc_name],
        pca_df[pca_df['Group'] == 'PDZ_2'][pc_name],
        pca_df[pca_df['Group'] == 'PDZ_3'][pc_name]
    )
    print(f"ANOVA result for {pc_name}: {anova_result}")
    
    # If the ANOVA is significant, perform Tukey's HSD test
    if anova_result.pvalue < 0.05:  # Adjust threshold as necessary
        tukey = pairwise_tukeyhsd(endog=pca_df[pc_name], groups=pca_df['Group'], alpha=0.05)
        
        # Print the Tukey table with formatted p-values
        print(f"\nTukey HSD for {pc_name}:")
        for result in tukey.summary().data[1:]:
            group1, group2, meandiff, p_adj, lower, upper, reject = result
            # Ensure very small p-values are formatted in scientific notation
            if p_adj < 1e-10:
                pval_str = f"< 1e-10"  # Handle extremely small p-values that may round to zero
            else:
                pval_str = f"{p_adj:.10f}" if p_adj < 1e-4 else f"{p_adj:.4f}"
            print(f"{group1} vs {group2} -> Meandiff: {meandiff:.4f}, p-adj: {pval_str}, 95% CI: [{lower:.4f}, {upper:.4f}], Reject: {reject}")



#%%
#Tukey HSD for PC2 vs PC7
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# Perform Tukey's HSD test for PC2 and PC7
tukey_pc2 = pairwise_tukeyhsd(endog=pca_df['PC2'], groups=pca_df['Group'], alpha=0.05)
tukey_pc7 = pairwise_tukeyhsd(endog=pca_df['PC7'], groups=pca_df['Group'], alpha=0.05)

# Print the results
print("Tukey HSD Results for PC2:")
print(tukey_pc2)

print("\nTukey HSD Results for PC7:")
print(tukey_pc7)

# %%
# Create the organism type column (Single-celled vs Multi-celled)
single_celled_key = ['SR3', 'MB3', 'CO1', 'CO2', 'SR2', 'MB1', 'SR1']
multi_celled_key = ['BI2', 'BI1', 'BI3', 'ML1', 'ML2', 'HS', 'OM', 'OL', 'AQ', 'CE', 'DM']

# Assuming 'Unnamed: 0' column contains the organism labels (adjust if necessary)
data['Organism_Type'] = data['Unnamed: 0'].apply(lambda x: 'Single-celled' if any(sc in x for sc in single_celled_key) 
                                                 else 'Multi-celled' if any(mc in x for mc in multi_celled_key) 
                                                 else 'Unknown')

# Filter out any 'Unknown' entries if they exist
filtered_data = data[data['Organism_Type'] != 'Unknown']

# Reuse PCA DataFrame (with principal components)
pca_df['Organism_Type'] = filtered_data['Organism_Type']

# %%
from scipy.stats import ttest_ind

# Welch's T-test (default is unequal variance)
for i in range(1, 11):  # Testing the first 10 PCs
    pc_name = f'PC{i}'
    t_test_result = ttest_ind(
        pca_df[pca_df['Organism_Type'] == 'Single-celled'][pc_name],
        pca_df[pca_df['Organism_Type'] == 'Multi-celled'][pc_name],
        equal_var=False  # Welch's t-test
    )
    print(f"T-test result for {pc_name}: {t_test_result}")


# %%
#Visualize data
import matplotlib.pyplot as plt
import seaborn as sns

# PC1 vs PC4
plt.figure(figsize=(10, 6))
sns.scatterplot(x=pca_df['PC1'], y=pca_df['PC4'], hue=pca_df['Organism_Type'], palette='Set1', s=100)
plt.title('Scatter plot of PC1 vs PC4 by Organism Type')
plt.xlabel('PC1')
plt.ylabel('PC4')
plt.legend(title='Organism Type')
plt.grid(True)
plt.show()

# PC1 vs PC6
plt.figure(figsize=(10, 6))
sns.scatterplot(x=pca_df['PC1'], y=pca_df['PC6'], hue=pca_df['Organism_Type'], palette='Set1', s=100)
plt.title('Scatter plot of PC1 vs PC6 by Organism Type')
plt.xlabel('PC1')
plt.ylabel('PC6')
plt.legend(title='Organism Type')
plt.grid(True)
plt.show()

# PC4 vs PC6
plt.figure(figsize=(10, 6))
sns.scatterplot(x=pca_df['PC4'], y=pca_df['PC6'], hue=pca_df['Organism_Type'], palette='Set1', s=100)
plt.title('Scatter plot of PC4 vs PC6 by Organism Type')
plt.xlabel('PC4')
plt.ylabel('PC6')
plt.legend(title='Organism Type')
plt.grid(True)
plt.show()


# %%

import matplotlib.pyplot as plt
import seaborn as sns

# Filter to significant PCs: PC1, PC4, and PC6
significant_pcs = ['PC1', 'PC4', 'PC6']

# Create scatter plots for each pair of significant PCs
for i, pc_x in enumerate(significant_pcs):
    for pc_y in significant_pcs[i+1:]:
        plt.figure(figsize=(10, 6))
        
        # Scatter plot with different colors for Single-celled and Multi-celled
        sns.scatterplot(x=pca_df[pc_x], y=pca_df[pc_y], hue=pca_df['Organism_Type'], palette='Set1', s=100)
        
        # Add labels and title
        plt.title(f'Scatter plot of {pc_x} vs {pc_y} by Organism Type')
        plt.xlabel(pc_x)
        plt.ylabel(pc_y)
        plt.legend(title='Organism Type')
        plt.grid(True)
        
        # Show the plot
        plt.show()

# %%

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
data = pd.read_csv('PDZs_tree2_embeddings.csv')

# Define which columns are embeddings
embedding_columns = [col for col in data.columns if 'embedding_' in col]

# Create a 'Group' column based on the 'Unnamed: 0' column
data['Group'] = data['Unnamed: 0'].apply(lambda x: 'PDZ_1' if 'PDZ_1' in x else 'PDZ_2' if 'PDZ_2' in x else 'PDZ_3')

# Standardize the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(data[embedding_columns])

# Perform PCA to reduce dimensionality to 10 components
pca = PCA(n_components=10)
principal_components = pca.fit_transform(scaled_data)

# Create a DataFrame for PCA results and add group labels
pca_df = pd.DataFrame(principal_components, columns=[f'PC{i}' for i in range(1, 11)])
pca_df['Group'] = data['Group']

# Create the organism type column (Single-celled vs Multi-celled)
single_celled_key = ['SR3', 'MB3', 'CO1', 'CO2', 'SR2', 'MB1', 'SR1']
multi_celled_key = ['BI2', 'BI1', 'BI3', 'ML1', 'ML2', 'HS', 'OM', 'OL', 'AQ', 'CE', 'DM']

# Assuming 'Unnamed: 0' column contains the organism labels (adjust if necessary)
data['Organism_Type'] = data['Unnamed: 0'].apply(lambda x: 'Single-celled' if any(sc in x for sc in single_celled_key) 
                                                 else 'Multi-celled' if any(mc in x for mc in multi_celled_key) 
                                                 else 'Unknown')

# Filter out any 'Unknown' entries if they exist
filtered_data = data[data['Organism_Type'] != 'Unknown']

# Reuse PCA DataFrame (with principal components)
pca_df['Organism_Type'] = filtered_data['Organism_Type']

# Marker styles for single-celled and multi-celled
marker_dict = {'Single-celled': 'o', 'Multi-celled': 's'}

# Create a scatter plot for PC6 vs PC7 with different marker shapes
plt.figure(figsize=(10, 6))
for organism_type, marker in marker_dict.items():
    sns.scatterplot(x=pca_df[pca_df['Organism_Type'] == organism_type]['PC6'],
                    y=pca_df[pca_df['Organism_Type'] == organism_type]['PC7'],
                    hue=pca_df[pca_df['Organism_Type'] == organism_type]['Group'],
                    style=pca_df[pca_df['Organism_Type'] == organism_type]['Organism_Type'],
                    markers=marker,
                    s=100, palette='Set1', legend='full')

# Add labels and title
plt.title('Scatter plot of PC6 vs PC7 by PDZ Group and Organism Type', fontsize=16)
plt.xlabel('PC6', fontsize=14)
plt.ylabel('PC7', fontsize=14)
plt.legend(title='Group and Organism Type', fontsize=12)
plt.grid(True)
plt.show()

# %%
