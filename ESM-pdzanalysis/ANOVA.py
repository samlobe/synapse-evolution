#%%

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import numpy as np
from scipy.stats import ttest_ind

#%%
# ANOVA and post hoc tukey HSD test 

data = pd.read_csv('PDZs_tree2_embeddings.csv')

# Define which columns are embeddings
embedding_columns = [col for col in data.columns if 'embedding_' in col]

# Create a 'Group' column based on the 'Unnamed: 0' column
data['Group'] = data['Unnamed: 0'].apply(lambda x: 'PDZ_1' if 'PDZ_1' in x else 'PDZ_2' if 'PDZ_2' in x else 'PDZ_3')

scaler = StandardScaler()
scaled_data = scaler.fit_transform(data[embedding_columns])

# Perform PCA 
pca = PCA(n_components=10)  
principal_components = pca.fit_transform(scaled_data)


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
    if anova_result.pvalue < 0.05: 
        tukey = pairwise_tukeyhsd(endog=pca_df[pc_name], groups=pca_df['Group'], alpha=0.05)
        
        # Print the Tukey table with formatted p-values
        print(f"\nTukey HSD for {pc_name}:")
        for result in tukey.summary().data[1:]:
            group1, group2, meandiff, p_adj, lower, upper, reject = result
            # small p-values are formatted in scientific notation
            if p_adj < 1e-10:
                pval_str = f"< 1e-10"  
            else:
                pval_str = f"{p_adj:.10f}" if p_adj < 1e-4 else f"{p_adj:.4f}"
            print(f"{group1} vs {group2} -> Meandiff: {meandiff:.4f}, p-adj: {pval_str}, 95% CI: [{lower:.4f}, {upper:.4f}], Reject: {reject}")


# %%
#T test for single vs multi celled organisms 
            
# Create the organism type column (Single-celled vs Multi-celled)
single_celled_key = ['SR3', 'MB3', 'CO1', 'CO2', 'SR2', 'MB1', 'SR1']
multi_celled_key = ['BI2', 'BI1', 'BI3', 'ML1', 'ML2', 'HS', 'OM', 'OL', 'AQ', 'CE', 'DM']

# Assuming 'Unnamed: 0' column contains the organism labels 
data['Organism_Type'] = data['Unnamed: 0'].apply(lambda x: 'Single-celled' if any(sc in x for sc in single_celled_key) 
                                                 else 'Multi-celled' if any(mc in x for mc in multi_celled_key) 
                                                 else 'Unknown')

# Filter out any 'Unknown' entries if they exist
filtered_data = data[data['Organism_Type'] != 'Unknown']

# Reuse PCA DataFrame (with principal components)
pca_df['Organism_Type'] = filtered_data['Organism_Type']

# Welch's T-test (default is unequal variance)
for i in range(1, 11):  
    pc_name = f'PC{i}'
    t_test_result = ttest_ind(
        pca_df[pca_df['Organism_Type'] == 'Single-celled'][pc_name],
        pca_df[pca_df['Organism_Type'] == 'Multi-celled'][pc_name],
        equal_var=False 
    )
    print(f"T-test result for {pc_name}: {t_test_result}")


