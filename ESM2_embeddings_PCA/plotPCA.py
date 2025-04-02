#%%
# Read the CSV 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from adjustText import adjust_text
import matplotlib.lines as mlines


# read the embeddings
embeddings_df = pd.read_csv('PDZs_tree2_embeddings.csv',index_col=0)

# PCA
pca = PCA(n_components=5)
pca.fit(embeddings_df)
embeddings_pca = pca.transform(embeddings_df)



#%%
# plot the cumulative variance

plt.plot(np.arange(5)+1,pca.explained_variance_ratio_.cumsum()*100)
plt.xlabel('Principal component',fontsize=15)
plt.ylabel('Cumulative explained variance',fontsize=15)
# xticks should be 1 indexed
plt.xticks(np.arange(1,6),fontsize=12)
# display the yticks as percentages with no decimal places
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(decimals=1))



#%% 
# table of variance

# PCA with 10 components
pca_10 = PCA(n_components=10)
pca_10.fit(embeddings_df)

# Get the explained variance for each PC (not cumulative)
explained_variance = pca_10.explained_variance_ratio_ * 100

# Create a table of variances for PCs 1 through 10
variance_table = pd.DataFrame({
    'Principal Component': [f'PC{i+1}' for i in range(10)],
    'Explained Variance (%)': explained_variance
})

# Output the table
print(variance_table)




#%%
# Plot the complete PCA with all points and labels 

indices1 = embeddings_df.index[embeddings_df.index.str.endswith('PDZ_1') | embeddings_df.index.str.endswith('PDZ1')].tolist()
indices2 = embeddings_df.index[embeddings_df.index.str.endswith('PDZ_2') | embeddings_df.index.str.endswith('PDZ2')].tolist()
indices3 = embeddings_df.index[embeddings_df.index.str.endswith('PDZ_3') | embeddings_df.index.str.endswith('PDZ3')].tolist()
indices4 = embeddings_df.index[embeddings_df.index.str.endswith('PDZ_0') | embeddings_df.index.str.endswith('PDZ0')].tolist()

# Convert index names to integer indices
int_indices1 = [embeddings_df.index.get_loc(name) for name in indices1]
int_indices2 = [embeddings_df.index.get_loc(name) for name in indices2]
int_indices3 = [embeddings_df.index.get_loc(name) for name in indices3]
int_indices4 = [embeddings_df.index.get_loc(name) for name in indices4]


# Plotting the PCA results
plt.scatter(embeddings_pca[int_indices1,0], embeddings_pca[int_indices1,1], color='red', label='PDZ-1')
plt.scatter(embeddings_pca[int_indices2,0], embeddings_pca[int_indices2,1], color='blue', label='PDZ-2')
plt.scatter(embeddings_pca[int_indices3,0], embeddings_pca[int_indices3,1], color='green', label='PDZ-3')
plt.scatter(embeddings_pca[int_indices4,0], embeddings_pca[int_indices4,1], color='grey', label='PDZ-0')


# simplify the PDZ names - look for the first underscore and take the substring before it
names1 = [name[:name.find('_')] for name in indices1]
names2 = [name[:name.find('_')] for name in indices2]
names3 = [name[:name.find('_')] for name in indices3]
names4 = [name[:name.find('_')] for name in indices4]

# Adding annotations and collecting them into 'texts'
from adjustText import adjust_text
texts = []
for i, txt in enumerate(names1):
    texts.append(plt.text(embeddings_pca[int_indices1[i],0], embeddings_pca[int_indices1[i],1], txt, fontsize=12, ha='center',color='red'))
for i, txt in enumerate(names2):
    texts.append(plt.text(embeddings_pca[int_indices2[i],0], embeddings_pca[int_indices2[i],1], txt, fontsize=12, ha='center',color='blue'))
for i, txt in enumerate(names3):
    texts.append(plt.text(embeddings_pca[int_indices3[i],0], embeddings_pca[int_indices3[i],1], txt, fontsize=12, ha='center',color='green'))
for i, txt in enumerate(names4):
    texts.append(plt.text(embeddings_pca[int_indices4[i],0], embeddings_pca[int_indices4[i],1], txt, fontsize=12, ha='center',color='grey'))

adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'))


plt.xlabel('PC1', fontsize=15)
plt.ylabel('PC2', fontsize=15)
plt.title("PCA of evolutionary PDZs' ESM embeddings", fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend()
plt.show()



#%%

#Cleaned up PCA with markers and shapes

# Define the categories based on your key
single_celled_key = ['SR3', 'MB3', 'CO1', 'CO2', 'SR2', 'MB1', 'SR1']
multi_celled_key = ['BI2', 'BI1', 'BI3', 'ML1', 'ML2', 'HS', 'OM', 'OL', 'AQ', 'CE', 'DM']
label_keys = ['CO', 'HS', 'HS1', 'AQ']

def get_shape_and_color(index):
    if any(index.startswith(key) for key in single_celled_key):
        shape = '^'  # Triangle for single-celled
    else:
        shape = 'o'  # Circle for multi-celled

    if index.endswith('PDZ_1') or index.endswith('PDZ1'):
        color = 'red'
    elif index.endswith('PDZ_2') or index.endswith('PDZ2'):
        color = 'blue'
    elif index.endswith('PDZ_3') or index.endswith('PDZ3'):
        color = 'green'
    elif index.endswith('PDZ_0') or index.endswith('PDZ0'):
        color = 'grey'
    else:
        color = 'black'  # Default color if none of the PDZ types match

    return shape, color

# Plotting all data points with shape and color distinction
texts = []
for index in embeddings_df.index:
    shape, color = get_shape_and_color(index)
    int_index = embeddings_df.index.get_loc(index)
    plt.scatter(embeddings_pca[int_index, 0], embeddings_pca[int_index, 1], color=color, marker=shape)
    
    # Add labels for specified categories, truncated before the first "_"
    if any(index.startswith(key) for key in label_keys):
        label = index.split('_')[0]  # Truncate label before the first "_"
        if label == 'AQ':
            texts.append(plt.text(embeddings_pca[int_index, 0] - 0.05, embeddings_pca[int_index, 1] - 0.35, label, fontsize=10, ha='center', color=color))
        else:
            texts.append(plt.text(embeddings_pca[int_index, 0] - 0.05, embeddings_pca[int_index, 1], label, fontsize=10, ha='center', color=color))

# Adjust texts to minimize overlaps and include lines connecting labels to points
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'))

# Custom legend handles
single_celled_handle = mlines.Line2D([], [], color='black', marker='^', markersize=8, linestyle='None', label='Single-celled')
multi_celled_handle = mlines.Line2D([], [], color='black', marker='o', markersize=8, linestyle='None', label='Multi-celled')
pdz1_handle = mlines.Line2D([], [], color='red', marker='o', markersize=8, linestyle='None', label='PDZ-1')
pdz2_handle = mlines.Line2D([], [], color='blue', marker='o', markersize=8, linestyle='None', label='PDZ-2')
pdz3_handle = mlines.Line2D([], [], color='green', marker='o', markersize=8, linestyle='None', label='PDZ-3')
pdz0_handle = mlines.Line2D([], [], color='grey', marker='o', markersize=8, linestyle='None', label='PDZ-0')

# Adjust legend position and size
plt.legend(handles=[single_celled_handle, multi_celled_handle, pdz1_handle, pdz2_handle, pdz3_handle, pdz0_handle],
           loc='upper right', fontsize=8, bbox_to_anchor=(0.98, 0.98), borderaxespad=0.1)

plt.xlabel('PC1', fontsize=15)
plt.ylabel('PC2', fontsize=15)
plt.title("PCA of evolutionary PDZs' ESM embeddings", fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()


