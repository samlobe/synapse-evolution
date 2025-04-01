#%%
# Read the CSV 
#!/Users/riyanilkant/opt/anaconda3/envs/openmm/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


# read the embeddings
embeddings_df = pd.read_csv('PDZs_tree2_embeddings.csv',index_col=0)

# PCA
pca = PCA(n_components=5)
pca.fit(embeddings_df)
embeddings_pca = pca.transform(embeddings_df)

# plot
# plt.scatter(embeddings_pca[:,0],embeddings_pca[:,1])
# plt.xlabel('PC1',fontsize=15)
# plt.ylabel('PC2',fontsize=15)
# plt.title("PCA of evolutionary PDZs' ESM embeddings",fontsize=15)
# plt.tick_params(axis='both', which='major', labelsize=12)
# plt.show()




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

#%% table of variance

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
# if the indices begin with "Node" then put them in a separate list
# indices_nodes = embeddings_df.index[embeddings_df.index.str.startswith('Node')].tolist()

# Convert index names to integer indices
int_indices1 = [embeddings_df.index.get_loc(name) for name in indices1]
int_indices2 = [embeddings_df.index.get_loc(name) for name in indices2]
int_indices3 = [embeddings_df.index.get_loc(name) for name in indices3]
int_indices4 = [embeddings_df.index.get_loc(name) for name in indices4]

# int_nodes = [embeddings_df.index.get_loc(name) for name in indices_nodes]   

# Plotting the PCA results
plt.scatter(embeddings_pca[int_indices1,0], embeddings_pca[int_indices1,1], color='red', label='PDZ-1')
plt.scatter(embeddings_pca[int_indices2,0], embeddings_pca[int_indices2,1], color='blue', label='PDZ-2')
plt.scatter(embeddings_pca[int_indices3,0], embeddings_pca[int_indices3,1], color='green', label='PDZ-3')
plt.scatter(embeddings_pca[int_indices4,0], embeddings_pca[int_indices4,1], color='grey', label='PDZ-0')

# plt.scatter(embeddings_pca[int_nodes,0], embeddings_pca[int_nodes,1], color='gray', label='Node')

# simplify the PDZ names - look for the first underscore and take the substring before it
names1 = [name[:name.find('_')] for name in indices1]
names2 = [name[:name.find('_')] for name in indices2]
names3 = [name[:name.find('_')] for name in indices3]
names4 = [name[:name.find('_')] for name in indices4]

# names_nodes = [name for name in indices_nodes]

# annotate the points with the PDZ names
# for i, txt in enumerate(names1):
#     plt.annotate(txt, (embeddings_pca[int_indices1[i],0], embeddings_pca[int_indices1[i],1]), fontsize=12, ha='center')
# for i, txt in enumerate(names2):
#     plt.annotate(txt, (embeddings_pca[int_indices2[i],0], embeddings_pca[int_indices2[i],1]), fontsize=12, ha='center')
# for i, txt in enumerate(names3): # center aligned
#     plt.annotate(txt, (embeddings_pca[int_indices3[i],0], embeddings_pca[int_indices3[i],1]), fontsize=12, ha='center')

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
# for i, txt in enumerate(names_nodes):
    # texts.append(plt.text(embeddings_pca[int_nodes[i],0], embeddings_pca[int_nodes[i],1], txt, fontsize=12, ha='center',color='gray'))

# Adjust texts to minimize overlaps
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'))


plt.xlabel('PC1', fontsize=15)
plt.ylabel('PC2', fontsize=15)
plt.title("PCA of evolutionary PDZs' ESM embeddings", fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend()
plt.show()

#%%

#Cleaned up PCA with markers and shapes

import matplotlib.pyplot as plt
from adjustText import adjust_text
import matplotlib.lines as mlines

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


#%%

#Cleaned up PDZ3 in pink with markers ahd shapes

import numpy as np
from scipy.stats import linregress
from adjustText import adjust_text
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# Example definitions (replace with your actual data)
# embeddings_pca = np.random.rand(100, 2)  # Replace with actual PCA embeddings
# int_indices_pdz3 = np.random.choice(range(100), 20, replace=False)  # Replace with actual PDZ3 indices

# Define categories based on your key
single_celled_key = ['SR3', 'MB3', 'CO1', 'CO2', 'SR2', 'MB1', 'SR1']
multi_celled_key = ['BI2', 'BI1', 'BI3', 'ML1', 'ML2', 'HS', 'OM', 'OL', 'AQ', 'CE', 'DM']
label_keys = ['HS', 'AQ', 'HS1', 'CO1', 'CO2']

def get_shape_and_color(index):
    if any(index.startswith(key) for key in single_celled_key):
        shape = '^'  # Triangle for single-celled
        color = 'purple'
    else:
        shape = 'o'  # Circle for multi-celled
        color = 'fuchsia'  # Bright fuchsia pink

    return shape, color

# Combine PDZ-3 data points
combined_data = embeddings_pca[int_indices_pdz3, :]

# Assuming `indices_pdz3` is correctly defined and matches the order of `int_indices_pdz3`
# If not, you may need to ensure `indices_pdz3` is in the correct order or mapping

# Perform linear regression on the combined data
slope, intercept, r_value, p_value, std_err = linregress(combined_data[:, 0], combined_data[:, 1])

# Plot PDZ-3 points with different shapes and colors
texts_pdz3 = []
for i, index in enumerate(indices_pdz3):
    shape, color = get_shape_and_color(index)
    
    plt.scatter(combined_data[i, 0], combined_data[i, 1], color=color, marker=shape)

    # Add labels for specified categories, truncated before the first "_"
    if any(index.startswith(key) for key in label_keys):
        label = index.split('_')[0]  # Truncate label before the first "_"
        texts_pdz3.append(plt.text(combined_data[i, 0], combined_data[i, 1] - 0.35, label, fontsize=10, ha='center', color=color))

# Adjust texts to minimize overlaps for PDZ-3
adjust_text(texts_pdz3, arrowprops=dict(arrowstyle='->', color='black'))

# Plot linear regression line
line, = plt.plot(combined_data[:, 0], slope * combined_data[:, 0] + intercept, color='black', label='Linear Regression')

# Add arrow at the end of the linear regression line
x_start, x_end = plt.xlim()
y_start, y_end = slope * x_start + intercept, slope * x_end + intercept
plt.annotate('', xy=(x_end, y_end), xytext=(x_start, y_start),
             arrowprops=dict(arrowstyle='->', color='black'))

# Custom legend handles
single_celled_handle = mlines.Line2D([], [], color='purple', marker='^', markersize=10, linestyle='None', label='Single-celled')
multi_celled_handle = mlines.Line2D([], [], color='fuchsia', marker='o', markersize=10, linestyle='None', label='Multi-celled')
line_handle = mlines.Line2D([], [], color='black', linestyle='-', label='Linear Regression')

# Create the custom legend without PDZ-3
plt.legend(handles=[single_celled_handle, multi_celled_handle, line_handle], loc='upper right', bbox_to_anchor=(1, 1))

# Adjust the layout to make space for the legend
plt.subplots_adjust(right=0.75)  # Adjust as necessary

plt.xlabel('PC1', fontsize=15)
plt.ylabel('PC2', fontsize=15)
plt.title("PCA of evolutionary PDZ3 ESM embeddings with Linear Regression", fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=12)

plt.show()

#%%

# %%
# UMAP plots


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from adjustText import adjust_text

# Read the CSV
embeddings_df = pd.read_csv('PDZs_tree2_embeddings.csv', index_col=0)

# Scale the data
scaler = StandardScaler()
embeddings_scaled = scaler.fit_transform(embeddings_df)

# PCA
pca = PCA(n_components=5)
embeddings_pca = pca.fit_transform(embeddings_scaled)

# UMAP on PCA-reduced data
umap_model = umap.UMAP(n_components=2, random_state=42)
embeddings_umap = umap_model.fit_transform(embeddings_pca)

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
plt.figure(figsize=(12, 10))
texts = []

for index in embeddings_df.index:
    shape, color = get_shape_and_color(index)
    int_index = embeddings_df.index.get_loc(index)
    plt.scatter(embeddings_umap[int_index, 0], embeddings_umap[int_index, 1], color=color, marker=shape, label=color if index == embeddings_df.index[0] else "")
    
    # Add labels for specified categories, truncated before the first "_"
    if any(index.startswith(key) for key in label_keys):
        label = index.split('_')[0]  # Truncate label before the first "_"
        if label == 'AQ':
            texts.append(plt.text(embeddings_umap[int_index, 0] - 0.05, embeddings_umap[int_index, 1] - 0.35, label, fontsize=10, ha='center', color=color))
        else:
            texts.append(plt.text(embeddings_umap[int_index, 0] - 0.05, embeddings_umap[int_index, 1], label, fontsize=10, ha='center', color=color))

# Adjust texts to minimize overlaps and include lines connecting labels to points
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'))

# Custom legend handles
single_celled_handle = plt.Line2D([], [], color='black', marker='^', markersize=8, linestyle='None', label='Single-celled')
multi_celled_handle = plt.Line2D([], [], color='black', marker='o', markersize=8, linestyle='None', label='Multi-celled')
pdz1_handle = plt.Line2D([], [], color='red', marker='o', markersize=8, linestyle='None', label='PDZ-1')
pdz2_handle = plt.Line2D([], [], color='blue', marker='o', markersize=8, linestyle='None', label='PDZ-2')
pdz3_handle = plt.Line2D([], [], color='green', marker='o', markersize=8, linestyle='None', label='PDZ-3')
pdz0_handle = plt.Line2D([], [], color='grey', marker='o', markersize=8, linestyle='None', label='PDZ-0')

# Adjust legend position and size
plt.legend(handles=[single_celled_handle, multi_celled_handle, pdz1_handle, pdz2_handle, pdz3_handle, pdz0_handle],
           loc='upper right', fontsize=8, bbox_to_anchor=(0.98, 0.98), borderaxespad=0.1)

plt.xlabel('UMAP1', fontsize=15)
plt.ylabel('UMAP2', fontsize=15)
plt.title("UMAP projection of evolutionary PDZs' ESM embeddings", fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.show()

# %%

#Final confidence ellipses graphs 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA

# Function to plot an ellipse
def plot_ellipse(mean, cov, ax, color, label, alpha=0.2):
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = eigenvalues.argsort()[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    theta = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
    width, height = 2 * np.sqrt(eigenvalues)
    ellipse = Ellipse(xy=mean, width=width, height=height, angle=theta, color=color, alpha=alpha, label=label)
    ax.add_patch(ellipse)

# Your existing data and PCA transformation
embeddings_df = pd.read_csv('PDZs_tree2_embeddings.csv', index_col=0)
pca = PCA(n_components=2)
embeddings_pca = pca.fit_transform(embeddings_df)

# Define your PDZ categories
categories = {
    'PDZ-1': {'indices': embeddings_df.index.str.endswith(('PDZ_1', 'PDZ1')), 'color': 'red'},
    'PDZ-2': {'indices': embeddings_df.index.str.endswith(('PDZ_2', 'PDZ2')), 'color': 'blue'},
    'PDZ-3': {'indices': embeddings_df.index.str.endswith(('PDZ_3', 'PDZ3')), 'color': 'green'},
    #'PDZ-0': {'indices': embeddings_df.index.str.endswith(('PDZ_0', 'PDZ0')), 'color': 'brown'}
}

# Plotting
fig, ax = plt.subplots(figsize=(10, 8))

for label, props in categories.items():
    indices = props['indices']
    color = props['color']
    points = embeddings_pca[indices]
    
    # Scatter all points
    ax.scatter(points[:, 0], points[:, 1], color=color, label=label)
    
    if points.shape[0] > 1:  # Only calculate mean and covariance if there are enough points
        # Calculate the mean and covariance for the ellipse
        mean = np.mean(points, axis=0)
        cov = np.cov(points, rowvar=False)
        
        # Plot the ellipse for each group
        plot_ellipse(mean, cov, ax, color, label)
    else:
        print(f"Skipping {label} because it has less than 2 points.")

# Customize the plot
ax.set_xlabel('PC1', fontsize=15)
ax.set_ylabel('PC2', fontsize=15)
ax.set_title("PCA of evolutionary PDZs' ESM embeddings with Confidence Ellipses", fontsize=15)
ax.legend()

plt.show()

#%%

#%% CRIPT SCATTERPLOT 
import matplotlib.pyplot as plt

# Data
data = {
    'AQ': 0.65676222,
    'CO1': 0.621128631,
    'CO2': 0.675571253,
    'HS': 0.68269705,
    #'NODE3 x COCRIPT': 0.454594356,
    #'NODE3 x HS CRIPT': 0.450946341,
    #'NODE8 x COCRIPT': 0.476117407,
    #'NODE8 x HSCRIPT': 0.502145215,
    #'NODE11 x CO CRIPT': 0.532353828,
    #'NODE11 x HS CRIPT': 0.543203769
}

# Scatter plot
plt.figure(figsize=(10, 6))
for i, (key, value) in enumerate(data.items()):
    plt.scatter(value, 0, label=key, color=plt.cm.tab10(i))
plt.xlabel('LIS score')
plt.yticks([])
plt.legend(title='Data points', loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('LIS Scores of AF2 Simulations with CRIPT')
plt.show()

#%%

#Other PCs



# %%
