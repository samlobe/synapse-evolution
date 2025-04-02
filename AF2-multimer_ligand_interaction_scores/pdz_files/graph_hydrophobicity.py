

# %%
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Labels
categories = ["CO1", "CO2", "AQ", "HS"]

# Mean values
means = [0.248239815, 0.239948711, 0.243624887, 0.236148901]

# Standard deviations
std_devs = [0.002286059, 0.000627978, 0.006634613, 0.003865587]

# Define bar positions
x = np.arange(len(categories))

# Create figure
fig, ax = plt.subplots(figsize=(8, 6))

# Aesthetic adjustments
sns.set(style="white")
plt.figure(figsize=(8, 5))
colors = sns.color_palette("Blues", len(x))  # Professional colors

# Plot bars with error bars
bars = ax.bar(x, means, yerr=std_devs, capsize=5, color=colors, width=0.6)

# Add labels on bars
for bar, mean in zip(bars, means):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.009, 
            f"{mean:.3f}", ha='center', fontsize=8, fontweight='bold')

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=14, fontweight='bold')
ax.set_ylabel("Average Tetrahedral Water Fraction of GLGF Residues", fontsize=14, fontweight='bold')
ax.set_title("Comparison of Tetrahedral Water Fraction for PDZ Proteins", fontsize=16, fontweight='bold')

# Light grid lines
ax.yaxis.grid(True, linestyle='--', alpha=0.5)

# Show the plot
plt.show()

# %%
