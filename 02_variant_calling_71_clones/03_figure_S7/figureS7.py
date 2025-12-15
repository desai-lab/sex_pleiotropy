import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

# Load the data
file_path = 'GGE_clones_final.tsv'
data = pd.read_csv(file_path, sep='\t')

# Define clones of interest and their regimes
clones_of_interest = {
    'D04_A1.1_H06_1': 'Asexual',
    'D04_A1.1_H06_2': 'Asexual',
    'C04_A3_H02_1': 'Sexual',
    'C04_A3_H02_2': 'Sexual',
    'C02_A3_A03_1': 'Sexual',
    'C02_A3_A03_2': 'Sexual'
}

# Define colors for regimes
regime_colors = {'Asexual': '#4682B4', 'Sexual': '#EE7600'}

# Create figure with GridSpec for flexible subplot arrangement
fig = plt.figure(figsize=(18, 12))
gs = gridspec.GridSpec(3, 3, height_ratios=[1, 1, 1], width_ratios=[1, 1, 1])

# ---------- DEFINE COMMON BIN EDGES ----------
num_bins = 20
bin_edges = np.linspace(0, 1, num_bins + 1)  # Ensures equal bin widths from 0 to 1

# ---------- GLOBAL HISTOGRAM (Top - Centered in Middle Column) ----------
all_clones_alt_frequencies = []

# Identify all alt_counts and ref_counts columns in the dataset
alt_columns_all = [col for col in data.columns if col.endswith('_alt_counts')]
ref_columns_all = [col.replace('_alt_counts', '_ref_counts') for col in alt_columns_all]

# Compute alt allele frequencies for all clones
for alt_col, ref_col in zip(alt_columns_all, ref_columns_all):
    valid_rows = data[alt_col] > 0
    alt_freq = data.loc[valid_rows, alt_col] / (data.loc[valid_rows, alt_col] + data.loc[valid_rows, ref_col])
    all_clones_alt_frequencies.extend(alt_freq)

# Compute the average count per clone in each bin
hist_values, _ = np.histogram(all_clones_alt_frequencies, bins=bin_edges)
average_counts_per_clone = hist_values / len(alt_columns_all)  # Normalize by number of clones

# Create a single subplot for the global histogram in the middle column
global_ax = fig.add_subplot(gs[0, 1])
global_ax.bar(bin_edges[:-1], average_counts_per_clone, width=np.diff(bin_edges), color='grey', edgecolor='black', alpha=0.7, align='edge')
global_ax.set_title('All Clones', fontsize=14)
global_ax.set_xlabel('Alternate Allele Frequency')
global_ax.set_ylabel('Avg Count per Clone')
global_ax.grid(axis='y', linestyle='--', alpha=0.7)
global_ax.set_xlim(0, 1)  # Shared x-axis

# Empty spaces in top left and right to balance layout
fig.add_subplot(gs[0, 0]).axis('off')
fig.add_subplot(gs[0, 2]).axis('off')

# ---------- INDIVIDUAL CLONE HISTOGRAMS ----------
ordered_clones = [
    ('D04_A1.1_H06_1', 'D04_A1.1_H06_2'),
    ('C04_A3_H02_1', 'C04_A3_H02_2'),
    ('C02_A3_A03_1', 'C02_A3_A03_2')
]

all_clone_frequencies = []
for col, (clone_1, clone_2) in enumerate(ordered_clones):
    for row, clone in enumerate([clone_1, clone_2]):
        regime = clones_of_interest[clone]

        # Identify alt_counts and ref_counts columns **only for this specific clone**
        alt_columns_of_interest = [column for column in data.columns if clone in column and column.endswith('_alt_counts')]
        ref_columns_of_interest = [column.replace('_alt_counts', '_ref_counts') for column in alt_columns_of_interest]

        # Calculate alt frequencies for the specified clone
        alt_frequencies = []
        for alt_col, ref_col in zip(alt_columns_of_interest, ref_columns_of_interest):
            valid_rows = data[alt_col] > 0
            alt_freq = data.loc[valid_rows, alt_col] / (data.loc[valid_rows, alt_col] + data.loc[valid_rows, ref_col])
            alt_frequencies.extend(alt_freq)

        # Store frequencies for setting shared limits
        all_clone_frequencies.extend(alt_frequencies)

# Compute shared y-axis limits
y_min = 0
hist_values_all, _ = np.histogram(all_clone_frequencies, bins=bin_edges)
y_max = max(hist_values_all) * 1.3

# Re-plot individual histograms with shared x-axis and y-axis limits and updated titles
for col, (clone_1, clone_2) in enumerate(ordered_clones):
    for row, clone in enumerate([clone_1, clone_2]):
        regime = clones_of_interest[clone]

        # Assign fixed population numbers
        if regime == "Asexual":
            population_label = "Asexual Population 1"
        else:
            population_label = f"Sexual Population {col}"  # Assigns 1 to first column, 2 to second

        # Extract the clone number (1 or 2) from the end of the name
        clone_num = clone.split('_')[-1]

        # Format the title
        title = f"{population_label} Clone {clone_num}"

        # Identify alt_counts and ref_counts columns **only for this specific clone**
        alt_columns_of_interest = [column for column in data.columns if clone in column and column.endswith('_alt_counts')]
        ref_columns_of_interest = [column.replace('_alt_counts', '_ref_counts') for column in alt_columns_of_interest]

        # Calculate alt frequencies for the specified clone
        alt_frequencies = []
        for alt_col, ref_col in zip(alt_columns_of_interest, ref_columns_of_interest):
            valid_rows = data[alt_col] > 0
            alt_freq = data.loc[valid_rows, alt_col] / (data.loc[valid_rows, alt_col] + data.loc[valid_rows, ref_col])
            alt_frequencies.extend(alt_freq)

        # Assign subplot position using gridspec
        ax = fig.add_subplot(gs[row + 1, col])
        ax.hist(alt_frequencies, bins=bin_edges, color=regime_colors[regime], edgecolor='black', alpha=0.7)
        ax.set_title(title, fontsize=12)  # Set updated title
        ax.set_xlabel('Alternate Allele Frequency')
        ax.set_ylabel('Count')
        ax.grid(axis='y', linestyle='--', alpha=0.7)

        # Apply shared limits
        ax.set_xlim(0, 1)
        ax.set_ylim(y_min, y_max)

        # Remove x labels from the **upper row of clone histograms**
        if row == 0:
            ax.set_xlabel(None)
        # Remove y labels
        if col == 1 or col == 2:
            ax.set_ylabel(None)

# Also apply the shared y-axis limit to the global plot
global_ax.set_ylim(y_min, y_max)

# Add label 'a' for the global histogram
fig.text(0.33, 0.92, 'A', fontsize=30, fontweight='bold', ha='center', va='center')

# Add label 'b' for the rest of the plot (placed near the clone histograms)
fig.text(0.02, 0.62, 'B', fontsize=30, fontweight='bold', ha='center', va='center')

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.95])

# Show the final figure with all histograms
plt.savefig("figureS7.pdf", dpi=300)
plt.savefig("figureS7.jpg", dpi=300)
#plt.show()