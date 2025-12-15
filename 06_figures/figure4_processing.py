import pandas as pd
import re

## Figure out which genes are multi-hit or singletons

# Load the data
file_path = 'GGE_clones_final.tsv'  # Replace with your actual file path
data = pd.read_csv(file_path, sep='\t')

# Retain only `alt_counts` columns and `ANN_simpler` for mutation details
alt_count_columns = [col for col in data.columns if '_alt_counts' in col]
filtered_data = data[alt_count_columns + ['ANN_simpler']].copy()

# Extract the relevant mutation category from `ANN_simpler`
mutation_categories = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
filtered_data['Parsed_Category'] = filtered_data['ANN_simpler'].apply(
    lambda x: next((cat for cat in mutation_categories if cat in str(x)), 'unknown')
)

# Exclude synonymous and noncoding mutations
excluded_categories = ['synonymous', 'noncoding']
filtered_data = filtered_data[~filtered_data['Parsed_Category'].isin(excluded_categories)]

# Convert alt_counts columns to binary (0 if 0, 1 if nonzero)
filtered_data[alt_count_columns] = filtered_data[alt_count_columns].map(lambda x: 1 if x > 0 else 0)

# Identify populations by removing the clone suffix (e.g., _1, _2)
population_columns = [col.replace('_alt_counts', '') for col in alt_count_columns]
populations = list(set('_'.join(col.rsplit('_', 2)[:-2]) for col in population_columns))  # Unique populations

# Map alt_count columns to populations
pop_map = {}
for col in alt_count_columns:
    population = '_'.join(col.replace('_alt_counts', '').split('_')[:-2])
    if population not in pop_map:
        pop_map[population] = []
    pop_map[population].append(col)

# Create a dictionary to track genes per population
gene_population_map = {}

for idx, row in filtered_data.iterrows():
    # Extract the gene name from the 4th field of ANN_simpler
    gene = str(row['ANN_simpler']).split('|')[3] if len(str(row['ANN_simpler']).split('|')) > 3 else None
    if gene:  # Proceed only if a gene name is extracted
        for population, clones in pop_map.items():
            # If at least one clone in the population has the mutation
            if row[clones].sum() > 0:
                if gene not in gene_population_map:
                    gene_population_map[gene] = set()
                gene_population_map[gene].add(population)

# Separate genes into multi-hit and single-hit
multi_hit_genes = [gene for gene, pops in gene_population_map.items() if len(pops) > 1]
single_hit_genes = [gene for gene, pops in gene_population_map.items() if len(pops) == 1]

# Output results
print(f"Number of multi-hit genes: {len(multi_hit_genes)}")
print(f"Number of single-hit genes: {len(single_hit_genes)}")

# # Save the results to CSV files
# pd.DataFrame({'Multi-Hit Genes': multi_hit_genes}).to_csv('multi_hit_genes.csv', index=False)
# pd.DataFrame({'Single-Hit Genes': single_hit_genes}).to_csv('single_hit_genes.csv', index=False)

## For each clone, how many mutations are in multi-hit genes or singletons

# Step 1: Filter mutations corresponding to multi-hit and single-hit genes
multi_hit_genes_set = set(multi_hit_genes)
single_hit_genes_set = set(single_hit_genes)

# Create masks for multi-hit and single-hit mutations
multi_hit_mask = filtered_data['ANN_simpler'].apply(lambda x: str(x).split('|')[3] if len(str(x).split('|')) > 3 else None).isin(multi_hit_genes_set)
single_hit_mask = filtered_data['ANN_simpler'].apply(lambda x: str(x).split('|')[3] if len(str(x).split('|')) > 3 else None).isin(single_hit_genes_set)

# Filter data for multi-hit and single-hit mutations
multi_hit_mutations = filtered_data[multi_hit_mask]
single_hit_mutations = filtered_data[single_hit_mask]

# Step 2: Count multi-hit and single-hit mutations per clone
clone_multi_hit_counts = multi_hit_mutations[alt_count_columns].sum()
clone_single_hit_counts = single_hit_mutations[alt_count_columns].sum()

# Step 3: Combine counts into a single DataFrame
clone_summary = pd.DataFrame({
    'clone_id': [col.replace('_alt_counts', '') for col in alt_count_columns],
    'multi_count': clone_multi_hit_counts.values,
    'single_count': clone_single_hit_counts.values
})

# Step 4: Extract 'regime' (SEX/ASEX) based on the suffix and then remove it
clone_summary['regime'] = clone_summary['clone_id'].apply(
    lambda x: 'SEX' if re.search(r'_SEX$', x) else 'ASEX' if re.search(r'_ASEX$', x) else None
)

# Remove '_SEX' or '_ASEX' suffix from clone_id
clone_summary['clone_id'] = clone_summary['clone_id'].str.replace(r'(_SEX|_ASEX)$', '', regex=True)

# Step 5: Add the 'hit_tot' column (sum of single_count and multi_count)
clone_summary['hit_tot'] = clone_summary['single_count'] + clone_summary['multi_count']

# Step 6: Add the 'single_multi' column (single_count / multi_count), handle division by zero
clone_summary['single_multi'] = clone_summary.apply(
    lambda row: row['single_count'] / row['multi_count'] if row['multi_count'] > 0 else None, axis=1
)

# Optional: Save to a CSV file
clone_summary.to_csv('mut_dat2.csv', index=False)

print("Clone summary with stats saved to 'mut_dat2.csv'.")