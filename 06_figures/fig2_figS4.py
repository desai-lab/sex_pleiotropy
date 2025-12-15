import sys
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from collections import Counter
import seaborn as sns
from matplotlib_venn import venn2
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
from scipy.stats import shapiro, ttest_ind
from statsmodels.stats.multitest import multipletests
from scipy.stats import brunnermunzel

file_path = 'GGE_clones_final.tsv'
np.random.seed(24)

### FIGURE 3A: Stacked bar plot of fixed & non-fixed mutations across sexual & asexual clones

def generate_figure3a(file_path):

    COMBINE_FIXED_NONFIXED = True   # set False to go back to the original 4-section plot

    # Load the data
    #file_path = sys.argv[1]  # Replace with your actual file path
    data = pd.read_csv(file_path, sep='\t')

    # Retain only `alt_counts` columns and `ANN_simpler` for mutation details
    alt_count_columns = [col for col in data.columns if '_alt_counts' in col and not col.endswith(('_3_alt_counts', '_4_alt_counts'))]
    data = data[alt_count_columns + ['ANN_simpler']]

    # Extract the relevant mutation category from `ANN_simpler`
    mutation_categories = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
    data['Parsed_Category'] = data['ANN_simpler'].apply(lambda x: next((cat for cat in mutation_categories if cat in x), 'unknown'))

    # Convert `alt_counts` values to binary (0 if 0, 1 if nonzero)
    data[alt_count_columns] = data[alt_count_columns].apply(lambda x: x.map(lambda v: 1 if v > 0 else 0))

    # Identify base names and suffixes for filtering
    clone_columns = alt_count_columns
    sorted_clones = sorted(set(
        col.replace('_alt_counts', '') for col in clone_columns
    ))

    base_to_suffixes = defaultdict(list)
    for clone in sorted_clones:
        base = '_'.join(clone.split('_')[:-2])  # Base name before _SEX or _ASEX
        suffix = clone.split('_')[-2]  # Extract the suffix (_1, _2, etc.)
        base_to_suffixes[base].append(suffix)

    # Filter valid bases (all bases are valid as clones ending in _3 or _4 are excluded)
    valid_bases = [
        base for base, suffixes in base_to_suffixes.items()
        if '1' in suffixes and '2' in suffixes
    ]

    # Filter valid clone columns based on valid bases
    valid_clones = [
        col for col in alt_count_columns
        if any(base in col for base in valid_bases) and any(suffix in col for suffix in ['_1', '_2'])
    ]

    # Separate sexual and asexual clones
    sexual_clones = [col for col in valid_clones if '_SEX_' in col]
    asexual_clones = [col for col in valid_clones if '_ASEX_' in col]

    # Initialize lists to hold fixed and non-fixed mutations for sexual and asexual clones
    fixed_sexual_mutations = []
    nonfixed_sexual_mutations = []
    fixed_asexual_mutations = []
    nonfixed_asexual_mutations = []

    # Function to ensure consistent counts for fixed mutations
    def get_consistent_counts(base, clone_1, clone_2, data, fixed=True):
        if fixed:
            mutations = (data[clone_1] == 1) & (data[clone_2] == 1)
        else:
            mutations = (data[clone_1] == 1) ^ (data[clone_2] == 1)
        mutation_data = data.loc[mutations, ['Parsed_Category']]
        return mutation_data['Parsed_Category'].value_counts().to_dict()

    # Iterate over valid bases for sexual clones
    sexual_counts_fixed = {}
    sexual_counts_nonfixed = {}
    for base in valid_bases:
        clone_1 = f"{base}_1_SEX_alt_counts"
        clone_2 = f"{base}_2_SEX_alt_counts"
        if clone_1 in data.columns and clone_2 in data.columns:
            # Fixed mutations: consistent across both clones
            sexual_counts_fixed[base] = get_consistent_counts(base, clone_1, clone_2, data, fixed=True)
            # Non-fixed mutations: independent for each clone
            nonfixed_1 = (data[clone_1] == 1) & ~(data[clone_2] == 1)
            nonfixed_2 = (data[clone_2] == 1) & ~(data[clone_1] == 1)
            sexual_counts_nonfixed[clone_1] = data.loc[nonfixed_1, 'Parsed_Category'].value_counts().to_dict()
            sexual_counts_nonfixed[clone_2] = data.loc[nonfixed_2, 'Parsed_Category'].value_counts().to_dict()

    # Iterate over valid bases for asexual clones
    asexual_counts_fixed = {}
    asexual_counts_nonfixed = {}
    for base in valid_bases:
        clone_1 = f"{base}_1_ASEX_alt_counts"
        clone_2 = f"{base}_2_ASEX_alt_counts"
        if clone_1 in data.columns and clone_2 in data.columns:
            # Fixed mutations: consistent across both clones
            asexual_counts_fixed[base] = get_consistent_counts(base, clone_1, clone_2, data, fixed=True)
            # Non-fixed mutations: independent for each clone
            nonfixed_1 = (data[clone_1] == 1) & ~(data[clone_2] == 1)
            nonfixed_2 = (data[clone_2] == 1) & ~(data[clone_1] == 1)
            asexual_counts_nonfixed[clone_1] = data.loc[nonfixed_1, 'Parsed_Category'].value_counts().to_dict()
            asexual_counts_nonfixed[clone_2] = data.loc[nonfixed_2, 'Parsed_Category'].value_counts().to_dict()

    def _sum_dicts(*ds):
        c = Counter()
        for d in ds:
            c += Counter(d or {})
        return dict(c)

    combined_asexual = {}
    combined_sexual  = {}
    for base in valid_bases:
        # asexual
        a_c1 = f"{base}_1_ASEX_alt_counts"
        a_c2 = f"{base}_2_ASEX_alt_counts"
        combined_asexual[base] = _sum_dicts(
            asexual_counts_fixed.get(base, {}),
            asexual_counts_nonfixed.get(a_c1, {}),
            asexual_counts_nonfixed.get(a_c2, {})
        )
        # sexual
        s_c1 = f"{base}_1_SEX_alt_counts"
        s_c2 = f"{base}_2_SEX_alt_counts"
        combined_sexual[base] = _sum_dicts(
            sexual_counts_fixed.get(base, {}),
            sexual_counts_nonfixed.get(s_c1, {}),
            sexual_counts_nonfixed.get(s_c2, {})
        )

    # Per-clone combined counts (counts of all mutations seen in that clone)
    combined_asexual_clones = {}
    combined_sexual_clones  = {}

    for base in valid_bases:
        a1 = f"{base}_1_ASEX_alt_counts"; a2 = f"{base}_2_ASEX_alt_counts"
        s1 = f"{base}_1_SEX_alt_counts";  s2 = f"{base}_2_SEX_alt_counts"

        if a1 in data.columns:
            combined_asexual_clones[a1] = data.loc[data[a1] == 1, 'Parsed_Category'].value_counts().to_dict()
        if a2 in data.columns:
            combined_asexual_clones[a2] = data.loc[data[a2] == 1, 'Parsed_Category'].value_counts().to_dict()

        if s1 in data.columns:
            combined_sexual_clones[s1]  = data.loc[data[s1] == 1, 'Parsed_Category'].value_counts().to_dict()
        if s2 in data.columns:
            combined_sexual_clones[s2]  = data.loc[data[s2] == 1, 'Parsed_Category'].value_counts().to_dict()


    # Combine all categories into one plot with 4 sections on the x-axis
    fig, ax = plt.subplots(figsize=(10, 10))

    # Consistent color scheme for mutation categories and specific order
    colors = {
        'indel': '#B32316',
        'nonsense': '#3D0066',
        'missense': '#A369CB',
        'synonymous': '#EEC7FC',
        'noncoding': '#999999',
    }

    # Sort clones by descending total mutation counts within each group
    def sort_clones_by_counts(counts):
        return sorted(counts.keys(), key=lambda x: sum(counts[x].values()), reverse=True)

    # Define sections for the x-axis in the specified order

    # sections = [
    #     ("Fixed Asexual", asexual_counts_fixed),
    #     ("Fixed Sexual", sexual_counts_fixed),
    #     ("Non-Fixed Asexual", asexual_counts_nonfixed),
    #     ("Non-Fixed Sexual", sexual_counts_nonfixed),
    # ]

    if COMBINE_FIXED_NONFIXED:
        # one stacked bar per clone (fixed+nonfixed together)
        sections = [
            ("Combined Asexual Clones", combined_asexual_clones),
            ("Combined Sexual Clones",  combined_sexual_clones),
        ]
    else:
        sections = [
            ("Fixed Asexual", asexual_counts_fixed),
            ("Fixed Sexual", sexual_counts_fixed),
            ("Non-Fixed Asexual", asexual_counts_nonfixed),
            ("Non-Fixed Sexual", sexual_counts_nonfixed),
        ]



    x_labels = []
    x_positions = []
    current_x = 0
    section_starts = []
    section_ends = []

    for section_name, counts in sections:
        is_fixed_like = ("Fixed" in section_name)

        if is_fixed_like:
            # For fixed sections, only plot one bar per base (merge `_1` and `_2` clones)
            bases = sort_clones_by_counts(counts)
            section_positions = range(current_x, current_x + len(bases))
            x_positions.extend(section_positions)
            x_labels.extend([section_name] * len(bases))
            section_starts.append(current_x)
            section_ends.append(current_x + len(bases) - 1)

            bottom = [0] * len(bases)
            for cat in ['noncoding', 'synonymous', 'missense', 'nonsense', 'indel']:
                heights = [counts[base].get(cat, 0) for base in bases]
                ax.bar(section_positions, heights, bottom=bottom, label=cat if current_x == 0 else None, color=colors[cat])
                bottom = [sum(x) for x in zip(bottom, heights)]
        else:
            # For non-fixed sections, plot each clone individually
            clones = sort_clones_by_counts(counts)
            section_positions = range(current_x, current_x + len(clones))
            x_positions.extend(section_positions)
            x_labels.extend([section_name] * len(clones))
            section_starts.append(current_x)
            section_ends.append(current_x + len(clones) - 1)

            bottom = [0] * len(clones)
            for cat in ['noncoding', 'synonymous', 'missense', 'nonsense', 'indel']:
                heights = [counts[clone].get(cat, 0) for clone in clones]
                ax.bar(section_positions, heights, bottom=bottom, label=cat if current_x == 0 else None, color=colors[cat])
                bottom = [sum(x) for x in zip(bottom, heights)]

        current_x += (len(bases) + 2) if is_fixed_like else (len(clones) + 2)
        #current_x += len(bases) + 2 if "Fixed" in section_name else len(clones) + 2  # Add spacing between sections

    # Add section labels below the x-axis
    section_midpoints = [(start + end) // 2 for start, end in zip(section_starts, section_ends)]
    ax.set_xticks(section_midpoints)

    if COMBINE_FIXED_NONFIXED:
        ax.set_xticklabels(["Asexual Clones", "Sexual Clones"], fontsize=22, rotation=0)
    else:
        ax.text((section_starts[0] + section_ends[1]) // 2, -0.1, "Fixed Mutations",
                ha='center', fontsize=16, transform=ax.get_xaxis_transform())
        ax.text((section_starts[2] + section_ends[3]) // 2, -0.1, "Non-Fixed Mutations",
                ha='center', fontsize=16, transform=ax.get_xaxis_transform())
        ax.set_xticklabels(["Asexual Populations", "Sexual Populations", "Asexual Clones", "Sexual Clones"],
                        fontsize=14, rotation=0)


    # # Add lines and labels below x-axis sections
    # ax.text((section_starts[0] + section_ends[1]) // 2, -0.1, "Fixed Mutations", ha='center', fontsize=16, transform=ax.get_xaxis_transform())
    # ax.text((section_starts[2] + section_ends[3]) // 2, -0.1, "Non-Fixed Mutations", ha='center', fontsize=16, transform=ax.get_xaxis_transform())
    # ax.set_xticklabels(["Asexual Populations", "Sexual Populations", "Asexual Clones", "Sexual Clones"], fontsize=14, rotation=0)
    ax.tick_params(axis='x', length=0, pad=15)

    # Add legend in the correct order
    handles, labels = ax.get_legend_handles_labels()
    order = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
    ordered_handles = [handles[labels.index(cat)] for cat in order]
    ordered_labels = [cat.capitalize() for cat in order]
    ax.legend(ordered_handles, ordered_labels, loc='center', bbox_to_anchor=(0.375, 0.875), fontsize=16, frameon=False)

    # y axis
    ax.tick_params(axis='y', labelsize=12)
    ax.set_ylim(0, 145)
    ax.set_ylabel("Mutation Count", fontsize=26)

    # Remove top, bottom, and right borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.tight_layout()
    #plt.savefig("figure3a.pdf", dpi=300)
    #plt.savefig("figure3a.jpg", dpi=300)
    #plt.close(fig)

    ### ================== STATS ================== ###

    def compute_summary(label, counts_dict):
        total_counts = [sum(cat_counts.values()) for cat_counts in counts_dict.values()]
        mean_val = np.mean(total_counts)
        median_val = np.median(total_counts)
        return {
            'label': label,
            'counts': total_counts,
            'mean': mean_val,
            'median': median_val
        }

    def shapiro_test(label, counts):
        stat, p = shapiro(counts)
        return {
            'label': label,
            'statistic': stat,
            'pvalue': p
        }

    # Compute basic stats
    summary_stats = [
        compute_summary("Fixed Asexual", asexual_counts_fixed),
        compute_summary("Fixed Sexual", sexual_counts_fixed),
        compute_summary("Non-Fixed Asexual", asexual_counts_nonfixed),
        compute_summary("Non-Fixed Sexual", sexual_counts_nonfixed),
    ]

    # Compute Shapiro-Wilk test results
    shapiro_stats_raw = [
        shapiro_test(s['label'], s['counts']) for s in summary_stats
    ]

    # Extract p-values for BH correction
    raw_pvals = [s['pvalue'] for s in shapiro_stats_raw]

    # Apply Benjamini-Hochberg correction
    reject, pvals_corrected, _, _ = multipletests(raw_pvals, method='fdr_bh')

    # Combine into final list with corrected p-values
    shapiro_stats = []
    for original, corrected_p, is_rejected in zip(shapiro_stats_raw, pvals_corrected, reject):
        result = {
            'label': original['label'],
            'statistic': original['statistic'],
            'pvalue_raw': original['pvalue'],
            'pvalue_corrected': corrected_p,
            'rejected': is_rejected
        }
        shapiro_stats.append(result)

    # Run Brunner-Munzel tests
    def run_brunner_munzel(label, group1_counts, group2_counts):
        stat, p = brunnermunzel(group1_counts, group2_counts, alternative='two-sided')
        return {
            'label': label,
            'statistic': stat,
            'pvalue': p
        }

    # Extract mutation counts
    fixed_asexual_counts = [sum(d.values()) for d in asexual_counts_fixed.values()]
    fixed_sexual_counts = [sum(d.values()) for d in sexual_counts_fixed.values()]
    nonfixed_asexual_counts = [sum(d.values()) for d in asexual_counts_nonfixed.values()]
    nonfixed_sexual_counts = [sum(d.values()) for d in sexual_counts_nonfixed.values()]

    # Perform tests
    brunner_results_raw = [
        run_brunner_munzel("Fixed: Asexual vs Sexual", fixed_asexual_counts, fixed_sexual_counts),
        run_brunner_munzel("Non-Fixed: Asexual vs Sexual", nonfixed_asexual_counts, nonfixed_sexual_counts),
    ]

    # Apply BH correction to Brunner-Munzel p-values
    brunner_raw_pvals = [r['pvalue'] for r in brunner_results_raw]
    reject_brunner, corrected_brunner_pvals, _, _ = multipletests(brunner_raw_pvals, method='fdr_bh')

    # Combine results
    brunner_results = []
    for res, p_corr, reject in zip(brunner_results_raw, corrected_brunner_pvals, reject_brunner):
        brunner_results.append({
            'label': res['label'],
            'statistic': res['statistic'],
            'pvalue_raw': res['pvalue'],
            'pvalue_corrected': p_corr,
            'rejected': reject
        })

    return fig, summary_stats, shapiro_stats, brunner_results

### FIGURE 3A inset: Venn diagram of mutational targets across sexual and asexual clones

def generate_figure3c(file_path):
    data = pd.read_csv(file_path, sep='\t')

    # Extract the relevant mutation category and ORF from `ANN_simpler`
    mutation_categories = ['indel', 'nonsense', 'missense', 'synonymous']
    data['Parsed_Category'] = data['ANN_simpler'].apply(lambda x: next((cat for cat in mutation_categories if cat in x), 'unknown'))

    functional_categories = ['indel', 'nonsense', 'missense']

    # Filter data to keep only functional mutation categories
    data = data[data['Parsed_Category'].isin(functional_categories)]

    # Extract ORFs from `ANN_simpler` (4th field when split by '|')
    data['ORF'] = data['ANN_simpler'].str.split('|').str[3]

    # Separate sexual and asexual clones
    sexual_clones = [col for col in data.columns if '_SEX_' in col and '_alt_counts' in col]
    asexual_clones = [col for col in data.columns if '_ASEX_' in col and '_alt_counts' in col]

    # Identify mutational targets for sexual clones
    sexual_targets = set(
        data.loc[(data[sexual_clones].sum(axis=1) > 0), 'ORF']
    )

    # Identify mutational targets for asexual clones
    asexual_targets = set(
        data.loc[(data[asexual_clones].sum(axis=1) > 0), 'ORF']
    )

    ## Finding multi-hit and single-hit genes

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

    # Separate sexual and asexual populations based on `SEX` and `ASEX` in population_columns
    sexual_populations = list(set(
        '_'.join(pop.rsplit('_', 2)[:-2]) for pop in population_columns if '_SEX' in pop
    ))
    asexual_populations = list(set(
        '_'.join(pop.rsplit('_', 2)[:-2]) for pop in population_columns if '_ASEX' in pop
    ))

    # Separate multi-hit genes by population type
    sexual_multi_hit_genes = set(gene for gene in multi_hit_genes if any(pop in sexual_populations for pop in gene_population_map[gene]))
    asexual_multi_hit_genes = set(gene for gene in multi_hit_genes if any(pop in asexual_populations for pop in gene_population_map[gene]))

    # Separate single-hit genes by population type
    sexual_single_hit_genes = set(gene for gene in single_hit_genes if any(pop in sexual_populations for pop in gene_population_map[gene]))
    asexual_single_hit_genes = set(gene for gene in single_hit_genes if any(pop in asexual_populations for pop in gene_population_map[gene]))

    ## Generate Venn diagrams

    # Define the data for the three Venn diagrams
    venn_data = {
        "All Genes": {
            "asexual": asexual_targets,
            "sexual": sexual_targets,
        },
        "Multi-Hit Genes": {
            "asexual": asexual_multi_hit_genes,
            "sexual": sexual_multi_hit_genes,
        },
        "Single-Hit Genes": {
            "asexual": asexual_single_hit_genes,
            "sexual": sexual_single_hit_genes,
        },
    }

    # Calculate areas based on unique and intersection counts
    data_areas = [
        (len(data['asexual'] - data['sexual']),  # Sexual-only
        len(data['sexual'] - data['asexual']),  # Asexual-only
        len(data['sexual'] & data['asexual']))  # Intersection
        for data in venn_data.values()
    ]
    max_area = max(sum(areas) for areas in data_areas)  # Find the maximum total area

    # Define a list of vertical limit multipliers, one for each Venn diagram
    vertical_multipliers = [0.5, 0.25, 0.45]

    # Function to tightly center-align the Venn diagram and remove vertical whitespace
    def center_align_venn(ax, asexual_only, sexual_only, intersection, venn, i):
        total_area = sexual_only + asexual_only + intersection
        scale_factor = np.sqrt(max_area / total_area) if total_area > 0 else 1
        vertical_multiplier = vertical_multipliers[i]

        # If there is no intersection, manually center-align the circles and labels
        if intersection == 0:
            circle1 = venn.get_patch_by_id('10')
            circle2 = venn.get_patch_by_id('01')
            if circle1:
                circle1.set_center((-0.3, 0))  # Adjust position of Sexual-only circle
            if circle2:
                circle2.set_center((0.6, 0))   # Adjust position of Asexual-only circle
            venn.get_label_by_id('10').set_position((-0.3, 0))  # Adjust Sexual-only label
            venn.get_label_by_id('01').set_position((0.6, 0))   # Adjust Asexual-only label

        # Tightly set axis limits
        ax.set_xlim(-0.9 * scale_factor, 0.9 * scale_factor)  # Reduced padding
        ax.set_ylim(-vertical_multiplier * scale_factor, vertical_multiplier * scale_factor)  # Adjusted vertical limits

    # Create a single figure with GridSpec for precise layout
    fig = plt.figure(figsize=(6, 10))  # Reduced figure height
    gs = GridSpec(3, 2, width_ratios=[0.1, 3.5], height_ratios=[1, 0.4, 0.7], figure=fig, hspace=0)

    # Colors for the Venn diagrams
    venn_colors = ('#4682B4', '#EE7600')

    # Loop through the data and plot each Venn diagram
    for i, ((title, data), (asexual_only, sexual_only, intersection)) in enumerate(zip(venn_data.items(), data_areas)):
        # Title on the left
        ax_title = fig.add_subplot(gs[i, 0])
        ax_title.axis("off")
        ax_title.text(1.75, 0.5, title, fontsize=20, rotation=90, va='center', ha='center')

        # Venn diagram on the right
        ax_venn = fig.add_subplot(gs[i, 1])
        venn = venn2(
            subsets=(asexual_only, sexual_only, intersection),
            set_labels=('Asexual', 'Sexual') if i == 0 else (None, None),
            set_colors=venn_colors,
            ax=ax_venn,
            alpha=0.8
        )

        # Increase the font size of the number labels
        for label_id in ['10', '01', '11']:
            label = venn.get_label_by_id(label_id)
            if label:  # Check if the label exists
                label.set_fontsize(20)  # Set font size to 16 (adjust as needed)

        # Customize the intersection color
        if venn.get_patch_by_id('11'):  # Check if intersection exists
            venn.get_patch_by_id('11').set_color('#C4CBD2')
            venn.get_patch_by_id('11').set_alpha(0.8)

        # Adjust font size for "Sexual" and "Asexual" labels (first Venn diagram only)
        if i == 0:
            if venn.set_labels[0]:  # Adjust font size for "Sexual" label
                venn.set_labels[0].set_fontsize(16)
            if venn.set_labels[1]:  # Adjust font size for "Asexual" label
                venn.set_labels[1].set_fontsize(16)

        # Center-align the Venn diagram and remove vertical whitespace
        center_align_venn(ax_venn, asexual_only, sexual_only, intersection, venn, i)

    plt.subplots_adjust(top=1, bottom=0, left=0, right=1, hspace=0, wspace=0)

    # Save the figure
    #plt.savefig("figure3c.pdf", dpi=300)
    #plt.savefig("figure3c.jpg", dpi=300)
    #plt.close(fig)
    return fig

### FIGURE 3B: Number of putatively functional mutations in sexual & asexual clones

def generate_figure3b(file_path):
    # Load the data
    #file_path = sys.argv[1]  # Replace with your actual file path
    data = pd.read_csv(file_path, sep='\t')

    # Retain only `alt_counts` columns and `ANN_simpler` for mutation details
    alt_count_columns = [col for col in data.columns if '_alt_counts' in col and not col.endswith(('_3_alt_counts', '_4_alt_counts'))]
    data = data[alt_count_columns + ['ANN_simpler']]

    # Extract the relevant mutation category from `ANN_simpler`
    mutation_categories = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
    data['Parsed_Category'] = data['ANN_simpler'].apply(lambda x: next((cat for cat in mutation_categories if cat in x), 'unknown'))

    # Convert `alt_counts` values to binary (0 if 0, 1 if nonzero)
    data[alt_count_columns] = data[alt_count_columns].apply(lambda x: x.map(lambda v: 1 if v > 0 else 0))

    valid_clones = [
        col for col in alt_count_columns
    ]

    # Uncomment this if you want to filter out unpaired clones or _3 / _4 clones
    # base_to_suffixes = defaultdict(list)
    # for clone in sorted_clones:
    #     base = '_'.join(clone.split('_')[:-2])  # Base name before _SEX or _ASEX
    #     suffix = clone.split('_')[-2]  # Extract the suffix (_1, _2, etc.)
    #     base_to_suffixes[base].append(suffix)

    # # Filter valid bases (populations with both _1 and _2 clones)
    # valid_bases = [
    #     base for base, suffixes in base_to_suffixes.items()
    #     if '1' in suffixes and '2' in suffixes
    # ]

    # # Filter valid clone columns based on valid bases
    # valid_clones = [
    #     col for col in alt_count_columns
    #     if any(base in col for base in valid_bases) and any(suffix in col for suffix in ['_1', '_2'])
    # ]

    # Separate sexual and asexual clones
    sexual_clones = [col for col in valid_clones if '_SEX_' in col]
    asexual_clones = [col for col in valid_clones if '_ASEX_' in col]

    # Functional mutation categories
    functional_categories = ['indel', 'nonsense', 'missense']

    # Filter data to keep only functional mutation categories
    data = data[data['Parsed_Category'].isin(functional_categories)]

    # Calculate total functional mutations per clone
    clone_functional_mutations = {clone: data[clone].sum() for clone in valid_clones}

    # Prepare data for plotting
    plot_data = []
    for clone in sexual_clones:
        plot_data.append({'Clone_Type': 'Sexual', 'Functional_Mutations': clone_functional_mutations[clone]})
    for clone in asexual_clones:
        plot_data.append({'Clone_Type': 'Asexual', 'Functional_Mutations': clone_functional_mutations[clone]})

    plot_df = pd.DataFrame(plot_data)

    # Plot the dot plot
    fig = plt.figure(figsize=(5, 10.5))
    ax = sns.stripplot(
        hue='Clone_Type',
        data=plot_df,
        x='Clone_Type',
        y='Functional_Mutations',
        jitter=True,
        palette={'Sexual': '#EE7600', 'Asexual': '#4682B4'},
        legend=False,
        size=24,
        edgecolor='black', linewidth=2.4,
        order=['Asexual', 'Sexual']
    )

    # Add median lines
    for clone_type in ['Asexual', 'Sexual']:
        median = plot_df[plot_df['Clone_Type'] == clone_type]['Functional_Mutations'].median()
        x_position = 0 if clone_type == 'Asexual' else 1
        plt.hlines(median, x_position - 0.3, x_position + 0.3, colors='black', linewidth=12, zorder=3)

    # Customize the plot
    plt.ylabel("Total Functional Mutations per Clone", fontsize=30, labelpad=10)
    plt.xlabel("")

    # Customize tick labels and ticks
    ax.tick_params(axis='x', labelsize=25, width=2, length=6)  # Customize x-axis ticks
    ax.tick_params(axis='y', labelsize=20, width=2, length=6)  # Customize y-axis ticks

    for spine in ax.spines.values():
        spine.set_linewidth(2)  # Set thickness of axis lines

    # Remove top, and right borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    #plt.savefig("figure3b.pdf", dpi=300)
    #plt.savefig("figure3b.jpg", dpi=300)
    #plt.close(fig)

    ### ================== STATS ================== ###

    # === Shapiro-Wilk Normality Tests ===
    asexual_values = plot_df[plot_df['Clone_Type'] == 'Asexual']['Functional_Mutations'].tolist()
    sexual_values = plot_df[plot_df['Clone_Type'] == 'Sexual']['Functional_Mutations'].tolist()

    def run_shapiro(label, values):
        stat, p = shapiro(values)
        return {
            'label': label,
            'statistic': stat,
            'pvalue': p,
            'n': len(values),
            'mean': np.mean(values),
            'median': np.median(values)
        }

    # Run tests
    shapiro_results = [
        run_shapiro("Asexual", asexual_values),
        run_shapiro("Sexual", sexual_values),
    ]

    # Extract raw p-values
    raw_pvals = [res['pvalue'] for res in shapiro_results]

    # Apply BH correction
    reject, corrected_pvals, _, _ = multipletests(raw_pvals, method='fdr_bh')

    # Append corrected p-values to the results
    for i, res in enumerate(shapiro_results):
        res['pvalue_corrected'] = corrected_pvals[i]
        res['rejected'] = reject[i]

    # === Welch's t-test ===
    t_stat, t_pval = ttest_ind(asexual_values, sexual_values, equal_var=False)
    welch_result = {
        'label': 'Asexual vs Sexual (Welch)',
        't_statistic': t_stat,
        'pvalue': t_pval,
        'n_asexual': len(asexual_values),
        'n_sexual': len(sexual_values),
        'mean_asexual': np.mean(asexual_values),
        'mean_sexual': np.mean(sexual_values)
    }

    return fig, shapiro_results, welch_result

### FIGURE 3C: Fraction of functional mutations in multi-hit genes

def generate_figure3d(file_path):
    data = pd.read_csv(file_path, sep='\t')

    # Retain only `alt_counts` columns and `ANN_simpler` for mutation details
    alt_count_columns = [col for col in data.columns if '_alt_counts' in col]
    filtered_data = data[alt_count_columns + ['ANN_simpler']].copy()

    # Extract mutation categories
    mutation_categories = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
    filtered_data['Parsed_Category'] = filtered_data['ANN_simpler'].apply(
        lambda x: next((cat for cat in mutation_categories if cat in str(x)), 'unknown')
    )

    # Exclude synonymous and noncoding mutations
    excluded_categories = ['synonymous', 'noncoding']
    filtered_data = filtered_data[~filtered_data['Parsed_Category'].isin(excluded_categories)]

    # Convert `alt_counts` columns to binary (1 if mutation is present, 0 if absent)
    filtered_data[alt_count_columns] = filtered_data[alt_count_columns].map(lambda x: 1 if x > 0 else 0)

    # Identify populations by removing clone suffix (_1, _2, etc.)
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
        gene = str(row['ANN_simpler']).split('|')[3] if len(str(row['ANN_simpler']).split('|')) > 3 else None
        if gene:
            for population, clones in pop_map.items():
                if row[clones].sum() > 0:  # If at least one clone in the population has the mutation
                    if gene not in gene_population_map:
                        gene_population_map[gene] = set()
                    gene_population_map[gene].add(population)

    # Identify **multi-hit** genes (genes mutated in multiple populations)
    multi_hit_genes = {gene for gene, pops in gene_population_map.items() if len(pops) > 1}

    # Compute the ratio of mutations in multi-hit genes for each clone
    clone_mutation_ratios = []
    for clone in alt_count_columns:
        total_mutations = filtered_data[clone].sum()
        multi_hit_mutations = filtered_data.loc[
            filtered_data['ANN_simpler'].apply(lambda x: str(x).split('|')[3] in multi_hit_genes), clone
        ].sum()
        ratio = multi_hit_mutations / total_mutations if total_mutations > 0 else 0
        clone_type = 'Sexual' if '_SEX_' in clone else 'Asexual'
        clone_mutation_ratios.append({'Clone_Type': clone_type, 'Mutation_Ratio': ratio})

    # Convert to DataFrame
    plot_df = pd.DataFrame(clone_mutation_ratios)

    # Create the dot plot
    fig = plt.figure(figsize=(5, 10.5))
    ax = sns.stripplot(
        hue='Clone_Type',
        data=plot_df,
        x='Clone_Type',
        y='Mutation_Ratio',
        jitter=True,
        palette={'Sexual': '#EE7600', 'Asexual': '#4682B4'},
        legend=False,
        size=24,
        edgecolor='black', linewidth=2.4,
        order=['Asexual', 'Sexual']
    )

    # Add median lines
    for clone_type in ['Asexual', 'Sexual']:
        median = plot_df[plot_df['Clone_Type'] == clone_type]['Mutation_Ratio'].median()
        x_position = 0 if clone_type == 'Asexual' else 1
        plt.hlines(median, x_position - 0.3, x_position + 0.3, colors='black', linewidth=12, zorder=3)

    # Customize the plot
    plt.ylabel("Fraction of Functional Mutations in Multi-Hit Genes", fontsize=26, labelpad=10)
    plt.xlabel("")

    # Customize tick labels and ticks
    ax.tick_params(axis='x', labelsize=25, width=2, length=6)
    ax.tick_params(axis='y', labelsize=20, width=2, length=6)

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    # Remove top and right borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    #plt.savefig("figure3d.pdf", dpi=300)
    #plt.savefig("figure3d.jpg", dpi=300)
    #plt.close(fig)

    ### ================== STATS ================== ###

    # === Shapiro-Wilk Normality Tests ===
    asexual_values = plot_df[plot_df['Clone_Type'] == 'Asexual']['Mutation_Ratio'].tolist()
    sexual_values = plot_df[plot_df['Clone_Type'] == 'Sexual']['Mutation_Ratio'].tolist()

    def run_shapiro(label, values):
        stat, p = shapiro(values)
        return {
            'label': label,
            'statistic': stat,
            'pvalue': p,
            'n': len(values),
            'mean': np.mean(values),
            'median': np.median(values)
        }

    # Run tests
    shapiro_results = [
        run_shapiro("Asexual", asexual_values),
        run_shapiro("Sexual", sexual_values),
    ]

    # Apply BH correction
    raw_pvals = [res['pvalue'] for res in shapiro_results]
    reject, corrected_pvals, _, _ = multipletests(raw_pvals, method='fdr_bh')

    # Add corrected p-values to the results
    for i, res in enumerate(shapiro_results):
        res['pvalue_corrected'] = corrected_pvals[i]
        res['rejected'] = reject[i]

    # === Welch's t-test ===
    t_stat, t_pval = ttest_ind(asexual_values, sexual_values, equal_var=False)
    welch_result = {
        'label': 'Asexual vs Sexual (Welch)',
        't_statistic': t_stat,
        'pvalue': t_pval,
        'n_asexual': len(asexual_values),
        'n_sexual': len(sexual_values),
        'mean_asexual': np.mean(asexual_values),
        'mean_sexual': np.mean(sexual_values)
    }

    return fig, shapiro_results, welch_result

# Generate individual figures
fig3a, summary_stats, shapiro_stats, brunner_results = generate_figure3a(file_path)
fig3b, shapiro_results_3b, welch_result_3b = generate_figure3b(file_path)
fig3c = generate_figure3c(file_path)
fig3d, shapiro_results_3d, welch_result_3d = generate_figure3d(file_path)

# Create master figure (2 cols: top A spans both; bottom has B | D)
master_fig = plt.figure(figsize=(10, 18))
gs = GridSpec(2, 2, height_ratios=[3, 2.5], figure=master_fig, hspace=0.075, wspace=0.02)

# --- Top: A spans both columns ---
ax_master_3a = master_fig.add_subplot(gs[0, :])
fig3a.canvas.draw()
ax_master_3a.imshow(fig3a.canvas.buffer_rgba())
ax_master_3a.axis('off')
#ax_master_3a.set_title("A", loc='left', fontsize=30, fontweight='bold')
ax_master_3a.text(0.01, 1, "A", transform=ax_master_3a.transAxes,
                  fontsize=30, fontweight='bold', va='top', ha='right')

# ----- Inset C inside A (top-right) -----
# (x0, y0, width, height) in axes coords; adjust to taste
INSET = [0.575, 0.35, 0.34*1.2, 0.52*1.2]  # top-right placement
inset_ax = ax_master_3a.inset_axes(INSET, zorder=5)
inset_ax.set_facecolor('white')
fig3c.canvas.draw()
inset_ax.imshow(fig3c.canvas.buffer_rgba())
inset_ax.axis('off')
# panel label for C
#inset_ax.text(-0.08, 1.02, "C", transform=inset_ax.transAxes,fontsize=30, fontweight='bold', va='bottom', ha='right')

# --- Bottom row: only B (left) and D (right) ---
# B
ax_master_3b = master_fig.add_subplot(gs[1, 0])
fig3b.canvas.draw()
ax_master_3b.imshow(fig3b.canvas.buffer_rgba())
ax_master_3b.axis('off')
ax_master_3b.text(-0.07, 1, "B", transform=ax_master_3b.transAxes,
                  fontsize=30, fontweight='bold', va='top', ha='right')

# D
ax_master_3d = master_fig.add_subplot(gs[1, 1])
fig3d.canvas.draw()
ax_master_3d.imshow(fig3d.canvas.buffer_rgba())
ax_master_3d.axis('off')
ax_master_3d.text(-0.07, 1, "C", transform=ax_master_3d.transAxes,
                  fontsize=30, fontweight='bold', va='top', ha='right')

ax_master_3b.set_position(ax_master_3b.get_position().translated(+0.02, 0))  # move right
# ax_master_3d.set_position(ax_master_3d.get_position().translated(-0.01, 0))  # move left

# Save or display
plt.savefig("figure2.pdf", dpi=300, bbox_inches='tight')
plt.savefig("figure2.jpg", dpi=300, bbox_inches='tight')
#plt.show()

### FIGURE S4: Stacked bar plot of fixed & non-fixed mutations across sexual & asexual clones

def generate_figureS4(file_path):
    # Load the data
    #file_path = sys.argv[1]  # Replace with your actual file path
    data = pd.read_csv(file_path, sep='\t')

    # Retain only `alt_counts` columns and `ANN_simpler` for mutation details
    alt_count_columns = [col for col in data.columns if '_alt_counts' in col and not col.endswith(('_3_alt_counts', '_4_alt_counts'))]
    data = data[alt_count_columns + ['ANN_simpler']]

    # Extract the relevant mutation category from `ANN_simpler`
    mutation_categories = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
    data['Parsed_Category'] = data['ANN_simpler'].apply(lambda x: next((cat for cat in mutation_categories if cat in x), 'unknown'))

    # Convert `alt_counts` values to binary (0 if 0, 1 if nonzero)
    data[alt_count_columns] = data[alt_count_columns].apply(lambda x: x.map(lambda v: 1 if v > 0 else 0))

    # Identify base names and suffixes for filtering
    clone_columns = alt_count_columns
    sorted_clones = sorted(set(
        col.replace('_alt_counts', '') for col in clone_columns
    ))

    base_to_suffixes = defaultdict(list)
    for clone in sorted_clones:
        base = '_'.join(clone.split('_')[:-2])  # Base name before _SEX or _ASEX
        suffix = clone.split('_')[-2]  # Extract the suffix (_1, _2, etc.)
        base_to_suffixes[base].append(suffix)

    # Filter valid bases (all bases are valid as clones ending in _3 or _4 are excluded)
    valid_bases = [
        base for base, suffixes in base_to_suffixes.items()
        if '1' in suffixes and '2' in suffixes
    ]

    # Filter valid clone columns based on valid bases
    valid_clones = [
        col for col in alt_count_columns
        if any(base in col for base in valid_bases) and any(suffix in col for suffix in ['_1', '_2'])
    ]

    # Separate sexual and asexual clones
    sexual_clones = [col for col in valid_clones if '_SEX_' in col]
    asexual_clones = [col for col in valid_clones if '_ASEX_' in col]

    # Initialize lists to hold fixed and non-fixed mutations for sexual and asexual clones
    fixed_sexual_mutations = []
    nonfixed_sexual_mutations = []
    fixed_asexual_mutations = []
    nonfixed_asexual_mutations = []

    # Function to ensure consistent counts for fixed mutations
    def get_consistent_counts(base, clone_1, clone_2, data, fixed=True):
        if fixed:
            mutations = (data[clone_1] == 1) & (data[clone_2] == 1)
        else:
            mutations = (data[clone_1] == 1) ^ (data[clone_2] == 1)
        mutation_data = data.loc[mutations, ['Parsed_Category']]
        return mutation_data['Parsed_Category'].value_counts().to_dict()

    # Iterate over valid bases for sexual clones
    sexual_counts_fixed = {}
    sexual_counts_nonfixed = {}
    for base in valid_bases:
        clone_1 = f"{base}_1_SEX_alt_counts"
        clone_2 = f"{base}_2_SEX_alt_counts"
        if clone_1 in data.columns and clone_2 in data.columns:
            # Fixed mutations: consistent across both clones
            sexual_counts_fixed[base] = get_consistent_counts(base, clone_1, clone_2, data, fixed=True)
            # Non-fixed mutations: independent for each clone
            nonfixed_1 = (data[clone_1] == 1) & ~(data[clone_2] == 1)
            nonfixed_2 = (data[clone_2] == 1) & ~(data[clone_1] == 1)
            sexual_counts_nonfixed[clone_1] = data.loc[nonfixed_1, 'Parsed_Category'].value_counts().to_dict()
            sexual_counts_nonfixed[clone_2] = data.loc[nonfixed_2, 'Parsed_Category'].value_counts().to_dict()

    # Iterate over valid bases for asexual clones
    asexual_counts_fixed = {}
    asexual_counts_nonfixed = {}
    for base in valid_bases:
        clone_1 = f"{base}_1_ASEX_alt_counts"
        clone_2 = f"{base}_2_ASEX_alt_counts"
        if clone_1 in data.columns and clone_2 in data.columns:
            # Fixed mutations: consistent across both clones
            asexual_counts_fixed[base] = get_consistent_counts(base, clone_1, clone_2, data, fixed=True)
            # Non-fixed mutations: independent for each clone
            nonfixed_1 = (data[clone_1] == 1) & ~(data[clone_2] == 1)
            nonfixed_2 = (data[clone_2] == 1) & ~(data[clone_1] == 1)
            asexual_counts_nonfixed[clone_1] = data.loc[nonfixed_1, 'Parsed_Category'].value_counts().to_dict()
            asexual_counts_nonfixed[clone_2] = data.loc[nonfixed_2, 'Parsed_Category'].value_counts().to_dict()

    # Combine all categories into one plot with 4 sections on the x-axis
    fig, ax = plt.subplots(figsize=(20, 10))

    # Consistent color scheme for mutation categories and specific order
    colors = {
        'indel': '#B32316',
        'nonsense': '#3D0066',
        'missense': '#A369CB',
        'synonymous': '#EEC7FC',
        'noncoding': '#999999',
    }

    # Sort clones by descending total mutation counts within each group
    def sort_clones_by_counts(counts):
        return sorted(counts.keys(), key=lambda x: sum(counts[x].values()), reverse=True)

    # Define sections for the x-axis in the specified order
    sections = [
        ("Fixed Asexual", asexual_counts_fixed),
        ("Fixed Sexual", sexual_counts_fixed),
        ("Non-Fixed Asexual", asexual_counts_nonfixed),
        ("Non-Fixed Sexual", sexual_counts_nonfixed),
    ]

    x_labels = []
    x_positions = []
    current_x = 0
    section_starts = []
    section_ends = []

    for section_name, counts in sections:
        if "Fixed" in section_name:
            # For fixed sections, only plot one bar per base (merge `_1` and `_2` clones)
            bases = sort_clones_by_counts(counts)
            section_positions = range(current_x, current_x + len(bases))
            x_positions.extend(section_positions)
            x_labels.extend([section_name] * len(bases))
            section_starts.append(current_x)
            section_ends.append(current_x + len(bases) - 1)

            bottom = [0] * len(bases)
            for cat in ['noncoding', 'synonymous', 'missense', 'nonsense', 'indel']:
                heights = [counts[base].get(cat, 0) for base in bases]
                ax.bar(section_positions, heights, bottom=bottom, label=cat if current_x == 0 else None, color=colors[cat])
                bottom = [sum(x) for x in zip(bottom, heights)]
        else:
            # For non-fixed sections, plot each clone individually
            clones = sort_clones_by_counts(counts)
            section_positions = range(current_x, current_x + len(clones))
            x_positions.extend(section_positions)
            x_labels.extend([section_name] * len(clones))
            section_starts.append(current_x)
            section_ends.append(current_x + len(clones) - 1)

            bottom = [0] * len(clones)
            for cat in ['noncoding', 'synonymous', 'missense', 'nonsense', 'indel']:
                heights = [counts[clone].get(cat, 0) for clone in clones]
                ax.bar(section_positions, heights, bottom=bottom, label=cat if current_x == 0 else None, color=colors[cat])
                bottom = [sum(x) for x in zip(bottom, heights)]

        current_x += len(bases) + 2 if "Fixed" in section_name else len(clones) + 2  # Add spacing between sections

    # Add section labels below the x-axis
    section_midpoints = [(start + end) // 2 for start, end in zip(section_starts, section_ends)]
    ax.set_xticks(section_midpoints)

    # Add lines and labels below x-axis sections
    ax.text((section_starts[0] + section_ends[1]) // 2, -0.1, "Fixed Mutations", ha='center', fontsize=16, transform=ax.get_xaxis_transform())
    ax.text((section_starts[2] + section_ends[3]) // 2, -0.1, "Non-Fixed Mutations", ha='center', fontsize=16, transform=ax.get_xaxis_transform())
    ax.set_xticklabels(["Asexual Populations", "Sexual Populations", "Asexual Clones", "Sexual Clones"], fontsize=14, rotation=0)
    ax.tick_params(axis='x', length=0, pad=15)

    # Add legend in the correct order
    handles, labels = ax.get_legend_handles_labels()
    order = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
    ordered_handles = [handles[labels.index(cat)] for cat in order]
    ordered_labels = [cat.capitalize() for cat in order]
    ax.legend(ordered_handles, ordered_labels, loc='center', bbox_to_anchor=(0.095, 0.85), fontsize=14, frameon=False)

    # y axis
    ax.tick_params(axis='y', labelsize=12)
    ax.set_ylim(0, 120)
    ax.set_ylabel("Mutation Count", fontsize=16)

    # Remove top, bottom, and right borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.tight_layout()
    plt.savefig("figureS4.pdf", dpi=300)
    plt.savefig("figureS4.jpg", dpi=300)
    #plt.close(fig)
    return fig

generate_figureS4(file_path)



##### STATS #####

### FIGURE 3A: Stacked bar plot of fixed & non-fixed mutations across sexual & asexual clones

# Print summary stats
print("\n--- Mutation Count Stats for Figure 3A ---")
for stat in summary_stats:
    print(f"{stat['label']}: mean = {stat['mean']:.2f}, median = {stat['median']:.2f}, n = {len(stat['counts'])}")
print("------------------------------------------\n")

# Print Shapiro-Wilk results
print("--- Shapiro-Wilk Normality Tests (with BH correction) ---")
for result in shapiro_stats:
    print(f"{result['label']}: W = {result['statistic']:.4f}, "
          f"raw p = {result['pvalue_raw']:.4f}, "
          f"BH-corrected p = {result['pvalue_corrected']:.4f}, "
          f"{'REJECT' if result['rejected'] else 'n.s.'}")
print("----------------------------------------------------------\n")

# Print Brunner-Munzel results
print("--- Brunner-Munzel Test Results (with BH correction) ---")
for result in brunner_results:
    print(f"{result['label']}: BM stat = {result['statistic']:.4f}, "
          f"raw p = {result['pvalue_raw']:.4f}, "
          f"BH-corrected p = {result['pvalue_corrected']:.4f}, "
          f"{'REJECT' if result['rejected'] else 'n.s.'}")
print("--------------------------------------------------------\n")

### FIGURE 3B: Number of putatively functional mutations in sexual & asexual clones

print("\n--- Shapiro-Wilk Normality Tests for Figure 3B (with BH correction) ---")
for res in shapiro_results_3b:
    print(f"{res['label']}: W = {res['statistic']:.4f}, "
          f"raw p = {res['pvalue']:.4f}, "
          f"BH-corrected p = {res['pvalue_corrected']:.4f}, "
          f"{'REJECT' if res['rejected'] else 'n.s.'} "
          f"(n = {res['n']}, mean = {res['mean']:.2f}, median = {res['median']:.2f})")
print("-----------------------------------------------------------------------\n")

print("--- Welch's t-test (Asexual vs Sexual) ---")
print(f"t = {welch_result_3b['t_statistic']:.4f}, p = {welch_result_3b['pvalue']:.4f}")
print(f"n (Asexual) = {welch_result_3b['n_asexual']}, mean = {welch_result_3b['mean_asexual']:.2f}")
print(f"n (Sexual)  = {welch_result_3b['n_sexual']}, mean = {welch_result_3b['mean_sexual']:.2f}")
print("------------------------------------------\n")

### FIGURE 3C: Fraction of functional mutations in multi-hit genes

# Print Shapiro-Wilk results
print("\n--- Shapiro-Wilk Normality Tests for Figure 3C (with BH correction) ---")
for res in shapiro_results_3d:
    print(f"{res['label']}: W = {res['statistic']:.4f}, "
          f"raw p = {res['pvalue']:.4f}, "
          f"BH-corrected p = {res['pvalue_corrected']:.4f}, "
          f"{'REJECT' if res['rejected'] else 'n.s.'} "
          f"(n = {res['n']}, mean = {res['mean']:.2f}, median = {res['median']:.2f})")
print("------------------------------------------------------------------------\n")

# Print Welch's t-test result
print("--- Welch's t-test (Asexual vs Sexual) ---")
print(f"t = {welch_result_3d['t_statistic']:.4f}, p = {welch_result_3d['pvalue']:.4f}")
print(f"n (Asexual) = {welch_result_3d['n_asexual']}, mean = {welch_result_3d['mean_asexual']:.2f}")
print(f"n (Sexual)  = {welch_result_3d['n_sexual']}, mean = {welch_result_3d['mean_sexual']:.2f}")
print("------------------------------------------\n")