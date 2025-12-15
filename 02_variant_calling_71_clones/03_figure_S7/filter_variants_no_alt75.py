import pandas as pd
import sys
import os

# Chromosome lengths from w303_gge_corr.fasta
chromo_lengths = {
    'chrI': 230207, 'chrII': 813124, 'chrIII': 317068, 'chrIV': 1531751, 'chrV': 576825,
    'chrVI': 270161, 'chrVII': 1090940, 'chrVIII': 562643, 'chrIX': 439888, 'chrX': 745751,
    'chrXI': 666816, 'chrXII': 1078715, 'chrXIII': 924431, 'chrXIV': 784333, 'chrXV': 1091291,
    'chrXVI': 948066
}

# Telomere lengths provided
telomere_lengths = {
    'TEL01L': 801, 'TEL01R': 808, 'TEL02L': 6608, 'TEL02R': 806, 'TEL03L': 1098, 'TEL03R': 838,
    'TEL04L': 904, 'TEL04R': 7309, 'TEL05L': 6473, 'TEL05R': 7276, 'TEL06L': 5530, 'TEL06R': 431,
    'TEL07L': 781, 'TEL07R': 7306, 'TEL08L': 5505, 'TEL08R': 6539, 'TEL09L': 7784, 'TEL09R': 821,
    'TEL10L': 7767, 'TEL10R': 850, 'TEL11L': 807, 'TEL11R': 913, 'TEL12L': 12085, 'TEL12R': 13897,
    'TEL13L': 6344, 'TEL13R': 891, 'TEL14L': 7428, 'TEL14R': 1056, 'TEL15L': 847, 'TEL15R': 7370,
    'TEL16L': 7223, 'TEL16R': 5615
}

# Mapping for Roman numerals to numeric codes
romans_to_numbers = {
    'I': '01', 'II': '02', 'III': '03', 'IV': '04', 'V': '05', 'VI': '06', 
    'VII': '07', 'VIII': '08', 'IX': '09', 'X': '10', 'XI': '11', 'XII': '12',
    'XIII': '13', 'XIV': '14', 'XV': '15', 'XVI': '16'
}

def generate_low_coverage_blocks(coverage_file, coverage_threshold=710): #average 10x coverage for 71 clones
    """Generate low coverage blocks from coverage data."""
    # Load the coverage data
    coverage_data = pd.read_csv(coverage_file, sep="\t", header=None, names=["Chromosome", "Position", "Coverage"])

    # Filter for regions with coverage below the threshold
    low_coverage = coverage_data[coverage_data["Coverage"] < coverage_threshold]

    # Group adjacent positions into blocks
    grouped_regions = []
    current_chrom = None
    start = None
    end = None

    for _, row in low_coverage.iterrows():
        chrom = row["Chromosome"]
        pos = row["Position"]

        if chrom != current_chrom or (end is not None and pos != end + 1):
            # Save the previous block if it exists
            if current_chrom is not None:
                grouped_regions.append([current_chrom, start, end])

            # Start a new block
            current_chrom = chrom
            start = pos
            end = pos
        else:
            # Extend the current block
            end = pos

    # Add the last block
    if current_chrom is not None:
        grouped_regions.append([current_chrom, start, end])

    # Convert to a DataFrame
    grouped_df = pd.DataFrame(grouped_regions, columns=["Chromosome", "Start", "End"])
    return grouped_df

def remove_mitochondrial_mutations(mutations):
    """Filter out mitochondrial mutations."""
    return mutations[mutations['CHROM'] != 'chrMito']

def telomere_filtering(mutations):
    """Filter out telomeric mutations."""
    def in_telomere(row):
        chromo = row['CHROM']
        pos = row['POS']
        if chromo in chromo_lengths:
            roman_num = chromo[3:].upper()  # Extract Roman numeral
            numeric_code = romans_to_numbers.get(roman_num, None)
            if numeric_code:
                left_edge = telomere_lengths[f'TEL{numeric_code}L']
                right_edge = chromo_lengths[chromo] - telomere_lengths[f'TEL{numeric_code}R']
                return pos < left_edge or pos > right_edge
        return False

    # Apply telomere filtering
    mutations['in_telomere'] = mutations.apply(in_telomere, axis=1)
    return mutations[~mutations['in_telomere']]

def multipopulation_filtering(mutations, meta_data_file):
    """Filter out mutations present in 2 or more populations."""
    # Load metadata
    meta_data = pd.read_csv(meta_data_file, sep="\t")

    # Identify `alt_counts` columns
    alt_columns = [col for col in mutations.columns if col.endswith("_alt_counts")]

    # Map `file_name` to `pop_id`
    file_to_pop = dict(zip(meta_data["file_name"], meta_data["pop_id"]))
    clone_to_population = {col: file_to_pop.get(col.replace("_alt_counts", ""), None) for col in alt_columns}

    # Group `alt_counts` columns by populations
    population_counts = {}
    for col, population in clone_to_population.items():
        if population:
            if population not in population_counts:
                population_counts[population] = []
            population_counts[population].append(col)

    # Count populations with non-zero alt_counts
    def count_nonzero_populations(row):
        nonzero_pops = 0
        for population, cols in population_counts.items():
            if any(row[col] > 0 for col in cols):  # Check if any clone in this population has alt_counts > 0
                nonzero_pops += 1
        return nonzero_pops

    # Add a column and filter
    mutations["pop_nonzero_count"] = mutations.apply(count_nonzero_populations, axis=1)
    return mutations[mutations["pop_nonzero_count"] < 2]

def low_coverage_filtering(mutations, low_coverage_blocks):
    """Filter out variants in low coverage regions."""
    # Convert low coverage blocks to a dictionary for efficient lookup
    low_coverage_dict = {}
    for _, row in low_coverage_blocks.iterrows():
        if row["Chromosome"] not in low_coverage_dict:
            low_coverage_dict[row["Chromosome"]] = []
        low_coverage_dict[row["Chromosome"]].append((row["Start"], row["End"]))

    # Function to check if a position overlaps with any range in the low coverage dictionary
    def is_low_coverage(chrom, pos):
        if chrom in low_coverage_dict:
            for start, end in low_coverage_dict[chrom]:
                if start <= pos <= end:
                    return True
        return False

    # Apply the filtering
    mutations["in_low_coverage"] = mutations.apply(
        lambda row: is_low_coverage(row["CHROM"], row["POS"]), axis=1
    )
    return mutations[~mutations["in_low_coverage"]]

def low_alt_counts_filtering(mutations):
    """Filter out variants with low alt_counts (<5) across every clone."""
    # Identify `alt_counts` columns
    alt_columns = [col for col in mutations.columns if col.endswith("_alt_counts")]

    # Function to check if all alt_counts are below 5 for this variant
    def all_low_alt_counts(row):
        return all((row[col] < 5 if not pd.isna(row[col]) else True) for col in alt_columns)

    # Set alt_counts < 5 to 0 across all rows
    for col in alt_columns:
        mutations.loc[mutations[col] < 5, col] = 0

    # Apply the filtering
    mutations["low_alternate_counts"] = mutations.apply(all_low_alt_counts, axis=1)
    return mutations[~mutations["low_alternate_counts"]]

def alt_freq_filtering(mutations):
    """Filter out variants where all valid alt frequencies are below 0.75."""
    # Identify alt_counts and ref_counts columns
    alt_columns = [col for col in mutations.columns if col.endswith("_alt_counts")]
    ref_columns = [col.replace("_alt_counts", "_ref_counts") for col in alt_columns]

    # Create a DataFrame to store alt frequencies
    alt_frequencies_df = pd.DataFrame(index=mutations.index)

    # Calculate alt frequencies
    for alt_col, ref_col in zip(alt_columns, ref_columns):
        alt_freq_col = alt_col.replace("_alt_counts", "_alt_freq")
        alt_frequencies_df[alt_freq_col] = 0.0  # Default to 0.0
        valid_rows = mutations[alt_col] > 0
        alt_frequencies_df.loc[valid_rows, alt_freq_col] = (
            mutations.loc[valid_rows, alt_col].astype(float) /
            (mutations.loc[valid_rows, alt_col] + mutations.loc[valid_rows, ref_col])
        )

    # Create a mask to filter rows where all valid alt frequencies are >= 0.75
    valid_mask = (alt_frequencies_df >= 0.75) | (alt_frequencies_df == 0.0)  # Keep rows with 0 frequencies
    rows_to_keep = valid_mask.all(axis=1)  # True if all conditions are met for a row

    # Filter the mutations DataFrame
    return mutations[rows_to_keep]


def main():
    # Check if the user provided an input file
    if len(sys.argv) != 2:
        print("Usage: python filter_variants.py <mutations_file>")
        sys.exit(1)

    # Input file path
    mutations_file = sys.argv[1]

    # Metadata file
    meta_data_file = "../../meta_data_for_mutations.txt"

    # Coverage file
    coverage_file = "../../01_reference_correction_71_clones/coverage.txt"

    # Ensure input file exists
    if not os.path.exists(mutations_file):
        print(f"Error: File '{mutations_file}' does not exist.")
        sys.exit(1)

    # Load the mutations file
    mutations = pd.read_csv(mutations_file, sep="\t")

    print(f"Total rows before any filtering: {mutations.shape[0]}")

    # Step 1: Remove mitochondrial mutations
    mutations = remove_mitochondrial_mutations(mutations)
    print(f"Rows after mitochondrial filtering: {mutations.shape[0]}")

    # Step 2: Apply telomere filtering
    mutations = telomere_filtering(mutations)
    print(f"Rows after telomere filtering: {mutations.shape[0]}")

    # Step 3: Apply multipopulation filtering
    mutations = multipopulation_filtering(mutations, meta_data_file)
    print(f"Rows after multipopulation filtering: {mutations.shape[0]}")

    # Step 4: Generate low coverage blocks and apply filtering
    low_coverage_blocks = generate_low_coverage_blocks(coverage_file)
    mutations = low_coverage_filtering(mutations, low_coverage_blocks)
    print(f"Rows after low coverage filtering: {mutations.shape[0]}")

    # # Step 5: Apply alt_freq filtering
    # mutations = alt_freq_filtering(mutations)
    # print(f"Rows after alt_freq filtering: {mutations.shape[0]}")
    
    # Step 6: Apply low alt_counts filtering
    mutations = low_alt_counts_filtering(mutations)
    print(f"Rows after low alt_counts filtering: {mutations.shape[0]}")

    # Save the filtered mutations
    output_file = os.path.splitext(mutations_file)[0] + "_filtered.tsv"
    mutations.to_csv(output_file, sep="\t", index=False)

    print(f"Filtered mutations saved to {output_file}")

if __name__ == "__main__":
    main()