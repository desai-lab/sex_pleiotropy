import pandas as pd
import numpy as np
import os
import sys

def group_variants_optimized(df, alt_columns, max_distance=25, freq_tolerance=5):
    """
    Group nearby variants within a population based on proximity and allele frequency similarity.

    Parameters:
        df (pd.DataFrame): Input DataFrame with mutation data.
        alt_columns (list): List of columns containing allele counts.
        max_distance (int): Maximum distance (in base pairs) for grouping.
        freq_tolerance (float): Tolerance for allele frequency similarity.

    Returns:
        pd.Series: A Pandas Series containing group IDs for each mutation.
    """
    # Backup original values to restore later
    alt_original = df[alt_columns].copy()

    # Convert alt_columns to numeric, ensuring no SettingWithCopyWarning
    df.loc[:, alt_columns] = df[alt_columns].apply(pd.to_numeric, errors="coerce")

    # Sort DataFrame by CHROM and POS
    df = df.sort_values(["CHROM", "POS"]).reset_index(drop=True)

    # Initialize group IDs and tracking variables
    group_ids = np.zeros(len(df), dtype=int)
    current_group = 1

    # Track the starting position of the current group
    start_index = 0

    for i in range(1, len(df)):
        # Check if the current row can be grouped with the starting row
        within_distance = (
            df.loc[i, "CHROM"] == df.loc[start_index, "CHROM"] and
            abs(df.loc[i, "POS"] - df.loc[start_index, "POS"]) <= max_distance
        )

        # Ensure valid numeric arrays for comparison
        try:
            row_values = df.loc[i, alt_columns].to_numpy(dtype=float)
            start_values = df.loc[start_index, alt_columns].to_numpy(dtype=float)
        except ValueError as e:
            print(f"Error converting row or start_values to numeric at index {i}: {e}")
            row_values = np.full(len(alt_columns), np.nan)
            start_values = np.full(len(alt_columns), np.nan)

        similar_frequencies = np.allclose(
            row_values,
            start_values,
            atol=freq_tolerance,
            equal_nan=True  # Treat NaN values as equal
        )

        if within_distance and similar_frequencies:
            # Continue grouping
            group_ids[i] = current_group
        else:
            # Start a new group
            current_group += 1
            group_ids[i] = current_group
            start_index = i  # Update the starting row for the new group

    # Restore original `alt_columns` values with blanks (NaN)
    df.loc[:, alt_columns] = alt_original

    return pd.Series(group_ids, index=df.index)

def main():
    # Check for the required input file argument
    if len(sys.argv) != 2:
        print("Usage: python group_variants.py <mutations_file>")
        return

    # Input file from the command line
    mutations_file = sys.argv[1]
    meta_data_file = "../../meta_data_for_mutations.txt"

    # Ensure input file exists
    if not os.path.exists(mutations_file):
        print(f"Error: File '{mutations_file}' does not exist.")
        return

    # Ensure metadata file exists
    if not os.path.exists(meta_data_file):
        print(f"Error: Metadata file '{meta_data_file}' does not exist.")
        return

    # Load the mutations and metadata files
    mutations = pd.read_csv(mutations_file, sep="\t")
    meta_data = pd.read_csv(meta_data_file, sep="\t")

    # Identify `alt_counts` columns in the mutations file
    alt_columns = [col for col in mutations.columns if col.endswith("_alt_counts")]

    # Map `file_name` to `pop_id` and create `clone_to_population`
    file_to_pop = dict(zip(meta_data["file_name"], meta_data["pop_id"]))
    clone_to_population = {col: file_to_pop.get(col.replace("_alt_counts", ""), None) for col in alt_columns}

    # Group `alt_counts` columns by populations
    population_counts = {}
    for col, population in clone_to_population.items():
        if population:
            if population not in population_counts:
                population_counts[population] = []
            population_counts[population].append(col)

    # Apply the optimized grouping function for each population
    mutations["mutation_group"] = -1  # Initialize a column for mutation groups
    for population, cols in population_counts.items():
        # Extract relevant data for this population
        population_data = mutations[cols + ["CHROM", "POS"]]

        # Group the variants
        group_ids = group_variants_optimized(population_data, cols)

        # Assign group IDs back to the main DataFrame
        mutations.loc[population_data.index, "mutation_group"] = group_ids

    # Save grouped mutations
    output_file = os.path.splitext(mutations_file)[0] + "_grouped.tsv"
    mutations.to_csv(output_file, sep="\t", index=False)

    print(f"Mutation grouping complete. Results saved to {output_file}")

if __name__ == "__main__":
    main()
