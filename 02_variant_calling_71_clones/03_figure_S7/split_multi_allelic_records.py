import pandas as pd
from collections import defaultdict
import sys
import os

def annotation_splitter(a):
    # input is full annotation column
    # output is dictionary from ALT alleles (strings) to lists of annotations:
    anns = defaultdict(list)
    for annotation in a.split(','):
        alt_text = annotation.split('|')[0]
        anns[alt_text].append(annotation)
    return anns

def process_file(input_file):
    # Extract base name without extension for output
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f"{base_name}_split.tsv"

    # Base columns and initialization
    base_cols = ['CHROM', 'POS', 'REF', 'QUAL']
    dfs = []

    # Getting column name info
    nd = pd.read_csv(input_file, delimiter='\t', dtype=str)  # Read all columns as strings
    other_cols = [i for i in nd.columns if i not in base_cols + ['ALT', 'ANN']]
    sample_names = [i.split('.AD')[0] for i in other_cols]
    all_cols = base_cols + ['ALT', 'ANN'] + other_cols
    new_col_names = base_cols + ['ALT', 'ANN'] + [i + '_ref_counts' for i in sample_names] + [i + '_alt_counts' for i in sample_names]

    # Process the file, splitting multi-allelic records
    new_rows = []

    for _, row in nd.iterrows():
        alt = row['ALT']
        if ',' not in alt:  # Only one alternate allele
            new_rows.append(
                list(row[base_cols + ['ALT', 'ANN']]) +
                [
                    int(row[c].split(',')[0]) if pd.notna(row[c]) else None  # Reference count
                    for c in other_cols
                ] +
                [
                    int(row[c].split(',')[1]) if pd.notna(row[c]) else None  # Alternate count
                    for c in other_cols
                ]
            )
        else:
            ad = annotation_splitter(row['ANN'])  # Dictionary from ALT alleles to annotations
            alt_split = alt.split(',')
            base_row = list(row[base_cols])
            for i, current_alt in enumerate(alt_split):
                current_ann = ','.join(ad[current_alt])
                new_rows.append(
                    base_row + [current_alt, current_ann] +
                    [
                        int(row[c].split(',')[0]) if pd.notna(row[c]) else None  # Reference count
                        for c in other_cols
                    ] +
                    [
                        int(row[c].split(',')[i + 1]) if pd.notna(row[c]) else None  # Alternate count
                        for c in other_cols
                    ]
                )

    dfs.append(pd.DataFrame(new_rows, columns=new_col_names))

    nd = pd.concat(dfs)  # Combine DataFrames if multiple (just one in this case)

    # Remove rows where ALT = '*'
    nd = nd[nd['ALT'] != '*']

    # Save the output
    nd.to_csv(output_file, sep='\t', index=False)
    print(f"Processed file saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python split_multi_allelic_records.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    process_file(input_file)
