import sys
import pandas as pd

def main():
    # Check if the input file is provided
    if len(sys.argv) != 2:
        print("Usage: python change_headers.py <input_tsv_file>")
        sys.exit(1)

    # Input and output files
    input_file = sys.argv[1]
    output_file = "GGE_clones_final.tsv"
    meta_data_file = "../../meta_data_for_mutations.txt"

    # Load the metadata
    meta_data = pd.read_csv(meta_data_file, sep='\t')
    file_to_clone = dict(zip(meta_data['file_name'], meta_data['clone_id']))
    clone_to_sex = dict(zip(meta_data['clone_id'], meta_data['SEX']))

    # Load the GGE_clones TSV
    gge_clones = pd.read_csv(input_file, sep='\t')

    # Update column names
    new_columns = []
    for col in gge_clones.columns:
        updated = False
        for file_name in file_to_clone.keys():
            if col.startswith(file_name):
                clone_id = file_to_clone[file_name]
                sex_type = clone_to_sex.get(clone_id, "UNKNOWN")  # Fallback to UNKNOWN if not found
                suffix = col[len(file_name):]  # Extract suffix
                # Rearrange to have clone_id, SEX, and suffix
                new_columns.append(f"{clone_id}_{sex_type}{suffix}")
                updated = True
                break
        if not updated:
            new_columns.append(col)

    gge_clones.columns = new_columns

    # Save the updated file
    gge_clones.to_csv(output_file, sep='\t', index=False)

    print(f"Updated TSV saved as '{output_file}'")

if __name__ == "__main__":
    main()
