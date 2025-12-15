import pandas as pd
import os
import argparse

# Required mapping for mutation types
mut_type_simple = {
    'missense': ['missense_variant'],
    'nonsense': ['stop_lost', 'stop_gained', 'start_lost'],
    'synonymous': ['synonymous_variant', 'stop_retained_variant', 'initiator_codon_variant'],
    'noncoding': ['upstream_gene_variant', 'splice_region_variant', 'intron_variant', 'splice_acceptor_variant', '', 
                  'intergenic_region', 'downstream_gene_variant'],
    'indel': ['conservative_inframe_insertion', 'conservative_inframe_deletion', 'disruptive_inframe_insertion', 'disruptive_inframe_deletion', 
           'frameshift_variant'],
}

# Create a mutation simplification dictionary
mu_simplify = dict()
for m in mut_type_simple:
    mu_simplify.update({i: m for i in mut_type_simple[m]})

# Function definitions for annotation parsing
def simplify_mutation_type(m):
    mut_type_list = [mu_simplify[i] for i in str(m).split('&')]
    mut_types_in_order = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
    mut_types_here = [mt for mt in mut_types_in_order if mt in mut_type_list]
    if len(mut_types_here) > 0:
        return mut_types_here[0]  # Returns the first significant mutation type (priority order)
    else:
        return 'NA'  # Default if no known type is found

def simple_ann(a):
    sa = a.split('|')
    
    # Determine the mutation type
    mut_type = simplify_mutation_type(sa[1])
    
    # Extract amino acid position unless the mutation type is 'noncoding'
    if mut_type != 'noncoding' and len(sa) > 13:
        aa_pos = sa[13]
    else:
        aa_pos = ''  # Exclude amino acid position for noncoding mutations

    # Extract specific amino acid change if present
    if len(sa) > 11 and sa[10]:
        aa_change = sa[10]  # Example: p.Tyr14Cys
    else:
        aa_change = ''
    
    # Construct the simplified annotation
    return '|'.join([str(i).replace(',', '_') for i in [sa[0], mut_type, sa[2], sa[3], sa[4], aa_pos, aa_change]])

def annotation_parser(row):
    # Simplifies annotations in the `ANN` column based on given rules
    annotations = row['ANN']
    anns = []
    for a in str(annotations).split(','):
        sa = a.split('|')
        if len(sa) > 1:
            anns.append(simple_ann(a))  # Process each annotation
    
    # Define the order of importance for impacts and mutation types
    impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    mut_type_order = ['indel', 'nonsense', 'missense', 'synonymous', 'noncoding']
    
    # Sort annotations by impact, then by mutation type
    anns_sorted = sorted(
        anns,
        key=lambda x: (
            impact_order.index(x.split('|')[2]) if x.split('|')[2] in impact_order else len(impact_order),
            mut_type_order.index(x.split('|')[1]) if x.split('|')[1] in mut_type_order else len(mut_type_order)
        )
    )
    
    # Return only the most significant annotation
    return anns_sorted[0] if anns_sorted else ''

# Command-line input handling
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simplify annotations in a TSV file.")
    parser.add_argument("input_file", help="Path to the input TSV file.")
    args = parser.parse_args()

    # Extract input and output file paths
    input_file = args.input_file
    base_name = os.path.splitext(input_file)[0]
    output_file = f"{base_name}_simplified.tsv"

    # Read the input TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Ensure the required columns are present
    if not all(col in df.columns for col in ['ALT', 'ANN']):
        raise ValueError("The required columns ['ALT', 'ANN'] are missing from the input file.")

    # Apply the annotation parser
    df['ANN_simpler'] = df.apply(annotation_parser, axis=1)

    # Save the updated DataFrame to a new file
    df.to_csv(output_file, sep='\t', index=False)

    print(f"Simplified annotations saved to {output_file}")
