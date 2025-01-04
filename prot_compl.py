"""
@author: Royal Truman

Find the molecular complexes in which the same protein is used.
First, extract the necessary data from the CORUM 5.0 database, https://mips.helmholtz-muenchen.de/corum/download
Download all annotated human protein complexes in *.txt format. Only the columns we need will be extracted from the
file downloaded, corum_humanComplexes.txt; specifically, the complexes and the proteins which comprise them.
This script identifies all te proteins in the dataset and for each identifies the complexes the protein is found in.
The list of complexes are output separated by |.
"""

from pandas import read_csv, errors

def process_corum_data(input_file, output_sma_file, output_complexes_file):
    """
    Processes CORUM data to extract human complexes and their subunit information.

    Args:
        input_file: Path to the input CORUM data file.
        output_sma_file: Path to the output file containing complex_id, complex_name, and subunits_uniprot_id for human complexes.
        output_complexes_file: Path to the output file containing subunits_uniprot_id and their associated complex_ids.
    """

    try:
        # 1. Read data and select human records
        df = read_csv(input_file, sep='\t')
        human_df = df[df['organism'] == 'Human']  # Remove bad recs from downloaded file

        # 2. Extract and save selected columns WITH HEADER
        sma_df = human_df[['complex_id', 'complex_name', 'subunits_uniprot_id']]
        sma_df.to_csv(output_sma_file, sep='\t', index=False, header=True)

        # 3 & 4. Process and output complex-protein associations
        complexes_dict = {}  # The dict keys will be the protein IDs, and the values lists of complex IDs
        with open(output_sma_file, 'r') as sma_file:
            next(sma_file)  # Skip header line which was now added.
            for line in sma_file:
                # Split each line into complex_id, complex_name (ignored using _), and subunits_uniprot_id). line.strip() removes leading/trailing whitespace
                complex_id, _, subunits_uniprot_id = line.strip().split('\t')
                # Handles cases where there are multiple subunits listed in the subunits_uniprot_id column, separated by semicolons.
                for subunit in subunits_uniprot_id.split(';'):
                    subunit = subunit.strip()
                    if subunit:  # Check if subunit is not empty
                        # If the current subunit is not already a key in the dict create a new entry
                        if subunit not in complexes_dict:
                            complexes_dict[subunit] = []
                        complexes_dict[subunit].append(complex_id)

        # Write output with headers
        with open(output_complexes_file, 'w') as complexes_file:
            complexes_file.write("subunits_uniprot_id\tcomplex_ids\n")
            # Sort the subunits alphabetically before writing them to the file.
            sorted_subunits = sorted(complexes_dict.keys())
            for subunit in sorted_subunits:
                # Join the list of complex IDs associated with the current protein into a single string separated |
                complex_ids = "|".join(complexes_dict[subunit])
                complexes_file.write(f"{subunit}\t{complex_ids}\n")

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
    except errors.ParserError:
        print(f"Error: Could not parse input file '{input_file}'. Check file format.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Example usage (same as before)
input_file = "corum_humanComplexes.txt"
output_sma_file = "corum_humanComplexes_sma.txt"
output_complexes_file = "Complexes_prot_found_in.txt"

process_corum_data(input_file, output_sma_file, output_complexes_file)
print("Processing complete. Check output files.")
