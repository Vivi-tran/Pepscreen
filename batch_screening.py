import os
import subprocess
import pandas as pd
import sys
import logging

# Clear existing handlers
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

# Configure logging
logging.basicConfig(
        filename = 'log.txt',
        filemode= 'w',
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
)
# Function to run bash commands
def run_command(command):
    process = subprocess.run(command, shell=True, text=True, capture_output=True)
    if process.returncode != 0:
        logging.info(f"Error occurred while running command: {command}")
        logging.info(process.stderr)
        sys.exit(1)
    return process.stdout

# Function to process peptide sequences and predict binding affinity
def process_peptides(input_csv, protein_sequence, output_csv):
    # Read input CSV
    df = pd.read_csv(input_csv)

    if 'id' not in df.columns or 'pep_seq' not in df.columns:
        logging.error("Input CSV must contain 'id' and 'pep_seq' columns.")
        sys.exit(1)

    # Prepare output list
    results = []

    for _, row in df.iterrows():
        peptide_id = row['id']
        peptide_sequence = row['pep_seq']

        logging.info(f"Processing peptide ID: {peptide_id}, Sequence: {peptide_sequence}")

        # Run feature extraction scripts
        run_command(f"python ./softwares/iupred2a/gene_intrinsic_disorder.py {protein_sequence} {peptide_sequence}")
        run_command(f"python ./codes/gen_env_padding_feature.py {protein_sequence} {peptide_sequence}")
        run_command(f"python ./codes/Physicochemical_characteristics.py {protein_sequence} {peptide_sequence}")
        run_command(f"python ./codes/SP_score.py {protein_sequence} {peptide_sequence}")

        # Run prediction script and capture the output
        prediction_output = run_command(f"python ./codes/predicted.py")
        # Process the prediction output
        try:
            # Split the output based on the lines, assuming there are two relevant lines
            output_lines = prediction_output.strip().splitlines()

            # Extract values from each line
            pred_pK_line = output_lines[0]
            pred_dG_line = output_lines[1]

            # Extract numeric values from the lines
            pred_pK = float(pred_pK_line.split(":")[1].strip())
            pred_dG = float(pred_dG_line.split(":")[1].strip())

            # Round the values
            pred_pK = round(pred_pK, 3)
            pred_dG = round(pred_dG, 3)
        except Exception as e:
            logging.error(f"Error parsing prediction output for peptide ID {peptide_id}: {prediction_output.strip()}")
            logging.error(f"Error details: {e}")
            sys.exit(1)

        results.append({
            'id': peptide_id,
            'peptide_sequence': peptide_sequence,
            'pred_pK': pred_pK,
            'pred_dG': pred_dG
        })

    # Save results to CSV
    output_df = pd.DataFrame(results)
    output_df.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")

# Main function
def main():
    if len(sys.argv) != 3:
        print("Usage: python process_peptides.py <input_csv> <output_csv>")
        logging.error("Usage: python process_peptides.py <input_csv> <output_csv>")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_csv = sys.argv[2]

    # Define protein sequence directly in the script
    protein_sequence = "SSPSEGLCPPGHHISEDGRDCISCKYGQDYSTHWNDLLFCLRCTRCDSGEVELSPCTTTRNTVCQCEEGTFREEDSPEMCRKCRTGCPRGMVKVGDCTPWSDIECVHKES"

    process_peptides(input_csv, protein_sequence, output_csv)

if __name__ == "__main__":
    main()
