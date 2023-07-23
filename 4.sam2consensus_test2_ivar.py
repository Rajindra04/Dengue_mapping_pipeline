import argparse
import glob
import os
import subprocess
import matplotlib.pyplot as plt

# Usage: python sam_to_consensus.py --input_dir input_dir --reference_fasta reference.fasta --output_dir output_dir 

def run_command(command):
    print("Running command:", command)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.wait()

    stdout = process.stdout.read().decode("utf-8")
    stderr = process.stderr.read().decode("utf-8")

    if process.returncode != 0:
        error_message = stderr.strip()
        raise Exception("Command execution failed with return code {}: {}".format(process.returncode, error_message))

    return stdout, stderr

def plot_coverage(coverage_file, output_dir):
    # Load coverage data from the file
    positions = []
    depths = []
    with open(coverage_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            position = int(columns[1])
            depth = int(columns[2])
            positions.append(position)
            depths.append(depth)

    # Create the coverage plot
    plt.plot(positions, depths)
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title('Read Coverage')
    plt.grid(True)

    # Save the coverage plot as an image file
    coverage_plot_file = os.path.join(output_dir, os.path.splitext(os.path.basename(coverage_file))[0] + '.png')
    plt.savefig(coverage_plot_file)
    plt.close()

    print("Coverage plot generated:", coverage_plot_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing BAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for consensus FASTA files.')
    args = parser.parse_args()

    input_dir = args.input_dir
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Find all BAM files in the input directory
    bam_files = glob.glob(os.path.join(input_dir, "*sorted.bam"))

    for bam_file in bam_files:
        print("Processing:", bam_file)

        # Get the sample name from the filename
        sample_name = os.path.splitext(os.path.basename(bam_file))[0]

        try:
            # Calculate read coverage using samtools depth
            coverage_file = f"{sample_name}_coverage.txt"
            coverage_qual_file = f"{sample_name}_qual.txt"
            coverage_command = f"samtools depth {bam_file} > {os.path.join(output_dir, coverage_file)} && samtools depth -Q 20 {bam_file} > {os.path.join(output_dir, coverage_qual_file)}"
            stdout, stderr = run_command(coverage_command)
            print("Command output:", stdout)
            print("Command error:", stderr)

            # Generate the coverage plot and save as an image file
            plot_coverage(os.path.join(output_dir, coverage_file), output_dir)

            # Generate pileup using BAM file and pipe it to ivar consensus
            pileup_command = f"samtools mpileup -d 1000 -A -Q 0 {bam_file} | ivar consensus -p {sample_name} -q 20 -t 0"
            stdout, stderr = run_command(pileup_command)
            print("Command output:", stdout)
            print("Command error:", stderr)

            # Move the consensus file to the output directory
            consensus_fasta = f"{sample_name}.fa"
            os.rename(consensus_fasta, os.path.join(output_dir, consensus_fasta))

            print("Consensus sequence, coverage file, and coverage plot generated:", os.path.join(output_dir, consensus_fasta), os.path.join(output_dir, coverage_file))

        except Exception as e:
            print("Error occurred during processing:", str(e))
            continue
