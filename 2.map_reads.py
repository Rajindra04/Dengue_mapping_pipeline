import argparse
import pandas as pd
import subprocess
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
##usage- python 2.map_reads.py --samplesheet samplesheet.tsv
def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception("Command execution failed with return code {}, stderr: {}".format(process.returncode, stderr.decode("utf-8")))
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def samplesheet_verify(samplesheet):
    try:
        samples = pd.read_csv(samplesheet, sep='\t')
        if not all(col in samples.columns for col in ['sample_name', 'read1', 'read2']):
            raise ValueError("Sample sheet must have columns 'sample_name', 'read1', 'read2'.")
    except Exception as e:
        raise ValueError("Error reading sample sheet: {}".format(e))
    return samples

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--samplesheet', type=str, help='Path to sample sheet with columns "sample_name", "read1", "read2".')
    args = parser.parse_args()

    samples = samplesheet_verify(args.samplesheet)
    base_dir = os.path.abspath(os.path.dirname(__file__))
    
        
    # Reference indexing with bwa-mem2
    index_command = f"bwa-mem2 index {base_dir}/references/denv1.fasta"
    print(f"Indexing reference denv1.fasta...")
    stdout, stderr = run_command(index_command)
    print(stdout)
    print(stderr)

    index_command = f"bwa-mem2 index {base_dir}/references/denv2.fasta"
    print(f"Indexing reference denv2.fasta...")
    stdout, stderr = run_command(index_command)
    print(stdout)
    print(stderr)

    index_command = f"bwa-mem2 index {base_dir}/references/denv3.fasta"
    print(f"Indexing reference denv3.fasta...")
    stdout, stderr = run_command(index_command)
    print(stdout)
    print(stderr)

    index_command = f"bwa-mem2 index {base_dir}/references/denv4.fasta"
    print(f"Indexing reference denv4.fasta...")
    stdout, stderr = run_command(index_command)
    print(stdout)
    print(stderr)
           
    for i, row in samples.iterrows():
        sample_name = row['sample_name']
        read1 = row['read1']
        read2 = row['read2']

        sample_folder = f"{sample_name}_output"
        os.makedirs(sample_folder, exist_ok=True)
        os.chdir(sample_folder)

        # Trimmomatic
        trimmomatic_command = f" fastp -i {read1} -I {read2} -o {sample_name}_trimmed_R1.fastq.gz -O {sample_name}_trimmed_R2.fastq.gz"
        print(f"Running Trimmomatic for sample {sample_name}...")
        stdout, stderr = run_command(trimmomatic_command)
        print(stdout)
        print(stderr)  

        # FastQC
        fastqc_command = f"fastqc {sample_name}_trimmed_R1.fastq.gz {sample_name}_trimmed_R2.fastq.gz --outdir=./"
        print(f"Running FastQC for sample {sample_name}...")
        stdout, stderr = run_command(fastqc_command)
        print(stdout)
        print(stderr)
        
        ## Mapping with bwa-mem2
        bwa_command = f"bwa-mem2 mem {base_dir}/references/denv1.fasta {sample_name}_trimmed_R1.fastq.gz {sample_name}_trimmed_R2.fastq.gz > {sample_name}_denv1_aln.sam"
        print(f"Running bwa-mem2 denv1 for sample {sample_name}...")
        stdout, stderr = run_command(bwa_command)
        print(stdout)
        print(stderr)

        # Mapping with bwa-mem2
        bwa_command = f"bwa-mem2 mem {base_dir}/references/denv2.fasta {sample_name}_trimmed_R1.fastq.gz {sample_name}_trimmed_R2.fastq.gz > {sample_name}_denv2_aln.sam"
        print(f"Running bwa-mem2 denv2 for sample {sample_name}...")
        stdout, stderr = run_command(bwa_command)
        print(stdout)
        print(stderr)
        # Mapping with bwa-mem2
        bwa_command = f"bwa-mem2 mem {base_dir}/references/denv3.fasta {sample_name}_trimmed_R1.fastq.gz {sample_name}_trimmed_R2.fastq.gz > {sample_name}_denv3_aln.sam"
        print(f"Running bwa-mem2 denv3 for sample {sample_name}...")
        stdout, stderr = run_command(bwa_command)
        print(stdout)
        print(stderr)
        
        # Mapping with bwa-mem2
        bwa_command = f"bwa-mem2 mem {base_dir}/references/denv4.fasta {sample_name}_trimmed_R1.fastq.gz {sample_name}_trimmed_R2.fastq.gz > {sample_name}_denv4_aln.sam"
        print(f"Running bwa-mem2 denv4 for sample {sample_name}...")
        stdout, stderr = run_command(bwa_command)
        print(stdout)
        print(stderr)
        
        # Move SAM files to base directory
        sam_files_dir = os.path.join(base_dir, "sam_files")
        
        os.system(f"mv {sample_name}*.sam {sam_files_dir}")

        print("SAM files moved to the designated directory.")
        
        # Move SAM files to base directory
        sam_files_dir = os.path.join(base_dir, "sam_files")
        os.makedirs(sam_files_dir, exist_ok=True)
        os.system(f"mv {sample_name}*.sam {sam_files_dir}")
        
        # returning to base_dir
        os.chdir(base_dir)
        

    logging.info("Processing completed!")
