# Single-End RNA-Seq Quantification Pipeline with Salmon

This Bash script automates the process of quantifying single-end RNA-seq data using Salmon. It downloads the human reference transcriptome, builds a Salmon index, and quantifies multiple samples processed by fastp.

## Features

- Downloads the human reference transcriptome (Ensembl GRCh38)
- Builds a Salmon index
- Processes multiple single-end RNA-seq samples in batch
- Uses Salmon for transcript quantification

## Prerequisites

- Bash shell
- Salmon (accessible in the system PATH)
- wget
- gunzip
- Internet connection (for downloading the reference transcriptome)

## Usage

```bash
./script_name.sh <FASTP_FILES_DIR> <SALMON_QUANT_DIR>
<FASTP_FILES_DIR>: Directory containing fastp-processed single-end FASTQ files
<SALMON_QUANT_DIR>: Output directory for Salmon quantification results
Input Directory Structure
The input directory (FASTP_FILES_DIR) should have this structure:

FASTP_FILES_DIR/
├── sample1/
│   └── sample1_cleaned.fastq
├── sample2/
│   └── sample2_cleaned.fastq
└── ...
Output Directory Structure
The output directory (SALMON_QUANT_DIR) will have this structure:

SALMON_QUANT_DIR/
├── sample1_quant/
│   ├── quant.sf
│   └── ...
├── sample2_quant/
│   ├── quant.sf
│   └── ...
└── ...

How it works
The script takes two command-line arguments: the input directory and the output directory.
It downloads the human reference transcriptome if not already present.
It uncompresses the downloaded reference transcriptome.
It builds a Salmon index if not already present.
For each sample in the input directory:
It extracts the sample name (first 12 characters).
It finds the corresponding FASTQ file.
It creates an output directory for the sample.
It runs Salmon quantification on the sample.
The script continues until all samples are processed.
Reference Transcriptome
The script uses the Homo sapiens GRCh38 cDNA reference from Ensembl (release 109). If a different reference is needed, update the REFERENCE_URL variable in the script.

Salmon Parameters
The script uses the following Salmon parameters:

-l A: Automatically determine the library type
-p 8: Use 8 threads for quantification
Adjust these parameters in the script if needed.

Sample Naming Convention
The script assumes that the sample name is contained within the first 12 characters of the directory name. Adjust the cut -c1-12 command if your naming convention differs.

Error Handling
The script checks for the existence of input sample directories and files.
It creates output directories as needed.
Note
Ensure that your input FASTQ files are single-end reads and follow the naming convention that matches the sample name (first 12 characters) followed by any additional characters and ending with .fastq.
