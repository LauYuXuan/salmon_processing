# Paired-end read Salmon Quantification Pipeline

This Bash script automates the process of quantifying RNA-seq data using Salmon. It downloads the reference transcriptome, builds a Salmon index, and quantifies multiple samples processed by fastp.

## Features

- Downloads the human reference transcriptome (Ensembl GRCh38)
- Builds a Salmon index
- Processes multiple samples in batch
- Uses Salmon for transcript quantification
- Handles paired-end RNA-seq data

## Prerequisites

- Bash shell
- Salmon (accessible in the system PATH)
- wget
- Internet connection (for downloading the reference transcriptome)

## Usage

```bash
./script_name.sh <FASTP_FILES_DIR> <SALMON_QUANT_DIR>
<FASTP_FILES_DIR>: Directory containing fastp-processed FASTQ files
<SALMON_QUANT_DIR>: Output directory for Salmon quantification results
Input Directory Structure
The input directory (FASTP_FILES_DIR) should have this structure:

FASTP_FILES_DIR/
├── sample1/
│   ├── sample1_cleaned_1.fastq
│   └── sample1_cleaned_2.fastq
├── sample2/
│   ├── sample2_cleaned_1.fastq
│   └── sample2_cleaned_2.fastq
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
The script checks if the output directory is provided.
It creates the output directory if it doesn't exist.
It downloads the human reference transcriptome if not already present.
It builds a Salmon index if not already present.
For each sample in the input directory:
It identifies the paired-end FASTQ files.
It creates an output directory for the sample.
It runs Salmon quantification on the sample.
The script continues until all samples are processed.
Reference Transcriptome
The script uses the Homo sapiens GRCh38 cDNA reference from Ensembl (release 109). If a different reference is needed, update the REFERENCE_FA and REFERENCE_URL variables in the script.

Salmon Parameters
The script uses the following Salmon parameters:

-l A: Automatically determine the library type
-p 8: Use 8 threads for quantification
Adjust these parameters in the script if needed.

Error Handling
The script checks if the output directory is provided.
It verifies the existence of input sample directories and files.
Note
Ensure that your input FASTQ files follow the naming convention <sample_name>_cleaned_1.fastq and <sample_name>_cleaned_2.fastq for paired-end data.
