#!/bin/bash

# Input directory name for fastp processed files
FASTP_FILES_DIR=$1

# Full path of the output quant directory
SALMON_QUANT_DIR=$2

# Check if the output directory was provided
if [ -z "${SALMON_QUANT_DIR}" ]; then
    echo "No output directory for Salmon quantification files provided."
    exit 1
fi

# Step 0: Create directories if they don't exist
mkdir -p "${SALMON_QUANT_DIR}"

# Step 1: Download reference transcriptome if not already present
REFERENCE_FA="Homo_sapiens.GRCh38.cdna.all.fa.gz"
REFERENCE_URL="https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/${REFERENCE_FA}"

if [ ! -f "${REFERENCE_FA}" ]; then
    echo "Downloading reference transcriptome..."
    wget -c "${REFERENCE_URL}" -O "${REFERENCE_FA}" &> download.log
fi

# Step 2: Build the Salmon index if not already present
INDEX_DIR="Homo_sapiens.GRCh38.cdna.all.salmon_index"

if [ ! -d "${INDEX_DIR}" ]; then
    echo "Building Salmon index..."
    salmon index -t "${REFERENCE_FA}" -i "${INDEX_DIR}"
fi

# Step 3: Quantify each sample with Salmon
for sample_dir in "${FASTP_FILES_DIR}"/*; do
    if [ -d "${sample_dir}" ]; then
        sample_name=$(basename "${sample_dir}")
        echo "Processing ${sample_name}..."

        # Path to read files
        read1="${sample_dir}/${sample_name}_cleaned_1.fastq"
        read2="${sample_dir}/${sample_name}_cleaned_2.fastq"

        # Output directory for Salmon quant
        quant_dir="${SALMON_QUANT_DIR}/${sample_name}_quant"

        # Make sure the output directory exists
        mkdir -p "${quant_dir}"

        # Run Salmon quant
        salmon quant -i "${INDEX_DIR}" -l A \
                     -1 "${read1}" \
                     -2 "${read2}" \
                     -p 8 \
                     -o "${quant_dir}"
    fi
done

echo "Salmon quantification complete."
