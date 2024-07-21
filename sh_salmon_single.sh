#!/bin/bash

# Set the first argument as the directory containing the cleaned single-end read files
FASTP_FILES_DIR="$1"
# Set the second argument as the directory where you want to store Salmon quantifications
SALMON_QUANT_DIR="$2"

# Set the URL for the reference fasta
REFERENCE_URL="https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
REFERENCE_FA="Homo_sapiens.GRCh38.cdna.all.fa.gz"

# Step 1: Download the reference fasta file if it doesn't exist
if [ ! -f "${REFERENCE_FA}" ]; then
    echo "Downloading reference transcriptome..."
    wget -c "${REFERENCE_URL}" -O "${REFERENCE_FA}" &> download.log
    # Uncompress the reference fasta file
    gunzip "${REFERENCE_FA}"
    # Update REFERENCE_FA to the uncompressed filename
    REFERENCE_FA="${REFERENCE_FA%.gz}"
fi

# Step 2: Build the Salmon index if not already present
INDEX_DIR="Homo_sapiens.GRCh38.cdna.all.salmon_index"

if [ ! -d "${INDEX_DIR}" ]; then
    echo "Building Salmon index..."
    salmon index -t "${REFERENCE_FA}" -i "${INDEX_DIR}"
fi

# Step 3: Quantify each sample with Salmon for single-end reads
for sample_dir in "${FASTP_FILES_DIR}"/*; do
    if [ -d "${sample_dir}" ]; then
        # Extracts the base name and then cuts at the first '.' or '-', using sed
        sample_name=$(basename "${sample_dir}" | cut -c1-12)

        # Find the first fastq file that matches the sample_name pattern
        read1=$(find "${sample_dir}" -maxdepth 1 -name "${sample_name}*.fastq" | head -n 1)

        # Output directory for Salmon quant
        quant_dir="${SALMON_QUANT_DIR}/${sample_name}_quant"

        # Make sure the output directory exists
        mkdir -p "${quant_dir}"

        # Run Salmon quant for single-end read
        salmon quant -i "${INDEX_DIR}" -l A \
                     -r "${read1}" \
                     -p 8 \
                     -o "${quant_dir}"
    fi
done

echo "Salmon quantification complete."
