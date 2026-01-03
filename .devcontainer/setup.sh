#!/bin/bash

# Update and install system dependencies
sudo apt-get update
sudo apt-get install -y samtools perl wget tabix libpq-dev

/opt/conda/bin/pip install psycopg2-binary hgvs pysam biopython

mkdir -p data

if [ ! -f data/Homo_sapiens_assembly19.fasta ]; then
    echo "Downloading hg19... this may take a moment."
    wget -P data/ https://storage.googleapis.com/gatk-legacy-bundles/b37/Homo_sapiens_assembly19.fasta
    wget -P data/ https://storage.googleapis.com/gatk-legacy-bundles/b37/Homo_sapiens_assembly19.fasta.fai
    wget -P data/ https://storage.googleapis.com/gatk-legacy-bundles/b37/Homo_sapiens_assembly19.dict
fi

# Update Nextflow
nextflow self-update
nextflow -version

# Annovar
chmod +x annovar/*.pl
echo "Before you can use annovar you must download the annovar database files"
