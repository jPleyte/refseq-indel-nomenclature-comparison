#!/bin/bash

# Update and install system dependencies
sudo apt-get update
sudo apt-get install -y samtools perl wget tabix libpq-dev awscli

/opt/conda/bin/pip install psycopg2-binary hgvs pysam biopython

mkdir -p data

if [ ! -f data/Homo_sapiens_assembly19.fasta ]; then
    echo "Downloading hg19... this may take a moment."
    
    aws s3 cp s3://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta data/ --no-sign-request
    aws s3 cp s3://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai data/ --no-sign-request
    aws s3 cp s3://broad-references/hg19/v0/Homo_sapiens_assembly19.dict data/ --no-sign-request    
fi

# Update Nextflow
nextflow self-update
nextflow -version

# Annovar
chmod +x annovar/*.pl
echo "Before you can use annovar you must download the annovar database files"
perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGeneWithVer annovar/humandb/
# perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGeneWithVerMRNA annovar/humandb/

# Customise the terminal command prompt
printf "export PS1='\\[\\e[3;36m\\]\${PWD#/workspaces/} ->\\[\\e[0m\\] '\n" >> $HOME/.bashrc
export PS1='\[\e[3;36m\]${PWD#/workspaces/} ->\[\e[0m\] '

NEXTFLOW_DIR="/workspaces/$(basename $CONTAINER_WORKSPACE_FOLDER)/nextflow"

cat << 'EOF' >> $HOME/.bashrc
cleanup() {
    cd $NEXTFLOW_DIR
    nextflow clean -f 
    rm -f $NEXTFLOW_DIR/.nextflow.log*
    rm -f $NEXTFLOW_DIR/results/*    
    echo "Removed nextflow logs and results"
}
EOF