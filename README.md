# refseq-indel-nomenclature-comparison

This project uses CIGAR strings to identify insertion or deletion alignment differences between hg19 and RefSeq transcript sequences. By simulating variants immediately downstream of these changes, the pipeline generates cDNA and protein nomenclature using both ANNOVAR and the hgvs/UTA Python package. Finally, the outputs are compared to determine if the two methods yield identical nomenclature for the same variant.

## Running the workflow

## Running the workflow in Github Codespaces

This repsoitory has a ``.devcontainer`` directory that makes it able to be deployed as a GitHub Codespace. When the Codespace is first launched the ``.devcontainer/setup.sh`` script will install necessary dependencies including the hg19 fasta and annovar databases. 

To run the workflow change to the nextflow directory and launch the worfklow:
```bash
nextflow run .
```

By default the workflow queries the UTA database and creates a list of variants. But if you have a list of variants you want to process you can place the the variants in ``/tmp/curated_gaps_and_variants.csv`` and then launch the wokflow with the ``-stub`` parameter. This causes the workflow to skip the first step that generates variants, and to use your list instead.

```bash
nextflow run . -stub
```

## Running the workflow locally 

To run the workflow on your local environment edit the ``nextflow/nextflow.config`` file and make the following changes to the ``mac`` profile:
* Set the ``params.fasta`` parameter to the full path to the hg19 fasta
* Set the ``params.uta_schema`` parameter to the schema in your UTA database
* Set the env ``ANNOVAR_HOME`` to directory where you installed annovar. 

Launch the workflow specifying the ``mac`` profile

```bash
run main.nf -profile mac
```
