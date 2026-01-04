# refseq-indel-nomenclature-comparison

This project uses CIGAR strings to identify insertion or deletion alignment differences between hg19 and RefSeq transcript sequences. By simulating variants immediately downstream of these changes, the pipeline generates cDNA and protein nomenclature using both ANNOVAR and the hgvs/UTA Python package. Finally, the outputs are compared to determine if the two methods yield identical nomenclature for the same variant.

