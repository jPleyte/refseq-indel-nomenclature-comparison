'''
Created on Jan 3, 2026

@author: pleyte
'''
import logging.config

import csv
import argparse
import hgvs.parser

import hgvs.dataproviders.uta
import hgvs.assemblymapper

from rinc.util import chromosome_map
from rinc.variant_transcript import VariantTranscript
from hgvs.transcriptmapper import TranscriptMapper
from rinc.util.log_config import LogConfig

ASSEMBLY_VERSION = "GRCh37"

class HgvsNomenclature(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
    
    def get_variants(self, variants_file: str):
        """
        Read the csv file that has the variants that will be processed. 
        """
        variants = [] 
        with open(variants_file, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                variants.append(VariantTranscript(row['chromosome'], 
                                                  int(row['position']),
                                                  row['reference'], 
                                                  row['alt'],
                                                  row['refseq_transcript']))
        
        self._logger.debug(f"Read {len(variants)} variants from {variants_file}")
        return variants

    def _get_exon(self, hdp: hgvs.dataproviders.uta.UTA_postgresql, var_c: hgvs.sequencevariant.SequenceVariant, refseq_chromosome: str) -> int:
        """
        Return exon for a transcript position
        UTA's exon start at zero.
        """
        tm = TranscriptMapper(hdp, var_c.ac, alt_ac=refseq_chromosome, alt_aln_method='splign')
        c_position = var_c.posedit.pos.start.base
        
        exon_ord = None
        note = None
        
        for exon in tm.tx_exons:
            if exon['tx_start_i'] <= c_position <= exon['tx_end_i']:
                exon_ord = exon['ord']
            
            if 'I' in exon['cigar'] or 'D' in exon['cigar']:
                note = f"The I or D is in exon {exon['ord']} with cigar {exon['cigar']}"
                
        if exon_ord is None: 
            self._logger.debug(f"Unable to find exon for {var_c}")
            
        return exon_ord, tm.strand, note
    
    def update_with_changes(self, variants: list[VariantTranscript]):
        """
        Update the variants with c. and p. changes using hgvs 
        """
        
        hp = hgvs.parser.Parser()
        hdp = hgvs.dataproviders.uta.connect()
        am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=ASSEMBLY_VERSION)

        try:
            for var in variants:
                try:
                    chromosome_map.get_refseq(var.chromosome)
                    # The g. is easy since we know it is a substitution
                    var_g = hp.parse_hgvs_variant(f"{chromosome_map.get_refseq(var.chromosome)}:g.{var.position}{var.reference}>{var.alt}")
                    var_c = am.g_to_c(var_g, var.cdna_transcript)
                    var_p = am.c_to_p(var_c)
                    
                    if var_p.posedit: 
                        var_p.posedit.uncertain = False
                    
                    var.g_dot = var_g.format().replace(var_g.ac + ':', '')
                    var.c_dot = var_c.format().replace(var_c.ac + ':', '')
                    var.protein_transcript = var_p.ac
                    var.p_dot1 = var_p.format(conf={"p_3_letter": False}).replace(var_p.ac + ':', '')
                    var.p_dot3 = var_p.format(conf={"p_3_letter": True}).replace(var_p.ac + ':', '')
                    
                    var.exon, var.strand, note = self._get_exon(hdp, var_c, chromosome_map.get_refseq(var.chromosome))
                    if note: 
                        var.notes.append(note)
                    
                    tx_info = hdp.get_tx_info(var_c.ac, chromosome_map.get_refseq(var.chromosome), 'splign')
                    var.gene = tx_info['hgnc']
                    
                except Exception as e:
                    self._logger.error(f"Error processing variant {var}: {e}")
                    raise
        finally:
            hdp.close()
            

    def write(self, out_filename, variant_transcripts: list[VariantTranscript]):
        headers = ['chromosome', 'position', 'reference', 'alt',
                   'cdna_transcript', 'hu.protein_transcript', 
                   'hu.exon', 'hu.strand', 'hu.gene', 
                   'hu.g_dot', 'hu.c_dot', 'hu.p_dot1', 'hu.p_dot3', 
                   'hu.note']
        
        with open(out_filename, 'w', newline='') as output:
            writer = csv.writer(output)
            writer.writerow(headers)
            
            for v in variant_transcripts:
                writer.writerow([v.chromosome, v.position, v.reference, v.alt, 
                                 v.cdna_transcript, v.protein_transcript, 
                                 v.exon, v.strand, v.gene,
                                 v.g_dot, v.c_dot, v.p_dot1, v.p_dot3,
                                 ";".join(v.notes)])
            
        self._logger.info(f"Wrote {len(variant_transcripts)} variant transcripts to {out_filename}")

def _parse_args():
    parser = argparse.ArgumentParser(description='Use hgvs to determine protein changes for a list of variants')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--variants", help="File with variants (csv)", required=True)
    parser.add_argument("--out", help="output file (csv)", required=True)
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)

    args = _parse_args()

    hn = HgvsNomenclature()

    # Read variants from csv 
    variants = hn.get_variants(args.variants)

    # Update variants with c. and p. changes 
    hn.update_with_changes(variants)

    hn.write(args.out, variants)
    
if __name__ == '__main__':
    main()        