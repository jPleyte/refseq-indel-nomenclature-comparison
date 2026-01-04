'''
Created on Jan 3, 2026

@author: pleyte
'''

import csv
import logging
import argparse
import hgvs.parser
import hgvs.variantmapper
import hgvs.dataproviders.uta
import hgvs.assemblymapper

from rinc.util import chromosome_map
from variant_transcript import VariantTranscript

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
    
    def get_variants(self, variants_file):
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
        return variants

    def update_with_changes(self, variants: list[VariantTranscript]):
        """
        Update the variants with c. and p. changes using hgvs 
        """

        hp = hgvs.parser.Parser()
        hdp = hgvs.dataproviders.uta.connect()
        am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=ASSEMBLY_VERSION)

        for var in variants:
            try:
                chromosome_map.get_refseq(var.chromosome)
                var_g = hp.parse_hgvs_variant(f"{chromosome_map.get_refseq(var.chromosome)}:g.{var.position}{var.reference}>{var.alt}")
                var_c = am.g_to_c(var_g, var.cdna_transcript)
                
                # Map to protein variant
                p_variant = am.c_to_p(var_c)
                
                var.p_dot = str(p_variant)
                var.c_dot = str(var_c)
            except Exception as e:
                self._logger.error(f"Error processing variant {var}: {e}")
                raise


def _parse_args():
    parser = argparse.ArgumentParser(description='Use hgvs to determine protein changes for a list of variants')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--variants", help="File with variants (csv)", required=True)
    parser.add_argument("--out", help="output file (csv)", required=True)
    args = parser.parse_args()
    return args    

def main():
    args = _parse_args()
    
    hn = HgvsNomenclature()
    
    # Read variants from csv 
    variants = hn.get_variants(args.variants)

    # Update variants with c. and p. changes 
    hn.update_with_changes(variants)
    
    hn.write(args.out, variants)
    
if __name__ == '__main__':
    main()        