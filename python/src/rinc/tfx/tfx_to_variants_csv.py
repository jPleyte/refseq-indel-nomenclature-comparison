'''
Read transcirpt nomenclature from a tfx json and write out just the variants to csv

Created on Jan 15, 2026

@author: pleyte
'''
from rinc.util.log_config import LogConfig
import argparse
import csv
import logging.config
import json
from rinc.variant_transcript import VariantTranscript
from rinc.io import variant_helper
from rinc.util import chromosome_map

class TfxToVariantsCsv(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
    
    def get_variant_transcripts(self, tfx_json_file: str, no_ccds=False) -> list[VariantTranscript]:
        """
        """
        with open(tfx_json_file, 'r') as f:
            data = json.load(f)
        
        self._logger.info(f"Read {len(data)} tfx from {tfx_json_file}")
        
        variant_transcripts = []
        
        for x in data['transcriptEffects']:
            if not x['cdnaTranscript']:
                self._logger.info(f"Skipping variant that doesn't have transcripts: {x['chromosome']}-{x['position']}-{x['reference']}-{x['alt']}-{x['cdnaTranscript']}")
                continue
            if no_ccds and x['cdnaTranscript'].startswith('CCDS'):
                self._logger.info(f"Skipping CCDS transcript: {x['chromosome']}-{x['position']}-{x['reference']}-{x['alt']}-{x['cdnaTranscript']}")
                continue
            v = VariantTranscript(x['chromosome'], x['position'], x['reference'], x['alt'], x['cdnaTranscript'])
            g_dot_part = x['sequenceVariant']
            refseq_chromosome = chromosome_map.get_refseq(v.chromosome)
            v.g_dot = f"{refseq_chromosome}:{g_dot_part}"
            
            variant_transcripts.append(v)
            
        return variant_transcripts

    def write(self, output_file: str, variants: list[VariantTranscript]):
        """
        Write variants to file using variant_helper module
        """
        variant_helper.write_variants(output_file, variants)
        self._logger.info(f"Wrote {len(variants)} variants to {output_file}")
        
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Read transcirpt nomenclature from a tfx json and write out just the variants to csv')

    parser.add_argument('-i', '--in',
                dest="input",
                help='Transcript Effects/tfx (json)',
                required=True)

    parser.add_argument('-o', '--out',
                dest="output",
                help='Variants (csv)',
                required=True)

    parser.add_argument("--no_ccds", action='store_true', help="Filter out all CCDS transcripts", required=False)
    parser.add_argument('--version', action='version', version='0.0.1')

    return parser.parse_args()

def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    logger = logging.getLogger("rinc.etl.tfx_to_variants_csv")

    args = _parse_args()

    t2v = TfxToVariantsCsv()
    variants = t2v.get_variant_transcripts(args.input, args.no_ccds)
    t2v.write(args.output, variants)

    logger.info(f"Wrote {len(variants)} variants to {args.output}")


if __name__ == '__main__':
    main()
