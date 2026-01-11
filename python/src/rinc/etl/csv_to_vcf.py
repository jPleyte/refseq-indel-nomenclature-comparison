'''
Convert csv to tab delimited SnpEff input file
Created on Jan 9, 2026

@author: pleyte
'''
import argparse
import csv
import logging.config
import pysam 

from rinc.util.log_config import LogConfig
from rinc.util import chromosome_map


VERSION = '0.0.1'

class CsvToVcf(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
    
    def read(self, in_file_csv):
        """
        Read the csv and return list of variants 
        """         
        variants = {}
        
        with open(in_file_csv, mode='r', newline='', encoding='utf-8') as csvfile:
            for row in csv.DictReader(csvfile):
                chromosome = row['chromosome']
                position = row['position']
                reference = row['reference']
                alt = row['alt']
                key = f"{chromosome}-{position}-{reference}-{alt}"
                if key not in variants:
                    variants[key] = row
        
        return variants.values()
    
    def write(self, variants: list, out_file_vcf):
        """
        Write list of variants to VCF file 
        """
        header = pysam.VariantHeader()
        header.add_line("##fileformat=VCFv4.2")        
        header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
        
        for x in chromosome_map.refseq_to_ncbi.values():            
            header.contigs.add(x)

        with pysam.VariantFile(out_file_vcf, "w", header=header) as vcf_out:
            for x in variants:
                # pysam expects zero based position, and then adds one
                start_pos = int(x['position']) - 1
                rec = vcf_out.new_record(
                    contig=x['chromosome'].replace('chr', ''),
                    start=start_pos, 
                    stop=start_pos + len(x['reference']),
                    alleles=(x['reference'], x['alt']),
                    id=".",
                    qual=None,
                    filter="PASS")
            
                vcf_out.write(rec)

def _parse_args():
    parser = argparse.ArgumentParser(description='Read variants from a csv formatted file and write out a VCF')

    parser.add_argument('--in',
                dest="input",
                help='csv list of variants',
                required=True)

    parser.add_argument('--out',
                dest="output",
                help='Vcf out file',
                required=True)

    parser.add_argument('--version', action='version', version=VERSION)

    return parser.parse_args()


def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    logger = logging.getLogger("rinc.etl.csv_to_vcf")

    args = _parse_args()

    c2v = CsvToVcf()
    variants = c2v.read(args.input)

    c2v.write(variants, args.output)

    logger.info(f"Wrote {len(variants)} variants to {args.output}")


if __name__ == '__main__':
    main()
