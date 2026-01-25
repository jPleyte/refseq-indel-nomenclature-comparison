'''
Created on Jan 22, 2026

@author: pleyte
'''
import argparse
import logging.config
import pandas as pd
from rinc.util.log_config import LogConfig
from collections import Counter
import gffutils

class WriteExonDetail(object):
    '''
    classdocs
    '''
    def __init__(self, gff_db_file, accession_index):
        self._logger = logging.getLogger(__name__)        
        self._transcripts = set()        

        self._transcript_lookup_counts = Counter()
        
        self._gff_db = gffutils.FeatureDB(gff_db_file)
        self._accession_index_df = pd.read_parquet(accession_index)
    
    def _get_accession_feature_id_map(self):
        accession_feature_id_map = {} 
        for mrna in self._gff_db.features_of_type('mRNA'):
            if 'refseq_accession' in mrna.attributes:
                acc = mrna.attributes['refseq_accession'][0]
                accession_feature_id_map[acc] = mrna.id
                
        for cds in self._gff_db.features_of_type('CDS'):
            if 'ccds_accession' in cds.attributes:
                acc = cds.attributes['ccds_accession'][0]
                accession_feature_id_map[acc] = cds.id
        
        return accession_feature_id_map
        
    def add_transcripts(self, file):
        """
        Add a datasource 
        """
        df = pd.read_csv(file, dtype=str)
        transcripts = set(df['cdna_transcript'])
        self._transcripts.update(transcripts)
        
        self._logger.info(f"Read {df.shape[0]} rows and {len(transcripts)} transcripts from {file}")
    
    def _get_exon_details(self, accession, feature):
        """
        """        
        if accession.startswith('NM'):
            transcript_feature = feature
        elif accession.startswith('CCDS'):
            parent_ids = feature.attributes.get('Parent', [])
            if not parent_ids:
                return []
            transcript_feature = self._gff_db[parent_ids[0]]
        else:
            raise ValueError(f"Unknow accession type: {accession}")

        return self._get_exon_details_mrna(accession, transcript_feature)
            
    def _get_exon_details_mrna(self, accession, feature):
        """
        Return the list of exons associated with a transcritp accession  
        """
        exons = list(self._gff_db.children(feature, featuretype='exon', order_by='start'))
        
        strand = feature.strand
        if strand == '-':
            exons.reverse()
        
        exon_details = []
        exon_number = 0
        for x in exons:
            exon_number += 1
            gap_attribute = x.attributes.get('Gap', [None])[0]
            if gap_attribute:
                self._transcript_lookup_counts['gaps'] += 1
                
            exon_details.append({'refseq_chrom': x.chrom,
                                 'accession': accession,
                                 'exon_number': exon_number,
                                 'start': x.start,
                                 'stop': x.stop,
                                 'strand': strand,
                                 'gene': x.attributes['gene'][0],
                                 'length': len(x)
                                 })

        return exon_details
            
    
    def _get_feature(self, accession):
        """
        """
        feat_id = self._accession_index_df['feature_id'].get(accession)
        if feat_id:
            return self._gff_db[feat_id]
        
        return None
    
    def get_tx_exon_details(self) -> list[dict]:
        """
        Iterate over the set of transcripts and look up exon coordinates and cigar strings
        """
        transcript_exon_details = []
        
        n = 0
        
        self._transcripts.clear()        
        self._transcripts.add("NM_001167672.3")
        self._transcripts.add("NM_015068.3")
        self._transcripts.add("NM_013386.4")
        
        for x in self._transcripts:
            feature = self._get_feature(x)
            if not feature:
                self._transcript_lookup_counts["tx_not_found"] += 1
                continue
        
            exon_details = self._get_exon_details(x, feature)
            self._transcript_lookup_counts["tx_found"] += 1
            transcript_exon_details.extend(exon_details)
            self._transcript_lookup_counts["tx_exons"] += len(exon_details)
            if not exon_details:
                self._transcript_lookup_counts["jdebug_no_exons"] += 1
            
            n += 1
            if n % 10000 == 0:
                self._logger.info(f"Processed {n}/{len(self._transcripts)}: {self._transcript_lookup_counts}")            
        
        return transcript_exon_details
    
    def write(self, out_file, tx_exon_details):
        """
        Write the exons to file
        """        
        df = pd.DataFrame(tx_exon_details)
        df.to_csv(out_file)
        self._logger.info(f"Wrote {df.shape[0]} exon rows for {len(self._transcripts)} transcritps to {out_file}")
        self._logger.info(f"Summary: {self._transcript_lookup_counts}")

def _parse_args():
    parser = argparse.ArgumentParser(description='Read multiple csvs and extract trancripts, write out exon coords for each transcript')

    parser.add_argument('--gff_db', help="The SQLite gff db by gffutils (db)", required=True)
    parser.add_argument('--accession_index', help="Transcript accession index for looking up exons in the gff db (parquet)", required=True)
    parser.add_argument('--transcripts', action='append', help="File containing trancripts (csv)", required=True)    
    parser.add_argument("--out", help="output file (csv)", dest="output", required=True)

    parser.add_argument("--version", action="version", version="0.0.1")

    return parser.parse_args()

    
def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    wed = WriteExonDetail(args.gff_db, args.accession_index)
    
    if args.transcripts:
        for t_file in args.transcripts:
            wed.add_transcripts(t_file)
    
    tx_exon_details = wed.get_tx_exon_details()
    wed.write(args.output, tx_exon_details)
    

if __name__ == '__main__':
    main()        