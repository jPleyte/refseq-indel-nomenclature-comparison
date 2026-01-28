'''
Helper methods for workign with gffutils  
- createDb: converts a gff to the SQLite database that is required for working with gffutils.
  - usage: python -m rinc.io.gffutils_helper createDb in_gff GCF_file.gff.gz --out_db GCF_file.gff.db 
- createAccessionIndex: Create an index linking every reffseq and ccds acession to their id in the gff db and 
    includes a field with the cigar strings for transcripts that have diffeerences from the reference genome. 
  - usage: python -m rinc.io.gffutils_helper createAccessionIndex --gff_db GCF_file.gff.db --out_parquet GCF_file_accession_index.parquet

The gff can be any but i am using GCF_000001405.25_GRCh37.p13_genomic.gff.gz from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/

Known bug:
* CCDS accessions have more than one entry
    
Usage:
 
Created on Jan 22, 2026

 
@author: pleyte
'''
import argparse
import gffutils
import pandas as pd
from collections import Counter

_index_counter = Counter()

def _create_gff_db(args):
    # Create the db file 
    gffutils.create_db(args.in_gff, 
                       args.out_db,                            
                       merge_strategy='create_unique',                       
                       force=True,
                       keep_order=True, 
                       verbose = True)
    print(f"Finished creating {args.out_db}")

def _get_accession_gap_map(gff_db) -> dict:
    """
    Find every feature of type 'cDNA_match' that indicates a transcript has a ref mismatch, indicated by a cigar string
    Returns a map from refseq accession to the cigar string for the whole transcript. 
    """
    gap_lookup = {}
    
    for match in gff_db.features_of_type('cDNA_match'):
        # Target attribute looks like ['NM_001167672.3 1403 17723 +']
        target = match.attributes.get('Target', [''])[0]
        gap = match.attributes.get('Gap', [None])[0]
        
        if target and gap:
            # Extract accession and target range: Parts --> ['NM_001167672.3', '1403', '17723', '+']
            parts = target.split(' ')
            acc = target.split(' ')[0]            
            t_start = int(parts[1])
            t_end = int(parts[2])
            
            gap_lookup[acc] = {
                'gap': gap,
                'target_start': t_start,
                'target_end': t_end
            }
            
            _index_counter['gaps'] += 1
    
    print(f"  Found {len(gap_lookup)} RefSeq transcripts with gap strings")
    return gap_lookup    
    
def _get_refseq_accession_gff_indexes(gff_db, gap_map_map):
    """
    Iterate over all features of type mRNA finding refseq trancripts (NM, and NP prefix)
    and return a list of mappings from accession to feature id in the gff.
    """
    n = 0
    data = []
    for mrna_feature in gff_db.features_of_type('mRNA'):        
        names = mrna_feature.attributes.get('Name', [])
        for name in names:
            if name.startswith(('NM_', 'NP_')):
                # Look up the gap in our pre-built dictionary Returns None if the accession isn't in the gap_lookup
                gm = gap_map_map.get(name, {})
            
                data.append({'accession': name, 
                             'feature_id': mrna_feature.id, 
                             'type': 'mRNA', 
                             'gap': gm.get('gap'),
                             'target_start': gm.get('target_start'),
                             'target_end': gm.get('target_end')
                             })
        n = n + 1
        if n % 50000 == 0:
            print(f"  Processed {n} mRNA features so far")
    
    print(f"  Found {len(data)} RefSeq transcripts")
    _index_counter['refseq_transcripts'] += len(data)
    return data
            

def _get_ccds_accession_gff_indexes(gff_db, gap_map_map):
    """
    Iterate over all features of type CDS finding protein transcripts that have a CCDS accession and return a list
    lf mappings from CCDS accession to feature id in the gff.
    """
    n = 0
    data = []
    for feature in gff_db.features_of_type('CDS'):
        dbxrefs = feature.attributes.get('Dbxref', [])
        for xref in dbxrefs:
            if xref.startswith('CCDS:'):
                # Get the CCDS accession
                ccds_id = xref.split(':')[-1]
                
                if not ccds_id:
                    raise ValueError(f"Unable to get CCDS from dbxref string {xref} at CDS feature id {feature.id}")
                
                # Get the CCDS's parent 
                parent_id = feature.attributes.get('Parent', [None])[0]
                current_gap_info = {}
                
                # 2. Look up the parent mRNA object to find its NM_ name
                if parent_id:                    
                    parent_mrna = gff_db[parent_id]
                    names = parent_mrna.attributes.get('Name', [])
                    
                    # 3. Check if any of the mRNA's names have a Gap in our map
                    for name in names:                        
                        if name in gap_map_map:
                            current_gap_info = gap_map_map[name]
                            _index_counter['ccds_with_gap'] += 1
                            break
                    
                data.append({
                    'accession': ccds_id,
                    'feature_id': feature.id,
                    'type': 'CDS',
                    'gap': current_gap_info.get('gap'),
                    'target_start': current_gap_info.get('target_start'),
                    'target_end': current_gap_info.get('target_end')
                })

        n = n + 1
        if n % 100000 == 0:
            print(f"  Processed {n} CDS features so far")
    
    print(f"  Found {len(data)} CCDS transcripts")
    _index_counter['ccds_transcripts'] += len(data)
    return data
    
def _create_accession_index(args):
    gff_db = gffutils.FeatureDB(args.gff_db)
    
    print("Pre-indexing Gap attributes from cDNA_match features for RefSeq transcripts...")
    gap_map_map = _get_accession_gap_map(gff_db)
    
    accession_gff_indexes = []
    print("Extracting mRNA RefSeq accessions...")
    accession_gff_indexes.extend(_get_refseq_accession_gff_indexes(gff_db, gap_map_map)) 
    
    print("Extracting CDS CCDS accessions...")
    accession_gff_indexes.extend(_get_ccds_accession_gff_indexes(gff_db, gap_map_map))
    
    print("Creating index dataframe...")
    df = (
        pd.DataFrame(accession_gff_indexes)
        .drop_duplicates(subset=['accession', 'gap', 'target_start', 'target_end'])
        .reindex(columns=['accession', 'feature_id', 'type', 'gap', 'target_start', 'target_end'])
        .astype({
            'target_start': 'Int64',
            'target_end': 'Int64',
            'gap': 'string'
        })
        .set_index('accession')
        .sort_index()
    )
    
    print(f"Saving dataframe {args.out_parquet}...")
    print(f"{_index_counter}")
    
    # Write out an binary, indexed version of the dataframe
    df.to_parquet(args.out_parquet)
    
    # Optionally, write out a csv to make it easier for human review
    if args.out_csv:        
        print(f"Saving csv dataframe: {args.out_csv}")
        df.to_csv(args.out_csv)
        
def get_feature_id(db, accession):
    if accession.startswith(('NM', 'NP')):
        for mrna in db.features_of_type('mRNA'):
            if 'refseq_accession' in mrna.attributes and mrna.attributes['refseq_accession'][0] == accession:
                # acc = transcript.attributes['refseq_accession'][0]
                #refseq_map[acc] = transcript.id
                return mrna.id

    elif accession.startswith("CCDS"):
        for cds in db.features_of_type('CDS'):
            if 'ccds_accession' in cds.attributes and cds.attributes['ccds_accession'][0] == accession:
                #acc = cds.attributes['ccds_accession'][0]
                #ccds_map[acc] = cds.id
                return cds.id
    
    else:
        return ValueError(f"Unknown accession prefix: {accession}")

def get_feature_by_accession(gff_db, accession_index_df, acc):
    try:
        # Index lookups in a sorted pandas index are O(log n) - extremely fast
        feat_id = accession_index_df.loc[acc, 'feature_id']
        return gff_db[feat_id]    
    except KeyError:
        return None
    
def main():
    parser = argparse.ArgumentParser(description='Convert a GFF to a db file usable by gffutils')

    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    gff_to_db_parser = subparsers.add_parser("createDb", help="Convert gff to db")
    gff_to_db_parser.set_defaults(func=lambda args: _create_gff_db(args))
    gff_to_db_parser.add_argument('--in_gff', help="Input file (gff)", required=True)    
    gff_to_db_parser.add_argument("--out_db", help="output file (db)", required=True)
    
    gff_to_db_parser = subparsers.add_parser("createAccessionIndex", help="Create an index of Feature ids to the records with refseq and ccds accessions")
    gff_to_db_parser.set_defaults(func=lambda args: _create_accession_index(args))
    gff_to_db_parser.add_argument("--gff_db", help="gffutils gff database (db)", required=True)
    gff_to_db_parser.add_argument("--out_parquet", help="Index file (parquet)", required=True)
    gff_to_db_parser.add_argument("--out_csv", help="Index file (csv)", required=False)
    
    args = parser.parse_args()
    
    if hasattr(args, 'func'):
        args.func(args)
        
    print("Done.")    

if __name__ == '__main__':
    main()
