'''
Takes any a list of files, extracts the cdna transcripts and write a distinct list to file

Created on Feb 7, 2026

@author: pleyte
'''
from rinc.util.log_config import LogConfig
import argparse
import logging.config
import pandas as pd

def get_transcripts(filenames: list):
    unique_transcripts = set()
    for f in filenames:
        df = pd.read_csv(f, usecols=['cdna_transcript'], dtype={'cdna_transcript': str})
        valid_values = df['cdna_transcript'].dropna().unique()
        print(f"Read {valid_values.shape[0]} unique transcripts from {f}")
        unique_transcripts.update(valid_values)
    
    return unique_transcripts

def write(output_filename: str, transcripts: list):
    """
    Write transcripts to file 
    """
    if transcripts:
        with open(output_filename, 'w') as f:
            for transcript in sorted(transcripts):
                f.write(f"{transcript}\n")
        print(f"Wrote {len(transcripts)} unique transcripts to {output_filename}")
    else:
        logging.warning("No transcripts found to write.")
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Extract the cdna transcripts from any number of files')

    parser.add_argument('--nomenclature', action='append', help="Filename with cdna_transcript field (csv)", required=True)    
    parser.add_argument("--out", help="output file (csv)", dest="output", required=True)
    parser.add_argument("--version", action="version", version="0.0.1")
    return parser.parse_args()

def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    args = _parse_args()
    transcripts = get_transcripts(args.nomenclature)
    write(args.output, transcripts)


if __name__ == '__main__':
    main()