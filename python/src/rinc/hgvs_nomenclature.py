'''
Created on Jan 3, 2026

@author: pleyte
'''

import logging
import argparse

class HgvsNomenclature(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
    
    def get_changes(self, in_file):
        """
        Read the csv file that has the variants that will be processed. 
        """

def _parse_args():
    parser = argparse.ArgumentParser(description='Use hgvs to determine protein changes for a list of variants')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--variants_file", help="File with variants (csv)", required=True)
    parser.add_argument("--out", help="output file (csv)", required=True)
    args = parser.parse_args()
    return args    

def main():
    args = _parse_args()
    
    hn = HgvsNomenclature()
    
    variant_transcripts = hn.get_changes(args.variants_file)
    
    hn.write(args.out, variant_transcripts)
    
if __name__ == '__main__':
    main()        