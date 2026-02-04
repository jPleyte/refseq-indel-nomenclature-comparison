'''
Generate variant-transcript batch files for submission to Mutalyzer.
We don't want to submit all variants to mutalzyer.  We just want the ones where 
the two input dataframes have the same c_dot values. 
 
Created on Jan 30, 2026

@author: pleyte
'''
import argparse
import logging.config
from rinc.util.log_config import LogConfig
import requests
import re
from rinc.util.pdot import PDot
from collections import Counter
from rinc.variant_transcript import VariantTranscript
from rinc.io import variant_helper
from rinc.util import chromosome_map

class MutalyzerBatch(object):
    '''
    Create batch files for submission to Mutalyzer 
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._p_dot_converter = PDot()
        self._api_response_counter = Counter()
        self._protein_response_counter = Counter()
        self._mutalizer_normalize_url = "https://mutalyzer.nl/api/normalize/{value}"
    
    def get_variant_transcripts(self, variants_file: str) -> list[VariantTranscript]:
        """        
        Receives a list of variants, and uses Mutalyzer's normalizer API to get transcripts and nomenclature.
        Returns a list of of variant-transcripts      
        """
        variants = variant_helper.get_variants(variants_file)
        
        n = 0
        variant_transcripts = []
        for v in variants:
            n = n + 1
            if n % 100 == 0:
                self._logger.info(f"Processed {n}/{len(variants)} variants")

            variant_transcripts.extend(self._get_variant_transcripts(v))

        self._logger.info(f"Read {len(variant_transcripts)} variants from {variants_file}")
        return variant_transcripts
    
    def _get_variant_transcripts(self, v: VariantTranscript):
        """
        Take a single variant and 
        1. Convert to g. and request all c. and transcripts for that g. from mutalyzer using their API 
        2. For every c. and transcript, request p. from mutalyzer
        
        Returns a list of all the transcripts for a single variant 
        """
        variant_transcripts = []

        simple_g_dot = self._get_simple_g_dot(v)
        transcript_cdots = self._get_transcripts(simple_g_dot)
        
        for t in transcript_cdots:
            cp_nomenclature = self._get_cdna_and_protein_nomenclature(t['cdna_transcript'], t['c_dot'])
            variant_transcripts.append(self._get_variant_transcript(v, t['g_dot'], cp_nomenclature))
        
        return variant_transcripts
            
            
    def _get_variant_transcript(self, v: VariantTranscript, g_dot: str, cdna_protein_values: dict):
        """
        Set the VariantTranscript fields  
        """
        v = VariantTranscript(v.chromosome, v.position, v.reference, v.alt, cdna_protein_values['cdna_transcript']) 
        v.g_dot = g_dot
        v.gene = cdna_protein_values['gene']
        v.c_dot = cdna_protein_values['c_dot']
        v.p_dot1 = cdna_protein_values['p_dot1']
        v.p_dot3 = cdna_protein_values['p_dot3']
        v.protein_transcript = cdna_protein_values['protein_transcript']
        v.additional_fields['code'] = cdna_protein_values['code']
        return v
        
    def _get_simple_g_dot(self, v: VariantTranscript) -> str:
        """
        Return a simple g dot string     
        """
        chromosome = chromosome_map.get_refseq(v.chromosome)
        return f"{chromosome}:g.{v.position}{v.reference}>{v.alt}"
    
    def _get_transcripts(self, simple_g_dot) -> list:
        """
        Use request a list of refseq transcripts from Mutalyzer API 
        """
        response, _ = self._get_api_response(simple_g_dot)
         
        transcript_cdots = []
        if response:
            corrected_g_dot = response['corrected_description']
            for c in response['equivalent_descriptions']['c']:
                refseq_transcript, c_dot = self._get_cdna_values(c['description'])
                if refseq_transcript and c_dot:
                    # jdebug self._cdna_response_counter['c_dot.has_value'] += 1
                    transcript_cdots.append({
                        'g_dot': corrected_g_dot,
                        'cdna_transcript': refseq_transcript,
                        'c_dot': c_dot
                    })                    
                else:                    
                    # jdebug if this happens, then consider logging the code 
                    #self._cdna_response_counter['c_dot.no_value'] += 1
                    raise ValueError(f"No c. for {simple_g_dot}: why not?")

        return transcript_cdots
                
    def _get_cdna_values(self, c_description: str):
        """
        Parse cDNA transcript and c. out of the c description string 
        """        
        regex_pattern = r"(?P<chromosome>NC_\d+\.\d+)\((?P<transcript>NM_\d+\.\d+)\):(?P<c_dot>c\.[^(\s]+)"
        match = re.search(regex_pattern, c_description)
        
        if match:
            d = match.groupdict()            
            return d['transcript'], d['c_dot']
        else:
            self._logger.warning(f"Unknown cDNA format: {c_description}")        
            return None, None 
        
    def _get_api_response(self, normalize_parameter):
        """
        Make API request toMutalyzer 
        """
        encoded_parameter = requests.utils.quote(normalize_parameter)
        url = self._mutalizer_normalize_url.format(value=encoded_parameter)
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            self._api_response_counter[f'{response.status_code}.success'] += 1
            return response.json(), None
        else:
            error_code = self._get_response_error_code(response.json())
            self._api_response_counter[f'{response.status_code}.{error_code}'] += 1
            self._logger.info(f"Unable to get normalizer result for {normalize_parameter}: {error_code}")
            return None, error_code

    def _get_response_error_code(self, mutalyzer_response: dict) -> str:
        """
        Return the error code from a mutalyzer reponse        
        """
        try:
            transcript_request = mutalyzer_response['custom']['corrected_description']
            errors = mutalyzer_response['custom']['errors']
            if len(errors) != 1:
                # jdebug: If this happens, then just concatenate error codes  
                raise ValueError(f"Received {len(errors)} errors when only one was expected for {transcript_request}")
            error_code = None
            for e in errors:
                error_code = e['code']
                self._logger.info(f"Received error code {error_code} for {transcript_request}: {e['details']}")
                return error_code
        except (KeyError, IndexError):
            self._logger.warning(f"Unable to process error response: {mutalyzer_response}")
            raise

    def _get_cdna_and_protein_nomenclature(self, cdna_transcript, c_dot):
        """
        Use Mutalyzer to request the p. nomenclature for a transcript 
        """
        transcript_c_dot = f"{cdna_transcript}:{c_dot}"
        response_data, error_code = self._get_api_response(transcript_c_dot)
        
        if response_data:
            info_code = self._get_response_success_code(response_data)
            cdna_accession, protein_accession, p_dot3_raw, gene = self._get_protein_values(response_data)
            
            if cdna_accession != cdna_transcript:
                raise ValueError(F"Response cDNA does not match request cDNA: {cdna_accession} != {cdna_transcript}")
            
            if protein_accession and p_dot3_raw:
                p_dot3 = self._p_dot_converter.get_remove_parenthesis(protein_accession, p_dot3_raw, is_three_letter=True)
                p_dot1 = self._p_dot_converter.get_p_dot1(protein_accession, p_dot3)
                self._protein_response_counter[f'success.{info_code}'] += 1
            else:
                self._logger.info(f"No protein values returned for {transcript_c_dot}")
                self._protein_response_counter[f'no_protein_data.{info_code}'] += 1
        else:
            # Mutalyzer returned an error message instead of protein nomenclature
            info_code = error_code
            p_dot1 = p_dot3 = protein_accession = gene = None

                
        return {
            'cdna_transcript': cdna_transcript,
            'c_dot': c_dot,
            'p_dot1': p_dot1,
            'p_dot3': p_dot3,
            'protein_transcript': protein_accession,
            'gene': gene,
            'code': info_code         
        }
    
    def _get_protein_values(self, data):
        """
        Extract protein transcript accession and p. from json data.         
        """
        regex_pattern = r"(?P<cdna_accession>NM_[^(\s]+)\((?P<protein_accession>NP_[^)]+)\):(?P<p_dot>p\.\([^)]+\))"
        match = re.search(regex_pattern, data['protein']['description'])
        if match:
            d = match.groupdict()
            return d['cdna_accession'], d['protein_accession'], d['p_dot'], data['gene_id']
        else:
            self._logger.warning(f"Unknown protein format: {data['protein']['description']}")        
            return None, None, None, None

    def _get_response_success_code(self, mutalyzer_response: dict) -> str:
        """
        Return the success code from Mutalyzer.
        Multiple codes are concatonated 
        """
        infos = mutalyzer_response['infos']
        return '_'.join([x['code'] for x in infos])
        
    def write(self, out_filename, variant_transcripts: list[VariantTranscript]):
        """
        Write variant transcripts to file  
        """
        variant_helper.write_variant_transcripts(out_filename, variant_transcripts, ['code'])
        self._logger.info(f"Wrote {len(variant_transcripts)} rows to {out_filename}")
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Write v')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--variants", help="list of variants in vcf coordinate format (csv)")
    parser.add_argument("--output", help="Nomenclature output file (csv)", required=True)
    
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    mb = MutalyzerBatch()
    
    variant_transcripts = mb.get_variant_transcripts(args.variants)
    mb.write(args.output, variant_transcripts)

if __name__ == '__main__':
    main()