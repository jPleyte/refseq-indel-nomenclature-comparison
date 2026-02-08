'''
Generate variant-transcript batch files for submission to Mutalyzer.
We don't want to submit all variants to mutalzyer.  We just want the ones where 
the two input dataframes have the same c_dot values. 
 
Bad variant code key
- EINTRONIC: "Intronic position `54-107` was used with a non intronic reference sequence `NM_001018837.2`."
- ESPLICESITE: Variant affects one or more splice sites.
- EOUTOFBOUNDARY: Position `*70901` is invalid; it is 70557 nucleotides after the end of the sequence (sequence length is 1006 nucleotides).

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
import os
import csv
import time
from _collections import deque

PRINT_REQUEST_STATS = True

class MutalyzerNomenclatureRequester(object):
    '''
    Create batch files for submission to Mutalyzer 
    '''
    def __init__(self, nomenclature_output_file: str, bad_variants_output_file: str, transcript_filter: str=None, rate_limit=20):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._file_nomenclature_output = nomenclature_output_file
        self._file_bad_variants_output = bad_variants_output_file
        self._file_handle_nomenclature = None
        self._file_handle_bad_variants = None
         
        # Optional cDNA transcript filter
        self._transcript_filter_filename = transcript_filter
        self._transcript_filter = None
        self._transcript_filter_counter = Counter()
        
        # If processing is interrupted, we can resume where we left off 
        self._already_processed_variants = set()
        self._bad_variants = set()
         
        self._p_dot_converter = PDot()
        self._variant_counter = Counter()
        self._api_response_counter = Counter()
        self._protein_response_counter = Counter()
        
        self._mutalizer_normalize_url = "https://mutalyzer.nl/api/normalize/{value}"
        
        self._nomenclature_field_names = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript', 'c_dot', 'exon', 'g_dot', 'gene', 'genomic_region_type', 'p_dot1', 'p_dot3', 'protein_transcript', 'protein_variant_type', 'code']
        self._bad_variant_field_names = ['chromosome', 'position', 'reference', 'alt', 'code']
        
        # Tracking requests per minute and latency
        self._request_timestamps = deque()
        self._latencies = deque(maxlen=100)
        self._per_variant_request_counter = None
        self._per_variant_transcript_filter_counter = None
        
        # Rate limit for api requests
        self._requests_per_minute = rate_limit
        self._minimum_delay = 60.0 / self._requests_per_minute if rate_limit else None
        self._last_request_time = 0.0
    
    @property
    def current_rpm(self) -> int:
        """Returns Requests Per Minute based on the last 60s."""
        return len(self._request_timestamps)
    
    @property
    def avg_latency(self) -> float:
        """Returns average response time in seconds."""
        if not self._latencies:
            return 0.0
        return sum(self._latencies) / len(self._latencies)
    
    def __enter__(self):
        # Look for a previous run 
        self._update_already_processed_variants()
        self._update_bad_variants()
        
        # Open file handles for nomenclature and bad variants 
        self._logger.debug(f"Openning nomenclature file {self._file_nomenclature_output}")
        self._file_handle_nomenclature = open(self._file_nomenclature_output, 'a', newline='', buffering=1)
        
        self._logger.debug(f"Openning bad variant file {self._file_bad_variants_output}")
        self._file_handle_bad_variants = open(self._file_bad_variants_output, 'a', newline='', buffering=1)
        
        # If a transcript filter was given create a set of transcripts  
        if self._transcript_filter_filename:
            with open(self._transcript_filter_filename, 'r') as f:
                # strip() removes the newline character \n
                self._transcript_filter = {val for line in f if (val := line.strip()) and val.startswith("NM")}
                
            self._logger.info(f"Created cDNA transcript filter with {len(self._transcript_filter)} distinct transcripts")
    
        # Use a session so we don't have to redo the ssl handshake for each request 
        self._session = requests.Session()
        
        return self  
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._file_handle_nomenclature:
            self._file_handle_nomenclature.close()
        
        if self._file_handle_bad_variants:
            self._file_handle_bad_variants.close()

        self._session.close()
        
        return False
    
    def _record_request(self, duration: float):
        """
        Record how long the last request took
        """
        now = time.time()
        self._request_timestamps.append(now)
        self._latencies.append(duration)
        
        # Remove timestamps older than 60 seconds to keep the window accurate
        while self._request_timestamps and self._request_timestamps[0] < now - 60:
            self._request_timestamps.popleft()
            
    def _update_already_processed_variants(self):
        """
        If a previous run of the module was stopped before all variants were complete, the progress can be loaded and we'll resume
        where we left off. This function reads the previously started nomenclature and adds all variants to the already_processed_variants
        list. 
        """
        if os.path.exists(self._file_nomenclature_output):
            with open(self._file_nomenclature_output, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    self._already_processed_variants.add(self._get_variant_identifier(row))
                    self._variant_counter['already_processed'] += 1
            
            self._logger.info(f"Resuming... {len(self._already_processed_variants)} variants already processed: {self._file_nomenclature_output}")
        else:
            self._logger.info(f"Previous nomenclature not detected")


    def _update_bad_variants(self):
        """
        If a previous run of the module was stopped before all variants were complete, the progress can be loaded and we'll resume
        where we left off. This function loads the previously started bad variants file.  
        """
        if os.path.exists(self._file_bad_variants_output):
            with open(self._file_bad_variants_output, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    self._bad_variants.add(self._get_variant_identifier(row))
            
            self._logger.info(f"Resuming... Read {len(self._bad_variants)} already known bad variants from: {self._file_bad_variants_output}")
        else:
            self._logger.info(f"Previous bad variants not detected")
        
    def _get_variant_identifier(self, v):
        if type(v) == dict:
            return f"{v['chromosome']}-{v['position']}-{v['reference']}-{v['alt']}"
        else:
            return f"{v.chromosome}-{v.position}-{v.reference}-{v.alt}"
        
    def fetch_and_write_variant_transcripts(self, variants_file: str) -> list[VariantTranscript]:
        """        
        Receives a list of variants, and uses Mutalyzer's normalizer API to get transcripts and nomenclature.
        Returns a list of of variant-transcripts      
        """
        variants = variant_helper.get_variants(variants_file)
        self._logger.info(f"Read {len(variants)} variants from {variants_file}")
        
        nomenclature_writer = csv.DictWriter(self._file_handle_nomenclature, fieldnames=self._nomenclature_field_names)
        bad_variants_writer = csv.DictWriter(self._file_handle_bad_variants, fieldnames=self._bad_variant_field_names)
        
        if os.path.getsize(self._file_nomenclature_output) == 0:
            nomenclature_writer.writeheader()
        if os.path.getsize(self._file_bad_variants_output) == 0:
            bad_variants_writer.writeheader()

        n = 0
        for v in variants:
            n = n + 1
            if n % 100 == 0:
                self._logger.info(f"Processed {n}/{len(variants)} variants")

            # Determine if we've already processed this variant 
            if self._get_variant_identifier(v) in self._already_processed_variants:
                self._variant_counter['already_processed'] += 1
                continue


            try:
                self._per_variant_request_counter = 0
                self._per_variant_transcript_filter_counter = 0

                start_time = time.perf_counter()                
                variant_transcripts, error_code = self._get_variant_transcripts(v)
                variant_duration = time.perf_counter() - start_time

                if error_code:
                    self.write_bad_variant(bad_variants_writer, v, error_code)
                    self._variant_counter['bad'] += 1
                else:
                    self._write_nomenclature_rows(nomenclature_writer, variant_transcripts)
                    self._variant_counter['new'] += 1

                    if PRINT_REQUEST_STATS:
                        print(f"Variant {n} involved {self._per_variant_request_counter} requests and took {variant_duration:.2f}s "
                              f"[RPM: {self.current_rpm}] "
                              f"[Avg Latency: {self.avg_latency:.2f}s] "
                              f"[Filter excluded: {self._per_variant_transcript_filter_counter} ]"
                        )
            except Exception as e:
                self._logger.warning(f"Unable to process variant {self._get_variant_identifier(v)}: {e}")
                self._variant_counter['error'] += 1
                self.write_bad_variant(bad_variants_writer, v, str(e))
        
        summary_lines = [
            "Completed processing:"
            f"  Variants: {self._variant_counter}"
            f"  Api responses: {self._api_response_counter}"
            f"  Protein results: {self._protein_response_counter}"
        ]
        
        if self._transcript_filter:
            summary_lines.append(f"  Transcript filter: {self._transcript_filter_counter}")
        
        self._logger.info("\n".join(summary_lines))

    def write_bad_variant(self, writer, v: VariantTranscript, error_code: str):
        """
        """
        writer.writerow({
            'chromosome': v.chromosome,
            'position': v.position,
            'reference': v.reference, 
            'alt': v.alt,
            'code': error_code
        })
    
    def _write_nomenclature_rows(self, writer, variant_transcripts: list[VariantTranscript]):
        """
        Append variant transcripts to the output file 
        """
        for v in variant_transcripts:            
            writer.writerow({
                'chromosome': v.chromosome,
                'position': v.position,
                'reference': v.reference, 
                'alt': v.alt,
                'cdna_transcript': v.cdna_transcript, 
                'c_dot': v.c_dot,
                'exon': v.exon,
                'g_dot': v.g_dot,
                'gene': v.gene, 
                'genomic_region_type': v.genomic_region_type, 
                'p_dot1': v.p_dot1, 
                'p_dot3': v.p_dot3, 
                'protein_transcript': v.protein_transcript, 
                'protein_variant_type': v.protein_variant_type,
                'code': v.additional_fields.get('code')
            })
        
    def _get_variant_transcripts(self, v: VariantTranscript) -> tuple[list[VariantTranscript], str]:
        """
        Take a single variant and -
        1. Convert to g. and request all c. and transcripts for that g. from mutalyzer using their API 
        2. For every c. and transcript, request p. from mutalyzer
        
        Returns a list of all the transcripts for a single variant
         or an error code if there was a problem.  
        """
        variant_transcripts = []
        
        simple_g_dot = self._get_simple_g_dot(v)
        transcript_cdots, error_code = self._get_transcripts(simple_g_dot)
        
        if error_code and transcript_cdots:
            raise ValueError(f"Error code should only exist when transcripts can not be found. See {v}")
        
        for t in transcript_cdots:
            cp_nomenclature = self._get_cdna_and_protein_nomenclature(chromosome_map.get_refseq(v.chromosome), 
                                                                      t['cdna_transcript'], 
                                                                      t['c_dot'])
            
            variant_transcripts.append(self._get_variant_transcript(v, t['g_dot'], cp_nomenclature))
        
        return variant_transcripts, error_code
            

    def _get_variant_transcript(self, v: VariantTranscript, g_dot: str, cdna_protein_values: dict) -> VariantTranscript:
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
        if len(v.reference) == 1:
            return f"{chromosome}:g.{v.position}{v.reference}>{v.alt}"
        else:
            start_pos = v.position
            end_pos = v.position + len(v.reference) - 1
            return f"{chromosome}:g.{start_pos}_{end_pos}{v.reference}>{v.alt}"
    
    def _get_transcripts(self, simple_g_dot) -> list:
        """
        Use request a list of refseq transcripts from Mutalyzer API 
        """
        response, error_code = self._get_api_response(simple_g_dot)
         
        transcript_cdots = []
        if response:
            corrected_g_dot = response['corrected_description']
            for c in response['equivalent_descriptions']['c']:
                refseq_transcript, c_dot = self._get_cdna_values(c['description'])
                if refseq_transcript and c_dot:
                    if self._transcript_filter is None or refseq_transcript in self._transcript_filter:
                        self._transcript_filter_counter['included'] += 1
                        transcript_cdots.append({
                            'g_dot': corrected_g_dot,
                            'cdna_transcript': refseq_transcript,
                            'c_dot': c_dot
                        })
                    elif self._transcript_filter is not None:
                        self._transcript_filter_counter['excluded'] += 1
                        self._per_variant_transcript_filter_counter += 1
                else:                    
                    raise ValueError(f"No c. for {simple_g_dot}: why not?")

        return transcript_cdots, error_code
                
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
        self._wait_for_rate_limit()
        
        self._per_variant_request_counter += 1
        start_time = time.perf_counter()
        
        encoded_parameter = requests.utils.quote(normalize_parameter)
        url = self._mutalizer_normalize_url.format(value=encoded_parameter)
        
        response = self._session.get(url, timeout=10)
        # response.raise_for_status()
        
        duration = time.perf_counter() - start_time
        self._record_request(duration)
        
        if response.status_code == 200:
            self._api_response_counter[f'{response.status_code}.success'] += 1
            return response.json(), None
        else:
            error_code = self._get_response_error_code(response.json())
            self._api_response_counter[f'{response.status_code}.{error_code}'] += 1
            self._logger.info(f"Unable to get normalizer result for {normalize_parameter}: {error_code}")
            return None, error_code

    def _wait_for_rate_limit(self):
        """
        Avoid annoying Mutalyzer by having a rate limit of no more than 20 requests per minute        
        """
        if self._requests_per_minute:
            now = time.time()
            time_since_last = now - self._last_request_time
            if time_since_last < self._minimum_delay:
                sleep_time = self._minimum_delay - time_since_last
                #self._logger.debug(f"Rate limiting: Sleeping for {sleep_time:.2f}s")
                time.sleep(sleep_time)
            
            # Update the timestamp AFTER the potential sleep
            self._last_request_time = time.time()
        
    def _get_response_error_code(self, mutalyzer_response: dict) -> str:
        """
        Return the error code from a mutalyzer reponse        
        """
        try:
            transcript_request = mutalyzer_response['custom']['corrected_description']
            errors = mutalyzer_response['custom']['errors']
            error_codes = []
            for e in errors:
                error_code = e['code']
                error_codes.append(error_code)
                self._logger.info(f"Received error code {error_code} for {transcript_request}: {e['details']}")

            return '_'.join(error_codes)
        except (KeyError, IndexError):
            self._logger.warning(f"Unable to process error response: {mutalyzer_response}. Keys not found in {mutalyzer_response}")
            raise

    def _get_c_dot_normalizer_request_string(self, refseq_chromosome: str, cdna_transcript: str, c_dot: str):
        """
        Return the string to send to Mutalyzer for the c. request. 
        - For an exonic c. like "c.123C>T" return a string like 'NM_001018837.2:c.123C>T'
        - For intronic or upstream the c. will be c.54-107del or c.100+5del in which case Mutalyzer needs the refseq chromosome in order to look at the referenc sequence. 
        """
        # is_intronic = bool(re.search(r'\d[+-]\d', c_dot))
        # if is_intronic: 
        return f"{refseq_chromosome}({cdna_transcript}):{c_dot}"
        # else:
            # return f"{cdna_transcript}:{c_dot}"
        
        
    def _get_cdna_and_protein_nomenclature(self, refseq_chromosome, cdna_transcript, c_dot):
        """
        Use Mutalyzer to request the p. nomenclature for a transcript 
        """
        request_string = self._get_c_dot_normalizer_request_string(refseq_chromosome, cdna_transcript, c_dot)
        response_data, error_code = self._get_api_response(request_string)
        
        if response_data:
            info_code = self._get_response_code(response_data)
            x_accession, protein_accession, p_dot3_raw, gene = self._get_protein_values(response_data)
            
            # Sanity check. not really neccessary
            if x_accession.startswith('NC') and x_accession != refseq_chromosome:
                raise ValueError(F"Response chromosome does not match request chromosome: {x_accession} != {refseq_chromosome}")
            elif x_accession.startswith('NM') and x_accession != cdna_transcript:
                raise ValueError(F"Response cDNA does not match request cDNA: {x_accession} != {cdna_transcript}")
            
            if protein_accession and p_dot3_raw:
                p_dot3 = self._p_dot_converter.get_remove_parenthesis(protein_accession, p_dot3_raw, is_three_letter=True)
                p_dot1 = self._p_dot_converter.get_p_dot1(protein_accession, p_dot3)
                self._protein_response_counter[f'success.{info_code}'] += 1
            else: 
                p_dot1 = p_dot3 = protein_accession = None
                self._logger.info(f"No protein values returned for {request_string} and code {info_code}")
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
        regex_pattern = r"(?P<accession>(?:NM|NC)_[^(\s]+)\((?P<protein_accession>NP_[^)]+)\):(?P<p_dot>p\..+)"
        
        if 'description' not in data['protein']:
            return data['corrected_model']['reference']['id'], None, None, data['gene_id']
            
        match = re.search(regex_pattern, data['protein']['description'])
        if match:
            d = match.groupdict()
            return d['accession'], d['protein_accession'], d['p_dot'], data['gene_id']
        else:
            raise ValueError(f"Unknown protein format: {data['protein']['description']}")

    def _get_response_code(self, mutalyzer_response: dict) -> str:
        """
        Return the success code from Mutalyzer.
        Multiple codes are concatonated 
        """
        if 'infos' in mutalyzer_response:
            infos = mutalyzer_response['infos']
            return '_'.join([x['code'] for x in infos])
        elif 'protein' in mutalyzer_response and 'errors' in mutalyzer_response['protein']:
            errors = mutalyzer_response['protein']['errors']
            return '_'.join([x['code'] for x in errors])
        else:
            return None
        
    def write(self, out_filename, variant_transcripts: list[VariantTranscript]):
        """
        Write variant transcripts to file  
        """
        variant_helper.write_variant_transcripts(out_filename, variant_transcripts, ['code'])
        self._logger.info(f"Wrote {len(variant_transcripts)} rows to {out_filename}")
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Send a list of variants to Mutalyzer using their web api and write out c. and p. nomenclature')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--variants", help="list of variants in vcf coordinate format (csv)")
    parser.add_argument("--nomenclature_output", help="Nomenclature output file (csv)", required=True)
    parser.add_argument("--bad_variants_output", help="Variants that could not be processed (csv)", required=True)
    parser.add_argument("--transcript_filter", help="Limit results to these cDNA transcripts (csv)", required=True)
    
    args = parser.parse_args()
    return args    

def main():
    # logging.config.dictConfig(LogConfig().stdout_config)
    logging.config.dictConfig(LogConfig().file_config)
    
    output = logging.root.handlers[0].baseFilename
    print(f"Log level={logging.root.getEffectiveLevel()}, output={output}")
    
    args = _parse_args()
    
    with MutalyzerNomenclatureRequester(args.nomenclature_output, args.bad_variants_output, transcript_filter=args.transcript_filter) as mnr:
        mnr.fetch_and_write_variant_transcripts(args.variants)

if __name__ == '__main__':
    main()