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

from rinc.util import chromosome_map, vcf_to_gdot
from rinc.variant_transcript import VariantTranscript
from hgvs.transcriptmapper import TranscriptMapper
from rinc.util.log_config import LogConfig
from hgvs.exceptions import HGVSDataNotAvailableError, HGVSInvalidIntervalError,\
    HGVSInvalidVariantError
from rinc.io import variant_helper
from rinc.util.tx_eff_pysam import PysamTxEff

ASSEMBLY_VERSION = "GRCh37"

class HgvsNomenclature(object):
    '''
    classdocs
    '''
    def __init__(self, fasta_file):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._pysqm_txeff = PysamTxEff(fasta_file)
    
    def _get_exon(self, hdp: hgvs.dataproviders.uta.UTA_postgresql, var_c: hgvs.sequencevariant.SequenceVariant, refseq_chromosome: str) -> int:
        """
        Return exon for a transcript position
        UTA's exon start at zero.
        """
        tm = TranscriptMapper(hdp, var_c.ac, alt_ac=refseq_chromosome, alt_aln_method='splign')

        c_position = var_c.posedit.pos.start.base
        c_position_0 = c_position - 1 
        
        exon_ord = None
        note = None
        
        for exon in tm.tx_exons:
            if exon['tx_start_i'] <= c_position_0 < exon['tx_end_i']:
                # Add one to make it one based since uta exons are zero based. 
                exon_ord = exon['ord'] + 1

            
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
        am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=ASSEMBLY_VERSION, replace_reference=False)

        try:
            for var in variants:
                try:
                    
                    var.g_dot, variant_type = vcf_to_gdot.get_gdot_plus(var.chromosome, var.position, var.reference, var.alt, self._pysqm_txeff)
                    var.additional_fields['variant_type'] = variant_type
                    
                    var_g = hp.parse_hgvs_variant(var.g_dot)
                    var_c = am.g_to_c(var_g, var.cdna_transcript)                        
                    var_p = am.c_to_p(var_c)
                    
                    # Don't wrap p. in parenthesis
                    if var_p.posedit: 
                        var_p.posedit.uncertain = False
                    
                    var.c_dot = var_c.format().replace(var_c.ac + ':', '')
                    var.protein_transcript = var_p.ac
                    var.p_dot1 = var_p.format(conf={"p_3_letter": False}).replace(var_p.ac + ':', '')
                    var.p_dot3 = var_p.format(conf={"p_3_letter": True}).replace(var_p.ac + ':', '')
                    
                    # Determine exon the way tfx does it
                    var.exon = self.get_exon_the_tfx_way(hdp, var)                    
                    exon_other, var.strand, note = self._get_exon(hdp, var_c, chromosome_map.get_refseq(var.chromosome))
                    var.additional_fields['exon_other'] = exon_other
                    
                    if note: 
                        var.notes.append(note)                    
                    
                    tx_info = hdp.get_tx_info(var_c.ac, chromosome_map.get_refseq(var.chromosome), 'splign')
                    var.gene = tx_info['hgnc']
                except HGVSDataNotAvailableError as e:
                    self._logger.warning(f"Unable to parse {var}, but it's ok, this happens when sequence not found in local SeqRepo: {e}")
                    var.notes.append("Encountered HGVSDataNotAvailableError.")                    
                except HGVSInvalidIntervalError as e:
                    self._logger.warning(f"Unable to parse {var}: {e}")
                    var.notes.append(f"Encountered HGVSInvalidIntervalError: {e}")
                except HGVSInvalidVariantError as e:
                    self._logger.warning(f"Invalid variant {var}: {e}")
                    var.notes.append(f"Encountered HGVSInvalidVariantError: {e}")
                except Exception as e:
                    self._logger.error(f"Error processing variant {var}: {e}")
                    raise
        finally:
            hdp.close()
            
    def get_exon_the_tfx_way(self, hdp, variant: VariantTranscript):
        """
        Determine the exon number using the tfx method
        """
        exons = hdp.get_tx_exons(variant.cdna_transcript,
                                 chromosome_map.get_refseq(variant.chromosome),
                                 'splign')
        # self.exon_list_hdp_get_tx_exons[f"{variant.chromosome}-{variant.position}-{variant.reference}-{variant.alt}-{variant.cdna_transcript}"] = exons
         
        return self._xxx_get_uta_exon(variant.cdna_transcript, variant.position, exons)
    
    
    

    def _xxx_get_uta_exon(self, primary_accession, genomic_position, exons):
        '''
        Return the exon at or closest to the given position on a transcript. If an exon is equidistant to two exon boundaries then the downstream/3' exon is chosen. 
        * primary_accession - The RefSeq transcript accession (e.g., NM_000100.1)
        * alternate_accession - The genomic accession which is the RefSeq chromosome (e.g., NC_000001.11), not the NCBI chromosome number (chr1).         
        * position - The c. position (eg 100 from c.100A>G)
        
        Returns first coding exon (Annovar style)
        '''
        # Convert position to zero based because UTA positions are zero-based, right-exclusive
        genomic_position_0 = genomic_position - 1

        # In order to return Annovar-style numbering which starts counting the first coding exon as Exon 1 is is necessary to
        # figure which exon UTA says is the first coding exon.
        coding_start_exon_ord = None

        # When the variant is within exon boundries this will be the exon number
        variant_exon_ord = None

        # For intronic variants we keep track of the nearest exon
        nearest_exon_ord = None
        nearest_exon_distance = None

        for exon in exons:

            # uta positions are zero-based
            genomic_start_0 = exon['alt_start_i']
            genomic_end_0 = exon['alt_end_i']

            # If position is between exon start and end then this is an exonic variant and we found our exon
            if genomic_start_0 <= genomic_position_0 <= genomic_end_0:
                variant_exon_ord = exon['ord']

            # Keep track of the first coding exon.
            if exon['tx_start_i'] == 0:
                coding_start_exon_ord = exon['ord']

            # Keep track of the nearest exon in case position is intronic.
            nearest_exon_ord, nearest_exon_distance = self._xxx_get_nearest_exon(genomic_position_0,
                                                                             genomic_start_0,
                                                                             genomic_end_0,
                                                                             exon['ord'],
                                                                             nearest_exon_ord,
                                                                             nearest_exon_distance)

        # Should not happen
        if coding_start_exon_ord is None:
            raise ValueError(f"Could not find the coding start exon (tx_start_i=0) for {primary_accession}")

        # Adjust the UTA exon number to be Annovar style first coding exon
        variant_exon_ord_adjusted = None
        if variant_exon_ord is not None:
            variant_exon_ord_adjusted = (variant_exon_ord - coding_start_exon_ord) + 1
        else:
            variant_exon_ord_adjusted = (nearest_exon_ord - coding_start_exon_ord) + 1

        return variant_exon_ord_adjusted
    
    def _xxx_get_nearest_exon(self, genomic_position, exon_start, exon_end, exon_ord, current_nearest_exon_ord, current_nearest_distance):
        '''
        Return the exon ordinal and distance that is closest to variant_position. Exon boundaries are left-inclusive, right exclusive. When the current 
        distance is the same distance as the distance to the exon then the downstream exon is chosen (see HGVS 3' rule).
        * genomic_position - the location from which to calculate distance    
        * exon_start, exon_end, exon_ord - the new exon to consider
        * current_nearest_exon_ord, current_nearest_distance - the current selection for closest exon     
        '''
        # To calculate distance to end subtract one because range is right-exclusive, meaning the last base in the exon is one less than the exon_end
        distance_to_end = abs(genomic_position - (exon_end - 1))
        if (current_nearest_exon_ord is None) or (distance_to_end < current_nearest_distance) or (distance_to_end == current_nearest_distance and exon_ord > current_nearest_exon_ord):
            current_nearest_exon_ord = exon_ord
            current_nearest_distance = distance_to_end

        distance_to_start = abs(genomic_position - exon_start)
        if (distance_to_start < current_nearest_distance) or (distance_to_start == current_nearest_distance and exon_ord > current_nearest_exon_ord):
            current_nearest_exon_ord = exon_ord
            current_nearest_distance = distance_to_start

        return current_nearest_exon_ord, current_nearest_distance
    
         
        
    def write(self, out_filename, variant_transcripts: list[VariantTranscript]):
        headers = ['chromosome', 'position', 'reference', 'alt',
                   'cdna_transcript', 'hu.protein_transcript', 
                   'hu.exon', 'hu.exon_other', 'hu.strand', 'hu.gene', 
                   'hu.g_dot', 'hu.c_dot', 'hu.p_dot1', 'hu.p_dot3', 'hu.variant_type',
                   'hu.note']
        
        with open(out_filename, 'w', newline='') as output:
            writer = csv.writer(output)
            writer.writerow(headers)

            for v in variant_transcripts:
                writer.writerow([v.chromosome, v.position, v.reference, v.alt, 
                                 v.cdna_transcript, v.protein_transcript, 
                                 v.exon, None if 'exon_other' not in v.additional_fields else v.additional_fields['exon_other'], 
                                 v.strand, v.gene,
                                 v.g_dot, v.c_dot, v.p_dot1, v.p_dot3, v.additional_fields['variant_type'],
                                 ";".join(v.notes)])
            
        self._logger.info(f"Wrote {len(variant_transcripts)} variant transcripts to {out_filename}")
    

def _parse_args():
    parser = argparse.ArgumentParser(description='Use hgvs to determine protein changes for a list of variant_helper')    
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--fasta", help="hg fasta", required=True)
    parser.add_argument("--variants", help="File with variants (csv)", required=True)
    parser.add_argument("--out", help="output file (csv)", required=True)
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)

    args = _parse_args()

    hn = HgvsNomenclature(args.fasta)

    # Read variants from csv 
    variants = variant_helper.get_variants(args.variants)
    
    logging.debug(f"Read {len(variants)} variants from {args.variants}")

    # Update variant_helper with c. and p. changes 
    hn.update_with_changes(variants)

    hn.write(args.out, variants)
    
if __name__ == '__main__':
    main()        