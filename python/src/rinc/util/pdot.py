'''
Convert one letter p. to three letter and vis versa
Created on Jan 10, 2026

@author: pleyte
'''
import hgvs.parser
import logging
from hgvs.exceptions import HGVSParseError
import re

aa3_to_1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 
    'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 
    'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*', 
    'Trp': 'W', 'Gln': 'Q' 
}

class PDot(object):
    '''
    Convert protein changes between one and three letter 
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._hp = hgvs.parser.Parser()
    
    def get_remove_parenthesis(self, transcript: str, p_dot: str, is_three_letter=True) -> str:
        """
        Remove the parenthesis indicating that the p. is "inferred" 
        """
        try:
            var_p = self._hp.parse_p_variant(f"{transcript}:{p_dot}")
            var_p.posedit.uncertain = False            
            p_dot = var_p.format(conf={"p_3_letter": is_three_letter})
            return p_dot.split(':')[1]
        except HGVSParseError as e:
            self._logger.warning(f"Error parsing {p_dot} for transcript {transcript}, will use backup method: {e}")

    def get_p_dot1(self, transcript: str, p_dot3: str) -> str:
        '''
        Convert 3-letter p. to 1-letter
        '''
        try:
            var_p = self._hp.parse_p_variant(f"{transcript}:{p_dot3}")
            p_dot1 = var_p.format(conf={"p_3_letter": False})
            return p_dot1.split(':')[1]
        except HGVSParseError as e:
            self._logger.warning(f"Error parsing {p_dot3} for transcript {transcript}, will use backup method: {e}")
        
        # Weird p. cause errors in hgvs so we resort to regex mapping        
        return self.get_map_pdot3_to_pdot1(p_dot3)
            

    def get_p_dot3(self, transcript, p_dot1: str) -> str:
        '''
        Convert 3-letter p. to 1-letter
        '''        
        var_p = self._hp.parse_p_variant(f"{transcript}:{p_dot1}")
        p_dot3 = var_p.format(conf={"p_3_letter": True})
        return p_dot3.split(':')[1]
    
    def get_map_pdot3_to_pdot1(self, p_dot3):
        """
        Use regex and an amno acid abbreviation map to convert three letter p. to one letter p.
        """
        pattern = re.compile('|'.join(aa3_to_1.keys()))
        return pattern.sub(lambda x: aa3_to_1[x.group()], p_dot3)
    
if __name__ == '__main__':
    pdot_converter = PDot()
    assert pdot_converter.get_map_pdot3_to_pdot1("p.TrpGln1973*") == "p.WQ1973*"
    