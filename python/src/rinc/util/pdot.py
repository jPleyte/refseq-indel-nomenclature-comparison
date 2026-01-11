'''
Convert one letter p. to three letter and vis versa
Created on Jan 10, 2026

@author: pleyte
'''
import hgvs.parser
import logging

class PDot(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._hp = hgvs.parser.Parser()
    
    def get_p_dot1(self, transcript, p_dot3: str) -> str:
        '''
        Convert 3-letter p. to 1-letter
        '''
        var_p = self._hp.parse_p_variant(f"{transcript}:{p_dot3}")
        p_dot1 = var_p.format(conf={"p_3_letter": False})
        return p_dot1.split(':')[1]
        
    def get_p_dot3(self, transcript, p_dot1: str) -> str:
        '''
        Convert 3-letter p. to 1-letter
        '''
        
        var_p = self._hp.parse_p_variant(f"{transcript}:{p_dot1}")
        p_dot3 = var_p.format(conf={"p_3_letter": True})
        return p_dot3.split(':')[1]
        