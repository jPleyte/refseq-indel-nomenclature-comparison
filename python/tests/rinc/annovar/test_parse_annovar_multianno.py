'''
Created on Jan 19, 2026

@author: pleyte
'''
import unittest
from rinc.annovar.parse_annovar_multianno import AnnovarMultianno

class TestAnnovarMultianno(unittest.TestCase):


    @classmethod
    def setUpClass(cls):
        cls._annm = AnnovarMultianno()
        
    def test_get_normalized_c_dot(self):
        self.assertEqual(self._annm._get_normalized_c_dot('c.C1307T'), 'c.1307C>T', "c dot not normalized correctly")
        self.assertEqual(self._annm._get_normalized_c_dot('c.1883T>C'), 'c.1883T>C', "c dot not normalized correctly")
    
    def test_get_normalized_p_dot1(self):
        self.assertEqual(self._annm._get_normalized_p_dot1('p.E1598X'), 'p.E1598*', "p dot not normalized correctly")
        self.assertEqual(self._annm._get_normalized_p_dot1('p.K297fsX26'), 'p.K297fs*26', "p dot not normalized correctly")
        self.assertEqual(self._annm._get_normalized_p_dot1('p.W815LfsX8'), 'p.W815Lfs*8', "p dot not normalized correctly")
        
        self.assertEqual(self._annm._get_normalized_p_dot1('p.G413Sfs*238'), 'p.G413Sfs*238', "p dot should not change with normalization")
        
        self.assertEqual(self._annm._get_normalized_p_dot1('p.*342LextX?'), 'p.*342Lext*?', "p dot not normalized correctly")
        self.assertEqual(self._annm._get_normalized_p_dot1('p.E273K'), 'p.E273K', "p dot not normalized correctly")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()