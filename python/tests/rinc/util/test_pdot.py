'''
Created on Jan 10, 2026

@author: pleyte
'''
import unittest
from rinc.util.pdot import PDot


class TestPDot(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._pd = PDot()
        
    def test_get_p_dot1(self):
        self.assertEqual(self._pd.get_p_dot1("NM_000352.3", "p.Arg175His"), "p.R175H", "Three letter not converted correctly")
        
    def test_get_p_dot3(self):
        self.assertEqual(self._pd.get_p_dot3("NM_000352.3", "p.R175H"), "p.Arg175His", "One letter not converted correctly")    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()