import unittest
import cgr
import os
from scipy.sparse import dok_matrix
import numpy

class CGRTest(unittest.TestCase):
    
    def test_guess_genome_size(self):
        self.assertEqual(cgr.estimate_genome_size("test_data/10kfile"),10240)
        
    def test_run_jellyfish(self):
        result = cgr.run_jellyfish("test_data/NC_012920.fasta", 11, "test_data/tmp.jf")
        os.remove("test_data/tmp.jf")
        self.assertEqual(result.returncode, 0)
        
    def test_get_kmer_list(self):
        cgr.run_jellyfish("test_data/NC_012920.fasta", 11, "test_data/tmp.jf")
        counts = []
        for i in cgr.get_kmer_list("test_data/tmp.jf"):
            counts.append(i)
        self.assertEqual(len(counts), 16438)
        self.assertEqual(len(counts[0]), 2)
        os.remove("test_data/tmp.jf")
        
    def test_get_coord(self):
        self.assertEqual(cgr.get_coord("A"), [0,0])
        self.assertEqual(cgr.get_coord("G"), [1,1])
        self.assertEqual(cgr.get_coord("C"), [1,0])
        self.assertEqual(cgr.get_coord("T"), [0,1])
        self.assertEqual(cgr.get_coord("CTGA"), [10,6])
    
    def test_get_k(self):
        self.assertEqual(cgr.get_k("test_data/five_mer.jf"), 5)
    
    def test_cgr_matrix(self):
        m = cgr.cgr_matrix("test_data/five_mer.jf")
        row1_col1 = m.toarray()[0][0]
        self.assertEqual(row1_col1, 70)
        
    def test_join_cgr(self):
        a = [ [1,1,1], [1,1,1] ]
        a = dok_matrix(a)
        b = [ [2,2,2], [2,2,2] ]
        b = dok_matrix(b)
        c = numpy.array( [ [  [1,2],  [1,2], [1,2]  ], [ [1,2], [1,2], [1,2] ] ] )
        d = cgr.join_cgr(a,b)
        self.assertEqual(d.tolist(), c.tolist() )
        
    def test_pad_two_channel_to_three(self):
        a =  numpy.array( [ [  [1,2],  [1,2], [1,2]  ], [ [1,2], [1,2], [1,2] ] ] )
        b = numpy.array( [ [  [1,2,0],  [1,2,0], [1,2,0]  ], [ [1,2,0], [1,2,0], [1,2,0] ] ])
        c = cgr.pad_two_channel_to_three(a)
        self.assertEqual(c.tolist(), b.tolist() )
    
if __name__ == '__main__':
    unittest.main()