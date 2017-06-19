import unittest
from rnapipe.samples import * 

class TestSeqFile(unittest.TestCase):
    '''Unit tests for SeqFile'''
    def test_basic_SE_filename(self):
        filename = "sample-1.fastq.gz"
        sf = SeqFile(filename)
        self.assertEqual(sf.name, "sample-1")
        self.assertEqual(sf.R1_or_R2, DEFAULT_METADATA["r"])
        self.assertEqual(sf.lane, DEFAULT_METADATA["lane"])
        self.assertEqual(sf.id, DEFAULT_METADATA["id"])
        self.assertEqual(sf.library, DEFAULT_METADATA["lb"])

    def test_basic_PE_filename(self):
        filename = "x1-rep1_R2.fastq.gz"
        sf = SeqFile(filename)
        self.assertEqual(sf.name, "x1-rep1")
        self.assertEqual(sf.R1_or_R2, "R2")
        self.assertEqual(sf.lane, DEFAULT_METADATA["lane"])
        self.assertEqual(sf.id, DEFAULT_METADATA["id"])
        self.assertEqual(sf.library, DEFAULT_METADATA["lb"])

    def test_mixed_filename(self):
        filename = "xxx_ID_123-A_L001_R1.fastq.gz"
        sf = SeqFile(filename)
        self.assertEqual(sf.name, "xxx")
        self.assertEqual(sf.R1_or_R2, "R1")
        self.assertEqual(sf.lane, "L001")
        self.assertEqual(sf.id, "123-A")
        self.assertEqual(sf.library, DEFAULT_METADATA["lb"])

    def test_full_filename(self):
        filename = "SM_sample_ID_aaa_LB_111_L008_R2.fastq.gz"
        sf = SeqFile(filename)
        self.assertEqual(sf.name, "sample")
        self.assertEqual(sf.R1_or_R2, "R2")
        self.assertEqual(sf.lane, "L008")
        self.assertEqual(sf.id, "aaa")
        self.assertEqual(sf.library, "111")

    def test_fail_sample_filename(self):
        filename = "sample_A_L001_R1.fastq.gz"
        sf = SeqFile(filename)
        self.assertEqual(sf.name, None)
        self.assertEqual(sf.R1_or_R2, None)
        self.assertEqual(sf.lane, None)
        self.assertEqual(sf.id, None)
        self.assertEqual(sf.library, None)

    def test_fail_tags_filename(self):
        filename = "sample-A_LB_123_ID_A_L001_R1.fastq.gz"
        sf = SeqFile(filename)
        self.assertEqual(sf.name, None)
        self.assertEqual(sf.R1_or_R2, None)
        self.assertEqual(sf.lane, None)
        self.assertEqual(sf.id, None)
        self.assertEqual(sf.library, None)

if __name__ == '__main__':
    unittest.main()

