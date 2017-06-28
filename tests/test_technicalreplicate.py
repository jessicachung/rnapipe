import unittest
from rnapipe.samples import * 

class TestTechincalReplicate(unittest.TestCase):
    '''Unit tests for TechincalReplicate class'''
    def test_basic_SE_filename(self):
        filename = "/dir/sample-1_R1.fastq.gz"
        tr = TechnicalReplicate(filename)
        self.assertEqual(tr.R1_path, "/dir/sample-1_R1.fastq.gz")
        self.assertEqual(tr.R2_path, None)
        self.assertEqual(tr.sample_name, "sample-1")
        self.assertEqual(tr.replicate_id, "sample-1")
        self.assertEqual(tr.lane, DEFAULT_METADATA["lane"])
        self.assertEqual(tr.id, DEFAULT_METADATA["id"])
        self.assertEqual(tr.library, DEFAULT_METADATA["lb"])
        self.assertEqual(tr.is_PE, False)
        self.assertEqual(tr.R1_trimmed, "sample-1_R1.trimmed.fastq.gz")
        self.assertEqual(tr.R2_trimmed, None)

    def test_basic_PE_filename(self):
        filename_1 = "/dir/x1-rep1_R1.fastq.gz"
        filename_2 = "/dir/x1-rep1_R2.fastq.gz"
        tr = TechnicalReplicate(filename_1, filename_2)
        self.assertEqual(tr.R1_path, "/dir/x1-rep1_R1.fastq.gz")
        self.assertEqual(tr.R2_path, "/dir/x1-rep1_R2.fastq.gz")
        self.assertEqual(tr.sample_name, "x1-rep1")
        self.assertEqual(tr.replicate_id, "x1-rep1")
        self.assertEqual(tr.lane, DEFAULT_METADATA["lane"])
        self.assertEqual(tr.id, DEFAULT_METADATA["id"])
        self.assertEqual(tr.library, DEFAULT_METADATA["lb"])
        self.assertEqual(tr.is_PE, True)
        self.assertEqual(tr.R1_trimmed, "x1-rep1_R1.trimmed.fastq.gz")
        self.assertEqual(tr.R2_trimmed, "x1-rep1_R2.trimmed.fastq.gz")

    def test_mixed_filename(self):
        filename = "seq/xxx_ID_123-A_L001_R1.fastq.gz"
        tr = TechnicalReplicate(filename)
        self.assertEqual(tr.R1_path, "seq/xxx_ID_123-A_L001_R1.fastq.gz")
        self.assertEqual(tr.R2_path, None)
        self.assertEqual(tr.sample_name, "xxx")
        self.assertEqual(tr.replicate_id, "xxx_ID_123-A_L001")
        self.assertEqual(tr.lane, "L001")
        self.assertEqual(tr.id, "123-A")
        self.assertEqual(tr.library, DEFAULT_METADATA["lb"])
        self.assertEqual(tr.is_PE, False)
        self.assertEqual(tr.R1_trimmed, "xxx_ID_123-A_L001_R1.trimmed.fastq.gz")
        self.assertEqual(tr.R2_trimmed, None)

    def test_full_filename(self):
        filename_1 = "data/SM_sample_ID_aaa_LB_111_L008_R1.fastq.gz"
        filename_2 = "data/SM_sample_ID_aaa_LB_111_L008_R2.fastq.gz"
        tr = TechnicalReplicate(filename_1, filename_2)
        self.assertEqual(tr.sample_name, "sample")
        self.assertEqual(tr.replicate_id, "SM_sample_ID_aaa_LB_111_L008")
        self.assertEqual(tr.lane, "L008")
        self.assertEqual(tr.id, "aaa")
        self.assertEqual(tr.library, "111")
        self.assertEqual(tr.is_PE, True)


    def test_fail_sample_filename(self):
        filename = "sample_A_L001_R1.fastq.gz"
        tr = TechnicalReplicate(filename)
        self.assertEqual(tr.sample_name, None)
        self.assertEqual(tr.lane, None)
        self.assertEqual(tr.id, None)
        self.assertEqual(tr.library, None)

    def test_fail_tags_filename(self):
        filename = "sample-A_LB_123_ID_A_L001_R1.fastq.gz"
        tr = TechnicalReplicate(filename)
        self.assertEqual(tr.sample_name, None)
        self.assertEqual(tr.lane, None)
        self.assertEqual(tr.id, None)
        self.assertEqual(tr.library, None)

if __name__ == '__main__':
    unittest.main()

