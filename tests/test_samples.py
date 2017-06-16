import unittest
from rnapipe.samples import * 

class TestSamples(unittest.TestCase):
    '''Unit tests for Samples'''
    def test_sample(self):
        f1 = SeqFile("sample.fastq.gz")
        samp = Sample(name="sample", condition="c1", covariates=[], files=[f1])
        self.assertEqual(samp.name, "sample")
        self.assertEqual(samp.condition, "c1")
        self.assertEqual(len(samp.files), 1)
        self.assertEqual(samp.technical_replicates, False)

    def test_PE_sample(self):
        f1_R1 = SeqFile("sample_ID_x_L001_R1.fastq.gz")
        f1_R2 = SeqFile("sample_ID_x_L001_R2.fastq.gz")
        samp = Sample(name="sample", condition="c1", covariates=[], files=[f1_R1, f1_R2])
        self.assertEqual(samp.name, "sample")
        self.assertEqual(samp.condition, "c1")
        self.assertEqual(len(samp.files), 2)
        self.assertEqual(samp.technical_replicates, False)

    def test_replicate_samples(self):
        f1 = SeqFile("sample_L001_R1.fastq.gz")
        f2 = SeqFile("sample_L002_R1.fastq.gz")
        samp = Sample(name="sample", condition="c1", covariates=[], files=[f1, f2])
        self.assertEqual(samp.name, "sample")
        self.assertEqual(samp.condition, "c1")
        self.assertEqual(len(samp.files), 2)
        self.assertEqual(samp.technical_replicates, True)


if __name__ == '__main__':
    unittest.main()

