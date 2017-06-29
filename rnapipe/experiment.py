'''

'''

import logging
from rnapipe.samples import Sample, TechnicalReplicate
from rnapipe.utils import read_csv, check_paired_files
from rnapipe.constants import *
import re
from glob import glob
from collections import defaultdict

class Experiment(object):
    '''
    '''
    def __init__(self, state):
        ''''''
        # Config options 
        self.paired_end = state.config.get_options("paired_end")
        self.trim_reads = state.config.get_options("trim_reads")
        self.alignment_method = state.config.get_options("alignment_method").lower()
        self.count_method = state.config.get_options("count_method").lower()

        # Check if index was provided in the config
        self.index_provided = any([(self.alignment_method == "star" and \
                                 state.config.get_options("star_index")),
                             (self.alignment_method == "hisat2" and \
                               state.config.get_options("hisat_index"))])
        logging.debug("Index provided: " + str(self.index_provided))

        # Get a dict of all samples plus their attributes grouped by condition
        self.samples = self.parse_samples(state)
        logging.debug("Sample dictionary: {}".format(self.samples))

        # Get a list of the samples
        self.sample_list = [s for l in self.samples.values() for s in l]

        # Check if there are multiple technical replicates per sample
        self.multiple_technical_replicates =  any(
                [s.multiple_technical_replicates for s in self.sample_list])
        
        # Get a flattened list of TechnicalReplicates
        tr_list = [s.technical_replicates for s in self.sample_list]
        tr_list = [tr for l in tr_list for tr in l]

        # Format TechnicalReplicates into a dictionary with ID as the key
        self.tr_dict = dict([(tr.replicate_id, tr) for tr in tr_list])

        # Get a list of R1 and R2 paths for inputs
        self.R1_files = [tr.R1_path for tr in tr_list]
        if self.paired_end:
            self.R2_files = [tr.R2_path for tr in tr_list]

        # Check
        self.check_paired_files()
        self.check_valid_alignment_method()

    def check_valid_alignment_method(self):
        '''Check if the alignment_method and count_method specifed in the 
        config file is a valid option'''
        if self.alignment_method not in ["star", "hisat2"]:
            logging.critical("Error: Invalid alignment_method in config file. " \
                    "Valid options are ['STAR', 'HISAT2'].")
            exit(EXIT_INVALID_ARGUMENT)
        if self.count_method not in ["featurecounts", "htseq-count"]:
            logging.critical("Error: Invalid count_method in config file. " \
                    "Valid options are ['featureCounts', 'HTSeq-count'].")
            exit(EXIT_INVALID_ARGUMENT)

    def check_paired_files(self):
        files_1 = [re.sub("_R1.fastq.gz$", "", x) for x in self.R1_files]
        files_2 = [re.sub("_R2.fastq.gz$", "", x) for x in self.R2_files]
        if self.paired_end and not (sorted(files_1) == sorted(files_2)):
            logging.critical("Incorrect pairings")
            sys.exit(99)   ## TODO
        return

    def get_all_fastq(self, seq_dir):
        '''Get a list of all fastq files from the specified seq_dir.'''
        seqfiles = glob("{}/*.fastq.gz".format(seq_dir))
        seqfiles.sort()
        if not seqfiles:
            logging.critical("Error: No fastq.gz files in seq_dir " \
                "{}".format(seq_dir))
        logging.debug("Files in seq_dir {}: {}".format(seq_dir, ",".join(seqfiles)))
        return seqfiles
    
    def parse_samples_csv(self, samples_csv):
        '''Return samples in a nested list'''
        sample_csv_list = read_csv(samples_csv)
        sample_list = []
        for line in sample_csv_list:
            info = [x.strip() for x in line.split(",")]
            if len(info) < 2:
                logging.critical("Incorrectly formatted samples.csv file.")
                sys.exit(EXIT_INVALID_CSV)
            if re.search("[+-]", ",".join(info[1:])):
                logging.critical("Error: Non-allowed characters (+,-) in " \
                    "sample CSV file condition or covariate column.")
                sys.exit(EXIT_INVALID_CSV)
            # Remove SM tag if included in CSV file
            name = re.sub("^SM_", "", info[0])
            sample_list.append(info)
        # Make sure all lines have the same number of fields
        if len(set([len(x) for x in sample_list])) != 1:
            logging.critical("The number of fields in the samples.csv file " \
                "is inconsistent.")
            sys.exit(EXIT_INVALID_CSV)
        return sample_list
    
    def get_technical_replicates(self, filenames, paired_end):
        ''''''
        # Get R1 files
        R1_files = [x for x in filenames if re.search("_R1.fastq.gz", x)]
        # For each R1 file, find if there's a matching R2 file and create 
        # a TechnicalReplicate object
        technical_replicates = []
        for f1 in R1_files:
            if paired_end:
                f2 = re.sub("_R1.fastq.gz$", "_R2.fastq.gz", f1)
                try:
                    filenames.remove(f2)
                    filenames.remove(f1)
                    tr = TechnicalReplicate(f1, f2)
                except ValueError:
                    logging.warning("Not using {} because it has no matching pair" \
                            "".format(f1))
            else:
                filenames.remove(f1)
                tr = TechnicalReplicate(f1)
            technical_replicates.append(tr)
        # Check filenames list is empty
        if filenames:
            logging.warning("The following fastq files in the seq directory" \
                    "cannot not be used for analysis due to SE/PE information. " \
                    "Check if your filenames end with '_R#.fastq.gz': {}".format(filenames))
        return technical_replicates
    
    def parse_samples(self, state):
        '''Parse sample information from the specified samples.csv file
        into a dictionary with conditions as keys.'''
        paired_end = state.config.get_options("paired_end")
        # Get all fastq files in directory
        seqfiles = self.get_all_fastq(state.config.get_options("seq_dir"))
        # Parse sample.csv file
        sample_csv_list = self.parse_samples_csv(state.config.get_options("samples_csv"))
        # Create TechnicalReplicate objects for seq files
        technical_replicates = self.get_technical_replicates(seqfiles, paired_end)
        # Create sample dictionary where each condition is a key
        sample_dict = defaultdict(list)
        for info in sample_csv_list:
            name = info[0]
            tr = [x for x in technical_replicates if x.sample_name == name]
            if tr:
                sample = Sample(name=name, condition=info[1], covariates=info[2:],
                                technical_replicates=tr)
                sample_dict[info[1]].append(sample)
                logging.debug(sample)
            else:
                logging.warning("Sample {} has no fastq files in seq_dir.").format(e.args[0])
        return sample_dict

    def get_alignment_structure(self):
        '''Return structure for ruffus inputs and outputs for alignment'''
        if self.alignment_method == "star":
            return
        return




#    def parse_conditions(comparison_csv, sample_dict):
#        '''Parse comparisons from the specified comparison.csv file if provided.
#        If not provided, use all pairwise comparisons.'''
#        if comparison_csv:
#            comparisons_csv_list = read_csv(comparison_csv)
#            comparisons_list = [tuple(x.split(",")) for x in comparisons_csv_list]
#        else:
#            logging.info("No comparison_csv file provided. Using all pairwise comparisons.")
#            comparisons_list = list(itertools.combinations(sample_dict.keys(), 2))
#        logging.debug(comparisons_list)
#        return comparisons_list
#    
#    def check_conditions(comparisons_list, sample_dict):
#        '''Check conditions file'''
#        sample_conditions = sample_dict.keys()
#        comparison_conditions = [c for comparison in comparisons_list for c in comparison]
#        comparison_conditions = unique_list(comparison_conditions)
#    
#        # Check all comparisons have two items
#        if not all([len(c) == 2 for c in comparisons_list]):
#            logging.error("Some comparisons do not have exactly 2 conditions. " \
#                    "Edit your comparison_csv file and rerun.")
#            sys.exit(EXIT_INVALID_CSV)
#    
#        # Check conditions are concordant
#        for c in sample_conditions:
#            if c not in comparison_conditions:
#                logging.warning("Condition {} is listed in samples_csv but has " \
#                        "no comparisons in comparisons_csv. Samples with " \
#                        "condition {} will not used in any comparison".format(c))
#        for c in comparison_conditions:
#            if c not in sample_conditions:
#                logging.critical("Condition {} is listed in comparisons_csv but " \
#                        "has no samples in samples_csv.".format(c))
#                sys.exit(EXIT_INVALID_CSV)
#    
#        # Check number of replicates per condition
#        for c in sample_dict:
#            if len(sample_dict[c]) == 1:
#                logging.warning("Condition {} has only 1 replicate. Downstream " \
#                        "analyses with this sample may fail.")
#        return

