'''
Sample information for RNA-seq analysis. Contains:
    - name: sample name
    - condition: condition for DE analysis
    - files: list of fastq files belonging to sample
    - covariates: list of covariates the sample has for use in R analysis
'''

import logging
import re
import sys
import itertools
from glob import glob
from collections import defaultdict

class Sample(object):
    '''Class for sample information'''
    def __init__(self, name, condition, covariates = None, files = None):
        self.name = name
        self.condition = condition
        self.files = files if files else []
        self.covariates = covariates if covariates else []
        self.technical_replicates = contains_technical_replicates(self.files)

    def info(self):
        str = "\n".join(["\t%s" % self.name,
                         "\t\tCondition:  %s" % self.condition,
                         "\t\tCovariates: %s" % self.covariates,
                         "\t\tTech repl:  %s" % self.technical_replicates,
                         "\t\tFiles:      %s" % \
                         "\n\t\t            ".join(self.files)])
        return str

def get_all_fastq(seq_dir):
    '''Get a list of all fastq files from the specified seq_dir.'''
    seqfiles = glob("{}/*.fastq.gz".format(seq_dir))
    seqfiles.sort()
    if not seqfiles:
        logging.critical("Error: No fastq.gz files in seq_dir " \
            "{}".format(seq_dir))
    logging.debug("Files in seq_dir {}: {}".format(seq_dir, ",".join(seqfiles)))
    return seqfiles

def group_fastq_on_sm(seqfiles):
    '''Group fastq files from the same sample (e.g. paired data, technical
    replicates) into a dictionary with the SM tag as the key.'''
    sm_dict = defaultdict(list)
    for filename in seqfiles:
        try:
            # Check that the file is named correctly
            sm = re.search('(.+\/)?SM_([A-Za-z0-9-.]+)_ID_[A-Za-z0-9-.]+_' \
                'LB_[A-Za-z0-9-.]+_L[0-9]+_R[12].fastq.gz', filename).group(2)
            sm_dict[sm].append(filename)
        except AttributeError:
            logging.warning("FASTQ file {} is not in the correct name format. " \
                  "File will not be included in analysis.".format(filename))
    if not sm_dict:
        logging.critical("Error: No fastq.gz files included in analysis")
    logging.debug("Sample file groups: {}".format(sm_dict))
    return sm_dict

def read_csv(filename):
    '''Read in CSV file and remove comments after '#' and empty lines'''
    csv_list = []
    with open(filename) as f:
        for line in f:
            line = line.partition('#')[0].strip()
            if line:
                csv_list.append(line)
    return csv_list

def parse_samples(seq_dir, samples_csv):
    '''Parse sample information from the specified samples.csv file
    into a dictionary with conditions as keys.'''
    # Get all fastq files in directory
    seqfiles = get_all_fastq(seq_dir)
    # Group files on sample name
    sm_dict = group_fastq_on_sm(seqfiles)
    # Parse sample.csv file
    sample_csv_list = read_csv(samples_csv)
    sample_dict = defaultdict(list)
    for line in sample_csv_list:
        try:
            info = [x.strip() for x in line.split(",")]
            if re.search("[+-]", ",".join(info[1:])):
                logging.critical("Error: Non-allowed characters (+,-) in " \
                    "sample CSV file condition or covariate column.")
                sys.exit(1)
            # Remove SM tag if included in CSV file
            name = re.sub("^SM_", "", info[0])
            condition = info[1]
            covariates = info[2:]
            sample = Sample(name=name, condition=condition,
                covariates=covariates, files=sm_dict[name])
            sample_dict[condition].append(sample)
        except IndexError:
            logging.critical("Incorrectly formatted samples.csv file.")
            sys.exit(1)
        except KeyError as e:
            logging.warning("Sample {} has no fastq files in seq_dir.").format(e.args[0])
    logging.debug("Sample dictionary: {}".format(sample_dict))
    return sample_dict

def parse_conditions(comparison_csv, sample_dict):
    '''Parse comparisons from the specified comparison.csv file if provided.
    If not provided, use all pairwise comparisons.'''
    if comparison_csv:
        comparisons_csv_list = read_csv(comparison_csv)
        comparisons_list = [tuple(x.split(",")) for x in comparisons_csv_list]
    else:
        logging.info("No comparison_csv file provided. Using all pairwise comparisons.")
        comparisons_list = list(itertools.combinations(sample_dict.keys(), 2))
    logging.debug(comparisons_list)
    return comparisons_list

def unique_list(x):
    '''Remove duplicate items from a list'''
    return list(set(x))

def check_conditions(comparisons_list, sample_dict):
    '''Check conditions file'''
    sample_conditions = sample_dict.keys()
    comparison_conditions = [c for comparison in comparisons_list for c in comparison]
    comparison_conditions = unique_list(comparison_conditions)

    # Check all comparisons have two items
    if not all([len(c) == 2 for c in comparisons_list]):
        logging.error("Some comparisons do not have exactly 2 conditions. " \
                "Edit your comparison_csv file and rerun.")
        sys.exit(1)

    # Check conditions are concordant
    for c in sample_conditions:
        if c not in comparison_conditions:
            logging.warning("Condition {} is listed in samples_csv but has " \
                    "no comparisons in comparisons_csv. Samples with " \
                    "condition {} will not used in any comparison".format(c))
    for c in comparison_conditions:
        if c not in sample_conditions:
            logging.critical("Condition {} is listed in comparisons_csv but " \
                    "has no samples in samples_csv.".format(c))
            sys.exit(1)

    # Check number of replicates per condition
    for c in sample_dict:
        if len(sample_dict[c]) == 1:
            logging.warning("Condition {} has only 1 replicate. Downstream " \
                    "analyses with this sample may fail.")
    return

def contains_technical_replicates(files):
    '''Check files and return True if there are technical replicates'''
    r1_match = [re.match(r".*_R1\.fastq\.gz$", f) for f in files]
    replicates = sum([x is not None for x in r1_match])
    return replicates > 1
