'''
Sample information for RNA-seq analysis. Contains:
    - name: sample name
    - condition: condition for DE analysis
    - files: list of fastq files belonging to sample
    - covariates: list of covariates the sample has for use in R analysis
'''

import logging
import os
import re
import sys
import itertools
from glob import glob
from collections import defaultdict

# TODO: Move this elsewhere
DEFAULT_METADATA = {"id": "id", "lb": "lb", "r": "R1", "lane": None}

class Sample(object):
    '''Class for sample information'''
    def __init__(self, name, condition, covariates = None, files = None):
        self.name = name
        self.condition = condition
        self.files = files if files else []
        self.covariates = covariates if covariates else []
        self.technical_replicates = self.contains_technical_replicates(self.files)

    def __repr__(self):
        str = "Sample(name={}, condition={}, files={}, covariates={}, technical_replicates={})".format(
                self.name, self.condition, self.files, self.covariates, self.technical_replicates)
        return str

    def info(self):
        "Pretty representation of the Sample object"
        filenames = [f.path for f in self.files]
        str = "\n".join(["\t%s" % self.name,
                         "\t\tCondition:  %s" % self.condition,
                         "\t\tCovariates: %s" % self.covariates,
                         "\t\tTech repl:  %s" % self.technical_replicates,
                         "\t\tFiles:      %s" % \
                         "\n\t\t            ".join(filenames)])
        return str

    def contains_technical_replicates(self, files):
        '''Check files and return True if there are technical replicates'''
        is_R1 = [f.R1_or_R2 == "R1" for f in files]
        replicates = sum(is_R1)
        return replicates > 1

class SeqFile(object):
    '''Class for sequencing file'''
    def __init__(self, path):
        self.path = path
        self.name = None
        self.R1_or_R2 = None
        self.lane = None
        self.id = None
        self.library = None

        # Parse metadata using filename
        self.parse_metadata()

    def __repr__(self):
        str = "SeqFile(name={}, R1_or_R2={}, lane={}, id={}, library={}, path={})".format(
                self.name, self.R1_or_R2, self.lane, self.id, self.library, 
                self.path)
        return str

    def get_metadata(self, match, re_group):
        '''Retreive metadata from regex match and parse. If no metadata
        provided, use defaults.'''
        result = match.group(re_group)
        if result:
            value = result.split("_")[-1]
        else:
            value = DEFAULT_METADATA[re_group]
        return value

    def parse_metadata(self):
        metadata_string = re.sub(".fastq(.gz)?$", "", os.path.basename(self.path))
        metadata_list = metadata_string.split("_")
        try:
            metadata_match = re.match("^(SM_)?(?P<sm>[A-Za-z0-9-]+)(?P<id>_ID_[A-Za-z0-9-]+)?(?P<lb>_LB_[A-Za-z0-9-]+)?(?P<lane>_L[0-9]+)?(?P<r>_R[1,2])?$", metadata_string)
            logging.debug(metadata_match.groups())
        except AttributeError:
            logging.warning("FASTQ file {} is not in the correct name format. " \
                  "File will not be included in analysis.".format(self.path))
            return
        # Set new values
        self.name = self.get_metadata(metadata_match, "sm")
        self.R1_or_R2 = self.get_metadata(metadata_match, "r")
        self.lane = self.get_metadata(metadata_match, "lane")
        self.id = self.get_metadata(metadata_match, "id")
        self.library = self.get_metadata(metadata_match, "lb")
        return


def get_all_fastq(seq_dir):
    '''Get a list of all fastq files from the specified seq_dir.'''
    seqfiles = glob("{}/*.fastq.gz".format(seq_dir))
    seqfiles.sort()
    if not seqfiles:
        logging.critical("Error: No fastq.gz files in seq_dir " \
            "{}".format(seq_dir))
    logging.debug("Files in seq_dir {}: {}".format(seq_dir, ",".join(seqfiles)))
    return seqfiles


def group_fastq_on_sample_name(seqfiles):
    '''Group fastq files from the same sample (e.g. paired data, technical
    replicates) into a dictionary with the SM tag as the key.'''
    sm_dict = defaultdict(list)
    for filename in seqfiles:
        logging.debug("# PARSING FILENAME: " + filename)
        f = SeqFile(filename)
        logging.debug(f)
        if f.name:
            sm_dict[f.name].append(f)
    if not sm_dict:
        logging.critical("Error: No fastq.gz files included in analysis")
    logging.debug("Sample file groups: {}".format(sm_dict))
    return sm_dict

def read_csv(filename):
    '''Read in CSV file and remove comments after '#' and empty lines'''
    csv_list = []
    with open(filename) as f:
        contents = f.read().split("\n")
    # Remove comments and empty lines
    for line in contents:
        line = line.partition('#')[0].strip()
        if line:
            csv_list.append(line)
    return csv_list

def parse_samples_csv(samples_csv):
    '''Return samples in a nested list'''
    sample_csv_list = read_csv(samples_csv)
    sample_list = []
    for line in sample_csv_list:
        info = [x.strip() for x in line.split(",")]
        if len(info) < 2:
            logging.critical("Incorrectly formatted samples.csv file.")
            sys.exit(1)
        if re.search("[+-]", ",".join(info[1:])):
            logging.critical("Error: Non-allowed characters (+,-) in " \
                "sample CSV file condition or covariate column.")
            sys.exit(1) # TODO: Change to predefined exit codes
        # Remove SM tag if included in CSV file
        name = re.sub("^SM_", "", info[0])
        sample_list.append(info)
    # Make sure all lines have the same number of fields
    if len(set([len(x) for x in sample_list])) != 1:
        logging.critical("The number of fields in the samples.csv file " \
            "is inconsistent.")
        sys.exit(1)
    return sample_list

def parse_samples(seq_dir, samples_csv):
    '''Parse sample information from the specified samples.csv file
    into a dictionary with conditions as keys.'''
    # Get all fastq files in directory
    seqfiles = get_all_fastq(seq_dir)
    # Group files on sample name
    seqfile_dict = group_fastq_on_sample_name(seqfiles)
    # Parse sample.csv file
    sample_csv_list = parse_samples_csv(samples_csv)
    # Create sample dictionary where each condition is a key
    sample_dict = defaultdict(list)
    for info in sample_csv_list:
        try:
            name = info[0]
            logging.debug(seqfile_dict[name])
            sample = Sample(name=name, condition=info[1],
                            covariates=info[2:], files=sm_dict[name])
            sample_dict[info[1]].append(sample)
            logging.debug(sample)
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

