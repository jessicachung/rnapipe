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
from rnapipe.constants import *

class Sample(object):
    '''Class for sample information'''
    def __init__(self, name, condition, covariates = None, technical_replicates = None):
        self.name = name
        self.condition = condition
        self.technical_replicates = technical_replicates
        self.files = self.get_files()
        self.covariates = covariates if covariates else []
        self.multiple_technical_replicates = len(technical_replicates) > 1

    def __repr__(self):
        str = "Sample(name={self.name}, condition={self.condition}, " \
              "technical_replicates={self.technical_replicates}, files=" \
              "{self.files}, covariates={self.covariates}, " \
              "multiple_technical_replicates={self.multiple_technical_replicates})".format(**locals())
        return str

    def info(self):
        "Pretty representation of the Sample object"
        str = "\n".join(["\t%s" % self.name,
                         "\t\tCondition:  %s" % self.condition,
                         "\t\tCovariates: %s" % self.covariates,
                         "\t\tFiles:      %s" % \
                         "\n\t\t            ".join(self.files)])
        return str
    
    def get_files(self):
        '''Get files from TechnicalReplicate objects'''
        files = [x.R1_path for x in self.technical_replicates]
        files += [x.R2_path for x in self.technical_replicates if x.is_PE]
        return(files)


class TechnicalReplicate(object):
    '''Class for each technical replicate'''
    def __init__(self, R1, R2=None):
        self.R1_path = R1
        self.R2_path = R2
        self.sample_name = None
        self.replicate_id = None
        self.lane = None
        self.id = None
        self.library = None
        self.R1_trimmed = re.sub(".fastq.gz$", ".trimmed.fastq.gz", os.path.basename(self.R1_path))
        if self.R1_path and self.R2_path:
            self.is_PE = True
            self.R2_trimmed = re.sub(".fastq.gz$", ".trimmed.fastq.gz", os.path.basename(self.R2_path))
        else:
            self.is_PE = False
            self.R2_trimmed = None

        # Parse metadata using filename
        self.parse_metadata()

    def __repr__(self):
        str = "TechnicalReplicate(sample_name={self.sample_name}, " \
              "R1_path={self.R1_path}, R2_path={self.R2_path}, " \
              "replicate_id={self.replicate_id}, lane={self.lane}, " \
              "id={self.id}, library={self.library}, is_PE={self.is_PE}, " \
              "R1_trimmed={self.R1_trimmed}, R2_trimmed={self.R2_trimmed}" \
              "".format(**locals())
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
        metadata_string = re.sub("_R1.fastq.gz$", "", os.path.basename(self.R1_path))
        metadata_list = metadata_string.split("_")
        try:
            metadata_match = re.match("^(SM_)?(?P<sm>[A-Za-z0-9-]+)" \
                    "(?P<id>_ID_[A-Za-z0-9-]+)?(?P<lb>_LB_[A-Za-z0-9-]+)?" \
                    "(?P<lane>_L[0-9]+)?$", metadata_string)
            logging.debug("Metadata groups: {}".format(metadata_match.groups()))
        except AttributeError:
            logging.warning("FASTQ file {} is not in the correct name format. " \
                  "File will not be included in analysis.".format(self.R1_path))
            return
        # Set new values
        self.sample_name = self.get_metadata(metadata_match, "sm")
        self.replicate_id = metadata_string
        self.lane = self.get_metadata(metadata_match, "lane")
        self.id = self.get_metadata(metadata_match, "id")
        self.library = self.get_metadata(metadata_match, "lb")
        return



