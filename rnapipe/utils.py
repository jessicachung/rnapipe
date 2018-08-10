'''
Various utility functions that don't have a sensible home elsewhere
'''
import os
import logging
from rnapipe.pipeline_base.utils import safe_make_dir

def path_list_join(dir, file_list):
    '''Join directory to a list of files'''
    return [os.path.join(dir, x) for x in file_list]

def create_empty_outputs(outputs):
    '''Create empty dummy files for testing purposes'''
    if isinstance(outputs, str):
        outputs = [outputs]
    for output_filename in outputs:
        with open(output_filename, "w"):
            pass

def get_output_paths(results_dir, default_paths):
    '''Join base result directory path to output paths'''
    output_paths = [(a, os.path.join(results_dir, b)) for a,b in 
            default_paths.items()]
    output_paths = dict(output_paths)
    return output_paths

def make_output_dirs(output_dict):
    '''Create directory for each value in the dictionary'''
    for dir in output_dict.values():
        safe_make_dir(dir)

def check_paired_files(R1, R2):
    files_1 = [re.sub("_R1.fastq.gz$", "", x) for x in R1]
    files_2 = [re.sub("_R2.fastq.gz$", "", x) for x in R2]
    if not (sorted(files_1) == sorted(files_2)):
        logging.critical("Incorrect pairings")
        sys.exit(99)   ## TODO
    return

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

# Copied from http://www.ruffus.org.uk/faq.html
def re_symlink(input_file, soft_link_name):
    '''
    Helper function: relinks soft symbolic link if necessary
    '''
    # Guard agains soft linking to oneself: Disastrous consequences of
    # deleting the original files!!
    if input_file == soft_link_name:
        logging.debug("Warning: No symbolic link made. You are using the " \
                      "original data directory as the working directory.")
        return
    # Soft link already exists: delete for relink?
    if os.path.lexists(soft_link_name):
        # do not delete or overwrite real (non-soft link) file
        if not os.path.islink(soft_link_name):
            raise Exception("%s exists and is not a link" % soft_link_name)
        try:
            os.unlink(soft_link_name)
        except:
            logging.debug("Can't unlink %s" % (soft_link_name))
    logging.debug("os.symlink(%s, %s)" % (input_file, soft_link_name))
    #
    #   symbolic link relative to original directory so that the entire path
    #       can be moved around with breaking everything
    #
    os.symlink( os.path.relpath(os.path.abspath(input_file),
                os.path.abspath(os.path.dirname(soft_link_name))), soft_link_name)
