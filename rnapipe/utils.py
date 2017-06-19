'''
Various utility functions that don't have a sensible home elsewhere
'''
import os
import logging

def path_list_join(dir, file_list):
    '''Join directory to a list of files'''
    return [os.path.join(dir, x) for x in file_list]

def run_picard(state, stage, args):
    mem = int(state.config.get_stage_options(stage, "mem"))
    return run_java(state, stage, PICARD_JAR, mem, args)

def run_trimmomatic(state, stage, args):
    mem = int(state.config.get_stage_options(stage, "mem"))
    return run_java(state, stage, TRIMMOMATIC_JAR, mem, args)

def create_empty_outputs(outputs):
    '''Create empty dummy files for testing purposes'''
    if isinstance(outputs, str):
        outputs = [outputs]
    for output_filename in outputs:
        with open(output_filename, "w"):
            pass

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
