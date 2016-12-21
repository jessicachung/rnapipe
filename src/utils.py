'''
Various utility functions that don't have a sensible home elsewhere
'''
from os import path

def path_list_join(dir, file_list):
    '''Join directory to a list of files'''
    return [path.join(dir, x) for x in file_list]

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
