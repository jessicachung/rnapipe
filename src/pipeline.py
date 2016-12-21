'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import PipelineStages
from pipeline_base.utils import safe_make_dir
from samples import parse_samples, contains_technical_replicates
from os import path
from utils import path_list_join
import error_codes
import logging

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.DEBUG)

def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name="rnapipe")
    # Stages are dependent on the state
    stages = PipelineStages(state)

    # Get a list of all samples plus their attributes (condition, covariates, etc.)
    samples = parse_samples(seq_dir=state.config.get_options("seq_dir"),
                            samples_csv=state.config.get_options("samples_csv"))
    # Get a list of the samples
    sample_list = [s for sample_list in samples.values() for s in sample_list]
    # Check if there are multiple technical replicates per sample
    technical_replicates =  any([s.technical_replicates for s in sample_list])
    # Get a list of paths to all the FASTQ files
    fastq_files = [s.files for s in sample_list]
    fastq_files = [item for sublist in fastq_files for item in sublist]

    # Print out samples
    sample_text = [s.info() for s in sample_list]
    logging.debug("Analysis samples:\n{}".format("\n".join(sample_text)))

    # Make directories
    # NOTE: Move this to somewhere else?
    def get_output_paths(state):
        results_dir = state.config.get_options("results_dir")
        analysis_name = state.config.get_options("pipeline_id")
        output_path = {
            "fastqc": "qc/fastqc",
            "post_trim_fastqc": "qc/post_trim_fastqc",
            "seq": "seq",
            "alignments": "alignments",
            "star_index": "reference/star/",
            "hisat_index": "reference/hisat/",
            "counts": "counts",
            "voom": "analysis/voom_{name}".format(name=analysis_name),
            "edgeR": "analysis/edgeR_{name}".format(name=analysis_name),
            "cuffdiff": "analysis/cuffdiff_{name}".format(name=analysis_name)
        }
        output_path = [(a, path.join(results_dir, b)) for a,b in output_path.items()]
        output_path = dict(output_path)
        return output_path

    output_dir = get_output_paths(state)
    print(output_dir)

    safe_make_dir(output_dir["alignments"])
    safe_make_dir(output_dir["seq"])

    # Config options
    paired_end = state.config.get_options("paired_end")
    trim_reads = state.config.get_options("trim_reads")
    alignment_method = state.config.get_options("alignment_method").lower()

    # Check if alignment_method is valid
    if alignment_method not in ["star", "hisat2"]:
        print("Error: Invalid alignment_method in config file. " \
              "Valid options are ['STAR', 'HISAT2'].")
        exit(error_codes.INVALID_ARGUMENT)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.do_nothing,
        name="original_fastqs",
        output=fastq_files)

    pipeline.originate(
        task_func=stages.do_nothing,
        name="reference_genome",
        output=state.config.get_options("reference_genome"))

    # Pre-trim FastQC
    pipeline.transform(
        task_func=stages.fastqc,
        name="fastqc",
        input=output_from("original_fastqs"),
        filter=suffix(".fastq.gz"),
        output="_fastqc.zip",
        output_dir=output_dir["fastqc"],
        extras=[output_dir["fastqc"]])

    # Trimmomatic
    if trim_reads and paired_end:
        pipeline.transform(
            task_func=stages.trim_reads_pe,
            name="trim_reads_pe",
            input=output_from("original_fastqs"),
            # Get R1 file and the corresponding R2 file
            filter=formatter(".+/(?P<sample>[a-zA-Z0-9_]+)_R1.fastq.gz"),
            add_inputs=add_inputs("{path[0]}/{sample[0]}_R2.fastq.gz"),
            output=path_list_join(output_dir["seq"],
                       ["{sample[0]}_R1.trimmed.fastq.gz",
                        "{sample[0]}_R1.unpaired.fastq.gz",
                        "{sample[0]}_R2.trimmed.fastq.gz",
                        "{sample[0]}_R2.unpaired.fastq.gz"]))
    elif trim_reads:
        pipeline.transform(
            task_func=stages.trim_reads_se,
            name="trim_reads_se",
            input=output_from("original_fastqs"),
            filter=formatter(".+/(?P<sample>[a-zA-Z0-9_]+)_R1.fastq.gz"),
            output=path.join(output_dir["seq"],
                             "{sample[0]}_R1.trimmed.fastq.gz"))

    # Post-trim FastQC
    if paired_end and trim_reads:
        pipeline.transform(
            task_func=stages.fastqc,
            name="post_trim_fastqc_pe",
            input=output_from("trim_reads_pe"),
            filter=formatter(".+/(?P<sample>[a-zA-Z0-9_]+)_R1.trimmed.fastq.gz"),
            output=path_list_join(output_dir["post_trim_fastqc"],
                       ["{sample[0]}_R1.trimmed_fastqc.gz",
                        "{sample[0]}_R2.trimmed_fastqc.gz"]),
            extras=["results/qc/post_trim_fastqc/"])
    elif trim_reads:
        pipeline.transform(
            task_func=stages.fastqc,
            name="post_trim_fastqc_se",
            input=output_from("trim_reads_se"),
            filter=suffix(".trimmed.fastq.gz"),
            output=".trimmed_fastqc.gz",
            output_dir=output_dir["post_trim_fastqc"],
            extras=[output_dir["post_trim_fastqc"]])

    # Create reference index for STAR
    if alignment_method == "star" and not \
            state.config.get_options("star_index"):
        pipeline.transform(
            task_func=stages.create_star_index,
            name="create_star_index",
            input=output_from("reference_genome"),
            filter=suffix(""),
            add_inputs=add_inputs(state.config.get_options("gene_ref")),
            output=path_list_join(output_dir["star_index"],
                       ["SA", "Genome", "genomeParameters.txt"]),
            extras=[output_dir["star_index"]])
    elif alignment_method == "star":
        # Don't create index if index is supplied
        # NOTE: May need to make the timestamp for the index more recent than
        #       the reference files. Or change to an originate task instead.
        output_dir["star_index"] = state.config.get_options("star_index")
        pipeline.transform(
            task_func=stages.do_nothing,
            name="create_star_index",
            input=output_from("reference_genome"),
            filter=suffix(""),
            add_inputs=add_inputs(state.config.get_options("gene_ref")),
            output=path_list_join(output_dir["star_index"],
                       ["SA", "Genome", "genomeParameters.txt"]))

    # Create reference index for HISAT2
    if alignment_method == "hisat" and not \
            state.config.get_options("hisat_index"):
        hisat_basename = path.join(output_dir["hisat_index"], "genome")
        pipeline.transform(
            task_func=stages.create_hisat_index,
            name="create_hisat_index",
            input=state.config.get_options("reference_genome"),
            add_inputs=add_inputs(state.config.get_options("gene_ref")),
            output=path_list_join(output_dir["hisat_index"],
                       ["genome.1.ht2", "genome.2.ht2"]),
            extras=[hisat_basename])
    elif alignment_method == "hisat":
        # Don't create index if index is supplied
        hisat_basename = state.config.get_options("hisat_index")
        output_dir["hisat_index"] = path.dirname(hisat_basename)
        prefix = path.basename(hisat_basename)
        pipeline.transform(
            task_func=stages.do_nothing,
            name="create_hisat_index",
            input=state.config.get_options("reference_genome"),
            add_inputs=add_inputs(state.config.get_options("gene_ref")),
            output=path_list_join(output_dir["hisat_index"],
                       ["{prefix}.1.ht2".format(prefix=prefix),
                        "{prefix}.2.ht2".format(prefix=prefix)]),
            extras=[hisat_basename])


    if paired_end and trim_reads:
        trim_task_name = "trim_reads_pe"
    elif not paired_end and trim_reads:
        trim_task_name = "trim_reads_se"
    elif paired_end:
        trim_task_name = "original_fastqs"

    # If there are technical replicates, each is mapped independently.
    # This is so each technical replicate maintains a separate read group.
    # NOTE: Should the read group metadata be stored elsewhere instead of in filenames?
    if alignment_method == "star":
        align_task_name = "star_align"
        (pipeline.transform(
            task_func=stages.star_align,
            name=align_task_name,
            input=output_from(trim_task_name),
            filter=formatter(".+/SM_(?P<sm>[a-zA-Z0-9-]+)_ID_" \
                             "(?P<id>[a-zA-Z0-9-]+)_LB_(?P<lb>[a-zA-Z0-9]+)" \
                             "_L(?P<ln>[0-9]+)_R1(.trimmed)?.fastq.gz"),
            output="%s/{sm[0]}/SM_{sm[0]}_ID_{id[0]}_LB_{lb[0]}" \
                   "_L{ln[0]}.star.bam" % output_dir["alignments"],
            extras=[output_dir["star_index"], "{sm[0]}", "{id[0]}", "{lb[0]}",
                    "{ln[0]}"])
        ).follows("create_star_index")

    if alignment_method == "hisat2":
        align_task_name = "hisat_align"
        (pipeline.transform(
            task_func=stages.hisat_align,
            name="hisat_align",
            input=output_from(trim_task_name),
            filter=formatter(".+/SM_(?P<sm>[a-zA-Z0-9-]+)_ID_" \
                             "(?P<id>[a-zA-Z0-9-]+)_LB_(?P<lb>[a-zA-Z0-9]+)" \
                             "_L(?P<ln>[0-9]+)_R1(.trimmed)?.fastq.gz"),
            output="%s/{sm[0]}/SM_{sm[0]}_ID_{id[0]}_LB_{lb[0]}" \
                   "_L{ln[0]}.hisat2.bam" % output_dir["alignments"],
            extras=["{sm[0]}", "{id[0]}", "{lb[0]}", "{ln[0]}"])
        ).follows("create_hisat_index")


    # Sort BAM by coordinates
    pipeline.transform(
        task_func=stages.sort_bam_by_coordinate,
        name="sort_bam_by_coordinate",
        input=output_from(align_task_name),
        filter=suffix(".bam"),
        output=".sorted.bam")

    # Merge files with the same sample name
    # if technical_replicates:
        # merge_task_name = "merge_bams"
    pipeline.collate(
        task_func=stages.merge_bams,
        name="merge_bams",
        input=output_from(align_task_name),
        filter=formatter(".+/SM_(?P<sm>[a-zA-Z0-9-]+)_ID_.+.bam"),
        output="%s/{sm[0]}/{sm[0]}.bam" % output_dir["alignments"])

    # Sort BAM by name for HTSeq-count
    pipeline.transform(
        task_func=stages.sort_bam_by_name,
        name="sort_bam_by_name",
        input=output_from("merge_bams"),
        filter=suffix(".bam"),
        output=".nameSorted.bam")

    return pipeline
