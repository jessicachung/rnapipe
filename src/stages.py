'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

from pipeline_base.utils import safe_make_dir, run_java
from pipeline_base.stages import Stages
from os import path
from utils import *

# Barcoo
# if options.cluster == "barcoo":
TRIMMOMATIC_JAR = "/usr/local/trimmomatic/0.36/trimmomatic-0.36.jar"
PICARD_JAR = "/usr/local/picard/2.3.0/picard.jar"
java_tmp = "-Djava.io.tmpdir=$TMPDIR"

# Snowy
TRIMMOMATIC_JAR = "/usr/local/easybuild/software/Trimmomatic/0.36/trimmomatic-0.36.jar"
PICARD_JAR = "/usr/local/easybuild/software/picard/2.3.0/picard.jar"
java_tmp = "-Djava.io.tmpdir=$TMPDIR"

# Local
java_tmp = ""


class PipelineStages(Stages):
    def __init__(self, *args, **kwargs):
        super(PipelineStages, self).__init__(*args, **kwargs)
        self.reference_genome = self.get_options("reference_genome")
        self.gene_ref = self.get_options("gene_ref")
        self.adapter_seq = self.get_options("adapter_seq")
        self.trimmomatic_parameters = self.get_options("trimmomatic_parameters")

    def do_nothing(self, *args):
        '''Do nothing'''
        pass

    def fastqc(self, input, outputs, fastqc_dir):
        '''Run FastQC on fastq files'''
        safe_make_dir(fastqc_dir)
        create_empty_outputs(outputs)
        # If multiple fastq inputs, join into a string
        if isinstance(input, list):
            fastq_input = " ".join(input)
        else:
            fastq_input = input
        command = "fastqc -o {fastqc_dir} -f fastq {fastq_input}".format(
                      fastqc_dir=fastqc_dir,
                      fastq_input=fastq_input)
        # run_stage(self.state, 'fastqc', command)

    def trim_reads_pe(self, inputs, outputs):
        '''Trim PE reads with Trimmomatic'''
        safe_make_dir(path.dirname(outputs[0]))
        create_empty_outputs(outputs)
        cores = self.get_stage_options("trim_reads", "cores")
        input_R1, input_R2 = inputs
        output_R1, unpaired_R1, output_R2, unpaired_R2 = outputs
        args = "PE -threads {n_threads} -phred33 {input_R1} {input_R2} " \
               "{output_R1} {unpaired_R1} {output_R2} {unpaired_R2} " \
               "ILLUMINACLIP:{adapter_seq}:2:30:10 {parameters}".format(
                   n_threads=cores, input_R1=input_R1, input_R2=input_R2,
                   output_R1=output_R1, unpaired_R1=unpaired_R1,
                   output_R2=output_R2, unpaired_R2=unpaired_R2,
                   adapter_seq=self.adapter_seq,
                   parameters=self.trimmomatic_parameters)
        # run_trimmomatic("trim_reads_pe", args)

    def trim_reads_se(self, input, output):
        '''Trim SE reads with Trimmomatic'''
        cores = self.get_stage_options("trim_reads", "cores")
        args = "SE -threads {n_threads} -phred33 {input_fastq} " \
               "{output_fastq} ILLUMINACLIP:{adapter_seq}:2:30:10 " \
               "{parameters}".format(n_threads=cores, input_fastq=input,
                   output_fastq=output, adapter_seq=self.adapter_seq,
                   parameters=parameters)
        # run_trimmomatic("trim_reads_se", args)

    def create_star_index(self, inputs, outputs, output_dir):
        '''Generate index for STAR'''
        safe_make_dir(output_dir)
        create_empty_outputs(outputs)
        genome_fa, gene_gtf = inputs
        cores = self.get_stage_options("align", "cores")
        command = "STAR --runThreadN {n_threads} --runMode genomeGenerate " \
                  "--genomeDir {output_dir} --genomeFastaFiles {genome_fa} " \
                  "--sjdbGTFfile {gene_gtf}".format(output_dir=output_dir,
                      genome_fa=genome_fa, gene_gtf=gene_gtf)
        # run_stage(self.state, "create_star_index", command)

    def create_hisat_index(self, inputs, outputs, hisat_basename):
        '''Generate index for HISAT2'''
        safe_make_dir(path.dirname(hisat_basename))
        create_empty_outputs(outputs)
        genome_fa, gene_gtf = inputs
        cores = self.get_stage_options("build_index", "cores")
        command = "hisat2-build -p {n_threads} {genome_fa} {hisat_basename}" \
                  "".format(n_threads=cores, genome_fa=genome_fa, gene_gtf=gene_gtf)
        # run_stage(self.state, "create_hisat_index", command)

    def star_align(self, inputs, outputs, ref_dir, sm, id, lb, ln):
        '''Align fastq files with STAR'''
        safe_make_dir(path.dirname(outputs))  # TODO: Change to list if multiple outputs
        create_empty_outputs(outputs)
        cores = self.get_stage_options("align", "cores")
        # If PE fastq inputs, join into a string
        if isinstance(inputs, list):
            fastq_input = " ".join(inputs)
        else:
            fastq_input = inputs
        print(fastq_input)
        command = "STAR --runThreadN {cores} --genomeDir {ref_dir} " \
                  "--readFilesIn {fastq_input} --readFilesCommand zcat" \
                  "".format(cores=cores, ref_dir=ref_dir, fastq_input=fastq_input)
        # run_stage(self.state, 'star_align', command)
        # TODO: Add stranded information
        #       Include unmapped reads in the file

    def hisat_align(self, inputs, outputs, ref_basename, sm, id, lb, ln):
        '''Align fastq files with HISAT2'''
        cores = self.get_stage_options("align", "cores")
        # If PE fastq inputs, use hisat -1 and -2 arguments, else use -U
        if isinstance(inputs, list):
            fastq_input = "-1 {fastq_R1} -2 {fastq_R2}".format(
                              fastq_R1=inputs[0], fastq_R2=inputs[1])
        else:
            fastq_input = "-U {fastq}".format(fastq=inputs)
        command = "hisat2 --dta-cufflinks --rg-id {id}_{ln} --rg SM:{sm} " \
                  "--rg LB:{lb} --rg PL:Illumina -x {ref_basename} " \
                  "{fastq_input}".format(id=id, ln=ln, sm=sm, lb=lb,
                      ref_basename=ref_basename, fastq_input=fastq_input)
        # run_stage(self.state, "hisat_align", command)

    def sort_bam_by_coordinate(self, input, output):
        '''Sort BAM file by coordinates'''
        create_empty_outputs(output)
        command = "samtools sort {input} {output} && samtools index {output}"
        # run_stage(self.state, "sort_bam_by_coordinate", command)

    def merge_bams(self, inputs, output):
        '''Merge multiple BAM files into one BAM file'''
        create_empty_outputs(output)
        print(inputs)
        print(output)
        if isinstance(inputs, list):
            command = "" # TODO: Merge files together with samtools
        else:
            pass # TODO: Make symbolic link instead of merging

    def sort_bam_by_name(self, input, output):
        '''Sort BAM file by name'''
        create_empty_outputs(output)
        command = ""
