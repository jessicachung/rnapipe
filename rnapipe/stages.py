'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

from rnapipe.pipeline_base.utils import safe_make_dir
from rnapipe.pipeline_base.stages import Stages
from rnapipe.pipeline_base.runner import run_stage
import os
from rnapipe.utils import *
import re


### TODO
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
# TRIMMOMATIC_JAR = "/mnt/transient_nfs/anaconda3/share/trimmomatic-0.36-3/trimmomatic.jar"
TRIMMOMATIC_JAR = "/mnt/galaxy/gvl/anaconda2/share/trimmomatic-0.36-5/trimmomatic.jar"

# Spartan
TRIMMOMATIC_JAR = "/usr/local/easybuild/software/Trimmomatic/0.36/trimmomatic-0.36.jar"

class PipelineStages(Stages):
    def __init__(self, state, experiment, *args, **kwargs):
        super(PipelineStages, self).__init__(state, *args, **kwargs)
        self.reference_genome = self.get_options("reference_genome")
        self.gene_ref = self.get_options("gene_ref")
        self.adapter_seq = self.get_options("adapter_seq")
        self.trimmomatic_parameters = self.get_options("trimmomatic_parameters")
        self.paired_end = self.get_options("paired_end")
        self.experiment = experiment

    def java_command(self, jar_path, mem_in_gb, command_args):
        '''Build a string for running a java command'''
        # Bit of room between Java's max heap memory and what was requested.
        # Allows for other Java memory usage, such as stack.
        java_mem = min(mem_in_gb - 2, 1)
        return 'java -Xmx{mem}g -jar {jar_path} {command_args}'.format(
            jar_path=jar_path, mem=java_mem, command_args=command_args)
    
    def run_java(self, state, stage, jar_path, mem, args):
        command = self.java_command(jar_path, mem, args)
        run_stage(state, stage, command)

    def run_picard(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, "mem"))
        return self.run_java(self.state, stage, PICARD_JAR, mem, args)
    
    def run_trimmomatic(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, "mem"))
        return self.run_java(self.state, stage, TRIMMOMATIC_JAR, mem, args)

    def do_nothing(self, *args):
        '''Do nothing'''
        pass

    def create_symlinks (self, input, output):
        '''
        Make symlinks in the results directory
        '''
        if isinstance(input, str):
            input = [input]
        if isinstance(output, str):
            output = [output]
        for src, dest in zip(input, output):
            re_symlink(src, dest)


    def fastqc(self, input, outputs, fastqc_dir):
        '''Run FastQC on fastq files'''
        safe_make_dir(fastqc_dir)
        # If multiple fastq inputs, join into a string
        if isinstance(input, tuple) or isinstance(input, list):
            input = " ".join(input)
        command = "fastqc -o {fastqc_dir} -f fastq {fastq_input}".format(
                      fastqc_dir=fastqc_dir,
                      fastq_input=input)
        run_stage(self.state, "fastqc", command)

    def trim_reads(self, inputs, outputs, unpaired_R1=None, unpaired_R2=None):
        '''Trim reads with Trimmomatic'''
        cores = self.get_stage_options("trim_reads", "cores")
        if self.paired_end:
            input_R1, input_R2 = inputs
            output_R1, output_R2 = outputs
            args = "PE -threads {n_threads} -phred33 {input_R1} {input_R2} " \
                   "{output_R1} {unpaired_R1} {output_R2} {unpaired_R2} " \
                   "ILLUMINACLIP:{adapter_seq}:2:30:10 {parameters}".format(
                       n_threads=cores, input_R1=input_R1, input_R2=input_R2,
                       output_R1=output_R1, unpaired_R1=unpaired_R1,
                       output_R2=output_R2, unpaired_R2=unpaired_R2,
                       adapter_seq=self.adapter_seq,
                       parameters=self.trimmomatic_parameters)
        else:
            args = "SE -threads {n_threads} -phred33 {input_fastq} " \
                   "{output_fastq} ILLUMINACLIP:{adapter_seq}:2:30:10 " \
                   "{parameters}".format(n_threads=cores, input_fastq=inputs,
                       output_fastq=outputs, adapter_seq=self.adapter_seq,
                       parameters=self.trimmomatic_parameters)
        self.run_trimmomatic("trim_reads", args)


    def create_star_index(self, inputs, outputs, output_dir):
        '''Generate index for STAR'''
        safe_make_dir(output_dir)
        genome_fa, gene_gtf = inputs
        cores = self.get_stage_options("align", "cores")
        command = "STAR --runThreadN {n_threads} --runMode genomeGenerate " \
                  "--genomeDir {output_dir} --genomeFastaFiles {genome_fa} " \
                  "--sjdbGTFfile {gene_gtf}".format(n_threads=cores,
                      output_dir=output_dir, genome_fa=genome_fa,
                      gene_gtf=gene_gtf)
        run_stage(self.state, "star", command)

    def create_hisat_index(self, inputs, outputs, hisat_basename):
        '''Generate index for HISAT2'''
        safe_make_dir(os.path.dirname(hisat_basename))
        genome_fa, gene_gtf = inputs
        cores = self.get_stage_options("build_index", "cores")
        command = "hisat2-build -p {n_threads} {genome_fa} {basename}" \
                  "".format(n_threads=cores, genome_fa=genome_fa, gene_gtf=gene_gtf,
                      basename=hisat_basename)
        run_stage(self.state, "hisat", command)

    def star_align(self, inputs, output, ref_dir, sample):
        '''Align fastq files with STAR'''
        output_dir = os.path.dirname(output)
        safe_make_dir(output_dir)
        #logging.debug(self.experiment.tr_dict[sample])
        cores = self.get_stage_options("align", "cores")
        # If PE fastq inputs, join into a string
        if self.paired_end:
            fastq_input = " ".join(inputs)
        else:
            fastq_input = inputs
        command = "STAR --runThreadN {cores} --genomeDir {ref_dir} " \
                  "--readFilesIn {fastq_input} --readFilesCommand zcat " \
                  "--outFileNamePrefix {output_dir}/{sample}.star. " \
                  "--outSAMtype BAM Unsorted " \
                  "--outSAMunmapped Within " \
                  "".format(cores=cores, ref_dir=ref_dir, fastq_input=fastq_input, 
                          output_dir=output_dir, sample=sample)
        run_stage(self.state, 'star', command)
        # TODO: Add stranded and RG info

    def hisat_align(self, inputs, output, ref_basename, sample):
        '''Align fastq files with HISAT2 and sort'''
        cores = self.get_stage_options("hisat", "cores")
        mem = "{}G".format(self.get_stage_options("hisat", "mem"))
        safe_make_dir(os.path.dirname(output))
        #logging.debug(self.experiment.tr_dict[sample])
        output_log = re.sub(".bam$", ".log", output)
        # If PE fastq inputs, use hisat -1 and -2 arguments, else use -U
        if self.paired_end:
            fastq_input = "-1 {fastq_R1} -2 {fastq_R2}".format(
                              fastq_R1=inputs[0], fastq_R2=inputs[1])
        else:
            fastq_input = "-U {fastq}".format(fastq=inputs)
        if self.experiment.stranded == "FR":
            stranded = "--rna-strandness FR"
        elif self.experiment.stranded == "RF":
            stranded = "--rna-strandness RF"
        else:
            stranded = ""
        # Get RG information
        info = self.experiment.tr_dict[sample]
        command = "hisat2 -p {n_threads} --dta {stranded} " \
                  "--rg-id {sm}_{id}_{ln} --rg SM:{sm} " \
                  "--rg LB:{lb} --rg PL:Illumina -x {ref_basename} " \
                  "{fastq_input} 2> {output_log} | samtools view -bS - > " \
                  "{output_bam} 2>> {output_log}" \
                  "".format(n_threads=cores, stranded=stranded, id=info.id, 
                          ln=info.lane, sm=info.sample_name, lb=info.library,
                          ref_basename=ref_basename, fastq_input=fastq_input,
                          output_bam=output, output_log=output_log)
        run_stage(self.state, "hisat", command)
        # TODO: Make sure unmapped reads are included in the bam file
        #       Is it a good idea for BAMs to be cufflinks compatible?

    def sort_bam_by_coordinate(self, input, outputs):
        '''Sort BAM file by coordinates and then index'''
        output_bam = outputs[0]
        # Provide a bit of room between memory requested and samtools max memory
        mem = max(int(self.get_stage_options("samtools", "mem")) - 2, 1)
        command = "samtools sort -m {mem}G {input} > {output} && " \
                  "samtools index {output}".format(mem=mem, output=output_bam,
                      input=input)
        run_stage(self.state, "samtools", command)

    def merge_bams(self, inputs, output):
        '''Merge multiple BAM files into one BAM file. Make a symlink if
        there's only one BAM file.'''
        if len(inputs) == 1:
            re_symlink(inputs[0], output)
        else:
            # Select only bams and not bais, which is the first item in the list
            bam_inputs = " ".join([x[0] for x in inputs])
            bam_output = output[0]
            command = "samtools merge {output} {bam_inputs} && samtools index {output}".format(
                          output=bam_output, bam_inputs=bam_inputs)
            run_stage(self.state, "samtools", command)

    def sort_bam_by_name(self, input, output):
        '''Sort BAM file by name'''
        # Provide a bit of room between memory requested and samtools max memory
        mem = max(int(self.get_stage_options("samtools", "mem")) - 2, 1)
        input_bam = input[0]
        command = "samtools sort -n -m {mem}G {input} > {output}".format(
                      mem=mem, input=input_bam, output=output)
        run_stage(self.state, "samtools", command)

    def htseq_count(self, input, output):
        '''Count features with HTSeq-count'''
        # YAML converts yes/no to true/false, so we need to convert it back
        if self.experiment.stranded == "FR":
            stranded = "yes"
        elif self.experiment.stranded == "RF":
            stranded = "reverse"
        else:
            stranded = "no"
        # command = "samtools view -h -F 4 {input} | " \
        #           "htseq-count --format=sam --mode=union --order=name " \
        #           "--stranded={stranded} - {gtf_file} > {output}".format(
        #                input=input, stranded=stranded,
        #                gtf_file=self.gene_ref, output=output)
        command = "htseq-count --format=bam --mode=union --order=name " \
                  "--stranded={stranded} {input} {gtf_file} > {output}".format(
                       input=input, stranded=stranded,
                       gtf_file=self.gene_ref, output=output)
        run_stage(self.state, "htseq_count", command)

    def featurecounts(self, input, output):
        '''Count features with featureCounts'''
        create_empty_outputs(output)
        command = ""
        # run_stage(self.state, "featurecounts", command)


#     def stringtie_assembly(self, inputs, output):
#         '''Assemble transcripts with stringtie'''
#         input = inputs[0]
#         # XXX: handle SE cases?
#         if self.experiment.stranded == "FR":
#             stranded = "--fr"
#         elif self.experiment.stranded == "RF":
#             stranded = "--rf"
#         else:
#             stranded = ""
#         cores = self.get_stage_options("stringtie", "cores")
#         command = "stringtie -p {cores} {stranded} -l STRG -G {gtf} " \
#                   "-o {output} {input}".format(cores=cores,
#                       stranded=stranded, gtf=self.gene_ref,
#                       output=output, input=input)
#         run_stage(self.state, "stringtie", command)

    def stringtie_estimates(self, inputs, outputs):
        '''Get expression estimates with stringtie'''
        input = inputs[0]
        output_gtf = outputs[0]
        cores = self.get_stage_options("stringtie", "cores")
        command = "stringtie -p {cores} -e -B -G {gtf} -o {output} " \
                  "{input}".format(cores=cores, gtf=self.gene_ref,
                      output=output_gtf, input=input)
        run_stage(self.state, "stringtie", command)

    def stringtie_prepDE(self, inputs, outputs):
        '''Create count matrices'''
        input_dir = os.path.dirname(os.path.dirname(inputs[0][0]))
        command = "prepDE.py -i {input_dir} -g {gene_matrix} " \
                  "-t {transcript_matrix} -l 75".format(
                      input_dir=input_dir, gene_matrix=outputs[0],
                      transcript_matrix=outputs[1])
        run_stage(self.state, "stringtie", command)

