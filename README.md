# rnapipe

A bioinformatics pipeline for RNA-seq analysis with support for
running pipeline stages on a distributed compute cluster.

rnapipe is based on the [Ruffus](http://www.ruffus.org.uk/index.html)
pipeline library.

# Dependencies

If you are working on the university's HPC cluster, Spartan, most of the
dependencies are available using the module system.

- [Python 3](https://www.python.org/downloads/)
- [DRMAA](http://apps.man.poznan.pl/trac/slurm-drmaa)
- [Samtools](http://www.htslib.org/download/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](http://multiqc.info/)
- If trimming reads:
  - [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- For alignment:
  - [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml), or
  - [STAR](https://github.com/alexdobin/STAR)
- For quantification:
  - [HTSeq](https://github.com/simon-anders/htseq), or
  - [featureCounts](http://subread.sourceforge.net/), or
  - [StringTie](https://github.com/gpertea/stringtie)

# Installation

If you are working on Spartan, you can install rnapipe with the
following commands.

Create a Python virtual environment for rnapipe and activate it:
```
module load Python/3.6.4-intel-2018.u4
python3 -m venv rnapipe
source rnapipe/bin/activate
```

Install Ruffus and rnapipe:
```
pip install -U git+https://github.com/jessicachung/ruffus
pip install -U git+https://github.com/jessicachung/rnapipe
```

Install any dependencies. You can either install them normally or in a conda
environment. For example:
```
module load Anaconda3/5.3.1
conda create -n rna_seq
source activate rna_seq
conda install -c bioconda hisat2==2.0.5
```


# Usage

It is recommended you run the pipeline in a terminal multiplexer such as
[Screen](https://en.wikipedia.org/wiki/GNU_Screen) or
[Tmux](https://en.wikipedia.org/wiki/Tmux).

### Set DRMAA path

If you're using the pipeline on a cluster with a job scheduling system, you'll
need to tell the cluster where your DRMAA library is with `DRMAA_LIBRARY_PATH`.

```
# On Spartan
DRMAA_LIBRARY_PATH=/usr/local/slurm_drmaa/1.0.7-GCC-4.9.2/lib/libdrmaa.so
export DRMAA_LIBRARY_PATH
```

### Edit filenames of FASTQ files

rnapipe requires the input FASTQ files to be named in a certain format. It is
recommended you create symlinks of your FASTQ files with new names and keep the
original files unchanged.

The filenames should be in the format:
```
SM_xxx_ID_xxx_LB_xxx_Lxxx_Rx.fastq.gz
```
where:
- SM: Sample name
- ID: ID
- LB: Library
- L: Lane
- R: 1 (for forward) or 2 (for reverse)

Note: technical replicates must have the same sample name (SM) and will
be merged together after alignment.

Filenames like:
```
xxx_R1.fastq.gz
```
will also work, but no metadata will be included in the alignments. If your
data is single-end, `R1` is still needed in the filenames.

### Configuration files

Running the pipeline requires a configuration file. You can find an example
config file in [example/pipeline.config](example/pipeline.config).

Create a file called `samples.csv` containing a list of your sample names (SM)
and group. For example:
```
sample1,case
sample2,case
sample3,case
sample4,control
sample5,control
sample6,control
```

### Dry run

Print out steps the pipeline will run with `-n`. Specify the target task the
pipeline will complete with `--target_tasks`. Valid target tasks are listed
in the Stages section below.

```
rnapipe \
    --config pipeline.config \
    --verbose 3 \
    --target_tasks multiqc_fastqc \
    -n
```

### Running the pipeline

Run the pipeline, allowing multiple jobs to be submitted to the queue with
`--jobs` and `--use_threads`.

```
rnapipe \
    --config pipeline.config \
    --verbose 3 \
    --log_file rnapipe.log \
    --target_tasks multiqc_fastqc \
    --jobs 4 \
    --use_threads
```

### Help message

```
rnapipe --help
```

# Stages

```
create_star_index
create_hisat_index
fastqc
trim_reads
post_trim_fastqc
star_align
hisat_align
sort_bam_by_coordinate
merge_bams
sort_bam_by_name
htseq_count
featurecounts
stringtie_estimates
stringtie_prepDE
```

