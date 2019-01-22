# Default output directories
OUTPUT_PATHS = {
    "fastqc": "qc/fastqc",
    "post_trim_fastqc": "qc/post_trim_fastqc",
    "seq": "seq",
    "alignments": "alignments",
    "star_index": "reference/star",
    "hisat_index": "reference/hisat",
    "counts": "counts",
    #"stringtie_assembly": "stringtie_assemblies",
    "stringtie_estimates": "stringtie_estimates",
    "analysis": "analysis"
}

# Default metadata for sequencing files
DEFAULT_METADATA = {"id": "id", "lb": "library", "r": "R1", "lane": "L000"}

# Error codes
EXIT_INVALID_ARGUMENT = 2
EXIT_INVALID_CSV = 3
