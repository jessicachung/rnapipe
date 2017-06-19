# Default output directories
OUTPUT_PATHS = {
    "fastqc": "qc/fastqc",
    "post_trim_fastqc": "qc/post_trim_fastqc",
    "seq": "seq",
    "alignments": "alignments",
    "star_index": "reference/star",
    "hisat_index": "reference/hisat",
    "counts": "counts",
    "analysis": "analysis"
}

# Default metadata for sequencing files
DEFAULT_METADATA = {"id": "id", "lb": "lb", "r": "R1", "lane": None}

# Error codes
EXIT_INVALID_ARGUMENT = 2
EXIT_INVALID_CSV = 3
