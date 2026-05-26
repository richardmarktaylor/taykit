import os
from pathlib import Path

TOOLS = [
    {
        "command": "ancestry",
        "help": "Run 1000 Genomes 30x ancestry PCA and ADMIXTURE pipeline",
        "module": "taykit.tools.ancestry",
        "description": "Prepare a 1000 Genomes 30x GRCh38 ancestry reference and run ancestry inference for a sample file.",
        "arguments": [
            {
                "name": "--prepare",
                "kwargs": {
                    "action": "store_true",
                    "help": "Download tools/reference data and build reusable reference files, then exit.",
                },
            },
            {
                "name": "--sample-file",
                "kwargs": {
                    "type": Path,
                    "help": "Input sample file: .txt, .vcf, .gvcf, .vcf.gz, or .gvcf.gz",
                },
            },
            {
                "name": "--force",
                "kwargs": {
                    "action": "store_true",
                    "help": "Force rebuilding/redownloading where possible.",
                },
            },
            {
                "name": "--threads",
                "kwargs": {
                    "type": int,
                    "default": os.cpu_count() or 1,
                    "help": "Number of CPU threads to use.",
                },
            },
            {
                "name": "--admixture-k",
                "kwargs": {
                    "type": int,
                    "default": 5,
                    "help": "Number of ADMIXTURE ancestral clusters to model.",
                },
            },
        ],
    },
    {
        "command": "opus",
        "help": "Generate OPUS reports from raw DNA files",
        "module": "taykit.tools.opus",
        "description": "Generate OPUS reports from raw DNA files.",
        "arguments": [
            {
                "name": "genetic_file",
                "kwargs": {"nargs": "?", "help": "Path to the raw DNA file"},
            },
            {
                "name": "--output-path",
                "kwargs": {"help": "Folder to write output reports into"},
            },
        ],
    },
    {
        "command": "merge",
        "help": "Merge two or more raw DNA files",
        "module": "taykit.tools.merge",
        "description": "Merge two or more raw genetic data files into one distinct rsID/genotype dataset.",
        "arguments": [
            {
                "name": "genetic_files",
                "kwargs": {"nargs": "+", "help": "Two or more raw DNA files to merge"},
            },
            {
                "name": "--output-path",
                "kwargs": {"help": "Folder to write output files into"},
            },
            {
                "name": "--report-type",
                "kwargs": {
                    "choices": ["html", "json", "txt", "xml"],
                    "help": "Optional report format",
                },
            },
            {
                "name": "--output-format",
                "kwargs": {
                    "default": "txt",
                    "choices": [
                        "txt",
                        "tsv",
                        "csv",
                        "json",
                        "jsonl",
                        "vcf",
                        "vcfgz",
                        "parquet",
                    ],
                    "help": "Output format",
                },
            },
        ],
    },
    {
        "command": "liftover",
        "help": "Convert raw DNA files from GRCh37 to GRCh38",
        "module": "taykit.tools.liftover",
        "description": "Convert supported raw genetic files from GRCh37 coordinates to GRCh38 coordinates.",
        "arguments": [
            {
                "name": "genetic_file",
                "kwargs": {
                    "nargs": "?",
                    "help": "Raw DNA file to convert from GRCh37 to GRCh38",
                },
            },
            {
                "name": "--output-path",
                "kwargs": {"help": "Folder to write output files into"},
            },
            {
                "name": "--output-format",
                "kwargs": {
                    "default": "txt",
                    "choices": [
                        "txt",
                        "tsv",
                        "csv",
                        "json",
                        "jsonl",
                        "vcf",
                        "vcfgz",
                        "parquet",
                    ],
                    "help": "Output format",
                },
            },
            {
                "name": "--report-type",
                "kwargs": {
                    "choices": ["html", "json", "txt", "xml"],
                    "help": "Optional report format",
                },
            },
            {
                "name": "--assume-grch37",
                "kwargs": {
                    "action": "store_true",
                    "help": "Allow liftover when build cannot be detected",
                },
            },
        ],
    },
    {
        "command": "impute",
        "help": "Impute GRCh38 raw DNA files using 1000 Genomes 30x",
        "module": "taykit.tools.impute",
        "description": "Run strict GRCh38 genotype imputation using 1000 Genomes 30x, Beagle, ShapeIT4 and optional IMPUTE5.",
        "arguments": [
            {
                "name": "input_file",
                "kwargs": {
                    "nargs": "?",
                    "help": "Path to raw DNA input file .txt, .tsv, or .csv",
                },
            },
            {
                "name": "output_file",
                "kwargs": {"nargs": "?", "help": "Path to final output file"},
            },
            {
                "name": "--wizard",
                "kwargs": {
                    "action": "store_true",
                    "help": "Launch an interactive imputation command builder.",
                },
            },
            {
                "name": "--output-format",
                "kwargs": {"choices": ["tsv", "vcf", "vcf.gz"], "default": "tsv"},
            },
            {"name": "--min-quality", "kwargs": {"type": float, "default": 0.8}},
            {"name": "--min-dr2", "kwargs": {"type": float, "dest": "min_quality"}},
            {
                "name": "--impute-chr-x",
                "kwargs": {"action": "store_true", "default": False},
            },
            {
                "name": "--download-references",
                "kwargs": {
                    "action": "store_true",
                    "default": True,
                    "help": "Download missing reference files. This is on by default.",
                },
            },
            {
                "name": "--no-download-references",
                "kwargs": {
                    "action": "store_false",
                    "dest": "download_references",
                    "help": "Do not download missing reference files.",
                },
            },
            {
                "name": "--threads",
                "kwargs": {"type": int, "default": os.cpu_count() or 1},
            },
            {
                "name": "--phasing-tool",
                "kwargs": {"choices": ["shapeit4", "beagle"], "default": "shapeit4"},
            },
            {"name": "--require-shapeit4", "kwargs": {"action": "store_true"}},
            {
                "name": "--imputation-mode",
                "kwargs": {"choices": ["beagle", "impute5", "both"], "default": "both"},
            },
            {"name": "--allow-position-fallback", "kwargs": {"action": "store_true"}},
            {"name": "--keep-ambiguous-snps", "kwargs": {"action": "store_true"}},
            {"name": "--low-freq-max-maf", "kwargs": {"type": float, "default": 0.05}},
            {"name": "--beagle-memory", "kwargs": {"default": "96g"}},
            {"name": "--beagle-window-cm", "kwargs": {"default": "40.0"}},
            {"name": "--beagle-overlap-cm", "kwargs": {"default": "4.0"}},
            {"name": "--beagle-window-markers", "kwargs": {"default": "2000000"}},
            {"name": "--beagle-imp-states", "kwargs": {"default": "4000"}},
            {"name": "--beagle-phase-states", "kwargs": {"default": "800"}},
            {"name": "--beagle-iterations", "kwargs": {"default": "12"}},
            {
                "name": "--beagle-em",
                "kwargs": {"choices": ["true", "false"], "default": "true"},
            },
            {"name": "--beagle-ne", "kwargs": {"default": ""}},
            {"name": "--beagle-err", "kwargs": {"default": ""}},
            {"name": "--keep-temp", "kwargs": {"action": "store_true"}},
        ],
    },
]
