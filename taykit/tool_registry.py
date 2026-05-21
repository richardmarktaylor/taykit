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
                    "default": 32,
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
]
