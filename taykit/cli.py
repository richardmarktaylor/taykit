#!/usr/bin/env python3

import argparse

from taykit.tools import merge, opus


def main():
    parser = argparse.ArgumentParser(
        prog="taykit",
        description="Taylor bioinformatics command-line toolkit",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<tool>",
    )

    opus_parser = subparsers.add_parser(
        "opus",
        help="Generate OPUS reports from raw DNA files",
    )
    opus_parser.add_argument("genetic_file", nargs="?")
    opus_parser.add_argument("--output-path")

    merge_parser = subparsers.add_parser(
        "merge",
        help="Merge two or more raw DNA files",
    )
    merge_parser.add_argument("genetic_files", nargs="+")
    merge_parser.add_argument("--output-path")
    merge_parser.add_argument(
        "--report-type",
        choices=["html", "json", "txt", "xml"],
    )
    merge_parser.add_argument(
        "--output-format",
        default="txt",
        choices=["txt", "tsv", "csv", "json", "jsonl", "vcf", "vcfgz", "parquet"],
    )

    args = parser.parse_args()

    if args.command == "opus":
        opus.main(args)

    elif args.command == "merge":
        merge.main(args)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
