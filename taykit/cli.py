#!/usr/bin/env python3

import argparse
from taykit.tools import opus


def main():
    parser = argparse.ArgumentParser(
        prog="taykit",
        description="Taylor bioinformatics command-line toolkit"
    )

    subparsers = parser.add_subparsers(dest="command")

    opus_parser = subparsers.add_parser(
        "opus",
        help="Generate OPUS reports from raw DNA files"
    )

    opus_parser.add_argument("genetic_file", nargs="?")
    opus_parser.add_argument("--output-path")

    args = parser.parse_args()

    if args.command == "opus":
        opus.main(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
