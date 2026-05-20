#!/usr/bin/env python3

import argparse
from taykit.tool_registry import TOOLS


def main():
    parser = argparse.ArgumentParser(
        prog="taykit",
        description="Taylor bioinformatics command-line toolkit",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<tool>",
    )

    tool_handlers = {}

    for tool in TOOLS:
        tool_parser = subparsers.add_parser(
            tool.COMMAND,
            help=tool.HELP,
            description=tool.DESCRIPTION,
            epilog=tool.EPILOG,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )

        tool.configure_parser(tool_parser)
        tool_handlers[tool.COMMAND] = tool.main

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return

    tool_handlers[args.command](args)


if __name__ == "__main__":
    main()
