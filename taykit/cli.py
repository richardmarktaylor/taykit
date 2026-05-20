#!/usr/bin/env python3

import argparse
import importlib

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

    for tool in TOOLS:
        subparsers.add_parser(
            tool["command"],
            help=tool["help"],
        )

    args, remaining_args = parser.parse_known_args()

    if not args.command:
        parser.print_help()
        return

    selected_tool = None

    for tool in TOOLS:
        if tool["command"] == args.command:
            selected_tool = tool
            break

    if selected_tool is None:
        parser.print_help()
        return

    module = importlib.import_module(selected_tool["module"])

    tool_parser = argparse.ArgumentParser(
        prog=f"taykit {module.COMMAND}",
        description=module.DESCRIPTION,
        epilog=module.EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    module.configure_parser(tool_parser)

    tool_args = tool_parser.parse_args(remaining_args)

    module.main(tool_args)


if __name__ == "__main__":
    main()
