#!/usr/bin/env python3

import argparse
import importlib
import sys

from taykit.tool_registry import TOOLS


def _pyinstaller_hidden_imports():
    """
    This function is never called at runtime.

    It exists so PyInstaller can see and bundle the dynamically loaded tool
    modules when building the one-file executable.
    """

    import taykit.tools.ancestry  # noqa: F401
    import taykit.tools.impute  # noqa: F401
    import taykit.tools.liftover  # noqa: F401
    import taykit.tools.merge  # noqa: F401
    import taykit.tools.opus  # noqa: F401


def find_tool(command):
    for tool in TOOLS:
        if tool["command"] == command:
            return tool

    return None


def build_tool_parser(tool):
    parser = argparse.ArgumentParser(
        prog=f"taykit {tool['command']}",
        description=tool["description"],
        epilog=tool.get("epilog"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    for argument in tool.get("arguments", []):
        parser.add_argument(argument["name"], **argument["kwargs"])

    return parser


def print_main_help():
    print("Taylor bioinformatics command-line toolkit")
    print()
    print("Usage:")
    print("  taykit <tool> <input_file> [output_file] [options]")
    print()
    print("Available tools:")

    for tool in TOOLS:
        print(f"  {tool['command']:<12} {tool['help']}")

    print()
    print("Examples:")
    print("  taykit opus input.txt")
    print("  taykit merge input1.txt input2.txt")
    print("  taykit liftover input.txt output.txt")
    print("  taykit ancestry input.txt")
    print("  taykit impute input.txt output.vcf.gz --output-format vcf.gz")
    print("  taykit impute --wizard")
    print()
    print("Help:")
    print("  taykit <tool> --help")


def main():
    if len(sys.argv) == 1 or sys.argv[1] in ["-h", "--help"]:
        print_main_help()
        return

    command = sys.argv[1]
    remaining_args = sys.argv[2:]

    tool = find_tool(command)

    if tool is None:
        print(f"ERROR: Unknown taykit tool: {command}")
        print()
        print_main_help()
        sys.exit(1)

    tool_parser = build_tool_parser(tool)

    if not remaining_args or remaining_args[0] in ["-h", "--help"]:
        tool_parser.print_help()
        return

    tool_args = tool_parser.parse_args(remaining_args)

    module = importlib.import_module(tool["module"])
    module.main(tool_args)


if __name__ == "__main__":
    main()
