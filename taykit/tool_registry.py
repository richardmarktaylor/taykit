from taykit.tools import liftover, merge, opus

TOOLS = [
    {
        "command": "opus",
        "help": "Generate OPUS reports from raw DNA files",
        "module": "taykit.tools.opus",
    },
    {
        "command": "merge",
        "help": "Merge two or more raw DNA files",
        "module": "taykit.tools.merge",
    },
    {
        "command": "liftover",
        "help": "Convert raw DNA files from GRCh37 to GRCh38",
        "module": "taykit.tools.liftover",
    },
]
