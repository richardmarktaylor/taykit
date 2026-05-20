#!/usr/bin/env python3

import argparse
import csv
import gzip
import html
import json
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

COMMAND = "merge"

HELP = "Merge two or more raw DNA files"

DESCRIPTION = (
    "Merge two or more raw genetic data files into one distinct rsID/genotype dataset."
)

EPILOG = """
Examples:
  taykit merge file1.txt file2.txt
  taykit merge file1.txt file2.txt --output-path ./merged
  taykit merge file1.txt file2.txt --report-type html
  taykit merge file1.txt file2.txt --output-format json --report-type json
  taykit merge file1.txt file2.txt file3.txt --output-format vcfgz --report-type html

What this tool does:
  - Accepts two or more supported raw DNA files.
  - Detects the source/provider format automatically.
  - Detects the genome build where possible.
  - Rejects mixed GRCh37/GRCh38 merges.
  - Creates a non-conflicting merged output.
  - Creates a separate conflicting output.
  - Optionally creates a report in HTML, JSON, TXT, or XML.

Supported input providers:
  - AncestryDNA
  - SelfDecode
  - MyHeritage
  - 23andMe-style tab-delimited raw data

Output files:
  - merged.<format>
  - conflicting.<format>
  - report.<report-type>, only if --report-type is supplied

Default behaviour:
  - Output folder: same folder as the first input file
  - Output format: txt
  - Report: none

Conflict logic:
  - If an rsID appears in multiple files with the same genotype, it is merged.
  - If an rsID appears in multiple files with different genotypes, it is written to conflicting output only.
  - Genotypes such as AG and GA are treated as equivalent.
"""


SUPPORTED_OUTPUT_FORMATS = [
    "txt",
    "tsv",
    "csv",
    "json",
    "jsonl",
    "vcf",
    "vcfgz",
    "parquet",
]

SUPPORTED_REPORT_TYPES = [
    "html",
    "json",
    "txt",
    "xml",
]


LIFTOVER_MESSAGE = """
ERROR: Input files use different genome builds.

You cannot safely merge GRCh37 / Build 37 and GRCh38 / Build 38 files directly.

Please run liftover first so all files use the same genome build, then run this merge tool again.

No merge files were created.
"""


def build_parser():
    parser = argparse.ArgumentParser(
        prog="taykit merge",
        description="Merge two or more raw genetic data files into a distinct rsID/genotype dataset.",
        epilog="""
Examples:
  taykit merge file1.txt file2.txt
  taykit merge file1.txt file2.txt --output-path ./merged
  taykit merge file1.txt file2.txt --report-type html
  taykit merge file1.txt file2.txt --output-format json --report-type json
  taykit merge file1.txt file2.txt file3.txt --output-format vcfgz --report-type html

What this tool does:
  - Accepts two or more supported raw DNA files.
  - Detects the source/provider format automatically.
  - Detects the genome build where possible.
  - Rejects mixed GRCh37/GRCh38 merges.
  - Creates a non-conflicting merged output.
  - Creates a separate conflicting output.
  - Optionally creates a report in HTML, JSON, TXT, or XML.

Supported input providers:
  - AncestryDNA
  - SelfDecode
  - MyHeritage
  - 23andMe-style tab-delimited raw data

Output files:
  - merged.<format>
  - conflicting.<format>
  - report.<report-type>, only if --report-type is supplied

Default behaviour:
  - Output folder: same folder as the first input file
  - Output format: txt
  - Report: none

Conflict logic:
  - If an rsID appears in multiple files with the same genotype, it is merged.
  - If an rsID appears in multiple files with different genotypes, it is written to conflicting output only.
  - Genotypes such as AG and GA are treated as equivalent.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "genetic_files",
        nargs="+",
        help="Two or more raw DNA files to merge",
    )

    parser.add_argument(
        "--output-path",
        help="Folder to write merged output files into. Defaults to the same folder as the first input file.",
    )

    parser.add_argument(
        "--report-type",
        choices=SUPPORTED_REPORT_TYPES,
        help="Optional report format. Choices: html, json, txt, xml. Default: no report.",
    )

    parser.add_argument(
        "--output-format",
        default="txt",
        choices=SUPPORTED_OUTPUT_FORMATS,
        help="Merged output format. Choices: txt, tsv, csv, json, jsonl, vcf, vcfgz, parquet. Default: txt.",
    )

    return parser


def parse_args():
    parser = build_parser()
    args = parser.parse_args()

    if len(args.genetic_files) < 2:
        parser.print_help()
        print()
        print("ERROR: Please provide at least two genetic files to merge.")
        sys.exit(1)

    return args


def detect_provider_and_build(file_path: Path):
    header_lines = []

    with file_path.open("r", encoding="utf-8", errors="ignore") as file:
        for _ in range(60):
            line = file.readline()
            if not line:
                break
            header_lines.append(line.strip())

    joined = "\n".join(header_lines).lower()

    if "ancestrydna raw data download" in joined:
        return {
            "provider": "AncestryDNA",
            "build": "GRCh37",
        }

    if "source: selfdecode" in joined:
        if "assembly: grch38" in joined:
            return {
                "provider": "SelfDecode",
                "build": "GRCh38",
            }

        if "assembly: grch37" in joined or "build 37" in joined:
            return {
                "provider": "SelfDecode",
                "build": "GRCh37",
            }

        return {
            "provider": "SelfDecode",
            "build": "Unknown",
        }

    if "fileformat=myheritage" in joined:
        if "reference=build37" in joined:
            return {
                "provider": "MyHeritage",
                "build": "GRCh37",
            }

        if "reference=build38" in joined:
            return {
                "provider": "MyHeritage",
                "build": "GRCh38",
            }

        return {
            "provider": "MyHeritage",
            "build": "Unknown",
        }

    if "23andme" in joined:
        if "build 37" in joined:
            return {
                "provider": "23andMe",
                "build": "GRCh37",
            }

        if "build 38" in joined:
            return {
                "provider": "23andMe",
                "build": "GRCh38",
            }

        return {
            "provider": "23andMe",
            "build": "Unknown",
        }

    # Generic raw data detection
    for line in header_lines:
        lower = line.lower()

        if lower in [
            "rsid\tchromosome\tposition\tgenotype",
            "# rsid\tchromosome\tposition\tgenotype",
        ]:
            return {
                "provider": "GenericTab",
                "build": "Unknown",
            }

        if lower == "rsid\tchromosome\tposition\tallele1\tallele2":
            return {
                "provider": "GenericTabAlleles",
                "build": "Unknown",
            }

        if lower == "rsid,chromosome,position,result":
            return {
                "provider": "GenericCsv",
                "build": "Unknown",
            }

    return {
        "provider": "Unknown",
        "build": "Unknown",
    }


def normalise_rsid(rsid: str):
    rsid = rsid.strip()

    if rsid.startswith("GSA-"):
        rsid = rsid.replace("GSA-", "", 1)

    return rsid


def normalise_genotype(genotype: str):
    genotype = genotype.strip().upper()

    if not genotype:
        return ""

    genotype = genotype.replace("/", "")
    genotype = genotype.replace("|", "")
    genotype = genotype.replace(" ", "")

    if genotype in ["--", "00", "NN", "NOCALL", "NO_CALL"]:
        return ""

    if len(genotype) == 2 and all(base in "ACGT" for base in genotype):
        return "".join(sorted(genotype))

    return genotype


def load_genetic_data(file_path: Path, provider: str):
    records = {}

    if provider in ["SelfDecode", "23andMe", "GenericTab"]:
        with file_path.open("r", encoding="utf-8", errors="ignore") as file:
            for line in file:
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")

                if len(parts) < 4:
                    continue

                if parts[0].lower() == "rsid":
                    continue

                rsid = normalise_rsid(parts[0])
                chromosome = parts[1].strip()
                position = parts[2].strip()
                genotype = normalise_genotype(parts[3])

                if not rsid or not genotype:
                    continue

                records[rsid] = {
                    "rsid": rsid,
                    "chromosome": chromosome,
                    "position": position,
                    "genotype": genotype,
                }

    elif provider in ["AncestryDNA", "GenericTabAlleles"]:
        with file_path.open("r", encoding="utf-8", errors="ignore") as file:
            for line in file:
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")

                if len(parts) < 5:
                    continue

                if parts[0].lower() == "rsid":
                    continue

                rsid = normalise_rsid(parts[0])
                chromosome = parts[1].strip()
                position = parts[2].strip()
                genotype = normalise_genotype(parts[3] + parts[4])

                if not rsid or not genotype:
                    continue

                records[rsid] = {
                    "rsid": rsid,
                    "chromosome": chromosome,
                    "position": position,
                    "genotype": genotype,
                }

    elif provider in ["MyHeritage", "GenericCsv"]:
        with file_path.open("r", encoding="utf-8", errors="ignore", newline="") as file:
            reader = csv.reader(file)

            for row in reader:
                if not row:
                    continue

                if row[0].startswith("#"):
                    continue

                if row[0].upper() == "RSID":
                    continue

                if len(row) < 4:
                    continue

                rsid = normalise_rsid(row[0].replace('"', "").strip())
                chromosome = row[1].replace('"', "").strip()
                position = row[2].replace('"', "").strip()
                genotype = normalise_genotype(row[3].replace('"', "").strip())

                if not rsid or not genotype:
                    continue

                records[rsid] = {
                    "rsid": rsid,
                    "chromosome": chromosome,
                    "position": position,
                    "genotype": genotype,
                }

    return records


def validate_inputs(input_files):
    file_infos = []

    for file_path in input_files:
        info = detect_provider_and_build(file_path)
        file_infos.append(
            {
                "path": file_path,
                "provider": info["provider"],
                "build": info["build"],
            }
        )

    known_builds = {
        file_info["build"]
        for file_info in file_infos
        if file_info["build"] in ["GRCh37", "GRCh38"]
    }

    if len(known_builds) > 1:
        print(LIFTOVER_MESSAGE)
        sys.exit(1)

    unknown_providers = [
        str(file_info["path"])
        for file_info in file_infos
        if file_info["provider"] == "Unknown"
    ]

    if unknown_providers:
        print("ERROR: One or more files could not be recognised:")
        for path in unknown_providers:
            print(f"  - {path}")
        print()
        print("No merge files were created.")
        sys.exit(1)

    return file_infos


def merge_records(file_records):
    observed = defaultdict(list)

    for file_record in file_records:
        source_name = file_record["path"].name

        for rsid, record in file_record["records"].items():
            observed[rsid].append(
                {
                    "source": source_name,
                    "chromosome": record["chromosome"],
                    "position": record["position"],
                    "genotype": record["genotype"],
                }
            )

    merged = []
    conflicting = []

    for rsid in sorted(observed.keys()):
        entries = observed[rsid]
        genotypes = sorted({entry["genotype"] for entry in entries})

        chromosome = entries[0]["chromosome"]
        position = entries[0]["position"]

        if len(genotypes) == 1:
            merged.append(
                {
                    "rsid": rsid,
                    "chromosome": chromosome,
                    "position": position,
                    "genotype": genotypes[0],
                    "source_count": len(entries),
                    "sources": ", ".join(
                        sorted({entry["source"] for entry in entries})
                    ),
                }
            )
        else:
            conflicting.append(
                {
                    "rsid": rsid,
                    "chromosome": chromosome,
                    "position": position,
                    "genotypes": "; ".join(
                        f"{entry['source']}={entry['genotype']}"
                        for entry in sorted(entries, key=lambda item: item["source"])
                    ),
                    "source_count": len(entries),
                    "sources": ", ".join(
                        sorted({entry["source"] for entry in entries})
                    ),
                }
            )

    return merged, conflicting


def output_extension(output_format):
    if output_format == "vcfgz":
        return "vcf.gz"

    return output_format


def write_delimited(output_file: Path, rows, delimiter):
    with output_file.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=list(rows[0].keys())
            if rows
            else ["rsid", "chromosome", "position", "genotype"],
            delimiter=delimiter,
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(output_file: Path, rows):
    with output_file.open("w", encoding="utf-8") as file:
        json.dump(rows, file, indent=2)


def write_jsonl(output_file: Path, rows):
    with output_file.open("w", encoding="utf-8") as file:
        for row in rows:
            file.write(json.dumps(row) + "\n")


def write_vcf(output_file: Path, rows, gzip_output=False):
    opener = gzip.open if gzip_output else open
    mode = "wt" if gzip_output else "w"

    with opener(output_file, mode, encoding="utf-8") as file:
        file.write("##fileformat=VCFv4.2\n")
        file.write("##source=taykit-merge\n")
        file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for row in rows:
            chromosome = row.get("chromosome", ".") or "."
            position = row.get("position", ".") or "."
            rsid = row.get("rsid", ".") or "."
            genotype = row.get("genotype", ".") or "."

            file.write(
                f"{chromosome}\t"
                f"{position}\t"
                f"{rsid}\t"
                f"N\t"
                f"<GENOTYPE>\t"
                f".\t"
                f"PASS\t"
                f"GT_STRING={genotype}\t"
                f"GT\t"
                f"{genotype}\n"
            )


def write_parquet(output_file: Path, rows):
    try:
        import pandas as pd
    except ImportError:
        print("ERROR: Parquet output requires pandas and pyarrow.")
        print("Install them with:")
        print("  python -m pip install pandas pyarrow")
        sys.exit(1)

    dataframe = pd.DataFrame(rows)
    dataframe.to_parquet(output_file, index=False)


def write_output_files(output_dir: Path, merged, conflicting, output_format):
    extension = output_extension(output_format)

    merged_file = output_dir / f"merged.{extension}"
    conflicting_file = output_dir / f"conflicting.{extension}"

    if output_format == "txt":
        write_delimited(merged_file, merged, "\t")
        write_delimited(conflicting_file, conflicting, "\t")

    elif output_format == "tsv":
        write_delimited(merged_file, merged, "\t")
        write_delimited(conflicting_file, conflicting, "\t")

    elif output_format == "csv":
        write_delimited(merged_file, merged, ",")
        write_delimited(conflicting_file, conflicting, ",")

    elif output_format == "json":
        write_json(merged_file, merged)
        write_json(conflicting_file, conflicting)

    elif output_format == "jsonl":
        write_jsonl(merged_file, merged)
        write_jsonl(conflicting_file, conflicting)

    elif output_format == "vcf":
        write_vcf(merged_file, merged, gzip_output=False)
        write_vcf(conflicting_file, conflicting, gzip_output=False)

    elif output_format == "vcfgz":
        write_vcf(merged_file, merged, gzip_output=True)
        write_vcf(conflicting_file, conflicting, gzip_output=True)

    elif output_format == "parquet":
        write_parquet(merged_file, merged)
        write_parquet(conflicting_file, conflicting)

    return merged_file, conflicting_file


def build_summary(file_records, merged, conflicting, merged_file, conflicting_file):
    total_input_snps = sum(len(file_record["records"]) for file_record in file_records)

    unique_input_rsids = set()

    for file_record in file_records:
        unique_input_rsids.update(file_record["records"].keys())

    return {
        "input_files": [
            {
                "file": file_record["path"].name,
                "provider": file_record["provider"],
                "build": file_record["build"],
                "snp_count": len(file_record["records"]),
            }
            for file_record in file_records
        ],
        "total_input_snp_rows": total_input_snps,
        "unique_input_rsids": len(unique_input_rsids),
        "merged_non_conflicting_count": len(merged),
        "conflicting_count": len(conflicting),
        "merged_output": str(merged_file),
        "conflicting_output": str(conflicting_file),
    }


def write_html_report(output_file: Path, summary, conflicting):
    input_rows = ""

    for item in summary["input_files"]:
        input_rows += f"""
        <tr>
            <td>{html.escape(item["file"])}</td>
            <td>{html.escape(item["provider"])}</td>
            <td>{html.escape(item["build"])}</td>
            <td>{item["snp_count"]:,}</td>
        </tr>
        """

    conflict_rows = ""

    for item in conflicting:
        conflict_rows += f"""
        <tr>
            <td>{html.escape(item["rsid"])}</td>
            <td>{html.escape(item["chromosome"])}</td>
            <td>{html.escape(item["position"])}</td>
            <td>{html.escape(item["genotypes"])}</td>
            <td>{html.escape(item["sources"])}</td>
        </tr>
        """

    if not conflict_rows:
        conflict_rows = """
        <tr>
            <td colspan="5" class="text-success">No genotype conflicts found.</td>
        </tr>
        """

    html_output = f"""
<!doctype html>
<html>
<head>
    <meta charset="utf-8">
    <title>taykit merge report</title>
    <link
        href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.7/dist/css/bootstrap.min.css"
        rel="stylesheet"
    >
</head>
<body class="bg-light">

<div class="container py-5">

    <h1 class="mb-4">taykit merge report</h1>

    <div class="row g-3 mb-4">
        <div class="col-md-3">
            <div class="card shadow-sm">
                <div class="card-body">
                    <h6>Total input SNP rows</h6>
                    <div class="display-6">{summary["total_input_snp_rows"]:,}</div>
                </div>
            </div>
        </div>

        <div class="col-md-3">
            <div class="card shadow-sm">
                <div class="card-body">
                    <h6>Unique input rsIDs</h6>
                    <div class="display-6">{summary["unique_input_rsids"]:,}</div>
                </div>
            </div>
        </div>

        <div class="col-md-3">
            <div class="card shadow-sm">
                <div class="card-body">
                    <h6>Merged non-conflicting</h6>
                    <div class="display-6">{summary["merged_non_conflicting_count"]:,}</div>
                </div>
            </div>
        </div>

        <div class="col-md-3">
            <div class="card shadow-sm">
                <div class="card-body">
                    <h6>Conflicting rsIDs</h6>
                    <div class="display-6">{summary["conflicting_count"]:,}</div>
                </div>
            </div>
        </div>
    </div>

    <h2>Input files</h2>

    <table class="table table-bordered table-striped bg-white">
        <thead class="table-dark">
            <tr>
                <th>File</th>
                <th>Provider</th>
                <th>Build</th>
                <th>SNP count</th>
            </tr>
        </thead>
        <tbody>
            {input_rows}
        </tbody>
    </table>

    <h2 class="mt-5">Output files</h2>

    <div class="alert alert-primary">
        <strong>Merged:</strong> {html.escape(summary["merged_output"])}<br>
        <strong>Conflicting:</strong> {html.escape(summary["conflicting_output"])}
    </div>

    <h2 class="mt-5">Conflicting genotypes</h2>

    <table class="table table-bordered table-striped bg-white">
        <thead class="table-dark">
            <tr>
                <th>rsID</th>
                <th>Chromosome</th>
                <th>Position</th>
                <th>Genotypes</th>
                <th>Sources</th>
            </tr>
        </thead>
        <tbody>
            {conflict_rows}
        </tbody>
    </table>

</div>

</body>
</html>
"""

    with output_file.open("w", encoding="utf-8") as file:
        file.write(html_output)


def write_txt_report(output_file: Path, summary, conflicting):
    with output_file.open("w", encoding="utf-8") as file:
        file.write("taykit merge report\n")
        file.write("=" * 80 + "\n\n")

        file.write("Input files\n")
        file.write("-" * 80 + "\n")

        for item in summary["input_files"]:
            file.write(
                f"{item['file']} | "
                f"{item['provider']} | "
                f"{item['build']} | "
                f"{item['snp_count']:,} SNPs\n"
            )

        file.write("\nSummary\n")
        file.write("-" * 80 + "\n")
        file.write(f"Total input SNP rows: {summary['total_input_snp_rows']:,}\n")
        file.write(f"Unique input rsIDs: {summary['unique_input_rsids']:,}\n")
        file.write(
            f"Merged non-conflicting count: {summary['merged_non_conflicting_count']:,}\n"
        )
        file.write(f"Conflicting count: {summary['conflicting_count']:,}\n")
        file.write(f"Merged output: {summary['merged_output']}\n")
        file.write(f"Conflicting output: {summary['conflicting_output']}\n")

        file.write("\nConflicting genotypes\n")
        file.write("-" * 80 + "\n")

        if not conflicting:
            file.write("No genotype conflicts found.\n")
        else:
            for item in conflicting:
                file.write(
                    f"{item['rsid']} | "
                    f"{item['chromosome']} | "
                    f"{item['position']} | "
                    f"{item['genotypes']}\n"
                )


def write_xml_report(output_file: Path, summary, conflicting):
    root = ET.Element("taykit_merge_report")

    input_files = ET.SubElement(root, "input_files")

    for item in summary["input_files"]:
        input_file = ET.SubElement(input_files, "input_file")
        for key, value in item.items():
            child = ET.SubElement(input_file, key)
            child.text = str(value)

    summary_element = ET.SubElement(root, "summary")

    for key, value in summary.items():
        if key == "input_files":
            continue

        child = ET.SubElement(summary_element, key)
        child.text = str(value)

    conflicts = ET.SubElement(root, "conflicting_genotypes")

    for item in conflicting:
        conflict = ET.SubElement(conflicts, "conflict")
        for key, value in item.items():
            child = ET.SubElement(conflict, key)
            child.text = str(value)

    tree = ET.ElementTree(root)
    tree.write(output_file, encoding="utf-8", xml_declaration=True)


def write_report(output_dir: Path, report_type, summary, conflicting):
    if not report_type:
        return None

    report_file = output_dir / f"report.{report_type}"

    if report_type == "html":
        write_html_report(report_file, summary, conflicting)

    elif report_type == "json":
        with report_file.open("w", encoding="utf-8") as file:
            json.dump(
                {
                    "summary": summary,
                    "conflicting": conflicting,
                },
                file,
                indent=2,
            )

    elif report_type == "txt":
        write_txt_report(report_file, summary, conflicting)

    elif report_type == "xml":
        write_xml_report(report_file, summary, conflicting)

    return report_file


def run_merge(args):
    input_files = [Path(path).expanduser().resolve() for path in args.genetic_files]

    for input_file in input_files:
        if not input_file.exists():
            print(f"ERROR: File not found: {input_file}")
            sys.exit(1)

    file_infos = validate_inputs(input_files)

    file_records = []

    for file_info in file_infos:
        records = load_genetic_data(file_info["path"], file_info["provider"])

        file_records.append(
            {
                "path": file_info["path"],
                "provider": file_info["provider"],
                "build": file_info["build"],
                "records": records,
            }
        )

        print(
            f"Loaded {len(records):,} SNPs from "
            f"{file_info['path'].name} "
            f"({file_info['provider']}, {file_info['build']})"
        )

    if args.output_path:
        output_dir = Path(args.output_path).expanduser().resolve()
    else:
        output_dir = input_files[0].parent

    output_dir.mkdir(parents=True, exist_ok=True)

    merged, conflicting = merge_records(file_records)

    merged_file, conflicting_file = write_output_files(
        output_dir,
        merged,
        conflicting,
        args.output_format,
    )

    summary = build_summary(
        file_records,
        merged,
        conflicting,
        merged_file,
        conflicting_file,
    )

    report_file = write_report(
        output_dir,
        args.report_type,
        summary,
        conflicting,
    )

    print()
    print(f"Merged non-conflicting SNPs : {len(merged):,}")
    print(f"Conflicting SNPs            : {len(conflicting):,}")
    print(f"Merged output               : {merged_file}")
    print(f"Conflicting output          : {conflicting_file}")

    if report_file:
        print(f"Report                      : {report_file}")
    else:
        print(
            "Report                      : not created. Use --report-type html/json/txt/xml to create one."
        )


def configure_parser(parser):
    parser.add_argument(
        "genetic_files",
        nargs="+",
        help="Two or more raw DNA files to merge",
    )

    parser.add_argument(
        "--output-path",
        help="Folder to write merged output files into. Defaults to the same folder as the first input file.",
    )

    parser.add_argument(
        "--report-type",
        choices=SUPPORTED_REPORT_TYPES,
        help="Optional report format. Choices: html, json, txt, xml. Default: no report.",
    )

    parser.add_argument(
        "--output-format",
        default="txt",
        choices=SUPPORTED_OUTPUT_FORMATS,
        help="Merged output format. Choices: txt, tsv, csv, json, jsonl, vcf, vcfgz, parquet. Default: txt.",
    )


def main(args=None):
    if args is None:
        args = parse_args()

    if not getattr(args, "genetic_files", None) or len(args.genetic_files) < 2:
        parser = build_parser()
        parser.print_help()
        print()
        print("ERROR: Please provide at least two genetic files to merge.")
        sys.exit(1)

    run_merge(args)


if __name__ == "__main__":
    main()
