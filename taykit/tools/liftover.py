#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import html
import json
import sys
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, Tuple


COMMAND = "liftover"
HELP = "Convert raw DNA files from GRCh37 / Build 37 to GRCh38 / Build 38"
DESCRIPTION = "Convert supported raw genetic files from GRCh37 coordinates to GRCh38 coordinates."

EPILOG = """
Examples:
  taykit liftover ancestry.txt
  taykit liftover ancestry.txt --output-path ./converted
  taykit liftover ancestry.txt --output-format csv
  taykit liftover ancestry.txt --output-format parquet --report-type html

Supported input providers:
  - AncestryDNA
  - SelfDecode
  - MyHeritage
  - 23andMe-style tab-delimited raw data

Output files:
  - <input_name>_GRCh38.<format>
  - <input_name>_failed_liftover.<format>
  - report.<report-type>, only if --report-type is supplied

Default behaviour:
  - Output folder: same folder as input file
  - Output format: txt
  - Report: none

Notes:
  - This tool only performs GRCh37 to GRCh38 liftover.
  - Files already detected as GRCh38 are rejected.
  - Unknown-build files are rejected unless --assume-grch37 is supplied.
  - VCF output is a simplified VCF-like representation because consumer DNA files do not contain full REF/ALT context.
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

CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_FILENAME = "hg19ToHg38.over.chain.gz"

TAYKIT_CACHE_DIR = Path.home() / ".taykit"
CHAIN_PATH = TAYKIT_CACHE_DIR / CHAIN_FILENAME


def configure_parser(parser):
    parser.add_argument(
        "genetic_file",
        nargs="?",
        help="Raw DNA file to convert from GRCh37 to GRCh38",
    )

    parser.add_argument(
        "--output-path",
        help="Folder to write output files into. Defaults to the same folder as the input file.",
    )

    parser.add_argument(
        "--output-format",
        default="txt",
        choices=SUPPORTED_OUTPUT_FORMATS,
        help="Output format. Choices: txt, tsv, csv, json, jsonl, vcf, vcfgz, parquet. Default: txt.",
    )

    parser.add_argument(
        "--report-type",
        choices=SUPPORTED_REPORT_TYPES,
        help="Optional report format. Choices: html, json, txt, xml. Default: no report.",
    )

    parser.add_argument(
        "--assume-grch37",
        action="store_true",
        help="Allow liftover when the file format is recognised but the build cannot be detected.",
    )


def detect_provider_and_build(file_path: Path):
    header_lines = []

    with file_path.open("r", encoding="utf-8", errors="ignore") as file:
        for _ in range(80):
            line = file.readline()
            if not line:
                break
            header_lines.append(line.strip())

    joined = "\n".join(header_lines).lower()

    if "ancestrydna raw data download" in joined:
        return {"provider": "AncestryDNA", "build": "GRCh37"}

    if "source: selfdecode" in joined:
        if "assembly: grch38" in joined:
            return {"provider": "SelfDecode", "build": "GRCh38"}
        if "assembly: grch37" in joined or "build 37" in joined:
            return {"provider": "SelfDecode", "build": "GRCh37"}
        return {"provider": "SelfDecode", "build": "Unknown"}

    if "fileformat=myheritage" in joined:
        if "reference=build37" in joined:
            return {"provider": "MyHeritage", "build": "GRCh37"}
        if "reference=build38" in joined:
            return {"provider": "MyHeritage", "build": "GRCh38"}
        return {"provider": "MyHeritage", "build": "Unknown"}

    if "23andme" in joined:
        if "build 37" in joined:
            return {"provider": "23andMe", "build": "GRCh37"}
        if "build 38" in joined:
            return {"provider": "23andMe", "build": "GRCh38"}
        return {"provider": "23andMe", "build": "Unknown"}

    for line in header_lines:
        lower = line.lower()

        if lower in [
            "rsid\tchromosome\tposition\tgenotype",
            "# rsid\tchromosome\tposition\tgenotype",
        ]:
            return {"provider": "GenericTab", "build": "Unknown"}

        if lower == "rsid\tchromosome\tposition\tallele1\tallele2":
            return {"provider": "GenericTabAlleles", "build": "Unknown"}

        if lower == "rsid,chromosome,position,result":
            return {"provider": "GenericCsv", "build": "Unknown"}

    return {"provider": "Unknown", "build": "Unknown"}


def normalise_rsid(rsid: str):
    rsid = rsid.strip()

    if rsid.startswith("GSA-"):
        rsid = rsid.replace("GSA-", "", 1)

    return rsid


def normalise_genotype(genotype: str):
    genotype = genotype.strip().upper()
    genotype = genotype.replace("/", "")
    genotype = genotype.replace("|", "")
    genotype = genotype.replace(" ", "")

    if genotype in ["--", "00", "NN", "NOCALL", "NO_CALL"]:
        return ""

    return genotype


def chromosome_to_ucsc(chromosome: str) -> Optional[str]:
    chromosome = chromosome.strip().upper()

    if chromosome.startswith("CHR"):
        chromosome = chromosome[3:]

    if chromosome in {str(i) for i in range(1, 23)}:
        return f"chr{chromosome}"

    if chromosome in ["X", "23", "25"]:
        return "chrX"

    if chromosome in ["Y", "24"]:
        return "chrY"

    if chromosome in ["MT", "M", "26"]:
        return "chrM"

    return None


def ucsc_to_plain_chromosome(chromosome: str) -> Optional[str]:
    chromosome = chromosome.strip()

    if chromosome.startswith("chr"):
        chromosome = chromosome[3:]

    chromosome = chromosome.upper()

    if chromosome in {str(i) for i in range(1, 23)}:
        return chromosome

    if chromosome == "X":
        return "X"

    if chromosome == "Y":
        return "Y"

    if chromosome in ["M", "MT"]:
        return "MT"

    return None


def download_file(url: str, destination: Path):
    destination.parent.mkdir(parents=True, exist_ok=True)

    print(f"Downloading UCSC liftover chain file:")
    print(f"  {url}")

    with urllib.request.urlopen(url) as response, destination.open("wb") as output_file:
        output_file.write(response.read())

    print(f"Saved chain file:")
    print(f"  {destination}")


def verify_gzip_file(path: Path):
    try:
        with gzip.open(path, "rb") as file:
            file.read(64)
        return True
    except Exception:
        return False


def ensure_chain_file():
    if CHAIN_PATH.exists() and verify_gzip_file(CHAIN_PATH):
        return CHAIN_PATH

    if CHAIN_PATH.exists():
        CHAIN_PATH.unlink()

    download_file(CHAIN_URL, CHAIN_PATH)

    if not verify_gzip_file(CHAIN_PATH):
        raise RuntimeError(f"Downloaded chain file is not valid: {CHAIN_PATH}")

    return CHAIN_PATH


def get_liftover():
    try:
        from pyliftover import LiftOver
    except ImportError:
        print("ERROR: pyliftover is not installed.")
        print("Rebuild taykit after adding pyliftover to build.sh:")
        print("  python -m pip install pyliftover")
        sys.exit(1)

    chain_file = ensure_chain_file()
    return LiftOver(str(chain_file))


def lift_position(liftover, chromosome: str, position_1based: int) -> Optional[Tuple[str, int]]:
    ucsc_chromosome = chromosome_to_ucsc(chromosome)

    if not ucsc_chromosome:
        return None

    if position_1based <= 0:
        return None

    position_0based = position_1based - 1

    try:
        results = liftover.convert_coordinate(ucsc_chromosome, position_0based)
    except Exception:
        return None

    if not results:
        return None

    new_chromosome, new_position_0based, _strand, _score = results[0]

    plain_chromosome = ucsc_to_plain_chromosome(new_chromosome)

    if not plain_chromosome:
        return None

    return plain_chromosome, int(new_position_0based) + 1


def load_records(file_path: Path, provider: str):
    records = []

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

                records.append({
                    "rsid": normalise_rsid(parts[0]),
                    "chromosome": parts[1].strip(),
                    "position": parts[2].strip(),
                    "genotype": normalise_genotype(parts[3]),
                })

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

                records.append({
                    "rsid": normalise_rsid(parts[0]),
                    "chromosome": parts[1].strip(),
                    "position": parts[2].strip(),
                    "genotype": normalise_genotype(parts[3] + parts[4]),
                })

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

                records.append({
                    "rsid": normalise_rsid(row[0].replace('"', "").strip()),
                    "chromosome": row[1].replace('"', "").strip(),
                    "position": row[2].replace('"', "").strip(),
                    "genotype": normalise_genotype(row[3].replace('"', "").strip()),
                })

    return [
        record for record in records
        if record["rsid"] and record["chromosome"] and record["position"] and record["genotype"]
    ]


def output_extension(output_format):
    if output_format == "vcfgz":
        return "vcf.gz"

    return output_format


def write_delimited(output_file: Path, rows, delimiter):
    fieldnames = ["rsid", "chromosome", "position", "genotype"]

    with output_file.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()

        for row in rows:
            writer.writerow({
                "rsid": row.get("rsid", ""),
                "chromosome": row.get("chromosome", ""),
                "position": row.get("position", ""),
                "genotype": row.get("genotype", ""),
            })


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
        file.write("##source=taykit-liftover\n")
        file.write("##reference=GRCh38\n")
        file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for row in rows:
            file.write(
                f"{row.get('chromosome', '.')}\t"
                f"{row.get('position', '.')}\t"
                f"{row.get('rsid', '.')}\t"
                f"N\t"
                f"<GENOTYPE>\t"
                f".\t"
                f"PASS\t"
                f"GT_STRING={row.get('genotype', '.')}\t"
                f"GT\t"
                f"{row.get('genotype', '.')}\n"
            )


def write_parquet(output_file: Path, rows):
    try:
        import pandas as pd
    except ImportError:
        print("ERROR: Parquet output requires pandas and pyarrow.")
        sys.exit(1)

    dataframe = pd.DataFrame(rows)
    dataframe.to_parquet(output_file, index=False)


def write_rows(output_file: Path, rows, output_format):
    if output_format in ["txt", "tsv"]:
        write_delimited(output_file, rows, "\t")
    elif output_format == "csv":
        write_delimited(output_file, rows, ",")
    elif output_format == "json":
        write_json(output_file, rows)
    elif output_format == "jsonl":
        write_jsonl(output_file, rows)
    elif output_format == "vcf":
        write_vcf(output_file, rows, gzip_output=False)
    elif output_format == "vcfgz":
        write_vcf(output_file, rows, gzip_output=True)
    elif output_format == "parquet":
        write_parquet(output_file, rows)


def write_html_report(report_file: Path, summary, failed_rows):
    failed_table_rows = ""

    for row in failed_rows[:1000]:
        failed_table_rows += f"""
        <tr>
            <td>{html.escape(row.get("reason", ""))}</td>
            <td>{html.escape(row.get("rsid", ""))}</td>
            <td>{html.escape(row.get("chromosome", ""))}</td>
            <td>{html.escape(str(row.get("position", "")))}</td>
            <td>{html.escape(row.get("genotype", ""))}</td>
        </tr>
        """

    if not failed_table_rows:
        failed_table_rows = """
        <tr>
            <td colspan="5" class="text-success">No failed liftover records.</td>
        </tr>
        """

    html_output = f"""
<!doctype html>
<html>
<head>
    <meta charset="utf-8">
    <title>taykit liftover report</title>
    <link
        href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.7/dist/css/bootstrap.min.css"
        rel="stylesheet"
    >
</head>
<body class="bg-light">

<div class="container py-5">

    <h1 class="mb-4">taykit liftover report</h1>

    <div class="alert alert-primary">
        <strong>Input:</strong> {html.escape(summary["input_file"])}<br>
        <strong>Provider:</strong> {html.escape(summary["provider"])}<br>
        <strong>Build:</strong> GRCh37 to GRCh38
    </div>

    <div class="row g-3 mb-4">
        <div class="col-md-3">
            <div class="card shadow-sm">
                <div class="card-body">
                    <h6>Total records</h6>
                    <div class="display-6">{summary["total_records"]:,}</div>
                </div>
            </div>
        </div>

        <div class="col-md-3">
            <div class="card shadow-sm">
                <div class="card-body">
                    <h6>Converted</h6>
                    <div class="display-6">{summary["converted_records"]:,}</div>
                </div>
            </div>
        </div>

        <div class="col-md-3">
            <div class="card shadow-sm">
                <div class="card-body">
                    <h6>Failed</h6>
                    <div class="display-6">{summary["failed_records"]:,}</div>
                </div>
            </div>
        </div>
    </div>

    <h2>Output files</h2>

    <div class="alert alert-secondary">
        <strong>Converted:</strong> {html.escape(summary["converted_output"])}<br>
        <strong>Failures:</strong> {html.escape(summary["failed_output"])}
    </div>

    <h2 class="mt-5">Failed records</h2>

    <table class="table table-bordered table-striped bg-white">
        <thead class="table-dark">
            <tr>
                <th>Reason</th>
                <th>rsID</th>
                <th>Chromosome</th>
                <th>Position</th>
                <th>Genotype</th>
            </tr>
        </thead>
        <tbody>
            {failed_table_rows}
        </tbody>
    </table>

</div>

</body>
</html>
"""

    report_file.write_text(html_output, encoding="utf-8")


def write_txt_report(report_file: Path, summary, failed_rows):
    with report_file.open("w", encoding="utf-8") as file:
        file.write("taykit liftover report\n")
        file.write("=" * 80 + "\n\n")

        for key, value in summary.items():
            file.write(f"{key}: {value}\n")

        file.write("\nFailed records\n")
        file.write("-" * 80 + "\n")

        if not failed_rows:
            file.write("No failed records.\n")
        else:
            for row in failed_rows:
                file.write(
                    f"{row.get('reason', '')} | "
                    f"{row.get('rsid', '')} | "
                    f"{row.get('chromosome', '')} | "
                    f"{row.get('position', '')} | "
                    f"{row.get('genotype', '')}\n"
                )


def write_xml_report(report_file: Path, summary, failed_rows):
    root = ET.Element("taykit_liftover_report")

    summary_element = ET.SubElement(root, "summary")

    for key, value in summary.items():
        child = ET.SubElement(summary_element, key)
        child.text = str(value)

    failures = ET.SubElement(root, "failed_records")

    for row in failed_rows:
        failure = ET.SubElement(failures, "failed_record")
        for key, value in row.items():
            child = ET.SubElement(failure, key)
            child.text = str(value)

    tree = ET.ElementTree(root)
    tree.write(report_file, encoding="utf-8", xml_declaration=True)


def write_report(output_dir: Path, report_type, summary, failed_rows):
    if not report_type:
        return None

    report_file = output_dir / f"report.{report_type}"

    if report_type == "html":
        write_html_report(report_file, summary, failed_rows)
    elif report_type == "json":
        with report_file.open("w", encoding="utf-8") as file:
            json.dump({"summary": summary, "failed_records": failed_rows}, file, indent=2)
    elif report_type == "txt":
        write_txt_report(report_file, summary, failed_rows)
    elif report_type == "xml":
        write_xml_report(report_file, summary, failed_rows)

    return report_file


def run_liftover(args):
    if not args.genetic_file:
        print("ERROR: Please provide a genetic file.")
        print()
        print("Example:")
        print("  taykit liftover ancestry.txt --report-type html")
        sys.exit(1)

    input_file = Path(args.genetic_file).expanduser().resolve()

    if not input_file.exists():
        print(f"ERROR: File not found: {input_file}")
        sys.exit(1)

    info = detect_provider_and_build(input_file)
    provider = info["provider"]
    build = info["build"]

    if provider == "Unknown":
        print("ERROR: Unsupported or unknown DNA format.")
        sys.exit(1)

    if build == "GRCh38":
        print("ERROR: This file appears to already be GRCh38.")
        print("No liftover was performed.")
        sys.exit(1)

    if build == "Unknown" and not args.assume_grch37:
        print("ERROR: The file format was recognised, but the genome build could not be detected.")
        print("Use --assume-grch37 only if you are sure this file is GRCh37 / Build 37.")
        sys.exit(1)

    records = load_records(input_file, provider)

    print(f"Detected provider : {provider}")
    print(f"Detected build    : {build}")
    print(f"Loaded records    : {len(records):,}")

    liftover = get_liftover()

    converted_rows = []
    failed_rows = []

    for record in records:
        try:
            position = int(record["position"])
        except ValueError:
            failed_rows.append({**record, "reason": "invalid_position"})
            continue

        lifted = lift_position(liftover, record["chromosome"], position)

        if lifted is None:
            failed_rows.append({**record, "reason": "unmapped"})
            continue

        new_chromosome, new_position = lifted

        converted_rows.append({
            "rsid": record["rsid"],
            "chromosome": new_chromosome,
            "position": str(new_position),
            "genotype": record["genotype"],
        })

    if args.output_path:
        output_dir = Path(args.output_path).expanduser().resolve()
    else:
        output_dir = input_file.parent

    output_dir.mkdir(parents=True, exist_ok=True)

    extension = output_extension(args.output_format)

    converted_file = output_dir / f"{input_file.stem}_GRCh38.{extension}"
    failed_file = output_dir / f"{input_file.stem}_failed_liftover.{extension}"

    write_rows(converted_file, converted_rows, args.output_format)
    write_rows(failed_file, failed_rows, args.output_format)

    summary = {
        "input_file": str(input_file),
        "provider": provider,
        "source_build": build,
        "target_build": "GRCh38",
        "total_records": len(records),
        "converted_records": len(converted_rows),
        "failed_records": len(failed_rows),
        "converted_output": str(converted_file),
        "failed_output": str(failed_file),
    }

    report_file = write_report(output_dir, args.report_type, summary, failed_rows)

    print()
    print(f"Converted records : {len(converted_rows):,}")
    print(f"Failed records    : {len(failed_rows):,}")
    print(f"Converted output  : {converted_file}")
    print(f"Failures output   : {failed_file}")

    if report_file:
        print(f"Report            : {report_file}")
    else:
        print("Report            : not created. Use --report-type html/json/txt/xml to create one.")


def main(args=None):
    run_liftover(args)
