#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pysam

try:
    pysam.set_verbosity(0)
except AttributeError:
    pass

COMMAND = "impute"

HELP = "Impute GRCh38 raw DNA files using 1000 Genomes 30x"

DESCRIPTION = "Run strict GRCh38 genotype imputation using 1000 Genomes 30x, Beagle, ShapeIT4 and optional IMPUTE5."

EPILOG = """
Examples:
  taykit impute --wizard
  taykit impute sample.txt output.tsv
  taykit impute sample.txt output.vcf --output-format vcf
  taykit impute sample.txt output.vcf.gz --output-format vcf.gz
  taykit impute sample.txt output.tsv --download-references --threads 32

Notes:
  - Input must already be GRCh38.
  - Run taykit liftover first for GRCh37 files.
  - Reference files are stored in ~/.taykit/imputation.
  - Requires bcftools, bgzip, tabix, java, and optionally shapeit4/impute5.
"""


# =========================================================
# CONFIGURATION
# =========================================================

SCRIPT_DIR = Path.home() / ".taykit" / "imputation"
SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR = SCRIPT_DIR / "temporary_files"
THOUSAND_GENOMES_DIR = SCRIPT_DIR / "reference_1000_genome_30x_grch38"
CLEAN_REF_DIR = SCRIPT_DIR / "reference_1000_genome_30x_grch38_clean"
NCBI_DIR = SCRIPT_DIR / "reference_ncbi"

BCFTOOLS_CMD = "bcftools"
BGZIP_CMD = "bgzip"
TABIX_CMD = "tabix"
JAVA_CMD = "java"
SHAPEIT4_CMD = "shapeit4"
IMPUTE5_CMD = "impute5"
IMP5CONVERTER_CMD = "imp5Converter"

DEFAULT_IMPUTE_CHR_X = False

CPU_COUNT = os.cpu_count() or 1
DEFAULT_THREADS = max(1, CPU_COUNT - 2)
DEFAULT_BEAGLE_MEMORY = "96g"

DEFAULT_BEAGLE_WINDOW_CM = "40.0"
DEFAULT_BEAGLE_OVERLAP_CM = "4.0"
DEFAULT_BEAGLE_WINDOW_MARKERS = "2000000"
DEFAULT_BEAGLE_IMP_STATES = "4000"
DEFAULT_BEAGLE_PHASE_STATES = "800"
DEFAULT_BEAGLE_ITERATIONS = "12"
DEFAULT_BEAGLE_EM = "true"
DEFAULT_BEAGLE_NE = ""
DEFAULT_BEAGLE_ERR = ""

DEFAULT_MIN_QUALITY = 0.8
DEFAULT_LOW_FREQ_MAX_MAF = 0.05
DEFAULT_IMPUTATION_MODE = "beagle"
DEFAULT_PHASING_TOOL = "shapeit4"

RECOMB_MAP_DIR = (
    SCRIPT_DIR
    / "reference_recombination_maps"
    / "plink.GRCh38.map"
    / "no_chr_in_chrom_field"
)
RECOMB_MAP_TEMPLATE = "plink.chr{chrom}.GRCh38.map"

BEAGLE_JAR = SCRIPT_DIR / "reference_beagle" / "beagle.27Feb25.75f.jar"

LOG_FILE = SCRIPT_DIR / "imputation_log_grch38_30x.txt"

DBSNP_GRCH38_CHR_TSV = NCBI_DIR / "dbsnp_GRCh38_chrpos_rsid_snps.chr.tsv.gz"
USE_DBSNP_RSID_LOOKUP = True

AUTOSOMES = [str(c) for c in range(1, 23)]
INPUT_CHROMS = AUTOSOMES + ["X", "Y", "MT"]

REFERENCE_BASE_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20201028_3202_phased"
)

NC_TO_CHR_GRCH38 = {
    "NC_000001.11": "1",
    "NC_000002.12": "2",
    "NC_000003.12": "3",
    "NC_000004.12": "4",
    "NC_000005.10": "5",
    "NC_000006.12": "6",
    "NC_000007.14": "7",
    "NC_000008.11": "8",
    "NC_000009.12": "9",
    "NC_000010.11": "10",
    "NC_000011.10": "11",
    "NC_000012.12": "12",
    "NC_000013.11": "13",
    "NC_000014.9": "14",
    "NC_000015.10": "15",
    "NC_000016.10": "16",
    "NC_000017.11": "17",
    "NC_000018.10": "18",
    "NC_000019.10": "19",
    "NC_000020.11": "20",
    "NC_000021.9": "21",
    "NC_000022.11": "22",
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "MT",
}


# =========================================================
# LOGGING
# =========================================================


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            try:
                stream.write(data)
            except Exception:
                pass

    def flush(self):
        for stream in self.streams:
            try:
                stream.flush()
            except Exception:
                pass


def log(message: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {message}", flush=True)


def _run(cmd, *, check=True, shell=False, capture_output=False, text=False):
    printable = " ".join(str(x) for x in cmd) if isinstance(cmd, list) else str(cmd)
    log(f"[cmd] {printable}")
    return subprocess.run(
        cmd, check=check, shell=shell, capture_output=capture_output, text=text
    )


# =========================================================
# STATS / QC DATA CLASSES
# =========================================================


@dataclass
class PrepareStats:
    chrom: str
    input_variants: int = 0
    rsid_matched: int = 0
    position_fallback_matched: int = 0
    rejected_no_reference_match: int = 0
    rejected_not_biallelic_snp: int = 0
    rejected_ambiguous_snp: int = 0
    rejected_genotype_not_mappable: int = 0
    rejected_duplicate_marker: int = 0
    written_to_prepared_vcf: int = 0

    def log_summary(self) -> None:
        log(
            f"chr{self.chrom} prepare QC: input={self.input_variants:,} "
            f"written={self.written_to_prepared_vcf:,} rsid_match={self.rsid_matched:,} "
            f"pos_fallback={self.position_fallback_matched:,} no_ref={self.rejected_no_reference_match:,} "
            f"non_biallelic={self.rejected_not_biallelic_snp:,} ambiguous={self.rejected_ambiguous_snp:,} "
            f"gt_unmapped={self.rejected_genotype_not_mappable:,} dup={self.rejected_duplicate_marker:,}"
        )


@dataclass
class EngineResult:
    engine: str
    vcf_path: Optional[Path]


# =========================================================
# GENERAL HELPERS
# =========================================================


def chrom_sort_key(chrom: str) -> Tuple[int, int]:
    if chrom.isdigit():
        return (0, int(chrom))
    if chrom == "X":
        return (1, 23)
    if chrom == "Y":
        return (2, 24)
    if chrom == "MT":
        return (3, 25)
    return (9, 999)


def homebrew_available() -> bool:
    return shutil.which("brew") is not None


def install_with_homebrew(package_name: str) -> None:
    if not homebrew_available():
        raise RuntimeError(
            f"Required package is missing: {package_name}. Homebrew is not installed."
        )

    log(f"Installing missing dependency with Homebrew: {package_name}")
    _run(["brew", "install", package_name], check=True)


def require_or_install_command(command: str, brew_package: str) -> None:
    if shutil.which(command):
        return

    install_with_homebrew(brew_package)

    if not shutil.which(command):
        raise RuntimeError(
            f"Tried to install {brew_package}, but command is still missing: {command}"
        )


def ensure_tools_available(args) -> None:
    require_or_install_command(BGZIP_CMD, "htslib")
    require_or_install_command(TABIX_CMD, "htslib")
    require_or_install_command(BCFTOOLS_CMD, "bcftools")
    require_or_install_command(JAVA_CMD, "openjdk")

    if args.phasing_tool == "shapeit4":
        require_or_install_command(SHAPEIT4_CMD, "shapeit4")

    if args.imputation_mode in {"impute5", "both"}:
        require_or_install_command(IMPUTE5_CMD, "impute5")


def ensure_directories() -> None:
    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    THOUSAND_GENOMES_DIR.mkdir(parents=True, exist_ok=True)
    CLEAN_REF_DIR.mkdir(parents=True, exist_ok=True)
    NCBI_DIR.mkdir(parents=True, exist_ok=True)


def download_file(url: str, destination: Path, *, force: bool = False) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
    temp_path = destination.with_suffix(destination.suffix + ".download")
    log(f"Downloading {url}")
    urllib.request.urlretrieve(url, temp_path)
    temp_path.replace(destination)
    log(f"Downloaded {destination}")


def ensure_index(vcf_path: Path) -> None:
    tbi_path = Path(str(vcf_path) + ".tbi")
    csi_path = Path(str(vcf_path) + ".csi")
    if tbi_path.exists() or csi_path.exists():
        return
    log(f"No index for {vcf_path.name}, running tabix ...")
    _run([TABIX_CMD, "-p", "vcf", str(vcf_path)], check=True)


def normalise_chr(chrom: str) -> str:
    value = chrom.strip().replace('"', "")
    if not value:
        return value
    value = value.replace("CHR", "chr").replace("Chr", "chr")
    if value.startswith("chr"):
        value = value[3:]
    value_upper = value.upper()
    if value_upper in NC_TO_CHR_GRCH38:
        return NC_TO_CHR_GRCH38[value_upper]
    if value_upper in {"M", "MT", "MTR", "MITO", "MITOCHONDRIAL"}:
        return "MT"
    if value_upper == "23":
        return "X"
    if value_upper == "24":
        return "Y"
    return value_upper


def normalise_gt_string(gt: str) -> str:
    value = (gt or "").strip().upper().replace('"', "")
    value = value.replace("/", "").replace("|", "")
    if value in {"", "--", "0", "00"}:
        return "--"
    return value


def is_biallelic_snp_record(rec) -> bool:
    if len(rec.ref) != 1 or rec.ref not in "ACGT":
        return False
    if not rec.alts or len(rec.alts) != 1:
        return False
    alt = rec.alts[0]
    return len(alt) == 1 and alt in "ACGT"


def is_ambiguous_snp(ref: str, alt: str) -> bool:
    return {ref, alt} in [{"A", "T"}, {"C", "G"}]


def detect_map_contig_style(map_path: Path) -> str:
    with map_path.open("r") as f_in:
        for line in f_in:
            line = line.strip()
            if line:
                return "chr" if line.split()[0].startswith("chr") else "bare"
    raise RuntimeError(f"Could not detect map contig style from empty map: {map_path}")


def get_or_create_recomb_map_for_style(chrom: str, target_style: str) -> Path:
    source_map = RECOMB_MAP_DIR / RECOMB_MAP_TEMPLATE.format(chrom=chrom)
    if not source_map.exists():
        raise FileNotFoundError(f"Recombination map not found: {source_map}")
    source_style = detect_map_contig_style(source_map)
    if source_style == target_style:
        return source_map

    styled_dir = RECOMB_MAP_DIR.parent / f"{target_style}_in_chrom_field"
    styled_dir.mkdir(parents=True, exist_ok=True)
    out_map = styled_dir / RECOMB_MAP_TEMPLATE.format(chrom=chrom)
    if out_map.exists():
        return out_map

    log(
        f"Normalising recombination map style from {source_style} to {target_style} for chr{chrom}"
    )
    with source_map.open("r") as f_in, out_map.open("w") as f_out:
        for line in f_in:
            parts = line.split()
            if len(parts) < 4:
                continue
            chrom_field = normalise_chr(parts[0])
            parts[0] = f"chr{chrom_field}" if target_style == "chr" else chrom_field
            f_out.write("\t".join(parts) + "\n")
    return out_map


def get_reference_vcf_filename(chrom: str) -> str:
    if chrom == "X":
        return "CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz"
    if chrom in AUTOSOMES:
        return f"CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
    return ""


def ensure_1000g_reference_downloaded(chrom: str, *, download: bool) -> None:
    filename = get_reference_vcf_filename(chrom)
    if not filename:
        return
    vcf_path = THOUSAND_GENOMES_DIR / filename
    tbi_path = Path(str(vcf_path) + ".tbi")
    if vcf_path.exists() and tbi_path.exists():
        return
    if not download:
        raise FileNotFoundError(
            f"Reference file missing: {vcf_path}. Run with --download-references to fetch it."
        )
    download_file(f"{REFERENCE_BASE_URL}/{filename}", vcf_path)
    download_file(f"{REFERENCE_BASE_URL}/{filename}.tbi", tbi_path)


def get_raw_ref_vcf_path(chrom: str) -> Path:
    raw_name = get_reference_vcf_filename(chrom)
    if not raw_name:
        raise FileNotFoundError(
            f"No 30x reference filename defined for chromosome {chrom}"
        )
    raw_path = THOUSAND_GENOMES_DIR / raw_name
    if not raw_path.exists():
        raise FileNotFoundError(
            f"Raw reference VCF for chromosome {chrom} not found: {raw_path}"
        )
    return raw_path


def get_clean_ref_vcf_path(chrom: str) -> Path:
    raw_path = get_raw_ref_vcf_path(chrom)
    clean_name = raw_path.name.replace(".vcf.gz", ".strict.biallelic.dedup.vcf.gz")
    clean_path = CLEAN_REF_DIR / clean_name
    if clean_path.exists() and Path(str(clean_path) + ".tbi").exists():
        return clean_path

    log(f"[ref] Preparing strict biallelic reference for chr{chrom} → {clean_path}")
    proc = subprocess.run(
        [
            BCFTOOLS_CMD,
            "norm",
            "-m",
            "-any",
            "-d",
            "all",
            "-v",
            "snps",
            "-Oz",
            "-o",
            str(clean_path),
            str(raw_path),
        ]
    )
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools norm failed for {raw_path}")
    _run([TABIX_CMD, "-f", "-p", "vcf", str(clean_path)], check=True)
    return clean_path


def detect_reference_contig_name(ref_vcf_path: Path, logical_chrom: str) -> str:
    with pysam.VariantFile(str(ref_vcf_path)) as vcf_in:
        contigs = set(vcf_in.header.contigs.keys())
    candidates = [logical_chrom, f"chr{logical_chrom}"]
    if logical_chrom == "MT":
        candidates.extend(["M", "chrM", "chrMT"])
    for candidate in candidates:
        if candidate in contigs:
            return candidate
    raise RuntimeError(
        f"Could not find contig for chromosome {logical_chrom} in {ref_vcf_path.name}"
    )


def detect_contig_style_from_ref_contig(contig: str) -> str:
    return "chr" if contig.startswith("chr") else "bare"


def build_reference_contig_map(
    impute_chr_x: bool, *, download_references: bool
) -> Dict[str, str]:
    contig_map: Dict[str, str] = {}
    chroms = AUTOSOMES + (["X"] if impute_chr_x else [])
    for chrom in chroms:
        ensure_1000g_reference_downloaded(chrom, download=download_references)
        ref_path = get_clean_ref_vcf_path(chrom)
        contig_map[chrom] = detect_reference_contig_name(ref_path, chrom)
    return contig_map


def detect_tabix_tsv_contig_style(tsv_gz_path: Path) -> str:
    with gzip.open(tsv_gz_path, "rt") as f_in:
        for line in f_in:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            return "chr" if line.split("\t", 1)[0].startswith("chr") else "bare"
    raise RuntimeError(f"Could not detect contig style from empty TSV: {tsv_gz_path}")


def convert_contig_style(chrom: str, style: str) -> str:
    bare = normalise_chr(chrom)
    if style == "chr":
        return "chrMT" if bare == "MT" else f"chr{bare}"
    return bare


def get_or_create_dbsnp_tsv_for_style(target_style: str) -> Path:
    source_tsv = DBSNP_GRCH38_CHR_TSV
    if not source_tsv.exists():
        raise FileNotFoundError(
            f"dbSNP TSV not found: {source_tsv}. This script expects a downloaded/local TSV, not an API call."
        )
    if not Path(str(source_tsv) + ".tbi").exists():
        raise FileNotFoundError(f"dbSNP TSV index not found: {source_tsv}.tbi")
    source_style = detect_tabix_tsv_contig_style(source_tsv)
    if source_style == target_style:
        return source_tsv

    out_tsv = (
        NCBI_DIR / f"dbsnp_GRCh38_chrpos_rsid_snps.{target_style}.normalised.tsv.gz"
    )
    out_tbi = Path(str(out_tsv) + ".tbi")
    if out_tsv.exists() and out_tbi.exists():
        return out_tsv

    log(
        f"Normalising dbSNP contig style from {source_style} to {target_style} → {out_tsv.name}"
    )
    temp_plain = TEMP_DIR / f"{out_tsv.stem}.tmp.tsv"
    with gzip.open(source_tsv, "rt") as f_in, temp_plain.open("w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")
        for line in f_in:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                f_out.write(line + "\n")
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            parts[0] = convert_contig_style(parts[0], target_style)
            writer.writerow(parts)

    with out_tsv.open("wb") as f_out:
        proc_bgzip = subprocess.run([BGZIP_CMD, "-c", str(temp_plain)], stdout=f_out)
    if proc_bgzip.returncode != 0:
        raise RuntimeError(f"bgzip failed while creating {out_tsv}")
    _run([TABIX_CMD, "-f", "-s", "1", "-b", "2", "-e", "2", str(out_tsv)], check=True)
    temp_plain.unlink(missing_ok=True)
    return out_tsv


def validate_reference_files(impute_chr_x: bool, *, download_references: bool) -> None:
    chroms = AUTOSOMES + (["X"] if impute_chr_x else [])
    missing = []
    for chrom in chroms:
        try:
            ensure_1000g_reference_downloaded(chrom, download=download_references)
            get_raw_ref_vcf_path(chrom)
        except FileNotFoundError as e:
            missing.append(str(e))
        source_map = RECOMB_MAP_DIR / RECOMB_MAP_TEMPLATE.format(chrom=chrom)
        if not source_map.exists():
            missing.append(f"Recombination map not found: {source_map}")
    if missing:
        raise FileNotFoundError(
            "Missing required reference files:\n" + "\n".join(missing)
        )


# =========================================================
# INPUT LOADERS
# =========================================================


def build_variants_dict() -> Dict[str, List[Tuple[str, int, str]]]:
    return {chrom: [] for chrom in INPUT_CHROMS}


def resolve_header_indexes(fieldnames: Iterable[str]) -> Tuple[int, int, int, int]:
    raw_names = list(fieldnames)
    names = [f.strip().upper() for f in raw_names]

    def find_one(options: List[str]) -> int:
        for opt in options:
            if opt in names:
                return names.index(opt)
        raise RuntimeError(
            f"Could not find any of header names {options}. Found: {raw_names}"
        )

    idx_rsid = find_one(["RSID", "RSID#", "RSID ID", "SNP", "MARKER", "ID"])
    idx_chr = find_one(["CHROMOSOME", "CHROM", "CHR"])
    idx_pos = find_one(["POSITION", "POS"])
    idx_gt = find_one(["GENOTYPE", "RESULT", "ALLELES", "GENOTYPE CALL"])
    return idx_rsid, idx_chr, idx_pos, idx_gt


def load_delimited_variants(
    input_path: Path, delimiter: str
) -> Dict[str, List[Tuple[str, int, str]]]:
    variants_by_chrom = build_variants_dict()
    with input_path.open(newline="") as f_in:
        reader = csv.reader(f_in, delimiter=delimiter)
        header_indexes = None
        for row in reader:
            if not row:
                continue
            first = (row[0] or "").strip()
            if first.startswith("#"):
                candidate = delimiter.join(row).lstrip("#").strip()
                if not candidate:
                    continue
                cols = [c.strip() for c in candidate.split(delimiter)]
                try:
                    header_indexes = resolve_header_indexes(cols)
                    break
                except RuntimeError:
                    continue
            try:
                header_indexes = resolve_header_indexes(row)
                break
            except RuntimeError:
                continue

        if header_indexes is None:
            raise RuntimeError(f"Could not find usable header in {input_path.name}")

        idx_rsid, idx_chr, idx_pos, idx_gt = header_indexes
        for row in reader:
            if not row or (row[0] or "").strip().startswith("#"):
                continue
            if len(row) <= max(idx_rsid, idx_chr, idx_pos, idx_gt):
                continue
            rsid = (row[idx_rsid] or "").strip()
            chrom = normalise_chr((row[idx_chr] or "").strip())
            pos_str = (row[idx_pos] or "").strip()
            gt_str = normalise_gt_string((row[idx_gt] or "").strip())
            if not chrom or chrom not in variants_by_chrom or gt_str == "--":
                continue
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            variants_by_chrom[chrom].append((rsid, pos, gt_str))

    for chrom in variants_by_chrom:
        variants_by_chrom[chrom].sort(key=lambda x: x[1])
    return variants_by_chrom


def load_input_variants(input_path: Path) -> Dict[str, List[Tuple[str, int, str]]]:
    suffix = input_path.suffix.lower()
    if suffix in {".txt", ".tsv"}:
        log(f"Detected TXT/TSV input: {input_path.name}")
        return load_delimited_variants(input_path, "\t")
    if suffix == ".csv":
        log(f"Detected CSV input: {input_path.name}")
        return load_delimited_variants(input_path, ",")
    raise RuntimeError(
        f"Unrecognised input extension '{suffix}'. Expected .txt, .tsv, or .csv"
    )


# =========================================================
# GENOTYPE / VARIANT MATCHING
# =========================================================


def map_genotype_to_ref_alt(
    my_gt: str, ref: str, alts: Tuple[str, ...] | None
) -> Optional[str]:
    my_gt = normalise_gt_string(my_gt)
    if my_gt == "--":
        return None
    if len(my_gt) == 1 and my_gt in "ACGT":
        my_gt = my_gt * 2
    if len(my_gt) != 2:
        return None
    a1, a2 = my_gt[0], my_gt[1]
    if a1 not in "ACGT" or a2 not in "ACGT":
        return None
    alleles = [ref] + list(alts or [])
    idxs = []
    for allele in (a1, a2):
        if allele not in alleles:
            return None
        idxs.append(str(alleles.index(allele)))
    return "/".join(idxs)


def find_reference_variant_by_rsid_or_position(
    ref_vcf,
    contig: str,
    pos: int,
    rsid: str,
    input_gt: str,
    stats: PrepareStats,
    *,
    allow_position_fallback: bool,
    keep_ambiguous: bool,
):
    candidates = list(ref_vcf.fetch(contig, pos - 1, pos))
    if not candidates:
        stats.rejected_no_reference_match += 1
        return None

    if rsid:
        rsid_candidates = [
            rec
            for rec in candidates
            if rec.id == rsid or (rec.id and rsid in rec.id.split(";"))
        ]
        if rsid_candidates:
            candidates = rsid_candidates
            stats.rsid_matched += 1
        elif not allow_position_fallback:
            stats.rejected_no_reference_match += 1
            return None
        else:
            stats.position_fallback_matched += 1
    else:
        stats.position_fallback_matched += 1

    for rec in candidates:
        if not is_biallelic_snp_record(rec):
            continue
        alt = rec.alts[0]
        if is_ambiguous_snp(rec.ref, alt) and not keep_ambiguous:
            stats.rejected_ambiguous_snp += 1
            continue
        if map_genotype_to_ref_alt(input_gt, rec.ref, rec.alts or ()) is None:
            continue
        return rec

    if any(not is_biallelic_snp_record(rec) for rec in candidates):
        stats.rejected_not_biallelic_snp += 1
    else:
        stats.rejected_genotype_not_mappable += 1
    return None


def gt_to_result_string(rec, sample) -> str:
    gt = sample.get("GT")
    if not gt or any(a is None for a in gt):
        return "--"
    alleles = []
    for a in gt:
        try:
            base = rec.alleles[a]
        except Exception:
            return "--"
        if base is None or len(base) != 1 or base not in {"A", "C", "G", "T"}:
            return "--"
        alleles.append(base)
    if len(alleles) == 1:
        return alleles[0]
    if len(alleles) != 2:
        return "--"
    alleles.sort()
    return "".join(alleles)


# =========================================================
# PREPARED TARGET VCF
# =========================================================


def create_prepared_vcf_for_chromosome(
    chrom: str,
    ref_contig: str,
    variants: List[Tuple[str, int, str]],
    sample_id: str,
    *,
    allow_position_fallback: bool,
    keep_ambiguous: bool,
) -> Tuple[Optional[Path], PrepareStats]:
    stats = PrepareStats(chrom=chrom, input_variants=len(variants))
    if not variants:
        return None, stats

    ref_path = get_clean_ref_vcf_path(chrom)
    ref_vcf = pysam.VariantFile(str(ref_path))
    out_path = TEMP_DIR / f"prepared.chr{chrom}.vcf"
    seen_markers = set()
    out_rows = []

    for rsid, pos, gt_str in variants:
        ref_var = find_reference_variant_by_rsid_or_position(
            ref_vcf,
            ref_contig,
            pos,
            rsid,
            gt_str,
            stats,
            allow_position_fallback=allow_position_fallback,
            keep_ambiguous=keep_ambiguous,
        )
        if ref_var is None:
            continue

        ref_allele = ref_var.ref
        alt_alleles = ref_var.alts or ()
        gt = map_genotype_to_ref_alt(gt_str, ref_allele, alt_alleles)
        if gt is None:
            stats.rejected_genotype_not_mappable += 1
            continue

        alt_str = ",".join(alt_alleles)
        marker_key = (ref_contig, ref_var.pos, ref_allele, alt_str)
        if marker_key in seen_markers:
            stats.rejected_duplicate_marker += 1
            continue
        seen_markers.add(marker_key)

        var_id = rsid if rsid else (ref_var.id or ".")
        out_rows.append(
            (
                ref_var.pos,
                f"{ref_contig}\t{ref_var.pos}\t{var_id}\t{ref_allele}\t{alt_str}\t.\tPASS\t.\tGT\t{gt}\n",
            )
        )
        stats.written_to_prepared_vcf += 1

    ref_vcf.close()
    stats.log_summary()

    if not out_rows:
        out_path.unlink(missing_ok=True)
        return None, stats

    with out_path.open("w") as f_out:
        f_out.write("##fileformat=VCFv4.2\n")
        f_out.write("##source=RawDNA_to_1000G30x_Strict_Prepared_VCF\n")
        f_out.write("##reference=GRCh38\n")
        f_out.write(f"##contig=<ID={ref_contig}>\n")
        f_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f_out.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n"
        )
        for _, line in sorted(out_rows, key=lambda x: x[0]):
            f_out.write(line)
    return out_path, stats


def bgzip_and_tabix_single_vcf(vcf_path: Path) -> Path:
    gz_path = Path(str(vcf_path) + ".gz")
    sort_tmp_dir = TEMP_DIR / "tmp_bcftools_sort"
    sort_tmp_dir.mkdir(parents=True, exist_ok=True)
    gz_path.unlink(missing_ok=True)
    _run(
        [
            BCFTOOLS_CMD,
            "sort",
            "-T",
            str(sort_tmp_dir),
            "-Oz",
            "-o",
            str(gz_path),
            str(vcf_path),
        ],
        check=True,
    )
    _run([TABIX_CMD, "-f", "-p", "vcf", str(gz_path)], check=True)
    vcf_path.unlink(missing_ok=True)
    return gz_path


# =========================================================
# PHASING
# =========================================================


def phase_with_shapeit4(
    chrom: str,
    prepared_vcf_gz: Path,
    ref_vcf_path: Path,
    recomb_map_path: Path,
    threads: int,
) -> Optional[Path]:
    out_bcf = TEMP_DIR / f"phased.shapeit4.chr{chrom}.bcf"
    out_vcf = TEMP_DIR / f"phased.shapeit4.chr{chrom}.vcf.gz"
    out_bcf.unlink(missing_ok=True)
    out_vcf.unlink(missing_ok=True)

    cmd = [
        SHAPEIT4_CMD,
        "--input",
        str(prepared_vcf_gz),
        "--map",
        str(recomb_map_path),
        "--reference",
        str(ref_vcf_path),
        "--region",
        chrom if not str(chrom).startswith("chr") else chrom,
        "--thread",
        str(threads),
        "--output",
        str(out_bcf),
    ]
    log(f"Running ShapeIT4 reference-based phasing for chr{chrom} ...")
    proc = subprocess.run(cmd)
    if proc.returncode != 0 or not out_bcf.exists():
        log(
            f"ShapeIT4 failed for chr{chrom}; falling back to Beagle internal phasing if requested."
        )
        return None

    _run([BCFTOOLS_CMD, "view", "-Oz", "-o", str(out_vcf), str(out_bcf)], check=True)
    _run([TABIX_CMD, "-f", "-p", "vcf", str(out_vcf)], check=True)
    return out_vcf


def get_target_vcf_for_imputation(
    chrom: str,
    prepared_vcf_gz: Path,
    ref_vcf_path: Path,
    recomb_map_path: Path,
    args,
) -> Path:
    if args.phasing_tool != "shapeit4":
        return prepared_vcf_gz
    phased = phase_with_shapeit4(
        chrom, prepared_vcf_gz, ref_vcf_path, recomb_map_path, args.threads
    )
    if phased is None:
        if args.require_shapeit4:
            raise RuntimeError(
                f"ShapeIT4 phasing failed for chr{chrom} and --require-shapeit4 was set"
            )
        return prepared_vcf_gz
    return phased


# =========================================================
# IMPUTATION ENGINES
# =========================================================


def run_beagle_for_chromosome(
    chrom: str,
    target_vcf_gz: Path,
    args,
) -> Optional[Path]:
    if not BEAGLE_JAR.exists():
        raise FileNotFoundError(f"Beagle JAR not found: {BEAGLE_JAR}")
    ref_vcf_path = get_clean_ref_vcf_path(chrom)
    ref_contig = detect_reference_contig_name(ref_vcf_path, chrom)
    ref_style = detect_contig_style_from_ref_contig(ref_contig)
    recomb_map_path = get_or_create_recomb_map_for_style(chrom, ref_style)
    out_prefix = TEMP_DIR / f"imputed.beagle.chr{chrom}"

    cmd = [
        JAVA_CMD,
        f"-Xmx{args.beagle_memory}",
        "-jar",
        str(BEAGLE_JAR),
        f"gt={target_vcf_gz}",
        f"ref={ref_vcf_path}",
        f"map={recomb_map_path}",
        f"out={out_prefix}",
        f"nthreads={args.threads}",
        f"window={args.beagle_window_cm}",
        f"overlap={args.beagle_overlap_cm}",
        f"window-markers={args.beagle_window_markers}",
        f"imp-states={args.beagle_imp_states}",
        f"phase-states={args.beagle_phase_states}",
        f"iterations={args.beagle_iterations}",
        f"em={args.beagle_em}",
    ]
    if args.beagle_ne:
        cmd.append(f"ne={args.beagle_ne}")
    if args.beagle_err:
        cmd.append(f"err={args.beagle_err}")

    log(f"Running Beagle imputation for chr{chrom} ...")
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        log(f"Beagle failed for chr{chrom} with exit code {proc.returncode}")
        return None
    out_vcf = TEMP_DIR / f"imputed.beagle.chr{chrom}.vcf.gz"
    if not out_vcf.exists():
        log(f"Beagle completed for chr{chrom} but output was not found")
        return None
    ensure_index(out_vcf)
    return out_vcf


def run_impute5_for_chromosome(
    chrom: str,
    target_vcf_gz: Path,
    args,
) -> Optional[Path]:
    ref_vcf_path = get_clean_ref_vcf_path(chrom)
    ref_contig = detect_reference_contig_name(ref_vcf_path, chrom)
    ref_style = detect_contig_style_from_ref_contig(ref_contig)
    recomb_map_path = get_or_create_recomb_map_for_style(chrom, ref_style)

    with pysam.VariantFile(str(ref_vcf_path)) as ref_vcf:
        length = ref_vcf.header.contigs[ref_contig].length
    if not length:
        log(f"Could not determine contig length for chr{chrom}; skipping IMPUTE5")
        return None

    region = f"{ref_contig}:1-{length}"
    out_bcf = TEMP_DIR / f"imputed.impute5.chr{chrom}.bcf"
    out_vcf = TEMP_DIR / f"imputed.impute5.chr{chrom}.vcf.gz"
    out_bcf.unlink(missing_ok=True)
    out_vcf.unlink(missing_ok=True)

    cmd = [
        IMPUTE5_CMD,
        "--h",
        str(ref_vcf_path),
        "--g",
        str(target_vcf_gz),
        "--m",
        str(recomb_map_path),
        "--r",
        region,
        "--o",
        str(out_bcf),
        "--threads",
        str(args.threads),
    ]
    log(f"Running IMPUTE5 imputation for chr{chrom} ...")
    proc = subprocess.run(cmd)
    if proc.returncode != 0 or not out_bcf.exists():
        log(f"IMPUTE5 failed for chr{chrom} with exit code {proc.returncode}")
        return None

    _run([BCFTOOLS_CMD, "view", "-Oz", "-o", str(out_vcf), str(out_bcf)], check=True)
    _run([TABIX_CMD, "-f", "-p", "vcf", str(out_vcf)], check=True)
    return out_vcf


def run_imputation_engines(chrom: str, target_vcf_gz: Path, args) -> List[EngineResult]:
    results: List[EngineResult] = []
    if args.imputation_mode in {"beagle", "both"}:
        results.append(
            EngineResult(
                "beagle", run_beagle_for_chromosome(chrom, target_vcf_gz, args)
            )
        )
    if args.imputation_mode in {"impute5", "both"}:
        results.append(
            EngineResult(
                "impute5", run_impute5_for_chromosome(chrom, target_vcf_gz, args)
            )
        )
    return results


# =========================================================
# RSID ANNOTATION
# =========================================================


def annotate_imputed_vcf_with_rsids(
    chrom: str, in_vcf: Path, dbsnp_tsv: Path, engine: str
) -> Path:
    if not USE_DBSNP_RSID_LOOKUP:
        return in_vcf
    ensure_index(in_vcf)
    out_vcf = TEMP_DIR / f"imputed.{engine}.chr{chrom}.rsid.vcf.gz"
    log(f"Annotating RSIDs for chr{chrom} {engine} → {out_vcf.name} ...")
    _run(
        [
            BCFTOOLS_CMD,
            "annotate",
            "--force",
            "-a",
            str(dbsnp_tsv),
            "-c",
            "CHROM,POS,ID",
            "-Oz",
            "-o",
            str(out_vcf),
            str(in_vcf),
        ],
        check=True,
    )
    _run([TABIX_CMD, "-f", "-p", "vcf", str(out_vcf)], check=True)
    return out_vcf


# =========================================================
# OUTPUT / QC
# =========================================================

FinalRow = Tuple[str, str, int, str, str, Optional[float], str, str, str, str, str]
# rsid, chrom, pos, result, source, quality, ref, alt, gt, engine, quality_metric


def extract_info_float(value) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (tuple, list)):
        vals = []
        for item in value:
            try:
                vals.append(float(item))
            except Exception:
                pass
        return max(vals) if vals else None
    try:
        return float(value)
    except Exception:
        return None


def extract_imputation_quality(rec) -> Tuple[Optional[float], str]:
    for key in ("DR2", "R2", "INFO"):
        value = extract_info_float(rec.info.get(key))
        if value is not None:
            return value, key
    return None, ""


def extract_maf_from_record(rec) -> Optional[float]:
    af = extract_info_float(rec.info.get("AF"))
    if af is None:
        af = extract_info_float(rec.info.get("MAF"))
        return af
    return min(af, 1.0 - af)


def collect_preserved_original_rows_for_chrom(
    chrom: str,
    variants_by_chrom: Dict[str, List[Tuple[str, int, str]]],
    final_rows: List[FinalRow],
    written_keys: Set[Tuple[str, int, str]],
) -> int:
    rows = variants_by_chrom.get(chrom, [])
    if not rows:
        return 0
    log(f"Preserving original input rows for chr{chrom} ...")
    count = 0
    for rsid, pos, gt_str in rows:
        gt_str = normalise_gt_string(gt_str)
        if gt_str == "--":
            continue
        key = (chrom, pos, rsid)
        if key in written_keys:
            continue
        written_keys.add(key)
        final_rows.append(
            (
                rsid,
                chrom,
                pos,
                gt_str,
                "original",
                None,
                "N",
                ".",
                "./.",
                "original",
                "",
            )
        )
        count += 1
    return count


def collect_imputed_rows_and_qc_for_chrom(
    chrom: str,
    rsid_vcf: Path,
    final_rows: List[FinalRow],
    written_keys: Set[Tuple[str, int, str]],
    min_quality: float,
    engine: str,
    *,
    low_freq_only: bool,
    low_freq_max_maf: float,
) -> Tuple[int, int, float, int, int]:
    written_count = 0
    quality_n = 0
    quality_sum = 0.0
    quality_ge_08 = 0
    quality_ge_03 = 0

    log(f"Processing {rsid_vcf.name} ...")
    try:
        vcf_in = pysam.VariantFile(str(rsid_vcf), "r")
    except Exception as e:
        log(f"WARNING: Could not open {rsid_vcf.name}: {e}")
        return written_count, quality_n, quality_sum, quality_ge_08, quality_ge_03

    if not vcf_in.header.samples:
        vcf_in.close()
        return written_count, quality_n, quality_sum, quality_ge_08, quality_ge_03

    sample_name = vcf_in.header.samples[0]
    for rec in vcf_in:
        quality, quality_metric = extract_imputation_quality(rec)
        if quality is not None:
            quality_n += 1
            quality_sum += quality
            if quality >= 0.8:
                quality_ge_08 += 1
            if quality >= 0.3:
                quality_ge_03 += 1

        if not is_biallelic_snp_record(rec):
            continue
        if is_ambiguous_snp(rec.ref, rec.alts[0]):
            continue
        if quality is None or quality < min_quality:
            continue

        if low_freq_only:
            maf = extract_maf_from_record(rec)
            if maf is None or maf > low_freq_max_maf:
                continue

        sample = rec.samples[sample_name]
        gt_tuple = sample.get("GT")
        if not gt_tuple or any(a is None for a in gt_tuple):
            continue

        sep = "|" if getattr(sample, "phased", False) else "/"
        gt_field = sep.join(str(a) for a in gt_tuple)
        result = gt_to_result_string(rec, sample)
        if result == "--":
            continue

        rsid = rec.id.split(";")[0] if rec.id and rec.id.startswith("rs") else ""
        logical_chrom = normalise_chr(str(rec.contig))
        key = (logical_chrom, rec.pos, rsid)
        if key in written_keys:
            continue

        written_keys.add(key)
        final_rows.append(
            (
                rsid,
                logical_chrom,
                rec.pos,
                result,
                "imputed",
                quality,
                rec.ref,
                ",".join(rec.alts or ()),
                gt_field,
                engine,
                quality_metric,
            )
        )
        written_count += 1

    vcf_in.close()
    return written_count, quality_n, quality_sum, quality_ge_08, quality_ge_03


def write_final_rows_as_tsv(output_path: Path, final_rows: List[FinalRow]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "RSID",
                "CHROMOSOME",
                "POSITION",
                "RESULT",
                "SOURCE",
                "QUALITY",
                "QUALITY_METRIC",
                "ENGINE",
            ]
        )
        for (
            rsid,
            chrom,
            pos,
            result,
            source,
            quality,
            _ref,
            _alt,
            _gt,
            engine,
            quality_metric,
        ) in sorted(final_rows, key=lambda x: (chrom_sort_key(x[1]), x[2], x[0] or "")):
            writer.writerow(
                [
                    rsid,
                    chrom,
                    pos,
                    result,
                    source,
                    "" if quality is None else f"{quality:.4f}",
                    quality_metric,
                    engine,
                ]
            )


def open_vcf_output(output_path: Path, compress: bool):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(output_path, "wt") if compress else output_path.open("w")


def write_vcf_header(handle, sample_id: str) -> None:
    handle.write("##fileformat=VCFv4.2\n")
    handle.write("##source=imputation_pipeline_grch38_strict_shapeit4_impute5\n")
    handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    handle.write(
        '##INFO=<ID=SOURCE,Number=1,Type=String,Description="original or imputed">\n'
    )
    handle.write(
        '##INFO=<ID=QUALITY,Number=1,Type=Float,Description="Imputation quality score from source engine">\n'
    )
    handle.write(
        '##INFO=<ID=QUALITY_METRIC,Number=1,Type=String,Description="DR2, R2 or INFO">\n'
    )
    handle.write(
        '##INFO=<ID=ENGINE,Number=1,Type=String,Description="original, beagle or impute5">\n'
    )
    for chrom in INPUT_CHROMS:
        handle.write(f"##contig=<ID={chrom}>\n")
    handle.write(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n"
    )


def write_final_rows_as_vcf(
    output_path: Path, final_rows: List[FinalRow], sample_id: str, compress: bool
) -> None:
    with open_vcf_output(output_path, compress) as handle:
        write_vcf_header(handle, sample_id)
        for (
            rsid,
            chrom,
            pos,
            _result,
            source,
            quality,
            ref,
            alt_field,
            gt_field,
            engine,
            quality_metric,
        ) in sorted(final_rows, key=lambda x: (chrom_sort_key(x[1]), x[2], x[0] or "")):
            info_parts = [f"SOURCE={source}", f"ENGINE={engine}"]
            if quality is not None:
                info_parts.append(f"QUALITY={quality:.4f}")
            if quality_metric:
                info_parts.append(f"QUALITY_METRIC={quality_metric}")
            handle.write(
                f"{chrom}\t{pos}\t{rsid if rsid else '.'}\t{ref if ref else 'N'}\t{alt_field if alt_field else '.'}\t.\tPASS\t{';'.join(info_parts)}\tGT\t{gt_field if gt_field else './.'}\n"
            )


def count_input_rows(variants_by_chrom: Dict[str, List[Tuple[str, int, str]]]) -> int:
    return sum(len(variants_by_chrom.get(chrom, [])) for chrom in INPUT_CHROMS)


def write_run_summary(
    summary_path: Path,
    input_file_path: Path,
    output_file: Path,
    output_format: str,
    input_row_count: int,
    total_rows_written: int,
    total_quality_n: int,
    total_quality_sum: float,
    total_quality_ge_08: int,
    total_quality_ge_03: int,
    min_quality: float,
    impute_chr_x: bool,
    args,
    chroms_processed: List[str],
    prepare_stats: List[PrepareStats],
) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    mean_quality = (total_quality_sum / total_quality_n) if total_quality_n else 0.0
    with summary_path.open("w") as f_out:
        f_out.write("Imputation Run Summary\n======================\n\n")
        f_out.write(f"Input file: {input_file_path}\n")
        f_out.write(f"Output file: {output_file}\n")
        f_out.write(f"Output format: {output_format}\n")
        f_out.write(f"Input rows: {input_row_count:,}\n")
        f_out.write(f"Output rows: {total_rows_written:,}\n")
        f_out.write(
            f"Rows added/retained difference: {total_rows_written - input_row_count:,}\n\n"
        )
        f_out.write("Run settings\n------------\n")
        f_out.write("Reference panel: 1000 Genomes 30x GRCh38\n")
        f_out.write(f"Phasing tool: {args.phasing_tool}\n")
        f_out.write(f"Imputation mode: {args.imputation_mode}\n")
        f_out.write(f"Minimum quality threshold: {min_quality}\n")
        f_out.write(f"Impute chrX: {impute_chr_x}\n")
        f_out.write(f"Threads: {args.threads}\n")
        f_out.write(f"Beagle memory: {args.beagle_memory}\n")
        f_out.write(f"Chromosomes processed: {', '.join(chroms_processed)}\n\n")
        f_out.write("Quality summary\n---------------\n")
        f_out.write(f"Variants with quality present: {total_quality_n:,}\n")
        f_out.write(f"Overall mean quality: {mean_quality:.4f}\n")
        f_out.write(f"Variants with quality >= 0.8: {total_quality_ge_08:,}\n")
        f_out.write(f"Variants with quality >= 0.3: {total_quality_ge_03:,}\n\n")
        f_out.write("Prepared VCF strict QC summary\n------------------------------\n")
        for stat in prepare_stats:
            f_out.write(
                f"chr{stat.chrom}: input={stat.input_variants:,}, written={stat.written_to_prepared_vcf:,}, "
                f"rsid_match={stat.rsid_matched:,}, pos_fallback={stat.position_fallback_matched:,}, "
                f"no_ref={stat.rejected_no_reference_match:,}, non_biallelic={stat.rejected_not_biallelic_snp:,}, "
                f"ambiguous={stat.rejected_ambiguous_snp:,}, gt_unmapped={stat.rejected_genotype_not_mappable:,}, "
                f"duplicates={stat.rejected_duplicate_marker:,}\n"
            )


def cleanup_chrom_intermediate_files(chrom: str) -> None:
    patterns = [
        f"prepared.chr{chrom}.vcf*",
        f"phased.shapeit4.chr{chrom}.*",
        f"imputed.beagle.chr{chrom}.*",
        f"imputed.impute5.chr{chrom}.*",
    ]
    for pattern in patterns:
        for path in TEMP_DIR.glob(pattern):
            try:
                log(f"Cleaning up {path}")
                path.unlink()
            except OSError as e:
                log(f"WARNING: could not delete {path}: {e}")


# =========================================================
# WIZARD HELPERS
# =========================================================


def ask_default(prompt: str, default: str) -> str:
    answer = input(f"{prompt} [{default}]: ").strip()
    return answer if answer else default


def ask_yes_no(prompt: str, default: str) -> bool:
    while True:
        answer = input(f"{prompt} [{default}]: ").strip().lower()

        if not answer:
            answer = default.lower()

        if answer in {"y", "yes"}:
            return True

        if answer in {"n", "no"}:
            return False

        print("Please answer yes or no.")


def command_available(command: str) -> bool:
    return shutil.which(command) is not None


def quote_command_item(item: str) -> str:
    if not item:
        return "''"
    if any(character.isspace() for character in item):
        return "'" + item.replace("'", "'\\''") + "'"
    return item


def run_wizard() -> None:
    print("GRCh38 imputation command builder")
    print("================================")
    print()

    input_file = ask_default("Input genotype file", "dna_val_cleansed.txt")
    output_file = ask_default("Output file", "dna_val_imputed_30x.txt")

    input_path = Path(input_file).expanduser()

    if not input_path.exists():
        raise SystemExit(f"Input file does not exist: {input_file}")

    output_format = ask_default("Output format: tsv, vcf, or vcf.gz", "tsv")

    if output_format not in {"tsv", "vcf", "vcf.gz"}:
        raise SystemExit(f"Invalid output format: {output_format}")

    min_quality = ask_default("Minimum imputation quality / DR2 threshold", "0.8")
    threads = ask_default("Threads", str(DEFAULT_THREADS))

    download_references = ask_yes_no(
        "Download missing 1000 Genomes reference files?",
        "yes",
    )

    impute_chr_x = ask_yes_no(
        "Attempt chrX imputation? Usually no unless sex-aware handling is ready",
        "no",
    )

    print()
    print("Phasing options:")
    print("  1) beagle   - safest if ShapeIT4 is not installed")
    print("  2) shapeit4 - reference-based pre-phasing, better but requires ShapeIT4")

    phasing_choice = ask_default("Choose phasing option number", "1")

    if phasing_choice == "1":
        phasing_tool = "beagle"
    elif phasing_choice == "2":
        phasing_tool = "shapeit4"
    else:
        raise SystemExit("Invalid phasing choice")

    require_shapeit4 = False

    if phasing_tool == "shapeit4" and not command_available("shapeit4"):
        print("ShapeIT4 selected but not found.")
        if ask_yes_no("Switch to Beagle internal phasing?", "yes"):
            phasing_tool = "beagle"
        else:
            raise SystemExit("ShapeIT4 is required but was not found.")
    elif phasing_tool == "shapeit4":
        require_shapeit4 = ask_yes_no(
            "Require ShapeIT4 to succeed rather than falling back to Beagle?",
            "no",
        )

    print()
    print("Imputation engine:")
    print("  1) beagle  - default, stable")
    print("  2) impute5 - optional, may help lower-frequency variants")
    print("  3) both    - Beagle first, then IMPUTE5 for extra low-frequency variants")

    imputation_choice = ask_default("Choose imputation option number", "1")

    if imputation_choice == "1":
        imputation_mode = "beagle"
    elif imputation_choice == "2":
        imputation_mode = "impute5"
    elif imputation_choice == "3":
        imputation_mode = "both"
    else:
        raise SystemExit("Invalid imputation choice")

    if imputation_mode in {"impute5", "both"} and not command_available("impute5"):
        print("IMPUTE5 selected but not found.")
        if ask_yes_no("Switch to Beagle-only imputation?", "yes"):
            imputation_mode = "beagle"
        else:
            raise SystemExit("IMPUTE5 is required but was not found.")

    low_freq_max_maf = ask_default(
        "Low-frequency MAF max for extra IMPUTE5 variants when using both",
        "0.05",
    )

    allow_position_fallback = ask_yes_no(
        "Allow position-only matching when rsID does not match? Usually no",
        "no",
    )

    keep_ambiguous_snps = ask_yes_no(
        "Keep strand-ambiguous A/T and C/G SNPs? Usually no",
        "no",
    )

    keep_temp = ask_yes_no("Keep temporary files for debugging?", "no")

    print()
    print("Beagle tuning options:")

    beagle_memory = ask_default("Beagle memory", DEFAULT_BEAGLE_MEMORY)
    beagle_window_cm = ask_default("Beagle window cM", DEFAULT_BEAGLE_WINDOW_CM)
    beagle_overlap_cm = ask_default("Beagle overlap cM", DEFAULT_BEAGLE_OVERLAP_CM)
    beagle_window_markers = ask_default(
        "Beagle window markers",
        DEFAULT_BEAGLE_WINDOW_MARKERS,
    )
    beagle_imp_states = ask_default(
        "Beagle imputation states",
        DEFAULT_BEAGLE_IMP_STATES,
    )
    beagle_phase_states = ask_default(
        "Beagle phasing states",
        DEFAULT_BEAGLE_PHASE_STATES,
    )
    beagle_iterations = ask_default(
        "Beagle iterations",
        DEFAULT_BEAGLE_ITERATIONS,
    )
    beagle_em = ask_default("Beagle EM true/false", DEFAULT_BEAGLE_EM)
    beagle_ne = ask_default("Beagle ne, blank to omit", "")
    beagle_err = ask_default("Beagle err, blank to omit", "")

    command = [
        "taykit",
        "impute",
        input_file,
        output_file,
        "--output-format",
        output_format,
        "--min-quality",
        min_quality,
        "--threads",
        threads,
        "--phasing-tool",
        phasing_tool,
        "--imputation-mode",
        imputation_mode,
        "--low-freq-max-maf",
        low_freq_max_maf,
        "--beagle-memory",
        beagle_memory,
        "--beagle-window-cm",
        beagle_window_cm,
        "--beagle-overlap-cm",
        beagle_overlap_cm,
        "--beagle-window-markers",
        beagle_window_markers,
        "--beagle-imp-states",
        beagle_imp_states,
        "--beagle-phase-states",
        beagle_phase_states,
        "--beagle-iterations",
        beagle_iterations,
        "--beagle-em",
        beagle_em,
    ]

    if download_references:
        command.append("--download-references")

    if impute_chr_x:
        command.append("--impute-chr-x")

    if require_shapeit4:
        command.append("--require-shapeit4")

    if allow_position_fallback:
        command.append("--allow-position-fallback")

    if keep_ambiguous_snps:
        command.append("--keep-ambiguous-snps")

    if keep_temp:
        command.append("--keep-temp")

    if beagle_ne:
        command.extend(["--beagle-ne", beagle_ne])

    if beagle_err:
        command.extend(["--beagle-err", beagle_err])

    print()
    print("Command to run:")
    print(" ".join(quote_command_item(item) for item in command))
    print()

    if ask_yes_no("Run this command now?", "yes"):
        wizard_args = argparse.Namespace(
            wizard=False,
            input_file=input_file,
            output_file=output_file,
            output_format=output_format,
            min_quality=float(min_quality),
            impute_chr_x=impute_chr_x,
            download_references=download_references,
            threads=int(threads),
            phasing_tool=phasing_tool,
            require_shapeit4=require_shapeit4,
            imputation_mode=imputation_mode,
            allow_position_fallback=allow_position_fallback,
            keep_ambiguous_snps=keep_ambiguous_snps,
            low_freq_max_maf=float(low_freq_max_maf),
            beagle_memory=beagle_memory,
            beagle_window_cm=beagle_window_cm,
            beagle_overlap_cm=beagle_overlap_cm,
            beagle_window_markers=beagle_window_markers,
            beagle_imp_states=beagle_imp_states,
            beagle_phase_states=beagle_phase_states,
            beagle_iterations=beagle_iterations,
            beagle_em=beagle_em,
            beagle_ne=beagle_ne,
            beagle_err=beagle_err,
            keep_temp=keep_temp,
        )

        main(wizard_args)
    else:
        print("Command not run.")


# =========================================================
# MAIN
# =========================================================


def configure_parser(parser):
    parser.add_argument(
        "input_file",
        nargs="?",
        help="Path to raw DNA input file .txt, .tsv, or .csv",
    )

    parser.add_argument(
        "output_file",
        nargs="?",
        help="Path to final output file",
    )

    parser.add_argument(
        "--wizard",
        action="store_true",
        help="Launch an interactive imputation command builder.",
    )

    parser.add_argument(
        "--output-format", choices=["tsv", "vcf", "vcf.gz"], default="tsv"
    )

    parser.add_argument(
        "--min-quality",
        type=float,
        default=DEFAULT_MIN_QUALITY,
        help="Minimum DR2/R2/INFO quality to export imputed variants",
    )

    parser.add_argument(
        "--min-dr2",
        type=float,
        dest="min_quality",
        help="Backward-compatible alias for --min-quality",
    )

    parser.add_argument(
        "--impute-chr-x", action="store_true", default=DEFAULT_IMPUTE_CHR_X
    )

    parser.add_argument(
        "--download-references",
        action="store_true",
        help="Download missing 1000 Genomes 30x reference VCFs and indexes",
    )

    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)

    parser.add_argument(
        "--phasing-tool", choices=["shapeit4", "beagle"], default=DEFAULT_PHASING_TOOL
    )

    parser.add_argument("--require-shapeit4", action="store_true")

    parser.add_argument(
        "--imputation-mode",
        choices=["beagle", "impute5", "both"],
        default=DEFAULT_IMPUTATION_MODE,
    )

    parser.add_argument("--allow-position-fallback", action="store_true")
    parser.add_argument("--keep-ambiguous-snps", action="store_true")
    parser.add_argument(
        "--low-freq-max-maf", type=float, default=DEFAULT_LOW_FREQ_MAX_MAF
    )
    parser.add_argument("--beagle-memory", default=DEFAULT_BEAGLE_MEMORY)
    parser.add_argument("--beagle-window-cm", default=DEFAULT_BEAGLE_WINDOW_CM)
    parser.add_argument("--beagle-overlap-cm", default=DEFAULT_BEAGLE_OVERLAP_CM)
    parser.add_argument(
        "--beagle-window-markers", default=DEFAULT_BEAGLE_WINDOW_MARKERS
    )
    parser.add_argument("--beagle-imp-states", default=DEFAULT_BEAGLE_IMP_STATES)
    parser.add_argument("--beagle-phase-states", default=DEFAULT_BEAGLE_PHASE_STATES)
    parser.add_argument("--beagle-iterations", default=DEFAULT_BEAGLE_ITERATIONS)
    parser.add_argument(
        "--beagle-em", choices=["true", "false"], default=DEFAULT_BEAGLE_EM
    )
    parser.add_argument("--beagle-ne", default=DEFAULT_BEAGLE_NE)
    parser.add_argument("--beagle-err", default=DEFAULT_BEAGLE_ERR)
    parser.add_argument("--keep-temp", action="store_true")


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser(
            description=DESCRIPTION,
            epilog=EPILOG,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        configure_parser(parser)
        args = parser.parse_args()

    if getattr(args, "wizard", False):
        run_wizard()
        return

    if not args.input_file or not args.output_file:
        raise SystemExit(
            "ERROR: Provide input_file and output_file, or run:\n"
            "  taykit impute --wizard"
        )

    input_file_path = Path(args.input_file).expanduser().resolve()
    output_file = Path(args.output_file).expanduser().resolve()
    summary_file = output_file.with_suffix(output_file.suffix + ".summary.txt")
    sample_id = input_file_path.stem
    impute_chr_x = bool(args.impute_chr_x)

    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    with LOG_FILE.open("a", buffering=1) as log_f:
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = _Tee(original_stdout, log_f)
        sys.stderr = _Tee(original_stderr, log_f)
        try:
            log("=" * 80)
            log(
                f"Imputation run started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
            )
            log(f"Working directory: {SCRIPT_DIR}")
            log(f"Input file: {input_file_path}")
            log(f"Output file: {output_file}")
            log(f"Output format: {args.output_format}")
            log(f"Reference panel: 1000 Genomes 30x GRCh38")
            log(f"Phasing tool: {args.phasing_tool}")
            log(f"Imputation mode: {args.imputation_mode}")
            log(
                f"Strict QC: biallelic SNPs only, no automatic strand complement, ambiguous SNPs removed unless requested"
            )
            log("=" * 80)

            ensure_directories()
            ensure_tools_available(args)

            if not input_file_path.exists():
                raise FileNotFoundError(f"Input file not found: {input_file_path}")
            if args.imputation_mode in {"beagle", "both"} and not BEAGLE_JAR.exists():
                raise FileNotFoundError(f"Beagle JAR not found: {BEAGLE_JAR}")

            validate_reference_files(
                impute_chr_x, download_references=args.download_references
            )

            log("Reading input file once and partitioning by chromosome ...")
            variants_by_chrom = load_input_variants(input_file_path)
            input_row_count = count_input_rows(variants_by_chrom)
            log(f"Input rows loaded: {input_row_count:,}")

            log("Resolving reference contig names ...")
            ref_contig_map = build_reference_contig_map(
                impute_chr_x, download_references=args.download_references
            )
            ref_style = detect_contig_style_from_ref_contig(ref_contig_map["1"])
            log(f"Detected reference contig style: {ref_style}")
            dbsnp_tsv_for_annotation = get_or_create_dbsnp_tsv_for_style(ref_style)
            log(f"Using dbSNP annotation TSV: {dbsnp_tsv_for_annotation}")

            total_rows_written = 0
            written_keys: Set[Tuple[str, int, str]] = set()
            final_rows: List[FinalRow] = []
            prepare_stats: List[PrepareStats] = []

            total_quality_n = 0
            total_quality_sum = 0.0
            total_quality_ge_08 = 0
            total_quality_ge_03 = 0

            chroms_to_process = AUTOSOMES + (["X"] if impute_chr_x else [])

            for chrom in chroms_to_process:
                log("")
                log(f"========== Chromosome {chrom} ==========")
                chrom_variants = variants_by_chrom.get(chrom, [])
                if not chrom_variants:
                    log(f"No input variants for chr{chrom}.")
                    continue

                prepared_vcf, stats = create_prepared_vcf_for_chromosome(
                    chrom=chrom,
                    ref_contig=ref_contig_map[chrom],
                    variants=chrom_variants,
                    sample_id=sample_id,
                    allow_position_fallback=args.allow_position_fallback,
                    keep_ambiguous=args.keep_ambiguous_snps,
                )
                prepare_stats.append(stats)

                if prepared_vcf is None:
                    total_rows_written += collect_preserved_original_rows_for_chrom(
                        chrom, variants_by_chrom, final_rows, written_keys
                    )
                    continue

                prepared_vcf_gz = bgzip_and_tabix_single_vcf(prepared_vcf)
                ref_vcf_path = get_clean_ref_vcf_path(chrom)
                ref_contig = detect_reference_contig_name(ref_vcf_path, chrom)
                ref_style = detect_contig_style_from_ref_contig(ref_contig)
                recomb_map_path = get_or_create_recomb_map_for_style(chrom, ref_style)
                target_vcf_gz = get_target_vcf_for_imputation(
                    chrom, prepared_vcf_gz, ref_vcf_path, recomb_map_path, args
                )

                engine_results = run_imputation_engines(chrom, target_vcf_gz, args)
                if not engine_results or all(
                    result.vcf_path is None for result in engine_results
                ):
                    total_rows_written += collect_preserved_original_rows_for_chrom(
                        chrom, variants_by_chrom, final_rows, written_keys
                    )
                    if not args.keep_temp:
                        cleanup_chrom_intermediate_files(chrom)
                    continue

                for result in engine_results:
                    if result.vcf_path is None:
                        continue
                    rsid_vcf = annotate_imputed_vcf_with_rsids(
                        chrom, result.vcf_path, dbsnp_tsv_for_annotation, result.engine
                    )
                    low_freq_only = (
                        args.imputation_mode == "both" and result.engine == "impute5"
                    )
                    (
                        written_imputed,
                        quality_n,
                        quality_sum,
                        quality_ge_08,
                        quality_ge_03,
                    ) = collect_imputed_rows_and_qc_for_chrom(
                        chrom=chrom,
                        rsid_vcf=rsid_vcf,
                        final_rows=final_rows,
                        written_keys=written_keys,
                        min_quality=args.min_quality,
                        engine=result.engine,
                        low_freq_only=low_freq_only,
                        low_freq_max_maf=args.low_freq_max_maf,
                    )
                    total_rows_written += written_imputed
                    total_quality_n += quality_n
                    total_quality_sum += quality_sum
                    total_quality_ge_08 += quality_ge_08
                    total_quality_ge_03 += quality_ge_03
                    mean_quality = (quality_sum / quality_n) if quality_n else 0.0
                    log(
                        f"chr{chrom} {result.engine}: quality_variants={quality_n:,} mean={mean_quality:.4f} "
                        f">=0.8={quality_ge_08:,} >=0.3={quality_ge_03:,} exported={written_imputed:,}"
                    )

                preserved_count = collect_preserved_original_rows_for_chrom(
                    chrom, variants_by_chrom, final_rows, written_keys
                )
                total_rows_written += preserved_count
                log(f"chr{chrom}: preserved_original_rows={preserved_count:,}")

                if not args.keep_temp:
                    cleanup_chrom_intermediate_files(chrom)

            for chrom in ["Y", "MT"]:
                preserved_count = collect_preserved_original_rows_for_chrom(
                    chrom, variants_by_chrom, final_rows, written_keys
                )
                total_rows_written += preserved_count
                log(f"chr{chrom}: preserved_original_rows={preserved_count:,}")
                if not args.keep_temp:
                    cleanup_chrom_intermediate_files(chrom)

            if args.output_format == "tsv":
                write_final_rows_as_tsv(output_file, final_rows)
            elif args.output_format == "vcf":
                write_final_rows_as_vcf(
                    output_file, final_rows, sample_id, compress=False
                )
            elif args.output_format == "vcf.gz":
                write_final_rows_as_vcf(
                    output_file, final_rows, sample_id, compress=True
                )
            else:
                raise RuntimeError(f"Unsupported output format: {args.output_format}")

            write_run_summary(
                summary_path=summary_file,
                input_file_path=input_file_path,
                output_file=output_file,
                output_format=args.output_format,
                input_row_count=input_row_count,
                total_rows_written=total_rows_written,
                total_quality_n=total_quality_n,
                total_quality_sum=total_quality_sum,
                total_quality_ge_08=total_quality_ge_08,
                total_quality_ge_03=total_quality_ge_03,
                min_quality=args.min_quality,
                impute_chr_x=impute_chr_x,
                args=args,
                chroms_processed=chroms_to_process + ["Y", "MT"],
                prepare_stats=prepare_stats,
            )
            log(f"Summary file written to {summary_file}")

            overall_mean_quality = (
                (total_quality_sum / total_quality_n) if total_quality_n else 0.0
            )
            log("")
            log("========== Overall QC ==========")
            log(f"Total variants with quality: {total_quality_n:,}")
            log(f"Overall mean quality: {overall_mean_quality:.4f}")
            log(f"Total quality>=0.8: {total_quality_ge_08:,}")
            log(f"Total quality>=0.3: {total_quality_ge_03:,}")
            log(f"Done. Wrote {total_rows_written:,} rows to {output_file}")
            log(
                f"Imputation run finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
            )
            log("=" * 80)

        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr


if __name__ == "__main__":
    main()
