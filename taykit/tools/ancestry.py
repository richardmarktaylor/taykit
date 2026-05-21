#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import platform
import shutil
import subprocess
import tarfile
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import requests
from tqdm import tqdm

COMMAND = "ancestry"

HELP = "Run 1000 Genomes 30x ancestry PCA and ADMIXTURE pipeline"

DESCRIPTION = "Prepare a 1000 Genomes 30x GRCh38 ancestry reference and run ancestry inference for a sample file."

EPILOG = """
Examples:
  taykit ancestry --prepare
  taykit ancestry --prepare --force --threads 32
  taykit ancestry --sample-file sample.txt
  taykit ancestry --sample-file sample.vcf.gz --threads 16
  taykit ancestry --sample-file sample.txt --admixture-k 5

Notes:
  - Requires bcftools, tabix/htslib, PLINK 1.9, PLINK2, and ADMIXTURE.
  - PLINK2 and ADMIXTURE can be downloaded by the tool.
  - bcftools and htslib should be installed with Homebrew:
      brew install bcftools htslib
  - Reference data is large and should not live inside the Homebrew install folder.
"""

PROJECT_DIR = Path.home() / ".taykit" / "ancestry"
PROJECT_DIR.mkdir(parents=True, exist_ok=True)

TOOLS_DIR = PROJECT_DIR / "tools"

REFERENCE_DIR = PROJECT_DIR / "reference_1000_genome_30x_grch38"
REFERENCE_CLEAN_DIR = PROJECT_DIR / "reference_1000_genome_30x_grch38_clean"

WORK_DIR = PROJECT_DIR / "ancestry_work"
REF_WORK_DIR = WORK_DIR / "reference"
SAMPLE_WORK_DIR = WORK_DIR / "samples"

PLINK2_PATH = TOOLS_DIR / "plink2"
ADMIXTURE_PATH = TOOLS_DIR / "admixture"

DEFAULT_THREADS = 32
DEFAULT_ADMIXTURE_K = 5

ADMIXTURE_SNP_TARGET = 300_000
MINOR_ALLELE_FREQUENCY = "0.05"
MISSINGNESS_LIMIT = "0.01"
LD_WINDOW_SIZE = "200"
LD_STEP_SIZE = "50"
LD_R2_LIMIT = "0.2"

PLINK2_MAC_ARM64_URL = (
    "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_arm64_20260425.zip"
)
PLINK2_MAC_X86_64_URL = (
    "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_20260425.zip"
)
PLINK2_LINUX_X86_64_URL = (
    "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260425.zip"
)

ADMIXTURE_MAC_URLS = [
    "https://dalexander.github.io/admixture/binaries/admixture_macosx-1.3.0.tar.gz",
    "https://github.com/NovembreLab/admixture/raw/master/releases/admixture_macosx-1.3.0.tar.gz",
]

BASE_1000G_30X_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20201028_3202_phased"
)

PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)

AUTOSOMES = [str(i) for i in range(1, 23)]

REF_PREFIX = REF_WORK_DIR / "1000g_30x_grch38_autosomes"
REF_PRUNED_PREFIX = REF_WORK_DIR / "1000g_30x_grch38_autosomes_pruned"
REF_PCA_PREFIX = REF_WORK_DIR / "1000g_30x_grch38_pca"

POPULATION_REGION = {
    "ACB": "African",
    "ASW": "African",
    "ESN": "African",
    "GWD": "African",
    "LWK": "African",
    "MSL": "African",
    "YRI": "African",
    "CLM": "American",
    "MXL": "American",
    "PEL": "American",
    "PUR": "American",
    "CDX": "East Asian",
    "CHB": "East Asian",
    "CHS": "East Asian",
    "JPT": "East Asian",
    "KHV": "East Asian",
    "CEU": "European",
    "FIN": "European",
    "GBR": "European",
    "IBS": "European",
    "TSI": "European",
    "BEB": "South Asian",
    "GIH": "South Asian",
    "ITU": "South Asian",
    "PJL": "South Asian",
    "STU": "South Asian",
}


def log(message: str) -> None:
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def run(
    cmd: List[str], cwd: Optional[Path] = None, allow_fail: bool = False
) -> subprocess.CompletedProcess:
    log("Running: " + " ".join(str(item) for item in cmd))
    result = subprocess.run(cmd, cwd=cwd, text=True)
    if result.returncode != 0 and not allow_fail:
        raise SystemExit(
            f"ERROR: Command failed with exit code {result.returncode}: {' '.join(cmd)}"
        )
    return result


def homebrew_available() -> bool:
    return shutil.which("brew") is not None


def install_with_homebrew(package_name: str) -> None:
    if not homebrew_available():
        raise SystemExit(
            f"ERROR: Required package is missing: {package_name}\n\n"
            "Homebrew is not installed, so TayKit cannot install it automatically.\n"
            "Install Homebrew first:\n"
            "  https://brew.sh\n"
        )

    log(f"Installing missing dependency with Homebrew: {package_name}")
    run(["brew", "install", package_name])


def require_command(name: str, brew_package: str) -> str:
    path = shutil.which(name)

    if path:
        return path

    install_with_homebrew(brew_package)

    path = shutil.which(name)

    if not path:
        raise SystemExit(
            f"ERROR: Tried to install {brew_package}, but command is still missing: {name}"
        )

    return path


def ensure_external_dependencies() -> None:
    require_command("bcftools", "bcftools")
    require_command("tabix", "htslib")
    require_command("gzip", "gzip")
    get_plink1_path()


def get_plink1_path() -> str:
    local_plink = TOOLS_DIR / "plink1" / "plink"

    if local_plink.exists():
        return str(local_plink)

    plink1_path = shutil.which("plink")

    if plink1_path:
        return plink1_path

    install_with_homebrew("plink")

    plink1_path = shutil.which("plink")

    if plink1_path:
        return plink1_path

    raise SystemExit(
        "ERROR: PLINK 1.9 is required for merging the reference panel with the sample.\n\n"
        "TayKit tried to install it using:\n"
        "  brew install plink\n\n"
        "but the plink command is still not available."
    )


def require_admixture() -> str:
    if ADMIXTURE_PATH.exists():
        return str(ADMIXTURE_PATH)

    raise SystemExit(
        "ERROR: ADMIXTURE is required because this script always runs ADMIXTURE.\n\n"
        f"Expected local path:\n{ADMIXTURE_PATH}\n\n"
        "Run:\n"
        "./run_ancestry.sh --prepare --force --threads 32\n"
    )


def ensure_dirs() -> None:
    for path in [
        TOOLS_DIR,
        REFERENCE_DIR,
        REFERENCE_CLEAN_DIR,
        WORK_DIR,
        REF_WORK_DIR,
        SAMPLE_WORK_DIR,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def delete_reference_files(*paths: Path) -> None:
    for path in paths:
        path.unlink(missing_ok=True)


def delete_prefix_files(prefix: Path) -> None:
    for path in prefix.parent.glob(prefix.name + ".*"):
        path.unlink(missing_ok=True)


def download_file(
    url: str, out_path: Path, force: bool = False, min_bytes: int = 1
) -> None:
    if out_path.exists() and not force:
        if out_path.stat().st_size >= min_bytes:
            log(f"Already exists, skipping: {out_path}")
            return
        log(f"Existing file looks too small, re-downloading: {out_path}")

    tmp_path = out_path.with_suffix(out_path.suffix + ".download")
    tmp_path.unlink(missing_ok=True)

    log(f"Downloading: {url}")

    with requests.get(url, stream=True, timeout=90) as response:
        response.raise_for_status()
        total = int(response.headers.get("content-length", 0))

        with open(tmp_path, "wb") as handle:
            with tqdm(
                total=total, unit="B", unit_scale=True, desc=out_path.name
            ) as progress:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
                        progress.update(len(chunk))

    if tmp_path.stat().st_size < min_bytes:
        tmp_path.unlink(missing_ok=True)
        raise SystemExit(f"ERROR: Downloaded file is suspiciously small: {out_path}")

    tmp_path.rename(out_path)


def verify_bgzip_vcf(path: Path) -> bool:
    if not path.exists():
        return False

    result = subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    return result.returncode == 0


def verify_tabix_index(vcf_path: Path, tbi_path: Path) -> bool:
    if not vcf_path.exists() or not tbi_path.exists():
        return False

    result = subprocess.run(
        ["tabix", "-H", str(vcf_path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    return result.returncode == 0


def ensure_valid_downloaded_vcf(
    chrom: str,
    vcf_url: str,
    tbi_url: str,
    raw_path: Path,
    tbi_path: Path,
    force: bool = False,
) -> None:
    download_file(vcf_url, raw_path, force=force, min_bytes=10_000_000)
    download_file(tbi_url, tbi_path, force=force, min_bytes=1_000)

    if verify_bgzip_vcf(raw_path) and verify_tabix_index(raw_path, tbi_path):
        return

    log(
        f"Chr{chrom} raw VCF or tabix index failed validation. Deleting and downloading again."
    )

    delete_reference_files(raw_path, tbi_path)

    download_file(vcf_url, raw_path, force=True, min_bytes=10_000_000)
    download_file(tbi_url, tbi_path, force=True, min_bytes=1_000)

    if not verify_bgzip_vcf(raw_path):
        raise SystemExit(
            f"ERROR: chr{chrom} VCF is corrupt after re-download: {raw_path}"
        )

    if not verify_tabix_index(raw_path, tbi_path):
        raise SystemExit(
            f"ERROR: chr{chrom} tabix index is invalid after re-download: {tbi_path}"
        )


def detect_plink2_url() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == "darwin" and machine in {"arm64", "aarch64"}:
        return PLINK2_MAC_ARM64_URL

    if system == "darwin":
        return PLINK2_MAC_X86_64_URL

    if system == "linux" and machine in {"x86_64", "amd64"}:
        return PLINK2_LINUX_X86_64_URL

    raise SystemExit(
        f"ERROR: Unsupported platform for automatic PLINK2 download: {system} {machine}"
    )


def ensure_plink2(force: bool = False) -> None:
    if PLINK2_PATH.exists() and not force:
        log(f"PLINK2 already exists: {PLINK2_PATH}")
        return

    zip_path = TOOLS_DIR / "plink2.zip"
    download_file(detect_plink2_url(), zip_path, force=True, min_bytes=500_000)

    with zipfile.ZipFile(zip_path, "r") as archive:
        archive.extractall(TOOLS_DIR)

    if not PLINK2_PATH.exists():
        found = [path for path in TOOLS_DIR.rglob("plink2") if path.is_file()]
        if found:
            shutil.copy(found[0], PLINK2_PATH)

    if not PLINK2_PATH.exists():
        raise SystemExit(
            "ERROR: PLINK2 was downloaded but the plink2 binary was not found."
        )

    PLINK2_PATH.chmod(0o755)
    log(f"PLINK2 ready: {PLINK2_PATH}")


def ensure_admixture(force: bool = False) -> bool:
    if ADMIXTURE_PATH.exists() and not force:
        log(f"ADMIXTURE already exists: {ADMIXTURE_PATH}")
        return True

    for url in ADMIXTURE_MAC_URLS:
        try:
            tar_path = TOOLS_DIR / "admixture.tar.gz"
            download_file(url, tar_path, force=True, min_bytes=100_000)

            with tarfile.open(tar_path, "r:gz") as archive:
                archive.extractall(TOOLS_DIR)

            found = [path for path in TOOLS_DIR.rglob("admixture") if path.is_file()]
            if found:
                shutil.copy(found[0], ADMIXTURE_PATH)
                ADMIXTURE_PATH.chmod(0o755)
                log(f"ADMIXTURE ready: {ADMIXTURE_PATH}")
                return True

        except Exception as error:
            log(f"ADMIXTURE download attempt failed: {error}")

    log("ADMIXTURE not installed.")
    return False


def get_30x_vcf_name(chrom: str) -> str:
    if chrom == "X":
        return "CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz"

    return f"CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"


def get_30x_clean_vcf_name(chrom: str) -> str:
    return get_30x_vcf_name(chrom).replace(".vcf.gz", ".dedup.vcf.gz")


def download_and_clean_1000g_30x(force: bool = False) -> None:
    require_command("bcftools", "bcftools")
    require_command("tabix", "htslib")
    require_command("gzip", "gzip")

    for chrom in AUTOSOMES:
        vcf_name = get_30x_vcf_name(chrom)
        vcf_url = f"{BASE_1000G_30X_URL}/{vcf_name}"
        tbi_url = f"{vcf_url}.tbi"

        raw_path = REFERENCE_DIR / vcf_name
        tbi_path = Path(str(raw_path) + ".tbi")

        clean_name = get_30x_clean_vcf_name(chrom)
        clean_path = REFERENCE_CLEAN_DIR / clean_name
        clean_tbi_path = Path(str(clean_path) + ".tbi")

        log(f"=== Chr {chrom} ===")

        ensure_valid_downloaded_vcf(
            chrom=chrom,
            vcf_url=vcf_url,
            tbi_url=tbi_url,
            raw_path=raw_path,
            tbi_path=tbi_path,
            force=force,
        )

        if clean_path.exists() and clean_tbi_path.exists() and not force:
            if verify_bgzip_vcf(clean_path) and verify_tabix_index(
                clean_path, clean_tbi_path
            ):
                log(f"Cleaned chr{chrom} already exists and is valid, skipping.")
                continue

            log(f"Cleaned chr{chrom} file is corrupt or invalid. Rebuilding it.")
            delete_reference_files(clean_path, clean_tbi_path)

        if force:
            delete_reference_files(clean_path, clean_tbi_path)

        log(f"Creating deduplicated VCF for chr{chrom}")
        run(
            [
                "bcftools",
                "norm",
                "-d",
                "all",
                "-Oz",
                "-o",
                str(clean_path),
                str(raw_path),
            ]
        )

        if not verify_bgzip_vcf(clean_path):
            delete_reference_files(clean_path, clean_tbi_path)
            raise SystemExit(f"ERROR: Cleaned chr{chrom} VCF failed gzip validation.")

        run(["tabix", "-f", "-p", "vcf", str(clean_path)])

        if not verify_tabix_index(clean_path, clean_tbi_path):
            raise SystemExit(f"ERROR: Cleaned chr{chrom} VCF failed tabix validation.")


def get_clean_autosome_vcfs() -> List[Path]:
    vcfs = []

    for chrom in AUTOSOMES:
        path = REFERENCE_CLEAN_DIR / get_30x_clean_vcf_name(chrom)
        if not path.exists():
            raise SystemExit(f"ERROR: Missing cleaned reference VCF: {path}")
        vcfs.append(path)

    return vcfs


def download_population_panel(force: bool = False) -> Path:
    panel_path = REF_WORK_DIR / "integrated_call_samples_v3.20130502.ALL.panel"
    download_file(PANEL_URL, panel_path, force=force, min_bytes=10_000)
    return panel_path


def read_population_panel(panel_path: Path) -> Dict[str, Dict[str, str]]:
    data: Dict[str, Dict[str, str]] = {}

    with open(panel_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            sample = row.get("sample")
            population = row.get("pop")
            super_population = row.get("super_pop")

            if not sample or not population:
                continue

            data[sample] = {
                "population": population,
                "super_population": super_population
                or POPULATION_REGION.get(population, "Unknown"),
                "region": POPULATION_REGION.get(
                    population, super_population or "Unknown"
                ),
            }

    return data


def build_reference_plink(force: bool = False, threads: int = DEFAULT_THREADS) -> None:
    if REF_PREFIX.with_suffix(".bed").exists() and not force:
        log(f"Reference PLINK files already exist: {REF_PREFIX}")
        return

    concatenated_vcf = REF_WORK_DIR / "1000g_30x_grch38_autosomes.dedup.vcf.gz"

    if not concatenated_vcf.exists() or force:
        vcfs = get_clean_autosome_vcfs()

        if force:
            delete_reference_files(
                concatenated_vcf, Path(str(concatenated_vcf) + ".tbi")
            )

        run(
            [
                "bcftools",
                "concat",
                "--threads",
                str(threads),
                "-Oz",
                "-o",
                str(concatenated_vcf),
                *[str(vcf) for vcf in vcfs],
            ]
        )

        if not verify_bgzip_vcf(concatenated_vcf):
            delete_reference_files(
                concatenated_vcf, Path(str(concatenated_vcf) + ".tbi")
            )
            raise SystemExit(
                "ERROR: Concatenated autosomal VCF failed gzip validation."
            )

        run(["tabix", "-f", "-p", "vcf", str(concatenated_vcf)])

    run(
        [
            str(PLINK2_PATH),
            "--vcf",
            str(concatenated_vcf),
            "--chr",
            "1-22",
            "--max-alleles",
            "2",
            "--snps-only",
            "just-acgt",
            "--set-all-var-ids",
            "@:#:$r:$a",
            "--threads",
            str(threads),
            "--make-bed",
            "--out",
            str(REF_PREFIX),
        ]
    )


def prune_reference(force: bool = False, threads: int = DEFAULT_THREADS) -> None:
    if REF_PRUNED_PREFIX.with_suffix(".bed").exists() and not force:
        log(f"Pruned reference already exists: {REF_PRUNED_PREFIX}")
        return

    common_prefix = REF_WORK_DIR / "1000g_30x_grch38_common_snps"
    prune_prefix = REF_WORK_DIR / "1000g_30x_grch38_prune"

    log("Removing old common/pruned/PCA/ADMIXTURE reference outputs.")
    for prefix in [common_prefix, prune_prefix, REF_PRUNED_PREFIX, REF_PCA_PREFIX]:
        delete_prefix_files(prefix)

    for q_file in REF_WORK_DIR.glob("*.Q"):
        q_file.unlink(missing_ok=True)

    for p_file in REF_WORK_DIR.glob("*.P"):
        p_file.unlink(missing_ok=True)

    log("Filtering reference to common, high-quality autosomal SNPs.")
    run(
        [
            str(PLINK2_PATH),
            "--bfile",
            str(REF_PREFIX),
            "--maf",
            MINOR_ALLELE_FREQUENCY,
            "--geno",
            MISSINGNESS_LIMIT,
            "--threads",
            str(threads),
            "--make-bed",
            "--out",
            str(common_prefix),
        ]
    )

    log("LD-pruning common SNPs.")
    run(
        [
            str(PLINK2_PATH),
            "--bfile",
            str(common_prefix),
            "--indep-pairwise",
            LD_WINDOW_SIZE,
            LD_STEP_SIZE,
            LD_R2_LIMIT,
            "--threads",
            str(threads),
            "--out",
            str(prune_prefix),
        ]
    )

    log(
        f"Building final PCA/ADMIXTURE SNP panel capped at approximately {ADMIXTURE_SNP_TARGET:,} SNPs."
    )
    run(
        [
            str(PLINK2_PATH),
            "--bfile",
            str(common_prefix),
            "--extract",
            str(prune_prefix.with_suffix(".prune.in")),
            "--thin-count",
            str(ADMIXTURE_SNP_TARGET),
            "--seed",
            "42",
            "--threads",
            str(threads),
            "--make-bed",
            "--out",
            str(REF_PRUNED_PREFIX),
        ]
    )


def run_reference_pca(force: bool = False, threads: int = DEFAULT_THREADS) -> None:
    if REF_PCA_PREFIX.with_suffix(".eigenvec").exists() and not force:
        log(f"Reference PCA already exists: {REF_PCA_PREFIX}")
        return

    run(
        [
            str(PLINK2_PATH),
            "--bfile",
            str(REF_PRUNED_PREFIX),
            "--pca",
            "20",
            "allele-wts",
            "--threads",
            str(threads),
            "--out",
            str(REF_PCA_PREFIX),
        ]
    )


def run_reference_admixture(
    force: bool = False,
    k_values: List[int] = [DEFAULT_ADMIXTURE_K],
    threads: int = DEFAULT_THREADS,
) -> None:
    require_admixture()

    for k_value in k_values:
        output_q = REF_WORK_DIR / f"1000g_30x_grch38_autosomes_pruned.{k_value}.Q"
        output_p = REF_WORK_DIR / f"1000g_30x_grch38_autosomes_pruned.{k_value}.P"

        if output_q.exists() and output_p.exists() and not force:
            log(f"ADMIXTURE reference K={k_value} already exists, skipping.")
            continue

        delete_reference_files(output_q, output_p)

        run(
            [
                str(ADMIXTURE_PATH),
                "--cv",
                f"-j{threads}",
                str(REF_PRUNED_PREFIX.with_suffix(".bed")),
                str(k_value),
            ],
            cwd=REF_WORK_DIR,
        )

        produced_q = REF_WORK_DIR / f"{REF_PRUNED_PREFIX.name}.{k_value}.Q"
        produced_p = REF_WORK_DIR / f"{REF_PRUNED_PREFIX.name}.{k_value}.P"

        if produced_q.exists():
            produced_q.rename(output_q)

        if produced_p.exists():
            produced_p.rename(output_p)


def prepare(force: bool = False, threads: int = DEFAULT_THREADS) -> None:
    ensure_dirs()
    ensure_external_dependencies()
    ensure_plink2(force=force)

    if not ensure_admixture(force=force):
        raise SystemExit("ERROR: ADMIXTURE could not be installed or found.")

    download_and_clean_1000g_30x(force=force)
    download_population_panel(force=force)

    build_reference_plink(force=force, threads=threads)
    prune_reference(force=force, threads=threads)
    run_reference_pca(force=force, threads=threads)
    run_reference_admixture(force=force, threads=threads)

    log("Prepare complete. 1000 Genomes 30x GRCh38 ancestry reference files are ready.")


def sample_stem(sample_file: Path) -> str:
    name = sample_file.name

    for suffix in [".vcf.gz", ".gvcf.gz", ".vcf", ".gvcf", ".txt"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]

    return sample_file.stem


def normalise_raw_genotype(value: str) -> str:
    return value.strip().upper().replace("/", "").replace("|", "").replace(" ", "")


def normalise_raw_rsid(value: str) -> str:
    value = value.strip().replace('"', "")

    if value.startswith("GSA-"):
        value = value.replace("GSA-", "", 1)

    return value


def convert_text_to_vcf(sample_file: Path, out_vcf: Path) -> Path:
    opener = gzip.open if sample_file.suffix == ".gz" else open

    rows = []
    header = None
    delimiter = None

    with opener(sample_file, "rt", newline="") as in_handle:
        for raw_line in in_handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith("##"):
                continue

            if line.startswith("#"):
                possible_header = line.lstrip("#").strip()

                if "rsid" in possible_header.lower():
                    header = possible_header
                    delimiter = "\t" if "\t" in possible_header else ","

                continue

            if header is None:
                lower_line = line.lower()

                if (
                    lower_line.startswith("rsid")
                    or lower_line.startswith("snp name")
                    or lower_line.startswith("snp,")
                    or lower_line.startswith("rsid,")
                ):
                    header = line
                    delimiter = "\t" if "\t" in line else ","
                    continue

            if header is None:
                continue

            if delimiter is None:
                delimiter = "\t" if "\t" in header else ","

            reader = csv.reader([line], delimiter=delimiter)
            parts = next(reader)

            header_reader = csv.reader([header], delimiter=delimiter)
            header_parts = [item.strip().lower() for item in next(header_reader)]

            row_map = {}

            for index, field_name in enumerate(header_parts):
                if index < len(parts):
                    row_map[field_name] = parts[index].strip().replace('"', "")

            def get_value(*names: str) -> str:
                for name in names:
                    if name.lower() in row_map:
                        return row_map[name.lower()]
                return ""

            rsid = normalise_raw_rsid(get_value("rsid", "snp", "snp name"))

            chromosome = (
                get_value("chromosome", "chrom", "chr")
                .strip()
                .replace("chr", "")
                .replace("CHR", "")
            )

            position = get_value("position", "pos", "chrposition").strip()

            genotype = normalise_raw_genotype(
                get_value("genotype", "result", "plus strand", "alleles")
            )

            if not genotype:
                allele1 = get_value("allele1", "allele 1").strip().upper()
                allele2 = get_value("allele2", "allele 2").strip().upper()
                genotype = normalise_raw_genotype(allele1 + allele2)

            if not rsid or not chromosome or not position or len(genotype) != 2:
                continue

            allele_one = genotype[0]
            allele_two = genotype[1]

            if allele_one not in "ACGT" or allele_two not in "ACGT":
                continue

            try:
                pos_int = int(position)
            except ValueError:
                continue

            chrom_sort = (
                23
                if chromosome == "X"
                else int(chromosome)
                if chromosome.isdigit()
                else 99
            )

            reference_allele = allele_one
            alternate_allele = allele_two if allele_two != reference_allele else "."
            genotype_value = "0/0" if allele_one == allele_two else "0/1"

            rows.append(
                (
                    chrom_sort,
                    pos_int,
                    chromosome,
                    position,
                    rsid,
                    reference_allele,
                    alternate_allele,
                    genotype_value,
                )
            )

    if not rows:
        raise SystemExit(
            "ERROR: No usable SNP rows were found.\n\n"
            "Supported text formats include:\n"
            "  - AncestryDNA: rsid chromosome position allele1 allele2\n"
            "  - SelfDecode: rsid chromosome position genotype\n"
            "  - 23andMe: rsid chromosome position genotype\n"
            "  - MyHeritage: RSID,CHROMOSOME,POSITION,RESULT\n"
            "  - TayKit/GSLS-style: SNP Name, Chr, Position, Plus Strand\n"
        )

    rows.sort(key=lambda item: (item[0], item[1], item[4]))

    sample_name = sample_stem(sample_file)

    with open(out_vcf, "w") as out_handle:
        out_handle.write("##fileformat=VCFv4.2\n")
        out_handle.write("##source=taykit_ancestry_text_converter\n")
        out_handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out_handle.write(sample_name + "\n")

        for (
            _,
            _,
            chromosome,
            position,
            rsid,
            reference_allele,
            alternate_allele,
            genotype_value,
        ) in rows:
            out_handle.write(
                f"{chromosome}\t{position}\t{rsid}\t"
                f"{reference_allele}\t{alternate_allele}\t.\tPASS\t.\tGT\t{genotype_value}\n"
            )

    return out_vcf


def normalise_sample_to_bed(
    sample_file: Path,
    sample_prefix: Path,
    threads: int = DEFAULT_THREADS,
) -> None:
    if not sample_file.exists():
        raise SystemExit(f"ERROR: Sample file does not exist: {sample_file}")

    suffixes = "".join(sample_file.suffixes).lower()
    working_vcf = sample_file

    if suffixes.endswith(".txt"):
        working_vcf = sample_prefix.with_suffix(".converted.vcf")
        convert_text_to_vcf(sample_file, working_vcf)

    if suffixes.endswith(".gvcf") or suffixes.endswith(".gvcf.gz"):
        log(
            "Warning: gVCF detected. The script treats it as VCF input and ignores reference blocks."
        )

    sorted_prefix = sample_prefix.parent / f"{sample_prefix.name}_sorted"

    run(
        [
            str(PLINK2_PATH),
            "--vcf",
            str(working_vcf),
            "--chr",
            "1-22",
            "--max-alleles",
            "2",
            "--snps-only",
            "just-acgt",
            "--set-all-var-ids",
            "@:#:$r:$a",
            "--threads",
            str(threads),
            "--make-pgen",
            "--sort-vars",
            "--out",
            str(sorted_prefix),
        ]
    )

    run(
        [
            str(PLINK2_PATH),
            "--pfile",
            str(sorted_prefix),
            "--threads",
            str(threads),
            "--make-bed",
            "--out",
            str(sample_prefix),
        ]
    )


def extract_sample_overlap(
    sample_prefix: Path,
    sample_overlap_prefix: Path,
    threads: int = DEFAULT_THREADS,
) -> None:
    prune_file = REF_WORK_DIR / "1000g_30x_grch38_prune.prune.in"

    run(
        [
            str(PLINK2_PATH),
            "--bfile",
            str(sample_prefix),
            "--extract",
            str(prune_file),
            "--threads",
            str(threads),
            "--make-bed",
            "--out",
            str(sample_overlap_prefix),
        ]
    )


def extract_reference_overlap_for_sample(
    sample_overlap_prefix: Path,
    reference_overlap_prefix: Path,
    threads: int = DEFAULT_THREADS,
) -> None:
    run(
        [
            str(PLINK2_PATH),
            "--bfile",
            str(REF_PRUNED_PREFIX),
            "--extract",
            str(sample_overlap_prefix.with_suffix(".bim")),
            "--threads",
            str(threads),
            "--make-bed",
            "--out",
            str(reference_overlap_prefix),
        ]
    )


def merge_reference_and_sample(
    reference_prefix: Path,
    sample_overlap_prefix: Path,
    merged_prefix: Path,
    threads: int = DEFAULT_THREADS,
) -> bool:
    plink1_path = get_plink1_path()

    result = run(
        [
            plink1_path,
            "--bfile",
            str(reference_prefix),
            "--bmerge",
            str(sample_overlap_prefix.with_suffix(".bed")),
            str(sample_overlap_prefix.with_suffix(".bim")),
            str(sample_overlap_prefix.with_suffix(".fam")),
            "--threads",
            str(threads),
            "--make-bed",
            "--out",
            str(merged_prefix),
        ],
        allow_fail=True,
    )

    if result.returncode == 0:
        return True

    missnp_candidates = [
        merged_prefix.with_suffix(".missnp"),
        Path(str(merged_prefix) + ".missnp"),
    ]

    missnp = next((path for path in missnp_candidates if path.exists()), None)

    if not missnp:
        return False

    flipped_prefix = (
        sample_overlap_prefix.parent / f"{sample_overlap_prefix.name}_flipped"
    )

    log("Merge failed because of allele mismatches. Trying strand flip using .missnp.")

    run(
        [
            plink1_path,
            "--bfile",
            str(sample_overlap_prefix),
            "--flip",
            str(missnp),
            "--threads",
            str(threads),
            "--make-bed",
            "--out",
            str(flipped_prefix),
        ]
    )

    second_result = run(
        [
            plink1_path,
            "--bfile",
            str(reference_prefix),
            "--bmerge",
            str(flipped_prefix.with_suffix(".bed")),
            str(flipped_prefix.with_suffix(".bim")),
            str(flipped_prefix.with_suffix(".fam")),
            "--threads",
            str(threads),
            "--make-bed",
            "--out",
            str(merged_prefix),
        ],
        allow_fail=True,
    )

    return second_result.returncode == 0


def run_joint_pca(
    merged_prefix: Path,
    pca_prefix: Path,
    threads: int = DEFAULT_THREADS,
) -> None:
    run(
        [
            str(PLINK2_PATH),
            "--bfile",
            str(merged_prefix),
            "--pca",
            "20",
            "--threads",
            str(threads),
            "--out",
            str(pca_prefix),
        ]
    )


def read_fam_iids(fam_path: Path) -> List[str]:
    with open(fam_path, "r") as handle:
        return [line.strip().split()[1] for line in handle if line.strip()]


def read_q_matrix(q_path: Path) -> List[List[float]]:
    rows: List[List[float]] = []

    with open(q_path, "r") as handle:
        for line in handle:
            if line.strip():
                rows.append([float(value) for value in line.strip().split()])

    return rows


def run_sample_admixture(
    merged_prefix: Path,
    sample_id: str,
    panel_path: Path,
    k_value: int = DEFAULT_ADMIXTURE_K,
    threads: int = DEFAULT_THREADS,
    force: bool = False,
) -> Dict:
    admixture_path = require_admixture()

    q_path = Path(str(merged_prefix) + f".{k_value}.Q")
    p_path = Path(str(merged_prefix) + f".{k_value}.P")

    if force:
        delete_reference_files(q_path, p_path)

    if not q_path.exists() or not p_path.exists():
        run(
            [
                admixture_path,
                "--cv",
                f"-j{threads}",
                str(merged_prefix.with_suffix(".bed")),
                str(k_value),
            ],
            cwd=merged_prefix.parent,
        )

    if not q_path.exists():
        raise SystemExit(f"ERROR: ADMIXTURE did not create Q file: {q_path}")

    fam_iids = read_fam_iids(merged_prefix.with_suffix(".fam"))
    q_rows = read_q_matrix(q_path)
    panel = read_population_panel(panel_path)

    if len(fam_iids) != len(q_rows):
        raise SystemExit(
            f"ERROR: ADMIXTURE Q row count does not match FAM sample count: "
            f"{len(q_rows)} Q rows vs {len(fam_iids)} FAM rows"
        )

    if sample_id not in fam_iids:
        raise SystemExit(f"ERROR: Sample ID not found in merged FAM: {sample_id}")

    sample_index = fam_iids.index(sample_id)
    sample_q = q_rows[sample_index]

    cluster_region_totals: List[Dict[str, float]] = [{} for _ in range(k_value)]

    for iid, q_values in zip(fam_iids, q_rows):
        if iid == sample_id:
            continue

        if iid not in panel:
            continue

        region = panel[iid]["region"]

        for cluster_index, value in enumerate(q_values):
            cluster_region_totals[cluster_index][region] = (
                cluster_region_totals[cluster_index].get(region, 0.0) + value
            )

    cluster_labels = []

    for cluster_index, totals in enumerate(cluster_region_totals):
        if not totals:
            cluster_labels.append(f"Cluster_{cluster_index + 1}")
            continue

        best_region = max(totals.items(), key=lambda item: item[1])[0]
        cluster_labels.append(best_region)

    cluster_percentages = {
        f"Cluster_{index + 1}_{cluster_labels[index]}": round(value * 100, 2)
        for index, value in enumerate(sample_q)
    }

    region_percentages: Dict[str, float] = {}

    for index, value in enumerate(sample_q):
        region = cluster_labels[index]
        region_percentages[region] = region_percentages.get(region, 0.0) + value * 100

    region_percentages = {
        region: round(value, 2)
        for region, value in sorted(
            region_percentages.items(), key=lambda item: item[1], reverse=True
        )
    }

    return {
        "method": f"ADMIXTURE K={k_value} on sample-overlap SNPs only",
        "k": k_value,
        "q_file": str(q_path),
        "p_file": str(p_path),
        "cluster_labels": {
            f"Cluster_{index + 1}": cluster_labels[index] for index in range(k_value)
        },
        "cluster_percentages": cluster_percentages,
        "region_percentages": region_percentages,
        "note": (
            "ADMIXTURE was run on the SNPs overlapping between the sample and the reference panel. "
            "This is faster, but low-overlap samples may produce less stable ADMIXTURE estimates."
        ),
    }


def read_eigenvec(eigenvec_path: Path) -> Dict[str, List[float]]:
    eigenvectors: Dict[str, List[float]] = {}

    with open(eigenvec_path, "r") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.strip().split()
            if len(parts) < 4:
                continue

            individual_id = parts[1]
            principal_components = [float(value) for value in parts[2:]]
            eigenvectors[individual_id] = principal_components

    return eigenvectors


def classify_by_nearest_reference(
    eigenvec_path: Path,
    panel_path: Path,
    sample_id: str,
    top_n: int = 20,
) -> Dict:
    eigenvectors = read_eigenvec(eigenvec_path)
    panel = read_population_panel(panel_path)

    if sample_id not in eigenvectors:
        raise SystemExit(f"ERROR: Could not find sample ID in PCA output: {sample_id}")

    sample_pcs = eigenvectors[sample_id][:10]
    distances = []

    for reference_sample, reference_pcs in eigenvectors.items():
        if reference_sample == sample_id:
            continue

        if reference_sample not in panel:
            continue

        reference_pcs = reference_pcs[:10]
        distance = (
            sum(
                (sample_value - reference_value) ** 2
                for sample_value, reference_value in zip(sample_pcs, reference_pcs)
            )
            ** 0.5
        )

        distances.append(
            {
                "sample": reference_sample,
                "distance": distance,
                "population": panel[reference_sample]["population"],
                "super_population": panel[reference_sample]["super_population"],
                "region": panel[reference_sample]["region"],
            }
        )

    distances.sort(key=lambda row: row["distance"])
    nearest = distances[:top_n]

    if not nearest:
        raise SystemExit("ERROR: No nearby 1000 Genomes reference samples were found.")

    population_counts: Dict[str, int] = {}
    region_counts: Dict[str, int] = {}

    for row in nearest:
        population_counts[row["population"]] = (
            population_counts.get(row["population"], 0) + 1
        )
        region_counts[row["region"]] = region_counts.get(row["region"], 0) + 1

    population_percentages = {
        population: round(count / len(nearest) * 100, 2)
        for population, count in sorted(
            population_counts.items(), key=lambda item: item[1], reverse=True
        )
    }

    region_percentages = {
        region: round(count / len(nearest) * 100, 2)
        for region, count in sorted(
            region_counts.items(), key=lambda item: item[1], reverse=True
        )
    }

    return {
        "method": "PLINK2 joint PCA against sample-overlap 1000 Genomes 30x reference",
        "principal_components_used": 10,
        "nearest_reference_samples_used": top_n,
        "region_percentages": region_percentages,
        "population_percentages": population_percentages,
        "nearest_reference_samples": nearest,
    }


def get_sample_id_from_fam(fam_path: Path) -> str:
    with open(fam_path, "r") as handle:
        rows = [line.strip().split() for line in handle if line.strip()]

    if not rows:
        raise SystemExit(
            "ERROR: No samples found after converting sample file to PLINK format."
        )

    if len(rows) != 1:
        log(f"Warning: sample file produced {len(rows)} samples. Using the first one.")

    return rows[0][1]


def get_variant_count_from_bim(bim_path: Path) -> int:
    with open(bim_path, "r") as handle:
        return sum(1 for line in handle if line.strip())


def process_sample(
    sample_file: Path,
    force: bool = False,
    threads: int = DEFAULT_THREADS,
    admixture_k: int = DEFAULT_ADMIXTURE_K,
) -> Path:
    ensure_dirs()
    ensure_external_dependencies()

    if not REF_PRUNED_PREFIX.with_suffix(".bed").exists():
        log("Reference files are missing. Running prepare first.")
        prepare(force=False, threads=threads)

    require_admixture()

    stem = sample_stem(sample_file)
    sample_dir = SAMPLE_WORK_DIR / stem
    sample_dir.mkdir(parents=True, exist_ok=True)

    sample_prefix = sample_dir / stem
    sample_overlap_prefix = sample_dir / f"{stem}_overlap"
    reference_overlap_prefix = sample_dir / f"{stem}_reference_overlap"
    merged_prefix = sample_dir / f"{stem}_merged_with_1000g_30x"
    pca_prefix = sample_dir / f"{stem}_joint_pca"

    output_json = PROJECT_DIR / f"{stem}-ancestry.json"

    if output_json.exists() and not force:
        log(f"Output already exists, skipping: {output_json}")
        return output_json

    normalise_sample_to_bed(sample_file, sample_prefix, threads=threads)
    sample_id = get_sample_id_from_fam(sample_prefix.with_suffix(".fam"))

    extract_sample_overlap(sample_prefix, sample_overlap_prefix, threads=threads)
    overlap_count = get_variant_count_from_bim(
        sample_overlap_prefix.with_suffix(".bim")
    )

    log(f"Sample/reference ancestry SNP overlap: {overlap_count:,} variants.")

    extract_reference_overlap_for_sample(
        sample_overlap_prefix=sample_overlap_prefix,
        reference_overlap_prefix=reference_overlap_prefix,
        threads=threads,
    )

    reference_overlap_count = get_variant_count_from_bim(
        reference_overlap_prefix.with_suffix(".bim")
    )
    log(
        f"Reduced reference SNP count for this sample: {reference_overlap_count:,} variants."
    )

    if reference_overlap_count != overlap_count:
        log(
            "Warning: sample overlap count and reference overlap count differ. "
            "This may indicate duplicate or incompatible variant IDs."
        )

    if not merge_reference_and_sample(
        reference_prefix=reference_overlap_prefix,
        sample_overlap_prefix=sample_overlap_prefix,
        merged_prefix=merged_prefix,
        threads=threads,
    ):
        raise SystemExit(
            "ERROR: Could not merge sample with reduced 1000 Genomes 30x reference.\n"
            "Likely causes: genome build mismatch, too few overlapping SNPs, strand mismatch, or incompatible variant IDs."
        )

    merged_count = get_variant_count_from_bim(merged_prefix.with_suffix(".bim"))
    log(f"Merged PCA/ADMIXTURE dataset SNP count: {merged_count:,} variants.")

    run_joint_pca(merged_prefix, pca_prefix, threads=threads)

    panel_path = REF_WORK_DIR / "integrated_call_samples_v3.20130502.ALL.panel"

    pca_result = classify_by_nearest_reference(
        eigenvec_path=pca_prefix.with_suffix(".eigenvec"),
        panel_path=panel_path,
        sample_id=sample_id,
        top_n=20,
    )

    admixture_result = run_sample_admixture(
        merged_prefix=merged_prefix,
        sample_id=sample_id,
        panel_path=panel_path,
        k_value=admixture_k,
        threads=threads,
        force=force,
    )

    result = {
        "sample_file": str(sample_file),
        "sample_name": stem,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "reference": {
            "name": "1000 Genomes 30x GRCh38 phased reference",
            "source": BASE_1000G_30X_URL,
            "chromosomes_used": "1-22",
            "x_used": False,
            "y_used": False,
            "note": "X was deliberately excluded from ancestry PCA/ADMIXTURE preparation.",
        },
        "sample_qc": {
            "sample_overlap_variants": overlap_count,
            "reference_overlap_variants": reference_overlap_count,
            "merged_variants_used": merged_count,
            "note": (
                "The reference was reduced to the customer's overlapping SNPs before PCA and ADMIXTURE "
                "to avoid running ancestry models across SNPs missing from the customer sample."
            ),
        },
        "plink2": {
            "path": str(PLINK2_PATH),
            "threads": threads,
        },
        "plink1": {
            "path": get_plink1_path(),
            "used_for": "sample-wise reference merge via --bmerge",
        },
        "admixture": {
            "available": True,
            "path": str(ADMIXTURE_PATH),
            "k": admixture_k,
            "snp_target": ADMIXTURE_SNP_TARGET,
            "minor_allele_frequency": MINOR_ALLELE_FREQUENCY,
            "missingness_limit": MISSINGNESS_LIMIT,
            "ld_pruning": {
                "window_size": LD_WINDOW_SIZE,
                "step_size": LD_STEP_SIZE,
                "r2_limit": LD_R2_LIMIT,
            },
            "threads": threads,
        },
        "ancestry": {
            "pca_nearest_reference": pca_result,
            "admixture": admixture_result,
        },
        "limitations": [
            "This is broad ancestry inference, not a clinical diagnosis.",
            "PCA percentages are nearest-reference approximations.",
            "ADMIXTURE percentages are model-based cluster proportions, not literal ethnic identity labels.",
            "ADMIXTURE cluster names are inferred from the dominant 1000 Genomes region in each cluster.",
            "Accuracy depends heavily on the overlap between the customer file and the prepared 1000 Genomes 30x SNP panel.",
            "Low SNP overlap can make ADMIXTURE estimates unstable.",
            "Small SNP-array files may produce unstable estimates.",
            "gVCF files are treated conservatively as VCF input; reference blocks are not used.",
            "1000 Genomes population labels are reference labels, not exhaustive ethnic identity labels.",
        ],
    }

    with open(output_json, "w") as handle:
        json.dump(result, handle, indent=2)

    log(f"Wrote: {output_json}")
    return output_json


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="PLINK2 + PLINK 1.9 + ADMIXTURE + 1000 Genomes 30x GRCh38 ancestry pipeline"
    )

    parser.add_argument(
        "--prepare",
        action="store_true",
        help="Download tools/reference data and build reusable reference files, then exit.",
    )

    parser.add_argument(
        "--sample-file",
        type=Path,
        help="Input sample file: .txt, .vcf, .gvcf, .vcf.gz, or .gvcf.gz",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Force rebuilding/redownloading where possible.",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="Number of CPU threads to use for PLINK2, PLINK 1.9, bcftools concat, and ADMIXTURE.",
    )

    parser.add_argument(
        "--admixture-k",
        type=int,
        default=DEFAULT_ADMIXTURE_K,
        help="Number of ADMIXTURE ancestral clusters to model.",
    )

    return parser.parse_args()


def configure_parser(parser):
    parser.add_argument(
        "--prepare",
        action="store_true",
        help="Download tools/reference data and build reusable reference files, then exit.",
    )

    parser.add_argument(
        "--sample-file",
        type=Path,
        help="Input sample file: .txt, .vcf, .gvcf, .vcf.gz, or .gvcf.gz",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Force rebuilding/redownloading where possible.",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="Number of CPU threads to use.",
    )

    parser.add_argument(
        "--admixture-k",
        type=int,
        default=DEFAULT_ADMIXTURE_K,
        help="Number of ADMIXTURE ancestral clusters to model.",
    )


def main(args=None) -> None:
    if args is None:
        parser = argparse.ArgumentParser(
            description=DESCRIPTION,
            epilog=EPILOG,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        configure_parser(parser)
        args = parser.parse_args()

    if args.prepare:
        prepare(force=args.force, threads=args.threads)
        return

    if not args.sample_file:
        raise SystemExit("ERROR: Provide either --prepare or --sample-file")

    ensure_dirs()
    ensure_plink2(force=False)
    process_sample(
        args.sample_file,
        force=args.force,
        threads=args.threads,
        admixture_k=args.admixture_k,
    )


if __name__ == "__main__":
    main()
