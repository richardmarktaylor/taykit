#!/bin/bash

set -e

# External tools required by taykit impute
if ! command -v brew >/dev/null 2>&1; then
  echo "ERROR: Homebrew is required to install external tools."
  exit 1
fi

brew install bcftools htslib openjdk

# IMPUTE5 may not be available as a standard Homebrew formula on every machine.
# If this fails, install IMPUTE5 manually and ensure `impute5` is on PATH.
if ! command -v impute5 >/dev/null 2>&1; then
  echo "WARNING: impute5 not found on PATH."
  echo "Install IMPUTE5 manually or add it to your taykit tools mirror."
fi

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install pyinstaller pandas pyarrow pyliftover requests tqdm pysam

rm -rf build dist taykit.spec

pyinstaller \
  --onefile \
  --name taykit \
  --collect-submodules taykit \
  taykit/cli.py

echo "Built executable: dist/taykit"
