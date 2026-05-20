#!/bin/bash

set -e

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install pyinstaller pandas pyarrow pyliftover

rm -rf build dist taykit.spec

pyinstaller \
  --onefile \
  --name taykit \
  --collect-submodules taykit \
  taykit/cli.py

echo "Built executable: dist/taykit"
