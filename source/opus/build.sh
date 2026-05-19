#!/bin/bash

set -e

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install pyinstaller

pyinstaller --onefile --name opus opus.py

echo "Built executable: dist/opus"
