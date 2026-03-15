#!/usr/bin/env bash
set -euo pipefail

python3 scripts/run_extension.py
python3 scripts/run_cognitive_comparison.py
python3 scripts/plot_figure6_comparison.py
