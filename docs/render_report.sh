#!/usr/bin/env bash
# Renders the pages of the analysis report (docs/*.html) from the output of the pipeline and copies the
# figures of the manuscript to docs/figures/.
# Run from the repository root, after `Rscript analysis/run_all.R`:
#     bash docs/render_report.sh
# The exploratory report (docs/exploratory.qmd) is not rendered here: it reruns the 2025 analysis.
set -euo pipefail
cd "$(dirname "$0")"
test -f ../analysis/output/results.json || { echo "Run 'Rscript analysis/run_all.R' first." >&2; exit 1; }
# the figures of the main text of the manuscript, as produced by the pipeline (PNG preview and 600 dpi TIFF)
mkdir -p figures
for f in Fig1 Fig2; do
  cp "../analysis/output/figures/$f.png" "../analysis/output/figures/$f.tif" figures/
done
for page in index data models results sensitivity diagnostics reproducibility; do
  quarto render "$page.qmd"
done
