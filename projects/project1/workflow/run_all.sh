#!/usr/bin/env bash
set -euo pipefail
workflow_dir="$(cd "$(dirname "$0")" && pwd)"
step=${1:?选择分析步骤：10_copykat|11_infercnv}
shift
case "$step" in
  10_copykat|11_infercnv) Rscript "$workflow_dir/$step.R" "$@" ;;
  *) printf '%s\n' '未知步骤' >&2; exit 2 ;;
esac
