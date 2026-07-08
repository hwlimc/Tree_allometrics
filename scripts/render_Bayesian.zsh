#!/bin/zsh
set -euo pipefail

script_dir="${0:A:h}"
project_root="$(cd "$script_dir/.." && pwd -P)"
rmd_file="$script_dir/02_developing_BEF_Allometrics_Bayes.Rmd"
output_dir="$script_dir"
output_file="BEF_Bayesian.md"

usage() {
  cat <<'EOF'
Usage:
  ./render_Bayesian.zsh [--open]

This renders:
  02_developing_BEF_Allometrics_Bayes.Rmd

The Markdown output is written to:
  BEF_Bayesian.md

Options:
  --open    Open the rendered Markdown in the default application after rendering.
  --help    Show this help message.
EOF
}

open_after=false

for arg in "$@"; do
  case "$arg" in
    --open)
      open_after=true
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $arg" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript was not found in PATH." >&2
  exit 1
fi

if [[ ! -f "$rmd_file" ]]; then
  echo "R Markdown file not found: $rmd_file" >&2
  exit 1
fi

mkdir -p "$output_dir"
cd "$project_root"

Rscript --vanilla - "$rmd_file" "$output_dir" "$output_file" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
rmd_file <- args[[1]]
output_dir <- args[[2]]
output_file <- args[[3]]

missing_pkgs <- c("rmarkdown", "knitr")[
  !vapply(c("rmarkdown", "knitr"), requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them with: install.packages(c(\"rmarkdown\", \"knitr\"))",
    call. = FALSE
  )
}

rmarkdown::render(
  input = rmd_file,
  output_format = rmarkdown::github_document(
    toc = TRUE,
    html_preview = FALSE,
    df_print = "default"
  ),
  output_dir = output_dir,
  output_file = output_file,
  intermediates_dir = output_dir,
  clean = TRUE,
  envir = new.env(parent = globalenv())
)
RSCRIPT

md_file="$output_dir/$output_file"
echo "Rendered Markdown:"
echo "$md_file"

if [[ "$open_after" == true ]]; then
  open "$md_file"
fi
