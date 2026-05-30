#!/bin/zsh

PROJECT_DIR="/Users/hyli0001/wrd/b/Dynamic_allometrics/"
SCRIPT="$PROJECT_DIR/original_files/cleaned_data/kr_biomass_dataset.R"
DATE=$(date +%Y%m%d_%H%M%S)
LOG="$PROJECT_DIR/scripts/logfiles/log_data_prep_${DATE}.txt"

cd "$PROJECT_DIR" || exit 1
echo "Running $SCRIPT..."
Rscript "$SCRIPT" > "$LOG" 2>&1

if [[ $? -eq 0 ]]; then
  echo "Success"
else
  echo "Failed. Check $LOG"
fi
