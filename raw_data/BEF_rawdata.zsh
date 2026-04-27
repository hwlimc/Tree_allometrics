#!/bin/zsh

SCRIPT="/Users/hyli0001/Library/CloudStorage/OneDrive-Sverigeslantbruksuniversitet/wrd/b/Dynamic_allometrics/original_files/raw_data/kr_biomass_dataset.R"
LOG="run.log"

echo "Running $SCRIPT..."
Rscript "$SCRIPT" > "$LOG" 2>&1

if [[ $? -eq 0 ]]; then
  echo "Success"
else
  echo "Failed. Check $LOG"
fi
