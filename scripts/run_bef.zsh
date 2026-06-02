#!/bin/zsh

# ==========================================
# Usage:
# ./run_bef.zsh . processed_data/plot_biomass.txt "ft1.forest_type,sp_code" 4 4 4000 0.99 15 rsd
# ./run_bef.zsh . processed_data/plot_biomass.txt "PFT,sp_code" 4 4 4000 0.99 15 rsd
# ==========================================

WD=${1:-$(pwd)}
DATA_FILE=${2:-processed_data/plot_biomass.txt}
HIERARCHY=${3:-PFT,sp_code}
CHAINS=${4:-4}
CORES=${5:-4}
ITER=${6:-4000}
ADAPT=${7:-0.99}
TREE=${8:-15}
XVAR=${9:-rsd}

DATE=$(date +%Y%m%d_%H%M%S)

cd "$WD" || {
    echo "Cannot access working directory:"
    echo "$WD"
    exit 1
}

mkdir -p scripts/logfiles


if [ ! -f "$DATA_FILE" ]; then
    echo "ERROR: Data file does not exist:"
    echo "$DATA_FILE"
    exit 1
fi

echo "Data file found:"
echo "$DATA_FILE"

echo "Data file size:"
ls -lh "$DATA_FILE"

echo "First line / column names:"
head -n 1 "$DATA_FILE" | tr '\t' '\n' | nl -ba

SAFE_HIERARCHY=$(echo "$HIERARCHY" | tr ',' '_' | sed 's/[^A-Za-z0-9_]/_/g')
SAFE_DATA=$(basename "$DATA_FILE" | sed 's/[^A-Za-z0-9_]/_/g')

LOG_FILE="scripts/logfiles/log_BEF_${SAFE_DATA}_${SAFE_HIERARCHY}_${XVAR}_${DATE}.txt"

echo "====================================="
echo "Working Dir   : $WD"
echo "Data file     : $DATA_FILE"
echo "Hierarchy     : $HIERARCHY"
echo "Chains        : $CHAINS"
echo "Cores         : $CORES"
echo "Iterations    : $ITER"
echo "Adapt Delta   : $ADAPT"
echo "Max TreeDepth : $TREE"
echo "X variable    : $XVAR"
echo "Log file      : $LOG_FILE"
echo "====================================="

nohup zsh -c "
Rscript scripts/run_bef_zsh.R \
  '$DATA_FILE' \
  '$HIERARCHY' \
  '$CHAINS' \
  '$CORES' \
  '$ITER' \
  '$ADAPT' \
  '$TREE' \
  '$XVAR' \
  2>&1 | while IFS= read -r line; do
    printf '[%s] %s\n' \"\$(date '+%Y-%m-%d %H:%M:%S')\" \"\$line\"
  done
" > "$LOG_FILE" 2>&1 &

echo ""
echo "Started."
echo "PID  : $!"
echo "LOG  : $WD/$LOG_FILE"
echo ""
echo "Monitor:"
echo "tail -f $WD/$LOG_FILE"