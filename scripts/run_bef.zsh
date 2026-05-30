#!/bin/zsh

# ==========================================
# Usage:
#
# ./run_bef.zsh \
#   . \
#   ft1.forest_type \
#   4 \
#   4000 \
#   4 \
#   0.99 \
#   15
#
# ==========================================

WD=${1:-$(pwd)}
GROUP=${2:-ft1.forest_type}
CHAINS=${3:-4}
ITER=${4:-4000}
CORES=${5:-4}
ADAPT=${6:-0.99}
TREE=${7:-15}

DATE=$(date +%Y%m%d_%H%M%S)

LOG_FILE="scripts/logfiles/log_BEF_${GROUP}_${DATE}.txt"

cd "$WD" || {
    echo "Cannot access working directory:"
    echo "$WD"
    exit 1
}

echo "====================================="
echo "Working Dir : $WD"
echo "Group       : $GROUP"
echo "Chains      : $CHAINS"
echo "Iterations  : $ITER"
echo "Cores       : $CORES"
echo "Adapt Delta : $ADAPT"
echo "Tree Depth  : $TREE"
echo "====================================="

nohup Rscript scripts/run_bef_zsh.R \
    "$GROUP" \
    "$CHAINS" \
    "$ITER" \
    "$CORES" \
    "$ADAPT" \
    "$TREE" \
    > "$LOG_FILE" 2>&1 &

echo ""
echo "Started."
echo "PID  : $!"
echo "LOG  : $WD/$LOG_FILE"
echo ""

echo "Monitor:"
echo "tail -f $WD/$LOG_FILE"
