#!/bin/zsh

unsetopt BG_NICE 2>/dev/null || true

# Usage:
# ./run_bef.zsh WD DATA_FILE HIERARCHY CHAINS CORES ITER ADAPT TREE XVAR MODEL K_DEPTH FAMILY SCALE_X SPLIT_COL SPLIT_VALUES DROP_SPLIT
#
# MODEL is a shape name from scripts/00_model_shapes.R, for example exp_decay,
# michaelis_menten, or linear. FAMILY can be gamma, lognormal/lnorm,
# student/tdis, or gaussian/normal/ndis.
#
# Combined model:
# ./run_bef.zsh . processed_data/plot_biomass.txt "PFT,sp_code" 4 4 4000 0.99 15 rsd exp_decay 1 gamma FALSE none all TRUE
#
# Separate species models:
# ./run_bef.zsh . processed_data/plot_biomass.txt none 4 4 4000 0.99 15 rsd exp_decay 0 gamma FALSE sp_code all TRUE

WD=${1:-$(pwd)}
DATA_FILE=${2:-processed_data/plot_biomass.txt}
HIERARCHY=${3:-PFT,sp_code}
CHAINS=${4:-4}
CORES=${5:-4}
ITER=${6:-4000}
ADAPT=${7:-0.99}
TREE=${8:-15}
XVAR=${9:-rsd}
MODEL=${10:-exp_decay}
K_DEPTH=${11:-1}
FAMILY=${12:-gamma}
SCALE_X=${13:-FALSE}
SPLIT_COL=${14:-none}
SPLIT_VALUES=${15:-all}
DROP_SPLIT=${16:-TRUE}

safe_name() {
	printf '%s' "$1" | tr ',' '-' | sed 's/[^A-Za-z0-9_-]/-/g; s/-\{1,\}/-/g; s/^-//; s/-$//; s/^$/none/'
}

print_field() {
	printf '%-16s %s\n' "$1:" "$2"
}

cd "$WD" || {
	echo "Cannot access working directory: $WD"
	exit 1
}

if [ ! -f "$DATA_FILE" ]; then
	echo "ERROR: Data file does not exist: $DATA_FILE"
	exit 1
fi

mkdir -p scripts/logfiles

DATE=$(date +%Y%m%d_%H%M%S)
LOG_FILE="scripts/logfiles/log_BEF_$(safe_name "$MODEL")_$(safe_name "$SPLIT_COL")_$(safe_name "$SPLIT_VALUES")_$(safe_name "$(basename "$DATA_FILE")")_$(safe_name "$HIERARCHY")_$(safe_name "$XVAR")_${DATE}.txt"

echo "Data file found: $DATA_FILE"
ls -lh "$DATA_FILE"
echo "Column names:"
head -n 1 "$DATA_FILE" | tr '\t' '\n' | nl -ba

echo "====================================="
print_field "Working Dir" "$WD"
print_field "Data file" "$DATA_FILE"
print_field "Hierarchy" "$HIERARCHY"
print_field "Chains" "$CHAINS"
print_field "Cores" "$CORES"
print_field "Iterations" "$ITER"
print_field "Adapt Delta" "$ADAPT"
print_field "Max TreeDepth" "$TREE"
print_field "X variable" "$XVAR"
print_field "Model" "$MODEL"
print_field "k depth" "$K_DEPTH"
print_field "Family" "$FAMILY"
print_field "Scale x" "$SCALE_X"
print_field "Split column" "$SPLIT_COL"
print_field "Split values" "$SPLIT_VALUES"
print_field "Drop split" "$DROP_SPLIT"
print_field "Log file" "$LOG_FILE"
echo "====================================="

cmd=(
	Rscript scripts/run_bef_zsh.R
	"$DATA_FILE"
	"$HIERARCHY"
	"$CHAINS"
	"$CORES"
	"$ITER"
	"$ADAPT"
	"$TREE"
	"$XVAR"
	"$MODEL"
	"$K_DEPTH"
	"$FAMILY"
	"$SCALE_X"
	"$SPLIT_COL"
	"$SPLIT_VALUES"
	"$DROP_SPLIT"
)

(
	"${cmd[@]}" 2>&1 | while IFS= read -r line; do
		printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$line"
	done
) > "$LOG_FILE" 2>&1 &

pid=$!
disown "$pid" 2>/dev/null || true

echo ""
echo "Started."
echo "PID  : $pid"
echo "LOG  : $WD/$LOG_FILE"
echo ""
echo "Monitor:"
echo "tail -f $WD/$LOG_FILE"
