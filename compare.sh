#!/bin/bash
# Launch the TASOPT aircraft comparison dashboard (Pluto notebook).
#
# Usage:
#   ./compare.sh                        # opens dashboard in browser
#   ./compare.sh ac1.jld2 ac2.jld2     # pre-set file paths via env vars
#
# The Pluto server binds to port 1234, which Codespaces will forward
# and open automatically.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NOTEBOOK="$SCRIPT_DIR/example/compare_aircraft.jl"

if [ ! -f "$NOTEBOOK" ]; then
  echo "Error: notebook not found at $NOTEBOOK"
  exit 1
fi

echo "Starting TASOPT Comparison Dashboard..."
echo "  Notebook: $NOTEBOOK"
echo "  Port:     1234"
echo ""
echo "Once Pluto loads, paste your .jld2 file paths into the text fields."
echo "Save aircraft with: TASOPT.quicksave_aircraft(ac, \"myfile.jld2\")"
echo ""

julia --project="$SCRIPT_DIR" -e "
import Pkg
Pkg.add(\"Pluto\")
import Pluto
Pluto.run(
    notebook=\"$NOTEBOOK\",
    port=1234,
    launch_browser=false,
    host=\"0.0.0.0\"
)
"
