#!/usr/bin/env bash
# regression_gate.sh — Regression baseline drift gate for TASOPT.jl
#
# Checks that no regression baseline files have been modified without an
# explicit BASELINE-REGEN marker in the authorizing commit message.
#
# ─────────────────────────────────────────────────────────────────────────────
# PROBLEM
#   Commit e83166a7 (tasopt-j9l.58) silently regenerated test/default_sized.jl
#   and test/fixtures/engine_sweep_baseline.toml as a side-effect of a
#   refactor commit.  The numeric deltas were real (a Pom bug fix was bundled
#   in), but the developer did not notice the baseline had changed.
#
# SOLUTION
#   Any commit that modifies a watched baseline file must contain the line
#
#     BASELINE-REGEN: tasopt-<id>
#
#   in its commit message, where <id> is the beads issue that documents WHY
#   the baseline changed.  This script fails the push if any modifying commit
#   lacks that marker.
#
#   Initial creation of a new baseline file does NOT require the marker
#   (there is nothing upstream to regress against).
#
# ─────────────────────────────────────────────────────────────────────────────
# USAGE
#   ./test/tools/regression_gate.sh [--upstream-ref=<gitref>] [--verbose]
#
#   --upstream-ref=REF   Ref to compare against (default: upstream/main)
#   --verbose            Print each checked commit even when it passes
#
# EXIT CODES
#   0   All baseline changes are authorized (or no changes found)
#   1   Unauthorized baseline drift detected; push is blocked
#
# OPT-OUT
#   In the commit message body (not subject line, but body is also fine):
#
#     BASELINE-REGEN: tasopt-abc123
#
#   The <id> must be a valid beads issue that explains the change.
#   Multiple baselines changed in one commit?  List them on separate lines:
#
#     BASELINE-REGEN: tasopt-abc123
#     BASELINE-REGEN: tasopt-def456
#
# ─────────────────────────────────────────────────────────────────────────────
# WIRING
#   Pre-push hook — install once per clone:
#     cp test/tools/regression_gate.sh .git/hooks/pre-push
#     chmod +x .git/hooks/pre-push
#
#   Or source from an existing pre-push hook:
#     ./test/tools/regression_gate.sh || exit 1
#
#   bd preflight — run manually before opening a PR:
#     ./test/tools/regression_gate.sh
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
UPSTREAM_REF="upstream/main"
VERBOSE=0

for arg in "$@"; do
  case "$arg" in
    --upstream-ref=*)  UPSTREAM_REF="${arg#--upstream-ref=}" ;;
    --verbose)         VERBOSE=1 ;;
    *)  echo "Unknown argument: $arg" >&2; exit 1 ;;
  esac
done

# ── Watched baseline files ─────────────────────────────────────────────────────
# Add new baseline files here as they are created.
BASELINE_FILES=(
  "test/default_sized.jl"
  "test/fixtures/engine_sweep_baseline.toml"
  "test/fixtures/engine_benchmark_baseline.toml"
)

# Marker pattern.  A commit message must contain this (case-sensitive) to be
# considered an approved baseline regen.
#   BASELINE-REGEN: tasopt-<alphanum>
MARKER_PATTERN='BASELINE-REGEN:[[:space:]]+tasopt-[a-zA-Z0-9._]+'

# ── Helpers ────────────────────────────────────────────────────────────────────
RED='\033[0;31m'
YEL='\033[0;33m'
GRN='\033[0;32m'
RST='\033[0m'

info()  { printf "${GRN}[gate]${RST}  %s\n"   "$*"; }
warn()  { printf "${YEL}[gate]${RST}  %s\n"   "$*"; }
error() { printf "${RED}[gate]${RST}  %s\n"   "$*" >&2; }

# ── Resolve upstream ref ───────────────────────────────────────────────────────
if ! git rev-parse --verify "$UPSTREAM_REF" >/dev/null 2>&1; then
  warn "Upstream ref '$UPSTREAM_REF' not found; skipping gate."
  exit 0
fi

merge_base=$(git merge-base HEAD "$UPSTREAM_REF" 2>/dev/null) || {
  warn "Cannot compute merge-base with '$UPSTREAM_REF'; skipping gate."
  exit 0
}

# ── Check each baseline file ───────────────────────────────────────────────────
violations=0

for file in "${BASELINE_FILES[@]}"; do
  # Find all commits between merge-base and HEAD that touched this file.
  commits=$(git log --pretty=format:"%H" "$merge_base..HEAD" -- "$file" 2>/dev/null)
  [[ -z "$commits" ]] && continue

  while IFS= read -r sha; do
    # Is this commit an ADDITION of the file (first-ever creation)?
    # If the file did not exist in the commit's parent, it's a pure add —
    # no upstream regressions are possible, so no marker required.
    parent="${sha}^"
    if ! git cat-file -e "${parent}:${file}" 2>/dev/null; then
      [[ $VERBOSE -eq 1 ]] && info "  $sha — ADD (no marker required): $file"
      continue
    fi

    # File existed in the parent → this is a modification.  Require marker.
    msg=$(git log -1 --format="%s%n%b" "$sha")
    short=$(git log -1 --format="%h %s" "$sha")

    if printf '%s' "$msg" | grep -qE "$MARKER_PATTERN"; then
      [[ $VERBOSE -eq 1 ]] && info "  $sha — approved: $file"
    else
      error "Unauthorized baseline modification:"
      error "  file:   $file"
      error "  commit: $short"
      error ""
      error "  The commit modifies a watched baseline but its message lacks"
      error "  an explicit approval marker.  To approve a baseline change,"
      error "  add a line to the commit message body:"
      error ""
      error "    BASELINE-REGEN: tasopt-<id>"
      error ""
      error "  where <id> is the beads issue documenting why the baseline changed."
      violations=$((violations + 1))
    fi
  done <<< "$commits"
done

# ── Verdict ────────────────────────────────────────────────────────────────────
if [[ $violations -eq 0 ]]; then
  info "Regression gate passed — no unauthorized baseline modifications."
  exit 0
else
  error ""
  error "Regression gate FAILED: $violations unauthorized baseline modification(s)."
  error "Fix by amending the offending commit(s) to include BASELINE-REGEN markers,"
  error "or by reverting the baseline changes if they were unintentional."
  exit 1
fi
