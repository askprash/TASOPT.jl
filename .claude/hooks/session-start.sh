#!/bin/bash
# session-start.sh
# Installs Julia via conda-forge so it's available in Claude Code web sessions.
# Idempotent – safe to run multiple times; after the first run the container
# state is cached, so subsequent sessions start instantly.

set -euo pipefail

# Only run in remote (web) sessions
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
    exit 0
fi

MINIFORGE_DIR="/opt/miniforge3"
JULIA_BIN="${MINIFORGE_DIR}/bin/julia"

echo "[session-start] Starting Julia environment setup..."

# ── 1. Install Miniforge (conda-forge) if not present ─────────────────────────
if [ ! -f "${MINIFORGE_DIR}/bin/conda" ]; then
    echo "[session-start] Installing Miniforge..."
    TMPSH="$(mktemp /tmp/miniforge-XXXXX.sh)"
    curl -fsSL -o "${TMPSH}" \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    bash "${TMPSH}" -b -p "${MINIFORGE_DIR}" > /tmp/miniforge-install.log 2>&1
    rm -f "${TMPSH}"
    echo "[session-start] Miniforge installed."
else
    echo "[session-start] Miniforge already present – skipping."
fi

# ── 2. Install Julia via conda-forge if not present ───────────────────────────
if [ ! -f "${JULIA_BIN}" ]; then
    echo "[session-start] Installing Julia via conda-forge..."
    "${MINIFORGE_DIR}/bin/conda" install -n base -c conda-forge julia -y \
        > /tmp/julia-conda-install.log 2>&1
    echo "[session-start] Julia installed: $("${JULIA_BIN}" --version)"
else
    echo "[session-start] Julia already present: $("${JULIA_BIN}" --version) – skipping."
fi

# ── 3. Persist PATH and JULIA_PKG_SERVER for the Claude session ───────────────
# JULIA_PKG_SERVER="" bypasses the blocked pkg.julialang.org endpoint
# and forces Julia to download packages directly from GitHub.
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
    echo "export PATH=\"${MINIFORGE_DIR}/bin:\$PATH\""   >> "${CLAUDE_ENV_FILE}"
    echo "export JULIA_PKG_SERVER=\"\""                   >> "${CLAUDE_ENV_FILE}"
fi
export PATH="${MINIFORGE_DIR}/bin:${PATH}"
export JULIA_PKG_SERVER=""

# ── 4. Bootstrap the Julia General registry (first run only) ──────────────────
if [ ! -d "${HOME}/.julia/registries/General" ]; then
    echo "[session-start] Adding Julia General registry from GitHub..."
    "${JULIA_BIN}" -e '
        import Pkg
        Pkg.Registry.add("General")
        println("[session-start] Registry added.")
    ' 2>&1 | grep -v "^┌ Warning\|^│\|^└"
fi

# ── 5. Instantiate TASOPT project packages ────────────────────────────────────
PROJECT_DIR="${CLAUDE_PROJECT_DIR:-/home/user/TASOPT.jl}"
if [ -f "${PROJECT_DIR}/Project.toml" ]; then
    echo "[session-start] Instantiating TASOPT packages (this may take a few minutes on first run)..."
    "${JULIA_BIN}" --project="${PROJECT_DIR}" -e '
        import Pkg
        Pkg.instantiate()
        println("[session-start] TASOPT packages ready.")
    ' 2>&1 | grep -v "^┌ Warning\|^│\|^└" | grep -v "^$"
fi

echo "[session-start] Julia environment ready. Julia: $("${JULIA_BIN}" --version)"
