#!/bin/bash
# Regenerate Doxygen HTML into docs/doxygen/html/
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

DOXYFILE="docs/Doxyfile_Cleaned"
LOG_DIR="docs/logs"
LOG_FILE="$LOG_DIR/doxygen_log.txt"
mkdir -p "$LOG_DIR"

if ! command -v doxygen >/dev/null 2>&1; then
  if [ -x "$HOME/.local/bin/doxygen" ]; then
    export PATH="$HOME/.local/bin:$PATH"
  else
    echo "Error: doxygen not found on PATH. Install doxygen or place a binary in ~/.local/bin/"
    exit 1
  fi
fi

echo "=== CFD Solver documentation update ==="
echo "Config: $DOXYFILE"
echo "INPUT: src include CUDA_KERNELS Basic_Function_Files Incompressible_Solver Metal_Kernels docs/DoxygenMainPage.cpp"
echo "Running: $(command -v doxygen) ($(doxygen --version))"

doxygen "$DOXYFILE" 2>&1 | tee "$LOG_FILE"

if [ -f "docs/doxygen/html/index.html" ]; then
  echo "✅ Documentation generated: docs/doxygen/html/index.html"
  file_count=$(grep -c "Parsing file" "$LOG_FILE" || true)
  echo "📄 Files processed: ${file_count:-0}"
  echo "📊 Log: $LOG_FILE"
else
  echo "❌ Documentation was not generated — see $LOG_FILE"
  tail -20 "$LOG_FILE" || true
  exit 1
fi
