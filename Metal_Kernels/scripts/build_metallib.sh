#!/usr/bin/env bash
# Compile CFDSolverKernels.metal to default.metallib in ../compiled/
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${SCRIPT_DIR}/.."
OUT="${ROOT}/compiled"
mkdir -p "${OUT}"
xcrun metal -c "${ROOT}/shaders/CFDSolverKernels.metal" -o "${OUT}/CFDSolverKernels.air"
xcrun metallib "${OUT}/CFDSolverKernels.air" -o "${OUT}/CFDSolverKernels.metallib"
echo "Wrote ${OUT}/CFDSolverKernels.metallib"
