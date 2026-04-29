#!/usr/bin/env bash
# Smoke test for the Eigen fetch. Header-only, so all we check is that
# fetch_eigen.sh populated the expected include tree.

set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
CACHE="${REPO_ROOT}/.wasm-cache/eigen"

fail() { echo "FAIL: $1" >&2; exit 1; }

[[ -f "${CACHE}/include/Eigen/Dense" ]]   || fail "missing ${CACHE}/include/Eigen/Dense"
[[ -f "${CACHE}/include/Eigen/QR" ]]      || fail "missing ${CACHE}/include/Eigen/QR"
[[ -f "${CACHE}/include/Eigen/SVD" ]]     || fail "missing ${CACHE}/include/Eigen/SVD"
[[ -f "${CACHE}/include/Eigen/Core" ]]    || fail "missing ${CACHE}/include/Eigen/Core"

echo "OK: Eigen cache looks valid at ${CACHE}"
