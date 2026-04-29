#!/usr/bin/env bash
# Fetch Eigen headers for the WASM build of scikit-bio-binaries.
#
# Eigen is header-only: there is no compile step. We simply download a
# pinned release tarball and extract the Eigen/ subtree into
# .wasm-cache/eigen/include/. Idempotent: cache hit returns immediately.
#
# Pinned to Eigen 3.4.0 (stable release, ABI-compatible with skbb usage).
#
# Inputs (env, all optional):
#   EIGEN_VERSION   release tag (default: 3.4.0)
#   FORCE_REFRESH=1 wipe the cache before fetching

set -euo pipefail

EIGEN_VERSION="${EIGEN_VERSION:-3.4.0}"
EIGEN_SHA256="${EIGEN_SHA256:-8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72}"
EIGEN_URL="${EIGEN_URL:-https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz}"

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
CACHE="${REPO_ROOT}/.wasm-cache/eigen"
INCLUDE_DIR="${CACHE}/include"

log() { echo "[fetch_eigen] $*"; }

if [[ "${FORCE_REFRESH:-0}" == "1" ]]; then
  log "FORCE_REFRESH=1 — clearing cache"
  rm -rf "${CACHE}"
fi

if [[ -f "${INCLUDE_DIR}/Eigen/Dense" ]]; then
  log "cache hit: ${INCLUDE_DIR}/Eigen/Dense already present — skipping"
  exit 0
fi

mkdir -p "${CACHE}"
TARBALL="${CACHE}/eigen-${EIGEN_VERSION}.tar.gz"

if [[ ! -f "${TARBALL}" ]]; then
  log "downloading ${EIGEN_URL}"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 3 --retry-delay 2 -o "${TARBALL}" "${EIGEN_URL}"
  elif command -v wget >/dev/null 2>&1; then
    wget -q -O "${TARBALL}" "${EIGEN_URL}"
  else
    echo "ERROR: neither curl nor wget available" >&2
    exit 1
  fi
fi

# Verify checksum. Skip verification if EIGEN_SHA256 is explicitly empty.
if [[ -n "${EIGEN_SHA256}" ]]; then
  log "verifying checksum"
  ACTUAL="$(sha256sum "${TARBALL}" | awk '{print $1}')"
  if [[ "${ACTUAL}" != "${EIGEN_SHA256}" ]]; then
    echo "ERROR: sha256 mismatch on ${TARBALL}" >&2
    echo "  expected: ${EIGEN_SHA256}" >&2
    echo "  actual:   ${ACTUAL}" >&2
    rm -f "${TARBALL}"
    exit 1
  fi
fi

log "extracting into ${INCLUDE_DIR}"
mkdir -p "${INCLUDE_DIR}"
TMPDIR_EXTRACT="$(mktemp -d)"
trap 'rm -rf "${TMPDIR_EXTRACT}"' EXIT
tar -xzf "${TARBALL}" -C "${TMPDIR_EXTRACT}"

# Tarball layout: eigen-<version>/Eigen/..., eigen-<version>/unsupported/...
EXTRACTED="${TMPDIR_EXTRACT}/eigen-${EIGEN_VERSION}"
if [[ ! -d "${EXTRACTED}/Eigen" ]]; then
  echo "ERROR: unexpected tarball layout under ${EXTRACTED}" >&2
  ls -la "${TMPDIR_EXTRACT}" >&2
  exit 1
fi

rm -rf "${INCLUDE_DIR}/Eigen" "${INCLUDE_DIR}/unsupported"
cp -r "${EXTRACTED}/Eigen" "${INCLUDE_DIR}/Eigen"
if [[ -d "${EXTRACTED}/unsupported" ]]; then
  cp -r "${EXTRACTED}/unsupported" "${INCLUDE_DIR}/unsupported"
fi

log "done: ${INCLUDE_DIR}/Eigen"
