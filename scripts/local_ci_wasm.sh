#!/usr/bin/env bash
# Local mirror of the build-and-test-wasm GitHub Actions job.
#
# Spins up a clean Ubuntu 24.04 container with only the packages the
# CI workflow installs, then runs the same commands the workflow runs.
# Catches the class of bug where a local conda env quietly provides
# headers (lapacke.h, cblas.h, ...) that the CI runner doesn't have.
#
# Usage:
#   scripts/local_ci_wasm.sh
#
# Prerequisites: docker daemon running, internet access for apt + emsdk.

set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

echo "[local_ci_wasm] mounting ${REPO_ROOT} -> /source (ro) in ubuntu:24.04"

# Mount the repo READ-ONLY at /source and have the container rsync it
# into a private /work directory. This keeps build artifacts container-
# local so the host filesystem is never polluted by docker-as-root files.
# -t for unbuffered stdout. --rm cleans up the container after exit.
docker run --rm -t \
  --network=host \
  -v "${REPO_ROOT}:/source:ro" \
  -e DEBIAN_FRONTEND=noninteractive \
  ubuntu:24.04 \
  bash -ec '
    set -eu
    export EMSDK_QUIET=1

    echo "=== container identity ==="
    id
    uname -a
    cat /etc/os-release | head -3

    # Copy the repo into a private writable location, excluding build
    # artifacts and the .wasm-cache (force re-download in the container
    # so we test the genuine bootstrap, not the host-cached state).
    mkdir -p /work
    cp -a /source/. /work/
    rm -rf /work/.wasm-cache
    cd /work
    echo "=== /work contents ==="
    ls -la
    echo "=========================="

    # Step 5 in main.yml — "Install native build deps for expected-value
    # generators". Keep this list in sync with the apt-get install line
    # in .github/workflows/main.yml; if they drift, this harness no
    # longer mirrors CI and will silently mask CI failures.
    #
    # Required packages (see workflow yaml comment for details):
    #   libopenblas-dev   — libopenblas.so (cblas + LAPACK Fortran)
    #   libblas-dev       — cblas.h at multiarch path
    #   liblapacke-dev    — lapacke.h + liblapacke.so
    apt-get update -qq
    apt-get install -y --no-install-recommends \
        ca-certificates curl git python3 xz-utils \
        g++ make libopenblas-dev libblas-dev liblapacke-dev

    # Step 3 in main.yml — emscripten-core/setup-emsdk@v16 with version 5.0.3.
    if [ ! -d /opt/emsdk ]; then
      git clone --depth 1 https://github.com/emscripten-core/emsdk.git /opt/emsdk
    fi
    /opt/emsdk/emsdk install 5.0.3
    /opt/emsdk/emsdk activate 5.0.3
    # shellcheck source=/dev/null
    . /opt/emsdk/emsdk_env.sh

    # Workaround for an emsdk + GNU make interaction:
    # emsdk_env.sh prepends /opt/emsdk to PATH, and that directory
    # contains a `node` SUBDIRECTORY (the toolchain root for the
    # bundled node). When make does a PATH lookup for "node" via its
    # direct-execve fast path, it finds that directory FIRST and tries
    # to exec it, which returns EACCES ("Permission denied"). Bash
    # filters out non-regular files during PATH lookup and avoids this,
    # which is why bash invocations work but make recipes do not.
    #
    # CI does not hit this because actions/setup-node@v4 prepends a
    # real node binary path. We mirror that by putting the actual node
    # bin directory before /opt/emsdk in PATH.
    export PATH="$(dirname "${EMSDK_NODE}"):${PATH}"
    node --version

    # Steps 7/8 — fetch + verify Eigen.
    scripts/fetch_eigen.sh
    scripts/test_eigen_wasm.sh

    # Step 9 — build the WASM archive.
    make wasm

    echo "=== node sanity inside script context ==="
    which node
    ls -la "$(which node)"
    node --version
    /bin/sh -c "node --version"
    echo "=========================================="

    # Steps 10/11 — run WASM tests with the env CI uses. Keep the
    # BLASLIB override here in sync with .github/workflows/main.yml.
    env NOGPU=1 BLASLIB="-llapacke -lopenblas" make wasm_test
    env NOGPU=1 BLASLIB="-llapacke -lopenblas" make wasm_api_test

    echo
    echo "[local_ci_wasm] PASSED — local CI mirror is green"
  '
