#!/usr/bin/env bash
#
# Builds the PRED-TMR WebAssembly module and emits predtmr.js + predtmr.wasm
# directly into web/ (the static portal directory).
#
# Requirements: the Emscripten SDK (emsdk). The script looks for an existing
# installation in $EMSDK, /tmp/emsdk, /opt/emsdk, $HOME/emsdk or ../emsdk.
# If none is found it fails with instructions; the GitHub Actions workflow
# installs emsdk automatically (see .github/workflows/pages.yml).
#
# Usage:
#   scripts/build-wasm.sh
#
# Environment overrides:
#   EMSDK   path to an existing emsdk checkout/installation.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
repo_root="$(cd "${project_dir}/.." && pwd)"

emcc_candidate() {
    command -v emcc >/dev/null 2>&1
}

setup_emsdk() {
    if emcc_candidate; then
        echo "==> emcc found on PATH"
        return 0
    fi

    local candidates=(
        "${EMSDK:-}"
        "/tmp/emsdk"
        "/opt/emsdk"
        "${HOME}/emsdk"
        "${repo_root}/emsdk"
    )

    for candidate in "${candidates[@]}"; do
        if [[ -n "${candidate}" && -f "${candidate}/emsdk_env.sh" ]]; then
            echo "==> Activating emsdk at ${candidate}"
            # shellcheck disable=SC1091
            source "${candidate}/emsdk_env.sh" >/dev/null 2>&1
            if emcc_candidate; then
                return 0
            fi
        fi
    done

    echo "ERROR: Emscripten SDK (emsdk) not found." >&2
    echo "Install it with:" >&2
    echo "  git clone https://github.com/emscripten-core/emsdk.git /tmp/emsdk" >&2
    echo "  cd /tmp/emsdk && ./emsdk install latest && ./emsdk activate latest" >&2
    echo "  EMSDK=/tmp/emsdk scripts/build-wasm.sh" >&2
    return 1
}

setup_emsdk

build_dir="${project_dir}/build-wasm"

echo "==> Configuring WASM build (emcmake)"
emcmake cmake -S "${project_dir}" -B "${build_dir}" \
    -DPREDTMR_BUILD_WASM=ON \
    -DPREDTMR_BUILD_TESTS=OFF \
    -DCMAKE_BUILD_TYPE=Release

echo "==> Building predtmr.js / predtmr.wasm"
cmake --build "${build_dir}" --target predtmr_wasm --parallel

echo "==> Artifacts:"
ls -la "${project_dir}/web/predtmr.js" "${project_dir}/web/predtmr.wasm"
