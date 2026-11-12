# PRED-TMR

[![build](https://github.com/claudepasquier/PRED-TMR/actions/workflows/build.yml/badge.svg)](https://github.com/claudepasquier/PRED-TMR/actions/workflows/build.yml)
[![pages](https://github.com/claudepasquier/PRED-TMR/actions/workflows/pages.yml/badge.svg)](https://github.com/claudepasquier/PRED-TMR/actions/workflows/pages.yml)

A modern C++17 re-implementation of the original **PRED-TMR** transmembrane
segment predictor, with a pure, validated kernel shared by three targets:

- a native **CLI** (`predtmr`),
- the **unit tests** (doctest, run through CTest),
- a **WebAssembly** module (`predtmr.js` + `predtmr.wasm`) powering a static
  web portal that runs entirely in the browser.

The kernel (`src/core/`) is a faithful port of the historical `MarkTM()`
segmentation and `rating()` routines.

## Access URL

The static WASM portal is deployed via GitHub Pages by the workflow
[`.github/workflows/pages.yml`](.github/workflows/pages.yml).

  ```
  https://claudepasquier.github.io/PRED-TMR/
  ```

- **Local portal preview** — build the WASM module then serve `web/`:

  ```bash
  scripts/build-wasm.sh
  python3 -m http.server 8000 --directory web
  # then open http://localhost:8000
  ```

- **Native CLI (local)** — the binary runs locally via `build/predtmr`
  (see "Native build" / "CLI usage" below).

## Layout

```
PRED-TMR/
├── CMakeLists.txt          # native + WASM targets (PREDTMR_BUILD_WASM)
├── include/predtmr/        # public headers (indicators, predictor, sequence, parser)
├── src/
│   ├── core/               # validated kernel (do not modify)
│   ├── parser/             # FASTA / SwissProt / PDB / plain parsing
│   ├── cli/                # native command-line front-end
│   └── wasm/bindings.cpp   # embind bindings (phase B)
├── tests/                  # doctest unit tests
├── scripts/build-wasm.sh   # one-shot Emscripten build
└── web/                    # static portal + generated predtmr.js / predtmr.wasm
```

## Native build

Requires CMake ≥ 3.20 and a C++17 compiler.

```bash
cmake -S . -B build -DPREDTMR_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## CLI usage

```bash
# FASTA from a file
./build/predtmr tests/data/sample.fasta

# raw sequence from stdin, per-residue output
echo "MNSVPLVLAAVALLALAAALAAAAPGGGAAA" | ./build/predtmr --format plain --output L
```

Options:

| Option | Description |
| --- | --- |
| `--format <auto\|plain\|fasta\|swissprot\|pdb>` | Input format (default: auto-detect) |
| `--output <text\|L>` / `--per-residue` | Human report or per-residue flags |
| `-h`, `--help` | Show help |

## WebAssembly build

The WASM module is produced by Emscripten with `--bind` (embind),
`MODULARIZE=1`, `EXPORT_NAME=createPredTMR` and `ALLOW_MEMORY_GROWTH=1`. Output
(`predtmr.js` + `predtmr.wasm`) is written directly into `web/`.

### Prerequisites

The [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html)
must be available. Install it with:

```bash
git clone https://github.com/emscripten-core/emsdk.git /tmp/emsdk
cd /tmp/emsdk && ./emsdk install latest && ./emsdk activate latest
```

`scripts/build-wasm.sh` auto-detects an existing SDK via `emcc` on `PATH`, then
`$EMSDK`, `/tmp/emsdk`, `/opt/emsdk`, `$HOME/emsdk` or `../emsdk`.

### Build

```bash
scripts/build-wasm.sh
# or, if emsdk lives elsewhere:
EMSDK=/path/to/emsdk scripts/build-wasm.sh
```

Under the hood this runs:

```bash
emcmake cmake -S . -B build-wasm \
    -DPREDTMR_BUILD_WASM=ON -DPREDTMR_BUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build-wasm --target predtmr_wasm --parallel
```

### JavaScript API

```js
const module = await createPredTMR();

// Parse FASTA / SwissProt / PDB / raw residues with the shared C++ parser.
const { name, sequence, length, format } =
    module.parseSequence(">p\nMNSVPLVLAAVALLALAAALAAAAPGGGAAA", "auto");

// Run the kernel.
const { segments, perResidue } = module.predict(sequence);
// segments: [{ start, end, score }, ...]   perResidue: [0|1, ...] (length N)
```

## Web portal (local preview)

After building the WASM module, serve `web/` over HTTP (browsers block
`fetch` of `.wasm` over `file://`):

```bash
python3 -m http.server 8000 --directory web
# open http://localhost:8000
```

Paste a raw sequence or a FASTA record, press **Predict**, and inspect the
highlighted sequence and the segment table. Ships a few built-in examples
(bacteriorhodopsin, a synthetic single-pass TM, and a soluble protein).