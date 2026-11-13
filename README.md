# PRED-TMR

[![build](https://github.com/claudepasquier/PRED-TMR/actions/workflows/build.yml/badge.svg)](https://github.com/claudepasquier/PRED-TMR/actions/workflows/build.yml)
[![pages](https://github.com/claudepasquier/PRED-TMR/actions/workflows/pages.yml/badge.svg)](https://github.com/claudepasquier/PRED-TMR/actions/workflows/pages.yml)
[![CI pipeline (Koda)](https://src.koda.cnrs.fr/uoa-biol/PRED-TMR/badges/main/pipeline.svg)](https://src.koda.cnrs.fr/uoa-biol/PRED-TMR/-/pipelines)
[![Pages deployment](https://img.shields.io/website?url=https%3A%2F%2Fclaudepasquier.github.io%2FPRED-TMR%2F&label=online%20app)](https://claudepasquier.github.io/PRED-TMR/)

**▶ Run the app online: <https://claudepasquier.github.io/PRED-TMR/>**

A statistical method for predicting transmembrane segments in protein
sequences, originally described in:

> Pasquier, C., Promponas, V.J., Palaios, G.A., Hamodrakas, J.S. and
> Hamodrakas, S.J. (1999). *A novel method for predicting transmembrane
> segments in proteins based on a statistical analysis of the SwissProt
> database: the PRED-TMR algorithm.* **Protein Engineering** 12(5), 381–385.

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

## History

- **[original: 1999]** PRED-TMR was developed at the **Faculty of Biology,
  Department of Cell Biology and Biophysics, University of Athens,
  Panepistimiopolis, Athens 15701, Greece**, by **Claude Pasquier**,
  **Vasilis J. Promponas**, **Georgios A. Palaios**, **Ioannis S. Hamodrakas**
  and **Stavros J. Hamodrakas** (see [`AUTHORS.md`](AUTHORS.md)).
- **[modernised: 2026]** The modern C++17 re-implementation (validated kernel,
  CLI, WebAssembly portal, CI) was produced at **Université Côte d'Azur, CNRS,
  I3S (Sophia Antipolis), France** by **Claude Pasquier**.

## Historical access

The original 1999 service was hosted at the University of Athens:

- Legacy version accessible and executable:
  `http://athina.biol.uoa.gr/PRED-TMR/`

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
├── docs/ALGORITHM.md       # full description of the PRED-TMR method
├── CITATION.cff            # machine-readable citation metadata
├── AUTHORS.md              # original authors + modern port credits
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

## How it works

PRED-TMR refines a standard hydrophobicity analysis with the detection of
favourable terminal ("start"/"end") configurations of transmembrane regions.
A statistically derived hydrophobicity scale (`H_i = ln(P_i)`, where `P_i` is
the propensity of residue type *i* to occur inside a transmembrane segment)
scores candidate segments, while decapeptide-based left/right indicators
delimit their termini. A length penalty (`e^(-|l-21|)`) favours segments of
biologically plausible length (94% of observed segments are 17–25 residues).
Segments are selected greedily by decreasing maximum score, rejecting
overlapping candidates. See [`docs/ALGORITHM.md`](docs/ALGORITHM.md) for the
full description and [`include/predtmr/`](include/predtmr/) for the
implementation.

## Citation

If you use PRED-TMR in your research, please cite:

> Pasquier, C., Promponas, V.J., Palaios, G.A., Hamodrakas, J.S. and
> Hamodrakas, S.J. (1999). A novel method for predicting transmembrane
> segments in proteins based on a statistical analysis of the SwissProt
> database: the PRED-TMR algorithm. *Protein Engineering*, 12(5), 381–385.

Machine-readable citation metadata is provided in
[`CITATION.cff`](CITATION.cff); author information in
[`AUTHORS.md`](AUTHORS.md).