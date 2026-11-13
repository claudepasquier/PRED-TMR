# Changelog

## 1.0.0 — 2026-11-12

First release of the modernised PRED-TMR (C++17 edition).

### Added
- `src/core/`: validated C++17 kernel, a faithful port of the historical
  `MarkTM()` segmentation and `rating()` routines, without altering the
  protein logic.
- `src/parser/`: FASTA / SwissProt / PDB / plain sequence parsing.
- Native command-line interface `predtmr` (`src/cli/`) with `--format` and
  `--output`/`--per-residue` options.
- doctest unit tests (`tests/`) run through CTest.
- WebAssembly module (`src/wasm/bindings.cpp`, embind, MODULARIZE) powering a
  static web portal (`web/`) that runs entirely in the browser.
- `scripts/build-wasm.sh`: one-shot Emscripten build with SDK auto-detection.
- `docs/ALGORITHM.md`: full description of the PRED-TMR method.
- `validate_mg.py`: validation script comparing predictions against the
  historical reference results.
- CI workflows: GitHub Actions (`build.yml`, `pages.yml`) and CNRS Koda
  GitLab (`.gitlab-ci.yml`) running the same native build + test gates.
- Zenodo metadata, CITATION.cff with the original 1999 Protein Engineering
  paper, GPL-2.0-or-later license, AUTHORS.md.

### Changed
- All computation for the web portal now runs client-side (WebAssembly); no
  server, no CGI, no user data stored — sequences never leave the browser.
