# Authors

## Original PRED-TMR algorithm and publication

PRED-TMR was developed at the Faculty of Biology, Department of Cell Biology
and Biophysics, University of Athens, Panepistimiopolis, Athens 15701, Greece.

- **Claude Pasquier**
- **Vasilis J. Promponas**
- **Georgios A. Palaios**
- **Ioannis S. Hamodrakas**
- **Stavros J. Hamodrakas** (corresponding author, `shamodr@atlas.uoa.gr`)

Published as: Pasquier, C., Promponas, V.J., Palaios, G.A., Hamodrakas, J.S.
and Hamodrakas, S.J. (1999). *A novel method for predicting transmembrane
segments in proteins based on a statistical analysis of the SwissProt
database: the PRED-TMR algorithm.* **Protein Engineering** 12(5), 381–385.

The original work was supported by the EEC-TMR "GENE-QUIZ" grant
ERBFMRXCT960019.

## Modern C++17 port

- **Claude Pasquier** — **[modernised: 2026]** Université Côte d'Azur, CNRS,
  I3S, Sophia Antipolis, France — lead of the modernised edition
  (Software Architecture, Code Modernisation, Writing — Original Draft).
  Original affiliation (**[original: 1999]**): University of Athens (see above).

The modern re-implementation in this repository (validated kernel, CLI,
WebAssembly bindings, web portal and CI) was produced as a faithful port of
the historical `MarkTM()` segmentation and `rating()` routines, without
altering the protein logic.
