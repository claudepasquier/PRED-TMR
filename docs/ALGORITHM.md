# The PRED-TMR algorithm

This document summarises the PRED-TMR method as published in:

> Pasquier, C., Promponas, V.J., Palaios, G.A., Hamodrakas, J.S. and
> Hamodrakas, S.J. (1999). *A novel method for predicting transmembrane
> segments in proteins based on a statistical analysis of the SwissProt
> database: the PRED-TMR algorithm.* **Protein Engineering** 12(5), 381–385.

The modern C++ kernel in `src/core/` is a faithful port of these routines.

## Overview

PRED-TMR predicts transmembrane (TM) segments using only information
contained in the protein sequence itself. It refines a standard
hydrophobicity analysis with the detection of favourable terminal
configurations ("edges": starts and ends) of transmembrane regions. This
allows the method both to discard highly hydrophobic stretches that are not
delimited by clear start/end configurations, and to confirm putative TM
segments that are not distinguishable by their hydrophobic composition alone.

## Statistical background

The method is based on a statistical study of the SwissProt database
(release 35). 9,392 transmembrane proteins were automatically extracted via
the `TRANSMEM` feature-table keyword and stored in the DB-TMR database,
containing for each segment: the accession code, organism classification,
segment length, direction when deducible, five flanking residues on each
side, and the segment residues themselves. Segments extending beyond
sequence ends or with unknown end-points were discarded before statistics
were computed.

### Segment length distribution

The 40,548 segments with reliable end-points have an average length of
21.30 residues (standard deviation 2.56). The distribution is sharper than a
Gaussian; 60% of segments are exactly 21 residues and 94% are 17–25 residues
long. It is approximated by `f(l) = e^(-|l - 21|)`.

### Transmembrane propensities

The propensity of residue type *i* to occur inside a TM segment is

```
P_i = F_i^TM / F_i
```

where `F_i^TM` and `F_i` are the frequencies of residue type *i* inside TM
segments and in the entire SwissProt database. Values > 1 indicate a
preference for the lipid-associated structure; values < 1 indicate
unfavourable residues.

### Statistically derived "hydrophobicity" scale

Following Sipos and von Heijne (1993), the propensities are translated into
a hydrophobicity scale:

```
H_i = ln(P_i)
```

The hydrophobicity of a segment from position *m* to *p* is the sum of the
`H` values of its residues:

```
H_m^p = Σ_{k=m..p} H(R_k)
```

## Terminal (edge) indicators

Favourable configurations are computed for decapeptides centred at TM
segment borders (five residues outside, five inside the membrane; positions
0–9 in the decapeptide). The propensity of residue type *i* at decapeptide
position *p* is `P_i^p = F_i^p / F_i`.

For the N-terminal ("left") side, the propensity that position *p* is the
first residue of the lipid-associated structure is:

```
P_left^p = Σ_{k=-5..4} ln[ P_{R(p+k)}^{1-k} ]
```

For the C-terminal ("right") side, the propensity that position *p* is the
first residue outside the TM region is:

```
P_right^p = Σ_{k=-5..4} ln[ P_{R(p+k)}^{4-k} ]
```

Because some decapeptides score high on both sides (e.g. `ILFVSTFFTM`:
P_left = 1.75, P_right = 2.61), the difference between the two propensities
matters more than their absolute values. They are combined into start and
end indicators:

```
             P_left^p + min(P_left^p, P_left^p - P_right^p)
LeftInd_p  = ----------------------------------------------
                                  2

             P_right^p + min(P_right^p, P_right^p - P_left^p)
RightInd_p = ------------------------------------------------
                                  2
```

The minimum prevents a small opposing propensity from contributing more
than the matching one.

## Scoring

A well-defined TM region should score well on all three parameters
(`LeftInd`, `RightInd`, `H`), but in practice weak indicators are often
compensated by excellent values of the others. A length penalty filters
segments of improbable length:

```
LP_l = e^(-|l - 21|)
```

After normalising the hydrophobicity to a decapeptide
(`NH_m^p = 10 · H_m^p / l`), the score of a segment from *m* to *p* is:

```
Score_m^p = e^(LeftInd_m) + e^(NH_m^p) + e^(RightInd_p) + 1 - LP_l
```

where `l = p - m + 1`.

## Prediction algorithm

1. For each position *m*, compute the maximum score `MScore_m =
   max_p(Score_m^p)` over end positions `p` in `[m+1, m+40]`, keeping only
   segments with positive indicators (`LeftInd > 0`, `RightInd > 0`) and a
   normalised hydrophobicity above the cut-off (`NH_m^p > 2`; the cut-off
   was optimised on a subset of 64 reliable proteins).
2. Memorise `MScore_m` and the corresponding end position for each *m*.
3. Select the highest `MScore_m` and mark the corresponding region as
   transmembrane.
4. Repeatedly select the next-highest `MScore_m` whose segment does not
   overlap a previously marked region, until no further non-overlapping
   candidates remain.

## Accuracy

On a test set of 101 non-homologous TM proteins with reliable topologies,
PRED-TMR obtained an agreement factor `Q = 88.83%`, correlation coefficient
`C = 0.80` and segment-match ratio `SM = 0.954`, comparing well with DAS,
PHDhtm, TOPPRED, SOSUI and tmPRED (best correlation coefficient of the six
methods; see Table II of the paper). Applied to all TM proteins of SwissProt
release 35, accuracy decreased only slightly (`Q = 86.14`, `C = 0.73`,
`SM = 0.889`).

## References

- Sipos, L. and von Heijne, G. (1993) *Eur. J. Biochem.*, 213, 1333–1340.
- von Heijne, G. (1992) *J. Mol. Biol.*, 225, 487–494.
- Rost, B., Casadio, R., Fariselli, P. and Sander, C. (1995) *Protein Sci.*,
  4, 521–533.
- Cserzo, K., Wallin, E., Simon, I., von Heijne, G. and Elofsson, A. (1997)
  *Protein Engng*, 10, 673–676.
