"""Validation script: translate ORFs from the M. genitalium genome (2026 version)
and run the modern PRED-TMR port on them, comparing predicted segments with the
historical reference PRED-TMR/Results/result_mg.txt.

The 2026 genome may differ from the 1998 version used to produce result_mg.txt,
so differences are expected; the goal is to quantify agreement.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENOME = ROOT / "PRED-TMR" / "GCA_040556925.1_ASM4055692v1_genomic.fna"
REFERENCE = ROOT / "PRED-TMR" / "Results" / "result_mg.txt"
BINARY = ROOT / "build" / "predtmr"

CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq: str) -> str:
    return "".join(
        CODON_TABLE.get(seq[i : i + 3], "X") for i in range(0, len(seq) - 2, 3)
    )


def extract_orfs(genome: str) -> list[tuple[str, str]]:
    """Extract ORFs (>= 100 codons) from all 6 frames, named MG-like by index."""
    orfs: list[tuple[str, str]] = []
    strands = [("+", genome), ("-", reverse_complement(genome))]
    counter = 0
    for strand, seq in strands:
        for frame in range(3):
            sub = seq[frame:]
            protein = translate(sub)
            for match in re.finditer(r"M[^*]{99,}", protein):
                counter += 1
                orfs.append((f"ORF{counter:04d}", match.group(0)))
    return orfs


def parse_reference(path: Path) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = {}
    for line in path.read_text().splitlines():
        match = re.match(r"^(MG\d+)\s+(\d+)\s+(.*)$", line)
        if not match:
            continue
        segments = re.findall(r"\((\d+)-(\d+)\)", match.group(3))
        result[match.group(1)] = [(int(a), int(b)) for a, b in segments]
    return result


def main() -> int:
    genome_lines = GENOME.read_text().splitlines()
    genome = "".join(line.strip() for line in genome_lines if not line.startswith(">"))
    orfs = extract_orfs(genome)
    print(f"Genome length: {len(genome)} bp, ORFs extracted: {len(orfs)}")

    fasta = "\n".join(f">{name}\n{seq}" for name, seq in orfs) + "\n"
    tmp = Path("/tmp/mg_orfs.fasta")
    tmp.write_text(fasta)

    proc = subprocess.run(
        [str(BINARY), str(tmp)], capture_output=True, text=True, check=True
    )
    predicted: dict[str, list[tuple[int, int]]] = {}
    current = ""
    for line in proc.stdout.splitlines():
        if line.startswith("Sequence: "):
            current = line.split("Sequence: ", 1)[1]
            predicted[current] = []
        elif line.startswith("  ") and "  " in line:
            parts = line.split()
            if len(parts) == 4 and parts[0].isdigit():
                predicted[current].append((int(parts[1]), int(parts[2])))

    reference = parse_reference(REFERENCE)
    print(f"Reference entries: {len(reference)}, predicted ORFs: {len(predicted)}")

    # Compare: for each reference entry, find best-matching predicted ORF by
    # overlap of segment positions. The 2026 genome uses different gene IDs,
    # so match by sequence length + segment overlap.
    exact = 0
    partial = 0
    matched = 0
    for mg_id, segments in reference.items():
        best = None
        best_score = -1.0
        for name, pred in predicted.items():
            if not pred:
                continue
            score = 0.0
            for seg in segments:
                for pseg in pred:
                    overlap = min(seg[1], pseg[1]) - max(seg[0], pseg[0]) + 1
                    if overlap > 0:
                        score += overlap / min(
                            seg[1] - seg[0] + 1, pseg[1] - pseg[0] + 1
                        )
            score /= max(len(segments), 1)
            if score > best_score:
                best_score = score
                best = (name, pred)
        if best is None:
            continue
        matched += 1
        if best_score >= 0.99:
            exact += 1
        elif best_score >= 0.5:
            partial += 1
        print(
            f"{mg_id}: best={best[0]} overlap={best_score:.2f} "
            f"ref={len(segments)} pred={len(best[1])}"
        )

    print(
        f"\nSummary: matched={matched}/{len(reference)} exact={exact} partial={partial}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
