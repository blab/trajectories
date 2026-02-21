# Data Format

This document explains the two trajectory formats produced by the pipeline: **forwards trajectories** (root-to-tip) and **pairwise trajectories** (tip-to-tip).

## Worked example

Consider a small phylogenetic tree with 3 tips (A, B, C) and 2 internal nodes (X, Y), where X is the root. Each node has a 10-nucleotide sequence. Branch labels show the Hamming distance (number of nucleotide differences) between parent and child:

```
         ┌──── 2 ─── A
   ┌─ 1 ─Y
───X     └─ 1 ─ B
   │
   └────── 3 ────── C
```

Node sequences:

```
Pos: 1  2  3  4  5  6  7  8  9  10
X:   A  T  C  G  A  T  C  G  A  T     (root)
Y:   A  T  C  A  A  T  C  G  A  T     1 change from X (pos 4: G→A)
A:   A  T  C  A  A  G  C  A  A  T     2 changes from Y (pos 6: T→G, pos 8: G→A)
B:   A  T  C  A  A  T  C  G  G  T     1 change from Y (pos 9: A→G)
C:   A  G  C  G  G  T  C  G  A  C     3 changes from X (pos 2: T→G, pos 5: A→G, pos 10: T→C)
```

## Forwards trajectories

A forwards trajectory traces the evolutionary path from root to a single tip. Each entry in the FASTA file records a node along the path and its cumulative Hamming distance from the root. The header format is `>{node_name}|{cumulative_hamming_distance}`. Intermediate nodes with zero branch distance are skipped.

### Example: trajectory for tip A

Path: X → Y → A

File `A.fasta`:
```
>X|0
ATCGATCGAT
>Y|1
ATCAATCGAT
>A|3
ATCAAGCAAT
```

- `X|0` — root, cumulative distance 0
- `Y|1` — 1 mutation accumulated (X→Y distance is 1)
- `A|3` — 3 mutations accumulated (1 from X→Y + 2 from Y→A)

### Example: trajectory for tip B

Path: X → Y → B

File `B.fasta`:
```
>X|0
ATCGATCGAT
>Y|1
ATCAATCGAT
>B|2
ATCAATCGGT
```

### Example: trajectory for tip C

Path: X → C

File `C.fasta`:
```
>X|0
ATCGATCGAT
>C|3
AGCGGTCGAC
```

Node Y does not appear because tip C descends directly from root X.

### File naming

Each forwards trajectory is saved as `{tip_name}.fasta` inside a sharded tar.zst archive (e.g. `forwards-train-000.tar.zst`).

## Pairwise trajectories

A pairwise trajectory contains exactly two tip sequences. The first tip is labeled with `|0` and the second with `|{hamming_distance}`, where the Hamming distance is computed directly between the two tip sequences (ignoring positions with gaps `-` or ambiguous bases `N`).

### Example: pair A and B

```
Pos: 1  2  3  4  5  6  7  8  9  10
A:   A  T  C  A  A  G  C  A  A  T
B:   A  T  C  A  A  T  C  G  G  T
                    *     *  *
```

Differences at positions 6, 8, 9 — Hamming distance of 3.

File `A__B.fasta`:
```
>A|0
ATCAAGCAAT
>B|3
ATCAATCGGT
```

### Example: pair A and C

```
Pos: 1  2  3  4  5  6  7  8  9  10
A:   A  T  C  A  A  G  C  A  A  T
C:   A  G  C  G  G  T  C  G  A  C
        *     *  *  *     *     *
```

Differences at positions 2, 4, 5, 6, 8, 10 — Hamming distance of 6.

File `A__C.fasta`:
```
>A|0
ATCAAGCAAT
>C|6
AGCGGTCGAC
```

Note that the pairwise Hamming distance between tips (computed directly from sequences) does not necessarily equal the sum of branch distances along the tree path connecting them, because mutations along the tree may involve reversions or convergent changes. For instance, A and C are separated by tree path A→Y→X→C with distance 2 + 1 + 3 = 6, which in this case happens to match the pairwise Hamming distance, but this is not guaranteed in general.

### File naming

Each pairwise trajectory is saved as `{tip1}__{tip2}.fasta` (double underscore separator) inside a sharded tar.zst archive (e.g. `pairwise-train-000.tar.zst`).
