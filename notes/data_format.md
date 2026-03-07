# Data Format

This document explains the two trajectory formats produced by the pipeline: **forwards trajectories** (root-to-tip) and **pairwise trajectories** (tip-to-tip).

## Worked example

Consider a small phylogenetic tree with 3 tips (A, B, C) and 2 internal nodes (X, Y), where X is the root. Each node has a 10-nucleotide sequence. Branch labels show the Hamming distance (number of nucleotide differences, ignoring positions where either sequence has a gap `-` or ambiguous base `N`) between parent and child:

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
A:   A  T  C  G  A  G  C  G  A  T     2 changes from Y (pos 4: A→G reversion, pos 6: T→G)
B:   A  T  C  A  A  T  C  G  G  T     1 change from Y (pos 9: A→G)
C:   A  G  C  G  G  T  C  G  A  C     3 changes from X (pos 2: T→G, pos 5: A→G, pos 10: T→C)
```

Note that the Y→A branch includes a **reversion**: position 4 mutated G→A on the X→Y branch, then reverted A→G on the Y→A branch, returning to the original root nucleotide.

## Forwards trajectories

A forwards trajectory traces the evolutionary path from root to a single tip. Each entry in the FASTA file records a node along the path. The header format is `>{node_name}|{branch_hamming_distance}|{direct_hamming_distance}`, where:

- **branch distance** is the Hamming distance from the previous emitted node (0 for the root)
- **direct distance** is the Hamming distance from the start (root) node

Intermediate nodes with zero branch distance are skipped; the tip is always emitted.

### Example: trajectory for tip A

Path: X → Y → A

File `A.fasta`:
```
>X|0|0
ATCGATCGAT
>Y|1|1
ATCAATCGAT
>A|2|1
ATCGAGCGAT
```

- `X|0|0` — root, branch distance 0, direct distance 0
- `Y|1|1` — 1 mutation from X, 1 mutation from root X
- `A|2|1` — 2 mutations from Y (pos 4 reversion, pos 6), but only 1 mutation from root X (pos 6). The direct distance (1) is less than the cumulative branch distance (1 + 2 = 3) because position 4 reverted to the root nucleotide.

### Example: trajectory for tip B

Path: X → Y → B

File `B.fasta`:
```
>X|0|0
ATCGATCGAT
>Y|1|1
ATCAATCGAT
>B|1|2
ATCAATCGGT
```

- `B|1|2` — 1 mutation from Y (pos 9), but 2 mutations from root X (pos 4, 9)

### Example: trajectory for tip C

Path: X → C

File `C.fasta`:
```
>X|0|0
ATCGATCGAT
>C|3|3
AGCGGTCGAC
```

Node Y does not appear because tip C descends directly from root X. For single-hop paths, branch distance and direct distance are always equal.

### File naming

Each forwards trajectory is saved as `{tip_name}.fasta` inside a sharded tar.zst archive (e.g. `forwards-train-000.tar.zst`).

## Pairwise trajectories

A pairwise trajectory contains exactly two tip sequences. Headers use the same three-field format as forwards trajectories: `>{name}|{branch_distance}|{direct_distance}`. The first tip gets `|0|0` and the second gets `|{hamming}|{hamming}`, where the Hamming distance is computed directly between the two tip sequences (ignoring positions with gaps `-` or ambiguous bases `N`). Branch and direct distance are always identical for pairwise trajectories since there are only two frames.

### Example: pair A and B

```
Pos: 1  2  3  4  5  6  7  8  9  10
A:   A  T  C  G  A  G  C  G  A  T
B:   A  T  C  A  A  T  C  G  G  T
              *     *        *
```

Differences at positions 4, 6, 9 — Hamming distance of 3.

File `A__B.fasta`:
```
>A|0|0
ATCGAGCGAT
>B|3|3
ATCAATCGGT
```

### Example: pair A and C

```
Pos: 1  2  3  4  5  6  7  8  9  10
A:   A  T  C  G  A  G  C  G  A  T
C:   A  G  C  G  G  T  C  G  A  C
        *        *  *           *
```

Differences at positions 2, 5, 6, 10 — Hamming distance of 4.

File `A__C.fasta`:
```
>A|0|0
ATCGAGCGAT
>C|4|4
AGCGGTCGAC
```

Note that the pairwise Hamming distance between tips (computed directly from sequences) does not necessarily equal the sum of branch distances along the tree path connecting them, because mutations along the tree may involve reversions or convergent changes. For instance, A and C are separated by tree path A→Y→X→C with branch distances 2 + 1 + 3 = 6, but the pairwise Hamming distance is only 4, because the reversion at position 4 (which mutated on X→Y then reverted on Y→A) does not contribute to the pairwise distance.

### File naming

Each pairwise trajectory is saved as `{tip1}__{tip2}.fasta` (double underscore separator) inside a sharded tar.zst archive (e.g. `pairwise-train-000.tar.zst`).
