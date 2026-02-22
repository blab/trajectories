# Mutation accuracy metric

## Definition

Given a source sequence and a target sequence, the **truth set** is the list of exact mutations between them (e.g., A3T, G7C), with **N = |truth set|**. A model's **predicted set** is its list of mutations, with **M = |predicted set|**. A mutation is **correct** if it appears in both sets — same position, same nucleotide.

$$\text{mutation\_accuracy} = \frac{\text{correct} - |M - N|}{N}$$

where correct = |truth ∩ predicted|.

**Interpretation**: When the model predicts exactly N mutations, this reduces to the proportion of correct mutations. When it predicts too few or too many, the |M − N| term penalizes miscalibration of the mutation count.

**Properties**:
- Perfect prediction → 1.0
- Copy source (no mutations) → −1.0
- Random N mutations (N ≪ L) → ≈ 0
- Range: [−∞, 1.0], though in practice bounded below by ≈ −1 for reasonable models

## Worked example

Source sequence (L = 10):

```
Position:  1  2  3  4  5  6  7  8  9  10
Source:    A  C  G  T  A  C  G  T  A  C
Target:    A  C  T  T  A  C  A  T  G  C
```

**Truth set** (N = 3): {G3T, G7A, A9G}

---

### 1. Perfect prediction

Predicted: {G3T, G7A, A9G} → M = 3, correct = 3

$$\frac{3 - |3 - 3|}{3} = \frac{3}{3} = \textbf{1.0}$$

### 2. Copy source (no mutations predicted)

Predicted: {} → M = 0, correct = 0

$$\frac{0 - |0 - 3|}{3} = \frac{-3}{3} = \textbf{-1.0}$$

### 3. Correct count, all at wrong sites

Predicted: {A1T, C2G, T4A} → M = 3, correct = 0

$$\frac{0 - |3 - 3|}{3} = \frac{0}{3} = \textbf{0.0}$$

### 4. Correct count, right sites but wrong nucleotides

Predicted: {G3A, G7T, A9C} → M = 3, correct = 0

$$\frac{0 - |3 - 3|}{3} = \frac{0}{3} = \textbf{0.0}$$

Cases 3 and 4 score identically: a wrong mutation is a wrong mutation regardless of whether it found a real mutation site.

### 5. Two correct, one wrong

Predicted: {G3T, G7A, T4A} → M = 3, correct = 2

$$\frac{2 - |3 - 3|}{3} = \frac{2}{3} = \textbf{0.67}$$

### 6. Two correct, stops early

Predicted: {G3T, G7A} → M = 2, correct = 2

$$\frac{2 - |2 - 3|}{3} = \frac{1}{3} = \textbf{0.33}$$

The model got everything it attempted right, but is penalized for under-predicting the mutation count.

### 7. Two correct, pads with two extra wrong mutations

Predicted: {G3T, G7A, A1T, C6G} → M = 4, correct = 2

$$\frac{2 - |4 - 3|}{3} = \frac{1}{3} = \textbf{0.33}$$

Same score as case 6: one extra mutation costs the same as one missing mutation. Both represent a calibration error of 1.

### 8. Predicts 6 mutations, none correct

Predicted: {A1T, C2G, T4A, A5G, C6T, T8A} → M = 6, correct = 0

$$\frac{0 - |6 - 3|}{3} = \frac{-3}{3} = \textbf{-1.0}$$

Same score as copy-source: the model made twice the expected number of mutations and got none right.

## Summary table

| Scenario | correct | M | excess/deficit | Score |
|---|---|---|---|---|
| Perfect | 3 | 3 | 0 | **1.0** |
| Copy source | 0 | 0 | 3 | **−1.0** |
| Wrong sites, right count | 0 | 3 | 0 | **0.0** |
| Right sites, wrong nucleotides | 0 | 3 | 0 | **0.0** |
| 2 correct, 1 wrong | 2 | 3 | 0 | **0.67** |
| 2 correct, stops early | 2 | 2 | 1 | **0.33** |
| 2 correct, 1 extra wrong | 2 | 4 | 1 | **0.33** |
| 6 wrong mutations | 0 | 6 | 3 | **−1.0** |
