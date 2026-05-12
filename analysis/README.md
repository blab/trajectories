# Analysis

## Comparing UShER and Viridian trees for spike-lg dataset

The `compare_trees.py` script validates tree quality by comparing reversion-to-reference rates between UShER and Viridian mutation-annotated trees. Both trees are subsampled to 10,000 tips and analyzed over the spike S1 region (positions 21563–23617).

### Trees compared

- **UShER**: https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2026/02/21/public-2026-02-21.all.masked.pb.gz (~8M sequences)
- **Viridian**: https://ndownloader.figshare.com/files/49691037 (~4.5M sequences, from [Figshare article 27194547](https://figshare.com/articles/dataset/Global_Viridian_SARS-CoV-2_phylogenetic_tree/27194547))

The Viridian tree is built from consensus sequences assembled with an amplicon-aware pipeline ([Hunt et al., Nature Methods 2025](https://www.nature.com/articles/s41592-025-02947-1)), which corrects systematic errors from tiled amplicon sequencing.

### Results

| Metric | UShER | Viridian |
|---|---|---|
| Branches with mutations | 558,627 | 242,307 |
| Total mutations | 621,865 | 263,346 |
| Mean mutations/branch | 1.11 | 1.09 |
| **Total reversions** | **69,246** | **5,379** |
| **Reversion rate (% of mutations)** | **11.14%** | **2.04%** |
| Branches with reversions | 60,040 | 5,040 |
| % branches with reversions | 10.75% | 2.08% |

The Viridian tree has **5.5x fewer reversions** to reference in the spike S1 region. Reversions are mutations where the derived base matches the reference, which typically indicate systematic consensus errors rather than true evolutionary events. The lower reversion rate confirms the Viridian tree produces cleaner training data for the trajectories pipeline.

Indel rates could not be compared because the UShER MAT protobuf format encodes only single-nucleotide substitutions, not indels.

### Running the comparison

```bash
# Download trees
curl -L -o /tmp/usher.pb.gz "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2026/02/21/public-2026-02-21.all.masked.pb.gz"
curl -L -o /tmp/viridian.pb.gz "https://ndownloader.figshare.com/files/49691037"

# Run comparison (requires matUtils)
python analysis/compare_trees.py \
  --usher /tmp/usher.pb.gz \
  --viridian /tmp/viridian.pb.gz \
  --subsample 10000 \
  --trim-begin 21563 --trim-end 23617
```
