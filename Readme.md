BITFAM (custom GRN fork)
========================

BITFAM infers transcription factor (TF) activity from scRNA-seq by factorizing expression into per-cell activities and TF–gene weights while respecting a supplied TF–target prior.

How BITFAM works
----------------
- Model: Bayesian matrix factorization $X \approx Z W^T$ with a binary TF–target mask; optimized via Stan variational inference.
- Inputs: log-normalized counts matrix (genes × cells) and a binary TF–target network (genes × TFs) with identical row order.
- Outputs: `BITFAM_activities()` for TF activities; `BITFAM_weights()` for TF–gene weights consistent with the prior.

What this fork changes
----------------------
This is a modified version of BITFAM (https://github.com/jaleesr/BITFAM).
Changes relative to upstream:

- Package renamed from `BITFAM` to `BITFAM2` (version 1.2.0).
- `BITFAM()` takes a user-supplied binary gene × TF prior network matrix
  (`network`) directly, replacing the `species`, `interseted_TF` and
  `scATAC_obj` arguments and the ChIP-seq target sets bundled in
  `inst/extdata`.
- Removed internal variable-gene subsetting; gene selection is now the
  caller's responsibility.
- Removed automatic upper-casing of gene symbols.
- Random seed exposed as the `seed` argument (previously hard-coded to 100).
- `ncores` given a default value; roxygen documentation rewritten.
- Added validation that `data` and `network` row names agree.
- `BITFAM_scATAC()` added to the package exports.

Installation
------------
```r
# Prereqs: a working C++ toolchain for rstan (see https://mc-stan.org/rstan/)
install.packages(c("rstan", "Seurat", "devtools"))

# Install this fork
devtools::install_github("brandonlukas/BITFAM")

# Load
library(BITFAM2)  # functions are named BITFAM_*
```

Quickstart
----------
```r
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(BITFAM2)

# 1) Prepare expression (log-normalized counts)
cells <- NormalizeData(CreateSeuratObject(counts = raw_counts))
genes <- VariableFeatures(cells)
data <- GetAssayData(cells)[genes, , drop = FALSE]

# 2) Prepare prior network (long -> wide; binary)
# long_network has columns: source (TF), target (gene)
network <- long_network %>%
  filter(target %in% genes) %>%
  add_count(source) %>%
  filter(n >= 10) %>%               # drop TFs with too few targets
  bind_rows(tibble(source = "ENSURE_ALL_GENES", target = genes)) %>%
  mutate(value = 1) %>%
  pivot_wider(id_cols = target, names_from = source, values_fill = 0) %>%
  select(-ENSURE_ALL_GENES) %>%
  column_to_rownames("target")

# 3) Align gene order in data and network
data <- data[rownames(network), , drop = FALSE]

# 4) Run BITFAM and extract activities
fit <- BITFAM(data = data, network = as.matrix(network), ncores = 4)
Z <- BITFAM_activities(fit)   # cells × TFs activities
W <- BITFAM_weights(fit)      # genes × TFs weights
```

Notes
-----
- The model densifies inputs; ensure the dense matrix fits in memory.
- If convergence is slow, increase `iter` and/or decrease `tol_rel_obj` in `BITFAM()`.

License and attribution
-----------------------
Licensed under the Apache License, Version 2.0; see `LICENSE`. Modifications
made in this fork are stated above and in `NOTICE`.

Original method:

> Gao, S. et al. Identifying noise-tolerant, statistically significant
> transcription factor activity from single-cell RNA-seq data.
> *Genome Research* **31**, 1296–1311 (2021).
