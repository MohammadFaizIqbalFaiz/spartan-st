# Command-line interface

Spartan includes a lightweight CLI:

```bash
spartan --help
```

Available subcommands include:

- `preprocess-imaging`
- `preprocess-sequencing`
- `build-graphs`
- `domains`
- `svg`
- `alpha-select`

Example:

```bash
spartan domains \
  --h5ad input.h5ad \
  --out output_spartan.h5ad \
  --spatial-coord generic \
  --spatial-neighborhood knn \
  --spatial-neighs 10 \
  --gene-neighs 15 \
  --alpha 0.80 \
  --beta1 0.10 \
  --beta2 0.40 \
  --resolution 1.0
```
