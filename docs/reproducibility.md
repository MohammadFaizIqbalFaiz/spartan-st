# Reproducibility

## General user environment

```bash
conda env create -f envs/environment.core.yml
conda activate spartan-core
pip install -e .
```

## Paper reproduction environment

```bash
conda env create -f envs/environment.paper.lock.yml
conda activate spartan-paperS
pip install -e .
```

The documented release is `v0.2.0`.
