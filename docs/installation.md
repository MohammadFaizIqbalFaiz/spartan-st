# Installation

## General package environment

```bash
conda env create -f envs/environment.core.yml
conda activate spartan-core
pip install -e .
```

## Paper reproducibility environment

```bash
conda env create -f envs/environment.paper.lock.yml
conda activate spartan-paperS
pip install -e .
```

The paper environment contains the pinned package versions used for the manuscript analyses.
