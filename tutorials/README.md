# Spartan tutorial notebooks

This directory contains reviewer-oriented and tutorial-style notebooks demonstrating Spartan across multiple spatial transcriptomics technologies and analysis tasks:

1. `ImagingDatasetSpartan.ipynb`: Detailed notebook that reproduces results for the MERFISH dataset and gives a detailed demonstration of the alpha-selection strategy.
2. `SequencingDatasetSpartan.ipynb`: Detailed notebook that reproduces results for the Stereo-seq dataset and gives a detailed demonstration of the alpha-selection strategy.
3. `VisiumHD Analysis Using Spartan.ipynb`: Detailed notebook that reproduces Spartan's Visium HD analysis results in Figure 5 and Figure 6 of the paper.
4. `SpartanSVGDiscovery.ipynb`: Detailed notebook on Spartan's SVG discovery results. It reproduces all SVG figures presented in Figure 4 of the paper.
5. `Generate_Benchmarking_and_Ablation_Study_Figures tutorial.ipynb`: Detailed notebook to generate ablation and benchmarking study panel figures in Figure 3 of the paper.

---

## Spartan parameter settings used for SDMBench and recent-method benchmarking

This section documents the Spartan parameter settings used for the spatial-domain benchmarking analyses. `Alpha` and `Resolution` refer to the Delaunay configuration, whereas `Alpha (KNN)` and `Resolution (KNN)` refer to the KNN configuration. The `Resolution` values correspond to the configurations at which the maximum metric values were observed. Median rows summarize the median metric values across retained configurations. Unless otherwise noted, `β1=0.26` and `β2=0.24`.

The primary manuscript benchmarking metrics are `ARI`, `NMI`, `HOM`, `COM`, and the composite `R` score. `CHAOS` is included here as an additional GitHub/reproducibility diagnostic for spatial coherence and domain compactness. It was computed during practical reproducibility analyses to provide extra transparency on spatial stability, but it was not used as a primary manuscript ranking metric.

Metric strings are reported as:

```text
ARI/NMI/HOM/COM/CHAOS
```

### MERFISH (SDMBench)

| Bregma / sample | α (Delaunay) | α (KNN) | Resolution (Delaunay) | Resolution (KNN) | Resolution grid | β1,β2 | Best metrics (Delaunay) | Best metrics (KNN) | Median metrics (Delaunay) | Median metrics (KNN) |
|---|---:|---:|---:|---:|---|---|---|---|---|---|
| 0.04 | 0.75 | 0.69 | 0.77 | 0.62 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.53/0.65/0.68/0.63/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.55/0.66/0.68/0.63/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.53/0.65/0.68/0.63/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.54/0.66/0.68/0.63/0.03 |
| 0.09 | 0.75 | 0.69 | 0.79 | 0.71 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.54/0.66/0.67/0.64/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.53/0.66/0.68/0.64/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.44/0.59/0.61/0.57/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.53/0.66/0.68/0.64/0.03 |
| 0.14 | 0.75 | 0.69 | 0.77 | 0.78 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.52/0.63/0.65/0.6/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.52/0.62/0.65/0.6/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.51/0.62/0.64/0.6/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.51/0.62/0.64/0.59/0.03 |
| 0.19 | 0.75 | 0.69 | 0.72 | 0.77 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.69/0.74/0.73/0.74/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.66/0.72/0.72/0.73/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.65/0.71/0.7/0.71/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.64/0.71/0.71/0.71/0.03 |
| 0.24 | 0.75 | 0.69 | 0.75 | 0.73 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.68/0.7/0.7/0.7/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.69/0.73/0.72/0.73/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.68/0.7/0.7/0.7/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.69/0.73/0.72/0.73/0.03 |

### osmFISH

| Sample | α (Delaunay) | α (KNN) | Resolution (Delaunay) | Resolution (KNN) | Resolution grid | β1,β2 | Best metrics (Delaunay) | Best metrics (KNN) | Median metrics (Delaunay) | Median metrics (KNN) |
|---|---:|---:|---:|---:|---|---|---|---|---|---|
| Sample 1 | 0.676 | 0.728 | 1.07 | 0.71 | [0.5,1.5] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.57/0.65/0.67/0.64/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.65/0.73/0.74/0.73/0.02 | ARI/NMI/HOM/COM/CHAOS = 0.57/0.65/0.67/0.64/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.64/0.72/0.73/0.71/0.02 |

### STARmap1k

| Sample | α (Delaunay) | α (KNN) | Resolution (Delaunay) | Resolution (KNN) | Resolution grid | β1,β2 | Best metrics (Delaunay) | Best metrics (KNN) | Median metrics (Delaunay) | Median metrics (KNN) |
|---|---:|---:|---:|---:|---|---|---|---|---|---|
| Sample 1 | 0.712 | 0.661 | 0.99 | 0.84 | [0.5,1.5] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.6/0.71/0.7/0.71/0.06 | ARI/NMI/HOM/COM/CHAOS = 0.57/0.68/0.67/0.69/0.06 | ARI/NMI/HOM/COM/CHAOS = 0.6/0.71/0.7/0.71/0.06 | ARI/NMI/HOM/COM/CHAOS = 0.57/0.68/0.67/0.69/0.06 |

### STARmap

| Sample | α (Delaunay) | α (KNN) | Resolution (Delaunay) | Resolution (KNN) | Resolution grid | β1,β2 | Best metrics (Delaunay) | Best metrics (KNN) | Median metrics (Delaunay) | Median metrics (KNN) |
|---|---:|---:|---:|---:|---|---|---|---|---|---|
| Sample 417 | 0.799 | 0.728 | 0.52 | 0.5 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.84/0.81/0.82/0.8/0.07 | ARI/NMI/HOM/COM/CHAOS = 0.84/0.81/0.82/0.8/0.07 | ARI/NMI/HOM/COM/CHAOS = 0.82/0.79/0.81/0.78/0.07 | ARI/NMI/HOM/COM/CHAOS = 0.8/0.79/0.8/0.78/0.07 |
| Sample 419 | 0.799 | 0.728 | 0.5 | 0.5 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.61/0.71/0.74/0.67/0.07 | ARI/NMI/HOM/COM/CHAOS = 0.62/0.69/0.72/0.65/0.07 | ARI/NMI/HOM/COM/CHAOS = 0.54/0.62/0.66/0.58/0.07 | ARI/NMI/HOM/COM/CHAOS = 0.62/0.69/0.72/0.65/0.07 |
| Sample 424 | 0.799 | 0.728 | 0.52 | 0.54 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.76/0.71/0.74/0.68/0.07 | ARI/NMI/HOM/COM/CHAOS = 0.8/0.73/0.75/0.71/0.07 | ARI/NMI/HOM/COM/CHAOS = 0.76/0.71/0.74/0.68/0.07 | ARI/NMI/HOM/COM/CHAOS = 0.8/0.73/0.75/0.71/0.07 |

### Stereo-seq

Stereo-seq was evaluated using the KNN configuration in this table.

| Sample | α (KNN) | Resolution (KNN) | Resolution grid | β1,β2 | Best metrics (KNN) | Median metrics (KNN) |
|---|---:|---:|---|---|---|---|
| E9.5E1S1 | 0.56 | 0.62 | [0.5,1.2] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.34/0.61/0.64/0.58/0.04 | ARI/NMI/HOM/COM/CHAOS = 0.34/0.61/0.64/0.58/0.04 |
| E9.5E2S1 | 0.56 | 1.08 | [0.5,1.2] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.47/0.65/0.66/0.64/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.46/0.64/0.65/0.63/0.05 |
| E9.5E2S2 | 0.56 | 1.19 | [0.5,1.2] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.57/0.7/0.7/0.7/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.56/0.7/0.7/0.7/0.05 |
| E9.5E2S3 | 0.56 | 0.9 | [0.5,1.2] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.58/0.69/0.7/0.68/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.58/0.68/0.69/0.67/0.05 |
| E9.5E2S4 | 0.56 | 1.09 | [0.5,1.2] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.5/0.62/0.65/0.6/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.49/0.62/0.64/0.6/0.05 |
| E10.5E1S1 | 0.56 | 0.8 | [0.5,1.2] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.48/0.62/0.66/0.59/0.02 | ARI/NMI/HOM/COM/CHAOS = 0.48/0.62/0.66/0.59/0.02 |
| E10.5E1S2 | 0.56 | 0.64 | [0.5,1.2] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.26/0.46/0.47/0.46/0.02 | ARI/NMI/HOM/COM/CHAOS = 0.26/0.46/0.47/0.46/0.02 |
| E10.5E1S3 | 0.56 | 0.8 | [0.5,1.2] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.39/0.61/0.63/0.6/0.02 | ARI/NMI/HOM/COM/CHAOS = 0.38/0.61/0.62/0.6/0.02 |
| E10.5E2S1 | 0.56 | 1.12 | [0.5,1.2] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.44/0.65/0.66/0.64/0.03 | ARI/NMI/HOM/COM/CHAOS = 0.44/0.65/0.66/0.64/0.03 |

### BaristaSeq

| Slice | α (Delaunay) | α (KNN) | Resolution (Delaunay) | Resolution (KNN) | Resolution grid | β1,β2 | Best metrics (Delaunay) | Best metrics (KNN) | Median metrics (Delaunay) | Median metrics (KNN) |
|---|---:|---:|---:|---:|---|---|---|---|---|---|
| Slice 1 | 0.793 | 0.723 | 0.79 | 0.79 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.66/0.72/0.7/0.74/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.64/0.72/0.72/0.73/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.66/0.72/0.7/0.74/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.64/0.72/0.72/0.73/0.05 |
| Slice 2 | 0.793 | 0.723 | 0.75 | 0.68 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.79/0.8/0.82/0.79/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.78/0.79/0.8/0.79/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.79/0.8/0.82/0.79/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.78/0.79/0.8/0.79/0.05 |
| Slice 3 | 0.793 | 0.723 | 0.66 | 0.64 | [0.5,0.9] | [0.26,0.24] | ARI/NMI/HOM/COM/CHAOS = 0.62/0.7/0.71/0.69/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.59/0.66/0.67/0.65/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.62/0.7/0.71/0.69/0.05 | ARI/NMI/HOM/COM/CHAOS = 0.59/0.66/0.67/0.65/0.05 |

### Vizgen MERFISH

For Vizgen MERFISH, the Delaunay and KNN configurations used different resolution ranges. The first range shown is for Delaunay and the second range is for KNN.

| Sample | α (Delaunay) | α (KNN) | Resolution (Delaunay) | Resolution (KNN) | Resolution grid (Delaunay) | Resolution grid (KNN) | β1,β2 | Best metrics (Delaunay) | Best metrics (KNN) |
|---|---:|---:|---:|---:|---|---|---|---|---|
| S1R1 | 0.85 | 0.82 | 2 | 1.7 | [2.0,4.0] | [1.0,2.0] | [0.20,0.30] | ARI/NMI/HOM/COM/CHAOS = 0.31/0.66/0.62/0.71/0.007 | ARI/NMI/HOM/COM/CHAOS = 0.31/0.66/0.61/0.71/0.007 |
| S1R2 | 0.85 | 0.82 | 2.2 | 1.3 | [2.0,4.0] | [1.0,2.0] | [0.20,0.30] | ARI/NMI/HOM/COM/CHAOS = 0.28/0.65/0.62/0.68/0.007 | ARI/NMI/HOM/COM/CHAOS = 0.31/0.65/0.58/0.74/0.007 |
| S2R1 | 0.85 | 0.82 | 2.2 | 1.9 | [2.0,4.0] | [1.0,2.0] | [0.20,0.30] | ARI/NMI/HOM/COM/CHAOS = 0.36/0.68/0.63/0.74/0.007 | ARI/NMI/HOM/COM/CHAOS = 0.36/0.68/0.63/0.74/0.007 |
| S2R3 | 0.85 | 0.82 | 2.1 | 1.9 | [2.0,4.0] | [1.0,2.0] | [0.20,0.30] | ARI/NMI/HOM/COM/CHAOS = 0.38/0.71/0.64/0.78/0.007 | ARI/NMI/HOM/COM/CHAOS = 0.38/0.7/0.64/0.78/0.007 |


---

## Notes on CHAOS

- `CHAOS` is reported only as an additional spatial-stability diagnostic in the GitHub/reproducibility tables.
- Lower `CHAOS` values indicate stronger spatial compactness/coherence of predicted domains.
- `CHAOS` was not used as a primary manuscript ranking metric and should be interpreted as complementary to `ARI`, `NMI`, `HOM`, `COM`, and the composite `R` score.
- Including `CHAOS` in these tables provides additional transparency because spatially compact domains are important for spatial transcriptomics, especially when label-matching metrics alone do not fully capture spatial fragmentation.

---

## Competing method parameter settings

This section documents the key parameter settings used for the recent competing methods included in the benchmarking analyses. Each method was run in its respective recommended software environment because several tools require mutually incompatible dependencies. Common preprocessing and harmonized input preparation are described in the manuscript Methods; the tables below record the method-specific parameters used for each spatial transcriptomics dataset.

The reported settings are provided for transparency and reviewer reproducibility. Ground-truth labels were used only for metric calculation and benchmarking, not during unsupervised clustering by Spartan or the competing methods unless required by the original method workflow.

---

## NichePCA settings

| Dataset | Sample / slice | Neighbors | Resolution |
|---|---|---:|---:|
| MERFISH | sample 1 | KNN = 11 | 0.14 |
| MERFISH | sample 2 | KNN = 11 | 0.14 |
| MERFISH | sample 3 | KNN = 11 | 0.14 |
| MERFISH | sample 4 | KNN = 11 | 0.14 |
| MERFISH | sample 5 | KNN = 11 | 0.14 |
| STARmap | sample 417 | KNN = 4 | 0.20 |
| STARmap | sample 419 | KNN = 4 | 0.20 |
| STARmap | sample 424 | KNN = 4 | 0.20 |
| STARmap* / STARmap1k | sample 1 | KNN = 5 | 0.20 |
| osmFISH | sample 1 | KNN = 11 | 0.13 |
| BaristaSeq | slice 1 | KNN = 7 | 0.30 |
| BaristaSeq | slice 2 | KNN = 7 | 0.30 |
| BaristaSeq | slice 3 | KNN = 7 | 0.30 |
| Vizgen MERFISH | sample 1 | KNN = 28 | 0.90 |
| Vizgen MERFISH | sample 2 | KNN = 28 | 1.20 |
| Vizgen MERFISH | sample 3 | KNN = 28 | 0.90 |
| Vizgen MERFISH | sample 4 | KNN = 28 | 1.20 |
| Visium HD | sample 1 | KNN = 12 | 0.13 |

---

## BANKSY settings

| Dataset | Sample / slice | Neighbors | Resolution |
|---|---|---:|---:|
| STARmap | sample 417 | 10 | 0.50 |
| STARmap | sample 419 | 10 | 0.50 |
| STARmap | sample 424 | 10 | 0.50 |
| STARmap* / STARmap1k | sample 1 | 10 | 0.90 |
| MERFISH | sample 1 | 19 | 0.35 |
| MERFISH | sample 2 | 19 | 0.35 |
| MERFISH | sample 3 | 19 | 0.35 |
| MERFISH | sample 4 | 19 | 0.40 |
| MERFISH | sample 5 | 19 | 0.42 |
| osmFISH | sample 1 | 19 | 0.90 |
| BaristaSeq | slice 1 | 20 | 0.80 |
| BaristaSeq | slice 2 | 20 | 0.70 |
| BaristaSeq | slice 3 | 20 | 0.90 |
| Vizgen MERFISH | sample 1 | 28 | 3.00 |
| Vizgen MERFISH | sample 2 | 28 | 3.00 |
| Vizgen MERFISH | sample 3 | 28 | 3.00 |
| Vizgen MERFISH | sample 4 | 28 | 3.00 |

---

## SpatialLeiden settings

SpatialLeiden was evaluated using both PCA-based and MULTISPATI PCA-based configurations where applicable. The table below reports the neighborhood setting, resolution search interval, and layer-ratio setting used for each dataset and sample.

### MERFISH

| Configuration | Sample / Bregma | Neighbors | Resolutions [latent, spatial] | Layer ratio |
|---|---:|---:|---|---:|
| KNN10 + PCA | 0.04 | 10 | [0.163, 1.000] | 1.8 |
| KNN10 + MULTISPATI PCA | 0.04 | 10 | [0.175, 0.950] | 1.8 |
| KNN10 + PCA | 0.09 | 10 | [0.100, 1.500] | 1.8 |
| KNN10 + MULTISPATI PCA | 0.09 | 10 | [0.213, 1.375] | 1.8 |
| KNN10 + PCA | 0.14 | 10 | [0.100, 1.200] | 1.8 |
| KNN10 + MULTISPATI PCA | 0.14 | 10 | [0.150, 1.100] | 1.8 |
| KNN10 + PCA | 0.19 | 10 | [0.200, 1.200] | 1.8 |
| KNN10 + MULTISPATI PCA | 0.19 | 10 | [0.200, 1.350] | 1.8 |
| KNN10 + PCA | 0.24 | 10 | [0.225, 1.400] | 1.8 |
| KNN10 + MULTISPATI PCA | 0.24 | 10 | [0.275, 1.250] | 1.8 |

### STARmap

| Configuration | Sample | Neighbors | Resolution [latent, spatial] | Layer ratio |
|---|---:|---:|---|---:|
| KNN10 + PCA | 417 | 10 | [0.500, 0.600] | 1.6 |
| KNN10 + MULTISPATI PCA | 417 | 10 | [0.450, 0.600] | 1.6 |
| KNN10 + PCA | 419 | 10 | [0.250, 0.900] | 1.6 |
| KNN10 + MULTISPATI PCA | 419 | 10 | [0.300, 0.700] | 1.6 |
| KNN10 + PCA | 424 | 10 | [0.600, 0.800] | 1.6 |
| KNN10 + MULTISPATI PCA | 424 | 10 | [0.400, 1.050] | 1.6 |

### BaristaSeq

| Configuration | Slice | Neighbors | Resolution [latent, spatial] | Layer ratio |
|---|---|---:|---|---:|
| KNN10 + PCA | slice 1 | 10 | [0.500, 1.300] | 1.8 |
| KNN10 + MULTISPATI PCA | slice 1 | 10 | [0.350, 1.500] | 1.8 |
| KNN10 + PCA | slice 2 | 10 | [0.513, 1.262] | 1.8 |
| KNN10 + MULTISPATI PCA | slice 2 | 10 | [0.350, 1.450] | 1.8 |
| KNN10 + PCA | slice 3 | 10 | [0.600, 0.800] | 1.8 |
| KNN10 + MULTISPATI PCA | slice 3 | 10 | [0.600, 1.100] | 1.8 |

### STARmap* / STARmap1k

| Configuration | Sample | Neighbors | Resolution [latent, spatial] | Layer ratio |
|---|---|---:|---|---:|
| KNN10 + PCA | sample 1 | 10 | [0.600, 1.450] | 1.4 |
| KNN10 + MULTISPATI PCA | sample 1 | 10 | [0.500, 1.200] | 1.4 |

### osmFISH

| Configuration | Sample | Neighbors | Resolution [latent, spatial] | Layer ratio |
|---|---|---:|---|---:|
| KNN10 + PCA | sample 1 | 10 | [0.875, 1.200] | 1.2 |
| KNN10 + MULTISPATI PCA | sample 1 | 10 | [0.875, 1.300] | 1.2 |

---

## General notes

- NichePCA, BANKSY, and SpatialLeiden were run in their respective method-specific environments.
- Resolution entries for each sample indicate the resolution that gave the highest ARI for that sample. These settings support reproduction of the benchmarking results presented in Figure 3 and the corresponding supplementary figures.
- Where possible, common preprocessing, PCA construction, spatial coordinates, and evaluation labels were harmonized across methods.
- SpatialLeiden was evaluated using both PCA and MULTISPATI PCA configurations where available.
- For NichePCA, the Visium HD parameter settings were used to generate the spatial domains presented in Figure 5B of the paper.
- Method-specific parameters were selected according to the benchmarking protocol described in the manuscript Methods.
- Ground-truth labels were used only for benchmarking metric calculation.
