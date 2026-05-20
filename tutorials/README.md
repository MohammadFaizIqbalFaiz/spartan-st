# Spartan tutorial notebooks

This directory contains reviewer-oriented and tutorial-style notebooks demonstrating Spartan across multiple spatial transcriptomics technologies and analysis tasks:
1. ImagingDatasetSpartan.ipynb : Detailed notebook that reproduces results for MERFISH dataset and gives detailed demonstration of the alpha-selection strategy.
2. SequencingDatasetSpartan.ipynb : Detailed notebook that reproduces results for Stereo-seq dataset and gives detailed demonstration of the alpha-selection strategy.
3. VisiumHD Analysis Using Spartan.ipynb : Detailed notebook that reproduces Spartan's Visium HD analysis results in Figure 5 and Figure 6 of the paper.
4. SpartanSVGDiscovery.ipynb: Detailed notebook on Spartan's SVG discovery results. It reproduces all SVG figures presented in Figure 4 of the paper. 

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

| Configuration | Sample | Neighbors | Resolution [latent, spatial]  | Layer ratio |
|---|---:|---:|---|---:|
| KNN10 + PCA | 417 | 10 | [0.500, 0.600] | 1.6 |
| KNN10 + MULTISPATI PCA | 417 | 10 | [0.450, 0.600] | 1.6 |
| KNN10 + PCA | 419 | 10 | [0.250, 0.900] | 1.6 |
| KNN10 + MULTISPATI PCA | 419 | 10 | [0.300, 0.700] | 1.6 |
| KNN10 + PCA | 424 | 10 | [0.600, 0.800] | 1.6 |
| KNN10 + MULTISPATI PCA | 424 | 10 | [0.400, 1.050] | 1.6 |

### BaristaSeq

| Configuration | Slice | Neighbors | Resolution [latent, spatial]  | Layer ratio |
|---|---|---:|---|---:|
| KNN10 + PCA | slice 1 | 10 | [0.500, 1.300] | 1.8 |
| KNN10 + MULTISPATI PCA | slice 1 | 10 | [0.350, 1.500] | 1.8 |
| KNN10 + PCA | slice 2 | 10 | [0.513, 1.262] | 1.8 |
| KNN10 + MULTISPATI PCA | slice 2 | 10 | [0.350, 1.450] | 1.8 |
| KNN10 + PCA | slice 3 | 10 | [0.600, 0.800] | 1.8 |
| KNN10 + MULTISPATI PCA | slice 3 | 10 | [0.600, 1.100] | 1.8 |

### STARmap* / STARmap1k

| Configuration | Sample | Neighbors | Resolution [latent, spatial]  | Layer ratio |
|---|---|---:|---|---:|
| KNN10 + PCA | sample 1 | 10 | [0.600, 1.450] | 1.4 |
| KNN10 + MULTISPATI PCA | sample 1 | 10 | [0.500, 1.200] | 1.4 |

### osmFISH

| Configuration | Sample | Neighbors | Resolution [latent, spatial]  | Layer ratio |
|---|---|---:|---|---:|
| KNN10 + PCA | sample 1 | 10 | [0.875, 1.200] | 1.2 |
| KNN10 + MULTISPATI PCA | sample 1 | 10 | [0.875, 1.300] | 1.2 |

---

## Notes

- NichePCA, BANKSY, and SpatialLeiden were run in their respective method-specific environments.
- Resolution entry against each sample is the best resolution that give the highest NMI for that sample. Useful to reproduce results presented in Fig 4A-B and Supplementary Figure 13A-D.
- Where possible, common preprocessing, PCA construction, spatial coordinates, and evaluation labels were harmonized across methods.
- SpatialLeiden was evaluated using both PCA and MULTISPATI PCA configurations where available.
- For NichePCA, the given Visium HD parameter settings were used to generate the spatial domains, presented in Figure 5B of the paper. 
- Method-specific parameters were selected according to the benchmarking protocol described in the manuscript Methods.
- Ground-truth labels were used only for benchmarking metric calculation.
