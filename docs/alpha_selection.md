# Alpha-selection and operating-regime analysis

These workflows are intended for benchmarking-oriented analysis when ground-truth annotations are available. Spartan spatial domain detection itself remains unsupervised.

```python
summary_df = sp.tl.initiate_alpha_selection(
    adata,
    lower_alpha=0.50,
    upper_alpha=0.90,
    step_alpha=0.01,
    lower_resolution=0.50,
    upper_resolution=2.00,
    step_resolution=0.05,
    lower_nlsas=40,
    upper_nlsas=70,
    n_jobs=4,
    config="lsg",
    seed=1,
    use_nLSAS=True,
)
```

Consensus alpha:

```python
alpha_summary = sp.tl.consensus_alpha(summary_list, top_k=50)
```
