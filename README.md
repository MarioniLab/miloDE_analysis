# miloDE_analysis
Folder contains scripts to generate and analyse data for miloDE

Specifically:

`./generate_data/` contains scripts to generate data for the analysis (e.g. neighbourhood assignments and miloDE estimations):

```
`./generate_data/datasets/`: get SCE dta for mouse embryo (WT, chimera Tal1- and combined).
`./generate_data/dev/augur_per_nhood/`: generate data for ad hoc neighborourhood selection using AUCs.
`./generate_data/dev/DE_detection_n_cells/`: generate data for estimating how DE detection depends on number of samples and replicates.
`./generate_data/dev/embedding_analysis/`: generate data for assesing how embedding choice affects neighbourhood 'homogeneity' and DE detection.
`./generate_data/dev/order_analysis/`: generate data for assesing how order choice affects neighbourhood 'homogeneity' and DE detection.
`./generate_data/simulations/`: generate data for simulations and how miloDE performs on them.
`./generate_data/lung_macrophages/`: generate data for IPF analysis.
```


`./analysis/` contains scripts to analyse data and plot figures:

```
`./analysis/dev/augur_per_nhood`: analyse ad hoc neighborourhood selection using AUCs.
`./analysis/dev/DE_detection_n_cells`: analyse how DE detection depends on number of samples and replicates.
`./analysis/dev/embedding_analysis`: analyse how embedding choice affects neighbourhood 'homogeneity' and DE detection.
`./analysis/dev/order_analysis`: analyse how order choice affects neighbourhood 'homogeneity' and DE detection.
`./analysis/dev/nhood_refinement`: analyse how neighbourhood refinement balances number of neighbourhoods and 'coverage'.
`./analysis/dev/nhood_size_estimate`: analyse how [order-k] affects neighbourhood sizes.
`./analysis/dev/cartoons`: generate cartoons for Fig. 1.
`./analysis/simulations`: analyse simulations (Fig. 2).
`./analysis/chimera_tal1`: analyse Tal1- mouse embryo chimera data.
`./analysis/lung_macrophages`: analyse IPF, macrophages.
```
