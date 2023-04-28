# Necessary libraries
import os,sys
import numpy as np
import pandas as pd
import scanpy as sc
import pandas as pd
import scvi
import torch
import argparse

#import argparse
device = torch.device("cuda")

# set pars
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
parser.add_argument("n_layers", type=int, help="n_layers")
parser.add_argument("counts_assay", help="counts_assay")
args = parser.parse_args()


root_dir = '/nfs/research/marioni/alsu/snc/'
adata = sc.read_h5ad(root_dir + 'data/sces/thymus/sce' + '.h5ad')
adata.obs['batch'] = adata.obs['batch'].astype("category")
adata.layers["counts_all_genefull"] = adata.X

adata.X = adata.layers[args.counts_assay]
adata.layers["counts"] = adata.X


# Func
def _train_model(adata_ref, labels_key='celltype', batch_col=None , n_layers=2):
    # assign batch if provided
    if batch_col is None:
        adata_ref.obs['batch'] = '1'
    else:
        adata_ref.obs['batch'] = adata_ref.obs[batch_col].copy()

    # set up model    
    scvi.model.SCVI.setup_anndata(adata_ref, batch_key="batch", layer='counts')

    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=n_layers,
    )
    
    vae = scvi.model.SCVI(adata, n_layers = n_layers, n_latent=30, gene_likelihood="nb")
    vae.train()

    # retrain with scanvi
    lvae = scvi.model.SCANVI.from_scvi_model(vae, adata=adata, labels_key=labels_key, unlabeled_category="Unknown")
    lvae.train(max_epochs=50)

    # assign embedding
    adata.obsm["X_scANVI"] = lvae.get_latent_representation()
    return(lvae)


# run
vae = _train_model(adata, labels_key='celltype', batch_col='batch', n_layers = 2)

# save embedding
obsm = adata.obsm['X_scVI']
rownames = adata.obs_vector('cell')
df = pd.DataFrame(obsm, index=rownames)
df.to_csv(root_dir + 'data/sces/thymus/scvi/' + args.counts_assay +'__n_' + str(args.n_layers) +'.csv', index=True)



