# Necessary libraries
import os,sys
import numpy as np
import pandas as pd
import scanpy as sc
import pandas as pd
import scvi
import torch
#import argparse
device = torch.device("cuda")

root_dir = '/nfs/research/marioni/alsu/hubmap_metaRef/data/'
adata = sc.read_h5ad(root_dir + 'sces/mouse_embryo/chimera/mouse_embryo_chimera_tal' + '.h5ad')
adata.obs['sample'] = adata.obs['sample'].astype("category")

# Func
def _train_model(adata_ref, batch_col=None):
    if batch_col is None:
        adata_ref.obs['batch'] = '1'
    else:
        adata_ref.obs['batch'] = adata_ref.obs[batch_col].copy()
    ## Select genes
    adata_ref.layers['counts'] = adata_ref.X
    # adata_ref = adata_ref[:,hvgs].copy()

    scvi.model.SCVI.setup_anndata(adata_ref, batch_key="batch", layer='counts')

    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
    )

    vae_ref = scvi.model.SCVI(
        adata_ref,
        **arches_params
    )

    vae_ref.train()
    adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()
    return(vae_ref)

def _fit_model(adata_query, vae_ref, batch_col=None):
    if batch_col is None:
        adata_query.obs['batch'] = '1'
    else:
        adata_query.obs['batch'] = adata_query.obs[batch_col].copy()
    adata_query_fit = adata_query[:, vae_ref.adata.var_names].copy()
    adata_query_fit.layers['counts'] = adata_query_fit.X
    vae_q = scvi.model.SCVI.load_query_data(
        adata_query_fit,
        vae_ref
        )

    vae_q.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))
    adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()
    return(vae_q)


# split data
adata_ref = adata[adata.obs['type']=='wt'].copy()
adata_query = adata[adata.obs['type']=='chimera'].copy()

# run
vae_ref = _train_model(adata_ref, batch_col='sample')
vae_q_fit = _fit_model(adata_query, vae_ref, batch_col='sample')

# save
adata_full = adata_query.concatenate(adata_ref, batch_key='dataset')
adata_full.layers['counts'] = adata_full.X.copy()
adata_full.obs['_scvi_labels'] = 0
#adata_full.obs['batch'] = adata_full.obs['batch'].astype("str")
adata_full.obsm["X_scVI"] = vae_q_fit.get_latent_representation(adata_full)
adata_full.write_h5ad(root_dir + 'sces/mouse_embryo/chimera/mouse_embryo_chimera_tal_w_scarches' + '.h5ad')

