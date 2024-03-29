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
parser.add_argument("hvg_file", help="hvg_file")
args = parser.parse_args()

# load adata
root_dir = '/nfs/research/marioni/alsu/hubmap_metaRef/'
adata = sc.read_h5ad(root_dir + 'data/sces/mouse_embryo/chimera_WT/mouse_embryo_wt_chimera' + '.h5ad')
adata.obs['sample'] = adata.obs['sample'].astype("category")
adata.layers["counts"] = adata.X

# get genes
df = pd.read_csv(root_dir + 'data/sces/mouse_embryo/chimera_WT/HVGs/'+ args.hvg_file + '.csv')
genes = df['HVG'].tolist()
adata = adata[:, adata.var_names.isin(genes)]


# Func
def _train_model(adata_ref, batch_col=None , n_layers=2):
    if batch_col is None:
        adata_ref.obs['batch'] = '1'
    else:
        adata_ref.obs['batch'] = adata_ref.obs[batch_col].copy()

    scvi.model.SCVI.setup_anndata(adata_ref, batch_key="batch", layer='counts')

    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=n_layers,
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


# split data: first sn to sce
adata_ref = adata[adata.obs['type']=='wt'].copy()
adata_query = adata[adata.obs['type']!='wt'].copy()

# run
vae_ref = _train_model(adata_ref, batch_col='sample', n_layers = 2)
vae_q_fit = _fit_model(adata_query, vae_ref, batch_col='sample')

# get_observation
adata_full = adata_query.concatenate(adata_ref, batch_key='dataset')
adata_full.layers['counts'] = adata_full.X.copy()
adata_full.obs['_scvi_labels'] = 0
#adata_full.obs['batch'] = adata_full.obs['batch'].astype("str")
adata_full.obsm["X_scVI"] = vae_q_fit.get_latent_representation(adata_full)

# save embedding
obsm = adata_full.obsm['X_scVI']
rownames = adata_full.obs_vector('cell')
df = pd.DataFrame(obsm, index=rownames)
df.to_csv(root_dir + 'data/sces/mouse_embryo/chimera_WT/scvi/scvi_supervised_' + args.hvg_file +'.csv', index=True)




