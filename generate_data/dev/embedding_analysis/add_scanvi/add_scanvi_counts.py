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
def _train_model(adata, labels_key='celltype', batch_col=None , n_layers=2):
    # assign batch if provided
    if batch_col is None:
        adata.obs['batch'] = '1'
    else:
        adata.obs['batch'] = adata.obs[batch_col].copy()

    # set up model    
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch", layer='counts')
    # train
    vae = scvi.model.SCVI(adata, n_layers = n_layers, n_latent=30, gene_likelihood="nb")
    vae.train()

    # retrain with scanvi
    lvae = scvi.model.SCANVI.from_scvi_model(vae, adata=adata, labels_key=labels_key, unlabeled_category="Unknown")
    lvae.train(max_epochs=50)

    # assign embedding
    adata.obsm["X_scANVI"] = lvae.get_latent_representation()
    return(lvae)


# run
lvae = _train_model(adata, labels_key='celltype', batch_col='sample', n_layers = 2)

# save embedding
obsm = adata.obsm['X_scANVI']
rownames = adata.obs_vector('cell')
df = pd.DataFrame(obsm, index=rownames)
df.to_csv(root_dir + 'data/sces/mouse_embryo/chimera_WT/scvi/scanvi_unsupervised_' + args.hvg_file +'.csv', index=True)



