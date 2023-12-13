import tables
import anndata
import pacmap
import cello
import os
import random
import sys
import time
import traceback
import warnings
import numpy as np
import pandas as pd
import networkx as nx
import scanpy as sc
import scrublet as scr
import decoupler as dc
import scipy.sparse as sp
from typing import Dict, Optional

import matplotlib
import matplotlib.pyplot as plt
PLOTS_DIR = os.path.join(snakemake.output["plots"])

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['figure.max_open_warning'] = 300

reseed = 42
random.seed(reseed)
np.random.seed(reseed)
n_cores=snakemake.threads

sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.figdir = PLOTS_DIR
sc.settings.set_figure_params(
    dpi=120,
    dpi_save=600,
    vector_friendly=True,
    format="pdf",
    transparent=True)
sc.settings.autoshow = False
sc.settings.autosave = True
sc.logging.print_versions()


def anndata_from_h5(file: str,
                    analyzed_barcodes_only: bool = True) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count matrix.
            True to load a limited set of barcodes: only those analyzed by the
            algorithm. This allows relevant latent variables to be loaded
            properly into adata.obs and adata.obsm, rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.

    """

    d = dict_from_h5(file)
    X = sp.csc_matrix((d.pop('data'), d.pop('indices'), d.pop('indptr')),
                      shape=d.pop('shape')).transpose().tocsr()

    # check and see if we have barcode index annotations, and if the file is filtered
    barcode_key = [k for k in d.keys() if (('barcode' in k) and ('ind' in k))]
    if len(barcode_key) > 0:
        max_barcode_ind = d[barcode_key[0]].max()
        filtered_file = (max_barcode_ind >= X.shape[0])
    else:
        filtered_file = True

    if analyzed_barcodes_only:
        if filtered_file:
            # filtered file being read, so we don't need to subset
            print('Assuming we are loading a "filtered" file that contains only cells.')
            pass
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        elif 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        else:
            print('Warning: analyzed_barcodes_only=True, but the key '
                  '"barcodes_analyzed_inds" or "barcode_indices_for_latents" '
                  'is missing from the h5 file. '
                  'Will output all barcodes, and proceed as if '
                  'analyzed_barcodes_only=False')

    # Construct the anndata object.
    adata = anndata.AnnData(X=X,
                            obs={'barcode': d.pop('barcodes').astype(str)},
                            var={'gene_name': (d.pop('gene_names') if 'gene_names' in d.keys()
                                               else d.pop('name')).astype(str)},
                            dtype=X.dtype)
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_name', inplace=True)

    # For CellRanger v2 legacy format, "gene_ids" was called "genes"... rename this
    if 'genes' in d.keys():
        d['id'] = d.pop('genes')

    # For purely aesthetic purposes, rename "id" to "gene_id"
    if 'id' in d.keys():
        d['gene_id'] = d.pop('id')

    # If genomes are empty, try to guess them based on gene_id
    if 'genome' in d.keys():
        if np.array([s.decode() == '' for s in d['genome']]).all():
            if '_' in d['gene_id'][0].decode():
                print('Genome field blank, so attempting to guess genomes based on gene_id prefixes')
                d['genome'] = np.array([s.decode().split('_')[0] for s in d['gene_id']], dtype=str)

    # Add other information to the anndata object in the appropriate slot.
    _fill_adata_slots_automatically(adata, d)

    # Add a special additional field to .var if it exists.
    if 'features_analyzed_inds' in adata.uns.keys():
        adata.var['cellbender_analyzed'] = [True if (i in adata.uns['features_analyzed_inds'])
                                            else False for i in range(adata.shape[1])]

    if analyzed_barcodes_only:
        for col in adata.obs.columns[adata.obs.columns.str.startswith('barcodes_analyzed')
                                     | adata.obs.columns.str.startswith('barcode_indices')]:
            try:
                del adata.obs[col]
            except Exception:
                pass
    else:
        # Add a special additional field to .obs if all barcodes are included.
        if 'barcodes_analyzed_inds' in adata.uns.keys():
            adata.obs['cellbender_analyzed'] = [True if (i in adata.uns['barcodes_analyzed_inds'])
                                                else False for i in range(adata.shape[0])]

    return adata


def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def _fill_adata_slots_automatically(adata, d):
    """Add other information to the adata object in the appropriate slot."""

    for key, value in d.items():
        try:
            if value is None:
                continue
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == adata.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == adata.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print('Unable to load data into AnnData: ', key, value, type(value))


def cellbender_to_scrublet(data_path, scr_out_path, cb_z_out_path, h5ad_out_path, expected_dblt, sample_name):
    # load the data
    z = []
    with tables.open_file(data_path) as f:
        print(f)  # display the structure of the h5 file
        z = f.root.matrix.latent_gene_encoding.read() # read latents
    np.savetxt(cb_z_out_path, z, delimiter=',') #     tiny_10x_pbmc_latent_gene_expression.csv

    adata = anndata_from_h5(data_path)
    adata.uns["name"] = sample_name
    adata.var_names_make_unique()
    # load the latent representation into a new slot called 'X_cellbender'
    adata.obsm['X_cellbender'] = z

    scrub = scr.Scrublet(adata.X, expected_doublet_rate = expected_dblt)
    adata.obs['doublet_score'], adata.obs['predicted_doublets'] =     scrub.scrub_doublets(min_counts=2, min_cells=3,     min_gene_variability_pctl=85, n_prin_comps=30)
    embedding = pacmap.PaCMAP(
        n_components=2,
        n_neighbors=None,
        MN_ratio=0.5,
        FP_ratio=2.0,
        apply_pca=False)
    adata.obsm['X_pacmap'] = embedding.fit_transform(
        adata.obsm['X_cellbender'],
        init="pca")
    sc.pl.embedding(
        adata,
        basis='X_pacmap',
        color='doublet_score',
        title='PaCMAP: Doublets score derived using Scrublet in {sample}'.format(sample=adata.uns["name"]),
        save="_doublet-score_{sample}.pdf".format(sample=adata.uns["name"]))
    sc.tl.tsne(adata,
        use_rep="X_cellbender",
        n_jobs=n_cores,
        random_state=reseed,
        perplexity=30,
        metric='euclidean')
    sc.pl.tsne(adata,
        color='doublet_score',
        title='tSNE: Doublets score derived using Scrublet in {sample}'.format(sample=adata.uns["name"]),
        save="_doublet-score_{sample}.pdf".format(sample=adata.uns["name"]))
    # Save results:
    adata.write(h5ad_out_path)
    pd.DataFrame(adata.obs).to_csv(scr_out_path, sep = '\t', header = True) # scrublet_calls.tsv



cellbender_to_scrublet(data_path=snakemake.input["filt_h5"], scr_out_path=snakemake.output["scrublet_calls"], cb_z_out_path=snakemake.output["dr"], h5ad_out_path=snakemake.output["h5ad"], expected_dblt=snakemake.params["expected_dblt"], sample_name=snakemake.params["sample_run_name"])