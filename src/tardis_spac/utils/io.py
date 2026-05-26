import itertools
import pandas as pd
import numpy as np

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from scipy.sparse import csr_matrix
from anndata import AnnData

def load_bin(
    gem_file: str,
    bin_size: int,
) -> AnnData:
    """
    Read BGI data and image file, and return an AnnData object.
    Parameters
    ----------
    gem_file
        The path of the BGI data file.
    bin_size
        The size of the bin.
    Returns
    -------
    Annotated data object with the following keys:
        - :attr:`anndata.AnnData.obsm` ``['spatial']`` - spatial spot coordinates.
    """
    bin_s = bin_size
    dat_file = gem_file
    
    dat = pd.read_csv(dat_file, delimiter="\t", comment="#")
    
    dat['x'] -= dat['x'].min()
    dat['y'] -= dat['y'].min()

    dat['xp'] = (dat['x'] // bin_s) * bin_s
    dat['yp'] = (dat['y'] // bin_s) * bin_s
    dat['xb'] = np.floor(dat['xp'] / bin_s + 1).astype(int)
    dat['yb'] = np.floor(dat['yp'] / bin_s + 1).astype(int)

    dat['bin_ID'] = max(dat['xb']) * (dat['yb'] - 1) + dat['xb']

    trans_x_xb = dat[['x', 'xb']].drop_duplicates()
    trans_x_xb = trans_x_xb.groupby('xb')['x'].apply(
        lambda x: int(np.floor(np.mean(x)))).reset_index()
    trans_y_yb = dat[['y', 'yb']].drop_duplicates()
    trans_y_yb = trans_y_yb.groupby('yb')['y'].apply(
        lambda y: int(np.floor(np.mean(y)))).reset_index()

    trans_matrix = pd.DataFrame(list(itertools.product(
        trans_x_xb['xb'], trans_y_yb['yb'])), columns=['xb', 'yb'])
    trans_matrix = pd.merge(trans_matrix, trans_x_xb, on='xb')
    trans_matrix = pd.merge(trans_matrix, trans_y_yb, on='yb')
    trans_matrix['bin_ID'] = max(
        trans_matrix['xb']) * (trans_matrix['yb'] - 1) + trans_matrix['xb']

    trans_matrix['in_tissue'] = 1

    tissue_positions = pd.DataFrame()
    # barcode is str, not number
    tissue_positions['barcodes'] = trans_matrix['bin_ID'].astype(str)
    tissue_positions['in_tissue'] = trans_matrix['in_tissue']
    tissue_positions['array_row'] = trans_matrix['yb']
    tissue_positions['array_col'] = trans_matrix['xb']
    tissue_positions['pxl_row_in_fullres'] = trans_matrix['y']
    tissue_positions['pxl_col_in_fullres'] = trans_matrix['x']
    tissue_positions.set_index('barcodes', inplace=True)

    if 'MIDCount' in dat.columns:
        dat = dat.groupby(['geneID', 'xb', 'yb'])[
            'MIDCount'].sum().reset_index()
        dat['bin_ID'] = max(dat['xb']) * (dat['yb'] - 1) + dat['xb']

        unique_genes = dat['geneID'].unique()
        unique_barcodes = dat['bin_ID'].unique()
        gene_hash = {gene: index for index, gene in enumerate(unique_genes)}
        barcodes_hash = {barcodes: index for index,
                         barcodes in enumerate(unique_barcodes)}
        dat['gene'] = dat['geneID'].map(gene_hash)
        dat['barcodes'] = dat['bin_ID'].map(barcodes_hash)

        counts = csr_matrix((dat['MIDCount'], (dat['barcodes'], dat['gene'])))

    else:
        dat = dat.groupby(['geneID', 'xb', 'yb'])[
            'MIDCounts'].sum().reset_index()
        dat['bin_ID'] = max(dat['xb']) * (dat['yb'] - 1) + dat['xb']

        unique_genes = dat['geneID'].unique()
        unique_barcodes = dat['bin_ID'].unique()
        gene_hash = {gene: index for index, gene in enumerate(unique_genes)}
        barcodes_hash = {barcodes: index for index,
                         barcodes in enumerate(unique_barcodes)}
        dat['gene'] = dat['geneID'].map(gene_hash)
        dat['barcodes'] = dat['bin_ID'].map(barcodes_hash)

        counts = csr_matrix((dat['MIDCounts'], (dat['barcodes'], dat['gene'])))
    adata = AnnData(counts)
    adata.var_names = list(gene_hash.keys())
    adata.obs_names = list(map(str, barcodes_hash.keys()))

    adata.obs = adata.obs.join(tissue_positions, how="left")
    adata.obsm['spatial'] = adata.obs[[
        'pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()
    adata.obs.drop(columns=['in_tissue', 'array_row', 'array_col',
                   'pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True,)

    return adata
