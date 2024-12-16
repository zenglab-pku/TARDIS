import pandas as pd
import numpy as np
import scanpy as sc

import seaborn as sns
import matplotlib.pyplot as plt

import itertools

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from typing import Optional
from scipy.sparse import csr_matrix # type: ignore
from anndata import AnnData

import cv2

from typing import (
    Optional
)

def qc_guide_bins(
    gem_path: str,
    guide_prefix: Optional[str] = None,
    fig_path: Optional[str] = None
) -> None:
    df = pd.read_csv(gem_path, header=0, index_col=None, sep='\t', comment='#')
    if df.columns[0] != 'geneID':
        df.set_index(df.columns[0])
    if guide_prefix != None: df = df[df.geneID.str.startswith(guide_prefix)]
    dup_df = df[df.duplicated(subset=['x', 'y'], keep=False)]
    dedup_df = dup_df.drop_duplicates(subset=['x', 'y'])
    single_df = df.drop_duplicates(subset=['x', 'y'], keep=False)
    x = ['Total CID', 'Singlet CID', 'Doublet CID', 'Dedup CID', 'Combined CID']
    y = [df.shape[0], single_df.shape[0], dup_df.shape[0], dedup_df.shape[0], dedup_df.shape[0] + single_df.shape[0]]
    plt.figure(figsize=(4, 3))
    sns.barplot(x=x, y=y, hue=[1, 2, 3, 4, 5], palette='tab20b', alpha=0.5, legend=False)
    for i in range(5):
        plt.text(x[i], y[i] + (max(y) / 200), y[i], ha='center')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Count')
    plt.yticks([])
    sns.despine()
    if fig_path != None:
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()

def filter_guide_reads(
    gem_path: str,
    guide_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    binarilize: bool = False,
    assign_pattern: Optional[str] = 'max',
    filter_threshold: Optional[int] = None
) -> pd.DataFrame:
    df = pd.read_csv(gem_path, header=0, index_col=None, sep='\t', comment='#')
    if df.columns[0] != 'geneID':
        df.set_index(df.columns[0], drop=True, inplace=True)
    if guide_prefix != None: df = df[df.geneID.str.startswith(guide_prefix)]
    if filter_threshold != None:
        df = df[df.MIDCount > filter_threshold]
    single_df = df.drop_duplicates(subset=['x', 'y'], keep=False)
    dup_df = df[df.duplicated(subset=['x', 'y'], keep=False)]
    if assign_pattern == 'max':
        max_counts = dup_df.groupby(['x', 'y'])['MIDCount'].transform('max')
        dedup_df = dup_df[dup_df['MIDCount'] == max_counts]
    elif assign_pattern == 'drop':
        dedup_df = pd.DataFrame(columns=dup_df.columns)
    elif assign_pattern == 'all':
        dedup_df = dup_df
    else:
        print('Error: Assign pattern invalid.')
        return
    output_df = pd.concat([single_df, dedup_df], axis=0)
    if binarilize is True:
        output_df['MIDCount'] = 1
        output_df['ExonCount'] = 1
    if output_path == None:
        return output_df
    else:
        output_df.to_csv(filter_threshold, index=False, header=True, sep='\t')
        return

def load_bin(
    gem_file: str,
    bin_size: int,
    library_id: str,
    image_file: Optional[str] = None,
) -> AnnData:
    """
    Read BGI data and image file, and return an AnnData object.
    Parameters
    ----------
    gem_file
        The path of the BGI data file.
    image_file
        The path of the image file.
    bin_size
        The size of the bin.
    library_id
        The library id.
    Returns
    -------
    Annotated data object with the following keys:

        - :attr:`anndata.AnnData.obsm` ``['spatial']`` - spatial spot coordinates.
        - :attr:`anndata.AnnData.uns` ``['spatial']['{library_id}']['images']`` - *hires* images.
        - :attr:`anndata.AnnData.uns` ``['spatial']['{library_id}']['scalefactors']`` - scale factors for the spots.
    """
    library = library_id
    dat_file = gem_file
    image = image_file
    bin_s = bin_size

    dat = pd.read_csv(dat_file, delimiter="\t", comment="#")
    
    dat['x'] -= dat['x'].min()
    dat['y'] -= dat['y'].min()

    width = dat['x'].max() + 1
    height = dat['y'].max() + 1

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
    if image_file:
        image = cv2.imread(image)
    else: 
        image = None
    if image is not None:

        spatial_key = "spatial"
        adata.uns[spatial_key] = {library: {}}
        adata.uns[spatial_key][library]["images"] = {}
        adata.uns[spatial_key][library]["images"] = {"hires": image}

        tissue_hires_scalef = max(image.shape[0]/width, image.shape[1]/height)

        spot_diameter = bin_s / tissue_hires_scalef
    
        adata.uns[spatial_key][library]["scalefactors"] = {
            "tissue_hires_scalef": tissue_hires_scalef,
            "spot_diameter_fullres": spot_diameter,
        }

    return adata

def combine_guide_replicates(gdata):
    """
    Usage:
        combine guide replicates in <gdata> to a single gene name
    Return:
        single gene name anndata
    """
    sgs = gdata.var_names.str.split('_', n=1).str[0]
    sgs_grouped = pd.DataFrame(gdata.X.toarray(), columns=gdata.var_names)
    sgs_grouped = sgs_grouped.groupby(sgs, axis=1).sum()

    cgdata = AnnData(sgs_grouped, obs=gdata.obs, var=pd.DataFrame(index=sgs_grouped.columns))
    cgdata.obsm['spatial'] = gdata.obsm['spatial']
    return cgdata

def nmf_clustering(
    adata: AnnData,
    n_components: int = 10,
    random_state: int = 42,
    max_iter: int = 1000,
    verbose: int = 0,
    n_top_genes: int = 2000,
) -> AnnData:
    from sklearn.decomposition import NMF

    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    hvg_data = adata[:, adata.var['highly_variable']]
    cnt_matrix = hvg_data.X.toarray()
    model = NMF(n_components=n_components, verbose=verbose, random_state=random_state, max_iter=max_iter)
    cnt_matrix_trans = model.fit_transform(cnt_matrix)
    hvg_data.obsm['X_nmf'] = cnt_matrix_trans
    hvg_data.uns['X_nmf_components'] = model.components_
    return hvg_data

def nmf_consensus(
    adata: AnnData,
    min_clusters: Optional[int] = 4,
    max_clusters: Optional[int] = 10,
    n_resamples: Optional[int] = 100,
    resample_frac: Optional[float] = 0.8,
    random_state: Optional[int] = 42,
    n_cluster_genes: Optional[int] = 50,
) -> AnnData:
    from tqdm import tqdm
    from scipy.stats import pearsonr

    cnt_matrix_trans = adata.obsm['X_nmf']
    corr_matrix = np.zeros((cnt_matrix_trans.shape[1], cnt_matrix_trans.shape[1]))
    for i in tqdm(range(cnt_matrix_trans.shape[1])):
        for j in range(i, cnt_matrix_trans.shape[1]):
            corr_matrix[i, j] = pearsonr(cnt_matrix_trans[:, i], cnt_matrix_trans[:, j])[0]
    for i in tqdm(range(corr_matrix.shape[1])):
        for j in range(i):
            corr_matrix[i, j] = corr_matrix[j, i]

    import consensusclustering as cc
    from sklearn.cluster import AgglomerativeClustering

    cc_model = cc.ConsensusClustering(
        AgglomerativeClustering(),
        min_clusters=min_clusters,
        max_clusters=max_clusters,
        n_resamples=n_resamples,
        resample_frac=resample_frac,
        k_param='n_clusters'
    )
    cc_model.fit(corr_matrix)
    k = cc_model.best_k()
    cc_model.plot_clustermap(k)

    model = AgglomerativeClustering(n_clusters=k)
    model.fit(corr_matrix)
    nmf_clusters = model.labels_

    for i in range(k):
        cluster_index = np.where(nmf_clusters == i)[0]
        genes = adata.var_names[np.argsort(adata.uns['X_nmf_components'][cluster_index, :]).flatten()[-n_cluster_genes:]].unique()
        sc.tl.score_genes(adata, genes, score_name=f'nmf_cluster_{i}')
        adata.obs[f'nmf_cluster_{i}'] = (adata.obs[f'nmf_cluster_{i}'] - np.mean(adata.obs[f'nmf_cluster_{i}'])) / np.std(adata.obs[f'nmf_cluster_{i}'])
    
    return adata
