from typing import Optional
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData
from scipy import sparse
from matplotlib import cm
import matplotlib as mpl

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

def plot_guide_gene_summary(adata):
    """
    Plots a summary scatterplot for guide data showing, for each guide/gene:
      - Number of spots detected (x-axis)
      - Total expression rank (y-axis, inverted)
      - Mean detection (color)
    Parameters
    ----------
    adata : AnnData
        AnnData object for guides or genes. .X must be cells x features.
    Returns
    -------
    None
    """
    X = adata.X
    if sparse.issparse(X):
        X = X.tocsc()
        expr_count = np.array((X > 0).sum(axis=0)).ravel()
        total_expr = np.array(X.sum(axis=0)).ravel()
        mean_expr = np.array(X.mean(axis=0)).ravel()
    else:
        expr_count = (X > 0).sum(axis=0)
        total_expr = X.sum(axis=0)
        mean_expr = X.mean(axis=0)

    gene_rank = total_expr.argsort().argsort() + 1
    gene_rank = len(total_expr) - gene_rank + 1

    log_mean_expr = np.log1p(mean_expr)
    norm = mpl.colors.LogNorm(
        vmin=log_mean_expr[log_mean_expr > 0].min(),
        vmax=log_mean_expr.max()
    )

    plt.figure(figsize=(5, 3))
    scatter = plt.scatter(
        expr_count, gene_rank,
        c=mean_expr, cmap='viridis', norm=norm,
        s=10, alpha=0.8
    )
    plt.xlabel('Number of spots detected')
    plt.ylabel('Gene total expression rank')
    plt.title('')
    plt.yticks([])
    plt.xscale('log')
    plt.colorbar(scatter, label='Mean detection')
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.show()

def remove_mito_ribo_hk_lnc_genes(adata, housekeeping_list="He2020Nature_mouseHK.txt"):
    """
    Usage:
        strip <adata> with genes names beginning with "Mt", "mt-", "Rp", "Gm" and ending with "Rik" or "Rik#";
        strip also with in house housekeeping gene list of mouse housekeeping genes
    Returns:
        clean anndata object without the genes above.
    """
    return_data = adata.copy()
    return_data.var["mt"] = return_data.var_names.str.startswith("Mt")
    return_data.var["mt-"] = return_data.var_names.str.startswith("mt-")
    return_data.var["gm"] = return_data.var_names.str.startswith("Gm")
    return_data.var["Rb"] = return_data.var_names.str.startswith("Rp")
    return_data.var["rik"] = [True if "Rik" in str else False for str in return_data.var_names]
    return_data = return_data[:, ~return_data.var["mt"]].copy()
    return_data = return_data[:, ~return_data.var["mt-"]]
    return_data = return_data[:, ~return_data.var["Rb"]]
    return_data = return_data[:, ~return_data.var["gm"]]
    return_data = return_data[:, ~return_data.var["rik"]]

    with open(housekeeping_list, 'r') as f:
        for line in f:
            hk_genes = line.split('\t')
            break
    return_data = return_data[:, [gene for gene in return_data.var_names if gene not in hk_genes]]
    del return_data.var
    return return_data

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
def _get_spatial_xy(adata: AnnData) -> pd.DataFrame:
    """Resolve (x, y) per obs for duplicate detection."""
    if "spatial" in adata.obsm:
        spatial = np.asarray(adata.obsm["spatial"])
        if spatial.shape[1] >= 3:
            y, x = spatial[:, 1], spatial[:, 2]
        elif spatial.shape[1] == 2:
            y, x = spatial[:, 0], spatial[:, 1]
        else:
            raise ValueError("obsm['spatial'] must have 2 or 3 columns.")
    elif {"pxl_row_in_fullres", "pxl_col_in_fullres"}.issubset(adata.obs.columns):
        y = adata.obs["pxl_row_in_fullres"].to_numpy()
        x = adata.obs["pxl_col_in_fullres"].to_numpy()
    elif {"array_row", "array_col"}.issubset(adata.obs.columns):
        y = adata.obs["array_row"].to_numpy()
        x = adata.obs["array_col"].to_numpy()
    else:
        raise ValueError(
            "Cannot infer spatial coordinates. Provide obsm['spatial'] "
            "or obs columns (pxl_row_in_fullres/pxl_col_in_fullres or array_row/array_col)."
        )
    return pd.DataFrame({"barcode": adata.obs_names.astype(str), "x": x, "y": y})
    
def _adata_to_guide_long(adata: AnnData) -> pd.DataFrame:
    """Convert guide AnnData to GEM-like long format."""
    X = adata.X
    if sparse.issparse(X):
        X = X.tocoo()
        if X.nnz == 0:
            df = pd.DataFrame(columns=["geneID", "barcode", "MIDCount"])
        else:
            df = pd.DataFrame({
                "geneID": adata.var_names[X.col],
                "barcode": adata.obs_names[X.row].astype(str),
                "MIDCount": X.data,
            })
    else:
        idx = np.where(X > 0)
        df = pd.DataFrame({
            "geneID": adata.var_names[idx[1]],
            "barcode": adata.obs_names[idx[0]].astype(str),
            "MIDCount": X[idx],
        })
    if df.empty:
        return df.assign(x=pd.Series(dtype=float), y=pd.Series(dtype=float))
    spatial = _get_spatial_xy(adata)
    df = df.merge(spatial, on="barcode", how="left")
    return df

def _long_to_adata(
    df: pd.DataFrame,
    template: AnnData,
    guide_names: pd.Index,
) -> AnnData:
    """Rebuild AnnData from filtered long table."""
    if df.empty:
        out = template[:, guide_names].copy()
        out.X = sparse.csr_matrix(out.shape, dtype=np.float32)
        return out
    barcodes = df["barcode"].astype(str).unique()
    gene_to_idx = {g: i for i, g in enumerate(guide_names)}
    barcode_to_idx = {b: i for i, b in enumerate(barcodes)}
    rows = df["barcode"].map(barcode_to_idx).to_numpy()
    cols = df["geneID"].map(gene_to_idx).to_numpy()
    data = df["MIDCount"].to_numpy()
    X = sparse.csr_matrix(
        (data, (rows, cols)),
        shape=(len(barcodes), len(guide_names)),
        dtype=np.float32,
    )
    obs = template.obs.loc[template.obs_names.astype(str).isin(barcodes)].copy()
    obs.index = obs.index.astype(str)
    obs = obs.loc[barcodes]
    out = AnnData(X=X, obs=obs, var=template.var.loc[guide_names].copy())
    out.var_names = guide_names.astype(str)
    out.obs_names = pd.Index(barcodes, dtype=str)
    for key in template.obsm:
        obs_idx = template.obs_names.astype(str).isin(barcodes)
        obs_order = template.obs_names.astype(str)[obs_idx]
        reorder = obs_order.get_indexer(barcodes)
        out.obsm[key] = template.obsm[key][obs_idx][reorder]
    for key in template.uns:
        out.uns[key] = template.uns[key]
    return out

def _filter_guide_long_df(
    df: pd.DataFrame,
    guide_prefix: Optional[str] = None,
    binarilize: bool = False,
    assign_pattern: str = "max",
    filter_threshold: Optional[int] = None,
) -> pd.DataFrame:
    """Core filtering logic shared by GEM and h5 paths."""
    if df.empty:
        return df
    if guide_prefix is not None:
        df = df[df["geneID"].str.startswith(guide_prefix)]
    if filter_threshold is not None:
        df = df[df["MIDCount"] > filter_threshold]
    if df.empty:
        return df
    single_df = df.drop_duplicates(subset=["x", "y"], keep=False)
    dup_df = df[df.duplicated(subset=["x", "y"], keep=False)]
    if assign_pattern == "max":
        max_counts = dup_df.groupby(["x", "y"])["MIDCount"].transform("max")
        dedup_df = dup_df[dup_df["MIDCount"] == max_counts]
    elif assign_pattern == "drop":
        dedup_df = dup_df.iloc[0:0]
    elif assign_pattern == "all":
        dedup_df = dup_df
    else:
        raise ValueError(f"Invalid assign_pattern: {assign_pattern!r}. Use 'max', 'drop', or 'all'.")
    output_df = pd.concat([single_df, dedup_df], axis=0)
    if binarilize:
        output_df = output_df.copy()
        output_df["MIDCount"] = 1
    return output_df

def filter_guide_reads_h5(
    adata_or_path: Union[str, Path, AnnData],
    guide_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    binarilize: bool = False,
    assign_pattern: Optional[str] = "max",
    filter_threshold: Optional[int] = None,
    copy: bool = True,
) -> Optional[AnnData]:
    """
    AnnData/h5ad version of filter_guide_reads.
    Applies the same filtering logic as filter_guide_reads on a guide count matrix.
    """
    if isinstance(adata_or_path, (str, Path)):
        adata = sc.read_h5ad(str(adata_or_path))
    else:
        adata = adata_or_path.copy() if copy else adata_or_path
    long_df = _adata_to_guide_long(adata)
    filtered_df = _filter_guide_long_df(
        long_df,
        guide_prefix=guide_prefix,
        binarilize=binarilize,
        assign_pattern=assign_pattern or "max",
        filter_threshold=filter_threshold,
    )
    if guide_prefix is not None:
        guide_names = adata.var_names[adata.var_names.str.startswith(guide_prefix)]
    else:
        guide_names = adata.var_names
    out = _long_to_adata(filtered_df, adata, guide_names)
    if output_path is None:
        return out
    out.write_h5ad(output_path)
    return None
    
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

def dbscan_density_region(
    adata,
    guide,
    eps=5,
    min_samples=20,
    label="0",
    mode="most",
    density_level=7,
    step=1,
    bandwidth=2,
    region_label_columns=("pxl_row_in_fullres", "pxl_col_in_fullres"),
):
    """
    Perform DBSCAN clustering on the specified guide and extract the main high-density region.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing spatial information and expression matrix.
    guide : str
        The name of the guide gene of interest.
    eps : float, default=5
        The eps parameter for the DBSCAN algorithm, determining the minimum distance between points to be considered neighbors.
    min_samples : int, default=20
        Minimum number of samples in a neighborhood for a point to be considered a core point in DBSCAN.
    label : str, default="0"
        Cluster label to extract the region for. Ignored if mode="most".
    mode : {"most", "label"}, default="most"
        Region extraction mode: "most" selects the largest cluster automatically, "label" uses the specified label.
    density_level : int, default=7
        Number of density levels used to select the density threshold.
    step : int, default=1
        Grid interval for density estimation; smaller values give finer resolution but are slower.
    bandwidth : float, default=2
        Bandwidth for KernelDensity.
    region_label_columns : tuple, default=("pxl_row_in_fullres", "pxl_col_in_fullres")
        Column names for the returned DataFrame.

    Returns
    -------
    guide_adata : AnnData
        Subset of the AnnData object containing only samples of the queried guide with an added DBSCAN cluster label.
    core_points_df : DataFrame
        Coordinates of points in the core high-density region, with columns as specified by region_label_columns.
    """
    from sklearn.cluster import DBSCAN
    from sklearn.neighbors import KernelDensity
    import numpy as np
    import pandas as pd

    guide_adata = adata[adata[:, guide].X.toarray().flatten() > 0].copy()
    guide_dots = guide_adata.obsm['spatial']

    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(guide_dots)
    cluster_labels = clustering.labels_
    print("Number of Clusters (excluding noise):",
          len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0))
    guide_adata.obs["dbscan_cluster"] = "Noise"
    guide_adata.obs.loc[:, "dbscan_cluster"] = cluster_labels.astype(str)

    labels = np.asarray(guide_adata.obs["dbscan_cluster"])
    if mode == "most":
        filtered_labels = labels[labels != "-1"]
        unique_labels, counts = np.unique(filtered_labels, return_counts=True)
        if len(counts) == 0:
            raise ValueError("No clusters (other than noise) found.")
        max_count_idx = np.argmax(counts)
        region_label = unique_labels[max_count_idx]
    else:
        region_label = str(label)

    guide_dots_region = np.asarray(guide_adata.obsm['spatial'][:, 1:], dtype=np.float32)
    sel_mask = (labels == region_label)
    if not np.any(sel_mask):
        raise ValueError(f"No points found for label {region_label}. Function will terminate.")
    cluster_points = guide_dots_region[sel_mask]

    x_min, x_max = cluster_points[:, 0].min() - 1, cluster_points[:, 0].max() + 1
    y_min, y_max = cluster_points[:, 1].min() - 1, cluster_points[:, 1].max() + 1
    x_range = np.arange(x_min, x_max + 1, step, dtype=np.float32)
    y_range = np.arange(y_min, y_max + 1, step, dtype=np.float32)
    X, Y = np.meshgrid(x_range, y_range)
    grid_coords = np.vstack([X.ravel(), Y.ravel()]).T

    kde = KernelDensity(bandwidth=bandwidth, kernel="gaussian")
    kde.fit(cluster_points)
    Z = np.exp(kde.score_samples(grid_coords))
    Z = Z.reshape(X.shape)
    density_levels = np.linspace(Z.min(), Z.max(), density_level)
    level_density = density_levels[2]
    core_points_mask = Z >= level_density
    core_points = np.column_stack((X[core_points_mask], Y[core_points_mask]))
    core_points_df = pd.DataFrame(core_points, columns=region_label_columns)
    return guide_adata, core_points_df
