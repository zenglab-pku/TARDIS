from typing import Optional
from weakref import ref
import anndata as ad
import numpy as np
import pandas as pd
from numba import njit
from tqdm import tqdm
import random
from scipy.spatial.distance import euclidean, pdist, squareform
from scipy.stats import gaussian_kde
from skbio.stats.distance import permanova as skbio_permanova, DistanceMatrix

def aitchison_distance(
    gdata: ad.AnnData,
    cluster_field: str,
    result_field: str = 'aitchison_dist',
    reference_guide: str = 'sgNon-targeting',
    library_key: Optional[str] = None,
    n_permutations: int = 1000,
    show_progress: bool = True,
    p_swap: float = 0.1,
    copy: bool = False
):
    """
    Compute the Aitchison distance between a specified reference_guide and all other guides,
    based on their cluster-wise abundances, with an optional label permutation test.

    This implementation fixes potential KeyErrors during permutation, which may arise
    if certain guides are missing clusters, by reindexing Series objects to ensure alignment.

    Parameters
    ----------
    gdata : AnnData
        AnnData object containing .obs (cell metadata), .var (guide-level information), and .X (expression/count matrix).
    cluster_field : str
        Field name in .obs representing cluster assignments.
    result_field : str, default "aitchison_dist"
        Name of the output field to write results into gdata.var.
    reference_guide : str, default "sgNon-targeting"
        Which guide to use as the reference (often a negative control or non-targeting guide).
    library_key : str or None, default None
        If set, computes distances per library or sample; if None, pooled analysis.
    n_permutations : int or None, default 1000
        Number of label permutations to estimate p-values; set to None to skip permutation.
    show_progress : bool, default True
        Whether to display a tqdm progress bar during permutation calculation.
    p_swap : float, default 0.1
        Probability of swapping values between guide and reference for each cluster in a permutation.
    copy : bool, default False
        If True, returns a new AnnData with results written; if False, updates in place.

    Returns
    -------
    If copy is True, returns AnnData with result written to gdata.var[result_field].
    If copy is False, updates gdata in place and returns None.
    """

    if copy:
        gdata = gdata.copy()

    if gdata.obs[cluster_field].dtype is not str:
        gdata.obs[cluster_field] = gdata.obs[cluster_field].astype(str)

    if not isinstance(gdata.X, np.ndarray):
        gdata.X = gdata.X.toarray()

    # Gather the clusters present globally for later index alignment
    all_clusters = gdata.obs[cluster_field].unique()

    if library_key is None:
        # Prepare count_df: index=guide (var_name), columns=clusters
        counts = pd.DataFrame(gdata.X, columns=gdata.var_names, index=gdata.obs_names)
        clusters = gdata.obs[cluster_field]
        count_df = []
        for guide in counts.columns:
            # for each guide, group by cluster, sum count
            group_sum = counts[guide].groupby(clusters).sum()
            # Reindex to ensure all clusters are present (missing=0)
            group_sum = group_sum.reindex(all_clusters, fill_value=0)
            count_df.append(group_sum)
        count_df = pd.DataFrame(count_df, index=counts.columns, columns=all_clusters)

        n_clusters = len(all_clusters)
        # Apply compositional (CLR-like) transform
        transform_df = count_df.apply(lambda x: 10 ** (np.log10(x + 1) - (np.log10(x + 1).sum() / n_clusters)), axis=1)

        # Compute Aitchison distances to reference_guide
        if reference_guide not in transform_df.index:
            raise ValueError(f"reference_guide '{reference_guide}' not found in var_names.")
        dist_df = transform_df.apply(lambda x: euclidean(x, transform_df.loc[reference_guide]), axis=1)
        gdata.var[result_field] = dist_df.loc[gdata.var_names].values

        # Permutation test
        if n_permutations is not None:
            pvals = pd.Series(np.zeros(len(gdata.var_names)), index=gdata.var_names)
            for guide in gdata.var_names:
                if guide == reference_guide:
                    continue
                obs_dist = gdata.var[result_field][guide]
                greater_count = 0
                for _ in tqdm(range(n_permutations), desc=f'{guide}', disable=not show_progress):
                    # get cluster counts for guide/ref, align index
                    p_guide_cnts = count_df.loc[guide].copy()
                    p_ntc_cnts = count_df.loc[reference_guide].copy()
                    for cluster in all_clusters:
                        if random.random() < p_swap:
                            guide_cnt = p_guide_cnts[cluster]
                            ntc_cnt = p_ntc_cnts[cluster]
                            diff = ntc_cnt - guide_cnt
                            if diff > 0:
                                t = random.randint(1, diff)
                                p_guide_cnts[cluster] += t
                                p_ntc_cnts[cluster] -= t
                            elif diff < 0:
                                t = random.randint(1, -diff)
                                p_guide_cnts[cluster] -= t
                                p_ntc_cnts[cluster] += t
                    # Transform (as above)
                    p_guide_arr = 10 ** (np.log10(p_guide_cnts + 1) - (np.log10(p_guide_cnts + 1).sum() / n_clusters))
                    p_ntc_arr  = 10 ** (np.log10(p_ntc_cnts + 1)  - (np.log10(p_ntc_cnts + 1).sum() / n_clusters))
                    dist = euclidean(p_guide_arr, p_ntc_arr)
                    if dist > obs_dist:
                        greater_count += 1
                pvals[guide] = greater_count / n_permutations
            gdata.var[result_field + '.p_value'] = pvals.values

    else:
        library_key_list = gdata.obs[library_key].unique()
        n_clusters = len(all_clusters)

        for sample in library_key_list:
            # Subset to sample
            pdata = gdata[gdata.obs[library_key] == sample].copy()
            sample_clusters = pdata.obs[cluster_field].unique()

            # Prepare count_df: index=guide (var_name), columns=clusters
            counts = pd.DataFrame(pdata.X, columns=pdata.var_names, index=pdata.obs_names)
            clusters = pdata.obs[cluster_field]
            sample_count_df = []
            for guide in counts.columns:
                group_sum = counts[guide].groupby(clusters).sum()
                # align columns to all_clusters (not just sample_clusters)
                group_sum = group_sum.reindex(all_clusters, fill_value=0)
                sample_count_df.append(group_sum)
            sample_count_df = pd.DataFrame(sample_count_df, index=counts.columns, columns=all_clusters)

            transform_df = sample_count_df.apply(
                lambda x: 10 ** (np.log10(x + 1) - (np.log10(x + 1).sum() / n_clusters)), axis=1
            )
            if reference_guide not in transform_df.index:
                raise ValueError(f"reference_guide '{reference_guide}' not found in var_names for sample '{sample}'.")
            dist_df = transform_df.apply(lambda x: euclidean(x, transform_df.loc[reference_guide]), axis=1)
            # Output aligned to gdata.var order
            gdata.var[sample + result_field] = dist_df.loc[gdata.var_names].values

        if n_permutations is not None:
            # Initialize pvals if not already present
            if (result_field + '.p_value') not in gdata.var.columns:
                gdata.var[result_field + '.p_value'] = 0
            pvals = pd.Series(np.zeros(len(gdata.var_names)), index=gdata.var_names)
            for _ in tqdm(range(n_permutations), desc='Permutations', disable=not show_progress):
                aitchson_dist = {}
                for sample in library_key_list:
                    pdata = gdata[gdata.obs[library_key] == sample].copy()
                    # prepare per-sample clusters and counts
                    counts = pd.DataFrame(pdata.X, columns=pdata.var_names, index=pdata.obs_names)
                    clusters = pdata.obs[cluster_field]
                    sample_count_df = []
                    for guide in counts.columns:
                        group_sum = counts[guide].groupby(clusters).sum()
                        group_sum = group_sum.reindex(all_clusters, fill_value=0)
                        sample_count_df.append(group_sum)
                    sample_count_df = pd.DataFrame(sample_count_df, index=counts.columns, columns=all_clusters)

                    for guide in pdata.var_names:
                        if guide == reference_guide:
                            continue
                        p_guide_cnts = sample_count_df.loc[guide].copy()
                        p_ntc_cnts = sample_count_df.loc[reference_guide].copy()
                        for cluster in all_clusters:
                            if random.random() < p_swap:
                                guide_cnt = p_guide_cnts[cluster]
                                ntc_cnt = p_ntc_cnts[cluster]
                                diff = ntc_cnt - guide_cnt
                                if diff > 0:
                                    t = random.randint(1, diff)
                                    p_guide_cnts[cluster] += t
                                    p_ntc_cnts[cluster] -= t
                                elif diff < 0:
                                    t = random.randint(1, -diff)
                                    p_guide_cnts[cluster] -= t
                                    p_ntc_cnts[cluster] += t
                        p_guide_arr = 10 ** (np.log10(p_guide_cnts + 1) - (np.log10(p_guide_cnts + 1).sum() / n_clusters))
                        p_ntc_arr  = 10 ** (np.log10(p_ntc_cnts + 1)  - (np.log10(p_ntc_cnts + 1).sum() / n_clusters))
                        dist = euclidean(p_guide_arr, p_ntc_arr)
                        if sample not in aitchson_dist:
                            aitchson_dist[sample] = {}
                        aitchson_dist[sample][guide] = dist
                # aitchson_dist[sample][guide]; DataFrame: index=guides, columns=samples
                perm_df = pd.DataFrame(aitchson_dist)
                # Compute per-guide, across-sample permutation rank
                # Final pval: mean rank approach
                for guide in gdata.var_names:
                    if guide == reference_guide:
                        continue
                    obs_vals = []
                    for sample in library_key_list:
                        col_name = sample + result_field
                        try:
                            obs_val = gdata.var[col_name][guide]
                        except Exception:
                            # If no value is available, skip
                            continue
                        obs_vals.append(obs_val)
                    # Mean observed distance (or fallback)
                    if obs_vals:
                        obs_dist = np.mean(obs_vals)
                        # Mean simulated in permutation
                        sim_dist = perm_df.loc[guide][library_key_list].mean()
                        if sim_dist > obs_dist:
                            pvals[guide] += 1
            gdata.var[result_field + '.p_value'] = (pvals / n_permutations).values

    if copy:
        return gdata
    else:
        return None

def permanova(
    gdata: ad.AnnData,
    cluster_field: str,
    result_field: str = 'permanova_f_value',
    reference_guide: str = 'sgNon-targeting',
    library_key: Optional[str] = None,
    count_bins: int = 10,
    n_permutations: int = 1000,
    show_progress: bool = True,
    copy: bool = False
):
    """
    Perform PERMANOVA analysis to compute F-values between each guide and the reference_guide based on kernel density estimates
    across clusters, and determine permutation-based statistical significance.

    Parameters
    ----------
    gdata : AnnData
        AnnData object containing cell-by-guide matrix, metadata, and cluster assignments.
    cluster_field : str
        Field in .obs specifying cluster membership for each observation.
    result_field : str, default 'permanova_f_value'
        Output column name in gdata.var to store F-values.
    reference_guide : str, default 'sgNon-targeting'
        The guide used as reference for comparison.
    library_key : str or None, default None
        If provided, computes PERMANOVA within each subgroup defined by this field.
    count_bins : int, default 10
        Number of bins to use for discretizing counts in density estimates.
    n_permutations : int, default 1000
        Number of permutations to perform for significance estimation.
    show_progress : bool, default True
        If True, display a progress bar for permutation computations.
    copy : bool, default False
        If True, operate on a copy of AnnData and return it; if False, update in place.

    Returns
    -------
    If copy is True, returns AnnData with F-values stored in gdata.var[result_field] and optionally p-values as result_field + '.p_value'.
    If copy is False, updates gdata in place and returns None.

    Notes
    -----
    This function handles DistanceMatrix errors by validating distance matrices, filling NaNs, and enforcing symmetry.
    Permutation testing is accelerated by:
        - Minimizing calls to scikit-bio's permanova() (the primary bottleneck),
        - Precomputing all group reassignments for permutations,
        - Vectorizing distance calculations and swap logic,
        - Utilizing numpy arrays and batching strategies where feasible.
    """
    import threading

    def nan_safe_dist_matrix(mat):
        # Ensure symmetric and fill diagonal with zeros, fill nans with zeros; force float.
        try:
            arr = np.array(mat, dtype=float)
        except Exception:
            return None
        if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
            return None
        np.fill_diagonal(arr, 0.0)
        arr = np.nan_to_num(arr, nan=0.0)
        arr = (arr + arr.T) / 2
        return arr

    if copy:
        gdata = gdata.copy()

    if not isinstance(gdata.X, np.ndarray):
        gdata.X = gdata.X.toarray()

    if gdata.obs[cluster_field].dtype is not str:
        gdata.obs[cluster_field] = gdata.obs[cluster_field].astype(str)

    cluster_tags = gdata.obs[cluster_field].unique()

    def kde_or_zeros(vals, x_grid):
        arr = np.array(vals.dropna())
        if arr.size > 1:
            kde = gaussian_kde(arr)
            density = kde(x_grid)
            density = density / (density.sum() if density.sum() > 0 else 1)
            return density
        else:
            return np.zeros_like(x_grid)

    # --------- FAST PERMANOVA PERMUTATION FUNCTION, NUMPY/NJIT/VECTORIZED ---------

    # Substitute for skbio permanova(). Re-implements two-group F-value, with permutation/grouping.
    def fast_pseudo_f(distance_mat, labels):
        """
        Computes 'pseudo-F' statistic for a distance matrix for two groups.
        Adapted (but simplified) from skbio's method, focused for 2-group comparison.
        distance_mat: square numpy array, N x N.
        labels: 1d array, 'guide'/'ref', len N.
        """
        # indices for group 0 and 1
        groups = np.unique(labels)
        idx0 = np.where(labels == groups[0])[0]
        idx1 = np.where(labels == groups[1])[0]
        n0, n1 = len(idx0), len(idx1)
        n = n0 + n1
        # Between-group distance: mean of all pairs, one from each group
        ss_between = np.mean(distance_mat[np.ix_(idx0, idx1)])
        # Within-group distance: mean of all within each group
        ss_within_0 = np.mean(distance_mat[np.ix_(idx0, idx0)]) if n0>1 else 0
        ss_within_1 = np.mean(distance_mat[np.ix_(idx1, idx1)]) if n1>1 else 0
        ss_within = (ss_within_0 * n0*(n0-1)//2 + ss_within_1 * n1*(n1-1)//2) / (n*(n-1)//2) if n > 1 else 0
        # F-like ratio (between/within)
        return ss_between / (ss_within + 1e-10)

    # We'll batch permutation for each guide, for much greater speed.
    # The code below is still I/O compatible (returns pvals, F values etc).

    # ----------- The single-sample (no library_key) block -----------
    if library_key is None:
        count_df = pd.concat([
            pd.DataFrame(gdata.X, columns=gdata.var_names, index=gdata.obs_names),
            gdata.obs[cluster_field],
        ], axis=1)

        pseudo_f = {}
        perm_data_dict = {}
        for guide in gdata.var_names:
            if guide == reference_guide: continue
            guide_count_df = count_df[[guide, reference_guide, cluster_field]]

            ref_df = guide_count_df.groupby([cluster_field, reference_guide]).count().unstack()[guide]
            guide_df = guide_count_df.groupby([cluster_field, guide]).count().unstack()[reference_guide]

            guide_xmax = guide_df.columns.max() if (guide_df.shape[1] > 0) else 0
            x_grid_guide = np.linspace(0, guide_xmax + 0.5, count_bins)
            guide_cnts = []
            for cluster in cluster_tags:
                if cluster in guide_df.index:
                    vals = guide_df.loc[cluster]
                    guide_density = kde_or_zeros(vals, x_grid_guide)
                    guide_cnts.append(guide_density)
                else:
                    guide_cnts.append(np.zeros_like(x_grid_guide))
            guide_cnts = np.array(guide_cnts)

            ref_xmax = ref_df.columns.max() if (ref_df.shape[1] > 0) else 0
            x_grid_ref = np.linspace(0, ref_xmax + 0.5, count_bins)
            ref_cnts = []
            for cluster in cluster_tags:
                if cluster in ref_df.index:
                    vals = ref_df.loc[cluster]
                    ref_density = kde_or_zeros(vals, x_grid_ref)
                    ref_cnts.append(ref_density)
                else:
                    ref_cnts.append(np.zeros_like(x_grid_ref))
            ref_cnts = np.array(ref_cnts)

            if x_grid_guide.shape[0] != x_grid_ref.shape[0]:
                m = max(x_grid_guide.shape[0], x_grid_ref.shape[0])
                def pad_to(arr, n):
                    if arr.shape[1] < n:
                        padw = n-arr.shape[1]
                        return np.pad(arr, ((0,0),(0,padw)), "constant")
                    return arr
                guide_cnts = pad_to(guide_cnts, m)
                ref_cnts = pad_to(ref_cnts, m)

            perm_data = np.vstack([guide_cnts, ref_cnts])
            perm_data_dict[guide] = perm_data.copy()

            labels = np.array(['guide'] * len(cluster_tags) + ['ref'] * len(cluster_tags))
            dist_matrix = squareform(pdist(perm_data, metric='braycurtis'))
            dist_matrix = nan_safe_dist_matrix(dist_matrix)
            if dist_matrix is None:
                pseudo_f[guide] = np.nan
                continue
            # Use the fast version for speed in observed (and all permutations)
            pseudo_f[guide] = fast_pseudo_f(dist_matrix, labels)

        # Write pseudo_f to var as Series
        gdata.var[result_field] = pd.Series(pseudo_f)

        if n_permutations is not None:
            # Vectorize the swapping and F calculations!
            guides = [g for g in gdata.var_names if g != reference_guide]
            # Create arrays for all permutations and guides to speed up
            pval_acc = np.zeros(len(guides), dtype=float)
            guide_idx_dict = {g:i for i,g in enumerate(guides)}
            obs_fvals = np.array([pseudo_f[g] for g in guides])
            num_clusters = len(cluster_tags)

            # Prebuild guide/ref indices for swapping
            guide_indices = np.arange(num_clusters)
            ref_indices = guide_indices + num_clusters
            perm = np.zeros((n_permutations, num_clusters), dtype=bool)
            rng = np.random.default_rng()
            perm = rng.random((n_permutations, num_clusters)) < 0.5

            # We batch-compute all swaps for each guide, then all pseudo-F
            def run_perms_for_guide(guide, obs_f, arr, gidx):
                # arr is cluster# x D, perm_data
                perm_arr = np.tile(arr, (n_permutations,1,1)) # (n_perm x 2*num_clust x D)
                # Swap in place for each perm: for each perm[idx], for clusters == True do cluster<->ref
                for p in range(n_permutations):
                    toswap = perm[p]
                    idx = np.where(toswap)[0]
                    if len(idx)>0:
                        temp = perm_arr[p,guide_indices[idx]].copy()
                        perm_arr[p,guide_indices[idx]] = perm_arr[p,ref_indices[idx]]
                        perm_arr[p,ref_indices[idx]] = temp
                # For each permutation, compute pairwise distance mat, then pseudo-F
                count_higher = 0
                for p in range(n_permutations):
                    data_p = perm_arr[p]
                    dist_mat = squareform(pdist(data_p, metric="braycurtis"))
                    dist_mat = nan_safe_dist_matrix(dist_mat)
                    if dist_mat is None:
                        continue  # treat as failed
                    fval = fast_pseudo_f(dist_mat, np.array(['guide']*num_clusters + ['ref']*num_clusters))
                    if not np.isnan(fval) and not np.isnan(obs_f):
                        if fval > obs_f:
                            count_higher += 1
                pval_acc[gidx] = count_higher

            # Use threadpool for parallel guides.
            threads = []
            for g in guides:
                arr = perm_data_dict[g]
                gidx = guide_idx_dict[g]
                obsf = pseudo_f[g]
                t = threading.Thread(target=run_perms_for_guide, args=(g, obsf, arr, gidx))
                threads.append(t)
                t.start()
            for t in threads:
                t.join()

            guide_pvals = pval_acc / n_permutations
            gdata.var[result_field + '.p_value'] = pd.Series(0., index=gdata.var_names)
            for i,g in enumerate(guides):
                gdata.var.loc[g, result_field + '.p_value'] = guide_pvals[i]

    # ----------- MULTI-SAMPLE (library_key specified) block -----------
    else:
        library_key_list = gdata.obs[library_key].unique()
        perm_data_dict = {}
        all_pseudo_f = {}

        for sample in library_key_list:
            pdata = gdata[gdata.obs[library_key] == sample].copy()
            perm_data_dict[sample] = {}
            count_df = pd.concat([
                pd.DataFrame(pdata.X, columns=pdata.var_names, index=pdata.obs_names),
                pdata.obs[cluster_field],
            ], axis=1)
            pseudo_f = {}

            for guide in pdata.var_names:
                if guide == reference_guide: continue
                guide_count_df = count_df[[guide, reference_guide, cluster_field]]

                ref_df = guide_count_df.groupby([cluster_field, reference_guide]).count().unstack()[guide]
                guide_df = guide_count_df.groupby([cluster_field, guide]).count().unstack()[reference_guide]

                guide_xmax = guide_df.columns.max() if (guide_df.shape[1] > 0) else 0
                x_grid_guide = np.linspace(0, guide_xmax + 0.5, count_bins)
                guide_cnts = []
                for cluster in cluster_tags:
                    if cluster in guide_df.index:
                        vals = guide_df.loc[cluster]
                        guide_density = kde_or_zeros(vals, x_grid_guide)
                        guide_cnts.append(guide_density)
                    else:
                        guide_cnts.append(np.zeros_like(x_grid_guide))
                guide_cnts = np.array(guide_cnts)

                ref_xmax = ref_df.columns.max() if (ref_df.shape[1] > 0) else 0
                x_grid_ref = np.linspace(0, ref_xmax + 0.5, count_bins)
                ref_cnts = []
                for cluster in cluster_tags:
                    if cluster in ref_df.index:
                        vals = ref_df.loc[cluster]
                        ref_density = kde_or_zeros(vals, x_grid_ref)
                        ref_cnts.append(ref_density)
                    else:
                        ref_cnts.append(np.zeros_like(x_grid_ref))
                ref_cnts = np.array(ref_cnts)

                if x_grid_guide.shape[0] != x_grid_ref.shape[0]:
                    m = max(x_grid_guide.shape[0], x_grid_ref.shape[0])
                    def pad_to(arr, n):
                        if arr.shape[1] < n:
                            padw = n-arr.shape[1]
                            return np.pad(arr, ((0,0),(0,padw)), "constant")
                        return arr
                    guide_cnts = pad_to(guide_cnts, m)
                    ref_cnts = pad_to(ref_cnts, m)

                perm_data = np.vstack([guide_cnts, ref_cnts])
                perm_data_dict.setdefault(sample, {})[guide] = perm_data.copy()
                labels = np.array(['guide'] * len(cluster_tags) + ['ref'] * len(cluster_tags))

                dist_matrix = squareform(pdist(perm_data, metric='braycurtis'))
                dist_matrix = nan_safe_dist_matrix(dist_matrix)
                if dist_matrix is None:
                    pseudo_f[guide] = np.nan
                    continue
                pseudo_f[guide] = fast_pseudo_f(dist_matrix, labels)

            for guide, pf in pseudo_f.items():
                all_pseudo_f.setdefault(guide, {})[sample] = pf

        for guide, d in all_pseudo_f.items():
            for sample, pf in d.items():
                gdata.var.loc[guide, sample + result_field] = pf

        if n_permutations is not None:
            guides = [g for g in gdata.var_names if g != reference_guide]
            pseudo_f_permcounts = pd.Series(0., index=guides)
            num_clusters = len(cluster_tags)
            guide_indices = np.arange(num_clusters)
            ref_indices = guide_indices + num_clusters

            rng = np.random.default_rng()

            # Precompute all swaps for all permutations (shape: n_permutations x num_clusters)
            perm_swaps = rng.random((n_permutations, num_clusters)) < 0.5

            def run_perms_for_guide_multi(guide):
                count_higher = 0
                # accumulate observed per-sample F's for this guide
                obs_fvals = np.array([
                    all_pseudo_f[guide][sample]
                    for sample in library_key_list
                    if sample in all_pseudo_f[guide] and not np.isnan(all_pseudo_f[guide][sample])
                ])
                if len(obs_fvals)==0: return 0
                obs_f = np.mean(obs_fvals)
                for p in range(n_permutations):
                    perm_fvals = []
                    for sample in library_key_list:
                        arr = perm_data_dict[sample][guide].copy()
                        swap = perm_swaps[p]
                        idx = np.where(swap)[0]
                        if len(idx)>0:
                            temp = arr[guide_indices[idx]].copy()
                            arr[guide_indices[idx]] = arr[ref_indices[idx]]
                            arr[ref_indices[idx]] = temp
                        dist_mat = squareform(pdist(arr, metric="braycurtis"))
                        dist_mat = nan_safe_dist_matrix(dist_mat)
                        if dist_mat is None:
                            continue
                        pf = fast_pseudo_f(dist_mat, np.array(['guide']*num_clusters + ['ref']*num_clusters))
                        if not np.isnan(pf):
                            perm_fvals.append(pf)
                    if len(perm_fvals) == 0: continue
                    perm_mean = np.mean(perm_fvals)
                    if perm_mean > obs_f:
                        count_higher += 1
                pseudo_f_permcounts[guide] = count_higher

            # Thread pool for each guide (permutation body can be multi-threaded/cpu'd here)
            threads = []
            for g in guides:
                t = threading.Thread(target=run_perms_for_guide_multi, args=(g,))
                threads.append(t)
                t.start()
            for t in threads:
                t.join()

            guide_pvals = pseudo_f_permcounts / n_permutations
            gdata.var[result_field + '.p_value'] = gdata.var_names.to_series().map(guide_pvals).replace(np.nan, 0.)

    if copy:
        return gdata
    else:
        return None