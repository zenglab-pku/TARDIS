from typing import Optional, Iterable, Union
import anndata as ad
import numpy as np
import pandas as pd
from tqdm import tqdm
import statsmodels.api as sm

import concurrent.futures

# --- Module-level helper functions for multiprocess pickling (must be at top-level for pickling) ---

def _perm_kl_worker(
    guide_idx,
    idx_ntc,
    n_permutations,
    guide_matrix,
    epx,
    kl_div
):
    """Permutation test for one guide, returns p."""
    rng = np.random.default_rng()  # use default random generator for independence
    perm_greater_count = 0

    for _ in range(n_permutations):
        # shuffle cells
        idx = np.arange(guide_matrix.shape[0])
        rng.shuffle(idx)
        permuted = guide_matrix[idx, :]
        perm_guide_vec = permuted[:, guide_idx]
        perm_ntc_vec = permuted[:, idx_ntc]
        perm_pk_sum = perm_guide_vec.sum()
        perm_ntc_sum = perm_ntc_vec.sum()
        if perm_pk_sum == 0 or perm_ntc_sum == 0:
            continue
        perm_pk_guide = perm_guide_vec / perm_pk_sum + epx
        perm_qk_guide = perm_ntc_vec / perm_ntc_sum + epx
        perm_kl_div = _kl_div(
            perm_pk_guide / perm_pk_guide.sum(),
            perm_qk_guide / perm_qk_guide.sum()
        )
        if perm_kl_div > kl_div:
            perm_greater_count += 1
    return perm_greater_count / n_permutations

def _perm_kl_arg_wrapper(args):
    # args = (guide_idx, kl_div, idx_ntc, n_permutations, guide_matrix, epx)
    guide_idx, kl_div, idx_ntc, n_permutations, guide_matrix, epx = args
    return _perm_kl_worker(
        guide_idx,
        idx_ntc,
        n_permutations,
        guide_matrix,
        epx,
        kl_div
    )

def _perm_multi_worker(args):
    seed, library_key_list, guide_data, all_guides, reference_guide, guide_list, library_key, result_field = args
    rng = np.random.default_rng(seed)
    perm_entropy = {}
    for library_key_value in library_key_list:
        perm_entropy[library_key_value] = {}
        obs_mask = (guide_data.obs[library_key] == library_key_value)
        present_guides = [g for g in all_guides if g in guide_data.var_names]
        pdata = guide_data[obs_mask, present_guides].copy()
        if reference_guide not in pdata.var_names:
            continue
        idx_ntc = list(pdata.var_names).index(reference_guide)
        guide_matrix = pdata.X.copy()
        qk_sum = pdata[:, reference_guide].X.sum()
        if qk_sum == 0:
            continue
        epx = 1 / qk_sum
        for guide in guide_list:
            if guide not in pdata.var_names:
                perm_entropy[library_key_value][guide] = np.nan
                continue
            idx_guide = list(pdata.var_names).index(guide)
            idx_arr = np.arange(guide_matrix.shape[0])
            rng.shuffle(idx_arr)
            permuted = guide_matrix[idx_arr, :]
            perm_guide_vec = permuted[:, idx_guide]
            perm_ntc_vec = permuted[:, idx_ntc]
            perm_pk_sum = perm_guide_vec.sum()
            perm_ntc_sum = perm_ntc_vec.sum()
            if perm_pk_sum == 0 or perm_ntc_sum == 0:
                perm_entropy[library_key_value][guide] = np.nan
                continue
            perm_guide = perm_guide_vec / perm_pk_sum + epx
            perm_ntc = perm_ntc_vec / perm_ntc_sum + epx
            perm_entropy[library_key_value][guide] = _kl_div(
                perm_guide / perm_guide.sum(), perm_ntc / perm_ntc.sum())
    import pandas as pd
    perm_df = pd.DataFrame(perm_entropy)
    # Ranking
    for col in perm_df.columns:
        non_nan_mask = perm_df[col].notna()
        perm_df.loc[non_nan_mask, col + '.rank'] = perm_df.loc[non_nan_mask, col].rank(method='min', ascending=False)
    rank_cols = perm_df.columns[perm_df.columns.str.endswith('.rank')]
    if len(rank_cols) > 0:
        perm_df['mean_rank'] = perm_df[rank_cols].mean(axis=1, skipna=True)
        perm_df['final_rank'] = perm_df['mean_rank'].rank(method='min', ascending=True)
    return perm_df['final_rank'] if 'final_rank' in perm_df else pd.Series(dtype=float)

def kl_divergence(
    guide_data: ad.AnnData,
    reference_guide: str = 'sgNon-targeting',
    result_field: str = 'kl_div',
    guide_list: Iterable = None,
    library_key: Optional[str] = None,
    n_permutations: Optional[int] = 1000,
    show_progress: bool = True,
    copy: bool = False,
    n_jobs: Optional[int] = None,
):
    """
    Compute the Kullback-Leibler (KL) divergence between each guide and a reference guide, 
    based on their overall count distributions, with an optional permutation test for significance.

    Parameters
    ----------
    guide_data : AnnData
        AnnData object containing .obs (cell metadata), .var (guide-level information), and .X (expression/count matrix).
    reference_guide : str, default "sgNon-targeting"
        The guide to use as the reference (typically a negative control).
    result_field : str, default "kl_div"
        Name of the output field to write KL divergences into guide_data.var.
    guide_list : Iterable or None, default None
        List of guides to calculate KL divergence for; if None, uses all guides except reference_guide.
    library_key : str or None, default None
        If set, compute KL divergence within each subgroup defined by this field; if None, pooled analysis.
    n_permutations : int or None, default 1000
        Number of label permutations to estimate p-values; set to None to skip permutation test.
    show_progress : bool, default True
        Whether to display a progress bar for permutation calculation.
    copy : bool, default False
        If True, return a copy of AnnData with results written; otherwise update in place.
    n_jobs : int or None, default None
        Number of parallel worker processes to use for permutation test (if enabled).

    Returns
    -------
    If copy is True, returns an AnnData object with KL divergences written to .var[result_field].
    Otherwise, returns None and updates guide_data in place.

    Notes
    -----
    This function computes the KL divergence for each guide versus the reference guide, handling zero or missing values robustly.
    If permutation testing is enabled, p-values are estimated using fast label shuffling in parallel.
    """

    if copy:
        guide_data = guide_data.copy()

    # If guide_list is not provided, use all variable names except the reference_guide
    if guide_list is None:
        guide_list = [x for x in guide_data.var_names if x != reference_guide]

    all_guides = list(set(guide_list) | {reference_guide})

    if not isinstance(guide_data.X, np.ndarray):
        guide_data.X = guide_data.X.toarray()

    if reference_guide not in guide_data.var_names:
        raise KeyError(f"Reference guide '{reference_guide}' not found in AnnData.var_names.")

    if library_key is None:
        # Keep all guides for slicing, so reference_guide is present
        pdata = guide_data[:, all_guides].copy()

        if result_field not in pdata.var.columns:
            pdata.var[result_field] = np.nan
        if n_permutations is not None and f"{result_field}.p_value" not in pdata.var.columns:
            pdata.var[f"{result_field}.p_value"] = 0.0

        qk_guide_vec = pdata[:, reference_guide].X.flatten()
        qk_sum = qk_guide_vec.sum()
        if qk_sum == 0:
            raise ValueError(f"Reference guide '{reference_guide}' sums to zero.")

        epx = 1 / qk_sum
        qk_guide = qk_guide_vec / qk_sum + epx

        # Calculate KL divergence for all guides
        for guide in guide_list:
            pk_guide_vec = pdata[:, guide].X.flatten()
            pk_sum = pk_guide_vec.sum()
            if pk_sum == 0:
                pdata.var.at[guide, result_field] = np.nan
                continue
            pk_guide = pk_guide_vec / pk_sum + epx
            kl_div = _kl_div(pk_guide / pk_guide.sum(), qk_guide / qk_guide.sum())
            pdata.var.at[guide, result_field] = kl_div

        # Parallelized permutation test for p-value
        if n_permutations is not None:
            for guide in guide_list:
                pdata.var.at[guide, f"{result_field}.p_value"] = 0.0
            idx_ntc = list(pdata.var_names).index(reference_guide)
            all_kl_divs = [pdata.var.at[guide, result_field] for guide in guide_list]
            guide_indices = [list(pdata.var_names).index(g) for g in guide_list]
            guide_matrix = pdata[:, all_guides].X

            # Only run permutation for valid (finite, non-nan) KL values!
            valid_guides = []
            valid_indices = []
            valid_kl_divs = []
            for i, g in enumerate(guide_list):
                kl = all_kl_divs[i]
                if np.isfinite(kl):
                    valid_guides.append(g)
                    valid_indices.append(guide_indices[i])
                    valid_kl_divs.append(kl)

            input_tuples = [
                (guide_idx, kl_div, idx_ntc, n_permutations, guide_matrix, epx)
                for guide_idx, kl_div in zip(valid_indices, valid_kl_divs)
            ]

            if n_jobs is None:
                import os
                n_jobs = min(8, (os.cpu_count() or 1))  # Cap at 8 by default

            # Multiprocessing: use process_map if available else ProcessPoolExecutor with tqdm progress
            if show_progress:
                try:
                    from tqdm.contrib.concurrent import process_map as _procmap
                except ImportError:
                    _procmap = None
            else:
                _procmap = None

            if _procmap is not None:
                perm_pvals = _procmap(
                    _perm_kl_arg_wrapper,
                    input_tuples,
                    max_workers=n_jobs,
                    desc="KL divergence permutations"
                )
            else:
                from tqdm import tqdm
                import concurrent.futures
                with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
                    perm_pvals = list(tqdm(
                        executor.map(_perm_kl_arg_wrapper, input_tuples),
                        total=len(valid_guides),
                        desc="KL divergence permutations",
                        disable=not show_progress
                    ))

            for g, p in zip(valid_guides, perm_pvals):
                pdata.var.at[g, f"{result_field}.p_value"] = p

        # Write back to guide_data.var (in case guide_data has >guide_list subset in var):
        for guide in guide_list:
            guide_data.var.at[guide, result_field] = pdata.var.at[guide, result_field]
            if n_permutations is not None:
                guide_data.var.at[guide, f"{result_field}.p_value"] = pdata.var.at[guide, f"{result_field}.p_value"]
    else:
        # Library key is provided, iterate over its unique values
        library_key_list = pd.unique(guide_data.obs[library_key])
        for library_key_value in library_key_list:
            obs_mask = (guide_data.obs[library_key] == library_key_value)
            present_guides = [g for g in all_guides if g in guide_data.var_names]
            if reference_guide not in guide_data.var_names:
                raise KeyError(f"Reference guide '{reference_guide}' not found in AnnData.var_names.")
            pdata = guide_data[obs_mask, present_guides].copy()
            result_col = f"{library_key_value}{result_field}"
            if result_col not in pdata.var.columns:
                pdata.var[result_col] = np.nan
            if f"{result_col}.rank" not in pdata.var.columns:
                pdata.var[f"{result_col}.rank"] = np.nan

            qk_guide_vec = pdata[:, reference_guide].X.flatten()
            qk_sum = qk_guide_vec.sum()
            if qk_sum == 0:
                pdata.var.loc[:, result_col] = np.nan
                continue
            epx = 1 / qk_sum
            qk_guide = qk_guide_vec / qk_sum + epx

            for guide in guide_list:
                if guide not in pdata.var_names:
                    pdata.var.at[guide, result_col] = np.nan
                    continue
                pk_guide_vec = pdata[:, guide].X.flatten()
                pk_sum = pk_guide_vec.sum()
                if pk_sum == 0:
                    pdata.var.at[guide, result_col] = np.nan
                    continue
                pk_guide = pk_guide_vec / pk_sum + epx
                kl_div = _kl_div(pk_guide / pk_guide.sum(), qk_guide / qk_guide.sum())
                pdata.var.at[guide, result_col] = kl_div

            non_nan_mask = pdata.var[result_col].notna()
            ranked = pdata.var.loc[non_nan_mask, result_col].rank(method='min', ascending=False)
            pdata.var.loc[non_nan_mask, f"{result_col}.rank"] = ranked
            for guide in guide_list:
                if guide in pdata.var.index:
                    guide_data.var.at[guide, f"{result_col}.rank"] = pdata.var.at[guide, f"{result_col}.rank"]

        rank_cols = guide_data.var.columns[guide_data.var.columns.str.endswith(result_field + '.rank')]
        if len(rank_cols) > 0:
            guide_data.var[result_field + '.mean_rank'] = guide_data.var[rank_cols].mean(axis=1, skipna=True)
            guide_data.var[result_field + '.final_rank'] = guide_data.var[result_field + '.mean_rank'].rank(method='min', ascending=True)

        # Permutation on multi-library (parallel over permutations)
        if n_permutations is not None:
            guide_data.var[result_field + '.p_value'] = np.zeros(len(guide_data.var))
            if n_jobs is None:
                import os
                n_jobs = min(8, (os.cpu_count() or 1))
            # Precompute the reference values for all guides
            final_ranks = guide_data.var[result_field + '.final_rank'].copy()

            # Stable set of seeds and argument tuples, include all needed context in args for multiproc
            if show_progress:
                try:
                    from tqdm import tqdm
                except ImportError:
                    tqdm = None
            seeds_list = [np.random.SeedSequence().entropy + i for i in range(n_permutations)]
            args_list = [
                (
                    seed,
                    library_key_list,
                    guide_data,
                    all_guides,
                    reference_guide,
                    guide_list,
                    library_key,
                    result_field,
                )
                for seed in seeds_list
            ]

            import concurrent.futures
            perm_results = []
            with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
                if show_progress and tqdm is not None:
                    for result in tqdm(executor.map(_perm_multi_worker, args_list),
                                       total=n_permutations, desc="Permutation", disable=not show_progress):
                        perm_results.append(result)
                else:
                    perm_results = list(executor.map(_perm_multi_worker, args_list))

            # Aggregate permutation nulls
            null_final_ranks = {guide: [] for guide in guide_list}
            for final_rank_series in perm_results:
                for guide in guide_list:
                    if guide in final_rank_series and np.isfinite(final_rank_series[guide]):
                        null_final_ranks[guide].append(final_rank_series[guide])

            # Now, compute p-value as fraction of permuted final_ranks < actual final_rank
            for guide in guide_list:
                true_rank = final_ranks.get(guide, np.nan)
                if not np.isfinite(true_rank):
                    continue
                pval = np.mean([r < true_rank for r in null_final_ranks[guide]]) if len(null_final_ranks[guide]) > 0 else np.nan
                guide_data.var.at[guide, result_field + '.p_value'] = pval

    if copy:
        return guide_data
    else:
        return None

import concurrent.futures
from functools import partial

def wasserstein_distance(
    guide_data: ad.AnnData,
    reference_guide: str = 'sgNon-targeting',
    result_field: str = 'w_dist',
    guide_list: Iterable = None,
    spatial_key: str = 'spatial',
    library_key: Optional[str] = None,
    n_permutations: Optional[int] = 1000,
    spatial_interval: float = 200,
    bin_width: Optional[Union[str, int]] = 'silverman',
    show_progress: bool = True,
    copy: bool = False,
    n_jobs: int = 4,
):
    """
    Compute the Wasserstein distance (Earth Mover's Distance) between spatial expression distributions
    of specified guides and a reference guide, optionally with permutation-based significance testing.

    Parameters
    ----------
    guide_data : AnnData
        AnnData object containing .obs (cell metadata), .var (guide-level metadata), and .X (expression matrix).
    reference_guide : str, default "sgNon-targeting"
        The guide to use as the reference distribution for distance calculation.
    result_field : str, default "w_dist"
        Name of the field in .var to which calculated Wasserstein distances will be written.
    guide_list : Iterable, optional
        List of guides to analyze. If None, will use all guides except the reference_guide.
    spatial_key : str, default "spatial"
        Key in .obsm giving spatial coordinates for each cell/sample.
    library_key : str or None, default None
        Field in .obs representing sample/library for stratified analysis; if None, all pooled.
    n_permutations : int or None, default 1000
        Number of permutations for label-shuffling based significance estimation.
    spatial_interval : float, default 200
        Bin interval (distance in spatial units) for discretizing the spatial distribution.
    bin_width : str or int, default "silverman"
        Bin width for kernel density estimation; if "silverman", uses Silverman's rule of thumb.
    show_progress : bool, default True
        If True, displays a progress bar for permutation computations.
    copy : bool, default False
        If True, returns a new AnnData object with results; otherwise, updates in place.
    n_jobs : int, default 4
        Number of parallel worker processes to use for permutation testing.

    Returns
    -------
    If copy is True, returns an AnnData object with Wasserstein distances written to .var[result_field].
    Otherwise, returns None and updates guide_data in place.

    Notes
    -----
    The function supports parallel computation of permutation-based null distributions for more efficient p-value computation.
    Wasserstein distances are computed between the distribution of observed guide expression and that of the reference guide,
    using spatial or other binned representations. Bin width can be auto-selected using "silverman" or set manually.
    """
    if copy:
        guide_data = guide_data.copy()
    if guide_list is None:
        guide_list = [x for x in guide_data.var_names if x != reference_guide]

    if not isinstance(guide_data.X, np.ndarray):
        guide_data.X = guide_data.X.toarray()

    guide_data = guide_data[:, guide_list].copy()

    def compute_grid_bins(x_min, x_max, y_min, y_max, spatial_interval):
        x_bins = max(int(np.floor((x_max - x_min) / spatial_interval)), 2)
        y_bins = max(int(np.floor((y_max - y_min) / spatial_interval)), 2)
        x_grid = np.linspace(x_min, x_max, x_bins)
        y_grid = np.linspace(y_min, y_max, y_bins)
        return x_grid, y_grid

    def _kde_for_guide(args):
        pdata, spatial_key, x_grid, y_grid, guide, bin_width, reference_guide = args
        if bin_width == 'silverman':
            cnt_grid = np.arange(
                pdata[:, guide].X.min(),
                pdata[:, guide].X.max(),
                max(1e-10, 1.06 * np.std(pdata[:, guide].X) * (len(pdata[:, guide].X) ** (-1 / 5)))
            ).round(2)
        else:
            cnt_grid = np.arange(pdata[:, guide].X.min(), pdata[:, guide].X.max(), bin_width)
        return (guide, _calculate_kde(pdata, spatial_key, x_grid, y_grid, cnt_grid, pdata[:, guide].X.flatten()))

    def _perm_worker(args):
        # args: (pdata, guide, idx_guide, idx_ntc, bin_width, x_grid, y_grid, reference_guide, current_wdist)
        pdata, guide, idx_guide, idx_ntc, bin_width, x_grid, y_grid, reference_guide, current_wdist = args
        perm_guide_vec, perm_ntc_vec = _permute_guide_bins(pdata[:, guide_list + [reference_guide]].X, idx_guide, idx_ntc)
        if bin_width == 'silverman':
            ntc_grid = np.arange(
                perm_ntc_vec.min(),
                perm_ntc_vec.max(),
                max(1e-10, 1.06 * np.std(perm_ntc_vec) * (len(perm_ntc_vec) ** (-1 / 5)))
            ).round(2)
            guide_grid = np.arange(
                perm_guide_vec.min(),
                perm_guide_vec.max(),
                max(1e-10, 1.06 * np.std(perm_guide_vec) * (len(perm_guide_vec) ** (-1 / 5)))
            ).round(2)
        else:
            ntc_grid = np.arange(perm_ntc_vec.min(), perm_ntc_vec.max(), bin_width).round(2)
            guide_grid = np.arange(perm_guide_vec.min(), perm_guide_vec.max(), bin_width).round(2)
        kde_ntc = _calculate_kde(pdata, spatial_key, x_grid, y_grid, ntc_grid, perm_ntc_vec)
        kde_guide = _calculate_kde(pdata, spatial_key, x_grid, y_grid, guide_grid, perm_guide_vec)
        w_dist = wasserstein_distance(kde_guide, kde_ntc)
        return w_dist > current_wdist

    if library_key is None:
        pdata = guide_data.copy()
        if result_field not in pdata.var.columns:
            pdata.var[result_field] = np.zeros(len(guide_list))
        x_min, x_max = pdata.obsm[spatial_key][:, 0].min(), pdata.obsm[spatial_key][:, 0].max()
        y_min, y_max = pdata.obsm[spatial_key][:, 1].min(), pdata.obsm[spatial_key][:, 1].max()
        x_grid, y_grid = compute_grid_bins(x_min, x_max, y_min, y_max, spatial_interval)
        # 加速KDE的计算
        kde_args = [
            (pdata, spatial_key, x_grid, y_grid, guide, bin_width, reference_guide)
            for guide in guide_list + [reference_guide]
        ]
        kde_matrix_dict = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as executor:
            for res in executor.map(_kde_for_guide, kde_args):
                kde_matrix_dict[res[0]] = res[1]
        for guide in guide_list:
            pdata.var[result_field][guide] = wasserstein_distance(kde_matrix_dict[guide], kde_matrix_dict[reference_guide])
        if n_permutations is not None:
            pdata.var[result_field + '.p_value'] = np.zeros(len(guide_list))
            idx_ntc = pdata.var_names.get_loc(reference_guide)
            perm_tasks = []
            for guide in guide_list:
                idx_guide = pdata.var_names.get_loc(guide)
                # accumulate tasks for n_permutations for this guide
                perm_tasks += [
                    (pdata, guide, idx_guide, idx_ntc, bin_width, x_grid, y_grid, reference_guide, pdata.var[result_field][guide])
                    for _ in range(n_permutations)
                ]
            # distribute tasks to threads/processes
            num_guides = len(guide_list)
            perm_results = [0 for _ in range(num_guides)]
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as executor:
                res_iter = executor.map(_perm_worker, perm_tasks)
                if show_progress:
                    from tqdm import tqdm
                    res_iter = tqdm(res_iter, total=len(perm_tasks), desc="Permutations")
                counting = [0 for _ in range(num_guides)]
                for idx, is_better in zip(range(len(perm_tasks)), res_iter):
                    guide_idx = guide_list.index(perm_tasks[idx][1])
                    counting[guide_idx] += is_better
            for i, guide in enumerate(guide_list):
                pdata.var[result_field + '.p_value'][guide] = counting[i] / n_permutations
        guide_data.var[result_field + '.p_value'] = pdata.var[result_field + '.p_value']
        guide_data.var[result_field] = pdata.var[result_field]
    else:
        library_key_list = guide_data.obs[library_key].unique()
        for library_key in library_key_list:
            pdata = guide_data[guide_data.obs[library_key] == library_key].copy()
            if library_key + result_field not in pdata.var.columns:
                pdata.var[library_key + result_field] = np.zeros(len(guide_list))
            x_min, x_max = pdata.obsm[spatial_key][:, 0].min(), pdata.obsm[spatial_key][:, 0].max()
            y_min, y_max = pdata.obsm[spatial_key][:, 1].min(), pdata.obsm[spatial_key][:, 1].max()
            x_grid, y_grid = compute_grid_bins(x_min, x_max, y_min, y_max, spatial_interval)
            kde_args = [
                (pdata, spatial_key, x_grid, y_grid, guide, bin_width, reference_guide)
                for guide in guide_list + [reference_guide]
            ]
            kde_matrix_dict = {}
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as executor:
                for res in executor.map(_kde_for_guide, kde_args):
                    kde_matrix_dict[res[0]] = res[1]
            for guide in guide_list:
                pdata.var[library_key + result_field][guide] = wasserstein_distance(kde_matrix_dict[guide], kde_matrix_dict[reference_guide])
            non_nan_mask = pdata.var[library_key + result_field].notna()
            pdata.var.loc[non_nan_mask, library_key + result_field + '.rank'] = pdata.var[library_key + result_field].rank(method='min', ascending=False)
            pdata.var.loc[non_nan_mask, library_key + result_field + '.rank'] = np.nan
            guide_data.var[library_key + result_field + '.rank'] = pdata.var[library_key + result_field + '.rank']
        rank_cols = guide_data.var.columns[guide_data.var.columns.str.endswith(result_field + '.rank')]
        guide_data.var[result_field + '.mean_rank'] = guide_data.var[rank_cols].mean(axis=1, skipna=True)
        guide_data.var[result_field + '.final_rank'] = guide_data.var[result_field + '.mean_rank'].rank(method='min', ascending=True)
        if n_permutations is not None:
            guide_data.var[result_field + '.p_value'] = np.zeros(len(guide_list))
            perm_tasks = []
            for _ in range(n_permutations):
                for library_key in library_key_list:
                    pdata = guide_data[guide_data.obs[library_key] == library_key].copy()
                    idx_ntc = pdata.var_names.get_loc(reference_guide)
                    for guide in guide_list + [reference_guide]:
                        idx_guide = pdata.var_names.get_loc(guide)
                        perm_tasks.append((pdata, guide, idx_guide, idx_ntc, bin_width, x_grid, y_grid, reference_guide, None))
            perm_dist_records = [{} for _ in range(len(library_key_list))]
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as executor:
                from tqdm import tqdm
                res_iter = executor.map(_perm_worker, perm_tasks)
                if show_progress:
                    res_iter = tqdm(res_iter, total=len(perm_tasks), desc="Permutations")

    if copy:
        return guide_data
    else:
        return None

def _calculate_kde(pdata, spatial_key, x_grid, y_grid, ntc_grid, cnt_vec):
    X, Y, CNT = np.meshgrid(x_grid, y_grid, ntc_grid)
    kde = sm.nonparametric.KDEMultivariate(data=[
        pdata.obsm[spatial_key][:, 0],
        pdata.obsm[spatial_key][:, 1],
        cnt_vec
    ], var_type='cco',
    bw='normal_reference')
    grid_coords = np.column_stack([X.ravel(), Y.ravel(), CNT.ravel()])
    kde_values = kde.pdf(grid_coords).reshape(X.shape)
    kde_matrix = np.zeros((kde_values.shape[0], kde_values.shape[1]))
    for i in range(kde_values.shape[2]):
        kde_matrix += kde_values[:, :, i] * ntc_grid[i]
    return kde_matrix

def _kl_div(pk, qk):
    return np.sum(np.where(pk != 0, pk * np.log(pk / qk), 0))

def _permute_guide_bins(guide_matrix, idx_guide, idx_ntc, swap_rate=0.5):
    t_bins = (guide_matrix[:, idx_guide] > 0) | (guide_matrix[:, idx_ntc] > 0)
    guide_vec = guide_matrix[t_bins, idx_guide].copy()
    ntc_vec = guide_matrix[t_bins, idx_ntc].copy()
    swap_mask = np.random.randn(len(guide_vec)) < 1 - swap_rate * 2
    temp = guide_vec[swap_mask]
    guide_vec[swap_mask] = ntc_vec[swap_mask]
    ntc_vec[swap_mask] = temp
    ret_guide_vec = guide_matrix[:, idx_guide].copy()
    ret_ntc_vec = guide_matrix[:, idx_ntc].copy()
    ret_guide_vec[t_bins] = guide_vec
    ret_ntc_vec[t_bins] = ntc_vec
    return ret_guide_vec, ret_ntc_vec
