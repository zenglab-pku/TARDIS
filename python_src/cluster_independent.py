from optparse import Option
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata

from scipy.spatial import distance_matrix
from scipy.stats import wasserstein_distance
from statsmodels.nonparametric.kernel_density import KDEMultivariate
from scipy.stats import entropy

from multiprocessing import Pool
from typing import Optional, Iterable

def rank_by_proportion_within_threhold(
    adata: anndata,
    thresholds: Optional[Iterable | int] = np.linspace(1, 501, 21),
    guide_list: Optional[Iterable | None] = None,
    control_guide: Optional[str] = 'sgNon-targeting',
    method: Optional[str] = 'bin',
    kernel_func: Optional[str] = 'log',
    kernel_scoring: Optional[Iterable | None] = None,
    return_fig: Optional[bool] = False
) -> None:
    proportion = {}
    if type(thresholds) is int:
        thresholds = np.linspace(1, thresholds, 21)
    if guide_list == None:
        guide_list = adata.var_names.tolist()
    if method == 'bin':
        for threshold in thresholds:
            proportion[threshold] = []
            for guide in guide_list:
                if guide == control_guide: continue
                guide_data = adata[adata[:, guide].X > 0]
                ntc_data = adata[adata[:, control_guide].X > 0]
                in_distance = ((distance_matrix(guide_data.obsm['spatial'], ntc_data.obsm['spatial']) < threshold).sum(axis=1) > 0).sum()
                proportion[threshold].append(in_distance / guide_data.shape[0])
    elif method == 'count':
        for threshold in thresholds:
            proportion[threshold] = []
            for guide in guide_list:
                if guide == control_guide: continue
                guide_data = adata[adata[:, guide].X > 0]
                ntc_data = adata[adata[:, control_guide].X > 0]
                in_distance = guide_data[(distance_matrix(guide_data.obsm['spatial'], ntc_data.obsm['spatial']) < threshold).sum(axis=1) > 0, guide].X.sum()
                proportion[threshold].append(in_distance / guide_data.shape[0])
    else:
        print('Error, method must be one of "bin" or "count"!')

    d_df = pd.DataFrame(proportion).melt()
    d_df.index = adata[:, (adata.var_names != control_guide) & np.isin(adata.var_names, guide_list)].var_names.tolist() * len(thresholds)
    d_df.reset_index(inplace=True)

    if kernel_scoring != None:
        sc = d_df.apply(lambda x: np.multiply(kernel_scoring, 1 - x), axis=1).sum(axis=1).sort_values()
    else:
        if kernel_func == 'log': 
            ln_scoring = np.array([np.log2(x) for x in np.linspace(1, np.exp2(10), len(thresholds))])
            sc = d_df.apply(lambda x: np.multiply(ln_scoring, 1 - x), axis=1).sum(axis=1).sort_values()
        elif kernel_func == 'linear':
            linear_scoring = np.linspace(0, 10, len(thresholds))
            sc = d_df.apply(lambda x: np.multiply(linear_scoring, 1 - x), axis=1).sum(axis=1).sort_values()
        elif kernel_func == 'exp':
            exp_scoring = np.array([np.power(2, x) for x in np.linspace(1e-32, np.log2(10), len(thresholds))])
            sc = d_df.apply(lambda x: np.multiply(exp_scoring, 1 - x), axis=1).sum(axis=1).sort_values()

    df = pd.DataFrame(sc, columns=['dist'])
    df['guide'] = df.index.str.split('_').str.get(0).tolist()
    df.reset_index(inplace=True)
    _, ax = plt.subplots(1, 1, figsize=(len(guide_list) // 3 * 10, 6))
    sns.barplot(df, y='ln_dist', x='guide', ax=ax, palette='rainbow', orient='v', alpha=0.6)
    plt.xticks(rotation=90)
    
    plt.gca().set_yscale('log')
    plt.xlabel('Gene')
    plt.ylabel('Distance score')
    plt.show()

    if return_fig: return ax
    else: return

def rank_by_kernel_estimated_distance(
    adata: anndata,
    control_guide: Optional[str] = 'sgNon-targeting',
    guide_list: Optional[Iterable] = None,
    n_permutation: Optional[int] = 50,
    n_process: Optional[int] = 8,
    return_fig: Optional[bool] = False,
    return_dataframe: Optional[bool] = True,
    sort_by_replicate: Optional[str | None] = '_'
) -> None:
    
    def compute_kde_statsmodels(data):
        kde = KDEMultivariate(data, var_type='ccu')
        return kde.pdf(data)

    def kde_worker(args):
        guide, control, spatial_coords, data, n_permutation = args
        if guide == control: return
        gene_a_expression = data[:, guide].X.flatten()
        gene_b_expression = data[:, control].X.flatten()
        data_a = np.vstack([spatial_coords[:, 0], spatial_coords[:, 1], gene_a_expression]).T
        data_b = np.vstack([spatial_coords[:, 0], spatial_coords[:, 1], gene_b_expression]).T
        Za = compute_kde_statsmodels(data_a)
        Zb = compute_kde_statsmodels(data_b)
        actual_distance = wasserstein_distance(Za, Zb)
        combined_expression = np.concatenate([gene_a_expression, gene_b_expression])
        permuted_distances = []
        np.random.seed(42)
        for _ in range(n_permutation):
            np.random.shuffle(combined_expression)
            permuted_a = combined_expression[:len(gene_a_expression)]
            permuted_b = combined_expression[len(gene_a_expression):]
            permuted_data_a = np.vstack([spatial_coords[:, 0], spatial_coords[:, 1], permuted_a]).T
            permuted_data_b = np.vstack([spatial_coords[:, 0], spatial_coords[:, 1], permuted_b]).T
            permuted_Za = compute_kde_statsmodels(permuted_data_a)
            permuted_Zb = compute_kde_statsmodels(permuted_data_b)
            permuted_distance = wasserstein_distance(permuted_Za, permuted_Zb)
            permuted_distances.append(permuted_distance)
        p_value = np.sum(np.array(permuted_distances) >= actual_distance) / len(permuted_distances)
        return guide, actual_distance, p_value
    
    if guide_list == None:
        guide_list = adata.var_names.tolist()
    if return_fig and return_dataframe:
        print('Error, can only return one of fig or dataframe!')
        return
    filtered_data = adata[adata[:, guide_list].X.sum(axis=1) > 0]
    spatial_coords = filtered_data.obsm['spatial']
    with Pool(processes=n_process) as pool:
        results = pool.map(kde_worker, [(guide, control_guide, spatial_coords, filtered_data, n_permutation) for guide in guide_list])
        d_df = pd.DataFrame(filter(None, results), columns=['guide', 'wd', 'p_value']).sort_values(by='wd', ascending=True)

    _, ax = plt.subplots(1, 1, figsize=(len(guide_list) // 4, 8))

    if sort_by_replicate is not None:
        d_df['replicate'] = d_df['guide'].str.split(sort_by_replicate).str.get(0)
        positions =np.arange(len(d_df)) * 2
        width = 0.5

        ax.bar(positions - width/2, np.array(d_df.iloc[:, -2]).flatten(), width=width, label='1')
        ax.bar(positions + width/2, np.array(d_df.iloc[:, -1]).flatten(), width=width, label='2')

        for idx, name in enumerate(d_df.index):
            ax.hlines(y=max([d_df.iloc[idx]['-log(wd)_1'], d_df.iloc[idx]['-log(wd)_2']]) + 0.02, xmin=2*idx-0.5, xmax=2*idx+0.5, colors='black')
            ax.text(x=2*idx, y=max([d_df.iloc[idx]['-log(wd)_1'], d_df.iloc[idx]['-log(wd)_2']]) + 0.03, s='*', ha='center')

        ax.set_xticks(positions)
        ax.set_xticklabels(d_df.index, rotation=90)
        ax.set_xlabel('Guide(Ordered by min in each pair)')
        ax.set_ylabel('-log(WD)')
        ax.set_ylim(np.min(d_df.iloc[:, -2:]) - 0.1, np.max(d_df.iloc[:, -2:]) + 0.1)
        plt.show()
    else:
        sns.barplot(d_df, x='guide', y='wd', hue='guide', legend=False, palette='rainbow')

        for idx, name in enumerate(d_df.index):
            if (d_df.loc[name, 'wd'] > 0) and (d_df.loc[name, 'wd'] > 0):
                ax.hlines(y=d_df.iloc[idx]['wd'] + 0.02, xmin=idx-0.25, xmax=idx+0.25, colors='black')
                ax.text(x=idx, y=d_df.iloc[idx]['wd'] + 0.03, s=('*' if d_df.p_value < 0.05 else '.'), ha='center')

        ax.set_xticklabels(d_df.index, rotation=90)
        ax.set_xlabel('Guide')
        ax.set_ylabel('WD')
        ax.set_ylim(np.min(d_df.loc[d_df.wd > 0, 'wd']) - 0.1, np.max(d_df.wd) + 0.1)
        plt.show()

    if return_dataframe: return d_df
    if return_fig: return ax
    return

def rank_by_relative_entropy(
    adata: anndata,
    reference_guide: Optional[str] = 'sum',
    control_guide: Optional[str] = 'sgNon-targeting',
    result_field: Optional[str] = 'KL distance',
    guide_list: Optional[Iterable] = None,
    n_top: Optional[int] = 50
) -> None:
    entro = {}
    if reference_guide == 'sum':
        qk_guide = adata.X.sum(axis=1).flatten() / adata.X.sum()
    elif reference_guide == 'ntc':
        qk_guide = adata[:, control_guide].X.flatten() / adata[:, control_guide].X.sum()
    else:
        print('Error, reference_guide must be one of "sum" or "ntc"!')
        return
    if guide_list != None and control_guide not in guide_list:
        guide_list = adata.var_names.tolist() + [control_guide]
    elif guide_list == None:
        guide_list = adata.var_names.tolist()

    for guide in guide_list:
        pk_guide = adata[:, guide].X.flatten() / adata[:, guide].X.sum()
        entro[guide] = entropy(pk_guide, qk_guide)
    
    df = pd.DataFrame(entro, index=[result_field]).T.sort_values(by=result_field, ascending=False).head(n_top).reset_index()

    adata.uns[result_field] = df
    return