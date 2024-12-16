import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import permanova, DistanceMatrix
from scipy.stats import chi2_contingency

from typing import Optional, Iterable

def rank_by_chi_square(
    adata: anndata,
    cluster_field: str,
    control_guide: Optional[str] = 'sgNon-targeting',
    guide_list: Optional[Iterable] = None,
    result_field: Optional[str] = 'Chi2 p-value'
) -> None:
    if cluster_field not in adata.obs.columns:
        print("Error, cluster field not in anndata obs, please perform clustering first!")
        return
    
    if guide_list != None and control_guide not in guide_list:
        guide_list = adata.var_names.tolist() + [control_guide]
    elif guide_list == None:
        guide_list = adata.var_names.tolist()
    if isinstance(adata.X, np.ndarray):
        c_df = pd.DataFrame(adata[:, guide_list].X, columns=guide_list)
    else:
        c_df = pd.DataFrame(adata[:, guide_list].X.toarray(), columns=guide_list)

    c_df[cluster_field] = adata.obs[cluster_field].tolist()
    c_df = c_df.groupby([cluster_field]).sum()

    c_df = c_df.div(c_df.sum(axis=1), axis=0) * c_df.loc[:, control_guide].sum()

    chi_dict = {}
    for guide in guide_list:
        if c_df.loc[:, guide].sum() == 0:
            chi_dict[guide] = 1
            continue
        if guide == control_guide: continue
        chi_dict[guide] = chi2_contingency(c_df.loc[:, [control_guide, guide]].T).pvalue
    
    pdf = pd.DataFrame(chi_dict, index=[result_field]).T.sort_values(by=result_field)

    adata.uns[result_field] = pdf
    return

def rank_by_cluster_permanova(
    adata: anndata,
    cluster_field: str,
    control_guide: Optional[str] = 'sgNon-targeting',
    n_permutation: Optional[int] = 999,
    guide_list: Optional[Iterable] = None,
    result_field: Optional[str] = 'PERMANOVA p-value'
) -> None:
    if cluster_field not in adata.obs.columns:
        print("Error, cluster field not in anndata obs, please perform clustering first!")
        return
    
    if guide_list != None and control_guide not in guide_list:
        guide_list = adata.var_names.tolist() + [control_guide]
    elif guide_list == None:
        guide_list = adata.var_names.tolist()

    if isinstance(adata.X, np.ndarray):
        c_df = pd.concat([pd.DataFrame(adata[:, guide_list].X, columns=guide_list, index=adata.obs_names), adata.obs[cluster_field]], axis=1)
    else:
        c_df = pd.concat([pd.DataFrame(adata[:, guide_list].X.toarray(), columns=guide_list, index=adata.obs_names), adata.obs[cluster_field]], axis=1)

    n_clusters = len(adata.obs[cluster_field].unique())
    p_value = {}
    for guide in guide_list:

        if guide == control_guide: continue
        g_df = c_df[[guide, control_guide, cluster_field]]

        df = pd.DataFrame(index=range(n_clusters))
        guide_cnts = np.array(pd.concat([df, g_df.groupby([cluster_field, guide]).count().unstack()[control_guide]], axis=0).fillna(0))
        ntc_cnts = np.array(pd.concat([df, g_df.groupby([cluster_field, control_guide]).count().unstack()[guide]], axis=0).fillna(0))
        print(df, g_df.groupby([cluster_field, control_guide]).count().unstack()[guide],
              df, g_df.groupby([cluster_field, guide]).count().unstack()[control_guide])
        print(pd.concat([df, g_df.groupby([cluster_field, control_guide]).count().unstack()[guide]], axis=0).fillna(0),
              pd.concat([df, g_df.groupby([cluster_field, guide]).count().unstack()[control_guide]], axis=0).fillna(0))
        data = np.vstack([guide_cnts, ntc_cnts])
        sample_ids = np.array([[name + str(c) for c in g_df[cluster_field].unique()] for name in ['guide_', 'control_']]).flatten()

        dist_matrix = squareform(pdist(data, metric='euclidean'))
        metadata = pd.DataFrame({
            'group': ['guide'] * n_clusters + ['control'] * n_clusters
        }, index=sample_ids)
        dm = DistanceMatrix(dist_matrix, ids=sample_ids)
        results = permanova(dm, metadata, column='group', permutations=n_permutation)
        p_value[guide] = results['p-value']
    
    p_df = pd.DataFrame(p_value.values(), index=p_value.keys()).sort_values(by=0, ascending=False)

    adata.uns[result_field] = p_df
    return
