import anndata
import cellcharter as cc
import scvi
import numpy as np
import squidpy as sq

from typing import Optional

def cluster_cellcharter(
    adata: anndata,
    batch_key: Optional[str | None],
    n_clusters: Optional[tuple | int] = 'auto',
    spatial_key: Optional[str] = 'spatial',
    cluster_field: Optional[str] = 'cluster',
    layer: Optional[str] = None,
    random_seed: Optional[int] = 42,
    n_nodes_hidden_layers: Optional[int] = 32,
    dim_latent_layers: Optional[int] = 10,
    n_hidden_layers: Optional[int] = 5,
    gene_likelihood_model: Optional[str] = 'poisson',
    latent_distribution: Optional[str] = 'normal',
    save_model: Optional[str | bool] = None,
    use_model: Optional[str | bool] = None,
    inplace: Optional[bool] = True
) -> None:
    
    if len(adata.var_names[adata.var_names.str.startswith('sg')]) != 0:
        print('Warning, clustering data with guide reads. Please check data.')
    scvi.settings.seed = random_seed
    if batch_key != None and batch_key not in adata.obs.columns:
        assert False, "Error, batch field not in adata obs, please check batch key settings."
    if not inplace:
        cluster_data = adata.copy()
    else:
        cluster_data = adata

    if use_model == False or use_model == None:
        scvi.model.SCVI.setup_anndata(cluster_data, batch_key=batch_key, layer=layer)
        print('Training scVI model, this may take a while...')
        model = scvi.model.SCVI(cluster_data, n_hidden=n_nodes_hidden_layers, n_latent=dim_latent_layers,
                                n_layers=n_hidden_layers, gene_likelihood=gene_likelihood_model,
                                latent_distribution=latent_distribution)
        model.train(early_stopping=True, enable_progress_bar=True)

        if save_model == True:
            print('Saving scVI model...')
            model.save("scvi.model", save_anndata=True, overwrite=True)
        elif save_model != False and save_model != None:
            print('Saving scVI model...')
            model.save(save_model, save_anndata=True, overwrite=True)   
    elif use_model == True:
        print('Loading scVI model...')
        model = scvi.model.SCVI.load("scvi.model")
    else:
        print('Loading scVI model...')
        model = scvi.model.SCVI.load(use_model)

    cluster_data.obsm["X_scVI"] = model.get_latent_representation(cluster_data).astype(np.float32)
    cluster_data.obs[batch_key] = cluster_data.obs[batch_key].astype("category")

    print('Calculating spatial neighbors...')
    sq.gr.spatial_neighbors(cluster_data, library_key=batch_key, coord_type="generic", delaunay=True, spatial_key=spatial_key)
    cc.gr.remove_long_links(cluster_data)
    cc.gr.aggregate_neighbors(cluster_data, n_layers=3, use_rep="X_scVI", out_key="X_cellcharter", sample_key=batch_key)

    print(f'Clustering with AutoK {n_clusters}...')
    if isinstance(n_clusters, tuple):
        autok = cc.tl.ClusterAutoK(n_clusters=(n_clusters[0], n_clusters[1]), max_runs=10, model_params=dict(random_state = random_seed))
    else:
        autok = cc.tl.Cluster(n_clusters=n_clusters)
    autok.fit(cluster_data, use_rep="X_cellcharter")
    cluster_data.obs[cluster_field] = autok.predict(cluster_data, use_rep="X_cellcharter")
    print(f'Almost done...')
    sq.pl.spatial_scatter(cluster_data, color=cluster_field, library_key=batch_key, size=50,
                          img=None, spatial_key=spatial_key, cmap="Set3", figsize=(15,15), shape=None)
    if not inplace:
        return cluster_data
    else: return