import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from typing import Optional

def plot_ranking_bar(
    adata: anndata,
    field_name: Optional[str] = None,
    figsize: Optional[tuple] = (5, 3),
    n_top: Optional[int] = 50
) -> None:
    pdf = adata.uns[field_name]
    pdf = pdf.sort_values(by=field_name, ascending=True).head(n_top).reset_index()
    _, ax = plt.subplots(1, 1, figsize=figsize)
    sns.barplot(pdf, x='index', y=field_name, palette='RdBu_r', ax=ax)
    plt.xticks(rotation=90, fontsize=8)
    plt.gca().set_yscale('log')
    plt.title(f'{field_name} Ranking')
    plt.xlabel('Guide')
    plt.ylabel(field_name)

    sns.despine()
    plt.show()

def plot_ranking_scatter(
    adata: anndata,
    field_name: Optional[str] = None,
    figsize: Optional[tuple] = (5, 3),
    n_top: Optional[int] = 50
) -> None:
    pdf = adata.uns[field_name]
    pdf = pdf.sort_values(by=field_name, ascending=True).head(n_top).reset_index()
    _, ax = plt.subplots(1, 1, figsize=figsize)
    sns.scatterplot(pdf, x='index', y=field_name, ax=ax, color='black', s=10)
    plt.xticks(rotation=90, fontsize=8)
    plt.title(f'{field_name} Ranking')
    plt.xlabel('Guide')
    plt.ylabel(field_name)

    plt.gca().invert_yaxis()
    sns.despine()
    plt.show()

def plot_ranking_hist(
    adata: anndata,
    field_name: Optional[str] = None,
    figsize: Optional[tuple] = (5, 3),
    n_top: Optional[int] = 50
) -> None:
    pdf = adata.uns[field_name]
    pdf['z_score'] = (pdf[field_name] - pdf[field_name].mean()) / pdf[field_name].std()

    _, ax = plt.subplots(1, 1, figsize=figsize)
    sns.kdeplot(pdf, x='z_score', ax=ax, fill=True, bw_adjust=0.5)
    plt.title(f'{field_name} Distribution')
    plt.xlabel('Z-score')
    plt.ylabel('Density')
    sns.despine()
    plt.show()