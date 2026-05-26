import matplotlib.pyplot as plt
import numpy as np

from statsmodels.nonparametric.kde import KDEUnivariate
import seaborn as sns
import matplotlib.colors as mcolors

def plot_top_kde(
    gdata,
    result_field: str,
    show_sgnt: bool = True,
    top_n: int = 2,
    violin: bool = True,
    sgnt_label: str = "sgNon-targeting",
    figsize: tuple = (5, 2.7),
    min_count: int = 0
):
    """
    Plot the KDE distribution of Aitchison distances, with options to highlight sgNon-targeting and the top N genes with the highest distances.

    Parameters
    ----------
    gdata : AnnData
        AnnData object with guide-level results in .var. Must contain the field specified by `result_field`.
    result_field : str
        Column name in gdata.var containing the distance metric (e.g., "aitchison_dist").
    show_sgnt : bool, default True
        Whether to highlight the sgNon-targeting guide in the plot.
    top_n : int, default 2
        Number of top guides (with highest distances) to highlight.
    violin : bool, default True
        Whether to draw 'violin' vertical lines for each data point in the lower panel.
    sgnt_label : str, default "sgNon-targeting"
        Name (index) of sgNon-targeting guide in gdata.var.
    figsize : tuple, default (5, 2.7)
        Figure size for matplotlib.
    min_count : int, default 0
        Minimum value in 'TotalCount' (if present in gdata.var) required for a guide to be included.

    Returns
    -------
    matplotlib.figure.Figure, matplotlib.axes.Axes
        The figure and axes containing the plot.
    """

    if (result_field not in gdata.var.columns) or (gdata.var[result_field].notna().sum() <= 0):
        raise ValueError(f"{result_field} not found in gdata.var or all values are NA/0")
        return

    var = gdata.var.copy()
    if 'TotalCount' in var.columns:
        var_sub = var[var['TotalCount'] > min_count]
    else:
        var_sub = var

    if len(var_sub) == 0 or result_field not in var_sub.columns:
        raise ValueError(f"{result_field} not found in valid gene list")
        return

    vals = var_sub[result_field].values
    kde = KDEUnivariate(vals)
    kde.fit()

    fig, ax = plt.subplots(2, 1, figsize=figsize, sharex=False, height_ratios=[3, 1], gridspec_kw={'hspace': 0})
    ax[0].set_ylim(0, max(kde.density)*1.15)

    ax[0].fill_between(kde.support, kde.density, alpha=1, color='gray')
    ax[0].plot(kde.support, kde.density, color='gray')

    lim = ax[0].get_xlim()
    ax[0].set_xlabel('')
    ax[0].set_xticks([])
    ax[0].set_ylabel('Density')
    sns.despine(ax=ax[0])

    if violin:
        for val in var_sub[result_field]:
            ax[1].vlines(x=val, ymin=0, ymax=1, color='gray', linestyles='-', alpha=0.3, linewidth=2)

    if show_sgnt and (sgnt_label in var_sub.index):
        sgnt_val = var_sub.loc[sgnt_label, result_field]
        color = 'red'
        yval = kde.density[np.abs(kde.support - sgnt_val).argmin()]
        ax[0].vlines(x=sgnt_val, ymin=0, ymax=yval,
                     color=color, linestyles='-', linewidth=3, alpha=1, zorder=10)
        ax[0].text(sgnt_val, yval, sgnt_label,
                   color=color, rotation=0, ha='left', va='bottom', fontsize=9, fontweight='bold', zorder=11)
        ax[1].vlines(x=sgnt_val, ymin=0, ymax=1, color=color, linestyles='-', linewidth=3, alpha=1, zorder=10)

    top_genes = (
        var_sub.loc[var_sub.index != sgnt_label, result_field]
        .sort_values(ascending=False)
        .head(top_n)
        .index
    )
    colors = plt.get_cmap("tab10")(np.linspace(0, 1, len(top_genes)))
    for i, gene in enumerate(top_genes):
        val = var_sub.loc[gene, result_field]
        color = colors[i % len(colors)]
        yval = kde.density[np.abs(kde.support - val).argmin()]
        ax[0].vlines(x=val, ymin=0, ymax=yval,
                     color=color, linestyles='--', alpha=1, linewidth=2)
        ax[0].text(val, yval, gene,
                   color=color, rotation=0, ha='left', va='bottom', fontsize=8, fontweight='bold')
        ax[1].vlines(x=val, ymin=0, ymax=1, color=color, linestyles='-', alpha=1, linewidth=3)

    ax[1].set_facecolor('whitesmoke')
    ax[1].set_ylim(0, 1)
    ax[1].set_xlim(lim)
    ax[1].set_xlabel(result_field.replace('_', ' ').capitalize())
    ax[1].set_ylabel('')
    ax[1].set_yticks([])

    plt.tight_layout()
    plt.show()

def plot_spatial_guides(
    gdata,
    scale_factor,
    image=None,
    figsize=(10, 10),
    palette='tab20b',
    s=3,
    edgecolor='none',
    legend=False,
    alpha=1.0,
    ax=None,
):
    """
    Plot spatial distribution of guides on a reference image.

    Parameters
    ----------
    gdata : AnnData object
        Gene/barcode matrix with .obsm['spatial'] coordinates and .var_names as guides.
    scale_factor : float
        Scale factor to apply to spatial coordinates, for matching image scale.
    image : PIL.Image or np.ndarray, optional
        Reference background (tissue) image to plot under spots.
    figsize : tuple, default (10, 10)
        Figure size for matplotlib.
    palette : str or sequence, default 'tab20b'
        Color palette for guides.
    s : float, default 3
        Point size.
    edgecolor : str, default 'none'
        Marker edge color.
    legend : bool, default False
        Whether to plot legend.
    alpha : float, default 1.0
        Marker transparency.
    ax : matplotlib.axes.Axes, optional
        Existing axis to plot on. If None, create new.

    Returns
    -------
    matplotlib.axes.Axes
        The axis with the plot.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    sparse_data = gdata[gdata.X.sum(axis=1) > 0]
    plot_df = pd.DataFrame({
        'x': sparse_data.obsm['spatial'][:, 2],
        'y': sparse_data.obsm['spatial'][:, 1],
        'guide': sparse_data.var_names[
            sparse_data.X.toarray().argmax(axis=1)
        ]
    })

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    if image is not None:
        ax.imshow(image)
    ax.axis('off')
    sns.scatterplot(
        x=plot_df.x * scale_factor,
        y=plot_df.y * scale_factor,
        hue=plot_df.guide,
        palette=palette,
        s=s,
        edgecolor=edgecolor,
        legend=legend,
        alpha=alpha,
        ax=ax,
    )
    ax.set_xlim(plot_df.x.min() * scale_factor, plot_df.x.max() * scale_factor)
    ax.set_ylim(plot_df.y.min() * scale_factor, plot_df.y.max() * scale_factor)
    return ax

def plot_ranking_scatter(gdata, NTC, result_field: str = 'dist'):
    """
    Plots ranked distances vs -log10(p) for all guides except NTC.

    Parameters
    ----------
    gdata : AnnData or compatible object with .var
        Annotated guide data containing result_field and result_field.p_value in .var.
    NTC : str
        Name of the negative control (NTC) guide to exclude from the plot.
    result_field : str
        Column name in gdata.var containing the distance metric (e.g., "aitchison_dist").
    """
    top_clones = [g for g in gdata.var_names if g != NTC]

    plot_data = gdata.var.dropna(subset=[result_field, f'{result_field}.p_value']).copy()
    plot_data = plot_data[plot_data.index.isin(top_clones)]
    plot_data[f'{result_field}.rank'] = plot_data[result_field].rank(ascending=False)
    plot_data['-log10(p)'] = -np.log10(plot_data[f'{result_field}.p_value'] + 1e-10)

    norm = mcolors.LogNorm(
        vmin=plot_data['-log10(p_val)'].replace(0, np.nan).min() if (plot_data['-log10(p_val)'] > 0).any() else 1e-3, 
        vmax=plot_data['-log10(p_val)'].max()
    )

    _, ax = plt.subplots(figsize=(5, 3))
    sns.scatterplot(
        plot_data,
        x=f'{result_field}',
        y=f'{result_field}.rank',
        hue='-log10(p_val)',
        palette='viridis',
        s=10,
        edgecolor="none",
        ax=ax,
        legend=False,
        hue_norm=norm,
    )
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='-log10(p_val)', pad=0.01)
    plt.gca().invert_yaxis()