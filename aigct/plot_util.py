import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def barchart(df: pd.DataFrame, x: str, y: str, x_label: str = '',
             y_label: str = '', ci: pd.DataFrame = None,
             filename: str = None, hatch: str = None,
             edgecolor=None,
             palette=None, errorbar=None, estimator: str = 'mean',
             title: str = '', title_fontsize: int = 20,
             y_label_fontsize: int = 15,
             x_label_fontsize: int = 15,
             label_size: int = 15, xtick_rotation: int = None,
             xtick_rotation_mode: str = None,
             file_dpi: int = 300, bbox_inches: str = 'tight'):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax = sns.barplot(df, x=x, y=y, palette=palette, hatch=hatch,
                     edgecolor=edgecolor, errorbar=errorbar,
                     estimator=estimator, hue=x)
    ax.tick_params('both', labelsize=label_size)
    ax.set_ylabel(y_label, fontsize=y_label_fontsize)
    ax.set_xlabel(x_label, fontsize=x_label_fontsize)
    if xtick_rotation:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=xtick_rotation,
                           ha="right", rotation_mode=xtick_rotation_mode)
    if ci is not None and len(ci) != 0:
        plt.errorbar(x=df[x], y=df[y], yerr=ci["CI"], fmt="none", c="k")
    plt.title(title, fontsize=title_fontsize)
    if filename is None:
        plt.show()
        return
    plt.savefig(filename, dpi=file_dpi, bbox_inches=bbox_inches)
    plt.savefig(filename.replace('.png', '.svg'), bbox_inches=bbox_inches)


def create_feature_palette(values: list) -> dict:
    if len(values) > 255:
        colors = sns.color_palette('husl', len(values))
    else:
        # Get a list of all available colors in Seaborn's color palette
        all_colors = sns.color_palette("husl", 256)

        # Randomly select colors from the list
        random_colors = np.random.choice(range(len(all_colors)),
                                         len(values), replace=False)
        colors = [all_colors[i] for i in random_colors]
    return dict(zip(values, colors))


def get_colors(num_colors: int) -> list:
    if num_colors > 255:
        colors = sns.color_palette('husl', num_colors)
    else:
        # Get a list of all available colors in Seaborn's color palette
        all_colors = sns.color_palette("husl", 256)

        # Randomly select colors from the list
        random_colors = np.random.choice(range(len(all_colors)),
                                         num_colors, replace=False)
        colors = [all_colors[i] for i in random_colors]
    return colors


def heatmap(df: pd.DataFrame, y_column: str = None, x_column: str = None,
            value_column: str = None, x_label: str = '', y_label: str = '',
            title: str = '', cmap: str = 'viridis', annot: bool = True,
            fmt: str = '.2f', linewidths: float = 0.5,
            title_fontsize: int = 20, label_fontsize: int = 15,
            tick_fontsize: int = 12, annot_fontsize: int = 10,
            figsize: tuple = (10, 8), filename: str = None,
            vmin: float = None, vmax: float = None,
            cbar: bool = True, cbar_kws: dict = None,
            square: bool = False, mask: pd.DataFrame = None,
            file_dpi: int = 300, bbox_inches: str = 'tight'):
    """
    Plot a heatmap using seaborn.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to plot as heatmap
    y_column : str, optional
        Column to use for y-axis categories. If provided along with x_column
        and value_column, the function will pivot the data before plotting.
    x_column : str, optional
        Column to use for x-axis categories. If provided along with y_column
        and value_column, the function will pivot the data before plotting.
    value_column : str, optional
        Column containing values to plot. If provided along with x_column and
        y_column, the function will pivot the data before plotting.
    x_label : str, optional
        Label for x-axis
    y_label : str, optional
        Label for y-axis
    title : str, optional
        Title of the plot
    cmap : str, optional
        Colormap to use
    annot : bool, optional
        Whether to annotate cells with values
    fmt : str, optional
        Format string for annotations
    linewidths : float, optional
        Width of lines between cells
    title_fontsize : int, optional
        Font size for title
    label_fontsize : int, optional
        Font size for axis labels
    tick_fontsize : int, optional
        Font size for tick labels
    annot_fontsize : int, optional
        Font size for annotations
    figsize : tuple, optional
        Figure size (width, height) in inches
    filename : str, optional
        If provided, save figure to this file
    vmin : float, optional
        Minimum value for colormap scaling
    vmax : float, optional
        Maximum value for colormap scaling
    cbar : bool, optional
        Whether to draw a colorbar
    cbar_kws : dict, optional
        Additional arguments for colorbar
    square : bool, optional
        Whether to make cells square-shaped
    mask : pd.DataFrame, optional
        Boolean DataFrame of same shape as df, True values will not be plotted
    file_dpi : int, optional
        DPI for saved figure
    bbox_inches : str, optional
        Bounding box in inches for saved figure
        
    Returns
    -------
    matplotlib.axes.Axes
        The matplotlib axes containing the plot
    """
    if (x_column is not None and y_column is not None and
            value_column is not None):
        df = df.pivot(index=y_column, columns=x_column, values=value_column)

    fig, ax = plt.subplots(figsize=figsize)

    # Set default colorbar keywords if not provided
    if cbar_kws is None:
        cbar_kws = {'shrink': 0.8}

    # Create the heatmap
    hm = sns.heatmap(df, ax=ax, cmap=cmap, annot=annot, fmt=fmt,
                     linewidths=linewidths, vmin=vmin, vmax=vmax,
                     cbar=cbar, cbar_kws=cbar_kws, square=square,
                     mask=mask, annot_kws={'size': annot_fontsize})

    # Set labels and title
    ax.set_xlabel(x_label, fontsize=label_fontsize)
    ax.set_ylabel(y_label, fontsize=label_fontsize)
    ax.set_title(title, fontsize=title_fontsize)

    # Set tick parameters
    ax.tick_params(axis='both', labelsize=tick_fontsize)

    # Rotate x-axis labels if needed
    plt.xticks(rotation=45, ha='right')

    # Tight layout
    plt.tight_layout()

    # Save figure if filename is provided
    if filename:
        plt.savefig(filename, dpi=file_dpi, bbox_inches=bbox_inches)
        # Also save as SVG for vector graphics
        plt.savefig(filename.replace('.png', '.svg'), bbox_inches=bbox_inches)

    return ax
