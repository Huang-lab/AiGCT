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
