"""
Create a PCA plot.
"""

__author__ = "Janek Sendrowski"
__contact__ = "j.sendrowski18@gmail.com"
__date__ = "2023-05-27"

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas_plink import read_plink
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

try:
    testing = False
    bed = snakemake.input[0]
    sample_names = snakemake.params.samples
    name_col = snakemake.params.get('name_col', 'name')
    add_names = snakemake.params.get('add_names', False)
    marker_size = snakemake.params.get('marker_size', 20)
    label_col = snakemake.params.get('label_col', None)
    label_dict = snakemake.params.get('label_dict', None)
    legend_title = snakemake.params.get('legend_title', None)
    legend_outside = snakemake.params.get('legend_outside', True)
    legend = snakemake.params.get('add_legend', True)
    legend_n_cols = snakemake.params.get('legend_n_cols', 3)
    legend_size = snakemake.params.get('legend_size', 6)
    subsample_size = snakemake.params.get('subsample_size', 0)
    seed = snakemake.params.get('seed', 0)
    cbar = snakemake.params.get('cbar', True)
    cmap = snakemake.params.get('cmap', None)
    out = snakemake.output[0]
except NameError:
    # testing
    testing = True
    bed = 'results/plink/Pongo_abelii/Pongo_abelii.bed'
    df = pd.read_csv(f"resources/metadata/Pongo_individuals.csv",sep="\t")
    sample_names = df[df["BAM_FOLDER"].str.contains("Pongo_abelii")].GVCF_ID.tolist()
    name_col = 'sample'
    add_names = False
    marker_size = 20
    label_col = 'sample'
    label_dict = {}
    legend_title = None
    legend = True
    legend_size = 6
    legend_n_cols = 3
    cbar = False
    legend_outside = True
    subsample_size = 10000
    seed = 0
    cmap = None
    out = "scratch/pca.png"

samples = pd.DataFrame({name_col: sample_names})


def load_data(
        bed: str,
        samples: pd.DataFrame,
        subsample_size: int,
        seed: int,
        name_col: str = 'name'
):
    """
    Load the data from the bed file and the samples file.

    :param bed: Bed file with the genotypes.
    :param samples: Dataframe with the samples.
    :param subsample_size: Subsample the SNPs to this size.
    :param seed: The seed for the random number generator. Use ``None`` for no seed.
    :param name_col: The column in the samples file that contains the names of the samples.
    :raises ValueError: If the samples in the samples file are not a subset of the samples in the bed file.
    """
    # load genotypes
    (bim, fam, genotypes) = read_plink(bed)

    # subsample SNPs if desired
    if subsample_size and genotypes.shape[0] > subsample_size:
        indices = np.random.default_rng(seed=seed).choice(range(genotypes.shape[0]), subsample_size, replace=False)
        genotypes = genotypes[indices]

    fam["sid"] = fam.fid
    fam.loc[fam.fid != fam.iid, "sid"] += "_" + fam.iid[fam.fid != fam.iid]

    # subset samples if desired
    samples = samples[samples[name_col].isin(fam.sid.values)].copy()

    # order samples by fam.sid
    order = {sample: i for i, sample in enumerate(fam.sid.values)}
    samples["ordering"] = samples[name_col].map(order)
    samples.sort_values(by="ordering", inplace=True)

    # create a map of sample name to index in fam.sid
    fam_order = {name: i for i, name in enumerate(fam.sid.values)}

    # get indices of samples in fam.sid order
    samples_indices = samples[name_col].map(fam_order).values

    # subset
    genotypes = genotypes[:, samples_indices]

    # Silence FutureWarning: The `numpy.may_share_memory` function is not implemented by Dask array.
    # You may want to use the da.map_blocks function or something similar to silence this warning.
    # Your code may stop working in a future release
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # impute missing genotype values
        genotypes = SimpleImputer(missing_values=np.nan, strategy="mean").fit_transform(genotypes)

    # transpose
    return np.array(genotypes).T, samples


def prepare_plot(
        x: np.ndarray,
        y: np.ndarray,
        label_dict: dict,
        samples: pd.DataFrame,
        label_col: str,
        marker_size: int = 20,
        cmap: str = None,
        legend: bool = True,
        legend_title: str = None,
        legend_outside: bool = True,
        legend_n_cols: int = 1,
        legend_size: int = 6,
        cbar: bool = None,
        name_col: str = 'name',
        add_names: bool = False,
) -> plt.Axes:
    """
    Prepare a PCA plot.

    :param x: The x values.
    :param y: The y values.
    :param label_dict: A dictionary that maps the labels to new labels.
    :param samples: Dataframe with the samples.
    :param label_col: The column in the samples file that contains the labels.
    :param marker_size: The size of the markers.
    :param cmap: The colormap.
    :param legend: Whether to add a legend.
    :param legend_title: The title of the legend.
    :param legend_outside: Whether to put the legend outside the plot.
    :param legend_n_cols: The number of columns in the legend.
    :param legend_size: The size of the legend.
    :param cbar: Whether to add a color bar.
    :param name_col: The column in the samples file that contains the names of the samples.
    :param add_names: Whether to add sample names to the plot.
    :return:
    """
    # recode labels if desired and replace missing values with placeholder
    if label_dict:
        samples['label_dict'] = [label_dict[x] if x in label_dict else '-' for x in samples[label_col]]

        label_col = 'label_dict'
    else:
        # replace missing values
        if label_col:
            samples.loc[samples[label_col].isna(), label_col] = '-'

    # plot the 2 components
    ax = sns.scatterplot(x=x, y=y, hue=samples[label_col] if label_col else None, s=marker_size, palette=cmap)

    # whether to include a legend
    if legend:
        # set legend title
        ax.get_legend().set_title(legend_title)

        # whether to place the legend outside the graph
        if legend_outside:
            plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=legend_n_cols, prop={'size': legend_size})
    else:
        if ax.get_legend():
            ax.get_legend().remove()

    # whether to show a color bar
    if cbar:
        norm = plt.Normalize(samples[label_col].min(), samples[label_col].max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        ax.figure.colorbar(sm)

    # add sample names if specified
    if add_names:
        eps = 0.15

        for i in range(samples.shape[0]):
            plt.text(x=x[i] + eps, y=y[i] + eps, s=samples[name_col][i], size=4, alpha=0.5)

    plt.gca().axis('square')

    plt.tight_layout(pad=5)

    return plt.gca()


genotypes, samples_used = load_data(
    bed=bed,
    samples=samples,
    subsample_size=subsample_size,
    seed=seed,
    name_col=name_col
)

pca = PCA(n_components=2)
pc = pca.fit_transform(genotypes)

prepare_plot(
    x=pc[:, 0],
    y=pc[:, 1],
    label_dict=label_dict,
    samples=samples_used,
    label_col=label_col,
    marker_size=marker_size,
    cmap=cmap,
    legend=legend,
    legend_title=legend_title,
    legend_outside=legend_outside,
    legend_n_cols=legend_n_cols,
    legend_size=legend_size,
    cbar=cbar,
    name_col=name_col,
    add_names=add_names
)

# add axis labels
v1 = pca.explained_variance_ratio_[0] * 100
v2 = pca.explained_variance_ratio_[1] * 100
plt.xlabel(f'{v1:.2f}%')
plt.ylabel(f'{v2:.2f}%')

plt.savefig(out, dpi=200, bbox_inches='tight', pad_inches=0.1)

if testing:
    plt.show()
