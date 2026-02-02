# lodei - local differential editing index calculation between two sets of RNA-seq samples.
# Copyright (C) 2024 Phillipp Torkler, Tilman Heise, Jessica Stelzl

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

mpl.rcParams['font.size'] = 7
mpl.rcParams['grid.color'] = '#b0b0b0'
mpl.rcParams['grid.linestyle'] = (0, (2, 3))
mpl.rcParams['grid.linewidth'] = 0.8
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['legend.markerscale'] = 1
mpl.rcParams['legend.fontsize'] = 6

lwd = 1.75


def extract_gene_name(attributes):
    for item in attributes.split(";"):
        if item.strip().startswith("gene_name="):
            return item.strip().split("=", 1)[1]
    return None


def filter_unique_genes(gff_df):
    if "gene_name" not in gff_df.columns:
        gff_df["gene_name"] = gff_df["attributes"].apply(extract_gene_name)
    # Search for duplicates by (seqid, gene_name)
    duplicated = gff_df.duplicated(subset=["seqid", "gene_name"], keep=False)
    n_total = len(gff_df)
    n_duplicates = duplicated.sum()
    n_unique = n_total - n_duplicates
    print(f"The GFF file has {n_duplicates} duplicated (seqid, gene_name) pairs.")
    print(f"For the analysis, only the {n_unique} unique (seqid, gene_name) entries are used.")
    gff_df = gff_df[~duplicated].copy()
    return gff_df


def merge_windows_with_gff(windows_df, gff_df):
    """
    This function merges a DataFrame of significant windows (`windows_df`) with gene annotation data (`gff_df`),
    ensuring only unique gene names are used for the analysis. If duplicated gene names are found in the GFF,
    they are removed.
    """
    merged = windows_df.merge(
        gff_df,
        left_on=["chrom", "name"], right_on=["seqid", "gene_name"],
        how="inner", suffixes=("", "_gene")
    )
    return merged


def compute_relative_positions(merged):
    """
    The function determines the relative location of the midpoint of each window in the corresponding gene, taking strand orientation into account.
    For genes on the '+' strand, the position is measured from the gene start; for genes on the '-' strand, it is measured from the gene end.
    """
    merged["mid"] = ((merged["wstart"] + merged["wend"]) / 2)
    merged["length"] = (merged["end"] - merged["start"] + 1)
    rel = ((merged["end"] - merged["mid"]).where(merged["strand"].eq("-"),
                                                 merged["mid"] - merged["start"]) / merged["length"]).clip(0, 1)
    return rel


def plot_metagene_hist(ax, rel_list, labels,
                       bins, wEI_list=None, norm_global=False):
    ys = []
    for idx, (rel, label) in enumerate(zip(rel_list, labels)):
        bin_values = [0] * (len(bins) - 1)
        wEI = wEI_list[idx]
        for j, v in enumerate(rel):
            for i in range(len(bins) - 1):
                if (v >= bins[i] and v < bins[i + 1]) or (i == len(bins) - 2 and v == bins[i + 1]):
                    bin_values[i] += wEI[j]
                    break
            # Calculate sum of wEI per bin
        ys.append(bin_values + [bin_values[-1]]) # + due to the step() function
        #if wEI is not None:
        #    ys.append(bin_values + [bin_values[-1]])
        #else:
        #    ys.append(np.zeros(len(bin_values)))
    # for global rescaling
    # print(ys)
    # print(f"len(ys[0]): {len(ys[0])}")
    # print(f"len(bins): {len(bins)}")
    # print(f"bins: {bins}")
    #x = []
    #for i in range(len(bins)-1):
    #    x.append(bins[i] + ( (bins[i+1] - bins[i])/2 ))
    max_y = 0
    for y in ys:
        ym = np.max(np.abs(y))
        if ym > max_y:
            max_y = ym
    for y, label in zip(ys, labels):
        if norm_global:
            y = np.array(y) / max_y
        else:
            y = np.array(y) / np.max(np.abs(y))
            # print(y)
        #ax.step(x, y, where='post', label=label, linewidth=lwd)
        ax.step(bins, y, where='post', label=label, linewidth=lwd)
    ax.set_xlabel("Relative position (5' → 3')")
    ax.set_ylabel("Relative signal")
    ax.set_ylim([-1.1, 1.1])
    x = np.round(np.arange(0, 1.1, .1), 1)
    ax.set_xticks(x, x)
    ax.legend(loc="best")
    ax.grid()


def make_plot(args):
    bins = int(args.get("num_bins"))

    gff_df = pd.read_csv(
        args["gff"], sep='\t', comment='#', header=None,
        names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    )
    gff_df = filter_unique_genes(gff_df)

    if args["labels"] is None:
        labels = [f"File {i}" for i in range(len(args["windows"]))]
    else:
        labels = args["labels"]

    rel_list = []
    wEI_list = []
    for txt_file in args["windows"]:
        windows_df = pd.read_csv(txt_file, sep="\t")
        merged = merge_windows_with_gff(windows_df, gff_df)

        rel = compute_relative_positions(merged)
        rel_list.append(rel)
        wEI_list.append(merged["wEI"].values)
    bin_edges = np.linspace(0, 1, bins + 1)

    fig, ax = plt.subplots()
    fig.set_size_inches(args["width"], args["height"])
    plot_metagene_hist(ax, rel_list, labels,
                       bin_edges, wEI_list, args["global"])
    fig.savefig(args["output"], facecolor='#FFFFFF', bbox_inches='tight', pad_inches=0.1)


def make_plot_heat(args):
    bins = int(args.get("num_bins"))

    gff_df = pd.read_csv(
        args["gff"], sep='\t', comment='#', header=None,
        names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    )
    gff_df = filter_unique_genes(gff_df)

    if args["labels"] is None:
        labels = [f"File {i}" for i in range(len(args["windows"]))]
    else:
        labels = args["labels"]

    rel_list = []
    wEI_list = []
    for txt_file in args["windows"]:
        windows_df = pd.read_csv(txt_file, sep="\t")
        merged = merge_windows_with_gff(windows_df, gff_df)
        # print(merged)
        rel = compute_relative_positions(merged)
        rel_list.append(rel)
        wEI_list.append(merged["wEI"].values)
        # print(rel)
    bin_edges = np.linspace(0, 1, bins + 1)

    fig = plt.figure()
    fig.set_size_inches(args["width"], args["height"])
    plot_metagene_heat(fig, rel_list, labels,
                       bin_edges, wEI_list, args["global"])
    fig.savefig(args["output"], facecolor='#FFFFFF', bbox_inches='tight', pad_inches=0.1)


def plot_metagene_heat(fig, rel_list, labels,
                       bins, wEI_list=None, norm_global=False):
    ys = []
    # print(f"bins: {bins}")
    for idx, (rel, label) in enumerate(zip(rel_list, labels)):
        bin_values = [0] * (len(bins) - 1)
        wEI = wEI_list[idx]
        for j, v in enumerate(rel):
            for i in range(len(bins) - 1):
                if (v >= bins[i] and v < bins[i + 1]):  # or (i == len(bins) - 2 and v == bins[i + 1]):
                    # if wEI is not None:
                    #    bin_values[i] += wEI[j]
                    bin_values[i] += wEI[j]
                    break
        # print(f"bin_values: {bin_values}")
        ys.append(bin_values)
        #if wEI is not None:
        #    ys.append(bin_values + [bin_values[-1]])
        #else:
        #    ys.append(np.zeros(len(bin_values)))
    # print(ys)
    max_y = 0
    for y in ys:
        ym = np.max(np.abs(y))
        if ym > max_y:
            max_y = ym
    
    mat = np.zeros((len(ys), len(ys[0])))
    for i, y in enumerate(ys):
        y = np.array(y)
        #print(y / np.max(np.abs(y)))
        if norm_global:
            mat[i, :] = y / max_y
        else:
            mat[i, :] = y / np.max(np.abs(y))

    ax_mat = fig.add_axes([0.25, 0.15, 0.65, 0.8])
    ax_mat.imshow(mat, vmin=-1, vmax=1, cmap="bwr", aspect="auto")
    ax_mat.set_yticks(np.arange(mat.shape[0]), labels)

    x = np.round(np.arange(0, 1.1, .1), 1)
    # print(f"shape: {mat.shape}")
    xrel = x * (mat.shape[1])
    xrel = np.array([int(xp * mat.shape[1]) - 0.5 for xp in x])
    ax_mat.set_xticks(xrel, x)

    ax_mat.set_xlabel("Relative position (5' → 3')")

    cax = fig.add_axes([.93, .15, 0.015, 0.25])

    cmap = plt.get_cmap('bwr')
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)

    cb = mpl.colorbar.ColorbarBase(
        cax,
        cmap=cmap,
        norm=norm,
        orientation='vertical',
        ticks=[-1, 0, 1]
    )
    cb.set_ticklabels(['-1', ' 0', ' 1'])
    cb.set_label('Intensity')
