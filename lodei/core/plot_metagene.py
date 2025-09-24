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

mpl.rcParams['font.size'] = 8
mpl.rcParams['grid.color'] = '#b0b0b0'
mpl.rcParams['grid.linestyle'] = (0, (2, 3))
mpl.rcParams['grid.linewidth'] = 0.8
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['legend.markerscale'] = 1
mpl.rcParams['legend.fontsize'] = 8

lwd = 2.5

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

def plot_metagene_hist(rel_list, labels, bins, output, title="", mean_signal=False, wEI_list=None, ylim=None):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for idx, (rel, label) in enumerate(zip(rel_list, labels)):
        bin_values = [0] * (len(bins) - 1)
        bin_counts = [0] * (len(bins) - 1)
        # If mean_signal: wEI_list must correspond to rel!
        wEI = wEI_list[idx] if mean_signal else None
        for j, v in enumerate(rel):
            for i in range(len(bins) - 1):
                if (v >= bins[i] and v < bins[i + 1]) or (i == len(bins) - 2 and v == bins[i + 1]):
                    bin_counts[i] += 1
                    if mean_signal and wEI is not None:
                        bin_values[i] += wEI[j]
                    break
        if mean_signal and wEI is not None:
            # Calculate sum of wEI per bin
            y = bin_values + [bin_values[-1]]
            ax.step(bins, y, where='post', label=label, linewidth=lwd)
        else:
            y = bin_counts + [bin_counts[-1]]
            ax.step(bins, y, where='post', label=label, linewidth=lwd)
    ax.set_xlabel("Relative gene position (5' → 3')")
    ax.set_ylabel("Summed wEI signal" if mean_signal else "Number of significant windows")
    ax.set_title(title)
    ax.legend(loc="best")
    ax.grid()
    if ylim is not None:
        ax.set_ylim(ylim)
    fig.savefig(output, facecolor='#FFFFFF', bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

def make_plot(args):
    gff_path = args["gff"]
    txt_files = args["windows"]
    bins = int(args.get("num_bins"))
    ylim = args.get("ylim", None)
    mean_signal = args.get("mean_signal", False)
    output = args["output"]

    gff_df = pd.read_csv(
        gff_path, sep='\t', comment='#', header=None,
        names=['seqid','source','type','start','end','score','strand','phase','attributes']
    )
    print(f"Loaded GFF file '{gff_path}' with {len(gff_df)} rows")
    gff_df = filter_unique_genes(gff_df)

    rel_list = []
    labels = []
    wEI_list = []
    for txt_file in txt_files:
        windows_df = pd.read_csv(txt_file, sep="\t")
        print(f"Loaded windows file '{txt_file}' with {len(windows_df)} rows")

        merged = merge_windows_with_gff(windows_df, gff_df)
        n_lost = len(windows_df) - len(merged)
        if n_lost > 0:
            print(f"Note: {n_lost} significant windows could not be matched because their gene name was duplicated in the GFF and thus removed from the analysis.")
        
        rel = compute_relative_positions(merged)
        rel_list.append(rel)
        labels.append(os.path.basename(txt_file))
        if mean_signal:
            # wEI must have the same order as rel!
            wEI_list.append(merged["wEI"].values)
    bin_edges = np.linspace(0, 1, bins + 1)
    plot_metagene_hist(rel_list, labels, bin_edges, output, mean_signal=mean_signal, wEI_list=wEI_list if mean_signal else None, ylim=ylim)