# lodei - local differential editing index calculation between two sets of RNA-seq samples.
# Copyright (C) 2024 Phillipp Torkler, Tilman Heise

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
import os
import pysam
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.patches as patches
import lodei.core.eifunctions as eif


def plot_single_signals(ax, df, pstart=0, s=None, e=None, grid_col="#C0C0C0",
                        grid_lwd=.5, grid_style=(0, (1, 5)), lwd=2, nox=False):
    # if s is None:
    #     s = 0
    # if e is None:
    #     e = df["x"][df.shape[0] - 1]
    s = 0
    e = df["x"][df.shape[0] - 1]
    sel = (df["x"] >= s) & (df["x"] <= e)
    
    columns = list(df.columns)
    sel_g0 = [0]
    for i in range(len(columns)):
        if columns[i][0:2] == "g0":
            sel_g0.append(i)
    sel_g1 = [0]
    for i in range(len(columns)):
        if columns[i][0:2] == "g1":
            sel_g1.append(i)
    sel_g0 = np.array(sel_g0)
    sel_g1 = np.array(sel_g1)
    # print(f"sel_g0: {sel_g0}")
    # print(f"sel_g1: {sel_g1}")
    g0 = df.iloc[:, sel_g0]
    g1 = df.iloc[:, sel_g1]
    
    cmap = mpl.colormaps["Blues"]
    cvalue = np.linspace(0, 1, num=(g0.shape[1] - 1) + 3)
    cvalue = cvalue[1:-1]
    for i in range(1, g0.shape[1]):
        # ax.plot(g0["x"][sel], g0.iloc[sel,i], c="blue")
        ax.plot(g0["x"][sel], g0.iloc[:, i][sel], c=cmap(cvalue[i]),
                label=g0.columns[i], zorder=10, linewidth=lwd)
    
    cmap = mpl.colormaps["Oranges"]
    cvalue = np.linspace(0, 1, num=(g1.shape[1] - 1) + 3)
    cvalue = cvalue[1:-1]
    for i in range(1, g1.shape[1]):
        ax.plot(g1["x"][sel], g1.iloc[:, i][sel], c=cmap(cvalue[i]),
                label=g1.columns[i], zorder=20, linewidth=lwd)
    # ax.grid(zorder=2)
    
    ylim = ax.get_ylim()
    ax.set_ylim(ylim)
    ax.set_xlim([s, e])
    ax.set_ylabel(r"$e^{A \rightarrow G}_{s,w}$")
    gridposx = ax.get_xticks()
    # print(gridpos)
    for x in gridposx:
        ax.plot([x, x], [-ylim[1], 2 * ylim[1]], c=grid_col,
                linestyle=grid_style, zorder=2, linewidth=grid_lwd)
    gridposy = ax.get_yticks()
    for y in gridposy:
        ax.plot([gridposx[0], gridposx[-1]], [y, y], c=grid_col,
                linestyle=grid_style, zorder=2, linewidth=grid_lwd)
    if nox:
        ax.get_xaxis().set_visible(False)
    else:
        xticks = ax.get_xticks()
        xlabels = [pstart]
        for i in range(1, len(xticks)):
            xlabels.append("+" + str(int(xticks[i])))
        ax.set_xticks(ax.get_xticks(), xlabels)
        
    ax.legend(loc="upper left")


def plot_group_signals(ax, df, s=None, e=None, diff=True, c="#484848", grid_col="#C0C0C0", grid_lwd=.5, grid_style=(0, (1, 5))):
    if s is None:
        s = 0
    if e is None:
        e = df["x"][df.shape[0] - 1]
    sel = (df["x"] >= s) & (df["x"] <= e)
    
    cmap = plt.cm.get_cmap("Blues")
    ax.plot(df["x"][sel], df["av_g0"][sel], c=cmap(.75), label=r"$EI_{G0}=$mean(G0)", zorder=10)
    cmap = plt.cm.get_cmap("Oranges")
    ax.plot(df["x"][sel], df["av_g1"][sel], c=cmap(.75), label=r"$EI_{G1}=$mean(G1)", zorder=10)
    ax.legend(loc="upper left")
    # ax.grid()
    
    ylim = ax.get_ylim()
    ax.set_ylim(ylim)
    ax.set_xlim([s, e])
    ax.set_ylabel("Local EI Signal")
    gridpos = ax.get_xticks()
    # print(gridpos)
    for x in gridpos:
        ax.plot([x, x], [-ylim[1], 2 * ylim[1]], c=grid_col, linestyle=grid_style, zorder=2, linewidth=grid_lwd)


def plot_average_signals(ax, df, pstart=0, s=None, e=None, grid_col="#C0C0C0",
                         grid_lwd=.5, grid_style=(0, (1, 5)), lwd=2, nox=False):
    # if s is None:
    #     s = 0
    # if e is None:
    #     e = df["x"][df.shape[0] - 1]
    s = 0
    e = df["x"][df.shape[0] - 1]
    sel = (df["x"] >= s) & (df["x"] <= e)
    
    ax.plot(df["x"][sel], df["av_g0"][sel], c="blue", label=r"$z^{A \rightarrow G}_{S}$", zorder=10, linewidth=lwd)
    ax.plot(df["x"][sel], df["av_g1"][sel], c="orange", label=r"$z^{A \rightarrow G}_{S'}$", zorder=10, linewidth=lwd)
    
    ylim = ax.get_ylim()
    ax.set_ylim(ylim)
    ax.set_xlim([s, e])
    ax.set_ylabel(r"$z^{A \rightarrow G}$")
    gridposx = ax.get_xticks()
    # print(gridpos)
    for x in gridposx:
        ax.plot([x, x], [-ylim[1], 2 * ylim[1]], c=grid_col, linestyle=grid_style, zorder=2, linewidth=grid_lwd)
    gridposy = ax.get_yticks()
    for y in gridposy:
        ax.plot([gridposx[0], gridposx[-1]], [y, y], c=grid_col, linestyle=grid_style, zorder=2, linewidth=grid_lwd)

    if nox:
        ax.get_xaxis().set_visible(False)
    else:
        xticks = ax.get_xticks()
        xlabels = [pstart]
        for i in range(1, len(xticks)):
            xlabels.append("+" + str(int(xticks[i])))
        ax.set_xticks(ax.get_xticks(), xlabels)
    ax.legend(loc="upper left")


def plot_group_differences(ax, df, pstart=0, s=None, e=None, diff=True, c="#484848",
                           grid_col="#C0C0C0", grid_lwd=.5, grid_style=(0, (1, 5)), lwd=2):
    # if s is None:
    #     s = 0
    #  if e is None:
    #     e = df["x"][df.shape[0] - 1]
    s = 0
    e = df["x"][df.shape[0] - 1]
    sel = (df["x"] >= s) & (df["x"] <= e)

    y = df["av_diff"].copy()
    y[np.isnan(y)] = 0
    ax.plot(df["x"][sel], y[sel], label=r"$\delta^{A \rightarrow G}_{S, S'}$", 
            zorder=10, c=c, linewidth=lwd)
    ax.legend(loc="upper left")
    # ax.grid()
    ylim = ax.get_ylim()
    ymax = np.max(np.abs(ylim))
    if ymax >= 30:
        ylim = [-ymax, ymax]
    else:
        ylim = [-30, 30]
        
    ax.set_ylim(ylim)
    ax.set_xlim([s, e])
    ax.set_ylabel(r"$\delta^{A \rightarrow G}_{S, S'}$")
    gridposx = ax.get_xticks()
    xticks = ax.get_xticks()
    xlabels = [pstart]
    for i in range(1, len(xticks)):
        xlabels.append("+" + str(int(xticks[i])))
    ax.set_xticks(ax.get_xticks(), xlabels)
    # print(gridpos)
    for x in gridposx:
        ax.plot([x, x], [-ylim[1], 2 * ylim[1]], c=grid_col,
                linestyle=grid_style, zorder=2, linewidth=grid_lwd)
    gridposy = ax.get_yticks()
    for y in gridposy:
        ax.plot([gridposx[0], gridposx[-1]], [y, y], c=grid_col,
                linestyle=grid_style, zorder=2, linewidth=grid_lwd)


def make_plots(args):
    # check if output directory exists and generate if not
    if os.path.exists(args["output"]) is False:
        os.makedirs(args["output"])
    # loading all required data
    fasta = pysam.FastaFile(args["fasta"])
    bedlike = pd.read_table(args["regions"], header=None, comment="#")
    
    # opening connections to BAM files
    bam_g0 = {}
    for s in args["group1"]:
        bam_g0[s] = pysam.AlignmentFile(s, "rb")

    bam_g1 = {}
    for s in args["group2"]:
        bam_g1[s] = pysam.AlignmentFile(s, "rb")
    
    total = bedlike.shape[0]
    flanking_space = 500
    grid_lwd = .75
    grid_style = (0, (1, 5))
    lwd = 1.9
    for i in range(total):
        chrom = str(bedlike.iloc[i, 0])
        start = bedlike.iloc[i, 1] - flanking_space
        end = bedlike.iloc[i, 2] + flanking_space
        strand = bedlike.iloc[i, 5]
        name = bedlike.iloc[i, 3]
        
        try:
            regions_g0, regions_g1, refseq = eif.fetch_regions(bam_g0, bam_g1, fasta, chrom, start, end)
            # print("Fetching regions Done")
            nuc_ref = "A"
            nuc_edit = "G"
            if strand == "-":
                nuc_ref = "T"
                nuc_edit = "C"
            
            data = eif.get_data(regions_g0, regions_g1, refseq, start, end, strand, 
                                nuc_ref=nuc_ref, nuc_edit=nuc_edit,
                                library=args["library"], debug=False)
            # print("get_data() Done")
            df = eif.calc_local_EI(data, hwidth=args["window_size"], step_size=args["step_size"],
                                   mincoverage=args["min_coverage"], rm_snp=args["rm_snps"])
            # print("calc_local_EI() Done")

            fig, ax = plt.subplots(3, 1)
            fig.subplots_adjust(hspace=0.0, wspace=0.05)
            fig.set_size_inches(10, 6)

            plot_single_signals(ax[0], df, grid_lwd=grid_lwd, grid_style=grid_style,
                                lwd=lwd, pstart=start, nox=True)
            plot_average_signals(ax[1], df, grid_lwd=grid_lwd, grid_style=grid_style,
                                 pstart=start, nox=True, lwd=lwd)
            plot_group_differences(ax[2], df, grid_lwd=grid_lwd, grid_style=grid_style,
                                   pstart=start, lwd=lwd)
            # plot_lEIregions(ax[2], regions)

            fig.savefig(f"{args['output']}/{i}_{name}_{start}_{end}_3.png", format="png",
                        dpi=300, facecolor="#FFFFFF", bbox_inches="tight", pad_inches=0.1)
            plt.close(fig)
            
            fig, ax = plt.subplots()
            fig.set_size_inches(10, 2.5)
            plot_group_differences(ax, df, grid_lwd=grid_lwd, grid_style=grid_style,
                                   pstart=start, lwd=lwd)
            fig.savefig(f"{args['output']}/{i}_{name}_{start}_{end}_2.png", format="png",
                        dpi=300, facecolor="#FFFFFF", bbox_inches="tight", pad_inches=0.1)
            plt.close(fig)
            
            fig, ax = plt.subplots()
            fig.set_size_inches(10, 2.5)
            plot_average_signals(ax, df, grid_lwd=grid_lwd, grid_style=grid_style,
                                 pstart=start, lwd=lwd)
            fig.savefig(f"{args['output']}/{i}_{name}_{start}_{end}_1.png", format="png",
                        dpi=300, facecolor="#FFFFFF", bbox_inches="tight", pad_inches=0.1)
            plt.close(fig)
            
            fig, ax = plt.subplots()
            fig.set_size_inches(10, 2.5)
            plot_single_signals(ax, df, grid_lwd=grid_lwd, grid_style=grid_style,
                                pstart=start, lwd=lwd)
            fig.savefig(f"{args['output']}/{i}_{name}_{start}_{end}_0.png", format="png",
                        dpi=300, facecolor="#FFFFFF", bbox_inches="tight", pad_inches=0.1)
            plt.close(fig)
            
        except Exception:
            print(f"Could not find regions for {name}: {chrom}, {start}, {end}")
        
    # closing all file connections
    for s in bam_g0:
        bam_g0[s].close()
    for s in bam_g1:
        bam_g1[s].close()
