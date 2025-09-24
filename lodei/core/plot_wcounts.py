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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

mpl.rcParams['font.size'] = 8
mpl.rcParams['grid.color'] = '#b0b0b0'
mpl.rcParams['grid.linestyle'] = (0, (2, 3))
mpl.rcParams['grid.linewidth'] = 0.8
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['legend.markerscale'] = 1
mpl.rcParams['legend.fontsize'] = 8

lwd = 2.5

# finds the min and max LoEI values of all mismatch pairs of an experiment.
# these limits are needed for gather_absolute_counts()
def find_signal_limits(signals):
    pairs = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
    min_signal = 0
    max_signal = 0
    for p in pairs:
        tmp_min = min(signals[p]["wEI"])
        tmp_max = max(signals[p]["wEI"])
        if tmp_min < min_signal:
            min_signal = tmp_min
        if tmp_max > max_signal:
            max_signal = tmp_max
    return min_signal, max_signal

# gathers how many windows are existing for a given range of LoDEI singals.
# So how many windows are < a th and how many are > a th.
def gather_absolute_counts(wei_signal, limits=[]):
    wei_values = np.round(np.arange(limits[0], limits[1]+1, .1), 1)
    counts = []
    for x in wei_values:
        if x > 0:
            counts.append(np.sum(wei_signal >= x)+1)
        elif x < 0:
            counts.append(np.sum(wei_signal <= x)+1)
        else:
            counts.append(1)
    return pd.DataFrame({"n": counts,
                         "limits": wei_values})

def plot_counts(ax, signals, title=""):
    pairs = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]

    handles = []
    spacer = Line2D([], [], linestyle='none', marker=None, label="")

    lpairs = [["AC", "AG", "AT"],
              ["CA", "CG", "CT"],
              ["GA", "GC", "GT"],
              ["TA", "TC", "TG"]]
    
    lcmaps = [plt.get_cmap("Reds"),
              plt.get_cmap("Blues"),
              plt.get_cmap("Greys"),
              plt.get_cmap("Greens")]

    for pairs, cmap in zip(lpairs, lcmaps):
        for i, p in enumerate(pairs):
            c = (i+1)*3
            line, = ax.plot(signals[p]["limits"], signals[p]["n"], 
                    label=p, c=cmap(c/11), linewidth=2)
            handles.append(line)
        handles.append(spacer)
    
    del handles[len(handles)-1]
    ax.set_yscale("log")
    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, loc="best")
    ax.set_title(title)
    ax.set_xlabel(r"$\delta^{X \rightarrow Y}$")
    ax.set_ylabel(r"#windows $\leq \delta^{X \rightarrow Y}$ for $\delta^{X \rightarrow Y} < 0$"+"\n"+r"#windows$\geq \delta^{X \rightarrow Y}$ for $\delta^{X \rightarrow Y} > 0$")
    ax.grid()
    ylim = ax.get_ylim()
    #print(ylim)
    ylim = (.8,ylim[1])
    ax.plot([0,0], ylim, c="#b0b0b0", zorder=4, linewidth=2.5)
    ax.set_ylim(ylim)

def make_plot(args):
    pairs = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]

    signals = {}
    for p in pairs:
        result_file = os.path.join(args["directory"], f"windows_{p}.txt")
        signals[p] = pd.read_table(result_file)

    pdata = {}
    
    min_s, max_s = find_signal_limits(signals)
    min_s = round(min_s) - 1
    max_s = round(max_s) + 1

    for p in pairs:
        pdata[p] = gather_absolute_counts(signals[p]["wEI"], limits=[min_s, max_s])



    fig, ax = plt.subplots()
    plot_counts(ax, pdata)
    fig.savefig(args["output"], facecolor='#FFFFFF', bbox_inches = 'tight', pad_inches = 0.1)

