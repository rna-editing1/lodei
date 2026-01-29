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
from matplotlib.lines import Line2D

mpl.rcParams['font.size'] = 8
mpl.rcParams['grid.color'] = '#b0b0b0'
mpl.rcParams['grid.linestyle'] = (0, (2, 3))
mpl.rcParams['grid.linewidth'] = 0.8
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['legend.markerscale'] = 1
mpl.rcParams['legend.fontsize'] = 8

lwd = 2.5


def get_ratios(signals, direction="-", max_value=80):
    pairs = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
    results = {"x": []}
    for p in pairs:
        results[p] = []

    if direction == "+":
        for i in range(1, max_value, 1):
            results["x"].append(i)
            tmp = []
            for p in pairs:
                r = int(np.sum((signals[p]["wEI"] >= i) & (signals[p]["wEI"] < i + 1)) + 1)
                results[p].append(r)
                tmp.append(r)
            m = float(np.median(tmp))
            for p in pairs:
                results[p][-1] = round(results[p][-1] / m, 2)

    elif direction == "-":
        for i in range(-max_value, 0, 1):
            results["x"].append(i)
            tmp = []
            for p in pairs:
                r = int(np.sum((signals[p]["wEI"] >= i) & (signals[p]["wEI"] < i + 1)) + 1)
                results[p].append(r)
                tmp.append(r)
            m = float(np.median(tmp))
            for p in pairs:
                results[p][-1] = round(results[p][-1] / m, 2)
    else:
        raise Exception("Unvalid direction")
    return results


def plot_ratios(ax, rpos, rneg, title="", y_max=1):
    pairs = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
    if y_max == 1:
        for p in pairs:
            y = np.max(np.log2(np.abs(rneg[p])))
            if y > y_max:
                y_max = y
            y = np.max(np.log2(np.abs(rpos[p])))
            if y > y_max:
                y_max = y

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
            c = (i + 1) * 3
            line, = ax.plot(rneg["x"], np.log2(rneg[p]),
                            label=p, c=cmap(c / 11), linewidth=2)
            line, = ax.plot(rpos["x"], np.log2(rpos[p]),
                            label=p, c=cmap(c / 11), linewidth=2)
            handles.append(line)
        handles.append(spacer)

    del handles[len(handles) - 1]

    labels = [h.get_label() for h in handles]
    ax.legend(handles, labels, loc="best")
    ax.set_xlabel(r"$\delta^{X \rightarrow Y}$")
    ax.set_ylabel("log2(r)")
    ax.grid()
    ax.set_ylim([-y_max - .5, y_max + .5])
    ax.set_title(title)


def make_plot(args):
    pairs = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]

    signals = {}
    for p in pairs:
        result_file = os.path.join(args["directory"], f"windows_{p}.txt")
        signals[p] = pd.read_table(result_file)

    ratios_neg = get_ratios(signals, direction="-", max_value=40)
    ratios_pos = get_ratios(signals, direction="+", max_value=40)

    fig, ax = plt.subplots()
    plot_ratios(ax, ratios_pos, ratios_neg, "", y_max=6)
    fig.savefig(args["output"], facecolor='#FFFFFF', bbox_inches='tight', pad_inches=0.1)
