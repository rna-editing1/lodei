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
import pysamstats
import pandas as pd
import numpy as np
import logging
import os
from datetime import datetime
import math


def get_nucleotides(ref, edit, strand, library="SR"):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    # single stranded. + genes have signal from reverse reads, - genes from forward reads
    if library == "SR":
        if strand == "+":
            return ref, edit, "_rev"
        elif strand == "-":
            return complement[ref], complement[edit], "_fwd"
    if library == "U":  # unstranded library
        if strand == "+":
            return ref, edit, ""
        elif strand == "-":
            return complement[ref], complement[edit], ""
    if library == "SF":
        if strand == "+":
            return ref, edit, "_fwd"
        elif strand == "-":
            return complement[ref], complement[edit], "_rev"
    # Return None if an unknown library-type is provided
    return None


def fetch_regions(bams_g0, bams_g1, fasta, chrom, start, end):
    logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Fetching data for {chrom},{start},{end}")
    regions_g0 = {}
    regions_g1 = {}
    take_ref = True
    ref = []
    for k, bam in bams_g0.items():
        region = pysamstats.load_variation_strand(
            bam, chrom=chrom,
            start=start, end=end, truncate=True,
            fafile=fasta, one_based=True, pad=True)
        regions_g0[k] = region
        if take_ref:
            ref = np.array(region["ref"], dtype="str")
            take_ref = False
    
    for k, bam in bams_g1.items():
        region = pysamstats.load_variation_strand(
            bam, chrom=chrom,
            start=start, end=end, truncate=True,
            fafile=fasta, one_based=True, pad=True)
        regions_g1[k] = region
    
    return regions_g0, regions_g1, ref


def get_data(regions_g0, regions_g1, ref, start, end, strand, nuc_ref, nuc_edit, library, debug=False):
    nuc, mut, st = get_nucleotides(nuc_ref, nuc_edit, strand, library=library)
    # Fetch data for samples of group 0
    # create empty first line
    mat_cov_0 = np.full(end - start, 0)
    mat_mut_0 = np.full(end - start, 0)

    for k, region in regions_g0.items():
        # print(f"getting data for {k}")
        # print(f'len matches: {region[f"matches{st}"]}')
        # add coverage and mutation to the corresponding matrix for each sample
        # mat_cov_0 = np.vstack((mat_cov_0, region[f"reads{st}"]))
        mat_cov_0 = np.vstack((mat_cov_0, region[f"matches{st}"] + region[f"mismatches{st}"]))
        # mat_cov_0 = np.vstack((mat_cov_0, region[f"reads_pp{st}"]))
        mat_mut_0 = np.vstack((mat_mut_0, region[f"{mut}{st}"]))
    # print(f"mat_cov_0.shape: {mat_cov_0.shape}")
    # print(f"mat_mut_0.shape: {mat_mut_0.shape}")
    # remove empty first line
    mat_cov_0 = mat_cov_0[1::, :]
    mat_mut_0 = mat_mut_0[1::, :]

    # group 1
    # create empty first line
    mat_cov_1 = np.full(end - start, 0)
    mat_mut_1 = np.full(end - start, 0)
    
    for k, region in regions_g1.items():
        # add coverage and mutation to the corresponding matrix for each sample
        mat_cov_1 = np.vstack((mat_cov_1, region[f"matches{st}"] + region[f"mismatches{st}"]))
        mat_mut_1 = np.vstack((mat_mut_1, region[f"{mut}{st}"]))
    
    # remove first line
    mat_cov_1 = mat_cov_1[1::, :]
    mat_mut_1 = mat_mut_1[1::, :]
    
    # print per position.... for debugging:
    if debug:
        for i in range(mat_cov_0.shape[1]):
            line = f"{start+i:>15} {ref[i]} "
            for j in range(mat_cov_0.shape[0]):
                line += f" {mat_mut_0[j,i]}/{mat_cov_0[j,i]}\t"
            print(f"{line}")
    
    # print(ref)
    # print(mat_cov_0)
    # print(mat_mut_0)
    # print("")
    # print(mat_cov_1)
    # print(mat_mut_1)
    
    dsignal = {"ref": ref, "cov0": mat_cov_0, "cov1": mat_cov_1,
               "mut0": mat_mut_0, "mut1": mat_mut_1, "nuc": nuc}
    return dsignal


# simple filter to remove potential SNPs.
def get_non_snp_positions(m, r=.8, n=0):
    # count the number of ratios >= r for each column of m
    esnp_counts = np.array([np.sum(m[:, i] > r) for i in range(m.shape[1])])
    columns = esnp_counts == n
    return columns


# gets the dict as an input
def get_EI(ref, c0, c1, m0, m1, nuc, mincoverage, debugstart=-1, rm_snp=False):
    # print(ref)
    # print(c0)
    # print(c1)
    # print("")
    # print(m0)
    # print(m1)
    
    # We are done, if no reference nucleotides are in the requested area:
    sel_ref = ref == nuc
    if np.sum(sel_ref) < 1:
        logging.debug("No reference nucleotide in the requested region.")
        return None
    
    # check min coverage for all columns of both matrices:
    sel_cov_0 = np.array([np.all(c0[:, i] >= mincoverage) for i in range(c0.shape[1])])
    sel_cov_1 = np.array([np.all(c1[:, i] >= mincoverage) for i in range(c1.shape[1])])
    
    # print(sel_ref)
    # print(sel_cov_0)
    # print(sel_cov_1)
    # print(np.where(sel_ref & sel_cov_0))
    # print(np.where(sel_ref & sel_cov_1))
    # print(ref[sel_ref & sel_cov_0])
    # print(c0[:, sel_ref & sel_cov_0])
    # print(c1[:, sel_ref & sel_cov_1])
    
    # only keep columns that fulfill all constraints:
    sel = sel_ref & sel_cov_0 & sel_cov_1
    # print(np.where(sel))
    # Again, if no positions fulfill all constraints, no EI can be calculated
    if np.sum(sel) < 1:
        logging.debug("No positions are left that fulfill all required constraints.")
        return None
    ref = ref[sel]
    c0 = c0[:, sel]
    c1 = c1[:, sel]
    m0 = m0[:, sel]
    m1 = m1[:, sel]
    
    # print for debugging:
    if debugstart >= 0:
        mpos = np.where(sel)
        for i in range(c0.shape[1]):
            line = f"{i:>4} {mpos[0][i]+debugstart:>10} {ref[i]} -> "
            # line = f"{i:>4} "
            for j in range(c0.shape[0]):
                line += f"{m0[j,i]}/{c0[j,i]}\t"
            line += " || "
            for j in range(c1.shape[0]):
                line += f"{m1[j,i]}/{c1[j,i]}\t"
            print(line)
    
    # print(c0)
    # print(m0)
    # print("-")
    # print(c1)
    # print(m1)
    
    # calculate EI for group 0
    ei_0 = []
    csn0 = []  # count #signal nucleotides
    ei_1 = []
    csn1 = []  # count #signal nucleotides
    rmat_0 = m0 / c0 * 100
    rmat_1 = m1 / c1 * 100

    sel_snps_0 = np.full(rmat_0.shape[1], True)
    sel_snps_1 = np.full(rmat_1.shape[1], True)
    sel_snps = np.full(rmat_0.shape[1], True)
    
    if rm_snp:
        sel_snps_0 = get_non_snp_positions(rmat_0, r=85)
        sel_snps_1 = get_non_snp_positions(rmat_1, r=85)
        sel_snps = sel_snps_0 & sel_snps_1  # select only those positions that are no SNPs in both groups

    if np.sum(sel_snps) < 1:
        logging.debug("No positions are left that fulfill all required constraints.")
        return None

    rmat_0 = rmat_0[:, sel_snps]
    for i in range(rmat_0.shape[0]):  # for a sample
        ei_0.append(np.sum(rmat_0[i, :]))
        csn0.append(np.sum(rmat_0[i, :] > 0))
    
    # calculate EI for group 1
    rmat_1 = rmat_1[:, sel_snps]
    for i in range(rmat_1.shape[0]):  # for a sample
        ei_1.append(np.sum(rmat_1[i, :]))
        csn1.append(np.sum(rmat_1[i, :] > 0))
    av1 = np.mean(ei_1)
    av0 = np.mean(ei_0)
    
    # get average number signal contributing positions
    if av1 - av0 > 0:
        csn = int(np.round(np.mean(csn1)))
        csn_float = np.round(np.mean(csn1), 2)
    else:
        csn = int(np.round(np.mean(csn0)))
        csn_float = np.round(np.mean(csn0), 2)
    return {"ei0": ei_0, "ei1": ei_1, "av0": av0, "av1": av1, "av_diff": av1 - av0, "cn": c0.shape[1], "csn": csn, "csn_float": csn_float}


def calc_local_EI(data, hwidth=75, step_size=25, mincoverage=5, rm_snp=False):
    length = len(data['ref']) + 1
    
    pos = []      # genomic coordinates
    av_diff = []  # averages between the two groups
    wstarts = []
    wends = []
    
    # lists that store EI values per sample. Only used for plotting for QC and visualization
    g0 = [[] for i in range(data["cov0"].shape[0])]
    g1 = [[] for i in range(data["cov1"].shape[0])]
    av_g0 = []
    av_g1 = []
    # null = []
    cn = []   # count reference nucleotides
    csn = []  # count signal contributing nucleotides
    csn_float = []  # count signal contributing nucleotides
    for i in range(hwidth, length - hwidth, step_size):
        pos.append(i)
        # print(f"{(i-hwidth)}:{i}:{(i+hwidth+1)}")
        # window = np.arange((i-hwidth), (i+hwidth+1),1)
        wstart = (i - hwidth)
        wstop = (i + hwidth + 1)
        wstarts.append(wstart)
        wends.append(wstop - 1)  # -1 because, the last position using a slice notation is not selected
        ei = get_EI(data["ref"][wstart:wstop],
                    data["cov0"][:, wstart:wstop],
                    data["cov1"][:, wstart:wstop],
                    data["mut0"][:, wstart:wstop],
                    data["mut1"][:, wstart:wstop],
                    data["nuc"], mincoverage, rm_snp=rm_snp)
        if ei is None:
            # av_diff.append(0)
            av_diff.append(np.nan)
            # add sample data
            av_g0.append(0)
            av_g1.append(0)
            for j in range(len(g0)):
                g0[j].append(0)
            for j in range(len(g1)):
                g1[j].append(0)
            # null.append(np.nan)
            cn.append(-1)
            csn.append(-1)
            csn_float.append(-1)
        else:
            av_diff.append(ei["av_diff"])
            
            # add sample data
            av_g0.append(ei["av0"])
            av_g1.append(ei["av1"])
            for j in range(len(ei["ei0"])):
                g0[j].append(np.round(ei["ei0"][j], 3))
            for j in range(len(ei["ei1"])):
                g1[j].append(np.round(ei["ei1"][j], 3))
            # null.append(ei["av_diff"])
            cn.append(ei["cn"])
            csn.append(ei["csn"])
            csn_float.append(ei["csn_float"])
        
    result = {"x": pos, "av_diff": np.array(av_diff)}
    result["av_g0"] = av_g0
    result["av_g1"] = av_g1
    for j in range(len(g0)):
        result[f"g0_{j}"] = g0[j]
    for j in range(len(g1)):
        result[f"g1_{j}"] = g1[j]
    # result["null"] = null
    result["cn"] = cn
    result["csn"] = csn
    result["csn_float"] = csn_float
    result["wstart"] = wstarts
    result["wend"] = wends
    return pd.DataFrame(result)


# chrom, start, name, strand, regions, hwsize, stream
def report_localei(df, chrom, start, name, strand, hwsize, stream):
    for i, row in df.iterrows():
        if not np.isnan(row["av_diff"]):
            line = f"{chrom}\t{start+df['wstart'][i]}\t{start+df['wend'][i]}\t"
            line += f"{name}\t{df['av_diff'][i]:.1f}\t{strand}\n"
            # line += f"{name}\t{df['av_diff'][i]:.1f}\t{strand}\t{df['cn'][i]}\t{df['csn'][i]}\t{df['csn_float'][i]}"  # \n"
            # line += f"\t{df['av_g0'][i]:.1f}\t{df['av_g1'][i]:.1f}\n"
            stream.write(line)
    stream.flush()


def load_signals(path):
    path = os.path.abspath(path)
    signals = {}
    pairs = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
    for p in pairs:
        tmp = pd.read_table(f"{path}/lei_windows_{p}.txt")
        # tmp.columns = ["chrom", "wstart", "wend", "name", "wEI", "strand", "cn", "csn", "csn_float", "av_0", "av_1"]
        tmp.columns = ["chrom", "wstart", "wend", "name", "wEI", "strand"]
        signals[p] = tmp
    return signals


def merge_windows(windows):
    chrom = []
    start = []
    stop = []
    strand = []
    name = []
    mEI = []
    tmp_chrom = windows.iloc[0, 0]
    tmp_start = windows.iloc[0, 1]
    tmp_stop = windows.iloc[0, 2]
    tmp_strand = windows.iloc[0, 5]
    tmp_name = windows.iloc[0, 3]
    tmp_mEI = [windows.iloc[0, 4]]
    for i in range(1, windows.shape[0]):
        if windows.iloc[i, 0] == tmp_chrom and windows.iloc[i, 1] <= tmp_stop and windows.iloc[i, 5] == tmp_strand and windows.iloc[i, 3] == tmp_name:
            tmp_stop = windows.iloc[i, 2]
        else:
            chrom.append(tmp_chrom)
            start.append(tmp_start)
            stop.append(tmp_stop)
            strand.append(tmp_strand)
            name.append(tmp_name)
            mEI.append(np.round(np.mean(tmp_mEI), 1))

            tmp_chrom = windows.iloc[i, 0]
            tmp_start = windows.iloc[i, 1]
            tmp_stop = windows.iloc[i, 2]
            tmp_strand = windows.iloc[i, 5]
            tmp_name = windows.iloc[i, 3]
            tmp_mEI = [windows.iloc[i, 4]]
    return pd.DataFrame({"chrom": chrom, "start": start,
                         "stop": stop, "name": name, "mEI": mEI, "strand": strand})


## returns the qvalues (neg and pos) for a given signal value
#def qvalue_for_signal_cutoff(signal, null, scutoff=1, correct=False):
#    s_pos = signal.loc[signal["wEI"] > 0, "wEI"]
#    s_neg = signal.loc[signal["wEI"] < 0, "wEI"]
#    
#    w_neg = 1
#    w_pos = 1
#    
#    if correct:
#        w_neg = np.sum(null["wEI"] < 0) / np.sum(s_neg < 0)
#        w_pos = np.sum(null["wEI"] > 0) / np.sum(s_pos > 0)
#    
#    a = int(np.sum(s_pos >= scutoff) * w_pos)
#    b = np.sum(null["wEI"] >= scutoff)  # FP
#    # print(f"pos: {b}/{a}")
#    if a <= 0:
#        if b > 0:
#            q_pos = 1
#        else:
#            q_pos = None
#    elif a <= 20:
#        q_pos = None
#    else:
#        q_pos = b / a
#    
#    a = int(np.sum(s_neg <= -scutoff) * w_neg)
#    b = np.sum(null["wEI"] <= -scutoff)
#    # print(f"neg: {b}/{a}")
#    if a <= 0:
#        if b > 0:
#            q_neg = 1
#        else:
#            q_neg = None
#    elif a <= 20:
#        q_neg = None
#    else:
#        q_neg = b / a
#    
#    if q_pos != None:
#        q_pos = np.min([1, q_pos])
#    
#    if q_neg != None:
#        q_neg = np.min([1, q_neg])
#    return q_neg, q_pos


## calculates the average q-value for a given mutation pair
## the average is computed based on the q-values for all other mutation pairs
## the function returns the average q-values for a list of signal values
#def get_qvalues(signals, values, index_signal="AG", correct=True):
#    values = np.round(values, 1)
#    q_pos = [1]*len(values)
#    q_neg = [1]*len(values)
#    pairs = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
#    if index_signal != "AG":
#        pairs = ["AC", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
#    bneg = 1.0
#    bpos = 1.0
#    qvdict = {}
#    for i, v in enumerate(values):
#        # print(f"{i}, {v}")
#        tmp_pos = []
#        tmp_neg = []
#        for index_null in pairs:
#            if index_signal != index_null:
#                qn, qp = qvalue_for_signal_cutoff(signals[index_signal], signals[index_null], scutoff=v, correct=correct)
#                if qn == None:
#                    tmp_neg.append(q_neg[-1])
#                else:
#                    tmp_neg.append(qn)
#            
#                if qp == None:
#                    tmp_pos.append(q_pos[-1])
#                else:
#                    tmp_pos.append(qp)
#        # print(tmp_neg)
#        # print(tmp_pos)
#        tmp = np.round(np.median(tmp_neg), 3)
#        if tmp < bneg:
#            bneg = tmp
#
#        tmp = np.round(np.median(tmp_pos), 3)
#        if tmp < bpos:
#            bpos = tmp
#        
#        q_neg[i] = bneg
#        q_pos[i] = bpos
#        qvdict[str(v)] = {"q_neg": bneg, "q_pos": bpos}
#        # print("")
#    qvdict["max"] = {"value": np.round(np.max(values), 1), "q_neg": bneg, "q_pos": bpos}
#    return pd.DataFrame({"values": values, "q_neg": q_neg, "q_pos": q_pos}), qvdict


## takes the output of get_qvalues and
## finds the signal cutoff for a given q-value
#def find_signal_cutoff_from_qvalue(qvdf, qcutoff=.1):
#    cut_neg = -math.inf
#    cut_pos = math.inf
#    for i, row in qvdf.iterrows():
#        if row["q_neg"] <= qcutoff:
#            cut_neg = -row["values"]
#            break
#    for i, row in qvdf.iterrows():
#        if row["q_pos"] <= qcutoff:
#            cut_pos = row["values"]
#            break
#    return cut_neg, cut_pos


#def add_qvalue_to_df(df, qvd):
#    tmp = np.ones(df.shape[0])
#    for i in range(df.shape[0]):
#        # print(f"{i}: {df.iloc[i,0]} | {df.iloc[i, 1]} | {df.iloc[i, 2]} | {df.iloc[i,3]} | {df.iloc[i,4]} | {df.iloc[i,5]} | ")
#        if df.iloc[i, 4] < 0:
#            if np.abs(df.iloc[i, 4]) > qvd["max"]["value"]:
#                tmp[i] = qvd["max"]["q_neg"]
#            else:
#                tmp[i] = qvd[str(np.abs(df.iloc[i, 4]))]["q_neg"]
#        elif df.iloc[i, 4] > 0:
#            if df.iloc[i, 4] > qvd["max"]["value"]:
#                tmp[i] = qvd["max"]["q_pos"]
#            else:
#                tmp[i] = qvd[str(df.iloc[i, 4])]["q_pos"]
#    df["q_value"] = tmp

def make_q_value_dict(signals, index_signal="AG"):
    # we start with the first x entries to get the first q-value
    # the q-value cannot decrease, but only increase from there on.
    s = np.array(signals[index_signal]["wEI"].copy())
    s.sort()
    th_neg = s[19]
    th_pos = s[len(s)-19]
    #print(f"th_neg: {th_neg}")
    #print(f"th_pos: {th_pos}")
    
    pairs = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
    FPs_neg = []
    FPs_pos = []
    for p in pairs:
        if p != index_signal:
            FPs_neg.append(np.sum(signals[p]["wEI"] <= th_neg))
            FPs_pos.append(np.sum(signals[p]["wEI"] >= th_pos))
        else:
            tp_fp_neg = np.sum(signals[p]["wEI"] <= th_neg)
            tp_fp_pos = np.sum(signals[p]["wEI"] >= th_pos)
    fp_neg = np.median(FPs_neg)
    fp_pos = np.median(FPs_pos)
    
    q_value_current_neg = min([1, round(fp_neg/tp_fp_neg, 4)])
    q_value_current_pos = min([1, round(fp_pos/tp_fp_pos, 4)])

    #print(f"neg: FP: {fp_neg:>10} | TP+FP: {tp_fp_neg} | q: {q_value_current_neg}")
    #print(f"pos: FP: {fp_pos:>10} | TP+FP: {tp_fp_pos} | q: {q_value_current_pos}")

    # now we have the starting q-values and can move from left to right (for negative values)
    # or from right to left (for positive values) to calculate all q-values
    qvdict = {}
    all_values = []
    all_qvalues = []
    values_neg = np.arange(th_neg, 0, 0.1)
    for value in values_neg:
        FPs_neg = []
        for p in pairs:
            if p != index_signal:
                FPs_neg.append(np.sum(signals[p]["wEI"] <= value))
            else:
                tp_fp_neg = np.sum(signals[p]["wEI"] <= value)
        fp_neg = np.median(FPs_neg)
        qvalue = round(fp_neg/tp_fp_neg, 4)
        if qvalue > q_value_current_neg:
            if qvalue > 1:
                qvalue = 1
            q_value_current_neg = qvalue
        qvdict[f"{round(value, 1)}"] = q_value_current_neg
        all_values.append(value)
        all_qvalues.append(q_value_current_neg)

    values_pos = np.arange(th_pos, 0, -0.1)
    for value in values_pos:
        FPs_pos = []
        for p in pairs:
            if p != index_signal:
                FPs_pos.append(np.sum(signals[p]["wEI"] >= value))
            else:
                tp_fp_neg = np.sum(signals[p]["wEI"] >= value)
        fp_pos = np.median(FPs_pos)
        qvalue = round(fp_pos/tp_fp_neg, 4)
        if qvalue > q_value_current_pos:
            if qvalue > 1:
                qvalue = 1
            q_value_current_pos = qvalue
        qvdict[f"{round(value, 1)}"] = q_value_current_pos
        all_values.append(value)
        all_qvalues.append(q_value_current_pos)
    return qvdict, pd.DataFrame({"values": all_values, "q_values": all_qvalues})

def qvalue2df(df, qvdict):
    values = [float(x) for x in list(qvdict.keys())]
    min_value = min(values)
    max_value = max(values)
    # print(f"{min_value}, {max_value}")
    q_value_neg_min_value = qvdict[str(min_value)]
    q_value_pos_max_value = qvdict[str(max_value)]
    # print(f"{q_value_neg_min_value}, {q_value_pos_max_value}")

    qvalues = [1]*len(df)
    for i in range(len(df)):
        wei = df.iloc[i, 4]
        if wei != 0:
            if wei > max_value:
                qvalues[i] = q_value_pos_max_value
            elif wei < min_value:
                qvalues[i] = q_value_neg_min_value
            else:
                qvalues[i] = qvdict[str(round(wei,1))]
    df["q_value"] = qvalues