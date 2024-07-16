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
import pysam
import pandas as pd
import numpy as np
import sys
import os
import time
import logging
import multiprocessing as mp
import lodei.core.eifunctions as eif
from datetime import datetime
from subprocess import Popen, PIPE


def find(args):
    if args["window_size"] < 10 or args["window_size"] > 50:
        raise Exception(f"Command line argument -w/--window_size needs to be an integer in between 10 and 50.")
    args["pairs"] = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
    args["step_size"] = 2*args["window_size"]+1
    path = os.path.abspath(args["output"])
    if args["self"] is False:
        # check if output directory exists and generate if not
        if os.path.exists(path) is False:
            os.makedirs(path)
        if os.path.exists(f"{path}/windows") is False:
            os.makedirs(f"{path}/windows")

    time_start = datetime.now()
    tstart = time.time()
    time_start_str1 = time_start.strftime("%Y-%m-%d %H:%M:%S")

    # Start single process version. Note, the multi-process version simply
    # splits the input data and calls this function again with one core per splitted input data
    if args["cores"] == 1:
        # create log file
        logging.basicConfig(filename=f"{path}/find_{args['subprocessid']}.log",
                            filemode='w', level=logging.INFO)

        print_info(args)
        logging.info(f"Start: {time_start_str1}\n")
        logging.info(f"Arguments: {args}")
        run(args)
    elif args["cores"] > 1:
        logging.basicConfig(filename=f"{path}/find_main.log",
                            filemode='w', level=logging.INFO)
        print_info(args)
        run_mp(args)
    else:
        raise Exception(f"Command line argument -c/--cores {args['cores']} provides an unvalid number of cores.")
    os.system(f"mv {path}/lei_windows* {path}/windows")
    logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Data collecition completed.")
    if args["self"] is False:
        logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Loading window data.")
        signals = eif.load_signals(f"{path}/windows")
        os.system(f"rm {path}/windows/*")  # remove old tables
        logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Loading window data completed.")

        ###############################################
        # averaged q-value version:
        logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: q-value calculation and window filtering.")
        values_signal = np.arange(0, 50, .1)
        qvalues = {}
        for p in args["pairs"]:
            logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Calculating q-values for {p}")
            qvalues[p], qd = eif.get_qvalues(signals, values_signal, index_signal=p, correct=False)
            # calculate q-values for windows and save data
            eif.add_qvalue_to_df(signals[p], qd)
            signals[p].to_csv(f"{path}/windows/windows_{p}.txt", header=True, index=False, sep="\t")

        for p in args["pairs"]:
            cutoff_neg, cutoff_pos = eif.find_signal_cutoff_from_qvalue(qvalues[p], qcutoff=.1)
            counts_neg = 0
            counts_pos = 0

            counts_neg = np.sum(signals[p]["wEI"] <= cutoff_neg)
            counts_pos = np.sum(signals[p]["wEI"] >= cutoff_pos)
            tmp = signals[p].loc[(signals[p]["wEI"] <= cutoff_neg) | (signals[p]["wEI"] >= cutoff_pos), :]
            # print(f'{p}: {cutoff_neg}, {cutoff_pos} | #neg: {counts_neg}, #pos: {counts_pos}, #total: {tmp.shape[0]}')
            logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: cutoffs {p}: {cutoff_neg}, {cutoff_pos} | #neg: {counts_neg}, #pos: {counts_pos}, #total: {tmp.shape[0]}")
            if tmp.shape[0] >= 1:
                tmp.to_csv(f"{path}/windows_qfiltered_{p}.txt", header=True, index=False, sep="\t")
            else:
                logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: No windows passed the q-value filter for {p}.")
        logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: q-value calculation and window filtering completed")

    os.system(f"rm -rf {path}/tmp")  # delete temporary files
    tstop = time.time()
    time_stop = datetime.now()
    time_stop_str1 = time_stop.strftime("%Y-%m-%d %H:%M:%S")
    logging.info("\n")
    logging.info(f"All done: {time_stop_str1}. Total time in hours: {np.round((tstop-tstart)/3600, 4)}\n")


def run(args):
    path = os.path.abspath(args["output"])
    logging.debug("Loading/Opening required file connections...")
    # open connection to save data
    connections_out = {}  # a dict that keeps the connections for each pair
    for p in args["pairs"]:
        if args["self"] is False:
            connections_out[p] = open(f"{path}/windows/lei_windows_{p}.txt", "w")
        else:
            connections_out[p] = open(f"{path}/lei_windows_{p}_{args['subprocessid']}.txt", "w")
    # loading all required data
    fasta = pysam.FastaFile(args["fasta"])
    gff = pd.read_table(args["gff"], header=None, comment="#")

    # opening connections to BAM files
    bam_g0 = {}
    for s in args["group1"]:
        logging.debug(f"Opening {s}")
        bam_g0[s] = pysam.AlignmentFile(s, "rb")

    bam_g1 = {}
    for s in args["group2"]:
        logging.debug(f"Opening {s}")
        bam_g1[s] = pysam.AlignmentFile(s, "rb")

    logging.debug("...finished opening file connections.")

    total = gff.shape[0]
    total_nucs = 0
    flanking_space = 500
    for i in range(total):
        chrom = str(gff.iloc[i, 0])
        start = gff.iloc[i, 3] - flanking_space
        end = gff.iloc[i, 4] + flanking_space
        strand = gff.iloc[i, 6]
        total_nucs += end - start
        # Try to get gene name from GTF. Currently only works for human gencode annotation.
        try:
            name = gff.iloc[i, 8].split("gene_name=")[1].split(";")[0]
        except Exception:
            name = str(i + 1)

        regions_g0, regions_g1, refseq = eif.fetch_regions(bam_g0, bam_g1, fasta, chrom, start, end)
        for p in args["pairs"]:
            # the data and local stuff needs to get done for each pair
            data = eif.get_data(regions_g0, regions_g1, refseq, start, end, strand,
                                nuc_ref=p[0], nuc_edit=p[1],
                                library=args["library"], debug=False)

            df = eif.calc_local_EI(data, hwidth=args["window_size"], step_size=args["step_size"],
                                   mincoverage=args["min_coverage"], rm_snp=args["rm_snps"])

            eif.report_localei(df, chrom, start, name, strand, args["window_size"], connections_out[p])

    # closing all file connections
    for p in args["pairs"]:
        connections_out[p].close()

    for s in bam_g0:
        bam_g0[s].close()
    for s in bam_g1:
        bam_g1[s].close()


def run_mp(args):
    # print("Run MP version!")
    logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Run MP version with {args['cores']} cores.")
    # check number of available cores and adjust if necessary
    if args["cores"] > mp.cpu_count():
        args["cores"] = mp.cpu_count()

    # check if output directory exists and generate if not
    path = os.path.abspath(args["output"])
    tmppath = path + "/tmp"
    # check if tmp directory already exist and remove all content from it.
    if os.path.exists(tmppath):
        print(f"Temporary directory {tmppath} already exists. Removing complete content in that directory!")
        lfiles = os.listdir(tmppath)
        for f in lfiles:
            os.remove(os.path.join(tmppath, f))
    else:
        os.makedirs(path + "/tmp")

    logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Splitting annotation for subprocesses")
    # split input data into chunks
    gff = pd.read_table(args["gff"], header=None, comment="#")
    total = gff.shape[0]
    if total < args["cores"]:
        args["cores"] = total
    chunk_size = int(total / args["cores"])
    start = 0
    fcount = 0
    for i in range(args["cores"]):
        if i < args["cores"] - 1:
            tmp = gff.iloc[start:(start + chunk_size), :]
        else:
            tmp = gff.iloc[start:total, :]
        tmp.to_csv(path + "/tmp/tmp_" + str(args["id"]) + "_" + str(fcount) + ".gff",
                   sep="\t", header=False,
                   index=False)
        start += chunk_size
        fcount += 1

    # start subprocesses for each tmp gff
    logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Starting subprocesses.")
    process_list = []
    for i in range(fcount):
        cmd = ["lodei", "find", "-a"]
        for s in args["group1"]:
            cmd.append(s)

        cmd.append("-b")
        for s in args["group2"]:
            cmd.append(s)

        cmd.append("-f")
        cmd.append(args["fasta"])

        cmd.append("-g")
        cmd.append(path + "/tmp/tmp_" + str(args["id"]) + "_" + str(i) + ".gff")

        cmd.append("-o")
        cmd.append(path + "/tmp")

        cmd.append("-i")
        cmd.append(str(args["id"]))

        cmd.append("-w")
        cmd.append(str(args["window_size"]))

        cmd.append("-m")
        cmd.append(str(args["min_coverage"]))

        if args["rm_snps"]:
            cmd.append("--rm_snps")

        cmd.append("--library")
        cmd.append(args["library"])

        cmd.append("--self")

        cmd.append("--subprocessid")
        cmd.append(str(i))

        process_list.append(Popen(cmd, stdout=PIPE, stderr=PIPE))

    # wait until jobs are finished
    if len(process_list) > args["cores"]:
        print("More process than available cores... Exiting.")
        sys.exit()
    # print(f"Number of processes started: {len(process_list)}")
    for p in process_list:
        p.wait()
    logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Subprocess completed.")
    logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Merge subprocess results.")
    # merge tmp files:
    for p in args["pairs"]:
        # merge temp results into final result file
        for i in range(fcount):
            os.system(f"cat {path}/tmp/lei_windows_{p}_{i}.txt >> {path}/lei_windows_{p}.txt")
    logging.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Merge subprocess results completed.")

    # copy log files into output directory:
    os.makedirs(path + "/logs")
    os.system("mv " + path + "/tmp/*.log " + path + "/logs")


def print_info(args):
    msg = "\nlodei Copyright (C) 2024.\n\n"
    msg += "This program comes with ABSOLUTELY NO WARRANTY.\n"
    msg += "This is free software (GNU GPLv3), and you are welcome to redistribute it under GPLv3 conditions.\n\n"
    msg += "Arguments:\n"
    msg += "group1:\n"
    for f in args["group1"]:
        msg += f"{f}\n"
    msg += "\n"
    msg += "group2:\n"
    for f in args["group2"]:
        msg += f"{f}\n"
    msg += "\n"
    msg += f"output      : {args['output']}\n"
    msg += f"gff         : {args['gff']}\n"
    msg += f"fasta       : {args['fasta']}\n"
    msg += f"cores       : {args['cores']}\n"
    msg += f"window_size : {args['window_size']}\n"
    msg += f"library type: {args['library']}\n"
    msg += f"min-coverage: {args['min_coverage']}\n"
    msg += f"remove SNPs : {args['rm_snps']}\n"
    print(msg)
    logging.info(msg)
