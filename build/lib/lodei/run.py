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
import argparse
import os
import sys
import lodei.core.find as find
import lodei.core.plot_regions as plotregions
import lodei.core.plot_wcounts as plotwindowcounts
import lodei.core.plot_metagene as plotmetagene
from lodei.core import convert

def main():
    parser = argparse.ArgumentParser(
        prog="lodei",
        description="""lodei is a collection of command line applications to
        calculate the local difference of nucleotide editing signals between two
        groups of samples.""",
        epilog="More information: https://github.com/rna-editing1/lodei")

    subparsers = parser.add_subparsers(
        title="Commands",
        description="""List of available subcommands.
        To get a detailed description of each of the
        commands use $ lodei <command> -h""",
        help="",
        metavar="")

    # --- find ---------
    parser_find = subparsers.add_parser(
        name="find",
        description="""'lodei find' detects differentially edited regions
        between two groups of samples for a given set of genomic locations.""",
        help="Finds differentially edited regions.")

    parser_find.add_argument(
        "-a", "--group1", type=str, nargs="+", required=True,
        help="Input BAM files for group 1 separated by white space.")

    parser_find.add_argument(
        "-b", "--group2", type=str, nargs="+", required=True,
        help="Input BAM files for group 2 separated by white space.")

    parser_find.add_argument(
        "-f", "--fasta", dest="fasta", type=str, required=True,
        help="""Fasta file. Must be the same fasta file used for producing
        BAM files of group1 and group2.""")

    parser_find.add_argument(
        "-g", "--gff", dest="gff", type=str, required=True,
        help="""Annotation in GFF file format. The program tries
                to find local editing regions
                within the given GFF annotations.""")

    parser_find.add_argument(
        '-o', '--output', dest='output', required=True,
        help="Output directory")

    parser_find.add_argument(
        '-c', '--cores', dest='cores', type=int, default=1, required=False,
        help="Number of cores used. [Default: 1]")

    parser_find.add_argument(
        '-i', '--id', dest='id', type=int, default=1, required=False,
        help="Optional run ID. [Default: 1]")

    parser_find.add_argument(
        '-w', '--window_size', dest='window_size', type=int, default=25, required=False,
        help="""Half-window size for editing index calculation. Must be between 10 and 50.
        The complete window size is: 2*w+1 [Default: 25]""")

    parser_find.add_argument(
        '-l', '--library', dest='library', type=str, default="SR", required=False,
        help="""Specifies the strandedness and sequencing type of your BAM files.
    Allowed values:
        SR   = single-end, reverse stranded,
        SF   = single-end, forward stranded,
        U    = single-end, unstranded,
        ISR = paired-end, reverse stranded,
        ISF = paired-end, forward stranded,
    When using paired-end options (ISR, ISF), LoDEI will split the BAM files into separate _1 and _2 files for processing.
    [Default: SR]"""
    )

    parser_find.add_argument(
        '-m', '--min_coverage', dest='min_coverage', type=int, default=5, required=False,
        help="""Minimum required read coverage at regions. Only positions with a coverage >= m
                in all samples of all groups are used for editing index calculation. [Default: 5]""")

    parser_find.add_argument(
        '-p', '--rm_snps', dest='rm_snps', action="store_true", default=False,
        help="""Simple heuristic to remove possible SNPs [Default: off]""")

    parser_find.add_argument(
        '--self', dest='self', action="store_true", default=False,
        help="Only used if the the program calls itself for multi processing reasons. [Default: False]")

    parser_find.add_argument(
        '--subprocessid', dest='subprocessid', type=int, default=0,
        help="For creating a separate log file when using -c > 1. [Default: 0]")

    parser_find.add_argument(
        '-v', '--verbose', dest='verbose', action="store_true", default=False,
        help="Verbose mode. [Default: off]")

    parser_find.set_defaults(func=find.find)

    # --- plot_regions ---------
    parser_plot_regions = subparsers.add_parser(
        name="plotregion",
        help="Generate plots of editing signals.",
        description="""'lodei plotregion' generates plots of editing signals for
        two groups of samples for a given set of genomic locations. Genomic locations
        are provided in BED-like format.""",)

    parser_plot_regions.add_argument(
        "-a", "--group1", type=str, nargs="+", required=True,
        help="Input BAM files for group 1 separated by white space.")

    parser_plot_regions.add_argument(
        "-b", "--group2", type=str, nargs="+", required=True,
        help="Input BAM files for group 2 separated by white space.")

    parser_plot_regions.add_argument(
        "-f", "--fasta", dest="fasta", type=str, required=True,
        help="""Fasta file. Must be the same fasta file used for producing
        the BAM files for group1 and group2.""")

    parser_plot_regions.add_argument(
        "-p", "--plotregions", dest="regions", type=str, required=True,
        help="""File in a BED-like format that specifies the regions to plot.
        First 6 columns need to be: chromosome, start, end, name, score and strand. """)

    parser_plot_regions.add_argument(
        '-o', '--output', dest='output', required=True,
        help="Output directory. Plots will be automatically named using information in the provided -p/--plotregions file.")

    parser_plot_regions.add_argument(
        '-w', '--window_size', dest='window_size', type=int, default=25, required=False,
        help="""Half-window size for editing index calculation. Must be between 10 and 50.
        The complete window size is: 2*w+1 [Default: 25]""")

    parser_plot_regions.add_argument(
        '-r', '--ref', dest='reference', type=str, default="A", required=False,
        help="Reference nucleotide examined for editing index calculation. [Default: A]")

    parser_plot_regions.add_argument(
        '-e', '--edit', dest='editing', type=str, default="G", required=False,
        help="""Specifies the mutated nucleotide for positions with reference
                nucleotides as set by -r to calculate the editing index. [Default: G]""")

    parser_plot_regions.add_argument(
        '-l', '--library', dest='library', type=str, default="SR", required=False,
        help="Specficy the sequencing library type. [Default: SR]")

    parser_plot_regions.add_argument(
        '-m', '--min_coverage', dest='min_coverage', type=int, default=5, required=False,
        help="""Minimum required read coverage at regions. Only positions with a coverage >= m
                in all samples of all groups are used for editing index calculation. [Default: 5]""")

    parser_plot_regions.add_argument(
        '--rm_snps', dest='rm_snps', action="store_true", default=False,
        help="""Simple heuristic to remove possible SNPs [Default: off]""")

    parser_plot_regions.set_defaults(func=plotregions.make_plots)

    # --- plot_windowcounts ---------
    parser_plot_wcounts = subparsers.add_parser(
        name="plotwindowcounts",
        help="Plots the number of detected windows for all mismatch pairs.",
        description="""'lodei plotwindowcounts' plots for each possible mismatch pair
        the number of detected windows below and above all signal thresholds.
        This plots provides additional feedback/QC about the detected signals.""",)

    parser_plot_wcounts.add_argument(
        "-d", "--directory", type=str, required=True,
        help="The output directory generated by LoDEI that contains all 'window_XY.txt' files")

    parser_plot_wcounts.add_argument(
        "-o", "--output", type=str, required=True,
        help="Output filename.")

    parser_plot_wcounts.set_defaults(func=plotwindowcounts.make_plot)

    # --- plot_metagene ---------
    parser_plot_metagene = subparsers.add_parser(
        name="plotmetagene",
        help="Visualizes the distribution of significant editing windows along gene length.",
        description=(
            "Plots the relative position of editing events across genes for one or more input files."
            "Each input file (tab-separated, e.g. LoDEI output) is shown as a separate line. "
            "Gene annotation is taken from a GFF3 file. Only genes with unique (chromosome, gene_name) "
            "combinations are considered; ambiguous entries are excluded and reported. "
            "The gene body is divided into bins of equal size, and the number of windows per bin is plotted."
        )
    )
    parser_plot_metagene.add_argument(
        "-w", "--windows", type=str, nargs="+", required=True,
        help="Path(s) to one or more tab-separated files containing significant editing windows. "
        "Each file will be shown as a separate line in the plot."
    )
    parser_plot_metagene.add_argument(
        "-g", "--gff", type=str, required=True,
        help="Path to GFF3 file containing gene annotations. Only unique (chromosome, gene_name) pairs are used."
    )
    parser_plot_metagene.add_argument(
        "-n", "--num-bins", type=int, default=5,
        help="Number of bins to divide each gene into (default: 5, i.e. 20%% steps along the gene body)."
    )
    parser_plot_metagene.add_argument(
        "--ylim", type=float, nargs=2, metavar=('YMIN', 'YMAX'),
        help="Set y-axis limits, e.g. --ylim -40 -20"
    )
    parser_plot_metagene.add_argument(
        "--wEI-signal", action="store_true",
        help="Plot summed wEI signal per bin instead of number of significant windows."
    )
    parser_plot_metagene.add_argument(
        "-o", "--output", required=True,
        help="Output filename for the plot (e.g. metagene.png)."
    )
    parser_plot_metagene.set_defaults(func=plotmetagene.make_plot)

    # --- convert ---------
    parser_convert = subparsers.add_parser(
        name="convert",
        description="""'lodei convert' provides file format conversion utilities
        for converting lodei output files to standard bioinformatics formats.""",
        help="Convert lodei output files to standard formats.")

    convert_subparsers = parser_convert.add_subparsers(
        title="Conversion subcommands",
        description="Available file format conversions.",
        help="",
        metavar="")

    # --- convert windows2bedgraph ---------
    parser_windows2bedgraph = convert_subparsers.add_parser(
        name="windows2bedgraph",
        description="""Convert lodei windows output from lodei find to bedGraph format.
        The output bedGraph file can be used for visualization in genome browsers
        like IGV or UCSC Genome Browser.""",
        help="Convert windows output to bedGraph format.")

    parser_windows2bedgraph.add_argument(
        "-i", "--input", dest="input", type=str, required=True,
        help="Input windows file from 'lodei find' output.")

    parser_windows2bedgraph.add_argument(
        "-o", "--output", dest="output", type=str, required=True,
        help="Output bedGraph file.")

    parser_windows2bedgraph.set_defaults(func=convert.windows2bedgraph)


    # --- version ---------
    parser_version = subparsers.add_parser(
        name="version",
        help="Show program version.",
        description="Prints the current version of lodei."
    )
    parser_version.set_defaults(func=print_version)

    # --- testrun ---------
    parser_testrun = subparsers.add_parser(
        name="testrun",
        help="""Performs a test run on a given data.""",
        description="""Performs a test run on a predefined test data set to check
        if the the program runs as expected.

        Download the test run data set from: http:...
        Unpack the dowloaded archive to a directory of your choice.
        Do not modify anything inside the unpacked folder! The testrun command
        uses hardcoded filenames and expects the test data as is after unpacking it.

        Run $ lodei testrun /path/to/unpacked/testdata
        """)
    parser_testrun.add_argument(
        'directory', type=str,
        help="Path to the unpacked test data.")
    parser_testrun.set_defaults(func=testrun)

    args = vars(parser.parse_args())

    if "func" in args:
        args["func"](args)
    else:
        parser.print_help()


def testrun(args):
    print("lodei testrun: start")
    if os.path.exists(args["directory"]) is False:
        print(f"{args['directory']} does not exist. Exiting.")
        sys.exit(1)
    # create a testoutput folder
    outpath = os.path.join(args["directory"], "testrun_output")
    if os.path.exists(outpath) is False:
        os.mkdir(outpath)

    testrunargs = {'group1': [os.path.join(args["directory"], "bam/s01.bam"),
                              os.path.join(args["directory"], "bam/s02.bam"),
                              os.path.join(args["directory"], "bam/s03.bam"),
                              os.path.join(args["directory"], "bam/s04.bam"),
                              os.path.join(args["directory"], "bam/s05.bam")],
                   'group2': [os.path.join(args["directory"], "bam/s06.bam"),
                              os.path.join(args["directory"], "bam/s07.bam"),
                              os.path.join(args["directory"], "bam/s08.bam"),
                              os.path.join(args["directory"], "bam/s09.bam"),
                              os.path.join(args["directory"], "bam/s10.bam"),],
                   'fasta': os.path.join(args["directory"], "annotation/genome.fa"),
                   'gff': os.path.join(args["directory"], "annotation/test_anno.gff3"),
                   'output': f'{outpath}',
                   'cores': 3, 'id': 1, 'window_size': 50,
                   'library': 'SR',
                   'min_coverage': 5, 'rm_snps': False, 'subprocessid': 0,
                   'self': False, 'verbose': True}
    print("lodei testrun: 'lodei find' start")
    find.find(testrunargs)
    print("lodei testrun: 'lodei find' finished")

    testrunargs = {'group1': [os.path.join(args["directory"], "bam/s01.bam"),
                              os.path.join(args["directory"], "bam/s02.bam"),
                              os.path.join(args["directory"], "bam/s03.bam"),
                              os.path.join(args["directory"], "bam/s04.bam"),
                              os.path.join(args["directory"], "bam/s05.bam")],
                   'group2': [os.path.join(args["directory"], "bam/s06.bam"),
                              os.path.join(args["directory"], "bam/s07.bam"),
                              os.path.join(args["directory"], "bam/s08.bam"),
                              os.path.join(args["directory"], "bam/s09.bam"),
                              os.path.join(args["directory"], "bam/s10.bam"),],
                   'fasta': os.path.join(args["directory"], "annotation/genome.fa"),
                   'regions': os.path.join(args["directory"], "annotation/example_plot_regions.txt"),
                   'output': os.path.join(outpath, "plotregions"),
                   'window_size': 50,
                   'library': 'SR',
                   'min_coverage': 5, "rm_snps": True}
    print("lodei testrun: 'lodei plotregion' start")
    plotregions.make_plots(testrunargs)
    print("lodei testrun: 'lodei plotregion' finished")
    print("lodei testrun: finished")


def print_version(args):
    print("lodei version 1.1.1")


if __name__ == "__main__":
    main()
