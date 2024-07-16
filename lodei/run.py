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
import argparse
import os
import sys
import lodei.core.find as find
import lodei.core.plot_regions as plotregions


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
        help="Specifies the sequencing library type. [Default: SR]")

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


if __name__ == "__main__":
    main()
