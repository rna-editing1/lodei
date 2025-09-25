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
import glob
import os
import logging
from subprocess import run
import sys
import shutil
import pandas as pd

logger = logging.getLogger("paired_end")
logger.setLevel(logging.INFO)

def is_paired_end(library_value):
    """
    Checks whether the provided library value indicates paired-end data.
    Allowed values for paired-end are: 'ISR' and 'ISF'.
    """
    paired_values = ['ISR', 'ISF']
    return library_value.upper() in paired_values

def setup_paired_logger(output_dir, log_filename="paired_processing.log"):
    '''
    Sets up a dedicated logger named 'paired_end' that writes log messages to a specified file.

    This function ensures that the logger is configured to write logs to a file located in the given
    output directory. It removes any existing handlers to prevent duplicate logging and sets up a
    new file handler with a specific log format.

    Parameters:
        output_dir (str): The directory where the log file will be created. If the directory does not exist, it will be created.
        log_filename (str, optional): The name of the log file. Defaults to "paired_processing.log".

    Returns:
        logging.Logger: A configured logger instance named 'paired_end'.
    '''

    # 1. create output directory if it does not exist
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 2. path to the log file
    log_path = os.path.join(output_dir, log_filename)
 
    # 4. delete existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.propagate = False  # supress console output

    # 5. configure the logger
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger

def create_temp_folder(output_dir):
    """
    Creates a temporary folder within the given output directory and,
    within that folder, creates two subdirectories 'r1' and 'r2' for storing
    the split BAM files.

    If the 'tmp_paired' folder already exists, it is removed entirely and then recreated.

    Args:
        output_dir (str): The base output directory.

    Returns:
        tuple: Absolute paths to the 'r1' and 'r2' subdirectories.
    """
    temp_dir = os.path.join(os.path.abspath(output_dir), "tmp_paired")
    r1_dir = os.path.join(temp_dir, "r1")
    r2_dir = os.path.join(temp_dir, "r2")

    if os.path.exists(temp_dir):
        logger.info(f"{temp_dir} already exists. Removing the existing directory.")
        shutil.rmtree(temp_dir)

    os.makedirs(temp_dir)    
    os.makedirs(r1_dir, exist_ok=True)
    os.makedirs(r2_dir, exist_ok=True)
    
    logger.info(f"Temporary directories created.")
    return r1_dir, r2_dir

def splitting_files(bam_files, r1_dir, r2_dir):
    """
    Splits the provided paired-end BAM files into two sorted and indexed BAM files using samtools.
    For each input BAM file, this function:
      1. Filters the BAM file for Read1 (using flag -f 64), sorts the output, writes the
         result to a new BAM file in r1_dir, and creates an index (.bai) for it.
      2. Filters the BAM file for Read2 (using flag -f 128), sorts the output, writes the
         result to a new BAM file in r2_dir, and creates an index (.bai) for it.

    The splitting rules for paired-end data (ISR and ISF) will later affect the downstream
    interpretation, but the splitting and indexing itself is done by these commands.

    Args:
        bam_files (list): List of input BAM file paths.
        r1_dir (str): Directory to store the sorted and indexed Read1 BAM files.
        r2_dir (str): Directory to store the sorted and indexed Read2 BAM files.

    Returns:
        tuple: Two lists containing the output paths of the sorted BAM files for R1 and for R2.
    """

    print("Prepare files for LoDEI paired-end processing...")
    logger.info("Splitting BAM files into R1 and R2...")

    new_bam_files_r1 = []
    new_bam_files_r2 = []

    for bam in bam_files:
        base = os.path.basename(bam)
        base_name, _ = os.path.splitext(base)

        r1_path = os.path.join(r1_dir, f"{base_name}_R1_sorted.bam")
        r2_path = os.path.join(r2_dir, f"{base_name}_R2_sorted.bam")
        
        logger.info(f"Splitting {bam} into:")
        logger.info(f"  Read1 -> {r1_path}")
        logger.info(f"  Read2 -> {r2_path}")
        
        cmd_r1 = f"samtools view -f 64 -b {bam} | samtools sort -o {r1_path}"
        result_r1 = run(cmd_r1, shell=True)
        if result_r1.returncode != 0:
            logger.error(f"Error splitting and sorting {bam} for Read1.")
            raise RuntimeError(f"Failed to process {bam} for Read1.")

        # Index für R1
        cmd_r1_index = f"samtools index {r1_path}"
        result_r1_index = run(cmd_r1_index, shell=True)
        if result_r1_index.returncode != 0:
            logger.error(f"Error indexing {r1_path}.")
            raise RuntimeError(f"Failed to index {r1_path}.")

        cmd_r2 = f"samtools view -f 128 -b {bam} | samtools sort -o {r2_path}"
        result_r2 = run(cmd_r2, shell=True)
        if result_r2.returncode != 0:
            logger.error(f"Error splitting and sorting {bam} for Read2.")
            raise RuntimeError(f"Failed to process {bam} for Read2")
        
        # Index für R2
        cmd_r2_index = f"samtools index {r2_path}"
        result_r2_index = run(cmd_r2_index, shell=True)
        if result_r2_index.returncode != 0:
            logger.error(f"Error indexing {r2_path}.")
            raise RuntimeError(f"Failed to index {r2_path}.")

        new_bam_files_r1.append(r1_path)
        new_bam_files_r2.append(r2_path)

        logger.info("Splitting completed.")

    return new_bam_files_r1, new_bam_files_r2


def process_single_end(args, files_r1_group1, files_r1_group2, files_r2_group1, files_r2_group2):
    """
    Calls the single-end processing workflow for the split BAM files.

    Args:
        args (dict): Arguments object with required settings.
        files_r1_group1 (list): List of R1 files for group1.
        files_r1_group2 (list): List of R1 files for group2.
        files_r2_group1 (list): List of R2 files for group1.
        files_r2_group2 (list): List of R2 files for group2.
    """
    logger.info("Starting single-end processing...")

    # Determine library type and assign R1/R2 roles
    if args["library"].upper() == "ISF":
        sf_group1 = files_r1_group1
        sr_group1 = files_r2_group1
        sf_group2 = files_r1_group2
        sr_group2 = files_r2_group2
    elif args["library"].upper() == "ISR":
        sr_group1 = files_r1_group1
        sf_group1 = files_r2_group1
        sr_group2 = files_r1_group2
        sf_group2 = files_r2_group2
    else:
        raise ValueError(f"Unsupported library type: {args['library']}")
    
    # Help function to run lodei for single-end processing
    def run_lodei(group1_files, group2_files, library_type, output_base_dir):
    
        # Define the output directory for the specific library type
        output_dir = os.path.join(output_base_dir, "tmp_paired", library_type)
        os.makedirs(output_dir, exist_ok=True)  # Ensure the directory exists

        # Build the command
        cmd = [
            "lodei", "find",
            "--group1", " ".join(group1_files),
            "--group2", " ".join(group2_files),
            "-f", args["fasta"],
            "-g", args["gff"],
            "-o", output_dir,
            "-c", str(args["cores"]),
            "--library", library_type,
            "--min_coverage", str(args["min_coverage"]),
            "--window_size", str(args["window_size"])
        ]

        if args["rm_snps"]:
            cmd.append("--rm_snps")

        # Log the command
        logger.info(f"Running command: {' '.join(cmd)}")
        os.system(" ".join(cmd))  # Execute the command
        # Execute the command
        # try:
        #     result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        #     logger.info(f"Single-end processing ({library_type}) completed successfully:\n{result.stdout}")
        # except subprocess.CalledProcessError as e:
        #     logger.error(f"Error during single-end processing ({library_type}):\n{e.stderr}")
        #     raise

    # Run the workflow for SF
    logger.info("Processing SF library...")
    run_lodei(sf_group1, sf_group2, "SF", args["output"])

    # Run the workflow for SR
    logger.info("Processing SR library...")
    run_lodei(sr_group1, sr_group2, "SR", args["output"])

    logger.info("Single-end processing completed.")

def merge_results(output_dir):
    """
    Merges paired window result files from SF and SR for all nucleotide pairs.

    For each nucleotide pair (e.g., AC, AG, AT, ...), this function finds the corresponding
    window result files from both SF and SR directories, merges them on window coordinates,
    averages the 'wEI' values, and writes the result to a new file in 'merged_windows'.

    Args:
        output_dir (str): Path to the main output directory containing 'tmp_paired/SF/windows'
                          and 'tmp_paired/SR/windows'.

    Output:
        Merged files are saved in {output_dir}/merged_windows/windows_{pair}.txt.
        Only windows present in both SF and SR files are included.
    """

    logger.info("Merging results from R1 and R2...")

    merged_dir = os.path.join(output_dir, "tmp_paired/merged_windows")
    os.makedirs(merged_dir, exist_ok=True)

    sf_pattern = os.path.join(output_dir, "tmp_paired", "SF", "windows", "windows_*.txt")
    sf_files = glob.glob(sf_pattern)

    if not sf_files:
        logger.info(f"No SF files found with pattern: {sf_pattern}")
        return

    for sf_path in sf_files:
        pair = os.path.basename(sf_path).split("_")[1].replace(".txt", "")
        sr_path = os.path.join(output_dir, "tmp_paired", "SR", "windows", f"windows_{pair}.txt")

        if not os.path.exists(sr_path):
            logger.info(f"SR file missing for {pair}, expected: {sr_path}")
            continue

        try:
            df_sf = pd.read_csv(sf_path, sep="\t")
            df_sr = pd.read_csv(sr_path, sep="\t")
        except Exception as e:
            logger.info(f"Error reading {pair}: {e}")
            continue

        keys = ["chrom", "wstart", "wend", "name", "strand"]
        # Merge on the key columns, keep only common rows
        try:
            df_merged = pd.merge(df_sf, df_sr, on=keys, suffixes=('_sf', '_sr'), how='inner')
        except Exception as e:
            logger.info(f"Error merging {pair}: {e}")
            continue

        # Calculate mean for wEI and round to two decimal places
        try:
            df_merged["wEI"] = ((df_merged["wEI_sf"] + df_merged["wEI_sr"]) / 2).round(1)
        except Exception as e:
            logger.info(f"Error calculating mean wEI for {pair}: {e}")
            continue

        # Final column order
        cols = ["chrom", "wstart", "wend", "name", "wEI", "strand"]
        df_final = df_merged[cols]

        out_path = os.path.join(merged_dir, f"lei_windows_{pair}.txt")
        try:
            df_final.to_csv(out_path, sep="\t", index=False)
            logger.info(f"Merged: {out_path} ({len(df_final)} rows)")
        except Exception as e:
            logger.info(f"Error writing for {pair}: {e}")

    logger.info("All pairs merged and saved.")

def clean_workspace(temp_dir):
    """
    Cleans up the temporary workspace:
    - Moves merged_windows one level up and renames it to windows.
    - Removes r1 and r2 subdirectories.
    - Renames tmp_paired to info_SE.
    Args:
        temp_dir (str): Path to the temporary directory (tmp_paired).
    """
    logger.info(f"Cleaning up temporary workspace: {temp_dir}")

    merged_dir = os.path.join(temp_dir, "merged_windows")
    parent_dir = os.path.dirname(temp_dir)
    target_windows = os.path.join(parent_dir, "windows")

    # Move merged_windows to parent and rename to windows
    if os.path.exists(merged_dir):
        if os.path.exists(target_windows):
            shutil.rmtree(target_windows)
        shutil.move(merged_dir, target_windows)

    # Remove r1 and r2 directories
    for sub in ["r1", "r2"]:
        sub_dir = os.path.join(temp_dir, sub)
        if os.path.exists(sub_dir):
            shutil.rmtree(sub_dir)

    # Rename tmp_paired to info_SE
    info_se_dir = os.path.join(parent_dir, "info_SE")
    if os.path.exists(info_se_dir):
        shutil.rmtree(info_se_dir)
    os.rename(temp_dir, info_se_dir)

def process_paired_end(args):
    """
    Main function for processing paired-end data.

    This function is called from find.py when the library parameter indicates paired-end data.
    It performs the following steps:
      1. Creates a temporary folder with subdirectories 'r1' and 'r2'.
      2. Splits the input BAM files into separate BAM files for Read1 and Read2 using samtools.
         The splitting rules are determined by the library type (ISR or ISF).
      3. Calls the existing single-end workflow for the split BAM files.
      4. Merges the results from the R1 and R2 analyses.
      5. Cleans up the temporary workspace.

    Args:
        args (dict): Arguments object with required settings (including 'output', 'group1', 'library', etc.)

    Returns:
    """

    # Create output directory and Log-File
    logger = setup_paired_logger(args["output"])
    logger.info("Paired-end processing started.")
    
    # 1. Create temporary directories.
    r1_dir, r2_dir = create_temp_folder(args["output"]) 

    # 2. Split the BAM files into separate files for R1 and R2.
    bam_files_r1_group1, bam_files_r2_group1 = splitting_files(args["group1"], r1_dir, r2_dir)
    bam_files_r1_group2, bam_files_r2_group2 = splitting_files(args["group2"], r1_dir, r2_dir)

    # 3. Process the split BAM files with the existing single-end workflow.

    process_single_end(args, bam_files_r1_group1, bam_files_r1_group2, bam_files_r2_group1, bam_files_r2_group2)

    # 4. Merge the results from R1 and R2.

    merge_results(args["output"])

    # 5. Clean up the temporary workspace.

    clean_workspace(os.path.join(args["output"], "tmp_paired"))

    logger.info("Paired-end processing completed.")
    return