
#!/usr/bin/env python

"""
This is a program for validating GeneLab Illumina metageonmics processed datasets.
"""

from enum import Flag
from http.client import PAYMENT_REQUIRED
from this import d
from dp_tools.core.check_model import ValidationProtocol, FlagCode, FlagEntry

import os
import sys
import argparse
import textwrap
import pandas as pd
import zipfile
import tarfile
import pathlib as path


parser = argparse.ArgumentParser(description="This program validates GeneLab Illumina metagenomics processed datasets. \
                                             Hard-coded variables that may need to be changed are near the top \
                                             of the script.")

GLDS_ID = input("GLDS ID (e.g. 'GLDS-276')")
sample_names_file = input("Path to file containing single-column of unique sample names: ")
additional_prefix = input("Add any expected additional filename prefix that was added to the files that describe multiple samples (default: \"\") Press 'return' if none: ")
single_ended = input("single ended? (y/n): ")
if single_ended == 'y':
    print("single ended")
    single_ended = True
elif single_ended == 'n':
    print("not single ended")
    single_ended = False
else: 
    single_ended = input("invalid input. single ended? (y/n): ")
    
### hard-coded stuff we might want to change ###
raw_reads_dir = "Raw_Sequence_Data/"
fastqc_dir = "FastQC_Outputs/"
filtered_reads_dir = "Filtered_Sequence_Data/"
assembly_based_dir = "Assembly-based_Processing/"
assemblies_dir = assembly_based_dir + "assemblies/"
genes_dir = assembly_based_dir + "predicted-genes/"
annotations_and_tax_dir = assembly_based_dir + "annotations-and-taxonomy/"
mapping_dir = assembly_based_dir + "read-mapping/"
combined_output_dir = assembly_based_dir + "combined-outputs/"
bins_dir = assembly_based_dir + "bins/"
MAGs_dir = assembly_based_dir + "MAGs/"
read_based_dir = "Read-based_Processing/"
logs_dir = additional_prefix + "processing_info/logs/"

expected_dirs = [raw_reads_dir, fastqc_dir, filtered_reads_dir, assembly_based_dir,
                 assemblies_dir, genes_dir, annotations_and_tax_dir, mapping_dir,
                 combined_output_dir, bins_dir, MAGs_dir, read_based_dir]

assembly_based_dirs = [assemblies_dir, genes_dir, annotations_and_tax_dir, mapping_dir,
                       combined_output_dir, bins_dir, MAGs_dir]

raw_suffix = "_HRremoved_raw.fastq.gz"
raw_R1_suffix = "_R1_HRremoved_raw.fastq.gz"
raw_R2_suffix = "_R2_HRremoved_raw.fastq.gz"
filtered_suffix = "_filtered.fastq.gz"
filtered_R1_suffix = "_R1_filtered.fastq.gz"
filtered_R2_suffix = "_R2_filtered.fastq.gz"

assembly_suffix = "-assembly.fasta"
mapping_dir_suffixes_all_have = [".bam", "-metabat-assembly-depth.tsv"]
mapping_info_suffix = "-mapping-info.txt"


expected_assembly_combined_outputs = [str(additional_prefix) + "Combined-contig-level-taxonomy-coverages-CPM.tsv",
                                      str(additional_prefix) + "Combined-gene-level-KO-function-coverages-CPM.tsv",
                                      str(additional_prefix) + "Combined-gene-level-taxonomy-coverages-CPM.tsv",
                                      str(additional_prefix) + "Combined-contig-level-taxonomy-coverages.tsv",
                                      str(additional_prefix) + "Combined-gene-level-KO-function-coverages.tsv",
                                      str(additional_prefix) + "Combined-gene-level-taxonomy-coverages.tsv"]

assembly_based_overview_table = assembly_based_dir + str(additional_prefix) + "Assembly-based-processing-overview.tsv"

expected_read_based_outputs = [str(additional_prefix) + "Gene-families-KO-cpm.tsv", 
                               str(additional_prefix) + "Gene-families-cpm.tsv", 
                               str(additional_prefix) + "Gene-families-grouped-by-taxa.tsv",
                               str(additional_prefix) + "Gene-families.tsv", 
                               str(additional_prefix) + "Metaphlan-taxonomy.tsv", 
                               str(additional_prefix) + "Pathway-abundances-cpm.tsv",
                               str(additional_prefix) + "Pathway-abundances-grouped-by-taxa.tsv", 
                               str(additional_prefix) + "Pathway-abundances.tsv",
                               str(additional_prefix) + "Pathway-coverages-grouped-by-taxa.tsv", 
                               str(additional_prefix) + "Pathway-coverages.tsv"]

expected_final_outputs_or_suffixes = [".fasta", "counts.tsv", "taxonomy.tsv", ".biom.zip", "taxonomy-and-counts.tsv", "read-count-tracking.tsv"]

processing_tar_file = str(additional_prefix) + "processing_info.tar"

expected_tar_contents = ["/Snakefile", "/config.yaml", "/envs", "/logs", "/scripts", "/unique-sample-IDs.txt"]

expected_log_file_suffixes = ["-CAT.log", "-assembly.log", "-bam-summarize-and-metabat.log", "-bowtie2-build.log", 
                              "-bbduk.log", "-kofamscan.log", "-pileup.log", "-prodigal.log", "-humann3-run.log"]

predicted_gene_file_suffixes = ["-genes.faa", "-genes.gff", "-genes.fasta"]
gene_fasta_suffixes = ["-genes.faa", "-genes.fasta"]
annotations_suffixes = ["-gene-coverage-annotation-and-tax.tsv", "-contig-coverage-and-tax.tsv"]

if additional_prefix == "":
    validation_log = str(GLDS_ID) + "-metagenomics-validation.log"
else:
    validation_log = str(GLDS_ID) + "_" + str(additional_prefix) + "metagenomics-validation.log"

################################################################################


################################################################################

### functions ###
def check_expected_directories(directory : path) -> FlagEntry :
    """ checks expected directories exist """

    if not os.path.isdir(directory):
        code = FlagCode.HALT
        message = "The directory '" + str(directory) + "' was expected but not found."
        #report_failure("The directory '" + str(directory) + "' was expected but not found.")
    else: 
        code = FlagCode.GREEN
        message = "The directory " + str(directory)+ " exists."
    return {"code": code, "message": message}


def check_for_file_and_contents(file : path):
    """ used by check_fastq_files and check_final_outputs functions """
    if os.path.exists(file): # path exists
        # now checking for content
        if not os.path.getsize(file) > 0: # file is empty
            code = FlagCode.HALT
            message = "The file " + str(file) + " exists, but the file is empty."
        else:
            code = FlagCode.GREEN # path exists and holds content
            message = "The file " + str(file) + " exists and holds content."
    else:
        code = FlagCode.HALT
        message = "The expected file '" + str(file) + "' does not exist."
    return {"code": code, "message": message} 


def check_raw_multiqc_outputs(file : path):
    """ makes sure all samples' read files are in the multiqc outputs """
    zip_file = zipfile.ZipFile(fastqc_dir + str(additional_prefix) + "raw_multiqc_data.zip")
    df = pd.read_csv(zip_file.open("multiqc_general_stats.txt"), sep = "\t", usecols = ["Sample"])
    file_prefixes_in_multiqc = df["Sample"].tolist()

    if not file in file_prefixes_in_multiqc:
        code = FlagCode.HALT
        message = "The raw multiqc output is missing the expected '" + file + "' entry."
    else:
        code = FlagCode.GREEN
        message = "The raw multiqc output contains the expected '" + file + "' entry."

    return {"code": code, "message": message} 


def check_filtered_multiqc_outputs(file : path):   
# checking filtered
    zip_file = zipfile.ZipFile(fastqc_dir + str(additional_prefix) + "filtered_multiqc_data.zip")
    df = pd.read_csv(zip_file.open("multiqc_general_stats.txt"), sep = "\t", usecols = ["Sample"])
    file_prefixes_in_multiqc = df["Sample"].tolist()

    if not file in file_prefixes_in_multiqc:
        code = FlagCode.HALT
        message = "The raw multiqc output is missing the expected '" + file + "' entry."
    else:
        code = FlagCode.GREEN
        message = "The raw multiqc output contains the expected '" + file + "' entry."

    return {"code": code, "message": message} 


def check_log_files(entry : path):

    ## filtered
    output_files_present = [f for f in os.listdir(filtered_reads_dir) if os.path.isfile(os.path.join(filtered_reads_dir, f))]

    if not any(output_file.endswith(entry) for output_file in output_files_present):
        code = FlagCode.HALT
        message = "An output file named or ending with '" + str(entry) + "' was expected but not found in " + str(filtered_reads_dir) + "."
    else:
        code = FlagCode.GREEN
        message = "The expected output file named or ending with " + str(entry) + " was found in " + str(filtered_reads_dir) + "."
    return {"code": code, "message": message} 
    

def check_general_fasta_format(file : path):
    
    if not os.path.getsize(file) > 0:
        code = FlagCode.HALT
        message = "The fasta file '" + str(file) + "' is empty but isn't expected to be."

    line_num = 0
    num_headers = 0
    num_seqs = 0

    with open(file) as in_file:

        for line in in_file:

            # keeping track of current line for reporting any problems
            line_num += 1

            if line.strip().startswith(">"):
                num_headers += 1
            else:
                num_seqs += 1

            if num_headers != num_seqs + 1 and num_headers != num_seqs:
                code = FlagCode.HALT
                message = "Fasta file '" + str(file) + "' does not seem to be formatted properly. Problem detected at line " + str(line_num) + "."
            else:
                code = FlagCode.GREEN
                message = "Fasta file " + str(file) + " is formatted properly."
    return {"code": code, "message": message}


def check_assembly_based_file(sample : path, file: path, failed_assemblies_list: path ,assembly = True):
    
    if not os.path.exists(file):
        code = FlagCode.HALT
        message = "The expected file '" + str(file) + "' does not exist."
    else: 
        if not os.path.getsize(file) > 0:
            # a sample can have no genes called even if the assembly produced contigs, so this is only throwing a warning if we are checking an assembly here
            if sample not in failed_assemblies_list and assembly == True:
                code = FlagCode.HALT
                message = "The file '" + str(file) + "' is empty, but that sample isn't noted in the 'Failed-assemblies.tsv' file as it should be if the assembly failed."
            elif sample in failed_assemblies_list and assembly == True:
                code = FlagCode.YELLOW
                message = "The file '" + str(file) + "' is exists and is noted in the 'Failed-assemblies.tsv' file as it should be if the assembly failed."
        else: 
            code = FlagCode.GREEN
            message = "The file '" + str(file) + "' exists and contains content."
    return {"code": code, "message": message}


    
def check_assembly_based_genes_file(sample : path, file: path, failed_assemblies_list :path, assembly = True):
    """ 
    separate function for working with expected output genes files, to handle
    cases where assemblies can succeed while there are still no gene calls

    just checks the file isn't empty if it exists
    """
    
    if os.path.exists(file) and sample not in failed_assemblies_list:

        if not os.path.getsize(file) > 0:
            code = FlagCode.HALT
            message = "The expected file '" + str(file) + "' exists, but appears to be empty when it shouldn't be."
        else: 
            code = FlagCode.GREEN
            message = "The expected file '" + str(file) + "' exists and holds content."
    else:
        code = FlagCode.HALT
        message = "The expected file '" + str(file) + "' is not in the failed assemblies list."
    return {"code": code, "message": message}


def check_assembly_summary(assembly_summary_path:path):
    # making sure assembly summary file is there

    if not os.path.exists(assembly_summary_path):
        code =  FlagCode.HALT
        message = "The assembly summary file, " + str(assembly_summary_path) + ", is expected but was not found."
    else: 
        code = FlagCode.GREEN
        message = "The assembly summary file, " + str(assembly_summary_path) + ", was found."
    return {"code": code, "message": message}

    
def check_bins_dir(output_files_present:path):
    """ makes sure outputs exist and checks formatting """

    ## bins_dir ##

    if output_files_present:

        output_fasta_bins = [filename for filename in output_files_present if filename.endswith(".fasta")]

        # checking for contents (checking fasta format not straightforward when there are softwraps, but don't want to remove them on these due to large contigs)
        for bin_file in output_fasta_bins:

            curr_file_path = bins_dir + bin_file

            if not os.path.getsize(curr_file_path) > 0:
                code=FlagCode.HALT
                message="The file '" + str(curr_file_path) + "' is empty, but shouldn't be there if that's the case."
            else:
                code=FlagCode.GREEN
                message="The file '" + str(curr_file_path) + "' holds content."

        # making sure summary table is there if there are any bins
        if len(output_fasta_bins) > 0:

            bins_summary_path = bins_dir + additional_prefix + "bins-overview.tsv"

            if not os.path.exists(bins_summary_path):
                code=FlagCode.HALT
                message="The bins summary file, " + str(bins_summary_path) + ", is expected but was not found."
            else:
                code=FlagCode.GREEN
                message="The file '" + str(bins_summary_path) + "' holds content."
    return {"code": code, "message": message}

def check_mags_dir(output_files_present:path):
    ## MAGs_dir ##

    if output_files_present:

        output_fasta_MAGs = [filename for filename in output_files_present if filename.endswith(".fasta")]

        # checking for contents (checking fasta format not straightforward when there are softwraps, but don't want to remove them on these due to large contigs)
        for MAG_file in output_fasta_MAGs:

            curr_file_path = MAGs_dir + MAG_file

            if not os.path.getsize(curr_file_path) > 0:
                code=FlagCode.HALT
                message= "The file '" + str(curr_file_path) + "' is empty, but shouldn't be there if that's the case."
            else:
                code=FlagCode.GREEN
                message="The file '" + str(curr_file_path) + "' holds content."

        # making sure summary table is there if there are any bins
        output_fasta_bins = [filename for filename in output_files_present if filename.endswith(".fasta")]
        if len(output_fasta_bins) > 0:

            MAGs_summary_path = MAGs_dir + additional_prefix + "MAGs-overview.tsv"

            if not os.path.exists(MAGs_summary_path):
                code = FlagCode.HALT
                message = "The MAGs summary file, " + str(MAGs_summary_path) + ", is expected but was not found."
            else:
                code=FlagCode.GREEN
                message="The file '" + str(MAGs_summary_path) + "' holds content."

    return {"code": code, "message": message}


def check_assembly_based_overview_table(expected_samples: path, overview_table_path:path):
    """ makes sure the output table exists and all input samples are in it """

    # now making sure all samples are in there
    # reading in table and getting sample IDs in list
    overview_tab = pd.read_csv(overview_table_path, sep = "\t")
    samples_in_tab = overview_tab['Sample_ID'].tolist()

    missing_sample_IDs = []

    for sample in expected_samples:
        if sample not in samples_in_tab:
            missing_sample_IDs.append(sample)

    if len(missing_sample_IDs) > 0:
        code=FlagCode.HALT
        message="The assembly overview table, '" + overview_table_path + "', doesn't have all the samples expected to be there."
    else:
        code=FlagCode.GREEN
        message="The assembly overview table, '" + overview_table_path + "', contains all the expected samples."
    return {"code": code, "message": message}
    


def check_processing_tar(entries:path,item:path):
    """ this makes sure a processing tar exists and contains what we expect """
    if entries[0] + item not in entries:
        code=FlagCode.HALT
        message = "The '" + str(processing_tar_file) + "' does not have '" + str(item) + "' as expected."
    else:
        code=FlagCode.GREEN
        message="The '" + str(processing_tar_file) + "' contains the expected '" + str(item) + "'."
    return {"code": code, "message": message}

def check_tar_log_files(samples:path):
    # checking log files
    for sample in samples:

        for suffix in expected_log_file_suffixes:

            target_log = logs_dir + sample + suffix

            if target_log not in entries:
                code= FlagCode.HALT
                message = "The '" + str(processing_tar_file) + "' does not have the '" + str(target_log) + "' log file as expected."
            else: 
                code = FlagCode.GREEN
                message = "The '" + str(processing_tar_file) + "' contains the expected '" + str(target_log) + "' log file."
    return {"code": code, "message": message}

################################################

vp = ValidationProtocol()


with vp.component_start(
    name="Directory Check",
    description="Check for directory existence from expected directory list",
):
        for directory in expected_dirs:
            with vp.component_start(
                name="directory",
                description="make sure the directory exists",
            ):
                with vp.payload(payloads=[{"directory": directory}]):
                    vp.add(check_expected_directories)

with open(sample_names_file) as f:
    sample_names = f.read().splitlines()


for sample in sample_names:
    with vp.component_start(
        name="File Check",
        description="make sure fastq and multiqc files exist and hold content",
    ):
        with vp.component_start(
            name="fastq file",
            description="make sure file exists and holds content",
        ):
            rawR_suffixes = [raw_R1_suffix,raw_R2_suffix]
            if single_ended == False:
                for R_suffix in rawR_suffixes:
                    with vp.payload(
                        payloads=[
                        {
                            "file" : raw_reads_dir + sample + R_suffix
                        }
                        ]
                    ):
                        vp.add(check_for_file_and_contents)
            else: 
                with vp.payload(
                        payloads=[
                        {
                            "file" : raw_reads_dir + sample + raw_suffix
                        }
                        ]
                    ):
                        vp.add(check_for_file_and_contents)

        with vp.component_start(
            name="raw multiqc file check",
            description="make sure raw multiqc files exist and hold content"
        ):
            R1_suffix = raw_R1_suffix.split(".")[0]
            R2_suffix = raw_R2_suffix.split(".")[0]
            rawR_suffixes = [R1_suffix,R2_suffix]
            if single_ended == False:
                suffix = raw_suffix.split(".")[0]
                for R_suffix in rawR_suffixes:
                    with vp.payload(
                        payloads=[
                            {
                                "file" : sample + R_suffix
                            }
                        ]
                    ):
                        vp.add(check_raw_multiqc_outputs)
            else:
                with vp.payload(
                    payloads=[
                        {
                            "file" : sample + raw_suffix
                        }
                    ]
                ):
                    vp.add(check_raw_multiqc_outputs)

        with vp.component_start(
            name="raw multiqc file check",
            description="make sure raw multiqc files exist and hold content"
        ):        
            R1_suffix = filtered_R1_suffix.split(".")[0]
            R2_suffix = filtered_R2_suffix.split(".")[0]
            filteredR_suffixes = [R1_suffix, R2_suffix]
            if single_ended == False:
                for filtered_R_suffix in filteredR_suffixes:
                        with vp.payload(
                            payloads=[
                            {
                                "file" : sample + filtered_R_suffix
                            }
                        ]
                    ):
                            vp.add(check_filtered_multiqc_outputs)
            else:
                with vp.payload(
                    payloads=[
                        {
                            "file" : sample + filtered_R_suffix
                        }
                    ]
                ):
                    vp.add(check_filtered_multiqc_outputs)

    failed_assemblies_list = []
    if os.path.exists(assemblies_dir + additional_prefix + "Failed-assemblies.tsv"):
        with open(assemblies_dir + additional_prefix + "Failed-assemblies.tsv") as failed_assemblies:
            for line in failed_assemblies:
                failed_assemblies_list.append(line.strip().split("\t")[0])
    successful_assemblies = list(set(sample_names) - set(failed_assemblies_list))
    curr_fasta_path = assemblies_dir + sample + assembly_suffix

    if sample not in failed_assemblies_list:
        file = assemblies_dir + sample + assembly_suffix
        with vp.component_start(
            name= "Assembly based file and fasta format check",
            description= "Makes sure output exists and are formatted properly"
        ):
            with vp.payload(
                payloads=[
                    {
                        "sample": sample, "file": curr_fasta_path, "failed_assemblies_list": failed_assemblies_list
                    }
                ]
            ):
                vp.add(check_assembly_based_file)

        with vp.component_start(
            name = "fasta format check",
            description= "Checking the format of fasta files"
        ):
            with vp.payload(
                payloads=[
                    {
                        "file" : curr_fasta_path
                    }
                ]
            ):
                vp.add(check_general_fasta_format)
    
    
    ## genes_dir ##

    for suffix in predicted_gene_file_suffixes:
        curr_file_path = genes_dir + sample + suffix
        with vp.component_start(
            name= "check assembly based genes files",
            description="Checking to see if assembly based gene files exist"
        ):
            with vp.payload(
                payloads=[
                    {
                        "sample":sample, "file": curr_file_path, "failed_assemblies_list": failed_assemblies_list, "assembly":False
                    }
                ]
            ):
# if any assemblies failed, these files won't exist for that assembly (they also may not exist if an assembly produced contigs too but not genes were called)
                vp.add(check_assembly_based_genes_file)
    
    for suffix in gene_fasta_suffixes:
        curr_file_path = genes_dir + sample + suffix
        if os.path.exists(curr_file_path) and os.path.getsize(curr_fasta_path) > 0:
            with vp.component_start(
            name = "fasta format check",
            description= "Checking the format of fasta files"
        ):
                with vp.payload(
                    payloads=[
                        {
                            "file" : curr_file_path
                        }
                    ]
                ):
                    vp.add(check_general_fasta_format)
    
    ## annotations_and_tax_dir ##
    for suffix in annotations_suffixes:
        curr_file_path = annotations_and_tax_dir + sample + suffix
        with vp.component_start(
            name="annotations and tax directory check",
            description="Checking for file existence and content"
        ):
            with vp.payload(
                payloads=[
                    {
                        "file":curr_file_path
                    }
                ]
            ):
                vp.add(check_for_file_and_contents)

## mapping_dir ##
    for suffix in mapping_dir_suffixes_all_have:
        curr_file_path = mapping_dir + sample + suffix

        # checking the file is present and not empty unless it is noted in the Failed-assemblies file
        if sample not in failed_assemblies_list:
            with vp.component_start(
                name="map directory check",
                description="checking map directory"
            ):
                with vp.payload(
                    payloads=[
                        {
                            "sample": sample, "file":curr_file_path, "failed_assemblies_list": failed_assemblies_list
                        }
                    ]
                ):
                    vp.add(check_assembly_based_file)
            

        # checking for mapping-info file for those that should have it
    if sample not in failed_assemblies_list:

        curr_file_path = mapping_dir + sample + mapping_info_suffix
        with vp.component_start(
                name="map directory check",
                description="checking map directory"
        ):
            with vp.payload(
                payloads=[
                    {
                        "sample":sample, "file":curr_file_path, "failed_assemblies_list":failed_assemblies_list
                    }
                ]
            ):
                vp.add(check_assembly_based_file)



assembly_summary_path = assemblies_dir + additional_prefix + "assembly-summaries.tsv"
with vp.component_start(
    name="assembly summary",
    description="Checking the assembly summary path"
):
    with vp.payload(
        payloads=[
            {
                "assembly_summary_path": assembly_summary_path
            }
        ]
    ):
        vp.add(check_assembly_summary) 


## combined_output_dir ##
for filename in expected_assembly_combined_outputs:

    curr_file_path = combined_output_dir + filename
    with vp.component_start(
        name="combined outputs directory check",
        description="Check if combined outputs directory exists and holds content"
    ):
        with vp.payload(
            payloads=[
                {
                    "file":curr_file_path
                }
            ]
        ):
            vp.add(check_for_file_and_contents)

# only if there were bins recovered
output_files_present = [f for f in os.listdir(bins_dir) if os.path.isfile(os.path.join(bins_dir, f))]
with vp.component_start(
    name="Check bins directory",
    description="Checking bins directory for existence and content"
):
    with vp.payload(
        payloads=[
            {
                "output_files_present":output_files_present
            }
        ]
    ):
        vp.add(check_bins_dir)

# only if there were MAGs recovered
output_files_present = [f for f in os.listdir(MAGs_dir) if os.path.isfile(os.path.join(MAGs_dir, f))]
with vp.component_start(
    name="Check mags directory",
    description="Checking mags directory for existence and content"
):
    with vp.payload(
        payloads=[
            {
                "output_files_present":output_files_present
            }
        ]
    ):
        vp.add(check_mags_dir)

# first making sure it exists and is not empty
with vp.component_start(
    name="Check assembly based overview table directory",
    description="Checking assembly based directory for existence and content"
):
    with vp.payload(
        payloads=[
            {
                "file":assembly_based_overview_table
            }
        ]
    ):
        vp.add(check_for_file_and_contents)

with vp.component_start(
    name="Check assembly based overview table",
    description="Checking assembly based overview table for expected samples"
):
    with vp.payload(
        payloads=[
            {
                "expected_samples":sample_names,"overview_table_path":assembly_based_overview_table
            }
        ]
    ):
        vp.add(check_assembly_based_overview_table)

for file in expected_read_based_outputs:
    with vp.component_start(
        name="Check read based outputs",
        description="Checking read based outputs for existence and content"
    ):
        with vp.payload(
            payloads=[
                {
                    "file":read_based_dir + file
                }
            ]
        ):
            vp.add(check_for_file_and_contents)

with vp.component_start(
    name="Check processing tar file",
    description="Checking processing tar file for existence and content"
):
    with vp.payload(
        payloads=[
            {
                "file":processing_tar_file
            }
        ]
    ):
        vp.add(check_for_file_and_contents)

with vp.component_start(
    name="Check processing tar file",
    description="Check tar log files for expected target log files"
):
    with vp.payload(
        payloads=[
            {
                "samples":sample_names
            }
        ]
    ):
        vp.add(check_tar_log_files)



with tarfile.open(processing_tar_file) as tar_obj:
    entries = tar_obj.getnames()

for item in expected_tar_contents:
    with vp.component_start(
        name="Check processing tar file",
        description="Checking processing tar file for expected tar contents"
    ):
        with vp.payload(
            payloads=[
                {
                    "entries":entries,"item":item
                }
            ]
        ):
            vp.add(check_processing_tar)

            

        

for entry in expected_final_outputs_or_suffixes:
    with vp.component_start(
        name= "Log Files Check",
        description= "Checking log files"
    ):
        with vp.payload(
            payloads=[
                {
            "entry": entry
                }
            ]
        ):
            vp.add(check_log_files)

            
 
'''

                 
'''
   
#########################################
# Running Protocol
#########################################

# Now printing the queued checks
print(vp.queued_checks())

print("********************************")
# Running the checks
vp.run()
#vp.run(skip_components = ["Tires"])

# And showing the results
print(vp.report()["flag_table"])

print("Individual Flag Row")
print(vp.report()["flag_table"].iloc[5].to_dict())


#########################################
# Exporting Protocol Results to File
#########################################

vp.report()["flag_table"].to_csv("VV_log.tsv", sep = "\t")
