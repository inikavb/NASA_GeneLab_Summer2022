#!/usr/bin/env python

"""
This is a program for validating GeneLab amplicon processed datasets.
"""
from concurrent.futures import process
from dis import code_info
from enum import Flag
from this import d
from typing import final
from dp_tools.core.check_model import ValidationProtocol, FlagCode, FlagEntry

import os
from re import T
import sys
import argparse
import textwrap
import pandas as pd
import zipfile
import tarfile
import pathlib as path


parser = argparse.ArgumentParser(description="This program validates GeneLab amplicon processed datasets. \
                                             Hard-coded variables that may need to be changed are near the top \
                                             of the script.")


GLDS_ID = input('GLDS ID (e.g. "GLDS-276")')
sample_names_file = input("Path to single-column file with unique sample names: ")
output_prefix = input("Output file prefix if there is one (press 'return' if none): ")
primers_already_trimmed = input("Primers trimmed prior to GeneLab processing? (y/n) ")
if primers_already_trimmed == 'y':
    print("trimmed primers")
    trimmed = True
elif primers_already_trimmed == 'n':
    print("untrimmed primers")
    primers_already_trimmed = False

single_ended = input("data are single-end sequencing? (y/n): ")
if single_ended == 'y':
    print("single ended")
    single_ended = True
elif single_ended == 'n':
    print("not single ended")
    single_ended = False



### hard-coded stuff we might want to change ###
raw_reads_dir = "Raw_Sequence_Data/"
fastqc_dir = "FastQC_Outputs/"
trimmed_reads_dir = "Trimmed_Sequence_Data/"
filtered_reads_dir = "Filtered_Sequence_Data/"
final_outputs_dir = "Final_Outputs/"

if not primers_already_trimmed:

        expected_dirs = [raw_reads_dir, fastqc_dir, trimmed_reads_dir, filtered_reads_dir, final_outputs_dir]

else:

    expected_dirs = [raw_reads_dir, fastqc_dir, filtered_reads_dir, final_outputs_dir]

raw_suffix = "_raw.fastq.gz"
raw_R1_suffix = "_R1_raw.fastq.gz"
raw_R2_suffix = "_R2_raw.fastq.gz"
primer_trimmed_suffix = "_trimmed.fastq.gz"
primer_trimmed_R1_suffix = "_R1_trimmed.fastq.gz"
primer_trimmed_R2_suffix = "_R2_trimmed.fastq.gz"
filtered_suffix = "_filtered.fastq.gz"
filtered_R1_suffix = "_R1_filtered.fastq.gz"
filtered_R2_suffix = "_R2_filtered.fastq.gz"

processing_tar_file = "processing_info.tar"

if output_prefix:
    processing_tar_file = output_prefix.rstrip("-").rstrip("_") + "_processing_info.tar"

expected_trimmed_outputs_or_suffixes = [str(output_prefix) + "cutadapt.log", str(output_prefix) + "trimmed-read-counts.tsv"]
expected_filtered_outputs_or_suffixes = ["filtered-read-counts.tsv"]
expected_final_outputs_or_suffixes = [".fasta", str(output_prefix) + "counts.tsv", str(output_prefix) + "taxonomy.tsv", ".biom.zip", str(output_prefix) + "taxonomy-and-counts.tsv", str(output_prefix) + "read-count-tracking.tsv"]

validation_log = str(GLDS_ID) + "_" + str(output_prefix) + "amplicon-validation.log"

######################################################################

def check_expected_directories(directory : path):
    """ checks expected directories exist """

    if not os.path.isdir(directory):
        code = FlagCode.HALT
        message = "The directory '" + str(directory) + "' was expected but not found."
        #report_failure("The directory '" + str(directory) + "' was expected but not found.")
    else: 
        code = FlagCode.GREEN
        message = "The directory " + str(directory)+ " exists."
    return {"code": code, "message": message}


def read_samples(file_path):
    """ reading unique sample names into list """

    with open(file_path) as f:
        sample_names = f.read().splitlines()

    return(sample_names)


def check_for_file_and_contents(file: path):
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


def check_raw_multiqc_outputs(sample : path):
    """ makes sure all samples' read files are in the multiqc outputs """

    # checking raw
    zip_file = zipfile.ZipFile(raw_multiqc_data_path)
    df = pd.read_csv(zip_file.open("multiqc_general_stats.txt"), sep = "\t", usecols = ["Sample"])
    file_prefixes_in_multiqc = df["Sample"].tolist()

    if not sample in file_prefixes_in_multiqc:
        code = FlagCode.HALT
        message = "The raw multiqc output is missing the expected '" + sample + "' entry."
    else:
        code = FlagCode.GREEN
        message = "The raw multiqc output contains the expected " + sample + " entry."
    return {"code": code, "message": message} 

def check_filtered_multiqc_outputs(sample : path):
    # checking filtered
    zip_file = zipfile.ZipFile(filt_multiqc_data_path)
    df = pd.read_csv(zip_file.open("multiqc_general_stats.txt"), sep = "\t", usecols = ["Sample"])
    file_prefixes_in_multiqc = df["Sample"].tolist()

    if not sample in file_prefixes_in_multiqc:
        code = FlagCode.HALT
        message = "The filtered multiqc output is missing the expected '" + sample + "' entry."
    else:
        code = FlagCode.GREEN
        message = "The filtered multiqc output contains the expected " + sample + " entry."
    return {"code": code, "message": message}



def check_general_fasta_format(file : path):

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



def check_intermediate_log_files_trimmed(entry : path):

    output_files_present = [f for f in os.listdir(trimmed_reads_dir) if os.path.isfile(os.path.join(trimmed_reads_dir, f))]
    if not any(output_file.endswith(entry) for output_file in output_files_present):
        code = FlagCode.HALT
        message = "An output file named or ending with '" + str(entry) + "' was expected but not found in " + str(trimmed_reads_dir) + "."
        
    else: 
        code = FlagCode.GREEN
        message = "The expected output file named or ending with " + str(entry) + " was found."

    return {"code": code, "message": message}

def check_intermediate_log_files_filtered(entry : path):
    
    output_files_present = [f for f in os.listdir(filtered_reads_dir) if os.path.isfile(os.path.join(filtered_reads_dir, f))]
    if not any(output_file.endswith(entry) for output_file in output_files_present):
        code = FlagCode.HALT
        message = "An output file named or ending with '" + str(entry) + "' was expected but not found in " + str(filtered_reads_dir) + "."
    else:
        code = FlagCode.GREEN
        message = "The expected output file named or ending with " + str(entry) + " was found."

    return {"code": code, "message": message}


def check_final_outputs(entry : path):
    """ makes sure outputs exist and checks formatting """

    if not any(output_file.endswith(entry) for output_file in output_files_present):
        code = FlagCode.HALT
        message = "An output file named or ending with '" + str(entry) + "' was expected but not found in " + str(final_outputs_dir) + "."
    else:
        code = FlagCode.GREEN
        message = "The expected output file named or ending with '" + str(entry) + "' was found in " + str(final_outputs_dir) + "."

    return {"code": code, "message": message}


def check_for_processing_tar(entries:path):
    """ this just makes sure a processing tar exists and at least has the Snakefile, as its contents can vary quite a bit """

    target_substring = str(output_prefix).rstrip("-") + "/Snakefile"

    base_target = "/Snakefile"

    if not any(target_substring.lower() in string.lower() for string in entries):

        if not any(base_target.lower() in string.lower() for string in entries):
            code = FlagCode.HALT
            message = "The '" + processing_tar_file + "' does not have a 'Snakefile' as expected."
        
    else:
        code = FlagCode.GREEN
        message = "The '" + processing_tar_file + "' contains the expected 'Snakefile.'"

    return {"code": code, "message": message}

######################################################################
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
            primer_trimmed_suffixes = [primer_trimmed_R1_suffix, primer_trimmed_R2_suffix]
            filtered_R_suffixes = [filtered_R1_suffix,filtered_R2_suffix]
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
                if primers_already_trimmed == False:
                    for primer_suffix in primer_trimmed_suffixes:
                        with vp.payload(
                        payloads=[
                        {
                            "file" : trimmed_reads_dir + sample + primer_suffix
                        }
                        ]
                    ):
                            vp.add(check_for_file_and_contents) 
                for filtered_R_suffix in filtered_R_suffixes:
                    with vp.payload(
                        payloads=[
                        {
                            "file" : filtered_reads_dir + sample + filtered_R_suffix
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
                if primers_already_trimmed == False:
                    for primer_suffix in primer_trimmed_suffixes:
                        with vp.payload(
                        payloads=[
                        {
                            "file" : trimmed_reads_dir + sample + primer_suffix
                        }
                        ]
                    ):
                            vp.add(check_for_file_and_contents) 
                for filtered_R_suffix in filtered_R_suffixes:
                    with vp.payload(
                    payloads=[
                    {
                        "file" : filtered_reads_dir + sample + filtered_R_suffix
                    }
                    ]
                ):
                        vp.add(check_for_file_and_contents)

# checking raw multiqc data
raw_multiqc_data_path = fastqc_dir + str(output_prefix) + "raw_multiqc_data.zip"
with vp.component_start(
        name="raw multiqc zip file check",
        description="make sure raw multiqc zip file exists and holds content"
    ):
        with vp.payload(
            payloads=[
                {
                    "file" : raw_multiqc_data_path
                }
            ]
        ):
            vp.add(check_for_file_and_contents)

for sample in sample_names:
    with vp.component_start(
        name="raw multiqc outputs check",
        description="make sure raw multiqc outputs exist and hold content"
    ):
        if not single_ended:
            R1_suffix = raw_R1_suffix.split(".")[0]
            R2_suffix = raw_R2_suffix.split(".")[0]
            raw_R_suffixes = [R1_suffix,R2_suffix]
            for raw_R_suffix in raw_R_suffixes:
                with vp.payload(
                    payloads=[
                        {
                            "sample":sample + raw_R_suffix
                        }
                    ]
                ):
                    vp.add(check_raw_multiqc_outputs)
            
        else: 
            suffix = raw_suffix.split(".")[0]
            with vp.payload(
                payloads=[
                    {
                        "sample" : sample + suffix
                    }
                ]
            ):
                vp.add(check_raw_multiqc_outputs)
            
filt_multiqc_data_path = fastqc_dir + str(output_prefix) + "filtered_multiqc_data.zip"
with vp.component_start(
        name="filtered multiqc zip file check",
        description="make sure filtered multiqc zip file exists and holds content"
    ):
        with vp.payload(
            payloads=[
                {
                    "file" : filt_multiqc_data_path
                }
            ]
        ):
            vp.add(check_for_file_and_contents)

for sample in sample_names:
    with vp.component_start(
        name="filtered multiqc outputs check",
        description="make sure filtered multiqc outputs exist and hold content"
    ):
        if not single_ended:
            R1_suffix = filtered_R1_suffix.split(".")[0]
            R2_suffix = filtered_R2_suffix.split(".")[0]
            filtered_R_suffixes = [R1_suffix,R2_suffix]
            for filtered_R_suffix in filtered_R_suffixes:
                with vp.payload(
                    payloads=[
                        {
                            "sample" : sample + filtered_R_suffix
                        }
                    ]
                ):
                    vp.add(check_filtered_multiqc_outputs)
            
        else: 
            suffix = filtered_suffix.split(".")[0]
            with vp.payload(
                payloads=[
                    {
                        "sample" : sample + suffix
                    }
                ]
            ):
                vp.add(check_filtered_multiqc_outputs)


## trimmed if needed
if not primers_already_trimmed:
    
    with vp.component_start(
        name="check trimmed intermediate log files",
        description="check for the expected trimmed intermediate log files"
    ):

        for entry in expected_trimmed_outputs_or_suffixes:
            with vp.payload(
                payloads=[
                    {
                        "entry" : entry
                    }
                ]
            ):
                vp.add(check_intermediate_log_files_trimmed)
## filtered
with vp.component_start(
        name="check filtered intermediate log files",
        description="check for the expected filtered intermediate log files"
    ):
    for entry in expected_filtered_outputs_or_suffixes:
        with vp.payload(
            payloads=[
                    {
                        "entry" : entry
                    }
                ]
        ):
            vp.add(check_intermediate_log_files_filtered)



# getting list of files in output dir
output_files_present = [f for f in os.listdir(final_outputs_dir) if os.path.isfile(os.path.join(final_outputs_dir, f))]

# making sure none of them are empty
with vp.component_start(
    name="check final output files",
    description="Checking final output files for existence and content"
):
    for output_file in output_files_present:
        with vp.payload(
            payloads=[
                {
                    "file" : final_outputs_dir + output_file
                }
            ]
        ):
            vp.add(check_for_file_and_contents)

with vp.component_start(
    name="check final output files",
    description="Checking final output files for existence and content"
):
    for entry in expected_final_outputs_or_suffixes:
        with vp.payload(
            payloads=[
            {
                "entry" : entry
            }
            ]
        ):
            vp.add(check_final_outputs)

# checking general fasta format is met
fasta_files_in_output_dir = [output_file for output_file in output_files_present if output_file.endswith(".fasta")]
with vp.component_start(
    name="check final outputs",
    description="Checking final output files for proper formatting "
):
    for fasta_file in fasta_files_in_output_dir:
        with vp.payload(
            payloads=[
                {
                    "file":final_outputs_dir + fasta_file
                }
            ]
        ): 
            vp.add(check_general_fasta_format)

with vp.component_start(
    name="check processing tar file",
    description="check processing tar file for existence and content"
):
    with vp.payload(
        payloads=[
            {
                "file":processing_tar_file
            }
        ]
    ):
        vp.add(check_for_file_and_contents)

with tarfile.open(processing_tar_file) as tar_obj:
    entries = tar_obj.getnames()
with vp.component_start(
    name="check processing tar file",
    description="check processing tar file for expected 'Snakefile'"
):
    with vp.payload(
        payloads=[
            {
                "entries":entries
            }
        ]
    ):
        vp.add(check_for_processing_tar)

##########################################
#  Running Protocol
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
