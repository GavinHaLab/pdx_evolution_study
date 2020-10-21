#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1 Mb Discordance Plotting Ichor Minimal Subclone Discordance
July 22nd, 2020
Anna Hoge
Ha Lab
Fred Hutchinson Cancer Research Center
"""

#imports
import os
import matplotlib.pyplot as plt
import re
import statistics
import glob
import math


#working directory
#os.chdir("/Volumes/fh")
os.chdir("/Volumes/ha_g")
#os.chdir("/Volumes/ha_g-1")
#os.chdir("/Volumes/ha_g-2")


#settings
"""
--SAMPLE_SIZE_THRESHOLDS - A list of integers to use as the minimum number of
fully analyzed pairings a cohort must have to be included in the plots.
--ADDITIONAL_DISCORDANCE_THRESHOLD - Two samples must have rounded ploidy-
adjusted copy number at a bin for that bin to be discordant, as well as having
a difference in unrounded ploidy-adjusted copy number of at least this number.
--USE_SUBCLONE_FILTER - Boolean indicating whether subclonal copy number calls
should be treated differently than clonal copy number calls.  If True, subclonal
copy number calls with cellular prevalence below SUCLONE_THRESHOLD will not
be discordant from either the loss or gain they have character of or from neutral
copy number.
--SUBCLONE_THRESHOLD - When USE_SUBCLONE_FILTER is set to True, subclones
with cellular prevalence above this number will not be treated differently from
clonal calls.
--USE_ROUNDED_CN - Boolean indicating whether to use ichorCNA's Corrected_Copy_Number
integer output (if True) or ichorCNA's logR_Copy_Number float output (if False)
as the copy number to be ploidy corrected.
--ICHOR_RESULTS_FILEPATH - Filepath to ichorCNA's results folder.
--USE_ALL_PAIRINGS - Boolean indicating whether all possible pairings between
a PT and its PDXs should be used, provided the PDXs are not of the same passage.
False would indicate only consecutive pairings should be used (e.g. given PT,
P1, and P3, pair PT-P1 and P1-P3 only, excluding PT-P3).
--PLOT_OUTPUT_PATH - Path to folder where output plots and files should go.
--MEDIANS_FILENAME - Name of output file listing cohort medians.
--PAIRING_DISCORDANCE_VALUES_FILENAME - Name of output file listing discordance
values for each analyzed pairing.
--WHITE_SAMPLES_FILENAME - Name of output file listing samples which had <5%
of bins with ploidy-adjusted CN != 2 and/or < 5% tumor fraction and were thus
not included in the discordance analysis.
"""
SAMPLE_SIZE_THRESHOLDS = [5]
ADDITIONAL_DISCORDANCE_THRESHOLD = 0.5
SUBCLONE_THRESHOLD = 0.9
USE_ALL_PAIRINGS = True
USE_ROUNDED_CN = True
USE_SUBCLONE_FILTER = True
ICHOR_RESULTS_FILEPATH = "projects/PDX/ichor/run_7_9_GHtest/results/"
PLOT_OUTPUT_PATH = "projects/PDX/production/ichor_minimal_subclone_discordance/plots/1mb_discordance/"
MEDIANS_FILENAME = "1mb_medians.txt"
PAIRING_DISCORDANCE_VALUES_FILENAME = "pairing_discordance_values.txt"
WHITE_SAMPLES_FILENAME = "white_samples.txt"


#global variables
"""
--FILEPATH_PAIRINGS_FILES - Path to folder containing pairings files for each
cohort.
--COHORT_NAMES_SELECT - List of cohort names to analyze.
"""
FILEPATH_PAIRINGS_FILES = "projects/PDX/pairings/"
COHORT_NAMES_SELECT = ["bcm_snp_breast", "europdx_wgs_brca", "europdx_wgs_crc", 
                "hci_snp_breast", "hci_wes_breast", "jax_snp_bladder", 
                "jax_snp_breast", "jax_snp_colorectal", "jax_snp_gbm", 
                "jax_snp_luad", "jax_snp_lusc", "jax_snp_melanoma", 
                "jax_snp_otherlung", "jax_snp_othertumor", "jax_snp_ovarian", 
                "jax_snp_sarcoma", "jax_snu_wes_gastric", "jax_wes", 
                "mdacc_wes_luad", "mdacc_wes_lusc", "mdacc_wes_otherlung", 
                "pdmr_wes_bladder", "pdmr_wes_colorectal", "pdmr_wes_headneck", 
                "pdmr_wes_lung", "pdmr_wes_othertumor", "pdmr_wes_pancreas", 
                "pdmr_wes_renal", "pdmr_wes_sarcoma", "pdmr_wes_skin", 
                "sibs_snp_hcc", "wistar_wes_melanoma", "wustl_wes_breast",
                "wustl_wes_pancreatic"]


#functions
def read_ichor_results(ichor_results_filepath, cohort_name, use_rounded_cn,
                       use_subclone_filter, subclone_threshold):
    """
    Given a filepath to ichor results organized into folders by cohort and a
    cohort name, returns a dictionary for the ichor results of the given
    cohort, where keys are sample names and values are ploidy-corrected copy
    number calls for each bin in the genome.
    """
    sample_info_dict = dict()
    cohort_results_folder = ichor_results_filepath + cohort_name
    sample_cn_files = glob.glob(cohort_results_folder + "/*/*.cna.seg")
    #for each sample
    for sample_cn_file in sample_cn_files:
        sample_param_file = sample_cn_file.replace(".cna.seg", ".params.txt")
        ploidy, subclone_fraction = extract_ploidy_and_subclone_fraction(sample_param_file)
        
        #get values that CN 1, 2, and 3 transform to with ploidy normalization,
        #used below for subclones
        if subclone_fraction != "NA":
            ploidy_adjusted_cn_1 = integer_round(1 / ploidy * 2)
            ploidy_adjusted_cn_2 = integer_round(2 / ploidy * 2)
            ploidy_adjusted_cn_3 = integer_round(3 / ploidy * 2)
        
        sample = sample_cn_file.split("/")[6]
        cn_list = []
        with open(sample_cn_file, "r") as f:
            f.readline()
            for line in f.readlines():
                line_parts = line.strip().split()
                #if there was no input data for that bin
                if line_parts[5] == "NA":
                    cn = "NA"
                #if there was input data for that bin
                else:
                    if use_rounded_cn:
                        #grab corrected copy number
                        cn = float(line_parts[7])
                    else:
                        #grab logR copy number
                        cn = float(line_parts[9])
                    #adjust CN if sublone filter is on, bin is subclonal,
                    #and subclone fraction is < given threshold
                    if use_subclone_filter:
                        is_subclone = float(line_parts[6])
                        #if the bin is subclonal and the subclone has less
                        #cellular prevalence than the set threshold
                        if is_subclone and (subclone_fraction < subclone_threshold):
                            #if the bin has less cellular prevalence than
                            #1 - the threshold, set the CN to 2
                            if subclone_fraction < (1 - subclone_threshold):
                                cn = 2.0
                            #otherwise, determine which ploidy-adjusted CNs the
                            #bin has character of
                            else:
                                #if it was a subclonal loss pre ploidy-correction
                                if integer_round(cn) == 1:
                                    #it now has character of everything from the ploidy-adjusted
                                    #value of what used to be CN 1 to the ploidy-adjusted
                                    #value of what used to be CN 2
                                    cn = str(ploidy_adjusted_cn_1) + "-" + str(ploidy_adjusted_cn_2)
                                #if it was a subclonal gain pre ploidy-correction
                                else:
                                    #it now has character of everything from the ploidy-adjusted
                                    #value of what used to be CN 2 to the ploidy-adjusted
                                    #value of what used to be CN 3
                                    cn = str(ploidy_adjusted_cn_2) + "-" + str(ploidy_adjusted_cn_3)
                #add CN to the CN list after adjusting for ploidy
                if isinstance(cn, str):
                    cn_list.append(cn)
                else:
                    ploidy_adjusted_cn = cn / ploidy * 2
                    cn_list.append(ploidy_adjusted_cn)
        sample_info_dict[sample] = cn_list
    return sample_info_dict


def integer_round(number):
    """
    Given a float, returns the integer it rounds to.  Python default rounding
    rounds 2.5 to 2, for example, favoring even numbers, while this function
    rounds 2.5 to 3.
    """
    if number - math.floor(number) < 0.5:
        return math.floor(number)
    else:
        return math.ceil(number)


def extract_ploidy_and_subclone_fraction(param_file):
    """
    Given a .params.txt file, returns the ploidy and subclone fraction reported
    therein.
    """
    with open(param_file, "r") as f:
        f.readline()
        info = f.readline().strip().split()
        ploidy = float(info[2])
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        subclone_info = f.readline().strip().split()
        subclone_fraction = subclone_info[-1]
    if subclone_fraction != "NA":
        subclone_fraction = float(subclone_fraction)
    return ploidy, subclone_fraction


def get_pairings_with_pt_list(pairings_files_filepath, cohort_name, use_all_pairings):
    """
    Given a pairings files filepath and a cohort name, reads the pairings file
    for that cohort and returns a list of pairings where each element in the
    list is a 3-element list itself specifying the first and second sample in
    the pairing and the PT those samples came from ("NA" if no PT sample).
    If use_all_pairings is given as True, will use the pairings files for all
    pairings as described in "settings" above, instead of only consecutive
    pairings.
    """
    if use_all_pairings and (cohort_name in ["jax_snp_otherlung", "pdmr_wes_bladder",
                                            "pdmr_wes_colorectal", "pdmr_wes_headneck",
                                            "pdmr_wes_lung", "pdmr_wes_pancreas",
                                            "pdmr_wes_renal", "pdmr_wes_sarcoma",
                                            "pdmr_wes_skin", "pdmr_wes_othertumor",
                                            "hci_snp_breast", "wistar_wes_melanoma",
                                            "sibs_snp_hcc", "europdx_wgs_crc",
                                            "europdx_wgs_brca"]):
        pairings_file = pairings_files_filepath + cohort_name + "_all_pairings_with_pt.txt"
    else:
        #the other cohorts only have 1 possible pairing scheme
        pairings_file = pairings_files_filepath + cohort_name + "_with_pt.txt"
    pairings_list = []
    with open(pairings_file, "r") as f:
        #remove header
        f.readline()
        lines = f.readlines()
        for line in lines:
            line_parts = line.strip().split()
            #if the line contains a pairing and not just a sample name from
            #the all_samples column
            if len(line_parts) > 1:
                pair1 = line_parts[0]
                pair2 = line_parts[1]
                pt = line_parts[2]
                pairings_list.append([pair1, pair2, pt])
    return pairings_list


def discordance_by_passage_comparison(cn_dict, pairings_list, additional_discordance_threshold,
                                      cohort_name, ichor_results_filepath,
                                      plot_output_path, pairing_discordance_values_filename,
                                      white_samples_filename):
    """
    Given a copy number dictionary and the corresponding pairings list and cohort
    name, returns two dictionaries where passage comparisons ("PT:P0", for
    example) are keys and values are lists of the number of bins discordant or
    fraction genome discordant for all passage comparisons in the pairings list.
    Additionally, writes the discordance values per pair to the given file, and
    the white samples to the given file.
    """
    all_white_samples = set()
    non_white_checked_samples = set()
    num_bins_discordant_by_passage_comparison_dict = dict()
    percent_genome_discordant_by_passage_comparison_dict = dict()
    
    with open(plot_output_path + pairing_discordance_values_filename, "a") as f:
        for pair in pairings_list:
            pair1 = pair[0]
            pair2 = pair[1]
            all_white_samples, non_white_checked_samples = update_all_white_samples_set(cn_dict,
                                                                                        all_white_samples,
                                                                                        non_white_checked_samples,
                                                                                        pair1,
                                                                                        cohort_name,
                                                                                        ichor_results_filepath)
            all_white_samples, non_white_checked_samples = update_all_white_samples_set(cn_dict,
                                                                                        all_white_samples,
                                                                                        non_white_checked_samples,
                                                                                        pair2,
                                                                                        cohort_name,
                                                                                        ichor_results_filepath)
            #if both samples are presumed to have adequate tumor fraction to analyze
            if (pair1 not in all_white_samples) and (pair2 not in all_white_samples):
                is_discordant_list = get_discordance_list(cn_dict, pair1, pair2,
                                                          additional_discordance_threshold)
                num_bins_discordant = is_discordant_list.count("1")
                #percent genome discordant is out of bins where both samples have values
                percent_genome_discordant = num_bins_discordant / (num_bins_discordant + is_discordant_list.count("0"))
                pair1_passage = get_passage(pair1)
                pair2_passage = get_passage(pair2)
                passage_comparison = pair1_passage + ":" + pair2_passage
                #add these values to the cohort's dictionary and write them
                #to the output file
                num_bins_discordant_by_passage_comparison_dict = update_dict(num_bins_discordant_by_passage_comparison_dict,
                                                                             passage_comparison,
                                                                             num_bins_discordant)
                percent_genome_discordant_by_passage_comparison_dict = update_dict(percent_genome_discordant_by_passage_comparison_dict,
                                                                                   passage_comparison,
                                                                                   percent_genome_discordant)
        
                f.write(cohort_name + "\t" + pair1 + "\t" + pair2 + "\t" + str(percent_genome_discordant) + "\n")
    with open(plot_output_path + white_samples_filename, "a") as f:
        for sample in all_white_samples:
            f.write(cohort_name + "\t" + sample + "\n")
    return num_bins_discordant_by_passage_comparison_dict, percent_genome_discordant_by_passage_comparison_dict

            
def update_dict(dictionary, key, value):
    """
    Given a dictionary, key, and value, appends the value to the list at
    dictionary[key] if key is already in the dictionary, and starts a list
    containing the value at dictionary[key] if key is new.
    """
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]
    return dictionary 


def get_passage(sample_name):
    """
    Given a string sample_name, returns a string with the passage information
    contained in the sample name (PT, P0, P1, etc...).
    """
    #look for passage information pattern in sample_name
    regex_results = re.match("([A-Z0-9a-z_-]+).(P[T0-9]+)", sample_name)
    #the passage information is the second element of the results
    passage = regex_results.groups()[1]
    return passage


def get_discordance_list(cn_dict, pair1, pair2, additional_discordance_threshold):
    """
    Given a copy number dictionary and 2 sample names, returns a list with
    a value for each bin in the two input dictionaries.  Values of "0" mean
    the two samples were not discordant at that location, whereas "1" means
    they were.  "NA" means there was no data for that bin in at least one of
    the two samples.
    """
    is_discordant_list = []
    #for each bin
    for i in range(len(cn_dict[pair1])):
        cn_pair1 = cn_dict[pair1][i]
        cn_pair2 = cn_dict[pair2][i]
        #if there is no data at that bin
        if cn_pair1 == "NA" or cn_pair2 == "NA":
            is_discordant_list.append("NA")
        #if at least one sample is subclonal
        elif isinstance(cn_pair1, str) or isinstance(cn_pair2, str):
            #if both samples are subclonal, call discordance
            if isinstance(cn_pair1, str) and isinstance(cn_pair2, str):
                discordance = call_subclone_subclone_discordance(cn_pair1, cn_pair2)
                is_discordant_list.append(discordance)
            #if sample 1 is subclonal and sample 2 is not, call discordance
            elif isinstance(cn_pair1, str):
                discordance = call_subclone_clone_discordance(cn_pair1, cn_pair2)
                is_discordant_list.append(discordance)
            #if sample 2 is subclonal and sample 1 is not, call discordance
            elif isinstance(cn_pair2, str):
                discordance = call_subclone_clone_discordance(cn_pair2, cn_pair1)
                is_discordant_list.append(discordance)
        #if neither sample is subclonal
        else:
            pair1_cn_call = call_cn_status(cn_pair1)
            pair2_cn_call = call_cn_status(cn_pair2)
            #if the calls are the same, they are not discordant
            if pair1_cn_call == pair2_cn_call:
                is_discordant_list.append("0")
            #if the calls are not the same
            else:
                #if they differ by less than the given threshold, they are not discordant
                if abs(cn_pair1 - cn_pair2) < additional_discordance_threshold:
                    is_discordant_list.append("0")
                #if they differ by more than or equal to the given threshold, they are discordant
                else:
                    is_discordant_list.append("1")
    return is_discordant_list


def call_subclone_subclone_discordance(subclonal_cn_call_1, subclonal_cn_call_2):
    """
    Given two subclonal calls, return "0" if they are not discordant (their
    ranges overlap) and "1" if they are discordant (their ranges do not overlap).
    """
    subclone_1_min = int(subclonal_cn_call_1[0])
    subclone_1_max = int(subclonal_cn_call_1[2])
    subclone_2_min = int(subclonal_cn_call_2[0])
    subclone_2_max = int(subclonal_cn_call_2[2])
    #if the ranges overlap,the they are not considered discordant
    if (subclone_1_min <= subclone_2_max) and (subclone_2_min <= subclone_1_max):
        return "0"
    else:
        return "1"


def call_subclone_clone_discordance(subclonal_cn_call, clonal_cn_call):
    """
    Given a subclonal CN call and a clonal CN call, returns "0" if they are
    not discordant (the range of the subclone call includes the clonal call)
    and "1" if they are discordant (the range of the subclone call does not
    include the clonal call).
    """
    subclone_min = int(subclonal_cn_call[0])
    subclone_max = int(subclonal_cn_call[2])
    rounded_clonal_cn_call = integer_round(clonal_cn_call)
    if (subclone_min <= rounded_clonal_cn_call) and (rounded_clonal_cn_call <= subclone_max):
        return "0"
    else:
        return "1"


def call_cn_status(copy_number):
    """
    Given a float copy_number, returns "neut" if copy_number rounds to 2,
    "loss" if it rounds to < 2, and "gain" if it rounds to > 2.
    """
    rounded_cn = integer_round(copy_number)
    if rounded_cn == 2:
        cn_status = "neut"
    elif rounded_cn < 2:
        cn_status = "loss"
    else:
        cn_status = "gain"
    return cn_status
        
    
def update_all_white_samples_set(cn_dict, all_white_samples,
                                 non_white_checked_samples, sample,
                                 cohort_name, ichor_results_filepath):
    """
    Given a dictionary where keys are patients and values are lists of copy numbers
    at different bins, a set of all_white_samples, a set of non-white checked
    samples, and a sample name string, returns an updated set with the
    new sample added to the appropriate set depending on if it is all-white or
    not, where all-white is defined as < 5% tumor fraction and/or <5% of bins
    with data altered from neutral.
    """
    #if sample has already been checked, don't recheck
    if (sample in all_white_samples) or (sample in non_white_checked_samples):
        return all_white_samples, non_white_checked_samples
    #if new sample, check to see if all-white
    else:
        sample_params_file =  ichor_results_filepath + cohort_name + "/" + sample + "/" + sample + ".params.txt"
        with open(sample_params_file, "r") as f:
            f.readline()
            info = f.readline()
            tumor_fraction = float(info.strip().split()[1])
        #if tfx < 5%, classify as all-white sample
        if tumor_fraction < 0.05:
            all_white_samples.add(sample)
        else:
            #also, if <5% bins are altered, classify as all-white sample
            cn_sample = cn_dict[sample]
            #remove NAs
            cn_sample_no_na = [cn for cn in cn_sample if cn != "NA"]
            #determine whether each CN in the list is neutral or, if subclonal,
            #has partial neutral character after ploidy normalization
            cn_sample_no_na_2_character = [has_cn_2_character(cn) for cn in cn_sample_no_na]
            #calculate percent of bins altered
            num_bins_neutral = cn_sample_no_na_2_character.count(True)
            num_bins_data = len(cn_sample_no_na)
            percent_altered = (num_bins_data - num_bins_neutral) / num_bins_data
            if percent_altered < 0.05:
                all_white_samples.add(sample)
            #if neither tfx nor altered bin thresholds capture sample, it is not all-white
            else:
                non_white_checked_samples.add(sample)
        return all_white_samples, non_white_checked_samples
    

def has_cn_2_character(cn_call):
    """
    Given a CN, returns True if the CN has at least partial neutral character,
    and False if not.  If the input is a float, this indicates clonal CN, and
    therefore True is returned iff the input rounds to 2.  If the input is a
    string of a range, e.g. "1-2", indicating subclonal CN, returns True if the
    range includes 2 and False if not.
    """
    #if the call is clonal, return True if CN = 2 and False if not
    if isinstance(cn_call, float):
        if integer_round(cn_call) == 2:
            return True
        else:
            return False
    #if the call is subclonal, return True if it has partial CN 2 character
    else:
        cn_call_min = int(cn_call[0])
        cn_call_max = int(cn_call[2])
        if (cn_call_min <= 2) and (2 <= cn_call_max):
            return True
        else:
            return False


def get_discordance_by_cohort(cohort_names, ichor_results_filepath, pairings_filepath,
                              additional_discordance_threshold,
                              use_all_pairings, use_rounded_cn, use_subclone_filter,
                              subclone_threshold, plot_output_path,
                              pairing_discordance_values_filename,
                              white_samples_filename):
    """
    Given a list of cohort names, a filepath to ichorCNA's results, and a
    path to the pairings files, returns num_bins_by_cohort_dict and
    percent_genome_by_cohort_dict.  The keys of these are cohort names.  The
    values are dictionaries where keys are passage comparisons and values are
    lists of discordance values.
    """
    num_bins_by_cohort_dict = dict()
    percent_genome_by_cohort_dict = dict()
    for cohort_name in cohort_names:
        print(cohort_name)
        #load in CN information and pairings list
        cn_dict = read_ichor_results(ichor_results_filepath,
                                          cohort_name, use_rounded_cn,
                                          use_subclone_filter,
                                          subclone_threshold)
        pairings_list = get_pairings_with_pt_list(pairings_filepath,
                                                  cohort_name, use_all_pairings)
        #get dictionaries of num_bins and percent_genome discordant, where
        #keys specify which passages are being compared and values are lists
        #of discordance values
        num_bins, percent_genome = discordance_by_passage_comparison(cn_dict,
                                                                     pairings_list,
                                                                     additional_discordance_threshold,
                                                                     cohort_name,
                                                                     ichor_results_filepath,
                                                                     plot_output_path,
                                                                     pairing_discordance_values_filename,
                                                                     white_samples_filename)
        num_bins_by_cohort_dict[cohort_name] = num_bins
        percent_genome_by_cohort_dict[cohort_name] = percent_genome
    return num_bins_by_cohort_dict, percent_genome_by_cohort_dict


def plot_discordance_by_cohort(discordance_dictionary, cohort_names, ylabel,
                               title, ylimit, plot_filename, sample_size_threshold,
                               medians_file_object):
    """
    Given a dictionary as described in get_discordance_by_cohort, generates
    a boxplot where each cohort is a box and boxes are sorted by median
    discordance.
    """
    cohort_medians = []
    cohort_means = []
    cohorts_with_data = []
    minimum = float("inf")
    maximum = float("-inf")
    
    for cohort_name in cohort_names:
        cohort_values = discordance_dictionary[cohort_name]
        if len(cohort_values) >= sample_size_threshold:
            median = statistics.median(cohort_values)
            mean = statistics.mean(cohort_values)
            cohorts_with_data.append(cohort_name)
            cohort_medians.append(median)
            cohort_means.append(mean)
            if min(cohort_values) < minimum:
                minimum = min(cohort_values)
            if max(cohort_values) > maximum:
                maximum = max(cohort_values)
    cohorts_sorted_by_median = [cohort for _, cohort in sorted(zip(cohort_medians, cohorts_with_data))]
    
    mean_of_means = statistics.mean(cohort_means)
    medians_file_object.write("#mean of means = " + str(round(mean_of_means, 5)) + "\n")
    medians_file_object.write("#minimum discordance = " + str(round(minimum, 5)) + "\n")
    medians_file_object.write("#maximum discordance = " + str(round(maximum, 5)) + "\n")
    
    sorted_cohort_medians = sorted(cohort_medians)
    for i in range(len(cohorts_sorted_by_median)):
        print(cohorts_sorted_by_median[i], sorted_cohort_medians[i])
        medians_file_object.write(cohorts_sorted_by_median[i] + " " + str(sorted_cohort_medians[i]) + "\n")
    
    median_of_medians = statistics.median(cohort_medians)
    plt.figure(figsize=(12, 6))
    plt.axhline(y = median_of_medians, color = "blue")
    plt.text(1, ylimit, "median of medians = " + str(round(median_of_medians, 5)),
             fontsize = 14)
    values_in_sorted_key_order = []
    for cohort_name in cohorts_sorted_by_median:
        values_in_sorted_key_order.append(discordance_dictionary[cohort_name])
        
    box_data = plt.boxplot(values_in_sorted_key_order, patch_artist = True)
    
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box_data[item], color="blue")
    plt.setp(box_data["boxes"], facecolor="springgreen")
    plt.setp(box_data["fliers"], markeredgecolor="springgreen")
    
    #setting axis labels
    plt.xticks(range(1, len(cohorts_sorted_by_median) + 1),
               cohorts_sorted_by_median, rotation = 90, fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.ylim(0 - ylimit * 0.1, ylimit * 1.1)
    plt.xlabel("Cohort", fontsize = 16)
    plt.ylabel(ylabel, fontsize = 16)
    plt.title(title, fontsize = 18)
    plt.savefig(plot_filename, bbox_inches = "tight")
    plt.show()


def break_up_pt_and_passage_comparisons(discordance_by_cohort_dict):
    """
    Given a discordance_by_cohort_dict, returns two dictionaries, one containing
    the data for passage comparisons that contain PT, and one for the comparisons
    between earlier and later passages.
    """
    pt_passage_dict = dict()
    early_late_dict = dict()
    for cohort in discordance_by_cohort_dict:
        discordance_dictionary = discordance_by_cohort_dict[cohort]
        pt_passage_discordance = []
        early_late_discordance = []
        for key in discordance_dictionary:
            if "PT" in key:
                pt_passage_discordance += discordance_dictionary[key]
            else:
                early_late_discordance += discordance_dictionary[key]
        pt_passage_dict[cohort] = pt_passage_discordance
        early_late_dict[cohort] =  early_late_discordance
    return pt_passage_dict, early_late_dict


def sorted_discordance_plots(ichor_results_filepath, pairings_filepath, cohort_names,
                             bin_size_name, bin_size_abbr, num_bins_ylimit,
                             plot_output_path, sample_size_thresholds, 
                             additional_discordance_threshold,
                             medians_filename, use_all_pairings,
                             use_rounded_cn, use_subclone_filter,
                             subclone_threshold, pairing_discordance_values_filename,
                             white_samples_filename):
    """
    Plots sorted discordance plots, both number of bins discordant
    and fraction genome discordant.
    """
    num_bins_by_cohort_dict, percent_genome_by_cohort_dict = get_discordance_by_cohort(cohort_names,
                                                                                       ichor_results_filepath,
                                                                                       pairings_filepath,
                                                                                       additional_discordance_threshold,
                                                                                       use_all_pairings,
                                                                                       use_rounded_cn,
                                                                                       use_subclone_filter,
                                                                                       subclone_threshold,
                                                                                       plot_output_path,
                                                                                       pairing_discordance_values_filename,
                                                                                       white_samples_filename)
    
    num_bins_pt_passage, num_bins_early_late = break_up_pt_and_passage_comparisons(num_bins_by_cohort_dict)
    percent_genome_pt_passage, percent_genome_early_late = break_up_pt_and_passage_comparisons(percent_genome_by_cohort_dict)
    
    f = open(plot_output_path + medians_filename, "a+")
    f.write("###Additional discordance threshold = {0} diff, subclone filtering = {1}, subclone threshold = {2}\n".format(additional_discordance_threshold, str(use_subclone_filter), str(subclone_threshold)))
    for sample_size_threshold in sample_size_thresholds:
        f.write("\n\n###Sample size threshold = {0}\n".format(sample_size_threshold))
        f.write("\n###Number of " + bin_size_name + " Discordant\n")
        f.write("###Discordance of PT vs. PDX Passage\n")
        plot_discordance_by_cohort(num_bins_pt_passage, cohort_names,
                                   "Number of " + bin_size_name + " Discordant",
                                   "Discordance of PT vs. PDX Passage",
                                   num_bins_ylimit,
                                   "{0}pt_vs_pdx_num_bins_{1}_{2}_diff_threshold_n{3}.png".format(plot_output_path, bin_size_abbr, str(additional_discordance_threshold), str(sample_size_threshold)),
                                   sample_size_threshold, f)
        f.write("\n###Number of " + bin_size_name + " Discordant\n")
        f.write("###Discordance of PDX Passage vs. Later PDX Passage\n")
        plot_discordance_by_cohort(num_bins_early_late, cohort_names,
                                   "Number of " + bin_size_name + " Discordant",
                                   "Discordance of PDX Passage vs. Later PDX Passage",
                                   num_bins_ylimit,
                                   "{0}early_vs_late_num_bins_{1}_{2}_diff_threshold_n{3}.png".format(plot_output_path, bin_size_abbr, str(additional_discordance_threshold), str(sample_size_threshold)),
                                   sample_size_threshold, f)
        f.write("\n###Fraction of " + bin_size_name + " Discordant\n")
        f.write("###Discordance of PT vs. PDX Passage\n")
        plot_discordance_by_cohort(percent_genome_pt_passage, cohort_names,
                                   "Fraction of " + bin_size_name + " Discordant",
                                   "Discordance of PT vs. PDX Passage", 1,
                                   "{0}pt_vs_pdx_percent_genome_{1}_{2}_diff_threshold_n{3}.png".format(plot_output_path, bin_size_abbr, str(additional_discordance_threshold), str(sample_size_threshold)),
                                   sample_size_threshold, f)
        f.write("\n###Fraction of " + bin_size_name + " Discordant\n")
        f.write("###Discordance of PDX Passage vs. Later PDX Passage\n")
        plot_discordance_by_cohort(percent_genome_early_late, cohort_names,
                                   "Fraction of " + bin_size_name + " Discordant",
                                   "Discordance of PDX Passage vs. Later PDX Passage", 1,
                                   "{0}early_vs_late_percent_genome_{1}_{2}_diff_threshold_n{3}.png".format(plot_output_path, bin_size_abbr, str(additional_discordance_threshold), str(sample_size_threshold)),
                                   sample_size_threshold, f)
    f.write("\n\n\n")
    f.close()


   
#run
if __name__ == "__main__":
    sorted_discordance_plots(ICHOR_RESULTS_FILEPATH, FILEPATH_PAIRINGS_FILES,
                             COHORT_NAMES_SELECT, "1 Mb Bins",
                            "1mb", 2500, PLOT_OUTPUT_PATH, SAMPLE_SIZE_THRESHOLDS,
                            ADDITIONAL_DISCORDANCE_THRESHOLD,
                            MEDIANS_FILENAME, USE_ALL_PAIRINGS,
                            USE_ROUNDED_CN, USE_SUBCLONE_FILTER,
                            SUBCLONE_THRESHOLD, PAIRING_DISCORDANCE_VALUES_FILENAME,
                            WHITE_SAMPLES_FILENAME)