#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Delta Groupings Ichor Minimal Subclone Discordance
July 22nd, 2020
Anna Hoge
Ha Lab
Fred Hutchinson Cancer Research Center
"""

#imports
import os
import matplotlib.pyplot as plt
import re
import random
import scipy.stats
import glob
import statistics
import math


#working directory
#os.chdir("/Volumes/fh")
os.chdir("/Volumes/ha_g")
#os.chdir("/Volumes/ha_g-1")
#os.chdir("/Volumes/ha_g-2")


#settings
"""
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
--PLOT_OUTPUT_PATH - Path to folder where output plots and files should go.
--PAIRING_DISCORDANCE_VALUES_FILENAME - Name of output file listing discordance
values for each analyzed pairing.
--WHITE_SAMPLES_FILENAME - Name of output file listing samples which had <5%
of bins with ploidy-adjusted CN != 2 and/or < 5% tumor fraction and were thus
not included in the discordance analysis.
"""
ADDITIONAL_DISCORDANCE_THRESHOLD = 0.5
PLOT_OUTPUT_PATH = "projects/PDX/production/ichor_minimal_subclone_discordance/plots/delta_groupings/"
ICHOR_RESULTS_FILEPATH = "projects/PDX/ichor/run_7_9_GHtest/results/"
SUBCLONE_THRESHOLD = 0.9
USE_ROUNDED_CN = True
USE_SUBCLONE_FILTER = True
P_VAL_FILENAME = "delta_groupings_p_vals.txt"
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


def get_pairings_with_pt_list(pairings_files_filepath, cohort_name):
    """
    Given a pairings files filepath and a cohort name, reads the pairings file
    for that cohort and returns a list of pairings where each element in the
    list is a 3-element list itself specifying the first and second sample in
    the pairing and the PT those samples came from ("NA" if no PT sample).
    """
    #these cohorts have multiple possible pairing schemes.  select correct one
    #for this analysis
    if "pdmr" in cohort_name:
        pairings_file = pairings_files_filepath + cohort_name + "_for_passage_number_change_with_pt2.txt"
    elif "sibs" in cohort_name or "wistar" in cohort_name:
        pairings_file = pairings_files_filepath + cohort_name + "_for_passage_number_change_with_pt.txt"
    else:
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


def discordance_by_delta(cn_dict, pairings_list, additional_discordance_threshold,
                         ichor_results_filepath, cohort_name, plot_output_path,
                         pairing_discordance_values_filename,
                         white_samples_filename):
    """
    Given a copy number dict and and pairings list, returns two dictionaries
    where deltas between passage comparisons are keys and values are lists of the
    number of bins discordant or fraction genome discordant for all passage
    comparisons of that delta.
    """
    with open(plot_output_path + pairing_discordance_values_filename, "a") as f:
        all_white_samples = set()
        non_white_checked_samples = set()
        num_bins_discordant_by_delta_dict = dict()
        percent_genome_discordant_by_delta_dict = dict()
        
        for pair in pairings_list:
            pair1 = pair[0]
            #not including PT pairings for this plot
            pair1_passage = get_passage(pair1)
            if pair1_passage != "PT":
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
                    pair2_passage = get_passage(pair2)
                    delta = int(pair2_passage[1:]) - int(pair1_passage[1:])
                    is_discordant_list = get_discordance_list(cn_dict, pair1, pair2,
                                                              additional_discordance_threshold)
                    num_bins_discordant = is_discordant_list.count("1")
                    #percent genome discordant is out of bins where both samples have CN values
                    percent_genome_discordant = num_bins_discordant / (num_bins_discordant + is_discordant_list.count("0"))
                    #add these values to the cohort's dictionary and write them
                    #to the output file
                    num_bins_discordant_by_delta_dict = update_dict(num_bins_discordant_by_delta_dict,
                                                                    delta,
                                                                    num_bins_discordant)
                    percent_genome_discordant_by_delta_dict = update_dict(percent_genome_discordant_by_delta_dict,
                                                                          delta,
                                                                          percent_genome_discordant)
                    
                    f.write(cohort_name + "\t" + pair1 + "\t" + pair2 + "\t" + str(delta) + "\t" + str(percent_genome_discordant) + "\n")
    with open(plot_output_path + white_samples_filename, "a") as f:
        for sample in all_white_samples:
            f.write(cohort_name + "\t" + sample + "\n")
    return num_bins_discordant_by_delta_dict, percent_genome_discordant_by_delta_dict

            
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
    

def union_dictionaries(dict1, dict2):
    """
    Given two dictionaries of lists dict1 and dict2, returns an updated version
    of dict1 that is the union of both inputs.
    """
    for key in dict2:
        if key in dict1:
            dict1[key] = dict1[key] + dict2[key]
        else:
            dict1[key] = dict2[key]
    return dict1


def discordance_by_passage_number_change_plot(percent_genome,
                                              plot_dots, plot_output_path,
                                              additional_discordance_threshold,
                                              grouping_name):
    """
    Given  a dictionary percent_genome, plots a violin plot of PDX passage
    differences vs the fraction of the genome discordant.
    """
    all_sorted_keys = sorted(percent_genome, key = lambda x: (x[1].isdigit(), x))
    #make list of values in sorted key order
    values_in_sorted_key_order = []
    sorted_keys = []
    for key in all_sorted_keys:
        values_in_sorted_key_order.append(percent_genome[key])
        sorted_keys.append(key)
    
    #violin plot
    violin_parts = plt.violinplot(values_in_sorted_key_order, showmedians = True)
    
    for i in range(len(values_in_sorted_key_order)):
        plt.text(i + 0.9, 1.0, "n=" + str(len(values_in_sorted_key_order[i])),
                 fontsize = 14)
    
    if plot_dots == True:
        #plot dots as scatter plot
        x_vals = []
        for i in range(len(values_in_sorted_key_order)):
            part_x_vals = [i + 1 for value in values_in_sorted_key_order[i]]
            x_vals += part_x_vals
        x_vals = [x + random.uniform(-0.05, 0.05) for x in x_vals]
        y_vals = [val for sublist in values_in_sorted_key_order for val in sublist]
        plt.scatter(x_vals, y_vals, alpha = 0.3, color = "blue")
    
    #make violin bodies green
    for violin in violin_parts["bodies"]:
        violin.set_facecolor("springgreen")
        violin.set_edgecolor("blue")
    #make violin lines blue
    for partname in ('cbars','cmins','cmaxes'):
        violin_part = violin_parts[partname]
        violin_part.set_edgecolor("blue")
        violin_part.set_linewidth(1)
    violin_parts["cmedians"].set_edgecolor("black")
    violin_parts["cmedians"].set_linewidth(2)
    
    #setting axis labels
    plt.xticks(range(1, len(sorted_keys) + 1), sorted_keys, fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.ylim(-0.1, 1.1)
    plt.xlabel("PDX Passage Difference", fontsize = 14)
    plt.ylabel("Fraction 1 Mb Bins Discordant", fontsize = 14)
    
    plot_filename = "{0}discordance_by_delta_1mb_{1}_diff_threshold_{2}.png".format(plot_output_path, str(additional_discordance_threshold), grouping_name)
    plt.savefig(plot_filename, bbox_inches = "tight")
            
    plt.show()


def group_dictionary_keys(dictionary, groupings):
    """
    Given a dictionary and a list of groupings where each element is a list of
    dictionary keys to be grouped together, returns grouped_dictionary.
    """
    grouped_dictionary = dict()
    for group in groupings:
        grouped_dictionary_key = ", ".join(str(number) for number in group)
        grouped_dictionary[grouped_dictionary_key] = []
        for number in group:
            grouped_dictionary[grouped_dictionary_key] += dictionary[number]
    return grouped_dictionary


def rank_sum_test(dictionary, key1, key2, p_val_file_object):
    """
    Given a dictionary of lists and two keys, conducts a rank
    sum test and prints the results as well as writes them to the given file.
    """
    test_stat, p_value = scipy.stats.mannwhitneyu(dictionary[key1],
                                                  dictionary[key2])
    median_key1 = round(statistics.median(dictionary[key1]), 4)
    median_key2 = round(statistics.median(dictionary[key2]), 4)
    if p_value < 0.05:
        print(key1 + " and " + key2 + " are different with p value = " + str(round(p_value, 4)))
        p_val_file_object.write(key1 + " (median = " + str(median_key1) + ") and " + key2 + " (median = " + str(median_key2) + ") are different with p value = " + str(round(p_value, 4)) + "\n")
    else:
        print("no difference between " + key1 + " and " + key2 + ", p value = " + str(round(p_value, 4)))
        p_val_file_object.write("no difference between " + key1 + " (median = " + str(median_key1) + ") and " + key2 + " (median = " + str(median_key2) + "), p value = " + str(round(p_value, 4)) + "\n")



#run   
if __name__ == "__main__":
    cohort_names = COHORT_NAMES_SELECT
    pairings_filepath = FILEPATH_PAIRINGS_FILES
    plot_output_path = PLOT_OUTPUT_PATH
    p_val_filepath = plot_output_path + P_VAL_FILENAME
    additional_discordance_threshold = ADDITIONAL_DISCORDANCE_THRESHOLD
    subclone_threshold = SUBCLONE_THRESHOLD
    use_rounded_cn = USE_ROUNDED_CN
    use_subclone_filter = USE_SUBCLONE_FILTER
    ichor_results_filepath = ICHOR_RESULTS_FILEPATH
    pairing_discordance_values_filename = PAIRING_DISCORDANCE_VALUES_FILENAME
    white_samples_filename = WHITE_SAMPLES_FILENAME
    
    with open(plot_output_path + pairing_discordance_values_filename, "w") as f:
        f.write("cohort\tsample1\tsample2\tdelta\tdiscordance\n")
    
    f = open(p_val_filepath, "w")
    f.write("###Additional discordance threshold = {0} diff, subclone filtering = {1}, subclone threshold = {2}\n".format(additional_discordance_threshold, str(use_subclone_filter), str(subclone_threshold)))
                
    all_cohorts_percent_genome_dict = dict()
    for cohort_name in cohort_names:
        print(cohort_name)
        #load in CN information and pairings list
        cn_dict = read_ichor_results(ichor_results_filepath, cohort_name, use_rounded_cn,
                   use_subclone_filter, subclone_threshold)
        pairings_list = get_pairings_with_pt_list(pairings_filepath,
                                                  cohort_name)
        #get dictionaries of num_bins and percent_genome discordant, where
        #keys specify passage delta and values are lists of discordance values
        num_bins, percent_genome = discordance_by_delta(cn_dict,
                                                        pairings_list,
                                                        additional_discordance_threshold,
                                                        ichor_results_filepath,
                                                        cohort_name,
                                                        plot_output_path,
                                                        pairing_discordance_values_filename,
                                                        white_samples_filename)
        all_cohorts_percent_genome_dict = union_dictionaries(all_cohorts_percent_genome_dict,
                                                             percent_genome)   
        
    #group passage deltas using 3 different grouping schemes
    grouped_dict1 = group_dictionary_keys(all_cohorts_percent_genome_dict, [[1, 2], [3, 4, 5, 6, 7, 8], [17, 18, 21]])
    grouped_dict2 = group_dictionary_keys(all_cohorts_percent_genome_dict, [[1, 2, 3], [4, 5, 6, 7, 8], [17, 18, 21]])
    grouped_dict3 = group_dictionary_keys(all_cohorts_percent_genome_dict, [[1, 2], [3, 4, 5], [6, 7, 8, 17, 18, 21]])
    
    #rank sum grouping1
    f.write("\n###Grouping 1\n")
    rank_sum_test(grouped_dict1, "1, 2", "3, 4, 5, 6, 7, 8", f)
    rank_sum_test(grouped_dict1, "3, 4, 5, 6, 7, 8", "17, 18, 21", f)
    rank_sum_test(grouped_dict1, "1, 2", "17, 18, 21", f)
    
    #rank sum grouping2
    f.write("\n###Grouping 2\n")
    rank_sum_test(grouped_dict2, "1, 2, 3", "4, 5, 6, 7, 8", f)
    rank_sum_test(grouped_dict2, "4, 5, 6, 7, 8", "17, 18, 21", f)
    rank_sum_test(grouped_dict2, "1, 2, 3", "17, 18, 21", f)
    
    #rank sum grouping3
    f.write("\n###Grouping 3\n")
    rank_sum_test(grouped_dict3, "1, 2", "3, 4, 5", f)
    rank_sum_test(grouped_dict3, "3, 4, 5", "6, 7, 8, 17, 18, 21", f)
    rank_sum_test(grouped_dict3, "1, 2", "6, 7, 8, 17, 18, 21", f)

    f.write("\n\n\n")
    
    #plot results
    discordance_by_passage_number_change_plot(grouped_dict1, False,
                                              plot_output_path,
                                              additional_discordance_threshold,
                                              "grouping1")
    discordance_by_passage_number_change_plot(grouped_dict2, False,
                                              plot_output_path,
                                              additional_discordance_threshold,
                                              "grouping2")
    discordance_by_passage_number_change_plot(grouped_dict3, False,
                                              plot_output_path,
                                              additional_discordance_threshold,
                                              "grouping3")

    f.close()
