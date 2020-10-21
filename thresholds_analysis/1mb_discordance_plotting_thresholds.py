#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1 Mb Discordance Plotting Thresholds
July 8th, 2020
Anna Hoge
Ha Lab
Fred Hutchinson Cancer Research Center
"""

#imports
import os
import matplotlib.pyplot as plt
import re
import statistics


#working directory
#os.chdir("/Volumes/fh")
os.chdir("/Volumes/ha_g")
#os.chdir("/Volumes/ha_g-1")
#os.chdir("/Volumes/ha_g-2")


#settings
"""
--SAMPLE_SIZE_THRESHOLDS - A list of integers to use as the minimum number of
fully analyzed pairings a cohort must have to be included in the plots.
--DISCORDANCE THRESHOLD PAIRS - A list of lists with the discordance thresholds
to use.  For example, [0.1, 0.4] means if logR_sample2 >= 0.3 and
logR_sample_1 < 0.3, for these samples to be discordant, logR_sample1 must
additionally be < 0.1 or be different by at least 0.4 from logR_sample2.
(Note that using >1 threshold pair will cause discordance values for each
threshold pair to be written to PAIRING_DISCORDANCE_VALUES_FILENAME in sequence.)
--USE_ALL_PAIRINGS - Boolean indicating whether all possible pairings between
a PT and its PDXs should be used, provided the PDXs are not of the same passage.
False would indicate only consecutive pairings should be used (e.g. given PT,
P1, and P3, pair PT-P1 and P1-P3 only, excluding PT-P3).
--PLOT_OUTPUT_PATH - Path to folder where output plots and files should go.
--MEDIANS_FILENAME - Name of output file listing cohort medians.
--PAIRING_DISCORDANCE_VALUES_FILENAME - Name of output file listing discordance
values for each analyzed pairing.
--WHITE_SAMPLES_FILENAME - Name of output file listing samples which had <5%
of bins with abs(logR) > 0.3 and were thus not included in the discordance
analysis.
"""
SAMPLE_SIZE_THRESHOLDS = [5]
DISCORDANCE_THRESHOLD_PAIRS = [[0.1, 0.3]]
USE_ALL_PAIRINGS = True
PLOT_OUTPUT_PATH = "projects/PDX/production/thresholds/plots/1mb_discordance/"
MEDIANS_FILENAME = "1mb_medians.txt"
PAIRING_DISCORDANCE_VALUES_FILENAME = "pairing_discordance_values.txt"
WHITE_SAMPLES_FILENAME = "white_samples.txt"


#global variables
"""
--MB_DATA_TABLES_FILEPATH - Path to folder containing tsv files for each cohort,
with 1 Mb bin logR values for each sample.
--FILEPATH_PAIRINGS_FILES - Path to folder containing pairings files for each
cohort.
--COHORT_NAMES_SELECT - List of cohort names to analyze.
"""
MB_DATA_TABLES_FILEPATH = "projects/PDX/1mb_severe_logr_and_cn_tables/"
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
def read_sample_info_file(filepath, cohort_name):
    """
    Given a filepath and a cohort name, returns a dictionary where keys are
    sample names and values are lists of logR values for each bin of that sample.
    """
    sample_info_dict = dict()
    sample_info_file = filepath + cohort_name + "_logr.tsv"
    with open(sample_info_file, "r") as f:
        for line in f.readlines():
            line_parts = line.strip().split()
            sample = line_parts[0].replace("\"", "")
            sample_info = line_parts[1:]
            sample_info_dict[sample] = sample_info
    return sample_info_dict


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


def discordance_by_passage_comparison(logr_dict, pairings_list, discordance_threshold1,
                                      discordance_threshold2, cohort_name,
                                      plot_output_path, pairing_discordance_values_filename,
                                      white_samples_filename):
    """
    Given a logR dict and pairings list, returns two dictionaries
    where passage comparisons ("PT:P0", for example) are keys and values are
    lists of the number of bins discordant or fraction genome discordant for all
    passage comparisons in the pairings list.  Additionally, writes the
    discordance values per pair to the given file, and the white samples to
    the given file.
    """
    all_white_samples = set()
    non_white_checked_samples = set()
    num_bins_discordant_by_passage_comparison_dict = dict()
    percent_genome_discordant_by_passage_comparison_dict = dict()
    
    with open(plot_output_path + pairing_discordance_values_filename, "a") as f:
        for pair in pairings_list:
            pair1 = pair[0]
            pair2 = pair[1]
            all_white_samples, non_white_checked_samples = update_all_white_samples_set(logr_dict,
                                                                                        all_white_samples,
                                                                                        non_white_checked_samples,
                                                                                        pair1)
            all_white_samples, non_white_checked_samples = update_all_white_samples_set(logr_dict,
                                                                                        all_white_samples,
                                                                                        non_white_checked_samples,
                                                                                        pair2)
            #if both samples are presumed to have adequate tumor fraction to analyze
            if (pair1 not in all_white_samples) and (pair2 not in all_white_samples):
                is_discordant_list = get_discordance_list(logr_dict, pair1, pair2,
                                                          discordance_threshold1,
                                                          discordance_threshold2)
                num_bins_discordant = is_discordant_list.count("1")
                #percent genome discordant is out of bins where both samples have logR values
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
    #add all-white samples to the all-white output file
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


def get_discordance_list(logr_dict, pair1, pair2, discordance_threshold1, discordance_threshold2):
    """
    Given a logR dict and 2 sample names, returns a list with
    a value for each bin in the two input dictionaries.  Values of "0" mean
    the two samples were not discordant at that location, whereas "1" means
    they were.  "NA" means there was no data for that bin in at least one of
    the two samples.
    """
    is_discordant_list = []
    #for each bin
    for i in range(len(logr_dict[pair1])):
        logr_pair1 = logr_dict[pair1][i]
        logr_pair2 = logr_dict[pair2][i]
        #if there is no data at that bin
        if logr_pair1 == "NA" or logr_pair2 == "NA":
            is_discordant_list.append("NA")
        #if pair1 is a gain and pair2 is not
        elif (float(logr_pair1) >= 0.3) and ((float(logr_pair2) < discordance_threshold1) or ((float(logr_pair2) < 0.3) and (float(logr_pair1) - float(logr_pair2) >= discordance_threshold2))):
            is_discordant_list.append("1")
        #if pair2 is a gain and pair1 is not
        elif (float(logr_pair2) >= 0.3) and ((float(logr_pair1) < discordance_threshold1) or ((float(logr_pair1) < 0.3) and (float(logr_pair2) - float(logr_pair1) >= discordance_threshold2))):
            is_discordant_list.append("1")
        #if pair1 is a loss and pair2 is not
        elif (float(logr_pair1) <= -0.3) and ((float(logr_pair2) > -discordance_threshold1) or ((float(logr_pair2) > -0.3) and (float(logr_pair1) - float(logr_pair2) <= -discordance_threshold2))):
            is_discordant_list.append("1")
        #if pair2 is a loss and pair1 is not
        elif (float(logr_pair2) <= -0.3) and ((float(logr_pair1) > -discordance_threshold1) or ((float(logr_pair1) > -0.3) and (float(logr_pair2) - float(logr_pair1) <= -discordance_threshold2))):
            is_discordant_list.append("1")
        #otherwise, they are not discordant
        else:
            is_discordant_list.append("0")
    return is_discordant_list

        
def update_all_white_samples_set(logr_dict, all_white_samples, non_white_checked_samples,
                                 sample):
    """
    Given a dictionary where keys are patients and values are lists of logR
    status at different bins, a set of all_white_samples, a set of non-white
    checked samples, and a sample name string, returns an updated set with the
    new sample added to the appropriate set depending on if it is all-white or not.
    """
    #if sample already checked, don't recheck
    if (sample in all_white_samples) or (sample in non_white_checked_samples):
        return all_white_samples, non_white_checked_samples
    #if new sample, check to see if all-white
    else:
        logr_sample = logr_dict[sample]
        abs_logr_sample_no_na = [abs(float(logr)) for logr in logr_sample if logr != "NA"]
        num_bins_above_threshold = len([logr for logr in abs_logr_sample_no_na if logr >= 0.3])
        num_bins_data = len(abs_logr_sample_no_na)
        percent_above_threshold = num_bins_above_threshold / num_bins_data
        #if >= 5% of the bins with data are gains or losses, it is not all-white
        if percent_above_threshold >= 0.05:
            non_white_checked_samples.add(sample)
        else:
            all_white_samples.add(sample)
        return all_white_samples, non_white_checked_samples


def get_discordance_by_cohort(cohort_names, data_tables_filepath, pairings_filepath,
                              discordance_threshold1, discordance_threshold2,
                              use_all_pairings, plot_output_path,
                              pairing_discordance_values_filename,
                              white_samples_filename):
    """
    Given a list of cohort names, a path to the logR data, and a
    path to the pairings files, returns num_bins_by_cohort_dict and
    percent_genome_by_cohort_dict.  The keys of these are cohort names.  The
    values are dictionaries where keys are passage comparisons and values are
    lists of discordance values.
    """
    num_bins_by_cohort_dict = dict()
    percent_genome_by_cohort_dict = dict()
    for cohort_name in cohort_names:
        #load in logR information and pairings list
        logr_dict = read_sample_info_file(data_tables_filepath, cohort_name)
        pairings_list = get_pairings_with_pt_list(pairings_filepath,
                                                  cohort_name, use_all_pairings)
        #get dictionaries of num_bins and percent_genome discordant, where
        #keys specify which passages are being compared and values are lists
        #of discordance values
        num_bins, percent_genome = discordance_by_passage_comparison(logr_dict,
                                                                     pairings_list,
                                                                     discordance_threshold1,
                                                                     discordance_threshold2,
                                                                     cohort_name,
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
    medians_file_object.write("#mean of means = " + str(round(mean_of_means, 4)) + "\n")
    medians_file_object.write("#minimum discordance = " + str(round(minimum, 4)) + "\n")
    medians_file_object.write("#maximum discordance = " + str(round(maximum, 4)) + "\n")
    
    sorted_cohort_medians = sorted(cohort_medians)
    for i in range(len(cohorts_sorted_by_median)):
        print(cohorts_sorted_by_median[i], sorted_cohort_medians[i])
        medians_file_object.write(cohorts_sorted_by_median[i] + " " + str(sorted_cohort_medians[i]) + "\n")
    
    median_of_medians = statistics.median(cohort_medians)
    plt.figure(figsize=(12, 6))
    plt.axhline(y = median_of_medians, color = "blue")
    plt.text(1, ylimit, "median of medians = " + str(round(median_of_medians, 4)),
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


def sorted_discordance_plots(data_tables_filepath, pairings_filepath, cohort_names,
                             bin_size_name, bin_size_abbr, num_bins_ylimit,
                             plot_output_path, sample_size_thresholds, discordance_threshold1,
                             discordance_threshold2, medians_filename, use_all_pairings,
                             pairing_discordance_values_filename,
                             white_samples_filename):
    """
    Plots sorted discordance plots, including number of bins discordant
    and fraction genome discordant.
    """
    num_bins_by_cohort_dict, percent_genome_by_cohort_dict = get_discordance_by_cohort(cohort_names,
                                                                                       data_tables_filepath,
                                                                                       pairings_filepath,
                                                                                       discordance_threshold1,
                                                                                       discordance_threshold2,
                                                                                       use_all_pairings,
                                                                                       plot_output_path,
                                                                                       pairing_discordance_values_filename,
                                                                                       white_samples_filename)
    
    num_bins_pt_passage, num_bins_early_late = break_up_pt_and_passage_comparisons(num_bins_by_cohort_dict)
    percent_genome_pt_passage, percent_genome_early_late = break_up_pt_and_passage_comparisons(percent_genome_by_cohort_dict)
    
    f = open(plot_output_path + medians_filename, "a+")
    f.write("###Discordance threshold = {0} or {1} diff\n".format(discordance_threshold1, discordance_threshold2))
    for sample_size_threshold in sample_size_thresholds:
        f.write("\n\n###Sample size threshold = {0}\n".format(sample_size_threshold))
        f.write("\n###Number of " + bin_size_name + " Discordant\n")
        f.write("###Discordance of PT vs. PDX Passage\n")
        plot_discordance_by_cohort(num_bins_pt_passage, cohort_names,
                                   "Number of " + bin_size_name + " Discordant",
                                   "Discordance of PT vs. PDX Passage",
                                   num_bins_ylimit,
                                   "{0}pt_vs_pdx_num_bins_{1}_{2}_or_{3}_diff_threshold_n{4}.png".format(plot_output_path, bin_size_abbr, str(discordance_threshold1), str(discordance_threshold2), str(sample_size_threshold)),
                                   sample_size_threshold, f)
        f.write("\n###Number of " + bin_size_name + " Discordant\n")
        f.write("###Discordance of PDX Passage vs. Later PDX Passage\n")
        plot_discordance_by_cohort(num_bins_early_late, cohort_names,
                                   "Number of " + bin_size_name + " Discordant",
                                   "Discordance of PDX Passage vs. Later PDX Passage",
                                   num_bins_ylimit,
                                   "{0}early_vs_late_num_bins_{1}_{2}_or_{3}_diff_threshold_n{4}.png".format(plot_output_path, bin_size_abbr, str(discordance_threshold1), str(discordance_threshold2), str(sample_size_threshold)),
                                   sample_size_threshold, f)
        f.write("\n###Fraction of " + bin_size_name + " Discordant\n")
        f.write("###Discordance of PT vs. PDX Passage\n")
        plot_discordance_by_cohort(percent_genome_pt_passage, cohort_names,
                                   "Fraction of " + bin_size_name + " Discordant",
                                   "Discordance of PT vs. PDX Passage", 1,
                                   "{0}pt_vs_pdx_percent_genome_{1}_{2}_or_{3}_diff_threshold_n{4}.png".format(plot_output_path, bin_size_abbr, str(discordance_threshold1), str(discordance_threshold2), str(sample_size_threshold)),
                                   sample_size_threshold, f)
        f.write("\n###Fraction of " + bin_size_name + " Discordant\n")
        f.write("###Discordance of PDX Passage vs. Later PDX Passage\n")
        plot_discordance_by_cohort(percent_genome_early_late, cohort_names,
                                   "Fraction of " + bin_size_name + " Discordant",
                                   "Discordance of PDX Passage vs. Later PDX Passage", 1,
                                   "{0}early_vs_late_percent_genome_{1}_{2}_or_{3}_diff_threshold_n{4}.png".format(plot_output_path, bin_size_abbr, str(discordance_threshold1), str(discordance_threshold2), str(sample_size_threshold)),
                                   sample_size_threshold, f)
    f.write("\n\n\n")
    f.close()


   
#run
if __name__ == "__main__":
    
    for discordance_threshold_pair in DISCORDANCE_THRESHOLD_PAIRS:
        discordance_threshold1 = discordance_threshold_pair[0]
        discordance_threshold2 = discordance_threshold_pair[1]
        print(discordance_threshold1, discordance_threshold2)
        
        sorted_discordance_plots(MB_DATA_TABLES_FILEPATH, FILEPATH_PAIRINGS_FILES,
                                 COHORT_NAMES_SELECT, "1 Mb Bins",
                                 "1mb", 2500, PLOT_OUTPUT_PATH, SAMPLE_SIZE_THRESHOLDS,
                                 discordance_threshold1, discordance_threshold2,
                                 MEDIANS_FILENAME, USE_ALL_PAIRINGS,
                                 PAIRING_DISCORDANCE_VALUES_FILENAME,
                                 WHITE_SAMPLES_FILENAME)