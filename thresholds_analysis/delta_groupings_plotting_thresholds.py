#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Delta Groupings Plotting Thresholds
July 8th, 2020
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
import statistics


#working directory
#os.chdir("/Volumes/fh")
os.chdir("/Volumes/ha_g")
#os.chdir("/Volumes/ha_g-1")
#os.chdir("/Volumes/ha_g-2")


#settings
"""
--DISCORDANCE THRESHOLD PAIRS - A list of lists with the discordance thresholds
to use.  For example, [0.1, 0.4] means if logR_sample2 >= 0.3 and
logR_sample_1 < 0.3, for these samples to be discordant, logR_sample1 must
additionally be < 0.1 or be different by at least 0.4 from logR_sample2.
(Note that using >1 threshold pair will cause discordance values for each
threshold pair to be written to PAIRING_DISCORDANCE_VALUES_FILENAME in sequence.)
--PLOT_OUTPUT_PATH - Path to folder where output plots and files should go.
--PAIRING_DISCORDANCE_VALUES_FILENAME - Name of output file listing discordance
values for each analyzed pairing.
--WHITE_SAMPLES_FILENAME - Name of output file listing samples which had <5%
of bins with abs(logR) > 0.3 and were thus not included in the discordance
analysis.
--P_VAL_FILENAME - Name of output to write p values for various statistical
tests to.
"""
DISCORDANCE_THRESHOLD_PAIRS = [[0.1, 0.3]]
PLOT_OUTPUT_PATH = "projects/PDX/production/thresholds/plots/delta_groupings/"
P_VAL_FILENAME = "delta_groupings_p_vals.txt"
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


def discordance_by_delta(logr_dict, pairings_list, discordance_threshold1,
                         discordance_threshold2, plot_output_path,
                         pairing_discordance_values_filename,
                         white_samples_filename):
    """
    Given a logR dict and and pairings list, returns two dictionaries
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
                    pair2_passage = get_passage(pair2)
                    delta = int(pair2_passage[1:]) - int(pair1_passage[1:])
                    is_discordant_list = get_discordance_list(logr_dict, pair1, pair2,
                                                              discordance_threshold1,
                                                              discordance_threshold2)
                    num_bins_discordant = is_discordant_list.count("1")
                    #percent genome discordant is out of bins where both samples have logR values
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
                                              discordance_threshold1,
                                              discordance_threshold2,
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
    
    plot_filename = "{0}discordance_by_delta_1mb_{1}_{2}_or_{3}_diff_threshold.png".format(plot_output_path, grouping_name, discordance_threshold1, discordance_threshold2)
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
    data_tables_filepath = MB_DATA_TABLES_FILEPATH
    pairings_filepath = FILEPATH_PAIRINGS_FILES
    plot_output_path = PLOT_OUTPUT_PATH
    p_val_filepath = plot_output_path + P_VAL_FILENAME
    pairing_discordance_values_filename = PAIRING_DISCORDANCE_VALUES_FILENAME
    white_samples_filename = WHITE_SAMPLES_FILENAME
    
    with open(plot_output_path + pairing_discordance_values_filename, "w") as f:
        f.write("cohort\tsample1\tsample2\tdelta\tdiscordance\n")
    
    
    f = open(p_val_filepath, "w")
    
    for discordance_threshold_pair in DISCORDANCE_THRESHOLD_PAIRS:
        discordance_threshold1 = discordance_threshold_pair[0]
        discordance_threshold2 = discordance_threshold_pair[1]
        print(discordance_threshold1, discordance_threshold2)
        f.write("###Discordance threshold = {0} or {1} diff\n".format(discordance_threshold1, discordance_threshold2))
                
        all_cohorts_percent_genome_dict = dict()
        for cohort_name in cohort_names:
            print(cohort_name)
            #load in logR information and pairings list
            logr_dict = read_sample_info_file(data_tables_filepath,
                                              cohort_name)
            pairings_list = get_pairings_with_pt_list(pairings_filepath,
                                                      cohort_name)
            #get dictionaries of num_bins and percent_genome discordant, where
            #keys specify passage delta and values are lists of discordance values
            num_bins, percent_genome = discordance_by_delta(logr_dict,
                                                            pairings_list,
                                                            discordance_threshold1,
                                                            discordance_threshold2,
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
                                                  plot_output_path, discordance_threshold1,
                                                  discordance_threshold2, "grouping1")
        discordance_by_passage_number_change_plot(grouped_dict2, False,
                                                  plot_output_path, discordance_threshold1,
                                                  discordance_threshold2, "grouping2")
        discordance_by_passage_number_change_plot(grouped_dict3, False,
                                                  plot_output_path, discordance_threshold1,
                                                  discordance_threshold2, "grouping3")
    
    f.close()
