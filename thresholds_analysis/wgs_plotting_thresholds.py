#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGS Plotting Thresholds
July 8th, 2020
Anna Hoge
Ha Lab
Fred Hutchinson Cancer Research Center
"""

#imports
import os
import matplotlib.pyplot as plt
import scipy.stats


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
--P_VAL_FILENAME - Name of output to write p values for various statistical
tests to.
--PAIRING_DISCORDANCE_VALUES_FILENAME - Name of output file listing discordance
values for each analyzed pairing.
--WHITE_SAMPLES_FILENAME - Name of output file listing samples which had <5%
of bins with abs(logR) > 0.3 and were thus not included in the discordance
analysis.
"""
DISCORDANCE_THRESHOLD_PAIRS = [[0.1, 0.3]]
PLOT_OUTPUT_PATH = "projects/PDX/production/thresholds/plots/wgs/"
P_VAL_FILENAME = "wgs_p_vals.txt"
PAIRING_DISCORDANCE_VALUES_FILENAME = "pairing_discordance_values.txt"
WHITE_SAMPLES_FILENAME = "white_samples.txt"


#global variables
"""
--MB_DATA_TABLES_FILEPATH - Path to folder containing tsv files for each cohort,
with 1 Mb bin logR values for each sample.
"""
MB_DATA_TABLES_FILEPATH = "projects/PDX/1mb_severe_logr_and_cn_tables/"
WGS_COHORT_NAMES = ["europdx_wgs_brca", "europdx_wgs_crc"]
BRCA_PAIRINGS = [['BC111_PT_WGS', 'BC111_P2_WGS', 'BC111_P5_WGS'],
 ['BC126_PT_WGS', 'BC126_P1_WGS', 'BC126_P5_WGS'],
 ['BC127_PT_WGS', 'BC127_P1_WGS', 'BC127_P5_WGS'],
 ['BC162_PT_WGS', 'BC162_P1_WGS', 'BC162_P5_WGS'],
 ['BC201_PT_WGS', 'BC201_P1_WGS', 'BC201_P5_WGS'],
 ['BC227_PT_WGS', 'BC227_P0_WGS', 'BC227_P6_WGS'],
 ['BC230_PT_WGS', 'BC230_P2_WGS', 'BC230_P6_WGS'],
 ['BC239_PT_WGS', 'BC239_P1_WGS', 'BC239_P5_WGS'],
 ['BC241_PT_WGS', 'BC241_P1_WGS', 'BC241_P5_WGS'],
 ['BC250_PT_WGS', 'BC250_P1_WGS', 'BC250_P4_WGS'],
 ['BC256_PT_WGS', 'BC256_P1_WGS', 'BC256_P5_WGS'],
 ['BC274_PT_WGS', 'BC274_P1_WGS', 'BC274_P5_WGS'],
 ['BC283_PT_WGS', 'BC283_P1_WGS', 'BC283_P5_WGS'],
 ['BC288_PT_WGS', 'BC288_P1_WGS', 'BC288_P5_WGS'],
 ['BC291_PT_WGS', 'BC291_P0_WGS', 'BC291_P5_WGS'],
 ['BC291A_PT_WGS', 'BC291A_P1_WGS', 'BC291A_P7_WGS'],
 ['BC297_PT_WGS', 'BC297_P1_WGS', 'BC297_P6_WGS'],
 ['BC302_PT_WGS', 'BC302_P1_WGS', 'BC302_P5_WGS'],
 ['BC325_PT_WGS', 'BC325_P1_WGS', 'BC325_P5_WGS'],
 ['BC326_PT_WGS', 'BC326_P1_WGS', 'BC326_P7_WGS'],
 ['BC332_PT_WGS', 'BC332_P1_WGS', 'BC332_P5_WGS'],
 ['BC336_PT_WGS', 'BC336_P1_WGS', 'BC336_P5_WGS'],
 ['BC343_PT_WGS', 'BC343_P1_WGS', 'BC343_P5_WGS'],
 ['BC408_PT_WGS', 'BC408_P0_WGS', 'BC408_P6_WGS'],
 ['BC426_PT_WGS', 'BC426_P1_WGS', 'BC426_P5_WGS'],
 ['BC430_PT_WGS', 'BC430_P0_WGS', 'BC430_P5_WGS'],
 ['BC453_PT_WGS', 'BC453_P1_WGS', 'BC453_P4_WGS'],
 ['BC457_PT_WGS', 'BC457_P1_WGS', 'BC457_P5_WGS'],
 ['BC464_PT_WGS', 'BC464_P1_WGS', 'BC464_P9_WGS'],
 ['BC474_PT_WGS', 'BC474_P0_WGS', 'BC474_P7_WGS'],
 ['BC509_PT_WGS', 'BC509_P1_WGS', 'BC509_P5_WGS'],
 ['BC513_PT_WGS', 'BC513_P1_WGS', 'BC513_P7_WGS'],
 ['BC519_PT_WGS', 'BC519_P1_WGS', 'BC519_P5_WGS'],
 ['BC531_PT_WGS', 'BC531_P1_WGS', 'BC531_P5_WGS'],
 ['BC666_PT_WGS', 'BC666_P1_WGS', 'BC666_P5_WGS'],
 ['BC701_PT_WGS', 'BC701_P1_WGS', 'BC701_P5_WGS'],
 ['BC715_PT_WGS', 'BC715_P1_WGS', 'BC715_P5_WGS'],
 ['BC794_PT_WGS', 'BC794_P1_WGS', 'BC794_P5_WGS'],
 ['BC801_PT_WGS', 'BC801_P1_WGS', 'BC801_P5_WGS'],
 ['BC878_PT_WGS', 'BC878_P0_WGS', 'BC878_P5_WGS'],
 ['BC899_PT_WGS', 'BC899_P1_WGS', 'BC899_P5_WGS'],
 ['BC989_PT_WGS', 'BC989_P1_WGS', 'BC989_P5_WGS'],
 ['BC1006_PT_WGS', 'BC1006_P1_WGS', 'BC1006_P5_WGS']]

CRC_PAIRINGS = [['CRC1_PT_WGS', 'CRC1_P0_WGS', 'CRC1_P3_WGS'],
 ['CRC2_PT_WGS', 'CRC2_P1_WGS', 'CRC2_P6_WGS'],
 ['CRC3_PT_WGS', 'CRC3_P1_WGS', 'CRC3_P5_WGS'],
 ['CRC4_PT_WGS', 'CRC4_P1_WGS', 'CRC4_P6_WGS'],
 ['CRC5_PT_WGS', 'CRC5_P1_WGS', 'CRC5_P5_WGS'],
 ['CRC7_PT_WGS', 'CRC7_P1_WGS', 'CRC7_P6_WGS'],
 ['CRC8_PT_WGS', 'CRC8_P1_WGS', 'CRC8_P6_WGS'],
 ['CRC9_PT_WGS', 'CRC9_P1_WGS', 'CRC9_P4_WGS'],
 ['CRC10_PT_WGS', 'CRC10_P0_WGS', 'CRC10_P3_WGS'],
 ['CRC11_PT_WGS', 'CRC11_P1_WGS', 'CRC11_P6_WGS'],
 ['CRC13_PT_WGS', 'CRC13_P1_WGS', 'CRC13_P4_WGS'],
 ['CRC14_PT_WGS', 'CRC14_P1_WGS', 'CRC14_P6_WGS'],
 ['CRC15_PT_WGS', 'CRC15_P1_WGS', 'CRC15_P4_WGS'],
 ['CRC16_PT_WGS', 'CRC16_P1_WGS', 'CRC16_P6_WGS'],
 ['CRC17_PT_WGS', 'CRC17_P1_WGS', 'CRC17_P6_WGS'],
 ['CRC18_PT_WGS', 'CRC18_P1_WGS', 'CRC18_P6_WGS'],
 ['CRC20_PT_WGS', 'CRC20_P1_WGS', 'CRC20_P5_WGS'],
 ['CRC21_PT_WGS', 'CRC21_P1_WGS', 'CRC21_P6_WGS'],
 ['CRC22_PT_WGS', 'CRC22_P1_WGS', 'CRC22_P5_WGS'],
 ['CRC23_PT_WGS', 'CRC23_P1_WGS', 'CRC23_P6_WGS'],
 ['CRC24_PT_WGS', 'CRC24_P1_WGS', 'CRC24_P6_WGS'],
 ['CRC25_PT_WGS', 'CRC25_P1_WGS', 'CRC25_P6_WGS'],
 ['CRC27_PT_WGS', 'CRC27_P1_WGS', 'CRC27_P6_WGS'],
 ['CRC28_PT_WGS', 'CRC28_P1_WGS', 'CRC28_P6_WGS'],
 ['CRC29_PT_WGS', 'CRC29_P1_WGS', 'CRC29_P5_WGS'],
 ['CRC30_PT_WGS', 'CRC30_P1_WGS', 'CRC30_P4_WGS'],
 ['CRC32_PT_WGS', 'CRC32_P1_WGS', 'CRC32_P6_WGS'],
 ['CRC33_PT_WGS', 'CRC33_P1_WGS', 'CRC33_P5_WGS'],
 ['CRC34_PT_WGS', 'CRC34_P0_WGS', 'CRC34_P3_WGS'],
 ['CRC35_PT_WGS', 'CRC35_P1_WGS', 'CRC35_P6_WGS'],
 ['CRC36_PT_WGS', 'CRC36_P0_WGS', 'CRC36_P2_WGS'],
 ['CRC37_PT_WGS', 'CRC37_P1_WGS', 'CRC37_P5_WGS'],
 ['CRC38_PT_WGS', 'CRC38_P1_WGS', 'CRC38_P6_WGS'],
 ['CRC39_PT_WGS', 'CRC39_P1_WGS', 'CRC39_P6_WGS'],
 ['CRC40_PT_WGS', 'CRC40_P0_WGS', 'CRC40_P2_WGS'],
 ['CRC41_PT_WGS', 'CRC41_P1_WGS', 'CRC41_P4_WGS'],
 ['CRC42_PT_WGS', 'CRC42_P1_WGS', 'CRC42_P5_WGS'],
 ['CRC43_PT_WGS', 'CRC43_P1_WGS', 'CRC43_P5_WGS'],
 ['CRC44_PT_WGS', 'CRC44_P1_WGS', 'CRC44_P6_WGS'],
 ['CRC45_PT_WGS', 'CRC45_P1_WGS', 'CRC45_P5_WGS'],
 ['CRC46_PT_WGS', 'CRC46_P1_WGS', 'CRC46_P6_WGS'],
 ['CRC47_PT_WGS', 'CRC47_P1_WGS', 'CRC47_P5_WGS'],
 ['CRC48_PT_WGS', 'CRC48_P1_WGS', 'CRC48_P6_WGS'],
 ['CRC49_PT_WGS', 'CRC49_P1_WGS', 'CRC49_P5_WGS'],
 ['CRC50_PT_WGS', 'CRC50_P0_WGS', 'CRC50_P2_WGS'],
 ['CRC51_PT_WGS', 'CRC51_P1_WGS', 'CRC51_P4_WGS'],
 ['CRC52_PT_WGS', 'CRC52_P1_WGS', 'CRC52_P6_WGS'],
 ['CRC53_PT_WGS', 'CRC53_P1_WGS', 'CRC53_P6_WGS'],
 ['CRC54_PT_WGS', 'CRC54_P1_WGS', 'CRC54_P6_WGS'],
 ['CRC55_PT_WGS', 'CRC55_P1_WGS', 'CRC55_P5_WGS'],
 ['CRC56_PT_WGS', 'CRC56_P1_WGS', 'CRC56_P5_WGS'],
 ['CRC57_PT_WGS', 'CRC57_P1_WGS', 'CRC57_P6_WGS'],
 ['CRC58_PT_WGS', 'CRC58_P1_WGS', 'CRC58_P6_WGS'],
 ['CRC59_PT_WGS', 'CRC59_P1_WGS', 'CRC59_P4_WGS'],
 ['CRC60_PT_WGS', 'CRC60_P1_WGS', 'CRC60_P4_WGS'],
 ['CRC61_PT_WGS', 'CRC61_P1_WGS', 'CRC61_P6_WGS'],
 ['CRC62_PT_WGS', 'CRC62_P1_WGS', 'CRC62_P6_WGS'],
 ['CRC63_PT_WGS', 'CRC63_P1_WGS', 'CRC63_P6_WGS'],
 ['CRC64_PT_WGS', 'CRC64_P0_WGS', 'CRC64_P4_WGS'],
 ['CRC65_PT_WGS', 'CRC65_P0_WGS', 'CRC65_P3_WGS'],
 ['CRC66_PT_WGS', 'CRC66_P1_WGS', 'CRC66_P5_WGS'],
 ['CRC67_PT_WGS', 'CRC67_P0_WGS', 'CRC67_P2_WGS'],
 ['CRC68_PT_WGS', 'CRC68_P1_WGS', 'CRC68_P6_WGS'],
 ['CRC69_PT_WGS', 'CRC69_P1_WGS', 'CRC69_P6_WGS'],
 ['CRC70_PT_WGS', 'CRC70_P1_WGS', 'CRC70_P4_WGS'],
 ['CRC71_PT_WGS', 'CRC71_P1_WGS', 'CRC71_P6_WGS'],
 ['CRC72_PT_WGS', 'CRC72_P1_WGS', 'CRC72_P4_WGS'],
 ['CRC73_PT_WGS', 'CRC73_P1_WGS', 'CRC73_P4_WGS'],
 ['CRC76_PT_WGS', 'CRC76_P1_WGS', 'CRC76_P4_WGS'],
 ['CRC77_PT_WGS', 'CRC77_P0_WGS', 'CRC77_P3_WGS'],
 ['CRC78_PT_WGS', 'CRC78_P1_WGS', 'CRC78_P4_WGS'],
 ['CRC79_PT_WGS', 'CRC79_P1_WGS', 'CRC79_P4_WGS'],
 ['CRC80_PT_WGS', 'CRC80_P1_WGS', 'CRC80_P6_WGS'],
 ['CRC81_PT_WGS', 'CRC81_P1_WGS', 'CRC81_P4_WGS'],
 ['CRC82_PT_WGS', 'CRC82_P1_WGS', 'CRC82_P4_WGS'],
 ['CRC83_PT_WGS', 'CRC83_P0_WGS', 'CRC83_P3_WGS'],
 ['CRC84_PT_WGS', 'CRC84_P0_WGS', 'CRC84_P3_WGS'],
 ['CRC85_PT_WGS', 'CRC85_P0_WGS', 'CRC85_P3_WGS'],
 ['CRC86_PT_WGS', 'CRC86_P1_WGS', 'CRC86_P6_WGS'],
 ['CRC87_PT_WGS', 'CRC87_P1_WGS', 'CRC87_P4_WGS'],
 ['CRC88_PT_WGS', 'CRC88_P0_WGS', 'CRC88_P3_WGS'],
 ['CRC89_PT_WGS', 'CRC89_P0_WGS', 'CRC89_P3_WGS'],
 ['CRC90_PT_WGS', 'CRC90_P0_WGS', 'CRC90_P3_WGS'],
 ['CRC91_PT_WGS', 'CRC91_P1_WGS', 'CRC91_P4_WGS'],
 ['CRC92_PT_WGS', 'CRC92_P1_WGS', 'CRC92_P4_WGS'],
 ['CRC93_PT_WGS', 'CRC93_P1_WGS', 'CRC93_P4_WGS'],
 ['CRC94_PT_WGS', 'CRC94_P1_WGS', 'CRC94_P3_WGS'],
 ['CRC95_PT_WGS', 'CRC95_P1_WGS', 'CRC95_P4_WGS'],
 ['CRC96_PT_WGS', 'CRC96_P1_WGS', 'CRC96_P4_WGS'],
 ['CRC97_PT_WGS', 'CRC97_P1_WGS', 'CRC97_P4_WGS'],
 ['CRC98_PT_WGS', 'CRC98_P0_WGS', 'CRC98_P3_WGS'],
 ['CRC99_PT_WGS', 'CRC99_P0_WGS', 'CRC99_P3_WGS'],
 ['CRC100_PT_WGS', 'CRC100_P0_WGS', 'CRC100_P3_WGS']]


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


def get_percent_genome_discordant(logr_dict, pair1, pair2, discordance_threshold1,
                                  discordance_threshold2):
    """
    Given a logR dictionary, 2 sample names, and discordance threshold,
    returns the fraction of the genome discordant between the two samples.
    """
    is_discordant_list = get_discordance_list(logr_dict, pair1, pair2,
                                              discordance_threshold1,
                                              discordance_threshold2)
    num_bins_discordant = is_discordant_list.count("1")
    percent_genome_discordant = num_bins_discordant / (num_bins_discordant + is_discordant_list.count("0"))
    return percent_genome_discordant


def is_white_sample(logr_dict, sample, cohort_abbr, plot_output_path, white_samples_filename):
    """
    Given a logR dict and a sample name, returns a Boolean indicating whether
    the sample has <5% altered bins.
    """
    logr_sample = logr_dict[sample]
    abs_logr_sample_no_na = [abs(float(logr)) for logr in logr_sample if logr != "NA"]
    num_bins_above_threshold = len([logr for logr in abs_logr_sample_no_na if logr >= 0.3])
    num_bins_data = len(abs_logr_sample_no_na)
    percent_above_threshold = num_bins_above_threshold / num_bins_data
    #if >= 5% of the bins with data are gains or losses, it is not all-white
    if percent_above_threshold >= 0.05:
        return False
    else:
        with open(plot_output_path + white_samples_filename, "a") as f:
            f.write(cohort_abbr + "\t" + sample + "\n")
        return True


def get_cohort_discordances(pairings_list, logr_dict,
                            discordance_threshold1,
                            discordance_threshold2,
                            cohort_abbr, plot_output_path,
                            pairing_discordance_values_filename,
                            white_samples_filename):
    """
    Given a pairings list and a logR dict for a cohort, as well as discordance
    threshold values, returns a dictionary "discordances" where keys are
    line IDs and values are lists of the discordance between [PT vs early PDX,
    PT vs later PDX, and early PDX vs later PDX] for each pairing in pairings_list.
    Also writes the all-white samples to the given filename.
    """
    discordances = dict()
    with open(plot_output_path + pairing_discordance_values_filename, "a") as f:
        for pairings in pairings_list:
            #check if any samples are all-white
            pt = pairings[0]
            early = pairings[1]
            late = pairings[2]
            pt_is_white = is_white_sample(logr_dict, pt, cohort_abbr,
                                          plot_output_path, white_samples_filename)
            early_is_white = is_white_sample(logr_dict, early, cohort_abbr,
                                             plot_output_path, white_samples_filename)
            late_is_white = is_white_sample(logr_dict, late, cohort_abbr,
                                            plot_output_path, white_samples_filename)
            #if none of them are all-white, analyze the line
            if (pt_is_white == False) and (early_is_white == False) and (late_is_white == False):
                sample = pt.split("_")[0]
                pt_early = get_percent_genome_discordant(logr_dict, pt, early,
                                                         discordance_threshold1,
                                                         discordance_threshold2)
                pt_late = get_percent_genome_discordant(logr_dict, pt, late,
                                                        discordance_threshold1,
                                                        discordance_threshold2)
                early_late = get_percent_genome_discordant(logr_dict, early, late,
                                                           discordance_threshold1,
                                                           discordance_threshold2)
                discordances[sample] = [pt_early, pt_late, early_late]
                f.write(cohort_abbr + "\t" + sample + "\t" + str(pt_early) + "\t" + str(pt_late) + "\t" + str(early_late) + "\n")
    return discordances


def wilcoxon(list1, list2, p_val_file_object, alternative_hypothesis):
    """
    Given two lists of values and an alternative hypothesis, computes the
    Wilcoxon signed-rank test and writes the results to the given file object.
    """
    test_stat, p_value = scipy.stats.wilcoxon(list1, list2,
                                                  alternative = alternative_hypothesis)
    rounded_p_value = str(round(p_value, 5))
    rounded_test_stat = str(round(test_stat, 5))
    if p_value < 0.05:
        print("different with p value = " + rounded_p_value)
        p_val_file_object.write("different with p value = {0}.  test stat = {1}.  n = {2}\n".format(rounded_p_value, rounded_test_stat, str(len(list1))))
    else:
        print("no difference, p value = " + rounded_p_value)
        p_val_file_object.write("no difference, p value = {0}.  test stat = {1}.  n = {2}\n".format(rounded_p_value, rounded_test_stat, str(len(list1))))


def wgs_plots(discordances, title, plot_output_path, group_abbr, p_val_file_object,
              discordance_threshold1, discordance_threshold2):
    """
    Given a dictionary of discordances as described in get_cohort_discordances(),
    plot those discordances and write the p values of statistical tests to the
    given file object.
    """
    pt_early_all = []
    pt_late_all = []
    early_late_all = []
    
    #with lines connecting samples
    for sample in discordances:
        pt_early = discordances[sample][0]
        pt_late = discordances[sample][1]
        early_late = discordances[sample][2]
        pt_early_all.append(pt_early)
        pt_late_all.append(pt_late)
        early_late_all.append(early_late)
        plt.plot([1, 2], [pt_early, pt_late], marker="o", color = "orange",
                 alpha = 0.2)
    plt.boxplot([pt_early_all, pt_late_all], labels = ["PT:early PDX", "PT:later PDX"])
    plt.xticks(fontsize = 14)
    plt.ylim(-0.1, 1)
    plt.ylabel("Fraction 1 Mb Bins Discordant", fontsize = 16)
    plt.title(title, fontsize = 16)
    plt.savefig("{0}wgs_comparison1_with_lines_{1}_or_{2}_diff_threshold_{3}.png".format(plot_output_path, discordance_threshold1, discordance_threshold2, group_abbr))
    plt.show()
    
    for sample in discordances:
        pt_early = discordances[sample][0]
        early_late = discordances[sample][2]
        plt.plot([1, 2], [pt_early, early_late], marker="o", color = "orange",
                 alpha = 0.2)
    plt.boxplot([pt_early_all, early_late_all], labels = ["PT:early PDX", "early PDX:later PDX"])
    plt.xticks(fontsize = 14)
    plt.ylim(-0.1, 1)
    plt.ylabel("Fraction 1 Mb Bins Discordant", fontsize = 16)
    plt.title(title, fontsize = 16)
    plt.savefig("{0}wgs_comparison2_with_lines_{1}_or_{2}_diff_threshold_{3}.png".format(plot_output_path, discordance_threshold1, discordance_threshold2, group_abbr))
    plt.show()
    
    #without lines connecting samples
    plt.boxplot([pt_early_all, pt_late_all], labels = ["PT:early PDX", "PT:later PDX"])
    plt.xticks(fontsize = 14)
    plt.ylim(-0.1, 1)
    plt.ylabel("Fraction 1 Mb Bins Discordant", fontsize = 16)
    plt.title(title, fontsize = 16)
    plt.savefig("{0}wgs_comparison1_{1}_or_{2}_diff_threshold_{3}.png".format(plot_output_path, discordance_threshold1, discordance_threshold2, group_abbr))
    plt.show()
    
    plt.boxplot([pt_early_all, early_late_all], labels = ["PT:early PDX", "early PDX:later PDX"])
    plt.xticks(fontsize = 14)
    plt.ylim(-0.1, 1)
    plt.ylabel("Fraction 1 Mb Bins Discordant", fontsize = 16)
    plt.title(title, fontsize = 16)
    plt.savefig("{0}wgs_comparison2_{1}_or_{2}_diff_threshold_{3}.png".format(plot_output_path, discordance_threshold1, discordance_threshold2, group_abbr))
    plt.show()
    
    p_val_file_object.write("PT:early PDX vs PT:later PDX (alternative = 'less')\n")
    wilcoxon(pt_early_all, pt_late_all, p_val_file_object, "less")
    p_val_file_object.write("\nPT:early PDX vs early PDX: later PDX (alternative = 'greater')\n")
    wilcoxon(pt_early_all, early_late_all, p_val_file_object, "greater")



#run
if __name__ == "__main__":
    p_val_file_object = open(PLOT_OUTPUT_PATH + P_VAL_FILENAME, "w")
    with open(PLOT_OUTPUT_PATH + PAIRING_DISCORDANCE_VALUES_FILENAME, "w") as f:
        f.write("cohort\tsample\tpt_early\tpt_late\tearly_late\n")
    for discordance_threshold_pair in DISCORDANCE_THRESHOLD_PAIRS:
        discordance_threshold1 = discordance_threshold_pair[0]
        discordance_threshold2 = discordance_threshold_pair[1]
        print(discordance_threshold1, discordance_threshold2)
        
        all_discordances = dict()
        p_val_file_object.write("###Discordance threshold = {0} or {1} diff\n\n\n".format(discordance_threshold1, discordance_threshold2))

        for cohort_name in WGS_COHORT_NAMES:
            if cohort_name == "europdx_wgs_brca":
                pairings_list = BRCA_PAIRINGS
                group_abbr = "brca"
            elif cohort_name == "europdx_wgs_crc":
                pairings_list = CRC_PAIRINGS
                group_abbr = "crc"
            logr_dict = read_sample_info_file(MB_DATA_TABLES_FILEPATH, cohort_name)
            discordances = get_cohort_discordances(pairings_list, logr_dict,
                                                   discordance_threshold1,
                                                   discordance_threshold2,
                                                   group_abbr, PLOT_OUTPUT_PATH,
                                                   PAIRING_DISCORDANCE_VALUES_FILENAME,
                                                   WHITE_SAMPLES_FILENAME)
            
            #analyze the given cohort
            p_val_file_object.write("\n###{0}\n".format(cohort_name))
            wgs_plots(discordances, cohort_name, PLOT_OUTPUT_PATH, group_abbr,
                      p_val_file_object, discordance_threshold1, discordance_threshold2)
            
            
            all_discordances.update(discordances)
            
        #analyze both cohorts together
        p_val_file_object.write("\n###{0}\n".format("both europdx wgs cohorts"))
        wgs_plots(all_discordances, "europdx wgs cohorts", PLOT_OUTPUT_PATH,
                  "both", p_val_file_object)
        p_val_file_object.write("\n\n\n")
    p_val_file_object.close()







