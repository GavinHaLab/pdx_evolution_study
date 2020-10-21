#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGS Plotting Ichor Minimal Subclone Discordance
July 22nd, 2020
Anna Hoge
Ha Lab
Fred Hutchinson Cancer Research Center
"""

#imports
import os
import matplotlib.pyplot as plt
import scipy.stats
import glob
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
--P_VAL_FILENAME - Name of output to write p values for various statistical
tests to.
--PAIRING_DISCORDANCE_VALUES_FILENAME - Name of output file listing discordance
values for each analyzed pairing.
--WHITE_SAMPLES_FILENAME - Name of output file listing samples which had <5%
of bins with ploidy-adjusted CN != 2 and/or < 5% tumor fraction and were thus
not included in the discordance analysis.
"""
PLOT_OUTPUT_PATH = "projects/PDX/production/ichor_minimal_subclone_discordance/plots/wgs/"
ICHOR_RESULTS_FILEPATH = "projects/PDX/ichor/run_7_9_GHtest/results/"
ADDITIONAL_DISCORDANCE_THRESHOLD = 0.5
SUBCLONE_THRESHOLD = 0.9
USE_ROUNDED_CN = True
USE_SUBCLONE_FILTER = True
P_VAL_FILENAME = "wgs_p_vals.txt"
PAIRING_DISCORDANCE_VALUES_FILENAME = "pairing_discordance_values.txt"
WHITE_SAMPLES_FILENAME = "white_samples.txt"


#global variables
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


def get_percent_genome_discordant(cn_dict, pair1, pair2,
                                  additional_discordance_threshold):
    """
    Given a CN dictionary, 2 sample names, and discordance threshold,
    returns the fraction of the genome discordant between the two samples.
    """
    is_discordant_list = get_discordance_list(cn_dict, pair1, pair2,
                                              additional_discordance_threshold)
    num_bins_discordant = is_discordant_list.count("1")
    percent_genome_discordant = num_bins_discordant / (num_bins_discordant + is_discordant_list.count("0"))
    return percent_genome_discordant
    
    
def is_white_sample(cn_dict, sample, cohort_name, ichor_results_filepath,
                    plot_output_path, white_samples_filename):
    """
    Given a CN dict and a sample name, returns a Boolean indicating whether
    the sample has < 5% tumor fraction and/or <5% of bins
    with data altered from neutral.
    """
    sample_params_file =  ichor_results_filepath + cohort_name + "/" + sample + "/" + sample + ".params.txt"
    with open(sample_params_file, "r") as f:
        f.readline()
        info = f.readline()
        tumor_fraction = float(info.strip().split()[1])
    #if tfx < 5%, classify as all-white sample
    if tumor_fraction < 0.05:
        #it is all-white
        with open(plot_output_path + white_samples_filename, "a") as f:
            f.write(cohort_name + "\t" + sample + "\n")
        return True
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
            #it is all-white
            with open(plot_output_path + white_samples_filename, "a") as f:
                f.write(cohort_name + "\t" + sample + "\n")
            return True
        #if neither tfx nor altered bin thresholds capture sample, it is not all-white
        else:
            return False
        
        
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


def get_cohort_discordances(pairings_list, cn_dict,
                            additional_discordance_threshold,
                            cohort_name, ichor_results_filepath,
                            plot_output_path,
                            pairing_discordance_values_filename,
                            white_samples_filename):
    """
    Given a pairings list and a CN dict for a cohort, as well as a discordance
    threshold value, returns a dictionary "discordances" where keys are
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
            pt_is_white = is_white_sample(cn_dict, pt, cohort_name,
                                          ichor_results_filepath,
                                          plot_output_path,
                                          white_samples_filename)
            early_is_white = is_white_sample(cn_dict, early, cohort_name,
                                             ichor_results_filepath,
                                             plot_output_path,
                                             white_samples_filename)
            late_is_white = is_white_sample(cn_dict, late, cohort_name,
                                            ichor_results_filepath,
                                            plot_output_path,
                                            white_samples_filename)
            #if none of them are all-white, analyze the line
            if (pt_is_white == False) and (early_is_white == False) and (late_is_white == False):
                sample = pt.split("_")[0]
                pt_early = get_percent_genome_discordant(cn_dict, pt, early,
                                                         additional_discordance_threshold)
                pt_late = get_percent_genome_discordant(cn_dict, pt, late,
                                                        additional_discordance_threshold)
                early_late = get_percent_genome_discordant(cn_dict, early, late,
                                                           additional_discordance_threshold)
                discordances[sample] = [pt_early, pt_late, early_late]
                f.write(cohort_name + "\t" + sample + "\t" + str(pt_early) + "\t" + str(pt_late) + "\t" + str(early_late) + "\n")
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
              additional_discordance_threshold):
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
    plt.savefig("{0}wgs_comparison1_with_lines_{1}_diff_threshold_{2}.png".format(plot_output_path, additional_discordance_threshold, group_abbr))
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
    plt.savefig("{0}wgs_comparison2_with_lines_{1}_diff_threshold_{2}.png".format(plot_output_path, additional_discordance_threshold, group_abbr))
    plt.show()
    
    #without lines connecting samples
    plt.boxplot([pt_early_all, pt_late_all], labels = ["PT:early PDX", "PT:later PDX"])
    plt.xticks(fontsize = 14)
    plt.ylim(-0.1, 1)
    plt.ylabel("Fraction 1 Mb Bins Discordant", fontsize = 16)
    plt.title(title, fontsize = 16)
    plt.savefig("{0}wgs_comparison1_{1}_diff_threshold_{2}.png".format(plot_output_path, additional_discordance_threshold, group_abbr))
    plt.show()
    
    plt.boxplot([pt_early_all, early_late_all], labels = ["PT:early PDX", "early PDX:later PDX"])
    plt.xticks(fontsize = 14)
    plt.ylim(-0.1, 1)
    plt.ylabel("Fraction 1 Mb Bins Discordant", fontsize = 16)
    plt.title(title, fontsize = 16)
    plt.savefig("{0}wgs_comparison2_{1}_diff_threshold_{2}.png".format(plot_output_path, additional_discordance_threshold, group_abbr))
    plt.show()
    
    p_val_file_object.write("PT:early PDX vs PT:later PDX (alternative = 'less')\n")
    wilcoxon(pt_early_all, pt_late_all, p_val_file_object, "less")
    p_val_file_object.write("\nPT:early PDX vs early PDX: later PDX (alternative = 'greater')\n")
    wilcoxon(pt_early_all, early_late_all, p_val_file_object, "greater")


#run
if __name__ == "__main__":
    ichor_results_filepath = ICHOR_RESULTS_FILEPATH
    subclone_threshold = SUBCLONE_THRESHOLD
    use_subclone_filter = USE_SUBCLONE_FILTER
    additional_discordance_threshold = ADDITIONAL_DISCORDANCE_THRESHOLD
    use_rounded_cn = USE_ROUNDED_CN
    plot_output_path = PLOT_OUTPUT_PATH
    p_val_filename = P_VAL_FILENAME
    pairing_discordance_values_filename = PAIRING_DISCORDANCE_VALUES_FILENAME
    white_samples_filename = WHITE_SAMPLES_FILENAME
    
    p_val_file_object = open(plot_output_path + p_val_filename, "w")
    all_discordances = dict()
    p_val_file_object.write("###Additional discordance threshold = {0} diff, subclone filtering = {1}, subclone threshold = {2}\n".format(additional_discordance_threshold, str(use_subclone_filter), str(subclone_threshold)))
    for cohort_name in WGS_COHORT_NAMES:
        if cohort_name == "europdx_wgs_brca":
            pairings_list = BRCA_PAIRINGS
            group_abbr = "brca"
        elif cohort_name == "europdx_wgs_crc":
            pairings_list = CRC_PAIRINGS
            group_abbr = "crc"
        cn_dict = read_ichor_results(ichor_results_filepath, cohort_name, use_rounded_cn,
                   use_subclone_filter, subclone_threshold)
        discordances = get_cohort_discordances(pairings_list, cn_dict,
                                               additional_discordance_threshold,
                                               cohort_name, ichor_results_filepath,
                                               plot_output_path,
                                               pairing_discordance_values_filename,
                                               white_samples_filename)
        
        #analyze the given cohort
        p_val_file_object.write("\n###{0}\n".format(cohort_name))
        wgs_plots(discordances, cohort_name, plot_output_path, group_abbr,
                  p_val_file_object, additional_discordance_threshold)
        
        all_discordances.update(discordances)
        
    #analyze both cohorts together
    p_val_file_object.write("\n###{0}\n".format("both europdx wgs cohorts"))
    wgs_plots(all_discordances, "europdx wgs cohorts", plot_output_path,
              "both", p_val_file_object, additional_discordance_threshold)
    p_val_file_object.write("\n\n\n")
    p_val_file_object.close()
    
    
    
    