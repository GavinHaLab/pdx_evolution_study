#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ploidy and Tfx Plotting Ichor
July 22nd, 2020
Anna Hoge
Ha Lab
Fred Hutchinson Cancer Research Center
"""

#imports
import glob
import os
import matplotlib.pyplot as plt


#working directory
os.chdir("/Volumes/ha_g/")


#settings
OUTPUT_DIR = "projects/PDX/production/ichor_minimal_subclone_discordance/plots/ploidy_and_tfx/"
ICHOR_RESULTS_FILEPATH = "projects/PDX/ichor/run_7_9_GHtest/results/"
ALL_WHITE_FILE = "projects/PDX/production/ichor_minimal_subclone_discordance/plots/1mb_discordance/white_samples.txt"
USE_ALL_PAIRINGS = True


#global variables
OUTPUT_FILENAME = "pairings_tfx_ploidy_comparisons.txt"
EXCL_WHITE_OUTPUT_FILENAME = "pairings_tfx_ploidy_comparisons_excl_all_white.txt"
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
def get_pairings_with_pt_list(pairings_files_filepath, cohort_name, use_all_pairings):
    """
    Given a pairings files filepath and a cohort name, reads the pairings file
    for that cohort and returns a list of pairings where each element in the
    list is a 3-element list itself specifying the first and second sample in
    the pairing and the PT those samples came from ("NA" if no PT sample).
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


def get_tfxs_and_ploidies(cohort_names, ichor_results_filepath):
    """
    Given a list of cohort names and a filepath to the ichor results for these
    cohorts, returns a dictionary with tumor fractions and a dictionary with
    ploidies.  Each dictionary's keys are cohort names and values are dictionaries
    themselves with sample name keys and values giving tumor fraction or ploidy.
    """
    tumor_fractions = dict()
    ploidies = dict()
    for cohort in cohort_names:
        tumor_fractions[cohort] = dict()
        ploidies[cohort] = dict()
        params_files = glob.glob(ichor_results_filepath + cohort + "/*/*.params.txt")
        for param_file in params_files:
            with open(param_file, "r") as f:
                f.readline()
                info = f.readline().strip().split()
                sample_name = info[0]
                ploidy = float(info[2])
                tumor_fraction = float(info[1])
                tumor_fractions[cohort][sample_name] = tumor_fraction
                ploidies[cohort][sample_name] = ploidy
    return tumor_fractions, ploidies


def separate_by_pt_pdx(sample_dictionary):
    """
    Given a dictionary specifying tumor fraction or ploidy values for each sample
    in each cohort, returns a list of the values from PT samples a list of the
    values from PDX samples
    """
    pt = []
    pdx = []
    for cohort in sample_dictionary:
        for sample in sample_dictionary[cohort]:
            if "PT" in sample:
                pt.append(sample_dictionary[cohort][sample])
            else:
                pdx.append(sample_dictionary[cohort][sample])
    return pt, pdx


def plot_ploidies_and_tumor_fractions(pt_tumor_fractions, pdx_tumor_fractions,
                                      pt_ploidies, pdx_ploidies, output_dir):
    """
    Given lists of pt_tumor_fractions, pdx_tumor_fractions, pt_ploidies, and
    pdx_ploidies, plots histograms of each and saves the results to the given
    output directory.
    """
    plt.hist(pt_tumor_fractions, bins = 50)
    plt.title("PT tumor fractions")
    plt.savefig(output_dir + "pt_tumor_fractions_hist.png")
    plt.show()
    plt.hist(pdx_tumor_fractions, bins = 50)
    plt.title("PDX tumor fractions")
    plt.savefig(output_dir + "pdx_tumor_fractions_hist.png")
    plt.show()
    plt.hist(pt_ploidies, bins = 50)
    plt.title("PT ploidies")
    plt.savefig(output_dir + "pt_ploidies_hist.png")
    plt.show()
    plt.hist(pdx_ploidies, bins = 50)
    plt.title("PDX ploidies")
    plt.savefig(output_dir + "pdx_ploidies_hist.png")
    plt.show()
    
    
def pairings_comparison(cohort_names, filepath_pairings_files,
                        use_all_pairings, tumor_fractions_dictionary,
                        ploidies_dictionary, output_dir, output_filename,
                        white_samples, excl_white_output_filename):
    """
    Given a list of cohort names, the filepath where the pairings files are stored,
    the option to use all pairings or direct lineage pairings, and dictionaries
    for tumor fractions and ploidies for each sample in each cohort, writes a file
    and draws a plot to the given output_dir comparing tumor fractions and ploidies
    for paired samples.
    """
    first_ploidies = []
    second_ploidies = []
    with open(output_dir + output_filename, "w") as f:
        with open(output_dir + excl_white_output_filename, "w") as excl_white_f:
            f.write("cohort\tsample1\tsample2\tsample1_tfx\tsample2_tfx\ttfx_diff\tsample1_ploidy\tsample2_ploidy\tploidy_diff\n")
            excl_white_f.write("cohort\tsample1\tsample2\tsample1_tfx\tsample2_tfx\ttfx_diff\tsample1_ploidy\tsample2_ploidy\tploidy_diff\n")
            for cohort in cohort_names:
                pairings = get_pairings_with_pt_list(filepath_pairings_files, cohort,
                                                     use_all_pairings)
                for pairing in pairings:
                    pair1 = pairing[0]
                    pair2 = pairing[1]
                    pair1_tumor_fraction = tumor_fractions_dictionary[cohort][pair1]
                    pair2_tumor_fraction = tumor_fractions_dictionary[cohort][pair2]
                    pair1_ploidy = ploidies_dictionary[cohort][pair1]
                    pair2_ploidy = ploidies_dictionary[cohort][pair2]
                    difference_in_tumor_fraction = round(abs(pair1_tumor_fraction - pair2_tumor_fraction), 3)
                    difference_in_ploidy = round(abs(pair1_ploidy - pair2_ploidy), 3)
                    f.write("\t".join((cohort, pair1, pair2, str(pair1_tumor_fraction),
                                      str(pair2_tumor_fraction),
                                      str(difference_in_tumor_fraction),
                                      str(pair1_ploidy), str(pair2_ploidy),
                                      str(difference_in_ploidy))) + "\n")
                    if (pair1 not in white_samples[cohort]) and (pair2 not in white_samples[cohort]):
                        excl_white_f.write("\t".join((cohort, pair1, pair2, str(pair1_tumor_fraction),
                                           str(pair2_tumor_fraction),
                                           str(difference_in_tumor_fraction),
                                           str(pair1_ploidy), str(pair2_ploidy),
                                           str(difference_in_ploidy))) + "\n")
                    first_ploidies.append(pair1_ploidy)
                    second_ploidies.append(pair2_ploidy)
    plt.scatter(first_ploidies, second_ploidies)
    plt.plot([0, 4], [0, 4], color = "red")
    plt.title("Ploidy comparisons of paired samples")
    plt.savefig(output_dir + "pairings_ploidy_comparison.png")
    plt.show()
    
    
def read_all_white_file(all_white_file, cohorts):
    """
    Given an file listing the all-white samples by cohort and a list of cohort
    names to include, return a dictionary where keys are cohorts and values
    are lists of all-white samples.
    """
    white_samples = dict()
    for cohort in cohorts:
        white_samples[cohort] = []
    with open(all_white_file, "r") as f:
        for line in f.readlines():
            line_parts = line.strip().split()
            cohort = line_parts[0]
            sample = line_parts[1]
            white_samples[cohort].append(sample)
    return white_samples


def write_list_to_file(input_list, filename):
    """
    Given an input_list and a filename, write each element of the input_list
    on a separate line to the file.
    """
    with open(filename, "w") as f:
        for element in input_list:
            f.write(str(element) + "\n")



#run
if __name__ == "__main__":
    cohort_names = COHORT_NAMES_SELECT
    ichor_results_filepath = ICHOR_RESULTS_FILEPATH
    filepath_pairings_files = FILEPATH_PAIRINGS_FILES
    use_all_pairings = USE_ALL_PAIRINGS
    output_dir = OUTPUT_DIR
    output_filename = OUTPUT_FILENAME
    excl_white_output_filename = EXCL_WHITE_OUTPUT_FILENAME
    all_white_file = ALL_WHITE_FILE
    
    #get tumor fractions and ploidies
    print("getting tumor fractions and ploidies")
    tumor_fractions, ploidies = get_tfxs_and_ploidies(cohort_names,
                                                      ichor_results_filepath)
    #separate into PT and PDXs
    print("separating by PT and PDX")
    pt_tumor_fractions, pdx_tumor_fractions = separate_by_pt_pdx(tumor_fractions)
    pt_ploidies, pdx_ploidies = separate_by_pt_pdx(ploidies)
    #plot ploidies and tumor fractions
    plot_ploidies_and_tumor_fractions(pt_tumor_fractions, pdx_tumor_fractions,
                                      pt_ploidies, pdx_ploidies, output_dir)
    write_list_to_file(pt_tumor_fractions, output_dir + "pt_tfxs.txt")
    write_list_to_file(pdx_tumor_fractions, output_dir + "pdx_tfxs.txt")
    write_list_to_file(pt_ploidies, output_dir + "pt_ploidies.txt")
    write_list_to_file(pdx_ploidies, output_dir + "pdx_ploidies.txt")

    print("comparing ploidy between pairs")
    white_samples = read_all_white_file(all_white_file, cohort_names)
    #plot and print to file differences in ploidies and tfxs between pairs
    pairings_comparison(cohort_names, filepath_pairings_files,
                        use_all_pairings, tumor_fractions,
                        ploidies, output_dir, output_filename,
                        white_samples, excl_white_output_filename)
    print("done")


      