#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PT Tfx vs Discordance Plotting Ichor
July 22nd, 2020
Anna Hoge
Ha Lab
Fred Hutchinson Cancer Research Center
"""

#imports
import os
import matplotlib.pyplot as plt
import re


#working directory
#os.chdir("/Volumes/fh")
os.chdir("/Volumes/ha_g")


#settings
ICHOR_DISCORDANCE_FILE = "projects/PDX/production/ichor_minimal_subclone_discordance/plots/1mb_discordance/pairing_discordance_values.txt"
THRESHOLD_DISCORDANCE_FILE = "projects/PDX/production/thresholds/plots/1mb_discordance/pairing_discordance_values.txt"
TFX_FILE = "projects/PDX/production/ichor_minimal_subclone_discordance/plots/ploidy_and_tfx/pairings_tfx_ploidy_comparisons.txt"
OUTPUT_DIR = "projects/PDX/production/ichor_minimal_subclone_discordance/plots/pt_tfx_vs_discordance/"


#functions
def get_pt_tfxs(tfx_file):
    """
    Given a file of sample pairings with their tumor fractions and ploidies,
    returns a dictionary pt_tfx where keys are sample names and values are
    tumor fractions.
    """
    pt_tfxs = dict()
    with open(tfx_file, "r") as f:
        f.readline()
        lines = f.readlines()
    for line in lines:
        line_parts = line.strip().split()
        sample = line_parts[1]
        if sample not in pt_tfxs:
            passage = get_passage(sample)
            if passage == "PT":
                tfx = float(line_parts[3])
                pt_tfxs[sample] = tfx
    return pt_tfxs
                

def get_pt_discordance_vs_tfx(discordance_file, pt_tfxs):
    """
    Given a list of pairing discordances and a dictionary pt_tfx as explained
    in get_pt_tfxs, returns a list discordances and a list tfxs, which give
    the discordances of each pairing in discordance_file that includes a PT
    and the tumor fraction of that PT, respectively.
    """
    discordances = []
    tfxs = []
    with open(discordance_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        line_parts = line.strip().split()
        sample = line_parts[1]
        passage = get_passage(sample)
        if passage == "PT":
            discordance = float(line_parts[3])
            tfx = pt_tfxs[sample]
            discordances.append(discordance)
            tfxs.append(tfx)
    return discordances, tfxs
            
                    
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


def plot_tfx_vs_discordance(tfxs, discordances, title):
    """
    Given a list tfxs, a list discordances, and a plot title, plots tfx vs
    discordance.
    """
    plt.scatter(tfxs, discordances)
    plt.title(title)
    plt.xlabel("PT tfx")
    plt.ylabel("Discordance values")
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)


def write_discordances_and_tfxs_to_file(discordances, tfxs, filename):
    with open(filename, "w") as f:
        f.write("discordance\ttfx\n")
        for i in range(len(discordances)):
            f.write(str(discordances[i]) + "\t" + str(tfxs[i]) + "\n")
            
            

#run
if __name__ == "__main__":
    ichor_discordance_file = ICHOR_DISCORDANCE_FILE
    threshold_discordance_file = THRESHOLD_DISCORDANCE_FILE
    tfx_file = TFX_FILE
    output_dir = OUTPUT_DIR

    pt_tfxs = get_pt_tfxs(tfx_file)
    ichor_discordances, ichor_tfxs = get_pt_discordance_vs_tfx(ichor_discordance_file, pt_tfxs)
    threshold_discordances, threshold_tfxs = get_pt_discordance_vs_tfx(threshold_discordance_file, pt_tfxs)
    plot_tfx_vs_discordance(ichor_tfxs, ichor_discordances, "ichorCNA")
    plt.savefig(output_dir + "pt_tfx_vs_discordance_ichor.png")
    plt.show()
    plot_tfx_vs_discordance(threshold_tfxs, threshold_discordances, "thresholds")
    plt.savefig(output_dir + "pt_tfx_vs_discordance_thresholds.png")
    plt.show()
    
    write_discordances_and_tfxs_to_file(ichor_discordances, ichor_tfxs, output_dir + "pt_tfx_vs_discordance_ichor.txt")
    write_discordances_and_tfxs_to_file(threshold_discordances, threshold_tfxs, output_dir + "pt_tfx_vs_discordance_thresholds.txt")

                
