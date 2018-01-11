#!/usr/bin/python

'''
This program takes two vcf files, one containing all SNPs in a population and
the other containing a set of fixed SNPs within a subset of the population. A
dataset is then built for a randomness test of dispersion of these fixed SNPs
within the total SNP set.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--all_SNPs',required=True,type=str,help='A vcf file of all snps observed in the population')
ap.add_argument('--fixed_SNPs',required=True,type=str,help='A vcf file of fixed snps within a subset of the population')

conf = ap.parse_args()

with open(conf.all_SNPs) as f:
    all_snp_lines = f.readlines()
with open(conf.fixed_SNPs) as f:
    fixed_snp_lines = f.readlines()


#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------

class Contig_SNPs(object):
    """A contig and associated SNPs found in a vcf file.

    Attributes:
        name: A string representing the contig name.
        snps: A list of integers representing bp locations of SNPs within the contigs
        fixed_snps: A list of integers representing bp locations of fixed SNPs from a subpopulation within the contigs
    """

    def __init__(self, contig, snps = [], fixed_snps = []):
        """Return a Contig_SNPs whose name is *name*, with no snps or fixed_snps"""
        self.name = contig
        self.snps = snps
        self.fixed_snps = fixed_snps

    def set_snps(self, snp_list):
        """Reset snps to a list of SNP locations"""
        self.snps = snp_list

    def set_fixed_snps(self, snp_list):
        """Reset fixed snps to a list of SNP locations"""
        self.fixed_snps = snp_list

    def add_snp(self, snp_location):
        """Add a SNP to snps"""
        self.snps.append(snp_location)

    def add_fixed_snp(self, snp_location):
        """Add a SNP to fixed_snps"""
        self.fixed_snps.append(snp_location)

    def get_randomness_string(self):
        """Determine if a fixed_snp is present for each snp in this object.
        Return a string that can be used for randomness testing."""
        fixed_set = Set(self.fixed_snps)
        variable_list = self.snps
        randomness_list = []
        for snp in variable_list:
            if snp in fixed_set:
                randomness_list.append('1')
            else:
                randomness_list.append('0')
        return "".join(randomness_list)



#-----------------------------------------------------
# Step 3
# Create an object for each contig containing associated SNPs
#-----------------------------------------------------

object_dict = {}
prev_contig = 'first'
object_list = []
for line in all_snp_lines:
    if line.startswith('#'):
        continue
    line = line.rstrip()
    split_line = line.split()
    contig = split_line[0]
    if prev_contig == 'first':
        contig_snps = []
    elif contig != prev_contig:
        contig_obj = Contig_SNPs(prev_contig)
        contig_obj.set_snps(contig_snps)
        object_dict[prev_contig] = contig_obj
        contig_snps = []
    prev_contig = contig
    snp_pos = split_line[1]
    contig_snps.append(snp_pos)
contig_obj = Contig_SNPs(prev_contig)
contig_obj.set_snps(contig_snps)
object_dict[prev_contig] = contig_obj
contig_snps = []


#-----------------------------------------------------
# Step 4
# Add SNP information into each contig detailing location
# of fixed SNPs
#-----------------------------------------------------

prev_contig = 'first'
object_list = []
for line in fixed_snp_lines:
    if line.startswith('#'):
        continue
    line = line.rstrip()
    split_line = line.split()
    contig = split_line[0]
    if prev_contig == 'first':
        contig_snps = []
    elif contig != prev_contig:
        contig_obj = object_dict[prev_contig]
        # contig_obj = contig_obj_list[0]
        contig_obj.set_fixed_snps(contig_snps)
        object_dict[prev_contig] = contig_obj
        contig_snps = []
    prev_contig = contig
    snp_pos = split_line[1]
    contig_snps.append(snp_pos)
contig_obj = object_dict[prev_contig]
contig_obj.set_fixed_snps(contig_snps)
object_dict[prev_contig] = contig_obj


#-----------------------------------------------------
# Step 5
# go through each SNP in a contig and print a 0 if
# it is variable and a 1 if it is fixed.
#-----------------------------------------------------

keys_list = object_dict.keys()
sorted_list = sorted(keys_list, key = lambda x: int(x.split('_')[1]))
# print sorted_list
for contig_name in sorted_list:
    contig_obj = object_dict[contig_name]
    # print contig_obj.name
    # print contig_obj.snps
    # print contig_obj.fixed_snps
    print "\t".join([contig_obj.name, contig_obj.get_randomness_string()])
