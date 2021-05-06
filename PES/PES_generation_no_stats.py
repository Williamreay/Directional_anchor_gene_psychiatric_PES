#!/usr/bin/env python3

"""
Created on Thurs 29 April 2020

Generating pharmagenic enrichment scores (PES)

Author: William Reay - william.reay@uon.edu.au - github: https://github.com/Williamreay
"""

import os
import sys
import argparse
from subprocess import Popen
import pdb

##Add command line flags for automation of PES automation protocol, set argsDict variable to be a dictionary containing the command line arguments

## Will calculate relevant statistics in the target cohort of interest


parser = argparse.ArgumentParser('Generation of PES profiles at a pre-specified P value')
parser.add_argument('--pathway_gwasfile', metavar=['filename'], type = str, help = 'GWAS summary statistics for variants only in genes part of the pathway')
parser.add_argument('--bfile', metavar='{prefix}', type = str, help = 'Target genotype files in plink binary format - .bed, .bim, .fam')
parser.add_argument('--PRsice2_binary', metavar=['filename'],type = str, help = 'Binary executable for PRsice2')
parser.add_argument('--phenotype', metavar='[str]', type = str, help = 'GWAS phenotype name')
parser.add_argument('--pathway', metavar='[str]', type = str, help = 'Pathway_name')
parser.add_argument('--A1', metavar='[str]', type = str, help = 'Name of the A1 allele column [effect allele] in the summary statistics')
parser.add_argument('--A2', metavar='[str]', type = str, help = 'Name of the A2 allele column [non-effect allele] in the summary statistics')
parser.add_argument('--CHR', metavar='[str]', type = str, help = 'Name of chromosome column in the summary statistics')
parser.add_argument('--BP', metavar='[str]', type = str, help = 'Name of SNP position column in the summary statistics')
parser.add_argument('--SNP', metavar='[str]', type = str, help = 'Name of SNP ID column in the summary statistics')
parser.add_argument('--p_threshold', metavar='[str]', type = str, help = 'P-value threshold for SNP inclusion')
parser.add_argument('--statistic', metavar='[str]', type = str, help = 'Statistic for PES calculation, i.e. OR or beta')
parser.add_argument('--P', metavar='[str]', type = str, help = 'Name of P-value column in the summary statistics')
parser.add_argument('--target_type', metavar='[str]', type = str, help = 'Whether target phenotype is binary')
inputFlags = parser.parse_args()

print(inputFlags)

## Calculate PES

Popen(inputFlags.PRsice2_binary + """ --A1 """ + inputFlags.A1 + """ --A2 """ + inputFlags.A2 + """ --chr """ + inputFlags.CHR + """ \
     --bp """ + inputFlags.BP + """ --stat """ + inputFlags.statistic + """ --pvalue """ + inputFlags.P + """ --snp """ + inputFlags.SNP +  """ \
     --base """ + inputFlags.pathway_gwasfile + """ --bar-levels """ + inputFlags.p_threshold + """ --fastscore --model add --no-regress --print-snp --perm 10000 --target \
     """ + inputFlags.bfile + """ --binary-target """ + inputFlags.target_type + """ --out """ + inputFlags.phenotype + """_""" + inputFlags.pathway, shell = True).wait()
