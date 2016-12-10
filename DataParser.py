#! usr/bin/env python2.7

from __future__ import  print_function

import sys, numpy
from collections import defaultdict
from decimal import *

getcontext().prec = 28

def deletionDataParser(file = ''):
    """Parser for handling deletion data. File to be extracted can be specified when function is called,
    maybe we want to specify something in the command line, or it can be read from sys.stdin

    Possible additions: determine whether heterozygous or homozygous and return that as well as the dictionary"""

    if file is '':
        f = sys.stdin
    else:
        f = open(file)

    inputData = dict() #key = (knockout, gene), value = gene expression level | knockout
                    #or wt, gene for wild type strain
    initial_wt = dict () #dictionary of wildtype expression levels
    initial_gene_exp = defaultdict(lambda: defaultdict(dict)) #dictionary of gene expression levels not including wt

    genes = list() #genes extracted from file

    for i, line in enumerate(f):
        if i == 0:
            genes = line.strip().replace('"', '').split('\t')
            genes.pop(0)
        elif line.rstrip():
            #print(line)
            geneLine = line.strip().split('\t')
            thisKnockout = geneLine.pop(0).replace('"', '')
            #print(geneLine)
            for i, gene in enumerate(genes):
                if thisKnockout == 'wt':
                    inputData[(thisKnockout, gene)] = float(geneLine[i])
                else:
                    inputData[(thisKnockout[:-5], gene)] = float(geneLine[i])

    for strain, x in inputData.items():
        knockout, gene = strain
        if knockout == 'wt':
            initial_wt[gene] = float(x)
        else:
            # gene = a, knockout = b
            initial_gene_exp[gene][knockout] = float(x)

    return inputData, initial_wt


