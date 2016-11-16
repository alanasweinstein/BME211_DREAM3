#! usr/bin/env python3

import sys

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
    genes = list() #genes extracted from file

    for i, line in enumerate(f):
        if i == 0:
            genes = line.strip().replace('"', '').split('\t')
            genes.pop(0)
        elif line.rstrip():
            #print(line)
            geneLine = line.strip().split('\t')
            thisKnockout = geneLine.pop(0).replace('"', '')
            print(geneLine)
            for i, gene in enumerate(genes):
                if thisKnockout == 'wt':
                    inputData[(thisKnockout, gene)] = geneLine[i]
                else:
                    inputData[(thisKnockout[:-5], gene)] = geneLine[i]

    return inputData, genes
