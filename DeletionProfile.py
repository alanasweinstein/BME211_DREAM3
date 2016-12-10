#!/usr/bin/env python2.7

"""
BME 211- DREAM 3 Assignment

Calculate the probabilities of a genetic network using deletion data

Data comes from DataParser as dictionary -> (knockout, gene): expression_val
- or "wt" instead of knockout

Estimate variance of each gene in the data (down a column)

Generate the CDF of the standard Gaussian distribution for each gene a (based on variance found)
-these are our initial reference points for wild type expression level (X(wt)) denoted by Xa(wt)

"""
from DataParser2 import *
import numpy
from collections import defaultdict
import math


def iterate(P, wt_exp, reestimate_wt, first = False):

    gene_exp = defaultdict(lambda: defaultdict(dict))

    #reformat P dictionary for gene expression levels
    for strain, x in P.items():
        knockout, gene = strain
        if knockout == 'wt':
            wt_exp[gene] = float(x)
        else:
            # gene = a, knockout = b
            gene_exp[gene][knockout] = float(x)

    var = 0

    # step 2
    for genePair, expLevel in P.items():
        # print(expLevel)
        b, a = genePair
        var += (expLevel - wt_exp[a]) ** 2

    var /= len(P)
    std = math.sqrt(var)

    new_wt_exp = {}

    #step3
    #calculate new wt expression levels for iteration
    for geneA, exp in wt_exp.items():
        setB = gene_exp[geneA].values()
        new_wt_exp[geneA] = (exp + sum(setB)) / (1.0 + len(setB))


    if first or not reestimate_wt:
        new_wt_exp = wt_exp


    gene_probs = {}  # all probabilities
    set_P = {}  # save to refine model: p_reg <0.95
    potential_regulators = list()  # if p_reg > 0.95


    for (b, a), x in P.items():
        if b not in ['wt', a]:
            if first:
                std = numpy.std(gene_exp[a].values())
            # gene_exp[a][b] = obs exp level of gene a in deletion strain of gene b
            # equation from Yip et al
            p_reg = 2 * norm.cdf(abs(float(x) - wt_exp[a]) / std) - 1
            # probabilities that gene b regulates gene a
            gene_probs[(b, a)] = p_reg

            if p_reg < 0.95:
                set_P[(b, a)] = float(x)
            else:
                potential_regulators.append((b, a))

    return set_P, new_wt_exp, std

#change whether to reestimate here
reestimate_wt = False
if reestimate_wt:
    print('reestimate WT')
else:
    print('Do no reestimate WT')


input_dict, initial_wt = deletionDataParser()

P, WT_exp, std = iterate(input_dict, initial_wt, reestimate_wt, first = True,)

while True:
    newP, newWT_exp, std = iterate(P, WT_exp, reestimate_wt)
    if newP == P:
        break
    else:
        P = newP
        if reestimate_wt:
            WT_exp = newWT_exp

        new_std = std

gene_probs = {}
potential_regulators = list()


#change the threshold here
threshold = 0.95
print('threshold:{}'.format(threshold))

for genePair, exp in input_dict.items():
    b, a = genePair
    #print(b)
    if b not in ['wt', a]:
        # gene_exp[a][b] = obs exp level of gene a in deletion strain of gene b
        # equation from Yip et al

        p_reg = 2 * norm.cdf(abs(float(exp) - WT_exp[a]) / new_std) - 1
        #pri/knt(p_reg)
        # probabilities that gene b regulates gene a
        gene_probs[(b, a)] = Decimal(p_reg)

        #print('here')

        if p_reg > threshold:
            signedEdge = None
            if exp < WT_exp[a]:
                signedEdge = '+'
            else:
                signedEdge = '-'

            potential_regulators.append((b, a, signedEdge))

truePairs = [('G2', 'G1'), ('G2', 'G3'), ('G3', 'G4'), ('G3', 'G5'), ('G3', 'G6'),  ('G3', 'G7'),
             ('G8', 'G7'), ('G8', 'G5'), ('G9', 'G5'), ('G9', 'G4'), ('G10', 'G7')]

realPositive = len(truePairs)

#precision = None #TP/(TP + FP)
#recall = None #TP/TP+FN
#specificity = None # TN/TN+FP

TP = 0.
FP = 0.
#FN = 0 # the ones we missed

for (b, a, sign) in sorted(potential_regulators, key=lambda x: x[0]):
    if (b,a) in truePairs:
        TP += 1.
    else:
        FP+= 1.

    print('{} {} regulates {} with prob {:.5f}  {}'.format(b, sign,a, gene_probs[(b, a)], (b,a) in truePairs))

FN = realPositive - TP
precision = TP/(TP+FP)
recall = TP/(TP+FN)

print('{} potential regulators found'.format(len(potential_regulators)))
print('precision: {}'.format(precision))
print('recall: {}'.format(recall))
