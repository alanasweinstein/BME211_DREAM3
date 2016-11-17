#!/usr/bin/env python2.7

'''
BME 211- DREAM 3 Assignment

Calculate the probabilities of a genetic network using deletion data


Data comes from DataParser as dictionary -> (knockout, gene): expression_val
	- or "wt" instead of knockout

Estimate variance of each gene in the data (down a column)

Generate the CDF of the standard Gaussian distribution for each gene a (based on variance found)
	-these are our initial reference points for wild type expression level (X(wt)) denoted by Xa(wt)

Genes 4 & 9 have a regulatory relationship: starting p=0.017
'''

import DataParser
import sys, numpy
from scipy.stats import norm
from collections import defaultdict

def del_profile(infile=sys.stdin):
	'''
	Calculate probabilities gene b regulates gene a
		- WT mean to start (will change in iterations)
		- loop through input_dict for (b,a):x
		- use gene_exp to calculate stdevs for gene a
		- compute probability of b -> and store in dict
	'''
	#input_dict = {(knockout, gene):expression_level, ('wt', gene):exp_level, ...}
	#(knockout, gene) = (b, a)
	input_dict, gene_list = DataParser.deletionDataParser(infile)
	
	wt_exp = {}
	gene_exp = defaultdict(lambda : defaultdict(dict))
	
	#reorder input_dict for calculating stdevs
	for strain, x in input_dict.items():
		knockout, gene = strain
		if knockout == 'wt':
			wt_exp[gene] = float(x)
		else:
			#gene = a, knockout = b
			gene_exp[gene][knockout] = float(x)
	
	wt_mean = numpy.mean(wt_exp.values())
	
	gene_probs = {} #all probabilities
	set_P = {} #save to refine model: p_reg > 0.05
	potential_regulators = {} #if p_reg < 0.05
	
	for (b, a), x in input_dict.items():
		if b != 'wt':
			a_std = numpy.std(gene_exp[a].values())
			#gene_exp[a][b] = obs exp level of gene a in deletion strain of gene b
			#equation from Yip et al
			p_reg = 2 * norm.cdf(abs(gene_exp[a][b] - wt_mean)/a_std) - 1
			#probabilities that gene b regulates gene a
			gene_probs[(b,a)] = p_reg
			
			if prob > 0.05:
				set_P[(b,a)] = x
			else:
				potential_regulators[b] = a
	
	#same format as input dict, but probabilities instead of expression values
	return gene_probs, set_P