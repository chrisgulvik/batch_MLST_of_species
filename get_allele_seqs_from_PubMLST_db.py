#!/usr/bin/env python


import argparse
import collections
import csv
import os
import re
import sys
import numpy as np
from Bio import SeqIO
from itertools import izip


def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-p', '--profiles', required=True, 
		help='tab-delimited file of variant profiles from PubMLST that define each ST')
	parser.add_argument('-s', '--ST', type=int, help='a ST to determine')
	parser.add_argument('-e', '--ext', default='tfa', help='file extension of FastA allele sequences file')
	parser.add_argument('-o', '--outdir', default=os.getcwd(), help='output directory')
	return parser.parse_args()

def get_locus_names(profiles):
	with open(profiles, 'r') as f:
		head = f.readline().rstrip()
		locus_names = [s for s in head.split('\t')]
		locus_names.pop(0)  #discard first item 'ST'
		if 'clonal_complex' in locus_names:
			locus_names.remove('clonal_complex')
		if 'mlst_clade' in locus_names:
			locus_names.remove('mlst_clade')
	print 'locus_names %s' % locus_names
	return locus_names

def locus_names_to_dict(locus_names, ext):
	'''
	Given a list of individual locus names,
	generate a dictionary of allele names and allele sequences.
	'''
	loci_dict = {}
	for locus in locus_names:
		locus_dict = locus + '_dict'
		with open(locus + '.' + ext, 'r') as f:
			locus_dict = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
		loci_dict.update(locus_dict)
	return loci_dict

def make_profiles_db(profiles):
	'''
	Creates a dictionary of profiles, 
	where ST numbers are keys and variant profiles are values
	'''
	ST_numbers  = np.genfromtxt(profiles, dtype=int, delimiter='\t', skip_header=1, usecols=0)
	ST_profiles = np.genfromtxt(profiles, dtype=int, delimiter='\t', skip_header=1, usecols=range(1,8))
	return dict(izip(ST_numbers,ST_profiles.tolist()))

def make_alleles_multifasta(locus_names_varnums, loci_dict, ST, outdir):
	'''
	Creates a multi-FastA file containing all allele sequences.
	'''
	out_series = []
	for allele in locus_names_varnums:
		out_series.append(loci_dict[allele])
	out_multifasta = open(os.path.join(outdir, 'ST-' + ST + '.fas'), 'w')
	SeqIO.write(out_series, out_multifasta, 'fasta')
	out_multifasta.close()

def make_concat_seq(ST, outdir):
	'''
	Given a ST-number multi-FastA file containing all allele sequences,
	produce a concatenated FastA file with the defline being the ST-number.
	'''
	individual_seqs = []
	with open(os.path.join(outdir, 'ST-' + ST + '.fas')) as out_multifasta:
		for line in out_multifasta:
			individual_seqs += re.findall('^[ATCG]+$', line)
		concat_seq = ''.join(individual_seqs)
	out_concatfasta = open(os.path.join(outdir, 'concat_ST-' + ST + '.fas'), 'w')
	out_concatfasta.write('>ST-' + ST + '\n' + concat_seq + '\n')  #FastA outfmt
	out_concatfasta.close()

def main():
	args = parseArgs()
	profiles = args.profiles
	ext = args.ext

	outdir = args.outdir
	if not os.path.exists(outdir): 
		os.mkdir(outdir)

	MLST_profiles = make_profiles_db(profiles)
	locus_names = get_locus_names(profiles)

	if args.ST is not None:
		ST = int(args.ST)
		num_STs_todo = 1
	else:
		print 'Found %s STs, and will produce FastA files for each...' % str(len(MLST_profiles))
		ST = len(MLST_profiles)
		num_STs_todo = len(MLST_profiles)

	while num_STs_todo > 0:
		profile = MLST_profiles[ST]
		locus_names_varnums = ['%s_%s' % t for t in zip(locus_names, profile)]
		# locus_names_varnums example: #['adk_7', 'atpA_1', 'dxr_1', 'glyA_10', 'recA_7', 'sodA_6', 'tpi_3']
		loci_dict = locus_names_to_dict(locus_names, ext)
		make_alleles_multifasta(locus_names_varnums, loci_dict, str(ST), outdir)
		make_concat_seq(str(ST), outdir)
		ST -= 1
		num_STs_todo -= 1

if __name__ == '__main__':
	main()
