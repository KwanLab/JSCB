#!/usr/bin/env python2

from __future__ import division
import argparse
import os
from Bio import SeqIO
import subprocess

def main():
	parser = argparse.ArgumentParser(description="Run the JS-CB algorithm on draft and complete genomes, to determine which genes are likely to be horizontally transferred. See Nucleic Acids Research 2007, 35, 4629-4639 and Open Biology 2017, 7, 170094.")
	parser.add_argument('-g','--genbank', help='Path to annotated genbank file of your genome.', required=True)
	parser.add_argument('-o','--output_dir', help='Path to directory to store output files.', default='.')

	args = vars(parser.parse_args())

	genbank_path = os.path.abspath(args['genbank'])
	output_dir = os.path.abspath(args['output_dir'])

	# Make output directory if it doesn't yet exist
	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)

	# First we need to make a concatenated version of the genbank file
	combined_record = None
	for seq_record in SeqIO.parse(genbank_path, 'genbank'):
		if combined_record is None:
			combined_record = seq_record
		else:
			combined_record = combined_record + seq_record
	combined_gbk_path = os.path.join(output_dir, 'combined.gbk')
	SeqIO.write(combined_record, combined_gbk_path, 'genbank')

	# Now we can run JSCB
	command_list = ['jscb.py', combined_gbk_path]
	subprocess.call(command_list)


main()
