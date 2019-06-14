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

	initial_dir = os.getcwd()

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
	os.chdir(output_dir)
	command_list = ['jscb.py', combined_gbk_path]
	subprocess.call(command_list)
	os.chdir(initial_dir)

	# Now we process the .tsv output file of JSCB and determine which genes are within the genomic islands
	# That have been identified, and which contigs they are on, etc.

	genomic_islands = dict() # dictionary keyed by GI name, holds sets of coordinate ranges
	jscb_output_path = os.path.join(output_dir, 'JSCB_output.tsv')
	with open(jscb_output_path) as jscb_output:
		for i,line in enumerate(jscb_output):
			if i > 0:
				line_list = line.rstrip().split('\t')
				gi_name = line_list[0]
				start_coordinate = int(line_list[1])
				end_coordinate = int(line_list[2])
				gi_set = set(range(start_coordinate, end_coordinate + 1))
				genomic_islands[gi_name] = gi_set

	# Now we collect the locus_tags of genes that fall into the genomic islands coordinate ranges
	genomic_island_genes = dict() # keyed by GI name, holds sets of locus_tags
	for gi_name in genomic_islands:
		genomic_island_genes[gi_name] = set()

	for seq_record in SeqIO.parse(combined_gbk_path, 'genbank'):
		for feature in seq_record.features:
			if feature.type == 'gene':
				# Work out if the gene overlaps with any genomic islands
				gene_coordinate_set = set(range(int(feature.location.start), int(feature.location.end) + 1))
				locus_tag = feature.qualifiers['locus_tag'][0] # Here we take just the first value
				for gi_name in genomic_islands:
					if len(gene_coordinate_set.intersection(genomic_islands[gi_name])) > 0:
						genomic_island_genes[gi_name].add(locus_tag)

	# Now we go through the original genbank file, noting the coordinates of genes and which contigs they
	# belong to. Note, we do not assume that genomic islands won't occur at contig boundaries because
	# the JSCB program did not see where these boundaries are.
	genomic_island_coordinates = dict() # Keyed by gi_name, then contig, holds two member list [lowest_coord, highest_coord] derived from component genes
	genomic_island_genes_by_contig = dict() # Keyed by gi_name, then contig, holds lists of locus_tags

	for gi_name in genomic_islands:
		genomic_island_coordinates[gi_name] = dict()
		genomic_island_genes_by_contig[gi_name] = dict()

	for seq_record in SeqIO.parse(genbank_path, 'genbank'):
		contig_name = seq_record.id
		for feature in seq_record.features:
			if feature.type == 'gene':
				locus_tag = feature.qualifiers['locus_tag'][0]
				for gi_name in genomic_island_genes:
					if locus_tag in genomic_island_genes[gi_name]:
						start_coordinate = int(feature.location.start)
						end_coordinate = int(feature.location.end)
						if contig_name in genomic_island_coordinates[gi_name]:
							if start_coordinate < genomic_island_coordinates[gi_name][contig_name][0]:
								genomic_island_coordinates[gi_name][contig_name][0] = start_coordinate
							if end_coordinate > genomic_island_coordinates[gi_name][contig_name][1]:
								genomic_island_coordinates[gi_name][contig_name][1] = end_coordinate
						else:
							genomic_island_coordinates[gi_name][contig_name] = [ start_coordinate, end_coordinate ]

						if contig_name in genomic_island_genes_by_contig[gi_name]:
							genomic_island_genes_by_contig[gi_name][contig_name].append(locus_tag)
						else:
							genomic_island_genes_by_contig[gi_name][contig_name] = [ locus_tag ]

	# Now we make a human-readable output table
	output_table_path = os.path.join(output_dir, 'hgt_genes_summary.tsv')
	output_table = open(output_table_path, 'w')
	output_table.write('GI_ID\tcontig\tstart\tend\tgenes\n')

	for gi_name in genomic_island_coordinates:
		if len(genomic_island_coordinates[gi_name]) > 1:
			# This means the GI spans multiple contigs
			# Because the order of contigs in the original file is likely random,
			# This means in reality there are two GIs
			GI_counter = 0
			for contig_name in genomic_island_coordinates[gi_name]:
				GI_counter += 1
				output_gi_name = gi_name + '_' + str(GI_counter)
				gene_string = ','.join(genomic_island_genes_by_contig[gi_name][contig_name])
				output_string = '\t'.join([ output_gi_name, contig_name, str(genomic_island_coordinates[gi_name][contig_name][0]), str(genomic_island_coordinates[gi_name][contig_name][1]), gene_string ])
				output_table.write(output_string + '\n')
		else:
			for contig_name in genomic_island_coordinates[gi_name]:
				gene_string = ','.join(genomic_island_genes_by_contig[gi_name][contig_name])
				output_string = '\t'.join([ gi_name, contig_name, str(genomic_island_coordinates[gi_name][contig_name][0]), str(genomic_island_coordinates[gi_name][contig_name][1]), gene_string ])
				output_table.write(output_string + '\n')

main()
