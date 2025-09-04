#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module computes a multiple sequence alignment (MSA) based on allele calling results.

Code documentation
------------------
"""


import os
import sys

import pandas as pd

try:
	from GetAlleles import get_alleles
	from utils import (
		constants as ct,
		mafft_wrapper as mw,
		file_operations as fo,
		fasta_operations as fao,
		iterables_manipulation as im,
		multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.GetAlleles import get_alleles
	from CHEWBBACA.utils import (
		constants as ct,
		mafft_wrapper as mw,
		file_operations as fo,
		fasta_operations as fao,
		iterables_manipulation as im,
		multiprocessing_operations as mo)


def concatenate_loci_alignments(sample, loci, sample_profile, fasta_index, output_directory):
	"""Concatenate the aligned sequences for a sample.

	Parameters
	----------
	sample : str
		Sample identifier.
	loci : list
		Loci identifiers.
	fasta_index : Bio.File._IndexedSeqFileDict
		Indexed FASTA file to get sequences from.
	output_directory : str
		Path to the output directory.

	Returns
	-------
	alignment_outfile : str
		Path to the FASTA file with the concatenated aligned sequences.
	"""
	alignment = ''
	for locus in loci:
		# Sequence headers include sample, locus, and allele IDs joined by '_'
		# Get allele ID
		allele_id = sample_profile[locus].tolist()[0]
		# Remove '*' from allele ID
		allele_id = allele_id.replace('*', '')
		# Get aligned sequence from index
		try:
			seqid = f'{sample}_{locus}_{allele_id}'
			alignment += str(fasta_index[seqid].seq)
		except Exception as e:
			seqid = f'{sample}_{locus}_0'
			alignment += str(fasta_index[seqid].seq)
	# Save alignment for sample
	alignment_outfile = fo.join_paths(output_directory,
									  [f'{sample}_cgMLST_alignment.fasta'])
	alignment_record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
											[sample, alignment])
	fo.write_lines([alignment_record], alignment_outfile)

	return alignment_outfile


def add_gaps(input_file, locus_id, output_file, gap_char, sample_ids):
	"""Add gap sequences to a MSA for samples that do not have an allele.

	Parameters
	----------
	input_file : str
		Path to the FASTA file with the aligned sequences.
	locus_id : str
		Locus identifier.
	gap_char : str
		Character to use to fill gaps.
	sample_ids : list
		Sample identifiers.

	Returns
	-------
	gapped_fasta : str
		Path to the FASTA file containing the updated MSA.
	"""
	records = fao.import_sequences(input_file)
	# Get list of samples where locus was identified
	dataset_sids = {k.split(f'_{locus_id}')[0]: k for k in records}
	# Get length of alignment to create gapped sequence
	# Just get length of the first record
	msa_len = len(records[list(records.keys())[0]])
	gapped_seq = gap_char * msa_len
	gapped_records = {}
	for sid in sample_ids:
		# Sample contains locus
		if sid in dataset_sids:
			gapped_records[dataset_sids[sid]] = records[dataset_sids[sid]]
		# Sample does not contain locus
		# Add gapped sequence
		else:
			gapped_seq_sid = f'{sid}_{locus_id}_0'
			gapped_records[gapped_seq_sid] = gapped_seq

	# Save Fasta file with gapped sequences
	outrecords = [f'>{k}\n{v}' for k, v in gapped_records.items()]
	fo.write_lines(outrecords, output_file)

	return output_file


def convert_msa_to_dna(input_file, dna_file, output_file, gap_char):
	"""
	"""
	protein_records = fao.import_sequences(input_file)
	dna_sequences = fao.import_sequences(dna_file)
	dna_records = {}
	for seqid, sequence in protein_records.items():
		# Check if it matches any record in the schema
		if seqid in dna_sequences:
			# Get allele sequence
			allele = dna_sequences[seqid]
			# Iterate over gapped protein sequence to create gapped DNA
			dna_index = 0
			gapped_dna = ''
			for i, char in enumerate(sequence):
				# Add codon if it is not a gap
				if char != gap_char:
					gapped_dna += allele[dna_index:dna_index+3]
					dna_index += 3
				# Add '---' if it is a gap
				elif char == gap_char:
					gapped_dna += gap_char * 3
			dna_records[seqid] = gapped_dna
		else:
			dna_records[seqid] = sequence * 3

	# Save DNA MSA
	dna_msa_recs = [f'>{k}\n{v}' for k, v in dna_records.items()]
	fo.write_lines(dna_msa_recs, output_file)

	return output_file


def determine_variable_positions(msa_files, output_files, ignore_missing=True):
	"""
	"""
	# Determine variable positions for protein MSAs
	variable_loci = []
	non_variable = []
	for fi, file in enumerate(msa_files):
		# Read MSA
		msa_records = fao.import_sequences(file)
		msa_sequences = list(msa_records.values())
		# Zip to pair chars in same positions
		zipped = list(zip(*msa_sequences))
		variable = []
		for i, pos in enumerate(zipped):
			# Get unique characters in position
			distinct = set(pos)
			# Optionally ignore gaps and missing data (N)
			if ignore_missing:
				if '-' in distinct:
					distinct.remove('-')
				if 'N' in distinct:
					distinct.remove('N')
			# If there are more than 1 unique character, it is a variable position
			if len(distinct) > 1:
				variable.append([i, pos])
		if len(variable) > 0:
			variable_pos = list(zip(*[i[1] for i in variable]))
			variable_pos = [''.join(i) for i in variable_pos]
			variable_records = [[allele_id, variable_pos[i]] for i, allele_id in enumerate(msa_records.keys())]
			variable_records = fao.fasta_lines(ct.FASTA_RECORD_TEMPLATE, variable_records)
			# Save variable positions
			fo.write_lines(variable_records, output_files[fi])
			variable_loci.append(output_files[fi])
		else:
			non_variable.append(output_files[fi])

	return variable_loci, non_variable


def create_full_msa(input_files, output_directory, loci_ids, sample_ids, results_alleles):
	"""
	"""
	concat = fo.join_paths(output_directory, ['concat'])
	# Concatenate all alignment files and index with BioPython
	fo.concatenate_files(input_files, concat)
	# Index file
	indexed_protein_concat = fao.index_fasta(concat)
	sample_MSA_outfiles = []
	for sid in sample_ids:
		# Get sample profile
		sample_profile = pd.read_csv(results_alleles,
										skiprows=range(1,sample_ids.index(sid)+1), # Only get header and sample profile
										nrows=1,
										delimiter='\t',
										header=0,
										dtype=str)
		alignment_file = concatenate_loci_alignments(sid,
														loci_ids,
														sample_profile,
														indexed_protein_concat,
														output_directory)
		sample_MSA_outfiles.append(alignment_file)

	return sample_MSA_outfiles


# Test
input_file = '/home/rmamede/test_chewie/features/ComputeMSA/spyogenes_data/results_alleles.tsv'
schema_directory = '/home/rmamede/test_chewie/features/ComputeMSA/spyogenes_data/spyogenes_wgMLST'
output_directory = '/home/rmamede/test_chewie/features/ComputeMSA/spyogenes_data/results_msa'
dna_msa = True
output_variable = True
gap_char = '-'
translation_table = 11
cpu_cores = 12
keep_locus_msa = False
only_locus_msa = False
def main(input_file, schema_directory, output_directory, dna_msa, output_variable, gap_char, translation_table, cpu_cores, keep_locus_msa, only_locus_msa):
	# Create output directory
	fo.create_directory(output_directory)
	# Get sample IDs
	sample_ids = fo.extract_column(input_file, delimiter='\t', column_index=0)

	# Create FASTA files with the alleles identified in the dataset
	# Call GetAlleles module
	# Removes the '*' in allele IDs
	_, dna_files, protein_files = get_alleles.main(input_file, schema_directory, None, output_directory, cpu_cores, False, True, translation_table)

	# Get FASTA files for loci identified in at least one sample
	if len(dna_files) == 0:
		sys.exit(ct.COMPUTEMSA_NO_ALLELES)

	# Run MAFFT to compute MSA
	print('\nRunning MAFFT to compute the MSA for each locus...')
	mafft_outdir = fo.join_paths(output_directory, ['mafft_results'])
	fo.create_directory(mafft_outdir)
	mafft_outfiles = [os.path.basename(file) for file in protein_files]
	mafft_outfiles = [fo.join_paths(mafft_outdir, [file]) for file in mafft_outfiles]
	mafft_inputs = [[file, mafft_outfiles[i]] for i, file in enumerate(protein_files)]
	common_args = []
	# Add common arguments to all sublists
	inputs = im.multiprocessing_inputs(mafft_inputs, common_args, mw.call_mafft)
	mafft_results = mo.map_async_parallelizer(inputs,
											  mo.function_helper,
											  cpu_cores,
											  show_progress=True)

	# Identify cases where MAFFT failed
	mafft_failed = [r[0] for r in mafft_results if r[1] is False]
	if len(mafft_failed) > 0:
		print(f'\nCould not determine MSA for {len(mafft_failed)} loci.')

	# Get files created by MAFFT
	mafft_success = [r[0] for r in mafft_results if r[1] is True]

	# Create folder to store gapped MSAs
	gapped_outdir = fo.join_paths(output_directory, ['gapped_MSAs'])
	fo.create_directory(gapped_outdir)

	# Add gap sequences when sample did not have an allele
	gapped_inputs = []
	for file in mafft_success:
		locus_id = fo.file_basename(file, False).split('_protein')[0]
		outfile_protein_gapped = fo.join_paths(gapped_outdir, [f'{locus_id}_protein_gapped.fasta'])
		gapped_inputs.append([file, locus_id, outfile_protein_gapped])

	common_args = [gap_char, sample_ids]
	gapped_inputs = im.multiprocessing_inputs(gapped_inputs, common_args, add_gaps)
	gapped_results = mo.map_async_parallelizer(gapped_inputs,
											   mo.function_helper,
											   cpu_cores,
											   show_progress=True)

	# Delete original file with ungapped MSA
	fo.delete_directory(mafft_outdir)

	# Convert protein MSAs to DNA MSAs
	if dna_msa:
		dna_inputs = []
		for i, file in enumerate(gapped_results):
			locus_id = fo.file_basename(file, False).split('_protein_gapped')[0]
			outfile_dna_gapped = fo.join_paths(gapped_outdir, [f'{locus_id}_dna_gapped.fasta'])
			dna_inputs.append([file, dna_files[i], outfile_dna_gapped])

		common_args = [gap_char]
		dna_inputs = im.multiprocessing_inputs(dna_inputs, common_args, convert_msa_to_dna)
		dna_results = mo.map_async_parallelizer(dna_inputs,
										        mo.function_helper,
												cpu_cores,
												show_progress=True)

	# Delete directories containing the FASTA files with DNA and protein sequences
	fo.delete_directory(os.path.dirname(dna_files[0]))
	fo.delete_directory(os.path.dirname(protein_files[0]))

	# Identify variable positions to get SNP MSA
	if output_variable:
		# Create folder to store variable position MSAs
		variable_outfolder = fo.join_paths(output_directory, ['variable_positions_MSAs'])
		fo.create_directory(variable_outfolder)
		# Define paths to output files
		variable_protein_outfiles = [fo.file_basename(file).replace('_gapped', '_variable') for file in gapped_results]
		variable_protein_outfiles = [fo.join_paths(variable_outfolder, [file]) for file in variable_protein_outfiles]
		variable_protein, non_variable_protein = determine_variable_positions(gapped_results, variable_protein_outfiles)

		# Determine variable positions for DNA MSAs
		if dna_msa:
			variable_dna_outfiles = [fo.file_basename(file).replace('_gapped', '_variable') for file in dna_results]
			variable_dna_outfiles = [fo.join_paths(variable_outfolder, [file]) for file in variable_dna_outfiles]
			variable_dna, non_variable_dna = determine_variable_positions(dna_results, variable_dna_outfiles)

	# User only wants the locus MSAs
	# Do not compute full MSAs
	if only_locus_msa:
		return

	# Create folder to store sample MSAs
	sample_protein_msas_folder = fo.join_paths(output_directory, ['sample_protein_MSAs'])
	fo.create_directory(sample_protein_msas_folder)

	# Create the full protein MSA
	print('\nCreating file with the full protein MSA...')
	# Get loci IDs for which a gapped protein MSA was created
	gapped_results_loci_ids = [fo.file_basename(file).split('_protein_gapped')[0] for file in gapped_results]
	sample_protein_MSA_outfiles = create_full_msa(gapped_results, sample_protein_msas_folder, gapped_results_loci_ids, sample_ids, input_file)
	# Concatenate sample protein alignments
	full_protein_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_PROTEIN_MSA])
	fo.concatenate_files(sample_protein_MSA_outfiles, full_protein_alignment)

	# Get length of the full alignment
	protein_msa_length = len(fo.read_lines(full_protein_alignment, strip=True, num_lines=2)[1])
	print(f'Protein MSA length: {protein_msa_length}')

	if dna_msa:
		sample_dna_msas_folder = fo.join_paths(output_directory, ['sample_dna_MSAs'])
		fo.create_directory(sample_dna_msas_folder)
		# Create the full DNA MSA
		print('Creating file with the full DNA MSA...')
		sample_dna_MSA_outfiles = create_full_msa(dna_results, sample_dna_msas_folder, gapped_results_loci_ids, sample_ids, input_file)
		# Concatenate sample dna alignments
		full_dna_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_DNA_CONCAT])
		fo.concatenate_files(sample_dna_MSA_outfiles, full_dna_alignment)

		# Get length of the full alignment
		dna_msa_length = len(fo.read_lines(full_dna_alignment, strip=True, num_lines=2)[1])
		print(f'Protein MSA length: {dna_msa_length}')

	# Delete folders with sample MSAs
	fo.delete_directory(sample_protein_MSA_outfiles)
	fo.delete_directory(sample_dna_msas_folder)

	# Create full MSAs with only variable positions if requested
	if output_variable:
		sample_protein_variable_msas_folder = fo.join_paths(output_directory, ['sample_protein_variable_MSAs'])
		fo.create_directory(sample_protein_variable_msas_folder)

		# Create the full protein MSA for the variable positions
		print('\nCreating file with the full protein MSA for the variable positions...')
		# Get loci IDs for which a gapped protein MSA was created and had variable positions
		variable_results_loci_ids = [fo.file_basename(file).split('_protein_variable')[0] for file in variable_protein]
		sample_protein_variable_MSA_outfiles = create_full_msa(variable_protein, sample_protein_variable_msas_folder, variable_results_loci_ids, sample_ids, input_file)
		# Concatenate sample protein alignments
		full_variable_protein_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_PROTEIN_MSA_VARIABLE])
		fo.concatenate_files(sample_protein_variable_MSA_outfiles, full_variable_protein_alignment)

		# Get length of the full alignment
		protein_variable_msa_length = len(fo.read_lines(full_variable_protein_alignment, strip=True, num_lines=2)[1])
		print(f'Protein variable MSA length: {protein_variable_msa_length}')

		if dna_msa:
			sample_dna_variable_msas_folder = fo.join_paths(output_directory, ['sample_dna_variable_MSAs'])
			fo.create_directory(sample_dna_variable_msas_folder)
			# Create the full DNA MSA
			print('Creating file with the full DNA MSA for the variable positions...')
			# Get loci IDs for which a gapped DNA MSA was created and had variable positions
			variable_results_loci_ids = [fo.file_basename(file).split('_dna_variable')[0] for file in variable_dna]
			sample_dna_variable_MSA_outfiles = create_full_msa(variable_dna, sample_dna_variable_msas_folder, variable_results_loci_ids, sample_ids, input_file)
			# Concatenate sample dna alignments
			full_variable_dna_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_DNA_MSA_VARIABLE])
			fo.concatenate_files(sample_dna_variable_MSA_outfiles, full_variable_dna_alignment)

			# Get length of the full alignment
			dna_variable_msa_length = len(fo.read_lines(full_variable_dna_alignment, strip=True, num_lines=2)[1])
			print(f'Protein MSA length: {dna_variable_msa_length}')
