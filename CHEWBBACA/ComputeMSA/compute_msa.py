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
from Bio import SeqIO

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
									  [f'{sample}.fasta'])
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


def determine_variable_positions(msa_file, output_file, gap_char, ambiguous_chars, gaps, ambiguous):
	"""
	"""
	# Read MSA
	msa_records = fao.import_sequences(msa_file)
	msa_sequences = list(msa_records.values())
	# Zip to pair chars in same positions
	zipped = list(zip(*msa_sequences))
	variable = []
	for i, pos in enumerate(zipped):
		# Get unique characters in position
		distinct = set(pos)

		# Ignore positions if any sequence has a gap or an ambiguous char
		if gaps == 'ignore' and gap_char in distinct:
			distinct.remove('-')
		if ambiguous == 'ignore' and any([c in distinct for c in ambiguous_chars]):
			for c in ambiguous_chars:
				if c in distinct:
					distinct.remove(c)

		# Exclude positions if any sequence has a gap or an ambiguous char
		if gaps == 'exclude' and gap_char in distinct:
			continue
		if ambiguous == 'exclude' and any([c in distinct for c in ambiguous_chars]):
			continue

		# If there are more than 1 unique character, it is a variable position
		if len(distinct) > 1:
			variable.append([i, pos])

	# Save variable positions
	# Create FASTA records with only variable positions
	# If no variable positions, do not create file
	if len(variable) > 0:
		variable_pos = list(zip(*[i[1] for i in variable]))
		variable_pos = [''.join(i) for i in variable_pos]
		variable_records = [[allele_id, variable_pos[i]] for i, allele_id in enumerate(msa_records.keys())]
		variable_records = fao.fasta_lines(ct.FASTA_RECORD_TEMPLATE, variable_records)
		fo.write_lines(variable_records, output_file)
		return [True, output_file]
	else:
		return [False, output_file]

#sample_id, sample_index, fasta_index, output_directory, loci_ids, profiles = inputs[0][:-1]
def create_full_msa(sample_id, sample_index, fasta_file, output_directory, loci_ids, profiles):
	"""
	"""
	# This tries to create an index for the FASTA file
	# but ends up loading an existing index created with SeqIO.index_db
	fasta_index = fao.index_fasta(fasta_file, True)
	# Get sample profile
	sample_profile = pd.read_csv(profiles,
								 skiprows=sample_index, # Only get header and sample profile
								 nrows=1,
								 delimiter='\t',
								 header=0,
								 dtype=str)

	# Remove all 'INF-' prefixes and '*' from identifiers
	# Replace special classifications by '0'
	masked_profile = sample_profile.apply(im.replace_chars)
	alignment_file = concatenate_loci_alignments(sample_id,
												 loci_ids,
												 masked_profile,
												 fasta_index,
												 output_directory)

	return alignment_file


# Test
# input_file = '/home/rmamede/test_chewie/features/ComputeMSA/spyogenes_data/results_alleles_shorter.tsv'
# schema_directory = '/home/rmamede/test_chewie/features/ComputeMSA/spyogenes_data/spyogenes_wgMLST'
# output_directory = '/home/rmamede/test_chewie/features/ComputeMSA/spyogenes_data/results'
# dna_msa = True
# output_variable = True
# translation_table = 11
# cpu_cores = 6
# only_loci_msas = False
# gaps = 'exclude'
# ambiguous = 'exclude'
# custom_mafft_params = None
# protein_input = False
# no_cleanup = False
def main(input_file, output_directory, schema_directory, dna_msa, output_variable,
		 translation_table, cpu_cores, only_loci_msas, gaps, ambiguous,
		 custom_mafft_params, protein_input, no_cleanup):
	# Create output directory
	fo.create_directory(output_directory)

	if os.path.isfile(input_file):
		# Get sample IDs
		sample_ids = fo.extract_column(input_file, delimiter='\t', column_index=0)
		# Create FASTA files with the alleles identified in the dataset
		# Call GetAlleles module
		# Removes the '*' in allele IDs
		_, dna_files, protein_files = get_alleles.main(input_file, schema_directory, None, output_directory, cpu_cores, False, True, translation_table)
	elif os.path.isdir(input_file):
		# Input is a directory with FASTA files
		# Copy input files to output directory
		files_to_copy = fo.listdir_fullpath(input_file, ct.FASTA_EXTENSIONS)
		copied_dir = fo.join_paths(output_directory, ['input_files'])
		fo.create_directory(copied_dir)
		copied_files = []
		for file in files_to_copy:
			fo.copy_file(file, copied_dir)
			copied_files.append(fo.join_paths(copied_dir, [os.path.basename(file)]))

		# Define list of files with DNA and protein sequences
		# If input files contain DNA sequences, translate to protein
		if protein_input is False:
			dna_files = copied_files
			translated_outdir = fo.join_paths(output_directory, ['translated'])
			fo.create_directory(translated_outdir)
			protein_files = [fao.translate_fasta(file, translated_outdir, translation_table) for file in copied_files]
			protein_files = [result[1] for result in protein_files]
		# If input files contain protein sequences
		else:
			dna_files = []
			protein_files = copied_files

	# Exit if no sequences were found
	if len(dna_files) == 0 and len(protein_files) == 0:
		sys.exit(ct.COMPUTEMSA_NO_ALLELES)

	# Run MAFFT to compute protein MSAs
	print('\nRunning MAFFT to compute the MSA for each locus...')
	mafft_outdir = fo.join_paths(output_directory, ['loci_MSAs'])
	mafft_protein_outdir = fo.join_paths(mafft_outdir, ['protein'])
	fo.create_directory(mafft_protein_outdir)
	mafft_outfiles = [os.path.basename(file) for file in protein_files]
	mafft_outfiles = [fo.join_paths(mafft_protein_outdir, [file]) for file in mafft_outfiles]
	mafft_inputs = [[file, mafft_outfiles[i]] for i, file in enumerate(protein_files)]
	common_args = custom_mafft_params.split() if custom_mafft_params is not None else [ct.MAFFT_DEFAULT_PARAMETERS[:-1]]
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
	if len(mafft_success) == 0:
		sys.exit(ct.NO_MSAS_CREATED)

	# Add gap sequences when samples did not have an allele identified
	# Only do this if input is a TSV file with allelic profiles
	# If it is a folder with FASTA files there is no information about samples missing loci
	if not os.path.isdir(input_file):
		# This will create new files with MSAs containing gap sequences (samples missing loci will show as gap sequences)
		gapped_inputs = []
		for file in mafft_success:
			locus_id = fo.file_basename(file, False).split('_protein')[0]
			gapped_inputs.append([file, locus_id, file])
		common_args = [ct.GAP_CHAR, sample_ids]
		gapped_inputs = im.multiprocessing_inputs(gapped_inputs, common_args, add_gaps)
		gapped_results = mo.map_async_parallelizer(gapped_inputs,
												mo.function_helper,
												cpu_cores,
												show_progress=True)
	# If the input was a directory with FASTA files, there is no need to add gap sequences
	else:
		gapped_results = mafft_success

	# Convert protein MSAs to DNA MSAs
	if dna_msa:
		dna_inputs = []
		if protein_input is False:
			aligned_dna_outdir = fo.join_paths(mafft_outdir, ['dna'])
			fo.create_directory(aligned_dna_outdir)
			for i, file in enumerate(gapped_results):
				locus_id = fo.file_basename(file, False).split('_protein')[0]
				outfile_dna = fo.join_paths(aligned_dna_outdir, [f'{locus_id}.fasta'])
				dna_inputs.append([file, dna_files[i], outfile_dna])

			# Convert to DNA MSA
			common_args = [ct.GAP_CHAR]
			dna_inputs = im.multiprocessing_inputs(dna_inputs, common_args, convert_msa_to_dna)
			dna_results = mo.map_async_parallelizer(dna_inputs,
													mo.function_helper,
													cpu_cores,
													show_progress=True)
		else:
			print('\nDNA MSAs were requested, but input files contain protein sequences.')

	# Delete directories containing the FASTA files with DNA and protein sequences
	if not no_cleanup:
		if protein_input is False:
			fo.delete_directory(os.path.dirname(dna_files[0]))
		fo.delete_directory(os.path.dirname(protein_files[0]))

	# Identify variable positions to get SNP MSA
	if output_variable:
		# Create folder to store variable position MSAs
		variable_protein_outfolder = fo.join_paths(mafft_outdir, ['variable_protein'])
		fo.create_directory(variable_protein_outfolder)
		# Define paths to output files
		variable_protein_outfiles = [fo.join_paths(variable_protein_outfolder, [fo.file_basename(file)]) for file in gapped_results]
		common_args = [ct.GAP_CHAR, ct.PROTEIN_AMBIGUOUS_CHARS, gaps, ambiguous]
		variable_inputs = im.multiprocessing_inputs([list(i) for i in zip(gapped_results, variable_protein_outfiles)], common_args, determine_variable_positions)
		variable_results = mo.map_async_parallelizer(variable_inputs,
													 mo.function_helper,
													 cpu_cores,
													 show_progress=True)

		# Collect variable and non-variable protein MSA files
		variable_protein = [r[1] for r in variable_results if r[0] is True]
		non_variable_protein = [r[1] for r in variable_results if r[0] is False]

		# Determine variable positions for DNA MSAs
		if dna_msa:
			if protein_input is False:
				# Create folder to store variable position MSAs
				variable_dna_outfolder = fo.join_paths(mafft_outdir, ['variable_dna'])
				fo.create_directory(variable_dna_outfolder)
				variable_dna_outfiles = [fo.join_paths(variable_dna_outfolder, [fo.file_basename(file)]) for file in dna_results]
				common_args = [ct.GAP_CHAR, ct.PROTEIN_AMBIGUOUS_CHARS, gaps, ambiguous]
				variable_inputs = im.multiprocessing_inputs([list(i) for i in zip(dna_results, variable_dna_outfiles)], common_args, determine_variable_positions)
				variable_results = mo.map_async_parallelizer(variable_inputs,
															 mo.function_helper,
															 cpu_cores,
															 show_progress=True)

				# Collect variable and non-variable DNA MSA files
				variable_dna = [r[1] for r in variable_results if r[0] is True]
				non_variable_dna = [r[1] for r in variable_results if r[0] is False]
			else:
				print('\nDNA MSAs were requested, but input files contain protein sequences.')

	# User only wants the loci MSAs or the input is a folder with FASTA files
	if only_loci_msas or os.path.isdir(input_file):
		print(f'MSAs for each input locus/file are available in {mafft_outdir}')
		sys.exit(0)

	# Create folder to store sample MSAs
	sample_msas_outdir = fo.join_paths(output_directory, ['sample_MSAs'])
	fo.create_directory(sample_msas_outdir)

	# Create the full protein MSA
	print('\nCreating file with the full protein MSA...')
	sample_protein_msas_outdir = fo.join_paths(sample_msas_outdir, ['protein'])
	fo.create_directory(sample_protein_msas_outdir)
	# Get loci IDs for which a gapped protein MSA was created
	gapped_results_loci_ids = [fo.file_basename(file).split('_protein')[0] for file in gapped_results]
	# Concatenate all loci MSAs and index with BioPython to retrieve sequences easily
	loci_msa_concat = fo.join_paths(sample_protein_msas_outdir, ['concat.fasta'])
	# Concatenate all alignment files and index with BioPython
	fo.concatenate_files(gapped_results, loci_msa_concat)

	# Index file with SeqIO.index_db to store record information as a file on disk
	# This allows to reload the index with multiprocessing
	# SeqIO.index creates the index in memory which cannot be shared between processes
	loci_msa_concat_index = fao.index_fasta(loci_msa_concat, True)

	# Pass the FASTA file used to create the index to each process and the function reloads the index for each process
	common_args = [loci_msa_concat, sample_protein_msas_outdir, gapped_results_loci_ids, input_file]
	inputs = [[sid, range(1,sample_ids.index(sid)+1)] for sid in sample_ids]
	inputs = im.multiprocessing_inputs(inputs, common_args, create_full_msa)
	sample_protein_MSA_outfiles = mo.map_async_parallelizer(inputs,
															mo.function_helper,
															cpu_cores,
															show_progress=True)

	# Concatenate sample protein alignments
	full_protein_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_PROTEIN_MSA])
	fo.concatenate_files(sample_protein_MSA_outfiles, full_protein_alignment)

	# Get length of the full alignment
	protein_msa_length = len(fo.read_lines(full_protein_alignment, strip=True, num_lines=2)[1])
	print(f'Protein MSA length: {protein_msa_length}')

	# Delete file indexed to get individual sequences
	fo.remove_files([loci_msa_concat])

	if not no_cleanup:
		fo.delete_directory(sample_protein_msas_outdir)

	if dna_msa:
		sample_dna_msas_outdir = fo.join_paths(sample_msas_outdir, ['dna'])
		fo.create_directory(sample_dna_msas_outdir)

		# Concatenate all loci MSAs and index with BioPython to retrieve sequences easily
		loci_msa_concat = fo.join_paths(sample_dna_msas_outdir, ['concat.fasta'])
		# Concatenate all alignment files and index with BioPython
		fo.concatenate_files(dna_results, loci_msa_concat)

		# Index file with SeqIO.index_db to store record information as a file on disk
		# This allows to reload the index with multiprocessing
		# SeqIO.index creates the index in memory which cannot be shared between processes
		loci_msa_concat_index = fao.index_fasta(loci_msa_concat, True)

		# Create the full DNA MSA
		print('Creating file with the full DNA MSA...')
		common_args = [loci_msa_concat, sample_dna_msas_outdir, gapped_results_loci_ids, input_file]
		inputs = [[sid, range(1,sample_ids.index(sid)+1)] for sid in sample_ids]
		inputs = im.multiprocessing_inputs(inputs, common_args, create_full_msa)
		sample_dna_MSA_outfiles = mo.map_async_parallelizer(inputs,
															mo.function_helper,
															cpu_cores,
															show_progress=True)

		# Concatenate sample dna alignments
		full_dna_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_DNA_MSA])
		fo.concatenate_files(sample_dna_MSA_outfiles, full_dna_alignment)

		# Get length of the full alignment
		dna_msa_length = len(fo.read_lines(full_dna_alignment, strip=True, num_lines=2)[1])
		print(f'DNA MSA length: {dna_msa_length}')

		# Delete folder with intermediate sample MSAs
		if not no_cleanup:
			fo.delete_directory(sample_dna_msas_outdir)

	# Create full MSAs with only variable positions if requested
	if output_variable:
		sample_protein_variable_msas_folder = fo.join_paths(sample_msas_outdir, ['variable_protein'])
		fo.create_directory(sample_protein_variable_msas_folder)

		# Get loci IDs for which a gapped protein MSA was created and had variable positions
		variable_results_loci_ids = [fo.file_basename(file).split('_protein')[0] for file in variable_protein]

		# Concatenate all loci MSAs and index with BioPython to retrieve sequences easily
		loci_msa_concat = fo.join_paths(sample_protein_variable_msas_folder, ['concat.fasta'])
		# Concatenate all alignment files and index with BioPython
		fo.concatenate_files(variable_protein, loci_msa_concat)

		# Index file with SeqIO.index_db to store record information as a file on disk
		# This allows to reload the index with multiprocessing
		# SeqIO.index creates the index in memory which cannot be shared between processes
		loci_msa_concat_index = fao.index_fasta(loci_msa_concat, True)

		# Create the full protein MSA for the variable positions
		print('\nCreating file with the full protein MSA for the variable positions...')
		common_args = [loci_msa_concat, sample_protein_variable_msas_folder, variable_results_loci_ids, input_file]
		inputs = [[sid, range(1,sample_ids.index(sid)+1)] for sid in sample_ids]
		inputs = im.multiprocessing_inputs(inputs, common_args, create_full_msa)
		sample_protein_variable_MSA_outfiles = mo.map_async_parallelizer(inputs,
															mo.function_helper,
															cpu_cores,
															show_progress=True)

		# Concatenate sample protein alignments
		full_variable_protein_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_PROTEIN_MSA_VARIABLE])
		fo.concatenate_files(sample_protein_variable_MSA_outfiles, full_variable_protein_alignment)

		# Get length of the full alignment
		protein_variable_msa_length = len(fo.read_lines(full_variable_protein_alignment, strip=True, num_lines=2)[1])
		print(f'Protein variable MSA length: {protein_variable_msa_length}')

		# Delete folder with intermediate sample protein variable MSAs
		if not no_cleanup:
			fo.delete_directory(sample_protein_variable_msas_folder)

		if dna_msa:
			sample_dna_variable_msas_folder = fo.join_paths(sample_msas_outdir, ['variable_dna'])
			fo.create_directory(sample_dna_variable_msas_folder)

			# Get loci IDs for which a gapped DNA MSA was created and had variable positions
			variable_results_loci_ids = [fo.file_basename(file, False) for file in variable_dna]

			# Concatenate all loci MSAs and index with BioPython to retrieve sequences easily
			loci_msa_concat = fo.join_paths(sample_dna_variable_msas_folder, ['concat.fasta'])
			# Concatenate all alignment files and index with BioPython
			fo.concatenate_files(variable_dna, loci_msa_concat)

			# Index file with SeqIO.index_db to store record information as a file on disk
			# This allows to reload the index with multiprocessing
			# SeqIO.index creates the index in memory which cannot be shared between processes
			loci_msa_concat_index = fao.index_fasta(loci_msa_concat, True)

			# Create the full DNA MSA for the variable positions
			print('\nCreating file with the full DNA MSA for the variable positions...')
			common_args = [loci_msa_concat, sample_dna_variable_msas_folder, variable_results_loci_ids, input_file]
			inputs = [[sid, range(1,sample_ids.index(sid)+1)] for sid in sample_ids]
			inputs = im.multiprocessing_inputs(inputs, common_args, create_full_msa)
			sample_dna_variable_MSA_outfiles = mo.map_async_parallelizer(inputs,
																mo.function_helper,
																cpu_cores,
																show_progress=True)

			# Concatenate sample dna alignments
			full_variable_dna_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_DNA_MSA_VARIABLE])
			fo.concatenate_files(sample_dna_variable_MSA_outfiles, full_variable_dna_alignment)

			# Get length of the full alignment
			dna_variable_msa_length = len(fo.read_lines(full_variable_dna_alignment, strip=True, num_lines=2)[1])
			print(f'Protein MSA length: {dna_variable_msa_length}')

			# Delete folder with intermediate sample DNA variable MSAs
			if not no_cleanup:
				fo.delete_directory(sample_dna_variable_msas_folder)

	# Delete folders with gapped MSAs for full and only variable positions
	if not no_cleanup:
		fo.delete_directory(mafft_outdir)
		fo.delete_directory(sample_msas_outdir)
