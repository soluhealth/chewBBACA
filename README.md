
[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/chewBBACA/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/chewbbaca)
[![Conda](https://img.shields.io/conda/dn/bioconda/chewbbaca?color=green)](https://anaconda.org/bioconda/chewbbaca)
[![chewBBACA](https://github.com/B-UMMI/chewBBACA/workflows/chewbbaca/badge.svg)](https://github.com/B-UMMI/chewBBACA/actions?query=workflow%3Achewbbaca)
[![Documentation Status](https://readthedocs.org/projects/chewbbaca/badge/?version=latest)](https://chewbbaca.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/github/license/B-UMMI/chewBBACA)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI:10.1099/mgen.0.000166](https://img.shields.io/badge/DOI-10.1099%2Fmgen.0.000166-blue)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166)

# chewBBACA

**chewBBACA** is a software suite for the creation and evaluation of core genome and whole genome MultiLocus Sequence 
Typing (cg/wgMLST) schemas and results. The "BBACA" stands for "BSR-Based Allele Calling Algorithm". BSR stands for 
BLAST Score Ratio as proposed by [Rasko DA et al.](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2). The "chew" part adds extra coolness to the name and could be thought of as "Comprehensive and Highly Efficient Workflow". chewBBACA allows to define the target loci in a schema based on multiple genomes (e.g. define target loci based on the distinct loci identified in a dataset of high-quality genomes for a species or lineage of interest) and performs allele calling to determine the allelic profiles of bacterial strains, easily scaling to thousands of genomes with modest computational resources. chewBBACA includes functionalities to annotate the schema loci, compute the set of loci that constitute the core genome for a given dataset, and generate interactive reports for schema and allele calling results evaluation to enable an intuitive analysis of the results in surveillance and outbreak detection settings or population studies. Pre-defined cg/wgMLST schemas can be downloaded from [Chewie-NS ](https://chewbbaca.online/) or adapted from other cg/wgMLST platforms.

### Check the [documentation](https://chewbbaca.readthedocs.io/en/latest/index.html) for implementation details and guidance on using chewBBACA.

## News

## 3.5.1 - 2026-01-06

chewBBACA no longer checks if input files have unique basename prefixes shorter than 30 characters. In the past, this was performed to ensure that sequence identifiers did not exceed the character limit (50 characters) enforced by BLAST when creating a database. The main changes to file name processing are the following:

- chewBBACA uses the file basename without the file extension as unique identifier (e.g. `GCF_008632635.1.fasta` is converted to `GCF_008632635.1`), instead of trying to determine the shortest unique prefix that can be used to identify each input file. It is still necessary for each file to have a unique identifier after the removal of the file extension (e.g. `GCF_008632635.1.fasta` and `GCF_008632635.1.fna` have different file extensions but the same identifier after removing the file extension, which is not allowed).
- The CreateSchema module uses the input file basenames without the file extension to define the identifiers for the loci in the created schemas (e.g. loci initially identified in the genomes `GCF_008632635.1.fasta` and `GCA_000006785.2_ASM678v2.fasta` are named as `GCF_008632635.1-proteinN.fasta` and `GCA_000006785.2_ASM678v2-proteinN.fasta`, respectively). We still recommend using short and unique file names without special characters (e.g.: `!@#?$^*()+`) for conciseness and to avoid potential issues.
- The AlleleCall module accepts and uses the new loci identifier format used by the CreateSchema module. The input genome or CDS files can also have basenames of any length as long as the basename without the file extension for each input file is unique. The output files created by the AlleleCall module use the full unique basenames (e.g. for the genome `GCA_000006785.2_ASM678v2.fasta`, the genome identifier used in the output files will be `GCA_000006785.2_ASM678v2`, instead of `GCA_000006785` used up until chewBBACA v3.5.0).
- The PrepExternalSchema module accepts schemas containing loci FASTA files with basenames longer than 30 characters.

Additionally, the CDS identifiers are converted to a different format (`lcl|SEQ1`, `lcl|SEQ2`...`lcl|SEQN`) before creating a BLAST database with `makeblastdb` and the `-parse_seqids` option to avoid issues related to some sequence identifiers being interpretd and modified (e.g. interpretd as PDB Chain IDs) when creating a database, resulting in errors when an identifier is modified and no longer matches the original identifier. This allowed to remove the check to verify that unique prefixes are not modified by BLAST during database creation.

Additional changes:

- Added the `--output-masked` option to the AlleleCall module to create a TSV file with the masked profiles (`INF-` prefixes are removed and the NIPH, NIPHEM, ASM, ALM, PLOT3, PLOT5, LOTSC, and PAMA classes are converted to `0`).

Check our [Changelog](https://github.com/B-UMMI/chewBBACA/blob/master/CHANGELOG.md) to learn about the latest changes.

## Citation

When using chewBBACA, please use the following citation:

> Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)
