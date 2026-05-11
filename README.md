# SOLU INSTRUCTIONS
## Deployment
Build image locally
`docker build . --tag solu/chewBBACA:$VERSION`
## Release
```
gcloud auth login
./release.sh -v <VERSION>
gcloud auth revoke
```

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

## 3.5.4 - 2026-04-22

This release adds the following bug fixes:

- The file with the cgMLST MSA was missing from the output directory created by the `AlleleCallEvaluator` module since the release of the `ComputeMSA` module (thanks to @victor5lm who reported this issue in https://github.com/B-UMMI/chewBBACA/issues/235). This change was unintended. The file with cgMLST MSA was re-added as one of the output files from the `AlleleCallEvaluator` module (`protein_msa.fasta` file). Simplified the conditions used to determine which steps to run.
- Made the definition of the sequence headers for the adapted schemas created by the `PrepExternalSchema` less rigid to better deal with unexpected headers. Added tests to validate schema adaptation for multiple external schema formats (e.g., [EnteroBase](https://enterobase.warwick.ac.uk/), [PubMLST](https://pubmlst.org/), and [Ridom](https://www.cgmlst.org/ncs)).
- Fixed the parsing of the file with loci annotations in the `LoadSchema` module.
- Fixed the creation of the paths to intermediate FASTA files used by BLASTp to determine representative alleles during schema adaptation. The paths were not properly created when the `SyncSchema` module called the `PrepExternalSchema` module using relative paths.
- Changed the index value used by the `select_highest_scores` function to sort and select the highest scoring BLASTp matches per target from `5` to `6`. The results were being sorted based on the length of the target sequences instead of the alignment raw score. This issue would not allow to identify the best scoring alignment in some cases, potentially leading to some sequences not being classified if a lower scoring alignment was selected. Fixing this issue results in a slight increase in the accuracy of the allele calling (thanks to @andreaderuvo for reporting this issue in https://github.com/B-UMMI/chewBBACA/issues/234).

Check our [Changelog](https://github.com/B-UMMI/chewBBACA/blob/master/CHANGELOG.md) to learn more about the latest changes.

## Citation

When using chewBBACA, please use the following citation:

> Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)
