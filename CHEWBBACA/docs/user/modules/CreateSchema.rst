CreateSchema - Create a gene-by-gene schema
===========================================

What is a Schema and how is it defined
::::::::::::::::::::::::::::::::::::::

A Schema is a pre-defined set of loci that is used in MLST analyses. Traditional MLST schemas relied in 7 loci that were internal fragments of housekeeping genes and each locus was defined by its amplification by a pair of primers yielding a fragment of a defined size.

In genomic analyses, schemas are a set of loci that are:

- Present in the majority of strains for core genome (cg) MLST schemas, typically a threshold of presence in 95% of the strains is used in schema creation. The assumption is that in each strain up to 5% of loci may not be identified due to sequencing coverage problems, assembly problems or other issues related to the use of draft genome assemblies.

- Present in at least one of the analyzed strains in the schema creation for pan genome/whole  genome (pg/wg) MLST schemas.

- Present in less than 95% of the strains for accessory genome (ag) MLST schemas.

It is important to consider that these definitions are always operational in nature, in the sense that the analyses are performed on a limited number of strains representing part of the biological diversity of a given species or genus and are always dependent on the definition of thresholds.  

In most cg/wg/pg/ag MLST schemas, contrary to MLST schemas, each locus corresponds to a coding sequence (CDS). However, depending on the allele calling algorithm, the alleles called for a given locus can be CDSs or best matches to existing CDSs without enforcing the need for the identified allele to be a CDS.  

In **chewBBACA**, schemas are composed of loci defined by CDSs and, by default, all the called alleles of a given locus are CDSs as defined by `Prodigal <https://github.com/hyattpd/Prodigal>`_ (for chewBBACA<=3.2.0) or `Pyrodigal <https://github.com/althonos/pyrodigal>`_ (for chewBBACA>=3.3.0) (it is also possible to provide FASTA files with CDSs and the ``--cds`` parameter to skip the gene prediction step). The use of Prodigal or Pyrodigal, instead of simply ensuring the presence of start and stop codons, adds an extra layer of confidence in identifying the most probable CDS for each allele. Because of this approach there may be variability in the size of the alleles identified by **chewBBACA** and by default a threshold of +/-20% of the mode of the size of the alleles of a given locus is used to identify a locus as present.

Create a wgMLST schema
::::::::::::::::::::::

Given a set of genome assemblies in FASTA format, chewBBACA offers the option to create a new schema by defining the distinct loci present in the genomes.

Basic Usage
-----------

::

	chewBBACA.py CreateSchema -i /path/to/InputAssembliesFolder -o /path/to/OutputFolder --n SchemaName --ptf /path/to/ProdigalTrainingFile

.. important::
	You should adjust the value passed to the ``--cpu`` parameter based on the specifications of your machine. chewBBACA will automatically adjust the value if it matches or exceeds the number of available CPU cores.

.. important::
	The use of a prodigal training file for schema creation is highly recommended. The training file is included in the newly created schema and is used to predict genes during the allele calling process.

Parameters
----------

::

    -i, --input-files           (Required) Path to the directory that contains the input FASTA files
                                or to a file with a list of full paths to FASTA files, one per line.

    -o, --output-directory      (Required) Output directory where the process will store intermediate
                                files and create the schema's directory.

    --n, --schema-name          (Optional) Name given to the schema folder (default: schema_seed).

    --ptf, --training-file      (Optional) Path to the Prodigal training file used by Pyrodigal to predict
                                genes and added to the schema folder (default: None).

    --bsr, --blast-score-ratio  (Optional) BLAST Score Ratio (BSR) value. The BSR is computed for each
                                BLASTp alignment and aligned sequences with a BSR >= than the defined
                                value are considered to be alleles of the same gene (default: 0.6).

    --l, --minimum-length       (Optional) Minimum sequence length value. Predicted coding sequences (CDSs)
                                shorter than this value are excluded (default: 201).

    --t, --translation-table    (Optional) Genetic code used to predict genes and to translate coding
                                sequences (CDSs) (default: 11).

    --st, --size-threshold      (Optional) Coding sequence (CDS) size variation threshold. Added to the
                                schema's config file to identify alleles with a size that deviates from
                                the locus length mode during the allele calling process (default: 0.2).

    --cpu, --cpu-cores          (Optional) Number of CPU cores that will be used to run the process (chewie
                                resets to a lower value if it is equal to or exceeds the total number of
                                available CPU cores)(default: 1).

    --b, --blast-path           (Optional) Path to the directory that contains the BLAST executables (default:
                                assumes BLAST executables were added to PATH).

    --pm, --prodigal-mode       (Optional) Prodigal running mode ("single" for finished genomes, reasonable
                                quality draft genomes and big viruses. "meta" for metagenomes, low quality
                                draft genomes, small viruses, and small plasmids) (default: single).

    --cds, --cds-input          (Optional) If provided, chewBBACA skips the gene prediction step and assumes
                                the input FASTA files contain coding sequences (default: False).
		
    --no-cleanup                (Optional) If provided, intermediate files generated during process execution
                                are not deleted at the end (default: False).

.. important::
  If you provide the ``--cds-input`` parameter, chewBBACA assumes that the input FASTA files contain CDSs and skips the gene prediction step. To avoid issues related with the format of the sequence headers, chewBBACA renames the sequence headers based on the unique basename prefix determined for each input file and on the order of the CDSs (e.g.: CDSs inside a file named ``GCF_000007125.1_ASM712v1_cds_from_genomic.fna`` are renamed to ``GCF_000007125-protein1``, ``GCF_000007125-protein2``, ..., ``GCF_000007125-proteinN``).

Outputs
-------

::

	OutputFolderName
	├── SchemaName
	│   ├── short
	│   │   ├── GenomeID_proteinN_short.fasta
	│   │   ├── ...
	│   │   └── GenomeID_proteinN_short.fasta
	│   ├── GenomeID_proteinN.fasta
	│   ├── ...
	│   ├── GenomeID_proteinN.fasta
	│   └── Training_file.trn
	├── invalid_cds.txt
	└── cds_coordinates.tsv

- One FASTA file per distinct locus identified in the schema creation process in the ``OutputFolderName/SchemaName`` directory. The name attributed to each FASTA file in the schema is based on the genome of origin of the representative allele chosen for that locus and on the order of gene prediction (e.g.: ``GCA-000167715-protein12.fasta``, first allele for the locus was identified in a genome assembly with the prefix ``GCA-000167715`` and the locus was the 12th gene predicted in that assembly).

- The ``OutputFolderName/SchemaName`` directory also contains a directory named ``short`` that includes FASTA files with the representative alleles for each locus.

- The training file passed to create the schema is also included in ``OutputFolderName/SchemaName`` and will be automatically detected during the allele calling process.

- The ``cds_coordinates.tsv`` file contains the coordinates (genome unique identifier, contig identifier, start position, stop position, protein identifier attributed by chewBBACA, and coding strand (chewBBACA<=3.2.0 assigns 1 to the forward strand and 0 to the reverse strand and chewBBACA>=3.3.0 assigns 1 and -1 to the forward and reverse strands, respectively)) of the CDSs identified in each genome. 

- The ``invalid_cds.txt`` file contains the list of alleles predicted that were excluded based on the minimum sequence size value and presence of ambiguous bases.

Workflow of the CreateSchema module
:::::::::::::::::::::::::::::::::::

.. image:: /_static/images/CreateSchema.png
   :width: 1000px
   :align: center

The CreateSchema module creates a schema seed based on a set of FASTA files with genome assemblies or CDSs. Brief description of the workflow:

	- If genome assemblies are given, the process starts by predicting CDSs for each genome using Pyrodigal. Alternatively, if FASTA files containing CDSs are provided (``--cds`` parameter), the process skips the gene prediction step.

	- The CDSs identified in the input files are deduplicated and translated (CDS translation identifies and excludes CDSs that contain ambiguous bases and with length below the theshold defined by the ``--l`` parameter), followed by a second deduplication step to determine the set of distinct translated CDSs. The schema creation process creates the same hash tables used in the allele calling process to store information about the distinct CDSs and distinct translated CDSs.

	- The distinct translated CDSs are clustered based on the proportion of minimizers (minimizers selected based on lexicographic order, k=5, w=5) shared with representative CDSs. The largest or one of the largest CDSs is selected as the first representative CDS. New representative CDSs are selected when CDSs share a low proportion (<0.2) of minimizers with any of the chosen representative CDSs.

	- Non-representative CDSs that share a proportion of minimizers ≥ 0.9 with the cluster representative are considered to correspond to the same locus and are excluded from the analysis.

	- The proportion of shared minimizers between non-representative CDSs is determined to exclude CDSs sharing a proportion of minimizers ≥ 0.9 with larger CDSs.

	- Intracluster and intercluster alignment with BLASTp enable identifying and excluding CDSs similar to representative or larger non-representative CDSs based on a BLAST Score Ratio (BSR) ≥ 0.6.

	- Each remaining distinct CDS is considered to be an allele of a distinct locus. The process ends by creating a schema seed, which includes one FASTA file containing a single representative allele per distinct locus identified in the analysis.
