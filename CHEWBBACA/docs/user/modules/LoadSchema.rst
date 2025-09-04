LoadSchema - Upload a schema to Chewie-NS
=========================================

The *LoadSchema* module enables the upload of schemas to Chewie-NS instances.

.. important::
    **You need to be a registered user with Contributor privileges to be able to upload schemas!** If you have registered at the `chewie-NS public website <https://chewbbaca.online/auth>`_ and want to contribute new schemas please contact us via e-mail at: imm-bioinfo@medicina.ulisboa.pt

To upload a schema to Chewie-NS it is required to provide:

- The **path** to the local schema.

  - The schema must have been created with chewBBACA v2.5.0 or a later version. If your schema was created with an older version, please adapt the schema with the ``PrepExternalSchema`` process or run the ``AlleleCall`` process to convert the schema to the latest chewBBACA version.

.. warning::
	**Only schemas that have been used with the same valid value per parameter can be uploaded (this restriction applies to the BLAST Score Ratio, Prodigal training file, minimum  sequence length, genetic code and sequence size variation threshold parameters).** Invalid or multiple values used in different allele calling runs for a single parameter can lead to inconsistent results; thus, it is strongly advised to always perform allele calling with the same set of parameters and refrain from altering the initial set of parameter values defined in the schema creation or adaptation processes.

- The **species ID** or **scientific name** of the species that the schema will be associated to.
  
	- To know the ID of a species you can:
		- Consult the `Overview <https://chewbbaca.online/stats>`_ table in the Chewie-NS public website.
		- Query the ``/species/list`` API endpoint through `Swagger UI <https://chewbbaca.online/api/NS/api/docs>`_.
		- Use a simple curl command: ``e.g.: curl -X GET "https://chewbbaca.online/NS/api/species/list" -H  "accept: application/json"``.
	
	Alternatively, you can use the `NSStats <https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/CHEWBBACA_NS/stats_requests.py>`_  module to get information about the list of species in Chewie-NS.

- A **name** for the schema.

  - The name should be short and concise. The name **must be unique** among the set of names for schemas of the same species (this means that the using the same name of an existing schema will lead to an error) and should not include spaces (e.g.: ``Project_cgMLST``, ``SRA_wgMLST``, ``Organization_cgMLST``).

- A **prefix** for the loci identifiers to facilitate the identification of the schema they belong to.

  - You may use the name of the schema as prefix to ensure prefix uniqueness for the loci of a schema.

Users may provide a description about the schema. The file with the description  will be sent to Chewie-NS and displayed in the schema's page in the website. Markdown syntax is  supported in order to allow greater customizability of the rendered description. For more information on the Markdown specification accepted by Chewie-NS please visit the `Github Flavored Markdown Specification page <https://github.github.com/gfm/>`_.

Sample description:

::

    # Whole-genome MLST schema for *Species name*

    This schema was created with [chewBBACA 2.5.0](https://github.com/B-UMMI/chewBBACA).

    ## Schema creation and validation

    A total of 100 *Species name* genomes were used to create this wgMLST schema.
    (Add more information that might be relevant to reproduce the process.
    For instance the assembly pipeline used or links to relevant external data repositories)

    ## Dataset

    All raw reads from SRA annotated as *Species name* were assembled into 100 genomes.
    (Add more information about dataset creation/collection. Include date of data download.
    A list of accession numbers or links to relevant external data repositories maybe useful)

    ## Citations

    (Add any relevant citations)

    For more information please access [external page](https://external/page)
    (If there is any external source with more information, link it here)

The process queries UniProt's SPARQL endpoint to retrieve annotations for the loci in the schema. The user that uploads the schema can provide a TSV file with annotations for some or all loci in the schema. The file with annotations must have the following columns:

	- First column: locus identifier (name of locus file without the ``.fasta`` extension).
	- Second column: UniProt protein name.
	- Third column: UniProt gene name.
	- Fourth column: URI for the UniProt entry.
	- Fifth column: user annotation (gene name commonly attributed by the user).
	- Sixth column: custom annotation (another term that the user might want to attribute).

However, no headers are necessary.

.. rst-class:: align-center

	+--------+-----------------------------------------+-----------------+------------------------------------------------+-------+---------------------------------------------------------------------------------+
	| locus1 | Penicillin-binding protein 1B, putative |     SAG0159     | https://www.uniprot.org/uniprotkb/Q8E240/entry | pbp1b | Penicillin-binding protein 1B                                                   |
	+--------+-----------------------------------------+-----------------+------------------------------------------------+-------+---------------------------------------------------------------------------------+
	| locus1 | Penicillin-binding protein 2A           |      pbp2A      | https://www.uniprot.org/uniprotkb/Q8DWZ3/entry | pbp2a | Penicillin-binding protein 2A                                                   |
	+--------+-----------------------------------------+-----------------+------------------------------------------------+-------+---------------------------------------------------------------------------------+
	| locus1 | Cell wall surface anchor family protein |     SAG1408     | https://www.uniprot.org/uniprotkb/Q8DYR5/entry | PI-2a | PI-2a ancillary protein 1                                                       |
	+--------+-----------------------------------------+-----------------+------------------------------------------------+-------+---------------------------------------------------------------------------------+
	| locus1 | Laminin-binding surface protein         |       lmb       | https://www.uniprot.org/uniprotkb/Q8DZ80/entry | lmb   | metal ABC transporter substrate-binding lipoprotein/laminin-binding adhesin Lmb |
	+--------+-----------------------------------------+-----------------+------------------------------------------------+-------+---------------------------------------------------------------------------------+	

It is not necessary to provide all annotation types for each locus nor for every locus. If no information is provided N/A will be automatically shown in the locus details page in Chewie-NS.

Basic Usage
-----------

To upload a schema for *Escherichia coli*, we could run one of the following commands:

- Providing the species ID:

::

	chewBBACA.py LoadSchema -i path/to/SchemaFolder -sp 9 -sn cgMLST_95 -lp cgMLST_95

- Providing the species name:

::

	chewBBACA.py LoadSchema -i path/to/SchemaFolder -sp "Escherichia coli" -sn cgMLST_95 -lp cgMLST_95

To upload a schema and provide a description and annotations:

::

    chewBBACA.py LoadSchema -i path/to/SchemaFolder -sp 9 -sn cgMLST_95 -lp cgMLST_95 --df description.txt --a annotations.tsv

To continue an upload that was interrupted or that aborted, we should provide the command used in the process that failed and add the ``--continue_up`` argument

::

    chewBBACA.py LoadSchema -i path/to/SchemaFolder -sp 9 -sn cgMLST_95 -lp cgMLST_95 --continue_up

.. important::
	If you cannot complete schema upload or if the information in the website is incorrect or missing, please contact us via e-mail: imm-bioinfo@medicina.ulisboa.pt

Parameters
----------

::

    -i, --schema-directory      (Required) Path to the directory of the schema to upload.

    -sp, --species-id           (Required) The integer identifier or name of the species that the
                                schema will be associated to in Chewie-NS.

    -sn, --schema-name          (Required) A brief and meaningful name that should help understand
                                the type and content of the schema.

    -lp, --loci-prefix          (Required) Prefix included in the name of each locus of the schema.

    --df, --description-file    (Optional) Path to a text file with a description about the schema.
                                Markdown syntax is supported in order to offer greater customizability
                                of the rendered description in the Frontend. Will default to the schema's
                                name if the user does not provide a valid path for a file (default: None).

    --a, --annotations          (Optional) Path to a TSV file with loci annotations. The first column has
                                loci identifiers (w/o .fasta extension), the second has UniProt protein names,
                                the third has UniProt gene names, the fourth has UniProt URIs, the fifth has
                                user annotations, and the sixth has custom annotations. (default: None).

    --cpu, --cpu-cores          (Optional) Number of CPU cores/threads that will be used to run the process
                                (chewie resets to a lower value if it is equal to or exceeds the total
                                number of available CPU cores/threads). This value is used to
                                accelerate the quality control step that validates schema alleles (default: 1).

    --ns, --nomenclature-server (Optional) The base URL for the Chewie-NS instance. The default value,
                                "main", will establish a connection to "https://chewbbaca.online/",
                                "tutorial" to "https://tutorial.chewbbaca.online/" and "local" to
                                "http://127.0.0.1:5000/NS/api/" (localhost). Users may also provide
                                the IP address to other Chewie-NS instances (default: main).

    --continue_up               (Optional) Check if the schema upload was interrupted and attempt
                                to continue upload (default: False).

Workflow of the LoadSchema module
:::::::::::::::::::::::::::::::::

.. image:: /_static/images/LoadSchema.png
   :width: 1200px
   :align: center

The LoadSchema module uploads local schemas to Chewie-NS. Brief description of the workflow:

	- The process starts by requesting the user credentials to ensure that the user has contributor privileges. Only contributors are allowed to upload schemas to Chewie-NS.

	- If the user is a contributor, the process checks if the species identifier provided by the user is valid and if the species is listed in Chewie-NS.

	- After this step, the process reads the schema’s configuration file to validate the schema parameter values and ensure that there is only a single value associated with each parameter. The initial validation steps are followed by the upload of the schema data to Chewie-NS.

	- The process reads the schema description, if the user provided one, or uses the schema name as description.

	- The alleles are translated and annotation terms for the loci are obtained through UniProt’s SPARQL endpoint. If the user provides custom loci annotations, the process reads the file provided by the user and adds the custom annotations to the loci annotation data to send to Chewie-NS.

	- After retrieving loci annotations, the process creates the schema in Chewie-NS by sending the schema’s parameter values and the list of file hashes to validate schema files uploaded in subsequent steps.

	- The loci are created and linked to the newly created schema by sending the loci identifiers and annotations to Chewie-NS.

	- The loci FASTA files are compressed and uploaded to Chewie-NS to add the allele sequences to the database and link them to the corresponding loci.

	- The last step in the process uploads the training file in the local schema and associates it to the newly created schema in Chewie-NS. After process completion, Chewie-NS will process the data that was sent to make the schema data and statistics available through the website and the API.
