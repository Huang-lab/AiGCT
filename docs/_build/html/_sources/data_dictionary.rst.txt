Data Dictionary
===============

We describe the logical data model of our database using standard relational 
database terminology such as table, column, and primary key.
The primary key is the collection of columns in a table that uniquely
identify a row in that table. Columns tagged with a (PK) symbol are part
of the primary key of the table.

Table: VARIANT_TASK
-------------------

Master list of tasks (categories of variant effect data).

Columns:

* CODE (PK) - Unique code for each task
* NAME
* SOURCE_TYPE
* DESCRIPTION

Table: VARIANT_DATA_SOURCE
--------------------------

Master list of sources of data.

Columns:

* CODE (PK) - Unique code for each source
* NAME
* DESCRIPTION

Table: VARIANT_EFFECT_SOURCE
----------------------------

Master list of sources of variant effect score data. Currently all score data
originates from a VEP.

Columns:

* CODE (PK) - Unique code for each source
* NAME
* SOURCE_TYPE - Currently, all values are VEP
* DESCRIPTION

Table: VARIANT
--------------

Master list of variants.

Columns:

* GENOME_ASSEMBLY (PK) - Currently, hg38 for all rows.
* CHROMOSOME (PK)
* POSITION (PK)
* REFERENCE_NUCLEOTIDE (PK)
* ALTERNATE_NUCLEOTIDE (PK)
* PRIOR_GENOME_ASSEMBLY
* PRIOR_CHROMOSOME
* PRIOR_POSITION
* PRIOR_PRIOR_GENOME_ASSEMBLY
* PRIOR_PRIOR_CHROMOSOME
* PRIOR_PRIOR_POSITION
* REFERENCE_AMINO_ACID
* ALTERNATE_AMINO_ACID
* AMINO_ACID_POSITION
* RS_DBSNP
* GENE_SYMBOL
* ENSEMBL_GENE_ID
* ENSEMBL_TRANSCRIPT_ID
* ENSEMBL_PROTEIN_ID
* ALLELE_FREQUENCY_SOURCE
* ALLELE_FREQUENCY

Table: VARIANT_EFFECT_LABEL
---------------------------

Variant labels organized by task.

Columns:

* TASK_CODE (PK)
* GENOME_ASSEMBLY (PK) - Currently, hg38 for all rows.
* CHROMOSOME (PK)
* POSITION (PK)
* REFERENCE_NUCLEOTIDE (PK)
* ALTERNATE_NUCLEOTIDE (PK)
* LABEL_SOURCE
* RAW_LABEL
* BINARY_LABEL

Table: VARIANT_EFFECT_SCORE
---------------------------

Variant effect prediction scores organized by task.

Columns:

* TASK_CODE (PK)
* GENOME_ASSEMBLY (PK) - Currently, hg38 for all rows.
* CHROMOSOME (PK)
* POSITION (PK)
* REFERENCE_NUCLEOTIDE (PK)
* ALTERNATE_NUCLEOTIDE (PK)
* SCORE_SOURCE
* RAW_SCORE
* RANK_SCORE - Normalized score, i.e. value between 0-1

Table: VARIANT_FILTER
---------------------

Table of system supplied named filters that can be used by the user to restrict the set
of variants to be retrieved or to be used in a benchmarking analysis. Filters are
defined in terms of a list of genes and/or variants to be included or excluded.

Columns:

* CODE (PK) - Unique code for each filter
* NAME
* INCLUDE_GENES - Y/N, include or exclude genes listed in VARIANT_FILTER_GENE table.
* INCLUDE_VARIANTS - Y/N, include or exclude variants listed in VARIANT_FILTER_VARIANT table.

Table: VARIANT_FILTER_GENE
--------------------------

Lists the set of genes to be included or excluded by a filter in the
VARIANT_FILTER table.

Columns:

* FILTER_CODE (PK) - Refers to a CODE in the variant_filter table.
* GENE_SYMBOL (PK)

Table: VARIANT_FILTER_VARIANT
-----------------------------

Lists the set of variants to be included or excluded by a filter in the
VARIANT_FILTER table.

Columns:

* FILTER_CODE (PK) - Refers to a CODE in the variant_filter table.
* GENOME_ASSEMBLY (PK)
* CHROMOSOME (PK)
* POSITION (PK)
* REFERENCE_NUCLEOTIDE (PK)
* ALTERNATE_NUCLEOTIDE (PK)




