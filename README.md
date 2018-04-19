# DeMinEr
DeMinEr is a tool to process ultra-deep sequencing data in order to identify rare nucleotide substitutions such as somatic hypermutations within a given genomic region (validated on IonProton libraies from 1kb PCR products). It was developed is to provide a comprehensive mutation analysis in 'AID off-target genes' but may be used for any region where mutation rate is of the same order of magnitude than the sequencing error rate (typical range:1E-5 - 1E-2).
Input data are expected to be nucleotide count tables in csv format (row: position *i*; column: A,C,G,T) stored in the *nucleotideCounts* directory. Such tables are typically obtained after bam-to-wig step as proposed by igvtools. A sequence of the genomic region (in fasta) must be found in the working directory.

## To process data with DeMiner:
  - clone or download repository,
  - unzip *nucleotideCounts* archive to get sample data, or create *nucleotideCounts* sub-directory and copy your own data.
  - open the *RunDeMiner* Jupyter notebook and run relevant cells. The python source files *DeMinEr.py* and *DeMinErReport.py* should be either in the working directory or in the PYTHONPATH.
