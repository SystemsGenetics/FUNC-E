# FUNC-E

FUNC-E is a Perl script used for functional enrichment of gene lists. It follows a similar approach to that of [DAVID] (https://david.ncifcrf.gov/) in that it performs enrichment analysis using a Fisher's test but then clusters enriched annotations using Kappa Statistics.  FUNC-E allows the user to provide their own annotation lists. It is fully executable on a UNIX command-line.

FUNC-E can be used for functional enrichment of a gene list, a probeset list (and can convert a probeset list into a gene list prior to enrichment). The gene or probeset list can contain multiple groupings (i.e. clusters or modules) and FUNC-E will perform functional enrichment on each group separately. 

# Usage
```
perl FUNC-E.pl [options]
```
For a complete listing of options, run FUNC-E with the -h flag

```
perl FUNC-E.pl -h
```

# Example 
The following example performs functional enrichment of a list of arabidopsis genes from the TAIR10 genome assembly and annotation.  It requires that four types of files using the arguments: --background, --query_list, --terms, and terms2features.  See the section below titled "Preparing Input Files" for imore information about these types of files.  Note that the --terms and --terms2features arguments can be provided as many times as there are files.  In this example, lists of terms (provided with the --terms toption) from  AraCyc, GO (Gene Ontology), IPR (InterPro), Pfam and PO (Plant Ontology) term lists have been prepared, as well as the mappings of these term lists to the genes.  These mappings are provided using the --terms2features option.  

```
perl FUNC-E.pl \
  --background arabidopsis_thaliana.TAIR10.genes.txt \
  --query_list modules.txt \
  --outprefix modules-enrichment \
  --terms AraCyc.terms.txt \
  --terms GO.terms.txt \
  --terms IPR.terms.txt \
  --terms Pfam.terms.txt \
  --terms PO.terms.txt \
  --terms2features arabidopsis_thaliana.TAIR10.genes2AraCyc.txt \
  --terms2features arabidopsis_thaliana.TAIR10.genes2GO.txt \
  --terms2features arabidopsis_thaliana.TAIR10.genes2IPR.txt \
  --terms2features arabidopsis_thaliana.TAIR10.genes2Pfam.txt \
  --terms2features arabidopsis_thaliana.TAIR10.genes2PO.txt  \
  --ecut 0.01  \
  --preset high 
```
Additionally, the --ecut option provides a p-value cutoff for enrichment, the --preset option provides a level of stringency for clustering of enriched terms, and the --outprefix provides the a prefix which is added to every output file created by this script.

# Installation
Before using FUNC-E the [R] (https://cran.r-project.org/) software must also be installed.  Addittionally, the following Perl modules must also be installed:

* Getopt::Long
* Text::NSP::Measures::2D::Fisher::right
* List::Util 
* Math::Complex
* Math::BigFloat
* Statistics::R

Once these prerequisites are met you can install FUNC-E by placing it anywhere in your path or calling the script from anywhere using the full path.

# Preparing Input Files

Before using KINC you must prepare your files.  You will need to prepare four files:
- A genomic background file containing the full list of genes.
- A query list of genes or probesets that will be analyzed for enrichment of function.
- One or more files containing a list of functional terms used for enrichment.
- One or more files that associates functional terms with genes.

Additionally, if you provide probesets rather than genes and you would like to convert probeset IDs to genes you must also prepare a file that maps probesets to genes.  Each of these files are described below.

## Genomic Background (--background option)
FUNC-E requires the complete list of genes. This serves as the genomic "background". This file should have a single column with each gene listed on a separate line. For example to first 10 lines of a TAIR10 background file would be:

```
AT1G01010
AT1G01020
AT1G01030
AT1G01040
AT1G01046
AT1G01050
AT1G01060
AT1G01070
AT1G01073
AT1G01080
```
## Query List (--query_list option)
The query list contains the gene or probeset list that will undergo functional enrichment.  This file allows you to specify multiple groups (i.e. modules) of genes for enrichment. The file should be tab-delimited.  The first column should contain the list of genes and the second column the group (or module) name name.  The group name allows for multiple groups of genes to be listed in the same file but enrichment  performed separately for each. The second column, however, may be left blank and only a single column of gene names can be provided.  

```
AT1G01010 Module1
AT1G01020 Module1
AT1G01030 Module1
AT1G01040 Module1
AT1G01046 Module2
AT1G01050 Module2
AT1G01060 Module2
AT1G01070 Module2
AT1G01073 Module2
AT1G01080 Module2
```

## Term List (--terms options)
A term list is a file that contains the terms used for enrichment.  This file should be a tab delimited file with three columns:  term category, term name and description.  The term name must be unique (e.g. term accession).  The term list can be contained in one large file with all terms from multiple vocabularies combined, or each vocabulary can be in separate files.  The following example combines terms from multiple vocabularies into a single list:
```
GO      GO:0000005      ribosomal chaperone activity
GO      GO:0000008      thioredoxin
IPR     IPR000002       Cdc20/Fizzy
IPR     IPR000003       Retinoid X receptor
IPR     IPR000005       Helix-turn-helix, AraC type
```
## Term Mapping List (--terms2features option)
A term mapping list maps the genes in the genomic background to terms in the term list.  This file should be tab delimited and consist of two columns:  locus/transcript ID, term name. The locus/transcript ID must be present in the genomics background file and the term name must be present in in the terms list file(s).  The following is example lines from a mapping list for rice genes:

```
LOC_Os01g01010  GO:0005097  
LOC_Os01g01010  GO:0005622  
LOC_Os01g01010  GO:0032313   
LOC_Os01g01030  GO:0005507  
```
## Other Files
If you provide probesets rather than genes and you want to convert probesets to genes for the enrichment, you must also provide a file that maps probesets to genes.  This file should be tab delimited and consist of three columns: probeset name, locus name, and mapping weight.  The maping weight provides a measure of the strength of the mapping with respect to redundancy and ambiguity.  Set the weight value to '1' if it is not needed.  If this file is not provided then it is assumed a 1:1 mapping between a term and gene. For example:

```
Os.32622.2.S1_x_at      LOC_Os01g01170  0.5
Os.33296.2.S1_at        LOC_Os01g01180  0.636363636363636
Os.33296.1.S1_at        LOC_Os01g01190  0.5
Os.33296.1.S1_x_at      LOC_Os01g01190  0.5    
```
