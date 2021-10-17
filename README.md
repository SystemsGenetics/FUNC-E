[![DOI](https://zenodo.org/badge/73501363.svg)](https://zenodo.org/badge/latestdoi/73501363)

# FUNC-E

FUNC-E is a Python package for functional enrichment analysis of gene lists. It follows a similar approach to that of [DAVID] (https://david.ncifcrf.gov/) in that it performs enrichment analysis using a Fisher's test but then clusters enriched annotations using Kappa Statistics.  FUNC-E provides the following benefits:

1. FUNC-E provides a command-line tool for inclusion in workflows.
2. FUNC-E provides an Application Programmers Interface (API) that can be used to incorporate functional enrichment into any Python script or application.
3. You can provide any list of vocabularies for which you have annotations (e.g., GO, KEGG, InterPro, Pfam, etc.)
4. FUNC-E is species agnostic.  You provide the gene names and the functional annotations.
5. FUNC-E provides both a command-line tool and API functions for creating term lists for the Gene Ontology, KEGG and InterPro.

# Installation
The Python version is currently in develop mode so you must manually clone and install it. Soon it will be available via pip.
```
git clone git@github.com:SystemsGenetics/FUNC-E.git
git checkout develop

cd FUNC_E
pip install .
```
# Preparing Input Files

Before using KINC you must prepare your files.  You will need to prepare four files:
- A genomic background file containing the full list of genes.
- A query list of genes that will be analyzed for enrichment of function.
- One or more files containing a list of functional terms used for enrichment.
- One or more files that associates functional terms with genes.

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
The query list contains the gene list that will undergo functional enrichment.  This file allows you to specify multiple groups (i.e. modules) of genes for enrichment. The file should be tab-delimited.  The first column should contain the list of genes and the second column the group (or module) name name.  The group name allows for multiple groups of genes to be listed in the same file but enrichment performed separately for each. The second column, however, may be left blank and only a single column of gene names can be provided.  

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

# Usage
## Command-line
### Generate Terms files
Bioinformatics tools such as [InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/), [Blast2GO](https://www.blast2go.com/) and [EnTAP](https://entap.readthedocs.io/en/latest/) (to name a few) provide the mapping of genes to controlled vocabulary terms, but creating the list of all terms in a vocabulary is still needed prior to enrichment. FUNC-E makes it easy to generate these for common vocabularies such as the Gene Ontology (GO), KEGG (KEGG) and InterPro (IPR).  

To generate a file of vocabulary terms from GO, KEGG and IPR use the following command:
```
FUNC-E-terms  --vocab GO KEGG IPR
```
This will create a file named `terms.tsv` ready for the format required by the `--terms` argument of `FUNC-E`

Alternatively, you can create separate files for each vocabulary:

```bash
FUNC-E-terms --outprefix KEGG --vocab KEGG
FUNC-E-terms --outprefix GO --vocab GO
FUNC-E-terms --outprefix IPR --vocab IPR
```

### Perform Functional Enrichment Analysis
FUNC-E provides the following usage instructions:

```
 FUNC-E [-h] --background BACKGROUND
     --query_list QUERY_LIST
     --ecut ECUT
     --terms TERMS [TERMS ...]
     --terms2features TERMS2FEATURES [TERMS2FEATURES ...]
     [--outprefix OUTPREFIX]
     [--module MODULE]
     [--vocab VOCAB [VOCAB ...]]
     [--similarity_threshold SIMILARITY_THRESHOLD]
     [--similarity_term_overlap SIMILARITY_TERM_OVERLAP]
     [--percent_similarity PERCENT_SIMILARITY]
     [--initial_group_membership INITIAL_GROUP_MEMBERSHIP]
     [--multiple_linkage_threshold MULTIPLE_LINKAGE_THRESHOLD]
     [--final_group_membership FINAL_GROUP_MEMBERSHIP]
     [--verbose VERBOSE]
```

For more detailed information about each argument please run the `FUNC-E -h` command.


#### Example
The following example performs functional enrichment of a list of arabidopsis genes from the TAIR10 genome assembly and annotation.  It requires that four types of files using the arguments: `--background`, `--query_list`, `--terms`, and `--terms2features`.  Note that the `--terms` and `--terms2features` arguments can be provided as many times as there are files.  In this example, lists of terms (provided with the `--terms` option) from  AraCyc, GO (Gene Ontology), IPR (InterPro), Pfam and PO (Plant Ontology) have been prepared, as well as the mappings of these term lists to the genes.  Genes to term mappings are provided using the `--terms2features` option.  

```bash
FUNC-E \
  --background arabidopsis_thaliana.TAIR10.genes.txt \
  --query_list modules.txt \
  --outprefix modules-enrichment \
  --terms AraCyc.terms.txt GO.terms.txt IPR.terms.txt \
          Pfam.terms.txt PO.terms.txt \
  --terms2features arabidopsis_thaliana.TAIR10.genes2AraCyc.txt \
                   arabidopsis_thaliana.TAIR10.genes2GO.txt \
                   arabidopsis_thaliana.TAIR10.genes2IPR.txt \
                   arabidopsis_thaliana.TAIR10.genes2Pfam.txt \
                   arabidopsis_thaliana.TAIR10.genes2PO.txt  \
  --ecut 0.01  
```
Additionally, the `--ecut` option provides a p-value cutoff for enrichment, and the --outprefix provides the a prefix which is added to every output file created by this script.



## Using the API
### Generate Terms files
Bioinformatics tools such as [InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/), [Blast2GO](https://www.blast2go.com/) and [EnTAP](https://entap.readthedocs.io/en/latest/) (to name a few) provide the mapping of genes to controlled vocabulary terms, but creating the list of all terms in a vocabulary is still needed prior to enrichment. FUNC-E makes it easy to generate these for common vocabularies such as the Gene Ontology (GO), KEGG (KEGG) and InterPro (IPR).  

To use the FUNC-E API to build a list of vocabularies, you must first import the package into your code:
```Python
from func_e.FUNC_E import FUNC_E
import func_e.vocabs.all as vocabs
```

To generate a Pandas DataFrame of vocabulary terms from GO, KEGG and IPR use the following function call:

```Python
terms = vocabs.getTerms(['GO', 'KEGG', 'IPR'])
```

### Perform Functional Enrichment Analysis
To perform functional enrichment using the FUNC-E API start by first instantiating a `FUNC_E` object.
```Python
fe = FUNC_E()
```

Next, you need to set the p-value cutoff for enrichment testing:
```Python
fe.setEnrichmentSettings({
    'ecut': 0.01
})
```

If you desire, you can change the clustering default settings as well:
```Python
fe.setClusteringSettings({
    'similarity_term_overlap': 3,
    'percent_similarity': 0.50,
    'initial_group_membership': 3,
    'multiple_linkage_threshold': 0.50,
    'final_group_membership':  3,
    'similarity_threshold': 0.5
})
```
The settings have the following meaning:

- `similarity_threshold`: This value is used to threshold the kappa scores. Pair-wise kappa scores are calculated for all genes. Kappa scores range between -1 to 1 and provide a measurement as to the similarity of annotations between two genes. Kappa scores greater than this value are considered meaningful and only those gene pairs with scores greater than this threshold are clustered. The default value if not specified is 0.35.
- `similarity_term_overlap`: Before kappa statistics are calculated two genes must share a specified number of terms. This parameter sets that minimum value. The default is 4.
- `percent_similarity`: Before clustering, seed groups are created, and when creating seed groups we want high quality groups. Therefore, the members of the seed groups must themselves share similarity with all other genes in the group greater or equal than the value specified by this parameter. The default is 0.50 (50 percent)
- `initial_group_membership`: When clustering, initial seed groups are created by grouping a gene with all other genes with which it has a significant (> similarity_threshold) kappa score. This parameter sets the minimum number of genes that must exist for a group to be considered a seed group. The default value is 4.
- `multiple_linkage_threshold`: After initial seed groups are formed an iterative process attempts to merge seed groups that have a specified percentage of genes in common. This parameter sets this percentage. The default is 0.50 (or seed groups must share 50 percent of genes to be merged).
- `final_group_membership`: This parameter sets the minimum number of terms in a cluster after all clustering. If the cluster has fewer terms it is thrown out. The default value is 4.


Next, FUNC_E can import the files needed for enrichment analysis. These are the same as the example files used in the command-line example above.
```Python
fe.importFiles({
    'background': 'arabidopsis_thaliana.TAIR10.genes.txt',
    'query': 'modules.txt',
    'terms2features': ['arabidopsis_thaliana.TAIR10.genes2AraCyc.txt',
                       'arabidopsis_thaliana.TAIR10.genes2GO.txt',
                       'arabidopsis_thaliana.TAIR10.genes2IPR.txt',
                       'arabidopsis_thaliana.TAIR10.genes2Pfam.txt',
                       'arabidopsis_thaliana.TAIR10.genes2PO.txt']
    'terms': ['IPR.terms.tsv', 'GO.terms.tsv', 'KEGG>terms.tsv']
})
```
Alternatively, you may have created the terms DataFrame using the `vocabs.getTerms()` function described above.  If so, you can leave out the `terms` argument in the `importFiles()` function call above and set the terms manually:

```Python
fe.setTerms(terms)
```
Now that FUNC-E has all of the necessary files and settings, you can perform functional enrichment:

```Python
fe.run()
```

If you only wish to perform enrichment analysis and not clustering you can provide the `cluster=False` argument:

```Python
fe.run(cluster=False)
```

If you want to limit enrichment to only a subset of modules and/or vocabularies you can provide the `modules` and `vocabs` arguments:

```Python
fe.run(modules=['module1', 'module2'], vocabs=['GO'])
```

Once completed you can access results using the following attributes of the `FUNC_E` object:

- `fe.enrichment`: a Pandas DataFrame containing the results of the enrichment test, including the p-value and bonferroni and bejamini corrected p-values
- `fe.clusters`: a Pandas DataFrame listing the the clusters that were identified contiaing the EASE score and geometric mean of p-values.
- `fe.cluster_terms`: a copy of the enrichment report, but with only clustered terms.

Finally, below are example commands to save results to a file:
```Python
fe.enrichment.sort_values('Fishers_pvalue').to_csv('FUNC-E.enriched_terms.tsv', sep="\t", index=None)

fe.clusters.to_csv('FUNC-E.clusters.tsv', sep="\t", index=None)

fe.cluster_terms.to_csv('FUNC-E.cluster_terms.tsv', sep="\t", index=None)
```
