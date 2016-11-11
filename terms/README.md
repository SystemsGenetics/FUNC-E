The shell scripts in this directory are meant to assist in creation of the terms list (needed for the --terms argument) for FUNC-E. Each script corresponds to a different vocabulary or database of controlled terms.  At the time that these scripts are added they are verified to work.  However, if the corresponding database changes file formats or URLs these scripts may break.

The following scripts are provided:

- AraCyc.sh:  used for creating a term list from TAIR's AraCyc database.
- GO.sh:  creates a term list from the Gene Ontology
- IPR.sh: creates a term list from InterPro
- KEGG.sh:  creates a term list of KEGG orthologs and pathways.
- PO.sh: creates a term list from the Plant Ontology
- Pfam.sh:  creates a list of terms from the Pfam database.
- RiceCyc.sh: creates a list of terms from the RiceCyc database

To use any of the terms from these databases your must have a mapping of the term names from the controlled vocaubulary/database mapped to the genes in your genomic background.  These mappings must be obtained from the site where your transcript or genome assembly is housed. If such mappings are not available it is possible to create them yourself using tools such as InterProScan, Blast2GO, and the KEGG Automatic Annotation Service (KEGG), for example.
