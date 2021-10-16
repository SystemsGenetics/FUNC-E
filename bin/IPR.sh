wget ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list
cat entry.list | awk -F"\t" '{print $1"\t"$2"\t"$4}' | grep IPR > IPR.terms.txt

