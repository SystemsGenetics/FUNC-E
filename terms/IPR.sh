wget ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list
perl -pi -e 's/^(IPR\d+)\s/IPR\t\1\t/' entry.list
cat entry.list | grep IPR | awk -F"\t" '{print $1"\t"$2"\t"$4}' > IPR.terms.txt

