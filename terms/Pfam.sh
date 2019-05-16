wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/pfamA.txt.gz
gunzip pfamA.txt.gz
cat pfamA.txt | awk -F"\t" '{print "Pfam\t"$2"\t"$5}' > Pfam.terms.txt
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/dead_family.txt.gz
gunzip dead_family.txt.gz
cat dead_family.txt | awk -F"\t" '{print "Pfam\t"$1"\t"$2" (deprecated): "$3}' > Pfam.deprecated-terms.txt
