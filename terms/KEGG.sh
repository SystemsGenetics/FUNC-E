
# Create a file of the KEGG orthologs
wget http://rest.kegg.jp/list/ko 
cat ko | awk -F"\t" '{print "KEGG\t"$1"\t"$2}' | perl -p -e 's/ko://' > KEGG_orthologs.terms.txt

# Create a file of the KEGG pathways
wget http://rest.kegg.jp/list/pathway
cat pathway | awk -F"\t" '{print "KEGG\t"$1"\t"$2}' | perl -p -e 's/path:map/ko/' > KEGG_pathways.terms.txt

# Create a file of the KEGG modules
wget http://rest.kegg.jp/list/md
cat md | awk -F"\t" '{print "KEGG\t"$1"\t"$2}' | perl -p -e 's/md://' > KEGG_modules.terms.txt

cat KEGG_orthologs.terms.txt KEGG_pathways.terms.txt KEGG_modules.terms.txt > KEGG.terms.txt


