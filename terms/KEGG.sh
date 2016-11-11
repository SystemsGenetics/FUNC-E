lynx "http://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir="

# first get the list of pathways
grep -P "^C" ko00001.keg > pathways.txt
perl -pi -e 's/^C\s+(\d+)\s(.*?)\s*\[.*/KEGG\tko\1\t\2/' pathways.txt 

# second get the list of othologs
egrep "K0|K1" ko00001.keg > orthologs.txt
perl -pi -e 's/^D\s+(K\d+)\s+.*?;\s+(.*?)/KEGG\t\1\t\2/' orthologs.txt

cat pathways.txt orthologs.txt > KEGG.terms.txt

