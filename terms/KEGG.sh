wget "http://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir=" -O ko00001.keg

# first get the list of pathways
grep -P "^C" ko00001.keg | grep "\[PATH:" | \
perl -p -e 's/^C\s+(\d+)\s(.*?)\s*\[PATH:(.*)\].*/KEGG\t\3\t\2/' > pathways.txt 

# second get the list of othologs
grep -P "^D" ko00001.keg | \
perl -p -e 's/^D\s+(K\d+)\s+.*?;\s+(.*?)/KEGG\t\1\t\2/' > orthologs.txt

cat pathways.txt orthologs.txt > KEGG.terms.txt

