# Retrieve the AraCyc Pathways for TAIR and convert to FUNC-E compatible terms file
wget https://www.arabidopsis.org/download_files/Pathways/Data_dumps/PMN9_September2014/aracyc_pathways.20140902
cat aracyc_pathways.20140902 | awk -F"\t" '{print "AraCyc\t"$1"\t"$2}' | sort -u > AraCyc.terms.txt

