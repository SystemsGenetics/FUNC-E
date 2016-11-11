cat ../../organisms/oryza_sativa/RiceCyc_v3.3/ricecyc_v3_3_pwy.txt | awk -F"\t" '{print "RiceCyc\t"$6"\t"$7}' | sort -u  > RiceCyc.terms.txt
