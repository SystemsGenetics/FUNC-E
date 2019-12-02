perl ../FUNC-E.pl \
  --query_list demo_query.txt \
  --background  oryza_sativa.MSU_v7_0.genes.txt \
  --terms GO.terms.txt \
  --terms IPR.terms.txt \
  --terms2features oryza_sativa.MSU_v7_0.genes2GO.txt \
  --terms2features oryza_sativa.MSU_v7_0.genes2IPR.txt \
  --outprefix demo_run_pl \
  --ecut 0.01 \
  --preset medium
