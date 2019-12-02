python3 ../FUNC-E.py \
  --query_list demo_query.txt \
  --background  oryza_sativa.MSU_v7_0.genes.txt \
  --terms GO.terms.txt IPR.terms.txt \
  --terms2features oryza_sativa.MSU_v7_0.genes2GO.txt oryza_sativa.MSU_v7_0.genes2IPR.txt \
  --outprefix demo_run_py \
  --ecut 0.01 \
  --preset medium \
  -v
