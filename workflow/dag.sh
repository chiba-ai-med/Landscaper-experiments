# DAG graph
snakemake -s workflow/preprocess.smk --rulegraph | dot -Tpng > plot/preprocess.png
snakemake -s workflow/landscaper.smk --rulegraph | dot -Tpng > plot/landscaper.png
snakemake -s workflow/rank_estimate.smk --rulegraph | dot -Tpng > plot/rank_estimate.png
snakemake -s workflow/random_walk.smk --rulegraph | dot -Tpng > plot/random_walk.png
snakemake -s workflow/vector_field.smk --rulegraph | dot -Tpng > plot/vector_field.png
snakemake -s workflow/fate_probability.smk --rulegraph | dot -Tpng > plot/fate_probability.png
snakemake -s workflow/graph_embedding.smk --rulegraph | dot -Tpng > plot/graph_embedding.png
snakemake -s workflow/discon_graph.smk --rulegraph | dot -Tpng > plot/discon_graph.png
snakemake -s workflow/gsea.smk --rulegraph | dot -Tpng > plot/gsea.png

# Oulhen's data
snakemake -s workflow/oulhen.smk --rulegraph | dot -Tpng > plot/oulhen.png
