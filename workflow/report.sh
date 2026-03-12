# HTML
mkdir -p report
snakemake -s workflow/preprocess.smk --report report/preprocess.html
snakemake -s workflow/landscaper.smk --report report/landscaper.html
snakemake -s workflow/rank_estimate.smk --report report/rank_estimate.html
snakemake -s workflow/random_walk.smk --report report/random_walk.html
snakemake -s workflow/vector_field.smk --report report/vector_field.html
snakemake -s workflow/fate_probability.smk --report report/fate_probability.html
snakemake -s workflow/graph_embedding.smk --report report/graph_embedding.html
snakemake -s workflow/discon_graph.smk --report report/discon_graph.html
snakemake -s workflow/gsea.smk --report report/gsea.html

# Oulhen's data
snakemake -s workflow/oulhen.smk --report report/oulhen.html
