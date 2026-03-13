# Landscaper-experiments

Energy Landscape Analysis workflows for single-cell transcriptomics using [Landscaper](https://github.com/chiba-ai-med/Landscaper).

This workflow consists of 10 workflows as follows:

- **workflow/preprocess.smk**: PCA, binarization, and sample stratification

![](https://github.com/chiba-ai-med/Landscaper-experiments/blob/main/plot/preprocess.png?raw=true)

- **workflow/landscaper.smk**: Landscaper execution, postprocessing, and DLGs computation

![](https://github.com/chiba-ai-med/Landscaper-experiments/blob/main/plot/landscaper.png?raw=true)

- **workflow/rank_estimate.smk**: PC rank selection by AUC cross-validation

- **workflow/random_walk.smk**: Transition probability matrices (Metropolis/Glauber) and random walk visualization

![](https://github.com/chiba-ai-med/Landscaper-experiments/blob/main/plot/random_walk.png?raw=true)

- **workflow/vector_field.smk**: Coarse-grained vector fields

![](https://github.com/chiba-ai-med/Landscaper-experiments/blob/main/plot/vector_field.png?raw=true)

- **workflow/graph_embedding.smk**: Graph embedding of state transition networks

![](https://github.com/chiba-ai-med/Landscaper-experiments/blob/main/plot/graph_embedding.png?raw=true)

- **workflow/discon_graph.smk**: Disconnectivity graphs and dendrograms

![](https://github.com/chiba-ai-med/Landscaper-experiments/blob/main/plot/discon_graph.png?raw=true)

- **workflow/fate_probability.smk**: Fate probability (absorption probabilities)

![](https://github.com/chiba-ai-med/Landscaper-experiments/blob/main/plot/fate_probability.png?raw=true)

- **workflow/gsea.smk**: GSEA for Basin and DLGs (waterfall plot, dot plot)

![](https://github.com/chiba-ai-med/Landscaper-experiments/blob/main/plot/gsea.png?raw=true)

- **workflow/oulhen.smk**: Oulhen dataset pipeline (encompasses most of the above, plus cross-species DLGs overlap analysis)

![](https://github.com/chiba-ai-med/Landscaper-experiments/blob/main/plot/oulhen.png?raw=true)

## Requirements

- Bash: GNU bash, version 4.2.46(1)-release (x86_64-redhat-linux-gnu)
- Snakemake: 6.5.3
- Singularity: 3.8.0

## How to reproduce this workflow

### In Local Machine

```
snakemake -s workflow/preprocess.smk -j 4 --use-singularity
snakemake -s workflow/landscaper.smk -j 4 --use-singularity
snakemake -s workflow/rank_estimate.smk -j 4 --use-singularity
snakemake -s workflow/random_walk.smk -j 4 --use-singularity
snakemake -s workflow/vector_field.smk -j 4 --use-singularity
snakemake -s workflow/graph_embedding.smk -j 4 --use-singularity
snakemake -s workflow/discon_graph.smk -j 4 --use-singularity
snakemake -s workflow/fate_probability.smk -j 4 --use-singularity
snakemake -s workflow/gsea.smk -j 4 --use-singularity
snakemake -s workflow/oulhen.smk -j 4 --use-singularity
```

### In Open Grid Engine

```
snakemake -s workflow/preprocess.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/landscaper.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/rank_estimate.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/random_walk.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/vector_field.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/graph_embedding.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/discon_graph.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/fate_probability.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/gsea.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/oulhen.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
```

### In Slurm

```
snakemake -s workflow/preprocess.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/landscaper.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/rank_estimate.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/random_walk.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/vector_field.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/graph_embedding.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/discon_graph.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/fate_probability.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/gsea.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/oulhen.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
```

## Authors

- Koki Tsuyuzaki
