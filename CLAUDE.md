# Landscaper-experiments

Landscaper（エネルギーランドスケープ解析）を用いた単一細胞トランスクリプトームのワークフロー集。

## プロジェクト構成

```
workflow/   Snakemakeワークフロー (.smk)
src/        R/shellスクリプト
data/       入力データ (Seuratオブジェクト, GO gene sets等)
output/     中間・最終出力
plot/       可視化結果
logs/       実行ログ
benchmarks/ 実行時間・メモリ記録
```

## ワークフロー一覧と依存関係

```
preprocess.smk          PCA・二値化・サンプル層別化
    ↓
landscaper.smk          Landscaper実行・後処理・DLGs算出
    ↓
    ├→ rank_estimate.smk      PC選択 (AUC cross-validation)
    ├→ random_walk.smk        遷移確率行列 (Metropolis/Glauber)
    │   ├→ vector_field.smk   粗視化ベクトル場
    │   ├→ graph_embedding.smk グラフ埋め込み
    │   └→ discon_graph.smk   非連結グラフ・デンドログラム
    ├→ fate_probability.smk   Fate Probability (吸収確率)
    └→ gsea.smk               GSEA (Basin / DLGs)

oulhen.smk              Oulhenデータ専用パイプライン (上記の大部分を内包)
```

各ワークフローは独立して実行: `snakemake -s workflow/<name>.smk`

## 共通設定・パラメータ

```python
PCS = list(range(3, 13))  # PC3〜PC12
SAMPLES = ['cont', 'DAPT', 'integrated']               # 基本3条件
# covありバージョン: 'cont_cov', 'DAPT_cov', 'integrated_cov'
OULHEN_SAMPLES = ['cont', 'DAPT', 'integrated']         # Oulhen用
GOS = ['bp', 'mf', 'cc']                                # GO ontology types
KFOLD = 5; SEED = 1                                     # PC AUC rank用
```

## Dockerコンテナ

| 用途 | イメージ |
|------|---------|
| Seurat解析全般 | `docker://koki/urchin_workflow_seurat:20251014` |
| Landscaper実行 | `docker://ghcr.io/chiba-ai-med/landscaper:pr-26` |
| GSEA | `docker://koki/urchin_workflow_gsea:20260113` |

## ルール定義の規約

```snakemake
rule rule_name:
    input:
        'output/{sample}/{pc}/input_file.tsv'
    output:
        'output/{sample}/{pc}/output_file.tsv'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/rule_name_{sample}_{pc}.txt'
    log:
        'logs/rule_name_{sample}_{pc}.log'
    shell:
        'src/rule_name.sh {input} {output} >& {log}'
```

- `resources: mem_mb=10000000` (10GB) がデフォルト。軽量処理は `1000000` (1GB)
- ログは必ず `>& {log}` でstdout+stderrをリダイレクト
- benchmarkとlogのファイル名にはワイルドカードを含める

## スクリプト規約

### シェルラッパー (src/*.sh)

```bash
#!/bin/bash
#$ -l nc=4
#$ -p -50
#$ -r yes
#$ -q node.q

#SBATCH -n 4
#SBATCH --nice=50
#SBATCH --requeue
#SBATCH -p node03-06
SLURM_RESTART_COUNT=2

Rscript src/script_name.R $@
```

- SGE (`#$`) と SLURM (`#SBATCH`) の両方のヘッダを記載
- `$@` で全引数をRscriptに転送
- 引数がハードコードされる場合もある (例: preprocess_landscaper_Oulhen.sh)

### Rスクリプト (src/*.R)

```r
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
# ...処理...
save(result, file=outfile)  # or write.table()
```

- 引数は `commandArgs(trailingOnly=TRUE)` で受け取る
- 共通関数は `source("src/Functions.R")` でロード
- Seurat v5対応: `GetAssayData()` の前に `JoinLayers()` が必要な場合がある

```r
if (inherits(seurat.integrated[["RNA"]], "Assay5")) {
  seurat.integrated[["RNA"]] <- JoinLayers(seurat.integrated[["RNA"]])
}
expr <- GetAssayData(seurat.integrated, assay="RNA", layer="data")
```

## 主要データオブジェクト

### Landscaper入力

| ファイル | 内容 |
|---------|------|
| `seurat.tsv` | cell x PC 行列 (float, ヘッダなし) |
| `group.tsv` | cell x 1 細胞型ラベル (ヘッダなし) |
| `cov.tsv` | cell x 1 共変量 (発生時間等、ヘッダなし、Oulhenでは不使用) |
| `BIN_DATA.tsv` | cell x PC 二値行列 ({-1, +1}) |

### Landscaper出力

| ファイル | 内容 |
|---------|------|
| `Allstates.tsv` | state x PC 二値パターン ({-1, +1}) |
| `E.tsv` | state x 1 エネルギー値 |
| `Basin.tsv` | Basin (安定吸引子) のstate ID |
| `Coordinate.tsv` | state x 2 ランドスケープ埋め込み座標 |
| `major_group.tsv` | state x 1 階層グループ |
| `igraph.RData` | 状態遷移ネットワーク (igraphオブジェクト) |
| `EnergyBarrier.tsv` | 状態間エネルギー障壁 |
| `SubGraph.tsv` | state → Basin所属 (吸引盆地ID) |
| `Allstates_major_group.tsv` | 二値パターン + state_id + major_group |

### 統合Seuratオブジェクト

`seurat_landscaper.RData` — Landscaper結果を統合したSeuratオブジェクト

- メタデータ: `state_idx`, `energy`, `is_basin`, `state_id`, `major_group`, `cell_time`, `celltype`, `celltype_colors`
- `@misc$landscaper`: Allstates, E, Basin, Coordinate, igraph等を格納したリスト
- 下流解析 (GSEA, DLGs, Fate Probability等) の入力

### 下流解析出力

| ファイル | 内容 |
|---------|------|
| `P_metropolis.tsv` / `P_glauber.tsv` | state x state 遷移確率行列 |
| `P_metropolis_fp.RData` / `P_glauber_fp.RData` | Fate Probability結果 |
| `dlgs.tsv` | DLGs (gene, DLG_score) |
| `*_gsea_*.RData` | fgsea結果 (res_list, gos) |
| `FINISH` | 複数出力ルールの完了マーカー |

## 出力ディレクトリ構造

```
output/{sample}/{pc}/                    中間出力
output/{comparison}/dlgs/{pc}/           DLGs (例: cont_DAPT)
plot/{sample}/{pc}/                      可視化
plot/{sample}/{pc}/Landscaper/           Landscaper直接出力
plot/{comparison}/dlgs/{pc}/             DLGs GSEA可視化
```

Oulhenデータの場合は `output/oulhen/`, `plot/oulhen/` 以下。

## 新しいワークフロー追加時の注意

1. **smkファイル冒頭**: `import re`, `min_version("6.5.3")`, PCS/SAMPLES定義, `container:` 指定
2. **wildcard_constraints**: `sample` と `pc` にはregex制約をつける
3. **Landscaper出力パス**: HPBaseデータは `plot/{sample}/{pc}/Landscaper/`、Oulhenは `plot/oulhen/{sample}/{pc}/Landscaper/`
4. **integrated先行**: cont/DAPTのLandscaperは `integrated` の `Coordinate.tsv` を共有座標として使用
5. **cov有無**: Oulhenデータにはcov.tsvがないため、postprocessスクリプトが異なる (`postprocess_landscaper_Oulhen.R`)
6. **FINISHマーカー**: 複数ファイルを出力するルールではFINISHファイルを使って完了を示す
