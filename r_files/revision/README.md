# Resubmission analysis layout

The reviewer-driven analyses are organized into four reproducible modules:

1. `revision_main/`: coordinate normalization and the revised pooled loop resource.
2. `ATAC_validation/`: matched-null validation of TSS-excluded ATAC overlap.
3. `downsampling/`: full-depth, 140M-contact, and 250M-contact sensitivity analyses.
4. `hrdp_genotype_diversity/`: SNP-based genetic diversity analysis of the ten strains.

Run `revision_main` first because its coordinate cache and output tables are
inputs to the ATAC and downsampling modules. The modules share
`r_files/funcs.R`, which is synchronized to the Google Drive `r_files` root.

Each module writes rebuildable files to its own `results/` directory. Large or
source-controlled inputs are stored under that module's `inputs/` directory in
the Google Drive copy and are not deleted by the post-commit synchronization.

## Software requirements

Use a recent R installation with `tidyverse`, `data.table`, `fs`, `digest`,
`scales`, `ggrepel`, `patchwork`, `GenomicRanges`, `GenomeInfoDb`, `IRanges`,
`S4Vectors`, `rtracklayer`, `AnnotationDbi`, `clusterProfiler`, and
`org.Rn.eg.db`. The genotype module additionally uses PLINK 2.0 and downloads
the official Apple Silicon binary when needed; on other platforms, set
`PLINK2_BIN` to a local PLINK 2.0 executable.

The Google Drive copy preserves the same `r_files/revision/` structure and
includes the shared `r_files/funcs.R`. Collaborators should download the four
module directories and `funcs.R` together so the scripts can resolve all
relative paths without using the original author's home directory.

---

## 3D Chromatin Loop vs Transcript Positional Relationships (5 Cases)

Evaluating whether gene transcripts reside inside, span across, or extend beyond chromatin loop boundaries (`loop_span: [x1, y2]`, `inter_anchor: (x2, y1)`).

```
[ 루프와 전사체의 위치 관계 5대 케이스 ]

┌─ 1. TSS가 앵커 [내부]에 위치할 때 (TSS inside anchor: x1~x2 or y1~y2)
│    ├── 【Case 1】 전사체가 루프 내부에 완전히 들어감 (Fully within loop span)
│    │              👉 22,390건 (39.30%)
│    └── 【Case 2】 전사체가 루프 경계를 넘어 바깥으로 삐져나감 (Spans loop boundary)
│                   👉 17,001건 (29.86%)
│
└─ 2. TSS가 앵커 [외부]에 위치할 때 (TSS outside anchor, TSS ±1kb 윈도우만 앵커에 닿음)
     ├─ 2-A. TSS가 루프 [바깥쪽] (x1 왼쪽 또는 y2 오른쪽)에 있을 때
     │    ├── 【Case 3】 전사체 몸통 전체가 루프 바깥에 머무름 (Outside loop span)
     │    │              👉 4,509건 (7.92%)
     │    └── 【Case 4】 전사체가 앵커를 뚫고 루프 안으로 진입하여 걸침 (Spans boundary into loop)
     │                   👉 3,714건 (6.52%)
     │
     └─ 2-B. TSS가 루프 [안쪽 공간] (x2와 y1 사이의 inter-anchor 영역)에 있을 때
          └── 【Case 5】 TSS는 안쪽 공간에 있고 전사체도 루프 내부에 머무름 (Within/spans loop)
                         👉 9,330건 (16.38%)
```

### Summary Statistics

* **전사체(Isoform) 레벨 (Ensembl 전사체 기준 총 56,944건)**:
  * **Case 1 (루프 완전 내부 포함)**: 22,390건 (39.30%)
  * **Case 2 (앵커 TSS $\rightarrow$ 루프 경계 돌파/걸침)**: 17,001건 (29.86%)
  * **Case 3 (외부 TSS $\rightarrow$ 전사체 몸통 전체 루프 바깥)**: 4,509건 (7.92%)
  * **Case 4 (외부 TSS $\rightarrow$ 전사체가 루프 안으로 진입하여 걸침)**: 3,714건 (6.52%)
  * **Case 5 (Inter-anchor TSS $\rightarrow$ 루프 내부 머무름)**: 9,330건 (16.38%)

* **유전자(Gene) 레벨 요약 (`df.direct.gene.assignment.position.flags`, 총 26,920건)**:
  * `all_transcripts_fully_within_loop_span` (모든 전사체가 루프 내부 완벽 포함): **14,101건 (52.38%)**
  * `overlap_without_full_transcript_containment` (일부/전체 전사체가 루프 경계에 걸침): **9,602건 (35.67%)**
  * `some_transcripts_fully_within_loop_span` (일부 전사체는 내부 포함, 일부는 걸침): **1,460건 (5.42%)**
  * `all_annotated_transcripts_outside_loop_span` (모든 전사체가 루프 바깥에 위치): **856건 (3.18%)**
  * `no_Ensembl_transcript_annotation` (Ensembl 전사체 정보 없는 EPD 전용 유전자): **901건 (3.35%)**

---

## Predicted CTCF Motif Annotations (`df.ctcf.evidence`)

CTCF sequence predictions are retained as an independent structural annotation. Motif counts do not select loops, define regulatory categories, determine direction, or imply in-vivo CTCF occupancy.

### CTCF Motif Dataset Summary

1. **FIMO 원본 예측 (Raw Predictions)**:
   * **`6,551,641개`** (100개의 다양한 CTCF 모티프 서열 모델로 쥐 게놈 rn7 전체를 스캔한 원본 결과)

2. **좌표 단위 비중복 정제 모티프 (`gr.ctcf.motif`)**:
   * **⭐️ `3,191,859개` (약 319만 개)**
   * 동일한 게놈 좌표에 여러 모티프 예측이 겹칠 경우, 가장 통계적 유의성이 높은(Minimum $p$-value) 최적의 대표 모티프를 선별하여 **3,191,859개의 고유 모티프 구간(GRanges)**으로 구성하여 사용.

3. **31,021개 루프 앵커에 오버랩된 모티프 총 수**:
   * **`2,775,067건`** (루프 앵커당 평균 수십 개씩 매핑됨)

### 31,021 Pooled Loops CTCF Classification (`predicted_ctcf_motif_annotation_class`)

| 분류 범주 | 루프 수 (`n`) | 비율 (`%`) | 구조적 의미 |
| :--- | :---: | :---: | :--- |
| **`predicted_motif_intervals_at_both_anchors`** | **29,980** | **96.6%** | 양쪽 앵커 모두 CTCF 모티프 보유 (전형적인 CTCF-매개 루프 구조) |
| **`predicted_motif_intervals_at_one_anchor`** | **929** | **3.0%** | 한쪽 앵커에만 CTCF 모티프 보유 |
| **`no_predicted_motif_interval_overlap`** | **112** | **0.4%** | 양쪽 앵커 모두 CTCF 모티프 없음 |

---
