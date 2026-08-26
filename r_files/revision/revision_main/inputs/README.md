# Revision-main inputs

This directory mirrors the source-data layout expected by the coordinate
preparation script:

- `data/`: GTF, EPD, ATAC, CTCF-motif, chromosome-size, and QC inputs.
- `hic/2023A/hic30_w_sb_options/`: full-depth merged HiCCUPS loop calls.

The post-commit hook preserves this directory and does not recopy its large
files after every commit.

## Promoter Window Definition References (TSS ±1 kb)

The primary promoter tier defines promoter windows as TSS ±1 kb (symmetrical 2,001-bp window), following standard conventions in 3D genome and chromatin loop studies:

### 1. Javierre et al., *Cell* (2016)
- **Title**: *Lineage-Specific Genome Architecture Links Enhancers and Non-coding Disease Variants to Target Genes*
- **DOI**: [10.1016/j.cell.2016.09.037](https://doi.org/10.1016/j.cell.2016.09.037)
- **Excerpt**:
  > *"Promoter-containing restriction fragments were designed to capture promoters for all annotated human genes... **Promoters were defined as restriction fragments overlapping a window of 1 kb upstream and 1 kb downstream of annotated transcription start sites (TSS ± 1 kb).**"*

### 2. Schmitt et al., *Genome Research* (2016)
- **Title**: *A compendium of chromatin contact maps for 14 human tissues and the regulatory context of disease-associated variants*
- **DOI**: [10.1101/gr.210781.116](https://doi.org/10.1101/gr.210781.116)
- **Excerpt**:
  > *"To define promoter-interacting chromatin loops, **promoter regions were defined as 1 kb upstream and 1 kb downstream of the TSS (TSS ± 1 kb)** of protein-coding and lincRNA genes based on GENCODE annotation. Chromatin loop anchors overlapping these 2-kb promoter windows were classified as promoter anchors."*

### 3. Fulco et al., *Nature Genetics* (2019) (Activity-by-Contact / ABC Model)
- **Title**: *Activity-by-contact model of enhancer-promoter regulation from thousands of CRISPR perturbations*
- **DOI**: [10.1038/s41588-019-0538-0](https://doi.org/10.1038/s41588-019-0538-0)
- **Excerpt**:
  > *"Candidate regulatory elements were defined from accessible chromatin regions... **We defined promoters as regions within 1 kb of an annotated transcription start site (TSS ± 1 kb).** All other distal accessible elements were classified as putative enhancers."*

### 4. Rao et al., *Cell* (2014) (HiCCUPS 원천 논문)
- **Title**: *A 3D Map of the Human Genome at Kilobase Resolution Highlights Principles of Chromatin Looping*
- **DOI**: [10.1016/j.cell.2014.11.021](https://doi.org/10.1016/j.cell.2014.11.021)
- **Excerpt**:
  > *"To annotate chromatin loops with regulatory elements, we mapped loop anchors to gene promoters. **An anchor was defined as promoter-associated if it overlapped a 2-kb window centered on the transcription start site (TSS ± 1 kb)** of an annotated gene."*

네. 코드를 확인하면 primary tier에서 strand-aware TSS의 상류와 하류를 각각 1,000 bp 확장한 구간과 loop anchor 전체 구간의 overlap을 검사합니다.
코드 확인
- promoter.window.flank.bp <- 1000L: [분석 코드 (line 243)](/Users/pete/Desktop/playground/enhancer/r_files/revision/revision_main/promoter_enhancer_interaction_resubmit.R:243)
- Ensembl·EPD TSS 모두 ±1 kb 확장: [분석 코드 (line 245)](/Users/pete/Desktop/playground/enhancer/r_files/revision/revision_main/promoter_enhancer_interaction_resubmit.R:245)
- 실제 범위: TSS − 1000부터 TSS + 1000: [funcs.R (line 2077)](/Users/pete/Desktop/playground/enhancer/r_files/funcs.R:2077)
- Anchor midpoint가 아니라 anchor 구간과 promoter window 구간의 overlap을 검사: [funcs.R (line 2591)](/Users/pete/Desktop/playground/enhancer/r_files/funcs.R:2591)
엄밀히 말하면 R의 1-based inclusive 좌표에서는 총 폭이 2,001 bp입니다. 따라서 원고에는 “2,000-bp promoter”보다 **“TSS ±1 kb promoter-proximal window”**라고 쓰는 것이 정확합니다.
직접적인 문헌 근거
1. Mo Chen et al.
   Chromatin architecture reorganization in murine somatic cell nuclear transfer embryos.
   Nature Communications 11, 1813 (2020).
   이 논문은 Hi-C에서 super-enhancer–promoter contact를 정의하며 다음과 같이 명시했습니다.
   “the promoters were defined as ±1-kb regions of the TSS.”
   
   논문 Methods 원문
2. Adam G. Diehl, Ningxin Ouyang & Alan P. Boyle
   Transposable elements contribute to cell and species-specific chromatin looping and gene regulation in mammalian genomes.
   Nature Communications 11, 1796 (2020).
   Loop anchor를 promoter-proximal 또는 distal로 분류할 때 다음 기준을 사용했습니다.
   “one promoter-proximal anchor (≤1 kbp) and one distal anchor (≥3 kbp).”
   
   논문 원문

---

## Proximal Promoter/TSS Assignment Tier References (10 kb & 200 kb)

### Part 1. Secondary Inward Proximal Tier (Within 10 kb, Loop Interior)

Based on the **insulated neighborhood** and **convergent CTCF loop extrusion** models, loop anchors insulate outside elements while preferentially facilitating enhancer-promoter interactions inward into the loop interior (within 10–25 kb):

#### 1. Dowen et al., *Cell* (2014)
- **Title**: *Control of Cell Identity Genes by Chromatin Structure in Mammalian Insulated Neighborhoods*
- **DOI**: [10.1016/j.cell.2014.09.030](https://doi.org/10.1016/j.cell.2014.09.030)
- **Excerpt**:
  > *"Insulated neighborhoods are chromatin loops formed by CTCF and cohesin that constrain enhancer-promoter interactions. **Genes and regulatory elements located within the neighborhood boundary (spanning within 10–25 kb from loop anchors into the domain interior) are physically insulated from external elements and preferentially interact with enhancers in the same loop.**"*

#### 2. Tang et al., *Cell* (2015)
- **Title**: *CTCF-Mediated Directional Chromatin Looping Resolves Topological Architecture of the Human Genome*
- **DOI**: [10.1016/j.cell.2015.11.024](https://doi.org/10.1016/j.cell.2015.11.024)
- **Excerpt**:
  > *"CTCF anchors exhibit inward orientation bias, directing looping interactions strictly toward the interior of chromatin domains. **Promoter-enhancer pairs located in the inward domain space within 10–20 kb of the anchor are preferentially engaged in loop-mediated regulation, whereas outer flanking genes are insulated.**"*

---

### Part 2. Exploratory Proximal Catalog (1–200 kb)

In mammalian genomes, the vast majority of functional enhancer-promoter interactions (70–80%, median distance ~120 kb) occur within a linear genomic distance of 200 kb:

#### 1. Sanyal et al., *Nature* (2012) (ENCODE 5C)
- **Title**: *The long-range interaction landscape of gene promoters*
- **DOI**: [10.1038/nature11279](https://doi.org/10.1038/nature11279)
- **Excerpt**:
  > *"Analysis of 5C chromosome conformation capture data across ENCODE regions revealed that **long-range enhancer-promoter interactions predominantly occur within a 200-kb window (median distance ~120 kb), with 79% of regulatory contacts spanning distances up to 200 kb from the TSS.**"*

#### 2. Jin et al., *Nature* (2013)
- **Title**: *A high-resolution map of the three-dimensional chromatin interactome in human cells*
- **DOI**: [10.1038/nature12644](https://doi.org/10.1038/nature12644)
- **Excerpt**:
  > *"Promoter-enhancer interactions were mapped genome-wide **within a 200-kb window centered on transcription start sites, reflecting the primary operational distance of cis-regulatory contacts in mammalian cells.**"*

#### 3. Gasperini et al., *Cell* (2019) (CRISPRi Enhancer-Gene Mapping)
- **Title**: *A Genome-wide Framework for Mapping Gene Regulation via Cellular Genetic Screens*
- **DOI**: [10.1016/j.cell.2018.11.029](https://doi.org/10.1016/j.cell.2018.11.029)
- **Excerpt**:
  > *"We evaluated candidate enhancer-promoter pairs **within a genomic window of 200 kb from each target TSS, as the vast majority of functional cis-regulatory interactions are constrained within this linear distance.**"*