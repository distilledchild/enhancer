# King 2025 Repeated H-MAGMA-only Gene Annotation

## Repeated Genes

   gene_name     hmagma_gene_id   CHR     START      STOP n_traits
      <char>             <char> <int>     <int>     <int>    <int>
1:    Klhl35 ENSRNOG00000017365     1 153767970 153773557        3
2:    Lrtomt ENSRNOG00000023434     1 156266655 156268704        3
3:    Rnf169 ENSRNOG00000026408     1 154254473 154311653        3
4:   Lamtor1 ENSRNOG00000020016     1 156272064 156291179        3
                                                                          traits
                                                                          <char>
1: crf_ny_incentive_value_index; crf_ny_lever_presses; pavca_ny_d5_response_bias
2: crf_ny_incentive_value_index; crf_ny_lever_presses; pavca_ny_d5_response_bias
3: crf_ny_incentive_value_index; crf_ny_lever_presses; pavca_ny_d5_response_bias
4: crf_ny_incentive_value_index; crf_ny_lever_presses; pavca_ny_d5_response_bias
   max_hmagma_Z min_hmagma_P min_hmagma_BH_FDR hmagma_NSNP_range
          <num>        <num>             <num>            <char>
1:       5.4582   2.4046e-08      5.636039e-05             22-22
2:       5.1112   1.6002e-07      1.762317e-04               8-8
3:       5.1097   1.6138e-07      1.762317e-04               9-9
4:       4.5858   2.2615e-06      1.236814e-03               5-5
   cmagma_NSNP_observed n_snp_cmagma_annotation
                 <char>                   <int>
1:                                            2
2:                                            2
3:                                            9
4:                                            1
   n_snp_hmagma_filtered_annotation added_snps_filtered_annotation
                              <int>                          <int>
1:                               39                             37
2:                               10                              8
3:                               20                             11
4:                                6                              5
   noise_snps_removed_annotation biological_priority
                           <int>              <char>
1:                            96         Exploratory
2:                            72         Exploratory
3:                            99         Exploratory
4:                            62            Moderate
                                                                                                                                                                              core_biology
                                                                                                                                                                                    <char>
1:                                                                                        Kelch-like protein; predicted CUL3-RING ubiquitin ligase substrate adaptor/proteostasis biology.
2: Readthrough/fusion-like leucine-rich transmembrane/O-methyltransferase locus; TOMT ortholog biology is linked to mechanotransduction and catecholamine methyltransferase-like function.
3:                                                                      Ubiquitin/chromatin DNA-damage-response regulator; limits/competes in RNF168/53BP1-type repair-factor recruitment.
4:                                                                        Lysosomal Ragulator anchor/adaptor upstream of mTORC1; links lysosomal nutrient signaling to MAPK/MTOR activity.
                                                                                                                                                                                    neuro_addiction_interpretation
                                                                                                                                                                                                            <char>
1:                                                                   No strong direct addiction evidence found yet; repeated H-MAGMA signal makes it a regulatory-contact candidate, not a mature biological lead.
2:                   Indirect candidate only. The catecholamine-methyltransferase-like naming is tempting, but published functional work is mainly auditory hair-cell mechanotransduction/deafness, not addiction.
3:                                                             No direct reward/addiction evidence found yet; could be relevant through chromatin/ubiquitin stress-response biology, but this remains speculative.
4: Best biological lead among the four because mTORC1 signaling is repeatedly implicated in neuronal plasticity and drug reward/cue learning; the gene-level claim still needs expression and locus-level support.
                                                              key_references
                                                                      <char>
1: NCBI Gene 283212; Human Protein Atlas KLHL35 brain single-cell expression
2:                                              PMID:18953341; PMID:28504928
3:                                           PMID:22733822; NCBI Gene 254225
4:                               PMID:31001086; PMID:29039413; PMID:24666346

## Interpretation

- `Lamtor1` is the strongest biological lead among the repeated genes because it sits in the lysosomal Ragulator/mTORC1 axis, which has published links to neuronal plasticity and drug reward/cue learning.
- `Klhl35`, `Lrtomt`, and `Rnf169` should be treated as exploratory regulatory-contact candidates until brain expression, pathway context, and independent phenotype/dataset recurrence are confirmed.
- All four genes are repeated H-MAGMA-only hits across the three selected King phenotypes at BH-FDR < 0.05; several also pass Bonferroni in CRF lever pressing or PavCA response bias.
- The local cMAGMA annotation contains a small number of SNPs for these genes, but those cMAGMA-local SNPs have no matching King GWAS p-values in the rn7 p-value files and no cMAGMA `genes.out` rows were produced. The observed gene-level signal therefore comes from strict H-MAGMA added SNPs.

## Output Files

- `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/king_1_3_repeated_genes_biology_annotation.tsv`
- `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/king_1_3_repeated_genes_trait_stats.tsv`
- `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/king_1_3_repeated_genes_snp_assignments.tsv`
- `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/king_1_3_repeated_genes_locus_figure_index.tsv`
- Figures: `/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/reports/repeated_gene_locus_figures`
