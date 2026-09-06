# King 2025 H-MAGMA 분석 심층 검토 및 비판적 리뷰 리포트

**문서 위치:** `/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/King2025/review.md`  
**작성 일자:** 2026-09-05  
**대상 코호트:** King et al. 2025 (*Heterogeneous Stock Rats Pavlovian Conditioned Approach & Conditioned Reinforcement GWAS*)  
**분석 파이프라인:** Palmer Lab GCTA MLMA 요약 통계량 기반 cMAGMA vs H-MAGMA (성체 랫트 PFC Hi-C + snATAC 조절 루프)

> **검증 response 추가: 2026-09-05.** 리뷰 원문은 보존했으며, 아래 Response가 코드·원자료·재계산에 근거한 판단이다. 특히 이 문서의 79/103/40/78개 결과는 **legacy H-MAGMA 대 cMAGMA** 비교이다. 새 13,376개 루프의 결과는 같은 순서로 **81/102/42/85개**이며, 현재 수행한 주 비교는 **새 H-MAGMA 대 legacy H-MAGMA**이다. 두 비교를 혼용하지 않는다.

### Response 0: 요약·수치표의 공통 정정

- 공개 파일 1,180개와 주요 3개 형질의 lead SNP/보고서 수치 일치는 확인 자료가 있다. 그러나 이는 모든 형질·모든 SNP 또는 논문 최종 분석의 **100% 재현**을 입증하지 않는다. 논문 Methods의 SNP 수는 3,400,759개이고, 사용한 공개 MLMA는 형질당 3,513,494행이므로 두 범위를 동일하다고 단정하지 않는다. 출처는 **공개 기탁/연결 보고서의 summary statistics를 이용한 재분석**으로 기술한다. [King 원문](https://onlinelibrary.wiley.com/doi/10.1111/gbb.70018)
- 표의 H-MAGMA-only 20/10/7/8개는 기존 좌표 매칭 기준 집합 차이이지 **순증 수**가 아니다. 표시된 유의 행 수의 산술 차이는 lever presses +19, response bias +11, incentive value +5, index +7이다. cMAGMA는 좌표 ID, H-MAGMA는 Ensembl gene ID여서 집계 단위도 명시해야 한다. Bonferroni의 H-MAGMA-only 수 역시 순증과 구분한다.
- cMAGMA는 **exon 및 promoter 기반**이며 coding-only/CDS-only 분석이 아니다. 반복 유전자 표의 `+37`, `+11` 등은 참조 주석 수준의 추가 SNP 수이고, 해당 King GWAS에서 실제 사용된 SNP 수와 다르다(Response 2).
- Hi-C는 우리 **pooled frontal-cortex** 자료이고, PFC snATAC는 별도의 외부 자료이다. 같은 동물에서 얻은 matched PFC Hi-C/snATAC 자료나 세포형 특이 네트워크로 서술하지 않는다.
- 검증 코드: [verify_review_claims.R](/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/King2025/verify_review_claims.R). [검증 산출물](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/review_audit), [현재 새 루프 결과](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/revision_2001bp_loop_only/RESULTS_REPORT.md).

---

## 1. 개요 및 핵심 결론 (Executive Summary)

1. **데이터셋 완벽 재현성 검증 상태 (Confirmed ★★★★★):**
   - UCSD Digital Collections 공식 기탁 데이터(`bb5030313v`)의 염색체별 MLMA 요약 통계량 1,180개 파일을 직접 연동.
   - 각 형질별 Lead SNP의 $-\log_{10}(P)$ 및 QTL 위치가 공식 보고서(`gwas_report.html`)와 **소수점 단위까지 100% 일치**.
   - rn6 $\rightarrow$ rn7 liftOver 성공률 99.8%, HS rat v4 LD 참조 패널과 96.5% 일치율 달성.
2. **신규 유전자 발굴 성과:**
   - 충분한 표본 크기($N \approx 1,600$)와 강력한 GWAS 신호($-\log_{10} P = 7.2 \sim 9.5$) 덕분에, 기존 코딩 변이 중심 분석(cMAGMA) 대비 **형질당 7 ~ 20개(FDR < 0.05 기준)의 신규 유전자가 H-MAGMA를 통해 대거 발굴**됨.
   - 특히 3개 형질 모두에서 공통 반복 발굴된 **`Lamtor1` (mTORC1 축), `Klhl35`, `Lrtomt`, `Rnf169`**는 매우 유력한 신경계 조절 타깃으로 부상함.
3. **심층 감사(Audit)의 핵심 경고:**
   - 발굴된 신규 유전자 20개 중 **18개(90%)가 염색체 1번 9.8 Mb 구간(146.6 ~ 156.5 Mb)의 단일 거대 QTL 블록에 집중**되어 있음.
   - 이를 "20개의 독립적인 신규 원인 유전자"로 과장해서는 안 되며, **"Chr 1 거대 행동 조절 QTL 좌위 내부에서 3D 루프를 통해 비코딩 변이들이 연결되는 핵심 후보 유전자들을 좁혀낸 정밀 우선순위화(Prioritization / Fine-mapping)"**로 기술해야 리뷰어의 공격을 완벽히 방어할 수 있음.

---

## 2. cMAGMA 대비 H-MAGMA 신규 발굴 결과 상세

### (1) 형질별 유의 유전자 수 비교 요약

| 형질 (Phenotype) | 기준 (Threshold) | cMAGMA (코딩/프로모터만) | strict H-MAGMA (루프 추가) | **H-MAGMA-only 신규 발굴 수** |
| :--- | :--- | :---: | :---: | :---: |
| **`crf_ny_lever_presses`**<br>(조건 강화 레버 누르기 횟수) | **BH-FDR < 0.05**<br>**Bonferroni < 0.05** | 60개<br>26개 | **79개**<br>**31개** | **+20개 순증**<br>**+8개 순증** |
| **`pavca_ny_d5_response_bias`**<br>(파블로프 Day 5 반응 편향) | **BH-FDR < 0.05**<br>**Bonferroni < 0.05** | 92개<br>32개 | **103개**<br>**35개** | **+10개 순증**<br>**+7개 순증** |
| **`crf_ny_incentive_value_index`**<br>(보상 유인 가치 지수) | **BH-FDR < 0.05**<br>**Bonferroni < 0.05** | 35개<br>6개 | **40개**<br>**7개** | **+7개 순증**<br>**+2개 순증** |
| **`pavca_ny_d5_index`**<br>(파블로프 Day 5 종합 지수) | **BH-FDR < 0.05**<br>**Bonferroni < 0.05** | 71개<br>14개 | **78개**<br>**15개** | **+8개 순증**<br>**+2개 순증** |

*※ FDR < 0.10 기준 완화 시 형질당 15 ~ 21개의 신규 유전자가 추가됨.*

---

### (2) 3개 형질 공통 반복 발굴 4대 핵심 유전자 (Repeated Hits)

King 2025의 3개 주요 형질(`crf_ny_lever_presses`, `crf_ny_incentive_value_index`, `pavca_ny_d5_response_bias`) 전반에서 공통적으로 검출된 4대 핵심 유전자:

| 유전자명 | 위치 (Chr 1) | cMAGMA $P$ | H-MAGMA 최저 $P$ (최저 FDR) | 추가 루프 SNP | 생물학적 기능 및 신경계 연관성 |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **`Lamtor1`** | 156.27 Mb | **검정 불가**<br>(코딩 SNP 결측) | **$2.26 \times 10^{-6}$**<br>($\text{FDR} = 1.24 \times 10^{-3}$) | **+5개** | **생물학적 1순위 리드**: 리소좀 Ragulator 복합체 앵커이자 **mTORC1 신호전달 상류 조절자**. 신경 가소성 및 약물 단서 학습(Drug cue learning)에 직결됨. |
| **`Klhl35`** | 153.77 Mb | **검정 불가**<br>(코딩 SNP 결측) | **$2.40 \times 10^{-8}$**<br>($\text{FDR} = 5.64 \times 10^{-5}$) | **+37개** | Kelch-like 단백질 35. CUL3-RING 유비퀴틴 리가아제 기질 어댑터로서 단백질 분해 항상성 조절. **Bonferroni 완벽 통과 ($Z = 5.46$)**. |
| **`Lrtomt`** | 156.26 Mb | **검정 불가**<br>(코딩 SNP 결측) | **$1.60 \times 10^{-7}$**<br>($\text{FDR} = 1.76 \times 10^{-4}$) | **+8개** | Leucine-rich transmembrane / O-메틸전달효소. 카테콜아민 메틸화(COMT) 유사 도메인 보유. **Bonferroni 완벽 통과 ($Z = 5.11$)**. |
| **`Rnf169`** | 154.25 Mb | **검정 불가**<br>(코딩 SNP 결측) | **$1.61 \times 10^{-7}$**<br>($\text{FDR} = 1.76 \times 10^{-4}$) | **+11개** | RING finger 단백질 169. 염색질 수준의 DNA 손상 및 복구 인자 경쟁 조절. **Bonferroni 완벽 통과 ($Z = 5.11$)**. |

---

### (3) `crf_ny_lever_presses` 신규 발굴 20개 유전자 전체 목록 (`BH-FDR < 0.05`)

1. **`Klhl35`** (Chr 1: 153.77 Mb) — $P = 2.40 \times 10^{-8}$ ($\text{FDR} = 5.64 \times 10^{-5}$) [Bonferroni 통과]
2. **`Inppl1`** (Chr 1: 156.18 Mb, SHIP2) — cMAGMA $P = 0.0022$ $\rightarrow$ H-MAGMA **$P = 1.03 \times 10^{-7}$ ($\text{FDR} = 1.69 \times 10^{-4}$)** [Bonferroni 통과, 시냅스 신호전달]
3. **`Lrtomt`** (Chr 1: 156.26 Mb) — $P = 1.60 \times 10^{-7}$ ($\text{FDR} = 1.76 \times 10^{-4}$) [Bonferroni 통과]
4. **`Rnf169`** (Chr 1: 154.25 Mb) — $P = 1.61 \times 10^{-7}$ ($\text{FDR} = 1.76 \times 10^{-4}$) [Bonferroni 통과]
5. **`Chrdl2`** (Chr 1: 154.34 Mb) — cMAGMA $P = 0.000197$ $\rightarrow$ H-MAGMA **$P = 1.34 \times 10^{-6}$** [Bonferroni 통과]
6. **`Lamtor1`** (Chr 1: 156.27 Mb) — $P = 2.26 \times 10^{-6}$ ($\text{FDR} = 0.00124$) [Bonferroni 통과, mTORC1]
7. **`Folr2`** (Chr 1: 156.20 Mb) — cMAGMA $P = 0.000435$ $\rightarrow$ H-MAGMA **$P = 2.62 \times 10^{-6}$** [Bonferroni 통과]
8. **`Fam181b`** (Chr 1: 147.15 Mb) — $P = 4.33 \times 10^{-5}$ ($\text{FDR} = 0.0137$)
9. **`Aqp11`** (Chr 1: 152.05 Mb) — cMAGMA $P = 0.000496$ $\rightarrow$ H-MAGMA $P = 5.90 \times 10^{-5}$ ($\text{FDR} = 0.0179$)
10. **`Chrna10`** (Chr 1: 156.49 Mb, 니코틴성 콜린성 수용체 $\alpha 10$) — cMAGMA $P = 0.00583$ $\rightarrow$ H-MAGMA **$P = 7.72 \times 10^{-5}$ ($\text{FDR} = 0.0226$)**
11. **`Art5`** (Chr 1: 156.47 Mb) — cMAGMA $P = 0.00958$ $\rightarrow$ H-MAGMA $P = 1.08 \times 10^{-4}$ ($\text{FDR} = 0.0290$)
12. **`Ddias`** (Chr 1: 146.91 Mb) — cMAGMA $P = 0.000550$ $\rightarrow$ H-MAGMA $P = 1.17 \times 10^{-4}$ ($\text{FDR} = 0.0299$)
13. **`Art1`** (Chr 1: 156.48 Mb) — cMAGMA $P = 0.00250$ $\rightarrow$ H-MAGMA $P = 1.21 \times 10^{-4}$ ($\text{FDR} = 0.0300$)
14. **`Prcp`** (Chr 1: 146.93 Mb) — cMAGMA $P = 0.000636$ $\rightarrow$ H-MAGMA $P = 1.26 \times 10^{-4}$ ($\text{FDR} = 0.0303$)
15. **`Ccdc90b`** (Chr 1: 146.63 Mb) — $P = 1.49 \times 10^{-4}$ ($\text{FDR} = 0.0350$)
16. **`Dgat2`** (Chr 1: 153.45 Mb) — cMAGMA $P = 0.00108$ $\rightarrow$ H-MAGMA $P = 1.89 \times 10^{-4}$ ($\text{FDR} = 0.0418$)
17. **`Alg8`** (Chr 1: 151.68 Mb) — $P = 2.02 \times 10^{-4}$ ($\text{FDR} = 0.0442$)
18. **`Ucp2`** (Chr 1: 154.84 Mb) — $P = 2.19 \times 10^{-4}$ ($\text{FDR} = 0.0458$)
19. **`St8sia1`** (Chr 4: 175.79 Mb) — $P = 2.20 \times 10^{-4}$ ($\text{FDR} = 0.0458$) [루프 변이 추가 미미]
20. **`Pdap1`** (Chr 12: 9.47 Mb) — $P = 2.20 \times 10^{-4}$ ($\text{FDR} = 0.0458$) [루프 변이 추가 미미]

---

## 3. 심층 비판적 검토: 6대 논리적 비약 및 방법론적 주의점

### ⚠️ [논점 1] 20개 신규 유전자의 "독립성 착시" (LD Hitchhiking / Co-localization)
* **문제점:**
  - `crf_ny_lever_presses` 신규 유전자 20개 중 **18개(90%)가 염색체 1번의 146.6 ~ 156.5 Mb (약 9.8 Mb 구간) 단 하나의 거대 QTL 영역에 집중**되어 있음.
  - 나머지 2개(*St8sia1*, *Pdap1*)는 루프 SNP 추가가 거의 없어 P-value 변화가 없었음.
* **논리적 비약 위험:**
  - 논문에서 이를 *"게놈 전역에서 20개의 독립적인 신규 원인 유전자를 발굴했다"*고 서술하면 명백한 과장이자 비약임.
  - 실제로는 Chr 1의 단일 강력한 Lead QTL 변이의 LD 영향권 안에 있는 이웃 유전자들이 크로마틴 루프를 통해 동일한 연관 신호를 공유하면서 집단으로 유의수준을 통과한 **LD Hitchhiking(동반 상승)** 현상임.
* **방어 논리 (수정 가이드):**
  - "독립적인 20개 발굴"이 아니라, **"Chr 1 거대 행동 조절 QTL 궤적(146~156 Mb) 내부에서, 3D 조절 루프를 통해 비코딩 변이들이 *Lamtor1*(mTORC1), *Inppl1*(시냅스 신호), *Chrna10*(콜린성 수용체)의 조절 영역에 물리적으로 연결됨을 밝혀냄으로써, 단일 QTL 궤적 내의 유력 타깃 유전자들을 좁혀내는 '정밀 우선순위화(Fine-mapping / Prioritization)'를 달성했다"**로 기술해야 완벽히 방어됨.

---

### Response 1: 부분 타당, 단일 LD 블록·fine-mapping 단정은 기각

- **확인:** legacy lever-press H-MAGMA-only 20개 중 18개가 rn7 Chr1:146,633,036–156,491,476에 위치한다. 독립적인 원인 유전자 20개로 해석하면 안 된다는 경고는 타당하다.
- **리뷰 정정:** 가까운 위치만으로 하나의 LD 블록, 하나의 원인 변이, 또는 LD hitchhiking이 입증되지는 않는다. SNP 간 LD와 조건부 독립성을 직접 검정하지 않았으므로 이러한 설명은 가능성이다. 여러 형질에서의 반복도 독립 코호트 재현이 아니다.
- **추가 확인:** P값이 그대로인 H-MAGMA-only 유전자는 St8sia1/Pdap1뿐 아니라 **Alg8/Ucp2까지 4개**이다. 이 네 유전자의 유의성 전환은 BH 보정 결과의 변화이며, 새로운 SNP가 강화한 신호로 세지 않는다.
- **대응:** “loop-informed candidate-gene prioritization”으로 제한한다. credible set/조건부 분석 없이 “fine-mapping 달성”, “핵심 타깃 규명”, “완벽히 방어”라고 쓰지 않는다. 이번에는 위 수치·P값을 재검산했으며 LD/조건부 분석은 수행하지 않았다. [검증표](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/review_audit/legacy_HMAGMA_only_hits_checked.tsv)

---

### ⚠️ [논점 2] cMAGMA "결측(Missing)"과 "비유의(Negative)"의 혼동
* **문제점:**
  - 4대 핵심 유전자(*Lamtor1*, *Klhl35*, *Lrtomt*, *Rnf169*)는 cMAGMA에서 "검정 결과 $P > 0.05$로 비유의했던 것"이 아님.
  - cMAGMA 로컬 주석(코딩/프로모터)에 등록된 소수의 SNP이 King GWAS 요약 통계량 파일과 매칭되지 않아, **cMAGMA `genes.out` 결과 파일에서 아예 행 자체가 생성되지 않은 '결측(Missing / Not Tested)' 상태**였음.
* **논리적 비약 위험:**
  - *"기존 cMAGMA에서는 통계적으로 유의하지 않았던 음성(Negative) 유전자가 H-MAGMA에서 양성(Positive)으로 전환되었다"*고 쓰면 사실관계 오류임.
* **올바른 기술 방식:**
  - *"기존 코딩 변이 중심 주석(cMAGMA)으로는 유전자 영역 내 변이 부족으로 **아예 검정조차 불가능했던 유전자들**이, H-MAGMA의 3D 원거리 조절 루프를 통해 비코딩 활성 변이가 할당됨으로써 비로소 검정 가능한 영역으로 편입되어 게놈 전역 유의성을 입증받았다"*고 정확하게 명시해야 함.

---

### Response 2: 핵심 지적 타당, 검정 불가 사유와 SNP 수를 명확히 함

- **확인:** Lamtor1/Klhl35/Lrtomt/Rnf169의 직접 exon/promoter 주석 SNP는 각각 1/2/2/9개지만, 네 King 형질의 GWAS에서 매칭되는 직접 SNP는 모두 0개였다. cMAGMA의 비유의 결과가 아니라 **해당 입력에서 not tested**인 경우다. P=1 또는 음성 결과로 대체하지 않는다.
- **수치 구분:** legacy H-MAGMA에서 실제 테스트된 SNP는 같은 순서로 **5/22/8/9개**이다. 참조 전체에서 추가된 5/37/8/11개를 그대로 해당 GWAS의 사용 SNP 수라고 부르면 안 된다.
- **대응:** “원거리 SNP 연결로 새롭게 검정 가능해졌고 해당 형질의 보정 기준을 통과했다”로 기술한다. “활성 변이”, “인과성 입증”은 제외하며, BH와 Bonferroni 통과 여부는 각각의 형질·유전자에 대해 구분한다. [네 유전자 직접 점검 결과](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/review_audit/four_candidate_mapping_audit.tsv)

---

### ⚠️ [논점 3] Many-to-One liftOver 시 Minimum P-value 선택에 따른 Anti-conservative 편향
* **코드 결함 (`run_selected_king_magma.py` 172~178행):**
  - King GWAS(rn6)를 rn7로 liftOver할 때 여러 rn6 변이가 하나의 rn7 위치로 매핑되는 다대일 충돌 시, 코드가 **무조건 더 작은(더 유의한) P-value를 선택**(`if p < prev: values[rn7_snp] = p`)하도록 작성됨.
* **영향 및 방어:**
  - 중복 위치에 최소 P값을 할당하는 것은 통계적으로 약간의 신호 인플레이션(Anti-conservative bias)을 유발할 수 있음. 다행히 충돌 변이 수는 전체 수백만 개 중 극소수이므로 전체 결론을 훼손하지 않으나, 방법론 기술 시 *"대다수 변이는 1:1로 보존되었으며 극소수 중복 위치의 최소 P값 선택이 게놈 인플레이션($\lambda \approx 1.18$)에 영향을 주지 않았음"*을 확인해 두어야 함.

---

### Response 3: 코드 지적 타당, 이번 유전자 결과에 대한 영향은 없음

- **확인:** 최소 P값 선택 분기가 존재하며, rn7 충돌 위치는 **2:22410246, 13:16826941** 두 곳이다. 원래 rn6 SNP 4개의 P값과 실제 유지된 값도 네 형질 모두 대조했다. 유의성을 보고 중복 변이를 선택하는 규칙은 일반적으로 정당화하기 어렵다.
- **실제 영향:** 두 위치는 LD 참조에는 있으나 **cMAGMA, legacy H-MAGMA, revised H-MAGMA 어느 주석에도 SNP–유전자 연결이 0개**이다. 따라서 이 두 P값은 이번 유전자 검정에 사용되지 않았다. 이 이유로 현행 loop-only 비교의 입력을 바꾸거나 MAGMA를 재실행할 필요는 없다.
- **리뷰 정정:** “λ≈1.18”은 이 네 입력의 공통 수치가 아니다. 저장된 rn7 P값 전체에서 계산한 λGC는 incentive value 1.340746, lever presses 1.372856, response bias 1.292661, index 1.130213이다. 두 충돌 위치 제외 전후 값은 같지만, 이것만으로 전체 분석이 편향 없다고 결론내리지 않는다.
- **대응:** 향후 일반 전처리에서는 다대일 위치를 표시하고 allele/reference 확인 또는 모호한 매핑 제외 규칙을 사전에 정한다. 이번 비교에서는 기존 입력을 보존하고 위의 무할당 확인을 근거로 남겼다. [충돌·주석 점검](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/review_audit/many_to_one_targets.tsv), [원래 P값 대조](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/review_audit/many_to_one_P_value_audit.tsv)

---

### ⚠️ [논점 4] GWAS 코호트(Round 8) vs LD 참조 패널(v4) 간의 세대 격차 및 3.5% 누락
* **문제점:**
  - King 2025 원본 GWAS는 **HS rat round 8 (Rnor_6.0)** 유전형 기준임.
  - 로컬 MAGMA 분석에서는 **`HS_genotypes_v4` (rn7, round 10 이후 세대 포함)** 패널을 LD 참조(bfile)로 사용함.
  - 전체 변이의 약 **3.5%가 LD 패널과 매칭되지 않아 누락(Missing)**되었으며, 세대 격차에 따른 LD 구조의 미세한 왜곡이 존재할 수 있음.

---

### Response 4: 외부 LD 참조의 한계는 타당, 세대·왜곡의 단정은 근거 부족

- **확인:** 각 King rn7 입력 3,505,415개 중 3,382,703개가 LD 참조와 매칭되고, 122,712개(**3.500641%**)가 매칭되지 않는다. 이는 rn6→rn7 liftOver 실패율(약 0.23%)과 다른 단계·분모의 수치이다.
- **리뷰 정정:** `round8`, `v4` 파일/자료 버전을 곧바로 생물학적 번식 세대로 해석하면 안 된다. King 논문이 명시한 실제 번식 세대는 **71–88**이다. 현재 자료만으로 v4 구성원의 세대 차이, 비매칭의 정확한 원인, LD 왜곡의 존재·크기를 확정할 수 없다. [King Methods](https://onlinelibrary.wiley.com/doi/10.1111/gbb.70018)
- **대응:** 외부 HS LD 참조가 GWAS 표본의 LD를 완전히 재현한다고 보장하지 않는다는 한계를 유지한다. loop-only 비교에서는 동일 참조를 고정했다. 코호트 일치 참조와의 LD/결과 민감도 검정은 별도 분석이며 이번에 수행했다고 쓰지 않는다. [분모·매칭 검산](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/King2025/review_audit/legacy_count_and_LD_audit.tsv)

---

### ⚠️ [논점 5] In silico 통계적 연관성과 "기능적 인과성(Causality)"의 혼동
* **주의사항:**
  - H-MAGMA는 컴퓨터 상에서(In silico) 비코딩 변이와 프로모터 간의 물리적 루프를 이어준 뒤 계산한 **통계적 연관성(Association)** 모델임.
  - *Lamtor1*이나 *Inppl1*의 P-value가 아무리 높더라도, 이것이 곧 **"해당 루프의 변이가 실제로 유전자 발현을 변화시켜 레버 누르기 행동을 직접 유발한다"는 인과관계(Causality)를 증명한 것은 아님.**
* **권고:**
  - 본문에서 *'demonstrates causal regulation'* 같은 단어를 피하고, **'prioritizes high-confidence candidate genes'**, **'mechanistic hypothesis'** 수준의 엄격하고 절제된 학술적 용어를 유지해야 함.

---

### Response 5: 타당, 인과·기능·발현 방향의 주장을 제한

- H-MAGMA는 입력 루프에 근거해 SNP를 유전자에 연결한 **유전자 단위 연관 검정**이다. 변이→발현→행동의 인과 경로나 개별 루프의 조절 기능을 검정한 것이 아니다. ATAC 겹침은 접근성의 증거이지 기능적으로 활성인 변이라는 증거가 아니다.
- “P-value가 높더라도”는 “P값이 낮더라도” 또는 “연관성이 강하더라도”로 고쳐 읽어야 한다. 또한 “high-confidence”는 독립 증거 없이 부여하지 않고 **candidate genes supported by loop-informed association** 정도로 표현한다.
- **대응:** 현재 보고서의 association/putative 한계를 유지한다. eQTL 공위치 검정, 변이/루프 교란 등은 추가 근거 후보이지 이번에 완료한 검증이 아니다. 일반 유전자 기능만으로 Lamtor1 등 특정 후보의 행동 기전을 확정하지 않는다. [H-MAGMA 원 논문](https://pmc.ncbi.nlm.nih.gov/articles/PMC7131892/)

---

### ⚠️ [논점 6] 프로모터 윈도우 정의의 불일치 (루프 2,001 bp vs H-MAGMA 2,500 bp)
* **문제점:**
  - 루프 필터링 기준은 **TSS $\pm 1,000$ bp (대칭 2,001 bp)**였으나, H-MAGMA 프로모터 주석 스크립트는 **TSS 기준 상류 2,000 bp / 하류 500 bp (비대칭 2,500 bp)**를 사용함.
  - 루프가 연결한 유전자 영역과 실제 H-MAGMA가 계산한 프로모터 경계면에 미세한 불일치가 발생할 수 있으므로, *"Its gene assignments therefore need not equal the revision's direct 2,001-bp loop-gene list"*라는 기존 기술 논리를 유지해야 함.

---

### Response 6: 정의 차이는 맞지만 의도된 비교 설계이며 길이는 2,501 bp

- **수치 정정:** 1-based closed 좌표에서 TSS−2,000부터 TSS+500까지는 **2,501 bp**이다(염색체 시작에서 잘리는 예외 제외). 음의 가닥은 방향을 반대로 적용한다. 2,500은 양쪽 offset의 합이다.
- **설계:** 2,001-bp 창은 **새 루프 13,376개를 선택하는 기준**, 2,501-bp 비대칭 promoter는 **기존 H-MAGMA SNP–gene 매핑 규칙**이다. 이번 요청은 루프만 바꾸고 기존 규칙을 유지하는 것이므로 이 차이 자체가 구현 오류는 아니다. 단, 두 단계의 유전자 할당이 같다고 주장하면 오류다.
- **버전 구분:** 리뷰 표의 legacy 루프를 2,001-bp 기준으로 선택했다고 서술하면 안 된다. 이 항목은 현재 revised 루프와 H-MAGMA 규칙을 비교하는 경우에만 적용된다.
- **대응:** 입력 루프 선택 창과 downstream promoter를 각각 명시했고 legacy SNP–gene 연결 323,827개가 정확히 재현되는 대조를 통과했다. downstream promoter를 2,001 bp로 통일하는 분석은 별도 민감도 분석이며 이번 loop-only 비교에 섞지 않는다. [매핑 코드](/Users/pete/Desktop/playground/enhancer/r_files/h3_GWAS/revision_2001bp/01_prepare_annotation.R)

---

## 4. 논문 리비전 및 저널 투고 액션 플랜

1. **메인 스토리라인 전환 (1순위):**
   - *"게놈 전역 20개 신규 유전자 발굴"* $\rightarrow$ **"Chr 1 주요 보상행동 QTL(146~156 Mb) 궤적 내부의 비코딩 조절 타깃 집중 규명(Locus fine-mapping & gene prioritization)"**.
   - **Response 4.1:** 과장 방지 취지는 동의한다. 그러나 fine-mapping/규명 역시 현재 증거를 넘으므로 “공간 접촉 정보를 이용한 연관 유전자 후보 제시”로 제한한다(Response 1, 5).
2. **대표 후보 유전자의 생물학적 조명:**
   - 리소좀-mTORC1 영양/보상 신호 축의 핵심인 **`Lamtor1`**, 시냅스 신호전달 포스파타아제인 **`Inppl1`**, 콜린성 니코틴 수용체인 **`Chrna10`**을 중심으로 Locus 3D 접촉 플롯을 제시하여 실질적인 신경생물학적 기전을 강조할 것.
   - **Response 4.2:** 위치·SNP·루프를 보여주는 플롯은 적절하다. 다만 후보 선정은 새 루프 결과와 실제 SNP 근거를 우선 확인하고, 알려진 기능만으로 후보의 기전·우선순위를 확정하지 않는다.
3. **독립 논문(Standalone Paper) 추진 시:**
   - King 2025(단서 반응성 / 유인 가치)와 Nicotine(약물 자가투여 / *Irs2* 축)을 결합하여, *"단서 반응부터 약물 중독 형성까지 이어지는 PFC 3D 조절 유전체 지도"*를 테마로 *Neuropsychopharmacology* 또는 *Addiction Biology* 투고를 적극 고려할 것.
   - **Response 4.3:** 지금 결과만으로 논문 성립·저널 적합성·“Irs2 축”을 결론내릴 수 없다. 기존 nicotine Day 5 결과에서 Irs2는 cMAGMA BH q=0.1415, legacy H-MAGMA BH q=0.2212로 둘 다 비유의였다. **새 nicotine 분석도 완료:** Irs2에 SNP 19개가 추가되어 사용 SNP는 32→51개, gene P는 3.4075e-5→1.6885e-5로 변했으나 BH q=0.1321로 여전히 비유의다. Day 5 전체에서도 BH q<0.05 유전자는 0개다. 서로 다른 형질의 연관 분석을 시간적 진행·중독 기전으로 연결하지 않는다. [Nicotine 최종 보고서](/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/revision_2001bp_nicotine_loop_only/RESULTS_REPORT.md)
