# Rat Prefrontal Cortex (PFC) 매칭 고품질 GWAS 후보 데이터셋 종합 평가 보고서

**작성 일시**: 2026-09-05  
**목적**: `King2025`처럼 통계적 검정력이 높고($N \ge 1,200 \sim 1,600+$, Lead SNP $-\log_{10} P \ge 7.5 \sim 9.5+$), 우리 연구의 핵심 조직인 **Adult Rat Prefrontal Cortex (PFC) 인핸서-프로모터 루프**와 생물학적으로 완벽히 부합하면서, **현실적으로 Summary Statistics 확보 및 재현이 가능한 데이터셋**을 발굴 및 평가함.

---

## 1. 벤치마크 기준: 왜 King2025는 성공했고, 다른 데이터셋은 아쉬웠는가?

| 평가 요소 | King2025 (대성공 ★★★★★) | Nicotine Day 5 (아쉬움 ★★★☆☆) | Lara2024 (아쉬움 ★★☆☆☆) |
| :--- | :--- | :--- | :--- |
| **샘플 수 ($N$)** | **$N \approx 1,600$** (충분한 검정력) | $N \approx 1,000 \sim 1,200$ | $N = 629$ (검정력 절대 부족) |
| **Lead SNP 시그널** | **$-\log_{10} P \approx 9.5$** (매우 강력한 QTL) | $-\log_{10} P \approx 5.0$ (피크 자체가 미약) | $-\log_{10} P \approx 6.0$ (중간 수준) |
| **뇌 부위 매칭** | **PavCA / Cue Lever Press** (PFC 보상회로) | Nicotine SA (중독) | Delay Discounting (충동성) |
| **H-MAGMA 결과** | **FDR < 0.05 신규 유전자 7~20개 발굴** | *Irs2* P값 개선되었으나 q < 0.05 미달 | 신규 FDR < 0.05 유전자 0개 |
| **데이터 가용성** | **UCSD 공개 아카이브 완비 (`bb15123938`)** | 자체 보유 요약 통계량 완비 | 완벽 재현 검증 통과 (`.mlma` 완비) |

> **핵심 교훈**:
> 아무리 생물학적 부위가 잘 맞아도 **① $N < 1,000$이거나 Lead SNP $-\log_{10} P < 7.0$이면 H-MAGMA에서 살아남는 신규 유전자가 나오지 않습니다.**
> 반대로 이론상 아무리 완벽해도 **② 공공 리포지토리에 완전한 Summary Statistics가 없거나, 비공개 공변량/샘플 필터링으로 '완벽 재현'이 불가능하면 논문에 탑재할 수 없습니다.**

---

## 2. 주요 후보 데이터셋별 심층 진단: 장점 vs 현실적 장애물

### [후보 1] Oxycodone Self-Administration GWAS (Palmer Lab / Carrette et al. / Chen Lab)
* **생물학적 매칭도**: **★★★★★ (최상)**
  - 우리가 구축한 E-P loop의 핵심 레퍼런스인 snATAC-seq/TSR 공공 데이터([GSE193757](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193757), Duttke & Telese et al.) 자체가 바로 **"Adult rat PFC under naive, cocaine, and oxycodone conditions"**에서 생성된 것입니다.
  - 즉, 우리 인핸서 루프의 분자생물학적 원천과 오피오이드 갈망/재발 신경망(mPFC-NAc 글루타메이트 투사 회로)이 100% 일치합니다.
* **통계적 검정력**: $N \approx 1,200 \sim 1,500+$, 주요 염색체(Chr 1, Chr 18 등) 게놈 전장 유의 수준($-\log_{10} P \ge 7.5 \sim 8.5$).
* **[현실적 장애물 및 팩트 체크]**:
  - **공공 Summary Statistics 부재**: 2021년 Carrette et al. (eNeuro) 논문은 대규모 랫 바이오뱅크 구축에 관한 논문이며, 옥시코돈 GWAS 전장 요약 통계량(`.mlma` 등)은 아직 C-GORD나 UCSD Library에 독립 다운로드 패키지로 공식 릴리즈되지 않았습니다.
  - **가용성 상태**: **In-House (연구실 내부 보유 여부에 전적으로 의존)**. 만약 Hao Chen 교수님 연구실 내부 네트워크나 Palmer 연구실 협력 드라이브에서 직접 `.mlma` 요약 통계량 파일을 전달받지 못한다면, 외부에서 즉시 다운로드하여 돌릴 수 없습니다.

---

### [후보 2] Gunturkun et al. 2022 (Frontiers in Psychiatry / Palmer Lab)
* **표현형**: Open Field Test (OFT 총 이동거리, 중앙 구역 시간), Novel Object Interaction (NOIT 신기성 탐색), Social Interaction (SIT)
* **통계적 검정력**:
  - 샘플 수: $N = 1,246$ (GeneNetwork 전체 축적 코호트는 $10,000+$)
  - Lead SNP: **`chr11:33359859` ($-\log_{10} P = 8.268$)**, **`chr10:94549701` ($-\log_{10} P = 7.286$)** 등 다수의 강력한 QTL 보유.
* **생물학적 매칭도**: **★★★★☆ (매우 우수)**
  - 신기성 탐색, 불안 조절, 사회적 접근/회피는 Medial Prefrontal Cortex (mPFC) - Basolateral Amygdala (BLA) - NAc 회로의 대표적 행동.
* **[현실적 장애물 및 팩트 체크]**:
  - **완벽 재현 차단 (Blocker)**:
    1. 논문에 명시된 DOI(`10.48810/P44W2` / `10.48810/P44W2Q`)가 공공에서 작동하지 않아 C-GORD 원천 데이터 접근 불가.
    2. 논문 분석군 1,246마리의 정확한 나이(Age) 공변량 및 공변량 선택 회귀 기준(>2% 분산 설명)이 비공개 상태라, GeneNetwork의 `trait_10424`와 로컬 지노타입으로 자체 GCTA를 돌렸을 때 논문 Table 2 수치와 미세한 불일치가 발생함.
  - **가용성 상태**: **Blocked for Exact Reproduction (근사 분석은 가능하나 논문 완벽 재현 검증은 저자 문의 필요)**.

---

### [후보 3] Cocaine Self-Administration 확장 형질 스크리닝 (UCSD `bb2334903t`, de Guglielmo et al. 2024 / George et al.)
* **배경**:
  - 로컬 `Cocaine2026/` 폴더에서 테스트했던 4개 기본 형질(`mean_to_01_03`, `pc1_lga`, `total_intake`, `shock_03`)은 신규 유전자가 나오지 않았습니다.
  - 그러나 다운로드된 UCSD 원본 아카이브(`bb2334903t`)에는 위 4개 외에도 **수십 개의 세부 행동 형질 `.mlma` 파일**이 포함되어 있습니다.
* **주목할 확장 형질**:
  - **Progressive Ratio (PR break point)**: 약물 갈망 및 동기화 수준 (PFC 보상 평가 회로 의존)
  - **Locomotor Sensitization**: 행동 감작 (VTA-PFC 도파민 가소성)
  - **Extinction & Reinstatement**: 소거 및 재발 (mPFC의 대표적 억제성 조절 행동)
* **생물학적 매칭도**: **★★★★★ (snATAC 원천 일치)**
  - Duttke et al. (GSE193757)의 핵심 조건인 코카인 노출 PFC와 직결.
* **통계적 검정력**: $N \approx 1,000 \sim 1,200$, 세부 형질 중 Lead SNP $-\log_{10} P \ge 7.0$ 이상인 형질 선별 필요.
* **[현실적 장점]**:
  - **이미 로컬 외부 SSD에 전 염색체 `.mlma` 파일들이 100% 다운로드되어 있음!**
  - 새로운 데이터를 내려받거나 재현 문제를 고민할 필요 없이, 기존 `.mlma` 파일들 중 최고 P-value 피크가 높은 형질을 스크리닝하여 즉시 H-MAGMA에 투입 가능.

---

### [후보 4] Kuhn et al. 2025 Heroin Vulnerability (Molecular Psychiatry / Kalivas & Palmer)
* **표현형**: 헤로인 자가투여 총 섭취량, Break Point, Escalation (12h), Nociception (Tail flick)
* **생물학적 매칭도**: **★★★★☆** (오피오이드 의존성 및 통각 민감도)
* **통계적 검정력**: $N = 874$, Lead SNP $-\log_{10} P \approx 6.5 \sim 7.2$.
* **[현실적 진단]**:
  - **완벽 재현 검증 성공**: 로컬 파이프라인(`kuhn_gwas_reproducibility_notes.md`)을 통해 논문의 Top SNP, 위치, P-value, beta가 100% 정확하게 재현됨.
  - **단점**: King2025에 비해 샘플 수($N=874$)가 적어, H-MAGMA 특이적 신규 유전자가 폭발적으로 나오지는 않는 한계 확인.

---

### [후보 5] 5-Choice Serial Reaction Time Task (5-CSRTT) Impulsivity GWAS (Mitchell / Palmer Lab)
* **표현형**: Premature Responding (행동 억제 실패 = 운동성 충동성), Omissions (주의 집중 결함)
* **생물학적 매칭도**: **★★★★★ (충동성/집행 기능, Prefrontal Cortex 전형적 기능)**
  - `Lara2024` 지연할인(Delay discounting, $N=629$)의 낮은 검정력을 완벽하게 대체할 수 있는 충동성 태스크.
* **통계적 검정력**: $N \approx 1,200 \sim 1,400+$, Lead SNP $-\log_{10} P \ge 7.0 \sim 8.0$.
* **[현실적 상태]**:
  - GeneNetwork `HSNIH-PalmerPublish` 또는 C-GORD에서 5-CSRTT 관련 형질(Premature responses)의 전장 summary-stats 추출 가능 여부 확인 필요.

---

## 3. 후보 데이터셋 종합 비교 매트릭스

| 데이터셋 | 표현형 (Phenotype) | 샘플 수 ($N$) | Lead SNP $-\log_{10} P$ | PFC 뇌 부위 적합도 | 실제 데이터 가용성 및 재현 상태 | 종합 추천도 |
| :--- | :--- | :--- | :--- | :--- | :--- | :---: |
| **King2025** | PavCA / Lever Press | $\approx 1,600$ | **$\approx 9.5$** | ★★★★★ | **완벽 재현 및 분석 완료 (`bb15123938`)** | **현재 1위 (기준)** |
| **Oxycodone SA** | 오피오이드 자가투여/갈망 | $\approx 1,500+$ | $\approx 7.5 - 8.5$ | **★★★★★ (snATAC 100% 일치)** | **공공 .mlma 미출시 (Chen 연구실 내부 확보 필요)** | **학술적 1위 / 가용성 보류** |
| **Gunturkun2022** | OFT 이동거리 / NOIT | $1,246$ | **$8.27$ (Chr 11)** | ★★★★☆ | **로컬 파이프라인 보유 / 단, DOI 미작동으로 완벽재현 보류** | **검정력 2위 / 재현 보류** |
| **Cocaine 확장** | 코카인 PR 갈망, 감작, 소거 | $\approx 1,100$ | 형질별 상이 | ★★★★★ | **외부 SSD에 전 염색체 `.mlma` 보유 (즉시 선별 가능)** | **현실적 대안 1위** |
| **5-CSRTT** | 충동성 (Premature response) | $\approx 1,400$ | $\approx 7.0 - 8.0$ | ★★★★★ | **GeneNetwork 데이터 존재 / 요약통계량 변환 필요** | **추천 후보 2위** |
| **Kuhn2025** | 헤로인 총 섭취량, Nociception | $874$ | $\approx 6.5 - 7.0$ | ★★★★☆ | **완벽 재현 완료 / 단, H-MAGMA 특이 유전자 부족** | **기확인 완료** |

---

## 4. 실행 로드맵 및 현실적 대안 제언

### [전략 A] 즉시 실행 가능한 현실적 최선책 (추가 다운로드/재현 장벽 없음)
1. **Cocaine2026 확장 형질 스크리닝**:
   - 이미 외부 SSD에 다운로드된 `Cocaine2026`의 모든 `.mlma` 파일들을 대상으로 `awk` 한 줄 스크립트로 **Lead SNP $-\log_{10} P \ge 7.0$ 이상인 형질(예: PR break point, Extinction 등)**을 빠르게 찾아냅니다.
   - 높은 시그널을 가진 형질이 발견되면 즉시 cMAGMA / H-MAGMA를 가동하여 신규 유전자를 확보할 수 있습니다.

### [전략 B] 학술적 최고봉(Oxycodone / Gunturkun) 데이터 획득 액션
1. **Oxycodone GWAS**:
   - 연구실(Hao Chen 교수님) 내부 서버나 Dropbox/공유 폴더에 보관된 **GCTA `.mlma` 최종 요약 통계량**이 존재하는지 확인 및 파일 수령. (수령 즉시 최고 설득력의 분석 가능)
2. **Gunturkun2022**:
   - `GWAS@ratgenes.org` 또는 저자(Palmer Lab)에게 1,246마리의 정확한 샘플 ID 목록 및 나이 공변량 데이터를 요청하여 완벽 재현 검증 락을 해제함.
