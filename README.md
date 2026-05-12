<div align="center">

<br/>

![Study Design](https://img.shields.io/badge/Design-Population--based%20Cohort-1e40af?style=flat-square)
![Method](https://img.shields.io/badge/Method-Marginal%20Structural%20Model-0891b2?style=flat-square)
![Data](https://img.shields.io/badge/Data-Hong%20Kong%20CDARS-059669?style=flat-square)
![Period](https://img.shields.io/badge/Period-2009--2023-6366f1?style=flat-square)
![Status](https://img.shields.io/badge/Status-Under%20Review-f59e0b?style=flat-square)

<p align="center">
  <b>Comparative safety profiles of biologic and targeted synthetic DMARDs <br/>
  in rheumatoid arthritis, using territory-wide electronic medical records <br/>
  and causal inference to address time-varying confounding.</b>
</p>

</div>

---

## 📖 Background

> Rheumatoid arthritis (RA) is a chronic autoimmune disease whose management increasingly relies on **biologic and targeted synthetic DMARDs (b/tsDMARDs)**. Although these drugs are effective, their comparative safety profiles remain unclear — randomised trials are too short and too small to detect rare adverse events, while conventional observational studies are vulnerable to **time-varying confounding** and **selection bias**.

This study addresses that evidence gap using a **territory-wide EMR database** and **marginal structural models (MSM)** to estimate real-world, head-to-head safety comparisons across b/tsDMARD classes.

---

## 🎯 Objectives

- Compare 5-year safety profiles across **four b/tsDMARD classes** in RA patients.
- Address time-varying confounding via **inverse probability weighting (IPW)**.
- Provide robust, clinically actionable evidence to support personalised therapy.

---

## 🧪 Drug Classes Compared

<table align="center">
<tr>
  <th align="center">Class</th>
  <th align="center">Agents</th>
</tr>
<tr>
  <td align="center"><b>TNFi</b></td>
  <td>adalimumab · certolizumab pegol · <b>etanercept</b> (ref) · golimumab · infliximab</td>
</tr>
<tr>
  <td align="center"><b>IL-6i</b></td>
  <td>sarilumab · tocilizumab</td>
</tr>
<tr>
  <td align="center"><b>Lymphocyte-targeting</b></td>
  <td>abatacept · rituximab</td>
</tr>
<tr>
  <td align="center"><b>JAKi</b></td>
  <td>baricitinib · tofacitinib · upadacitinib</td>
</tr>
</table>

---

## 🏥 Data Source

- **Database:** Hong Kong Clinical Data Analysis and Reporting System (**CDARS**)
- **Coverage:** ~80% of routine hospital admissions in Hong Kong (7.3M residents)
- **Study period:** 1 Jan 2009 – 31 Dec 2023
- **Population:** Adult RA patients (ICD-9-CM 714.0) initiating ≥1 b/tsDMARD, excluding other autoimmune indications and non-RA rituximab uses.

---

## 🔬 Methods

### Analytical Framework

```mermaid
flowchart LR
    A[Raw EMR<br/>CDARS] --> B[Cohort<br/>Construction]
    B --> C[Baseline IPTW]
    B --> D[Time-varying IPCW]
    C --> E[Stabilised IPW<br/>Pseudo-cohort]
    D --> E
    E --> F[Pooled Logistic<br/>Regression]
    F --> G[5-year Risk &<br/>Risk Difference]
    G --> H[Bootstrap CI<br/>+ Bonferroni]
```

### Key Methodological Features

<details>
<summary><b>📌 Marginal Structural Model (MSM)</b></summary>

Handles variables like **ESR, CRP, csDMARDs, glucocorticoids, NSAIDs, opioids** that act simultaneously as **mediators** and **time-varying confounders** — a setting in which conventional regression is biased.

</details>

<details>
<summary><b>📌 Stabilised Inverse Probability Weights</b></summary>

- **Baseline IPTW:** addresses confounding by indication
- **Time-varying IPCW:** addresses informative censoring from treatment switching
- Final weight = IPTW × stabilised IPCW, truncated at 1st/99th percentiles
- Mean stabilised weights: **0.87–1.02** across all arms
- Balance assessed via **Love plots** (SMD < 0.1 for most covariates)

</details>

<details>
<summary><b>📌 Outcome Modelling</b></summary>

- **Pooled logistic regression** with visit-by-agent interactions for non-proportional hazards
- **5-year absolute risks & risk differences** vs. etanercept (reference)
- **1,000 bootstrap replicates** for pointwise 95% percentile CIs
- **Bonferroni correction** across 9 outcomes → CI level = 1 − 0.05/9

</details>

<details>
<summary><b>📌 Competing Risk of Death</b></summary>

Discrete-time **multi-state model** with three states (event of interest / death / no event) to avoid bias from treating deaths as censored.

</details>

<details>
<summary><b>📌 Negative Control</b></summary>

**Fracture** risk is used as a negative control outcome — no plausible association with b/tsDMARD class is expected. Any detected difference would signal residual confounding.

</details>

<details>
<summary><b>📌 Missing Data</b></summary>

- ESR/CRP missing rates: **36%–44%** throughout follow-up
- **Multiple Imputation by Chained Equations (MICE)** applied
- Drug doses standardised by **WHO Defined Daily Dose (DDD)**

</details>

---

## 📊 Outcomes Examined

<div align="center">

| Cardiovascular | Infectious | Mental Health | Metabolic | Oncological | Other |
|:---:|:---:|:---:|:---:|:---:|:---:|
| MACE | Infection | Depression | Diabetes | Malignancy | Gastritis · Hospitalisation · Mortality |

*Plus **fracture** as negative control.*

</div>

---

## 📝 Citation

> *Manuscript under review. Full citation will be updated upon publication.*

---

## 🤝 Contact

For methodological questions or collaboration enquiries, please open an issue or contact the corresponding author.

<div align="center">

<sub>Built with ❤️ for reproducible clinical research</sub>

</div>
