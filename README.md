## MMF Outbreak Investigation
How  would housing or behavioral changes reduce the number of infections during the MMF outbreak at MMED ?
---
## Overview
We explore the epideomological modelling of Muizenburg Mathamatics Fever (MMF) outbreak at MMED. This repository contains data and R code for fitting a two-location SIR compartmental model to the MMF outbreak. The outbreak involved 34 cases across two locations (AIMS and Empire) over 4 days (15-19 June 2026).

---
Repository Structure
```
MMF/
├── data/
│   ├── MMF-Final+Locations+Doctors.csv   # cleaned outbreak data (34 cases)
│   ├── MMF-Final+Duplicates.csv          # raw data with duplicates
│   └── MMF-DoctorVisits.csv              # doctor visit records
│
├── scripts/
│   ├── Data_Cleaning.R                   # data cleaning and preparation
│   ├── VisualiseOriginalData.R           # exploratory visualisation
│   ├── SSIIRR_Full_Model.R              # main model: deterministic + MLE
│   └── Gillespie.R                       # stochastic Gillespie simulation
│
└── README.md
```
---

## Model Description

We fit a deterministic ODE model with two population groups:
Compartment	Description
S_A	Susceptible individuals at AIMS
S_E	Susceptible individuals at Empire
I_A	Infectious individuals at AIMS
I_E	Infectious individuals at Empire
R_A	Recovered individuals at AIMS
R_E	Recovered individuals at Empire
Transmission is time-varying:
Daytime (08:00–18:30): both groups mix — governed by `beta_M`
Night-time: within-group only — governed by `beta_NA` (AIMS) and `beta_NE` (Empire)

## Key Assumptions

Closed population: N = 42 ($N_A$ = 25, $N_E$ = 17) <br>
Infectious period = 24 hours → $γ_A$ = $γ_E$ = 1 <br>
3 index cases: Kimberley (Empire), Morgan and Mandie (AIMS) <br>
$R0$ ≈ $beta_M$ (since $γ = 1$) <br>

## Parameter Estimates

Parameter	Prior	MLE	95% CI <br>
$beta_M$ (daytime)	1.8	fitted	from Hessian <br>
$beta_NA$ (AIMS night)	1.5	fixed	<br>
$beta_NE$ (Empire night)	0	fixed	<br>
$gamma_A$	1	fixed	<br>
$gamma_E$	1	fixed	<br>

## How to Run

Clone the repository:
```bash
git clone https://github.com/vixenmaep/MMF.git
```
Set your working directory in R:
```r
setwd("path/to/MMF/MMF")
```
Install required packages:
```r
install.packages(c("deSolve", "ggplot2", "ellipse", "tidyr"))
```
Run the main model script:
```r
source("SSIIRR_Full_Model.R")
```
---
Authors
MMF 2026 Outbreak Investigation Group  
African Institute for Mathematical Sciences (AIMS), Muizenberg, RSA
---
Acknowledgements
Adapted from the ICI3D Fitting Tutorial  
(Bellan 2015, Blumberg 2025, Pearson 2026)
