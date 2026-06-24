## MMF Outbreak Investigation
How  would housing or behavioral changes reduce the number of infections during the MMF outbreak at MMED ?
---
## Overview
We explore the epideomological modelling of Muizenburg Mathamatics Fever (MMF) outbreak at MMED. This repository contains data and R code for fitting a two-location SIR compartmental model to the MMF outbreak. The outbreak involved 34 cases across two locations (AIMS and Empire) over 5 days (15-19 June 2026).

---
Repository Structure
```
MMF/
├── data/
│   ├── MMF-Final+Locations+Doctors.csv   # cleaned outbreak data (34 cases)
│   ├── MMF-Final+Duplicates.csv          # raw data with duplicate infections
│   └── MMF-DoctorVisits.csv              # doctor visit records
│
├── scripts/
│   ├── Data_Cleaning.R                   # data cleaning and preparation
│   ├── VisualiseOriginalData.R           # exploratory visualisation
|   ├── MLE_Fitting.R                     # Maximum Likelihood Estimation model fitting for transmission coefficients
│   ├── Deterministic_ODE_Base.R          # main model: deterministic + MLE
|   ├── Gillespie.R                       # stochastic Gillespie simulation
|   ├── Housing_intervention.R            # deterministic predictive models should some housing change or quarantine take place
|   ├── Behavioural_intervention.R        # deterministic predictive models should some behavioural change take place
|   └── Final_Networks.ipynb                     # Python code to express the social network within the MMED workshop
│
└── README.md
```
---

## Done Prior to Model and Prediction Building
# In Data_Cleaning.R:
Entries edited to be uniform,
Missing entries added,
Location of Housing added, 0 if infected individual lives in the AIMS building, and 1 if infected individual lives in the Empire building,
MMF-DoctorVisits.R file combined to the original data, 
Added infections having no effect on the model, ie the infections to already recovered individuals
# In MLE_Fitting.R:
Used MLE to estimate transmission coefficient for mixing groups during the day and within groups at night
## Begin Base Model Building
# In Deterministic_ODE_Base.R:
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
Plotted next to collected data

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
source("Deterministic_ODE_Base.R")
```
---
# In Gillespie.R:
Used the same model world and description as for the deterministic model
Modelled a random Gillespie Stochastic model
Found the average of 1000 random infectious curves over a 10 day period 
####### Someone needs to type in how to run stuff because I'm confused how the typing works

# In Housing_intervention.R:
Used the same model world and description as for the deterministic model
Modelled the infectious curve for the base model and 3 different housing interventions
These interventions include: 
Everyone mixing at the same rate throughout the day and night
Everyone isolating at night, ie no interactions at night for anyone
Everyone swapping which building they live in, ie the sizes of the populations swapped
These were modelled for the AIMS and Empire populations separately

# In Behavioural_intervention.R:
Used the same model world and description as for the deterministic model
Modelled the infectious curve for the base model and 1 different behavioural change
This behavioural change introduced the following probabilities into the transmission coefficients:
Multiplied by 0.03 if both the infectious and susceptible individuals observed the behavioural change
Multiplied by 0.14 if the susceptible individual observed the behavioural change and the infectious individual did not
Multiplied by 0.07 if the infectious individual observed the behavioural change and the susceptible individual did not
Multiplied by 1 if neither the infectious nor susceptible individuals observed the behavioural change

The behavioural change would also differ between night and day:
Probability of observing behavioural change during the day = 0.85
Probability of observing behavioural change during the day = 0.3

# In Final_Networks.ipynb
Showed different network setups of transmissions between specific individuals in the data


Authors
MMF 2026 Outbreak Investigation Group  
African Institute for Mathematical Sciences (AIMS), Muizenberg, RSA
---
Acknowledgements
Adapted from the ICI3D Fitting Tutorial  
(Bellan 2015, Blumberg 2025, Pearson 2026)
