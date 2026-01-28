# Mosquito Spatial Diffusion Modeling (SIT)

This project provides a comprehensive mathematical and numerical framework for the analysis, simulation, and parameter estimation of mosquito population control strategies, specifically the **Sterile Insect Technique (SIT)**. 

## Publication
Download our full paper: [Mosquito_Diffusion_Study.pdf](./Mosquito_diffusion_study.pdf)

## Project Overview
The core of this project is a deterministic spatial diffusion model describing the evolution of mosquito surface density $\rho(x, y, t)$. This study, conducted as part of a larger initiative led by **Facundo Muñoz, (CIRAD)**, aims to provide tools to optimize the release of sterile males to regulate vector-borne diseases.

The model is formulated as a **parabolic partial differential equation (PDE)**:
$$\frac{\partial \rho}{\partial t} = D \Delta \rho - \nu \rho - \rho \sum_{k=1}^{N_{traps}} f_k(x, y)$$

Where:
* **Diffusion ($D$)**: Movement and dispersion in the environment.
* **Mortality ($\nu$)**: Natural death rate of the population.
* **Capture ($f_k$)**: Impact of traps localized in the domain.


##  Key Features

### 1. Analytical Analysis
* Development of **analytical approximations** based on the truncation of the Duhamel series to understand early-stage diffusion.

### 2. Numerical Resolution & Benchmarking
The project implements and compares three major numerical methods to solve the diffusion-reaction equations:
* **Finite Difference Method (MDF)**: Implemented for its ease of use and computational speed.
* **Finite Element Method (MEF)**: Tested for handling complex geometries.
* **Finite Volume Method (MVF)**: Evaluated for its conservative properties.


### 3. Parameter Inference (Optimization)
A central part of this work is the **inverse problem**: estimating physical parameters from capture data. 
* **Algorithm**: Uses **L-BFGS-B** (Quasi-Newton method) to minimize the Mean Square Error (MSE).
* **Identifiability**: Addresses the critical "flat valley" identifiability issue between trap efficiency ($\gamma$) and trap radius ($R_{trap}$) by introducing functional constraints.
* **Performance**: Successfully converges in approximately 23 iterations with low error rates (~8-10%) for $D$ and $\nu$.



## Results Summary
The model confirms the duality between microscopic (individual Brownian motion) and macroscopic (PDE density) scales. The numerical framework proves robust enough to be applied to real-world SIT scenarios, providing physically consistent estimated values.

## Packages Used
python, numpy, pandas, matplotlib, scipy.optimize

## Author and Acknoledgment 
Aya Chaieb, Antonin Decouvelaere, Lucas Bourret, Baptiste Petiot, Mathieu Prioux, Facundo Muñoz

---
*This work is a continuation of research initiated by CIRAD regarding the calibration of probabilistic diffusion models.*