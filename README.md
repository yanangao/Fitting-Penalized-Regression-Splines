# Fitting-Penalized-Regression-Splines-
## Problem setup and goal
- Fit a penalized regression splines with p = 3 (i.e., cubic). the number of knots is 30, and they are equi-spaced within the domain of the data.
- Model assumption: y = f + e, where E(e) = 0 and COV(e) = σ^2 * I.
- Conduct a simulation study to compare four methods of choosing parameter λ: GCV, CV, Corrected AIC and L2 risk.
## Simulation study
The experimental set-up adopted in this project is the same as Lee(2003), which is designed to study the effects of (i) noise level, (ii) design density, (iii) degree of spatial variation and (iv) noise variance. 
## Citation 
Lee, T. C. (2003). Smoothing parameter selection for smoothing splines: a simulation study. Computational statistics & Data analysis, 42(1-2), 139-148.
