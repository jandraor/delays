
# Andrade & Duggan (2023)

This repository contains code for the paper:

[Jair Andrade](https://www.linkedin.com/in/jandraor/) and [Jim
Duggan](https://ie.linkedin.com/in/jduggan). *Anchoring the mean
generation time in the SEIR to mitigate biases in $\Re_0$ estimates due
to uncertainty in the distribution of the epidemiological delays*.

The analysis in this study can be reproduced by executing the files:

- **S1.rmd**
- **S2.rmd**
- **S3.rmd**
- **S4.rmd**
- **S5.rmd**
- **S6.rmd**
- **S7.rmd**

## Abstract

The basic reproduction number, $\Re_0$, is of paramount importance in
the study of infectious disease dynamics. Primarily, $\Re_0$ serves as
an indicator of the transmission potential of an emerging infectious
disease and the effort required to control the invading pathogen.
However, its estimates from compartmental models are strongly
conditioned by assumptions in the model structure, such as the
distribution of the latent and infectious periods (epidemiological
delays). To further complicate matters, models with dissimilar delay
structures produce equivalent incidence dynamics. Following a simulation
study, we reveal that the nature of such equivalency stems from a linear
relationship between $\Re_0$ and the mean generation time, along with
adjustments to other parameters in the model. Leveraging this knowledge,
we propose and successfully test an alternative parameterisation of the
SEIR model that produces accurate $\Re_0$ estimates regardless of the
distribution of the epidemiological delays, at the expense of biases in
other quantities deemed of lesser importance. We further explore this
approach’s robustness by testing various transmissibility levels,
generation times, and data fidelity (overdispersion). Finally, we apply
the proposed approach to data from the 1918 influenza pandemic. We
anticipate that this work will mitigate potentially large errors in the
estimation of $\Re_0$..

## Resources

- [Stan](https://mc-stan.org/)

- [readsdr](https://github.com/jandraor/readsdr)
