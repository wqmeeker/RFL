The purpose of this repo is to hold the codes needed to run and test
the Stan/rstan model for the random fatigue limit (RFL) model.

R, rstan, and Rtools (for windows) need to be installed.

Start an instance of R with the wd in this folder.

Send all of the commands in RunRFL.R to run the test.

As of 29 April 2022, the answers do not agree with the Bayesian analysis by Johnson_etal1999Technometrics discussion and maximum likelihood results from the PascualMeeker1999Technometrics paper (which do agree). 

The estimates of gamma and the distribution of gamma should be centerd around 220 (somewhat below the smallest level of stress with failures (270). But exp(mu_log_gamma) = exp(1.702) = 5.48. We suspect that the distribution of the random parameter gamma is not being specified correctly. We experimented with a Jacobian correction, but that did not work.
