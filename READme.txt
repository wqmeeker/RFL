The purpose of this repo is to hold the codes needed to run and test the Stan/rstan model for the random fatigue limit (RFL) model.

R, rstan, and Rtools (for windows) need to be installed.

Start an instance of R with the wd in this folder.

Send all of the commands in RunRFL.R to run the test.

 1 May 2022
After fixing my error in the prior specification for mu_log_gamma, I
am experimenting with specification of minimumly informative priors
for sigma_error and the parameters for the distribution of the gamma values.
