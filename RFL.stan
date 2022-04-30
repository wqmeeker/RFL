

//
//  Random Fatigue Limit (RFL) model where mu depends on stress
//    and there is a random fatigue limit that varies from unit to unit
//

data {
int<lower=1> N;//                   number of rows in the data matrix
int<lower=1, upper=2> M;//          number of columns in the tran_response matrix
matrix[N,M]  tran_response;//       transformed response (take log above)
int<lower=1, upper=4> cencode[N];// censor codes 1-failure, 2-right, 3-left, 4-interval
real weights[N];//                  case weights or counts
real stress[N];//                   stress variable
real log_s0_minus_gamma0;//         centering value for the log stress variable (can set to 0)
int<lower=2, upper=4> numdist; //   distribution code 2-weibull 4-lognormal
//
// parameter codes and values for prior
//
int<lower=0, upper=5> prior_codes[5]; // prior codes for the five parameters
//
//   all  Stan-level priors are for unrestricted parameters (only 0, 1 and 5 implemented so far)
//      0 flat 
//      1 normal (implied by lognormal for originally positive parameters)
//      5 location-scale-t (log-location-scale-t for originally positive parameters)
//
//    prior_parameters are location, scale, dof (dof needed only for l-s-t)
//      note the parameter order
//
real beta0star_prior_parameters[3];
real beta1_prior_parameters[3];
real log_sigma_error_prior_parameters[3]; 
real mu_log_gamma_prior_parameters[3];
real log_sigma_log_gamma_prior_parameters[3];  
}

parameters{
real beta0star;//	mu intercept after stress translation
real beta1;//		mu slope
real log_sigma_error;//	log error sigma	 
real mu_log_gamma;//	mu for log gamma	
real log_sigma_log_gamma;//	log sigma for log gamma	 
real Z_log_gamma[N];// Norm(0,1) log fatigue_limit standardized random fatigue limit
}

transformed parameters{
real beta0;//			mu intercept
real<lower=0> sigma_error;//	error term
real<lower=0> sigma_log_gamma;  //scale parameter for log gamma
beta0 =  beta0star - beta1*log_s0_minus_gamma0;
sigma_error = exp(log_sigma_error);
sigma_log_gamma = exp(log_sigma_log_gamma);
//print(", b0=", beta0, ", b1 = ", beta1,  ",  s_e= ", sigma_error,  ", mu-lg= ", mu_log_gamma,  ", s-lg= ", sigma_log_gamma);
}

model{
real mu; // log life location parameter depends on stress
real gamma;  // random fatigue limit for each observation
//
//   standardized log RFL is normal(0, 1) The following are equivalent
//
//   Z_log_gamma ~ normal(0, 1);
     target += normal_lpdf(Z_log_gamma | 0, 1);
//
for(n in 1:N){
//
//    compute gamma for the current observation
//
      gamma = exp(mu_log_gamma + Z_log_gamma[n]*sigma_log_gamma);
//
// fix log likelihood for values of gamma > stress[n]
//
     if(gamma >= stress[n]){
//     print("obs= ", n, " gamma too big: stress = ", stress[n],  ", gamma = ", gamma, ", cen = ", cencode[n],  ", b0=", beta0, ", b1 = ", beta1,  ",  s_e= ", sigma_error,  ", Z= ", Z_log_gamma[n], ", mu-lg= ", mu_log_gamma,  ", s-lg= ", sigma_log_gamma);
//
//   for a censored observation the likelihood of being censored is 1, 
//          so add nothing [log(1) = 0] for this observation and continue
//
       if(cencode[n]==2)continue;
       else{
//
// A failure (of any kind) is impossible so add in something close to -Inf and  [log(0) = -Inf] continue
//
      target += log(1.e-300);
	  continue;
	  }
     }
//
//    compute mu for the current observation
//
     mu = beta0 + beta1*log(stress[n] - gamma);
//
// for each observation, add in the likelihood contributions
//       according to the distribution and censoring type
//
// if(n==1)print("logN= ", tran_response[n, 1], ", S= ", stress[n], ", g= ", gamma, ", b0*= ", beta0star, ", b0= ", beta0, ", b1= ", beta1, ", mu= ", mu, ", sig-e= ", sigma_error,  ", mu-lg= ", mu_log_gamma,  ", sig-lg= ", sigma_log_gamma, ", log density before= ", target());
//
if(numdist==4){
//
//  lognormal
//
//  lognormal exact
   if (cencode[n]==1)
      target += weights[n]*normal_lpdf(tran_response[n, 1]| mu, sigma_error);
//  lognormal right censored
   if (cencode[n]==2)
      target += weights[n]*normal_lccdf(tran_response[n, 1]| mu, sigma_error);
//  lognormal left censored
   if (cencode[n]==3)
      target += weights[n]*normal_lcdf(tran_response[n, 1]| mu, sigma_error);
//  lognormal interval censored
   if (cencode[n]==4)
      target += weights[n]*log(normal_cdf(tran_response[n, 2], mu, sigma_error) -
                               normal_cdf(tran_response[n, 1], mu, sigma_error));
}
else if (numdist==2){
//
//  Weibull
//
//  Weibull exact
   if (cencode[n]==1)
      target += weights[n]*gumbel_lpdf(tran_response[n, 1]| mu, sigma_error);
//  Weibull right censored
   if (cencode[n]==2)
      target += weights[n]*gumbel_lccdf(tran_response[n, 1]| mu, sigma_error);
//  Weibull left censored
   if (cencode[n]==3)
      target += weights[n]*gumbel_lcdf(tran_response[n, 1]| mu, sigma_error);
//  Weibull interval censored
   if (cencode[n]==4)
      target += weights[n]*log(gumbel_cdf(tran_response[n, 2], mu, sigma_error) -
                            gumbel_cdf(tran_response[n, 1], mu,sigma_error));
}

else{
 reject("numdist incorrect value; found numdist=", numdist);
}
// if(n==1)print("logN= ", tran_response[n, 1], " S= ", stress[n], " g= ", gamma, " b0*= ", beta0star, " b0= ", beta0, " b1= ", beta1, " mu= ", mu,  " sig-e= ", sigma_error,  " mu-lg= ", mu_log_gamma,  " sig-lg= ", sigma_log_gamma, "log density after= ", target());			    
}
//
// specify prior distributions
//
//     priors for beta0star, beta1, log_sigma_error, mu_log_gamma, log_sigma_log_gamma
//
// if prior_codes==0 flat (default --- do not need to do anything)
// if prior_codes==1 normal  (for informative and weakly informative)
//        generally arises from a lognormal prior on originally positive parameters
// if prior_codes==5 loc-sc-t (for informative and weakly informative)
//        generally arises from a log-loc-sc-t prior on originally positive parameters
//
//
//     normal distribution parameters are mean and standard deviation
//     student_t parameters are dof, location, and scale
//
//    beta0star prior
//
if(prior_codes[1]==1)
beta0star ~ normal(beta0star_prior_parameters[1],  beta0star_prior_parameters[2]);
if(prior_codes[1]==5)
beta0star ~ student_t(beta0star_prior_parameters[3],  beta0star_prior_parameters[1],  beta0star_prior_parameters[2]);
//
//   beta1 prior
//
if(prior_codes[2]==1)
beta1 ~ normal(beta1_prior_parameters[1],  beta1_prior_parameters[2]);
if(prior_codes[2]==5)
beta1 ~ student_t(beta1_prior_parameters[3],  beta1_prior_parameters[1],  beta1_prior_parameters[2]);
//
//  log_sigma_error prior 
//
if(prior_codes[3]==1)
log_sigma_error ~ normal(log_sigma_error_prior_parameters[1],  log_sigma_error_prior_parameters[2]);
if(prior_codes[3]==5)
log_sigma_error ~ student_t(log_sigma_error_prior_parameters[3],  log_sigma_error_prior_parameters[1],  log_sigma_error_prior_parameters[2]);
//
//  mu_log_gamma prior 
//
if(prior_codes[4]==1)
mu_log_gamma ~ normal(mu_log_gamma_prior_parameters[1],  mu_log_gamma_prior_parameters[2]);
if(prior_codes[4]==5)
mu_log_gamma ~ student_t(mu_log_gamma_prior_parameters[3],  mu_log_gamma_prior_parameters[1],  mu_log_gamma_prior_parameters[2]);
//
//  log_sigma_log_gamma prior 
//
if(prior_codes[5]==1)
log_sigma_log_gamma ~ normal(log_sigma_log_gamma_prior_parameters[1],  log_sigma_log_gamma_prior_parameters[2]);
if(prior_codes[5]==5)
log_sigma_log_gamma ~ student_t(log_sigma_log_gamma_prior_parameters[3],  log_sigma_log_gamma_prior_parameters[1],  log_sigma_log_gamma_prior_parameters[2]);
}

