rangePrior <-
function(lower,upper,distribution,center.prob=0.99){

  tail.prob <- 1 - (1-center.prob)/2
switch(distribution,

       lognormal={
         mu <- (log(lower)+log(upper))/2
         sigma <- (log(upper)-log(lower))/(2*qnorm(tail.prob))},

       normal={
         mu <- (lower+upper)/2
         sigma <- (upper-lower)/(2*qnorm(tail.prob))},
       {stop("distribution",distribution,"not recognized")})
result <- c(mu=mu, sigma=sigma)
result
}
