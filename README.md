# Gibbs-Sampling
Code zum Gibbs Sampling

# Gibbs sampler for baysian linear regressions as general function --------


# Introduction ------------------------------------------------------------

# The aim of this procedure is to describe a function that solves any gibbs sampler problem
# with the following structure of known and unknown information in baysian linear regression.
# The procedure is cunstructed such, that burn in and actual length can be compared as good
# as different randomly chosen starting values.


# Enviroment settings -----------------------------------------------------

rm(list = ls())

# Algorithm ---------------------------------------------------------------

# The function is named gibbs_sampler. We need different specifications to apply it.
# At first we specify the number of repetitions (reps), the length (R), and the burn in (burn).
# We need a data matrix X, a variable of interest y, the a priori parameters for the
# beta distribution beta_prior and diag_prior and the a priori parameters for 
# the inverse gamma distribution a_prior and b_prior

gibbs_sampler <- function(reps, R, burn,
                          X, y,
                          beta_prior, diag_prior,
                          a_prior, b_prior){
  
  # software requirments
  if(!require("MASS")) install.packages("MASS")
  library(MASS)
  if(!require("invgamma")) install.packages("invgamma")
  library(invgamma)
  if(!require("coda")) install.packages("coda")
  library(coda)
  
  # storage for the final results is created
  # The structure is a list per repetitions, which includes a list with two dataframes.
  # The first dataframe stands for the burn in and the second for the actual length of
  # the gibbs_sampler
  results <- list()
  for (i in 1:reps) {
    results[[i]] <- list(data.frame(), data.frame())
  }
  
  # we define a counter of repetitions for the while loop
  rep_counter <- 1
  
  while(rep_counter <= reps) {
    
    # within each repetition we create some temporary storage for the results of beta and sigma2
    betas <- matrix(NA, nrow = burn + R, ncol = length(beta_prior))
    sigma2s <- matrix(NA, nrow = burn + R, ncol = 1)
    
    # we need to set starting values, they are results from random draws
    # and differ in each repetition.
    beta <- c(runif(ncol(X), -10, 10))
    sigma2 <- runif(1, 0, 10)
    
    
    # here the actual gibbs sampler algorithm starts
    
    for (i in 1:(burn + R)) {
      
      # at first draw two parameters needed for multivariate normal
      var_beta <- solve(1/sigma2 * (t(X) %*% X) + diag_prior)
      nu_beta <- var_beta %*% (1/sigma2 * (t(X) %*% y) + (diag_prior %*% beta_prior))
      # plug them in to get a vector beta from a multivariate normal
      beta <- mvrnorm(1, nu_beta, var_beta)
      
      # use beta and prior information to calculate the rate parameter
      # from a inverse gamma
      rate_par <- 0.5 * t(y - X %*% beta) %*% (y - X %*% beta) + b_prior
      # use information to draw sigma2 from a inverse gamma
      sigma2 <- rinvgamma(1, a_prior * (N/2), rate_par)
      
      # assign results to the temporary storage
      betas[i,] <- beta
      sigma2s[i] <- sigma2
    }
    
    # create a dataframe by result matrices betas and sigma2s
    created_names <- character()
    for (i in 1:ncol(X)) {
      created_names[i] <- paste("beta", i-1, sep = "")
    }
    created_names <- c(created_names, "sigma2")
    eval_df <- cbind(as.data.frame(betas), as.data.frame(sigma2s))
    names(eval_df) <- created_names
    
    # split data up according burn in period and actual length of gibbs sampler
    burn_in <- eval_df[(1:burn),] 
    gibbs_results <- eval_df[((burn + 1):(burn + R)),]
    
    # assign results to the final storage as described in the beginning
    results[[rep_counter]][[1]] <- burn_in
    results[[rep_counter]][[2]] <- gibbs_results
    
    # go to the next rep
    rep_counter <- rep_counter + 1
  }
  
  # only return the list of results
  return(results)
}

# save the function
save(gibbs_sampler, file = "gibbs_sampler.RData")




# Evaluation of a simulated dataframe -------------------------------------


# Construction of data ----------------------------------------------------
# in case we did not save it by now
if(!require("MASS")) install.packages("MASS")
library(MASS)
if(!require("invgamma")) install.packages("invgamma")
library(invgamma)
if(!require("coda")) install.packages("coda")
library(coda)



# number of observations in dataset
N <- 10000
# parameters, here a vektor mu and a matrix sigma needs to be defined
mu_of_X <- c(rep(0, 2))
sigma_of_X <- matrix(c(1,  0.3,
                       0.3,  1), 
                     nrow = 2, byrow = T)
# draw with above parameter from multivariate normal
X <- cbind(rep(1,N), mvrnorm(N, mu_of_X, sigma_of_X))

# define "true" values to compare gibbs sampler afterwards and finish the data generation
true_beta <- c(1.5, 0.5, -0.5)
true_sigma <- 0.5
y <- X %*% true_beta + rnorm(N, 0, true_sigma)

# prior distribution for beta is a multivariate normal distribution with
beta_prior <- rep(0, 3) 
diag_prior <- 0.01 * diag(3)
# prior distribution for sigma2 is a inverse gamma with
a_prior <- runif(1); b_prior <- runif(1)
# these priors are conjugate, important for analytical derivation
# they determine the parameters, especially sigma2 a lot



# Apply the gibbs_sampler - function --------------------------------------

results_list <- gibbs_sampler(5, 1000, 500, 
                              X, y, beta_prior, diag_prior, a_prior, b_prior)

system.time(results_list <- gibbs_sampler(5, 1000, 500, 
                                          X, y, beta_prior, diag_prior, a_prior, b_prior))

rm(list=setdiff(ls(), c("gibbs_sampler", "results_list", "result_list")))



# Analysis of dataframes ---------------------------------------------------
windows()
par(mfrow=c(2,2), byrow = T)
autocorr.plot(results_list[[1]][[1]][,2], main = "autocorr. b1 in burn-in", auto.layout = F, col ="red", lwd = 2)
autocorr.plot(results_list[[1]][[2]][,2], main = "autocorr. b1 after burn-in", auto.layout = F, col ="blue",lwd = 2)
autocorr.plot(results_list[[1]][[1]][,4], main = "autocorr. sigma2 in burn-in", auto.layout = F, col ="red", lwd = 2)
autocorr.plot(results_list[[1]][[2]][,4], main = "autocorr. sigma2 after burn-in", auto.layout = F, col ="blue", lwd = 2)

res_1_gibbs <- as.mcmc(results_list[[1]][[2]])
res_2_gibbs <- as.mcmc(results_list[[2]][[2]])

windows()
par(mfrow=c(4,2))
traceplot(res_1_gibbs[,2], main = "b1, starting val. 1", col = "blue")
traceplot(res_2_gibbs[,2], main = "b1, starting val. 2", col = "darkblue")
plot(density(results_list[[1]][[2]][,2]), xlab = "", main = "", col = "blue", lwd = 2)
plot(density(results_list[[2]][[2]][,2]), xlab = "", main = "", col = "darkblue", lwd = 2)
traceplot(res_1_gibbs[,4], main = "sigma2, starting val. 1", col = "blue")
traceplot(res_2_gibbs[,4], main = "sigma2, starting val. 2", col = "darkblue")
plot(density(results_list[[1]][[2]][,4]), xlab = "", main = "", col = "blue", lwd = 2)
plot(density(results_list[[2]][[2]][,4]), xlab = "", main = "", col = "darkblue", lwd = 2)
