surv <- read.csv("../data/surv.csv")

# Function to compute negative log-likelihood 
neg_log_likelihood <- function(params, likelihood_func, times, init, fin) {
  # Compute the likelihood for the given parameters
  likelihood <- likelihood_func(params, times, init, fin)
  
  # Handle edge cases
  if (any(is.na(likelihood)) || any(likelihood <= 0)) {
    return(Inf)  # Return Inf if likelihood is invalid
  }
  
  # Compute the log-likelihood
  log_likelihood <- sum(log(likelihood))
  
  # Return the negative log-likelihood (for minimisation)
  return(-log_likelihood)
}

# Marie's model likelihood function
marie_likelihood <- function(params, times, init, fin) {
  # Extract parameters
  theta1 <- params[1]
  theta2 <- params[2]
  s2 <- params[3]  # variance parameter
  
  # Ensure theta1 is non-negative
  if (theta1 < 0) return(rep(0, length(fin)))
  
  # Calculate the deterministic part of the model
  deterministic_part <- theta1 * exp(-theta2 * times)
  
  # Since the model is: fin = init * (deterministic_part + error)
  # We can rearrange to: fin/init = deterministic_part + error
  # So the error term is: fin/init - deterministic_part
  
  # Calculate likelihood using the correct error structure
  likelihoods <- dnorm(fin/init, mean = deterministic_part, sd = sqrt(s2))
  
  return(likelihoods)
}

# Function to fit Marie's model using maximum likelihood
fit_marie_model <- function(surv) {
  # Extract data
  times <- surv$times
  init <- surv$init
  fin <- surv$fin
  
  # Initial parameter guesses
  # Using the hint: we might want to use estimates from Question 1
  # For now, setting reasonable initial values as we don't have Question 1 stuff yet
  initial_theta1 <- 1.0  # Replace with estimate from Question 1 when available
  initial_theta2 <- 0.2  # Replace with estimate from Question 1 when available
  initial_s2 <- var(fin/init)  # Basic variance estimate
  
  initial_params <- c(initial_theta1, initial_theta2, initial_s2)
}