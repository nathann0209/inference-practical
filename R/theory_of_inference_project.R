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

  # Ensure parameters are valid
  if (theta1 < 0 || s2 <= 0) {
    return(rep(1e-10, length(fin)))  # Return small positive number instead of 0
  }

  # Calculate expected final counts
  expected_fin <- init * theta1 * exp(-theta2 * times)
  
  # Standard deviation scales with initial count
  sd_fin <- init * sqrt(s2)
  
  # Calculate likelihood using the correct error structure
  likelihoods <- dnorm(fin, mean = expected_fin, sd = sd_fin)
  
  return(likelihoods)
}

# Function to fit Marie's model using maximum likelihood
fit_marie_model <- function(surv) {
  # Extract data from surv
  times <- surv$times
  init <- surv$init
  fin <- surv$fin
  
  # Initial parameter guesses
  # (Replace initial_theta1 and initial_theta2 with estimates from Q1 when available)
  initial_theta1 <- 1.0
  initial_theta2 <- 0.2
  
  # Estimate initial s2 based on a rough fit of the model using initial guesses for parameters
  expected_fin <- init * initial_theta1 * exp(-initial_theta2 * times)
  residuals <- fin - expected_fin
  descaled_residuals <- residuals / init
  initial_s2 <- mean(descaled_residuals^2)
  
  initial_params <- c(initial_theta1, initial_theta2, initial_s2)

  # Print initial values for debugging
  cat("Initial params:", initial_params, "\n")


  # Minimize negative log-likelihood
  optim_result <- optim(
    par = initial_params,
    fn = neg_log_likelihood,
    likelihood_func = marie_likelihood,
    times = times,
    init = init,
    fin = fin,
    method = "BFGS",
    hessian = TRUE,
    control = list(trace = 1, maxit = 1000)  # Add tracing for debugging
  )

  # Check convergence
  if (optim_result$convergence != 0) {
    warning("Optimization did not converge. Code: ", optim_result$convergence)
  }

  return(optim_result)
}


plot_marie_fit <- function(surv, fit) {
  times <- surv$times
  init <- surv$init
  fin <- surv$fin
  
  # Extract the fitted parameters
  theta1_hat <- fit$par[1]
  theta2_hat <- fit$par[2]
  s2_hat <- fit$par[3]
  
  # Calculate expected final counts
  expected_fin <- init * theta1_hat * exp(-theta2_hat * times)
  
  # Calculate 95% prediction intervals (error scales with init)
  ci_multiplier <- 1.96
  sd_fin <- init * sqrt(s2_hat)
  lower_bound <- expected_fin - ci_multiplier * sd_fin
  upper_bound <- expected_fin + ci_multiplier * sd_fin

  # Plot observed final counts against times
  plot(times, fin, pch = 19, col = "blue",
       xlab = "Time (hours)", ylab = "Final Count (fin)",
       main = "Marie's Model Fit with 95% Prediction Intervals")

  # To draw a smooth curve, order the data by times
  ord <- order(times)
  lines(times[ord], expected_fin[ord], col = "red", lwd = 2)

  # Add the prediction interval as a shaded region
  polygon(c(times[ord], rev(times[ord])),
          c(lower_bound[ord], rev(upper_bound[ord])),
          col = rgb(1, 0, 0, 0.2), border = NA)

  legend("topright", legend = c("Observed", "Fitted", "95% Prediction Interval"),
         col = c("blue", "red", rgb(1, 0, 0, 0.2)),
         pch = c(19, NA, 15), lty = c(NA, 1, NA), bty = "n")
}


# Fit the model
marie_fit <- fit_marie_model(surv)

# Ensure the save directory exists
if (!dir.exists("../plots")) {
  dir.create("../plots")
}

# Create and save plot of the model fit with 95% prediction intervals
png(filename = "../plots/marie_fit.png", width = 800, height = 600)
plot_marie_fit(surv, marie_fit)
dev.off()


