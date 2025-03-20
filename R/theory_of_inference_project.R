surv <- read.csv("../data/surv.csv")

# Extract data from surv
times <- surv$times
init <- surv$init
fin <- surv$fin

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

# multiplicative_model's model likelihood function
multiplicative_model_likelihood <- function(params, times, init, fin) {
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
  
  # Calculate likelihood
  likelihoods <- dnorm(fin, mean = expected_fin, sd = sd_fin)
  
  return(likelihoods)
}

# log-linear model's likelihood function
log_linear_model_likelihood <- function(params, times, init, fin) {
  # Extract parameters
  phi1 <- params[1]
  phi2 <- params[2]
  sigma2 <- params[3]  # variance parameter
  
  # Ensure parameters are valid
  if (sigma2 <= 0) {
    # Return small positive likelihood if parameters are invalid
    return(rep(1e-10, length(fin)))
  }
  
  # Calculate mean of log(fin) under the model
  expected_log_fin <- log(init) + phi1 + phi2 * times
  
  # Standard deviation in log-scale
  sd_log_fin <- sqrt(sigma2)
  
  # Calculate likelihood using normal distribution on log(fin)
  likelihoods <- dnorm(log(fin), mean = expected_log_fin, sd = sd_log_fin)
  
  return(likelihoods)
}

optim_wrapper <- function(initial_params, likelihood_func) {
  # Perform optimisation
  optim_result <- optim(
    par = initial_params,
    fn = neg_log_likelihood,
    likelihood_func = likelihood_func,
    times = times,
    init = init,
    fin = fin,
    method = "BFGS",
    hessian = TRUE,
    control = list(trace = 1, maxit = 1000)  # Debug/tracing
  )
  
  # Check convergence
  if (optim_result$convergence != 0) {
    warning("Optimisation did not converge. Code: ", optim_result$convergence)
  }
  
  return(optim_result)
}
# Function to find mle estimates for multiplicative model
fit_multiplicative_model <- function() {
  
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
  optim_result <- optim_wrapper(initial_params, multiplicative_model_likelihood)

  return(optim_result)
}
# Function to find mle estimates for log linear model
fit_log_linear_model <- function() {
  
  # Initial parameter guesses
  initial_phi1 <- 0.0
  initial_phi2 <- 0.0
  
  # Estimate initial sigma^2 based on a rough fit of the model using initial guesses
  expected_log_fin <- log(init) + initial_phi1 + initial_phi2 * times
  residuals <- log(fin) - expected_log_fin
  initial_sigma2 <- mean(residuals^2)
  
  # Bundle into an initial parameter vector
  initial_params <- c(initial_phi1, initial_phi2, initial_sigma2)
  
  # Print initial values for debugging
  cat("Initial params:", initial_params, "\n")
  
  # Minimize negative log-likelihood
  optim_result <- optim_wrapper(initial_params, log_linear_model_likelihood)
  
  return(optim_result)
}

plot_model_fit <- function(param1, param2, variance_parameter,
                           model_name = "Model Name") {
  
  if (model_name == "Multiplicative Model") {
      predicted_ratio <- param1 * exp(-param2 * times)
      sd_ratio <- sqrt(variance_parameter)
  } else if (model_name == "Log Linear Model") {
      predicted_ratio <- exp(param1 + param2 * times)
      sd_ratio <- predicted_ratio * exp(0.5 * variance_parameter) * sqrt(exp(variance_parameter) - 1)
  } else {
      stop("Unsupported model type. Please use 'Multiplicative Model' or 'Log Linear Model'.")
  }
  
  
  # Compute observed ratio
  observed_ratio <- fin / init
  
  # Take sqrt(variance_parameter) as the standard deviation
  ci_multiplier <- 1.96
  
  # 95% prediction (or confidence) interval
  ci_upper <- predicted_ratio + ci_multiplier * sd_ratio
  ci_lower <- predicted_ratio - ci_multiplier * sd_ratio
  
  # Plot observed data
  plot(times, observed_ratio,
       pch = 19, col = "blue",
       xlab = "Time (hours)",
       ylab = "Relative Change (fin / init)",
       main = paste(model_name, "Fit: Relative Change vs Time"))
  
  
  # Order data for a smooth fitted line
  ord <- order(times)
  
  # Add fitted curve
  lines(times[ord], predicted_ratio[ord], col = "red", lwd = 2)
  
  # Add prediction/confidence interval as a shaded area
  polygon(
    x = c(times[ord], rev(times[ord])),
    y = c(ci_lower[ord], rev(ci_upper[ord])),
    col = rgb(1, 0, 0, 0.2),
    border = NA
  )
  
  # Add legend
  legend("topright",
         legend = c("Observed", "Fitted", "95% Interval"),
         col = c("blue", "red", rgb(1, 0, 0, 0.2)),
         pch = c(19, NA, 15),
         lty = c(NA, 1, NA),
         bty = "n")
}

# Function to analyse residuals 
analyse_residuals <- function(param1, param2, variance_parameter,
                                                   model_name = "Model Name") {
  # Calculate standardised residuals
  if (model_name == "Multiplicative Model") {
    expected_fin <- init *param1* exp(-param2* times)
    raw_residuals <- fin - expected_fin
    standardised_residuals <- raw_residuals / (init * sqrt(variance_parameter))
  } else if (model_name == "Log Linear Model") {
    expected_fin <- exp(log(init) + param1 + param2*times)
    raw_residuals <- log(fin) - log(expected_fin)
    standardised_residuals <- raw_residuals / sqrt(variance_parameter)
  } else {
    stop("Unsupported model type. Please use 'Multiplicative Model' or 'Log Linear Model'.")
  }
  
  # Create a panel of diagnostic plots
  par(mfrow = c(2, 2))
  
  # Plot 1: Residuals vs Fitted
  plot(expected_fin, standardised_residuals, 
       main = "Residuals vs Fitted",
       xlab = "Fitted Values", 
       ylab = "Standardised Residuals")
  abline(h = 0, col = "red", lty = 2)
  
  # Plot 2: Residuals vs Time
  plot(times, standardised_residuals, 
       main = "Residuals vs Time",
       xlab = "Time (hours)", 
       ylab = "Standardised Residuals")
  abline(h = 0, col = "red", lty = 2)
  
  # Plot 3: Normal Q-Q Plot
  qqnorm(standardised_residuals)
  qqline(standardised_residuals, col = "red")
  
  # Plot 4: Histogram of residuals with normal density overlay
  hist(standardised_residuals, 
       main = "Histogram of Residuals",
       xlab = "standardised Residuals", 
       probability = TRUE,
       breaks = 15)
  curve(dnorm(x, mean = mean(standardised_residuals), sd = sd(standardised_residuals)), 
        add = TRUE, col = "red", lwd = 2)
  
  # Reset the plotting parameters
  par(mfrow = c(1, 1))
  
  
  # Return residuals and test result for reporting
  return(list(
    residuals = raw_residuals,
    standardised_residuals = standardised_residuals,
    mean_residual = mean(standardised_residuals),
    sd_residual = sd(standardised_residuals)
  ))
}


# Fit multiplicative model
multiplicative_model_fit <- fit_multiplicative_model()

# Print multiplicative model mle estimates: 
theta1_hat <- multiplicative_model_fit$par[1]
theta2_hat <- multiplicative_model_fit$par[2]
s2_hat <- multiplicative_model_fit$par[3]
cat("theta1_hat: ", theta1_hat, "\n", "theta2_hat: ", theta2_hat, "\n", "s2_hat: ", s2_hat, "\n")

# Fit log linear model
log_linear_model_fit <- fit_log_linear_model()

# Print log linear model mle estimates: 
phi1_hat <- log_linear_model_fit$par[1]
phi2_hat <- log_linear_model_fit$par[2]
sigma2_hat <- log_linear_model_fit$par[3]
cat("phi1_hat: ", phi1_hat, "\n", "phi2_hat: ", phi2_hat, "\n", "sigma2_hat: ", sigma2_hat, "\n")


# Ensure the save directory exists
if (!dir.exists("../plots")) {
  dir.create("../plots")
}

# Create and save plot of the multiplicative model fit with 95% prediction intervals
png(filename = "../plots/multiplicative_model_fit.png", width = 800, height = 600)
plot_model_fit(theta1_hat, theta2_hat, s2_hat, "Multiplicative Model")
dev.off()

# Save the plot for multiplicative model residuals analysis
png(filename = "../plots/multiplicative_model_residuals.png", width = 800, height = 800)
residual_analysis_multiplicative <- analyse_residuals(theta1_hat, theta2_hat, s2_hat, "Multiplicative Model")
dev.off()

# Print summary statistics for report
cat("Residual analysis summary (multiplicative):\n")
cat("Mean of standardised residuals:", residual_analysis_multiplicative$mean_residual, "\n")
cat("SD of standardised residuals:", residual_analysis_multiplicative$sd_residual, "\n")

# Create and save plot of the log linear model fit with 95% prediction intervals
png(filename = "../plots/log_linear_model_fit.png", width = 800, height = 600)
plot_model_fit(phi1_hat, phi2_hat, sigma2_hat, "Log Linear Model")
dev.off()

# Save the plot for log linear model residuals analysis
png(filename = "../plots/log_linear_model_residuals.png", width = 800, height = 800)
residual_analysis_log_linear <- analyse_residuals(phi1_hat, phi2_hat, sigma2_hat, "Log Linear Model")
dev.off()

# Print summary statistics for report
cat("Residual analysis summary (log linear):\n")
cat("Mean of standardised residuals:", residual_analysis_log_linear$mean_residual, "\n")
cat("SD of standardised residuals:", residual_analysis_log_linear$sd_residual, "\n")
