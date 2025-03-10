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
  
  # Calculate likelihood
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

  # Order the data by times to draw a smooth curve
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

# Function to analyse residuals for Marie's model
analyse_marie_residuals <- function(surv, fit) {
  # Extract data and fitted parameters
  times <- surv$times
  init <- surv$init
  fin <- surv$fin
  
  theta1_hat <- fit$par[1]
  theta2_hat <- fit$par[2]
  s2_hat <- fit$par[3]
  
  # Calculate expected final counts
  expected_fin <- init * theta1_hat * exp(-theta2_hat * times)
  
  # Calculate raw residuals
  raw_residuals <- fin - expected_fin
  
  # Calculate standardized residuals - scaled by init to match error structure
  standardized_residuals <- raw_residuals / (init * sqrt(s2_hat))
  
  # Create a panel of diagnostic plots
  par(mfrow = c(2, 2))
  
  # Plot 1: Residuals vs Fitted
  plot(expected_fin, standardized_residuals, 
       main = "Residuals vs Fitted",
       xlab = "Fitted Values", 
       ylab = "Standardized Residuals")
  abline(h = 0, col = "red", lty = 2)
  
  # Plot 2: Residuals vs Time
  plot(times, standardized_residuals, 
       main = "Residuals vs Time",
       xlab = "Time (hours)", 
       ylab = "Standardized Residuals")
  abline(h = 0, col = "red", lty = 2)
  
  # Plot 3: Normal Q-Q Plot
  qqnorm(standardized_residuals)
  qqline(standardized_residuals, col = "red")
  
  # Plot 4: Histogram of residuals with normal density overlay
  hist(standardized_residuals, 
       main = "Histogram of Residuals",
       xlab = "Standardized Residuals", 
       probability = TRUE,
       breaks = 15)
  curve(dnorm(x, mean = mean(standardized_residuals), sd = sd(standardized_residuals)), 
        add = TRUE, col = "red", lwd = 2)
  
  # Reset the plotting parameters
  par(mfrow = c(1, 1))
  
  # Shapiro-Wilk test for normality
  shapiro_test <- shapiro.test(standardized_residuals)
  
  # Return residuals and test result for reporting
  return(list(
    residuals = raw_residuals,
    standardized_residuals = standardized_residuals,
    shapiro_test = shapiro_test,
    mean_residual = mean(standardized_residuals),
    sd_residual = sd(standardized_residuals)
  ))
}

# Fit the model
marie_fit <- fit_marie_model(surv)

# analyse residuals
residual_analysis <- analyse_marie_residuals(surv, marie_fit)

# Ensure the save directory exists
if (!dir.exists("../plots")) {
  dir.create("../plots")
}

# Create and save plot of the model fit with 95% prediction intervals
png(filename = "../plots/marie_fit.png", width = 800, height = 600)
plot_marie_fit(surv, marie_fit)
dev.off()

# Save the plot for residuals analysis
png(filename = "../plots/marie_residuals.png", width = 800, height = 800)
analyse_marie_residuals(surv, marie_fit)
dev.off()

# Print summary statistics for report
cat("Residual analysis summary:\n")
cat("Mean of standardized residuals:", residual_analysis$mean_residual, "\n")
cat("SD of standardized residuals:", residual_analysis$sd_residual, "\n")
cat("Shapiro-Wilk test p-value:", residual_analysis$shapiro_test$p.value, "\n")


