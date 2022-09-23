library(e1071)
library(dplyr)
library(data.table)
library(ggplot2)

#' Define a function to simulate the effects of violating multiple assumptions 
#' 
#' @param Means The mean values to test. Required.
#' @param Variances The variances to test. Required. 
#' @param Skews The skews to test. Required. 
#' @param ErrorThreshold How far off simulated values can be from measured values. 
#'     Default is 0.5.
#' @param Iterations The number of point to estimate in a distribution. Default is 10000.
#'
#' @return A large data.table of input summary statistics, actual simulated summary
#'    statistics, and the median RI score. 
ri_score_multiple_violations <- function(Means,
                                         Variances,
                                         Skews,
                                         ErrorThreshold = 0.5,
                                         Iterations = 10000) {
  
  # Create a dataframe of all possibilities to iterate through
  Parameters <- expand.grid(Means, Variances, Skews) 
  colnames(Parameters) <- c("Mean", "Variance", "Skew")
  
  # Set target (error threshold)
  targ <- ErrorThreshold
  
  # Define a function to calculate RI score
  RIscore <- function(RI_d, sigma = 3) {
    exp(-1 * (RI_d^2/(2 * sigma^2)))
  }
  
  # Get erlang distribution parameters
  get_erd <- function(variance, skew) {
    kappa <- (2/skew)^2
    lambda <- sqrt(kappa/variance)
    return(c(kappa, lambda))
  }
  
  # Generate simulation
  Simulation <- do.call(rbind, lapply(1:nrow(Parameters), function(row){
    
    # Pull out expected values
    MEAN <- Parameters$Mean[row]
    SKEW <- Parameters$Skew[row]
    VAR <- Parameters$Variance[row]
    
    message(paste0("Mean: ", MEAN, ", Skew: ", SKEW, ", Var: ", VAR))
    
    # Generate the target value
    if (MEAN == 0) {MEAN_TARGET = c(-targ, targ)} else {MEAN_TARGET = c(MEAN - MEAN*targ, MEAN + MEAN*targ)}
    if (SKEW == 0) {SKEW_TARGET = c(-targ, targ)} else {SKEW_TARGET = c(SKEW - SKEW*targ, SKEW + SKEW*targ)}
    VAR_TARGET = c(VAR - VAR*targ, VAR + VAR*targ)
    
    # Run rnorm is the skew is 0
    if (SKEW == 0) {
      
      # Generate simulated values
      values <- rnorm(10000, mean = MEAN, sd = sqrt(VAR))
      
      # Get simulated mean, skew, sd
      MEAN_SIM <- mean(values)
      SKEW_SIM <- skewness(values)
      VAR_SIM <- var(values)
      
      # Simulated values must be within a 10% window of the target values
      while(MEAN_SIM >= max(MEAN_TARGET) | MEAN_SIM <= min(MEAN_TARGET) |
            SKEW_SIM >= max(SKEW_TARGET) | SKEW_SIM <= min(SKEW_TARGET) |
            VAR_SIM >= max(VAR_TARGET) | VAR_SIM <= min(VAR_TARGET)) {
        
        values <- rnorm(10000, mean = MEAN, sd = sqrt(VAR))
        
        # Get simulated mean, skew, sd
        MEAN_SIM <- mean(values)
        SKEW_SIM <- skewness(values)
        VAR_SIM <- var(values)
        
      }
      
      # Once we have simulated values within the 10% of the targeted simulated value, #return results
      return(data.table(
        Mean = MEAN,
        SD = sqrt(VAR),
        Skew = SKEW,
        Mean_Sim = MEAN_SIM,
        SD_Sim = sqrt(VAR_SIM),
        Skew_Sim = SKEW_SIM,
        Median_RI_Score = median(RIscore(values))
      ))
      
    } else {
      
      # Get erlang distribution paramters
      params <- get_erd(VAR, SKEW)
      
      # Use gamma
      values <- rgamma(10000, shape = params[1], rate = params[2]) 
      if (SKEW < 0) {values <- values * -1}
      
      # Adjust mean
      if (round(mean(values)) != MEAN) {
        
        # Get difference
        diff <- MEAN - mean(values)
        values <- values + diff
        
      }
      
      # Get simulated mean, skew, var
      MEAN_SIM <- mean(values)
      SKEW_SIM <- skewness(values)
      VAR_SIM <- var(values)
      
      # Simulated values must be within a 50% window of the target values
      while(MEAN_SIM >= max(MEAN_TARGET) | MEAN_SIM <= min(MEAN_TARGET) |
            SKEW_SIM >= max(SKEW_TARGET) | SKEW_SIM <= min(SKEW_TARGET) |
            VAR_SIM >= max(VAR_TARGET) | VAR_SIM <= min(VAR_TARGET)) {
        
        # Use gamma
        values <- rgamma(10000, shape = params[1], rate = params[2]) 
        if (SKEW < 0) {values <- values * -1}
        
        if (round(mean(values)) != MEAN) {
          
          # Get difference
          diff <- MEAN - mean(values)
          values <- values + diff
          
        }
        
        # Get simulated mean, skew, sd
        MEAN_SIM <- mean(values)
        SKEW_SIM <- skewness(values)
        VAR_SIM <- var(values)
        
      }
      
      return(data.table(
        Mean = MEAN,
        SD = sqrt(VAR),
        Skew = SKEW,
        Mean_Sim = MEAN_SIM,
        SD_Sim = sqrt(VAR_SIM),
        Skew_Sim = SKEW_SIM,
        Median_RI_Score = median(RIscore(values))
      ))
      
    }
    
  })) %>% data.table()
  
  return(Simulation)
  
}

SimulationResults <- ri_score_multiple_violations(
  Means = seq(-5, 5, 1),
  Variances = c(0.5, 1, 2, 3, 4, 5, 6)^2,
  Skews = seq(-5, 5, 1)
)

ggplot(SimulationResults, aes(x = Mean, y = Mean_Sim)) + geom_point()
ggplot(SimulationResults, aes(x = SD, y = SD_Sim)) + geom_point()
ggplot(SimulationResults, aes(x = Skew, y = Skew_Sim)) + geom_point()
write.csv(SimulationResults, "~/Downloads/Multiple_Violations_Simulation.csv", quote = F, row.names = F)



