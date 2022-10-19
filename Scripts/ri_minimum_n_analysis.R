library(dplyr)
library(data.table)
library(MASS)
library(parallel)
library(doParallel)
library(foreach)

# Pull all true positives 
Annotations <- fread("~/Downloads/AnnoTP.txt")

# Subset to N > 30
MetaboliteCount <- table(Annotations$`Compound Name`, dnn = c("Compound.Name")) %>%
  data.frame() 
Anno30 <- Annotations[Annotations$`Compound Name` %in% 
                        MetaboliteCount[MetaboliteCount$Freq >= 30, "Compound.Name"],] %>%
  filter(`Compound Name` != "[PNNLMET0040] Impurity 001 [12.148]")
Compounds <- unique(Anno30$`Compound Name`)

# Build for loop
registerDoParallel(detectCores())

for (metab in Compounds) {
  
  message(metab)
  
  # Pull the distribution
  dist <- Anno30[Anno30$`Compound Name` == metab, "Retention Index"] %>% unlist()
  
  results <- foreach(n = 3:30, .combine = rbind) %dopar% {
    
    repeats <- foreach(iteration = 1:100, .combine = rbind) %dopar% {
    
      # Subset the distribution
      takeSample <- sample(dist, n)
     
      # Take all the distributional fits
      Normal <- tryCatch({fitdistr(takeSample, "normal", lower = c(0,0))}, error = function(e) {
        list("estimate" = c("mean" = NA, "sd" = NA), "sd" = c("mean" = NA, "sd" = NA))
      }) 
      LogNormal <- tryCatch({fitdistr(takeSample, "lognormal", lower = c(0,0))}, error = function(e) {
        list("estimate" = c("meanlog" = NA, "sdlog" = NA), "sd" = c("meanlog" = NA, "sdlog" = NA))
      }) 
      Logistic <- tryCatch({fitdistr(takeSample, "logistic", lower = c(0,0))}, error = function(e) {
        list("estimate" = c("location" = NA, "scale" = NA), "sd" = c("location" = NA, "scale" = NA))
      }) 
      Gamma <- tryCatch({fitdistr(takeSample, "normal", lower = c(0,0))}, error = function(e) {
        list("estimate" = c("shape" = NA, "rate" = NA), "sd" = c("shape" = NA, "rate" = NA))
      }) 
      
      # Return the data.table of fits 
      return(
        data.table(
          "Compound" = metab,
          "Iteration" = iteration,
          "Normal.Mean.Est" = Normal$estimate[[1]],
          "Normal.Mean.SE" = Normal$sd[[1]],
          "Normal.SD.Est" = Normal$estimate[[2]],
          "Normal.SD.SE" = Normal$sd[[2]],
          "LogNormal.LogMean.Est" = LogNormal$estimate[[1]],
          "LogNormal.LogMean.SE" = LogNormal$sd[[1]],
          "LogNormal.LogSD.Est" = LogNormal$estimate[[2]],
          "LogNormal.LogSD.SE" = LogNormal$sd[[2]],
          "Logistic.Location.Est" = Logistic$estimate[[1]],
          "Logistic.Location.SE" = Logistic$sd[[1]],
          "Logistic.Scale.Est" = Logistic$estimate[[2]],
          "Logistic.Scale.SE" = Logistic$sd[[2]],
          "Gamma.Shape.Est" = Gamma$estimate[[1]],
          "Gamma.Shape.SE" = Gamma$sd[[1]],
          "Gamma.Rate.Est" = Gamma$estimate[[2]],
          "Gamma.Rate.SE" = Gamma$sd[[2]]
        )
      )
                        
    }
    
    repeats$N <- n
    return(repeats)
    
  }
  
  fwrite(results,
    paste0("~/Git_Repos/metabolomics_retention_index_score/RI_Specific_Data/Minimum_N_Analysis/", metab, ".txt"),
    quote = F, row.names = F, sep = "\t"
  )

}
