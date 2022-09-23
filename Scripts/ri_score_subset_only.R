set.seed(450)

library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)
library(foreach)
library(doParallel)
library(trelliscopejs)
library(ggplot2)

######################
## Define functions ##
######################

#' Return the CoreMS annotation data formatted appropriately for holdout analysis
#'   with only the subsetted metabolites returned
#' 
#' @param InputFolder The folder with the CoreMS data 
#' @param ParametersFile Path to the parameters file with the distributional parameters
#' 
#' @return A data.table
format_data_subset_only <- function(InputFolder, ParametersFile) {
  
  #####################
  ## 1. Compile Data ##
  #####################
  
  # Read in the verified and cleaned annotation data
  CoreMS_Files <- list.files(InputFolder, full.names = T)
  
  message("Pulling all annotations...")
  
  # Pull all compounds
  Annotations <- do.call(rbind, lapply(CoreMS_Files, function(x) {
    
    # Read data 
    MQ <- fread(x)
    
    # Return all compounds 
    MQ %>% 
      dplyr::select(`Sample Name`, `Peak Index`, `Retention Index Score`, `Truth Annotation`,
             `Retention Index`, `Retention Index Ref`, `Compound Name`) %>%
      return()
    
  }))
  
  ######################################################################################
  ## 2. Remove Blank Compounds, Recalculate Rank, Remove Bins without a True Positive ##
  ######################################################################################
  
  message("Removing blanks...")
  
  # Remove blank compounds, recalculate all ranks, and remove bins without a true positive
  Pre_Holdout <- Annotations %>% 
    filter(`Compound Name` != "") %>%
    group_by(`Sample Name`, `Peak Index`) %>%
    mutate(`Retention Index Rank` = rank(-`Retention Index Score`)) %>%
    filter("True.Positive" %in% `Truth Annotation`) %>%
    ungroup()
  
  ####################################
  ## 3. Define True Positive Subset ##
  ####################################
  
  message("Subsetting down to true positives with n >= 30...")
  
  # Get n>=30 compounds
  Counts <- Pre_Holdout %>%
    filter(`Truth Annotation` == "True.Positive") %>%
    dplyr::select(`Compound Name`) %>%
    unlist() %>% 
    table(dnn = "Compound") %>% 
    data.frame() %>% 
    filter(Freq >= 30 & Compound != "[PNNLMET0040] Impurity 001 [12.148]") 
  
  Compounds <- Counts %>% dplyr::select(Compound) %>% unlist() %>% as.character()
  
  # Get TN Count
  TN_Counts <- Pre_Holdout %>%
    filter(`Truth Annotation` == "True.Negative") %>%
    dplyr::select(`Compound Name`) %>%
    unlist() %>% 
    table(dnn = "Compound") %>% 
    data.frame() %>% 
    filter(Compound != "[PNNLMET0040] Impurity 001 [12.148]" & Compound %in% Compounds) 
  TN_Compounds <- TN_Counts %>% dplyr::select(Compound) %>% unlist() %>% as.character()

  
  ##################################################################
  ## 4. Filter down to subset compounds only and recalculate rank ##
  ##################################################################
  
  NewSub <- Pre_Holdout %>%
    group_by(`Sample Name`, `Peak Index`) %>%
    filter(`Compound Name` %in% Compounds & "True.Positive" %in% `Truth Annotation`) %>%
    mutate(`Retention Index Rank` = rank(-`Retention Index Score`)) 
  
  # Read parameters file
  Parameters <- fread(ParametersFile)
  colnames(Parameters)[9] <- "Compound Name"
  
  # Merge Parameters
  NewSub <- merge(NewSub, Parameters, by = "Compound Name")

  return(NewSub)
  
}  

###############
## PULL DATA ##
###############

# Pull subset data 
subset_data <- format_data_subset_only(
  "~/Git_Repos/metabolomics_scoring_metrics/Data/CoreMS_Identification/",
  "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Parameter_Estimates.csv"
)

rankProp <- function(ranks, N) {
  round(length(ranks[ranks <= N & !is.na(ranks)]) / length(ranks), 4)
}

# Generate list of compounds in subset
Compounds <- subset_data %>% dplyr::select(`Compound Name`) %>% unique() %>% unlist()

######################################
## DEFINE FLEXIBLE HOLDOUT FUNCTION ##
######################################

holdout_analysis <- function(Data,
                             ScoreFunction,
                             Param1,
                             Param2,
                             EstimationFunction,
                             Outpath, 
                             NIterations = 50) {
  
  # Define the holdout function
  holdout_fun <- function(metabolite) { 
  
    # Get the frequency 
    Freq <- Data[Data$`Compound Name` == metabolite & Data$`Truth Annotation` == "True.Positive",] %>% nrow()
    
    # Pick a random number 1-10 to be the group that is held out
    group <- sample(1:10, 1)
    
    # Grab the local mean and standard deviation following hold out
    Dist <- Data %>%
      filter(`Compound Name` == metabolite & `Truth Annotation` == "True.Positive") %>%
      mutate(Group = sample((1:Freq %% 10) + 1)) %>%
      filter(Group != group) %>%
      dplyr::select(`Retention Index`) %>%
      unlist()
    NewParams <- EstimationFunction(Dist)
    
    # Create the final holdout results
    FinalData <- Data
    
    # If no parameters, return NA
    if (is.null(NewParams)) {
      FinalData[,c(Param1, Param2, "New RI Score", "New RI Rank")] <- NA
    } else {
      
      # Pull parameters
      FinalData[FinalData$`Compound Name` == metabolite, Param1] <- NewParams$estimate[1]
      FinalData[FinalData$`Compound Name` == metabolite, Param2] <- NewParams$estimate[2]
      
      FinalData <- FinalData %>%
        mutate(`New RI Score` = pmap(list(`Retention Index`, `Retention Index Ref`, !!sym(Param1), !!sym(Param2)), 
                                          ScoreFunction) %>% unlist()) %>%
        group_by(`Sample Name`, `Peak Index`) %>%
        mutate(`New RI Rank` = rank(-`New RI Score`)) %>%
        filter(`Compound Name` == metabolite & `Truth Annotation` != "") %>% 
        dplyr::select(`Compound Name`, `Sample Name`, `Peak Index`, `Retention Index Score`,
                      `Truth Annotation`, `Retention Index`, `Retention Index Ref`,
                      `Retention Index Rank`, !!sym(Param1), !!sym(Param2),
                      `New RI Score`, `New RI Rank`)
      
    }
    
    return(FinalData)
    
  }
  
  browser()
  
  # Use parallel computing 
  registerDoParallel(detectCores())
  for (i in 1:length(Compounds)) {
    message(paste("On compound:", Compounds[i]))
    results <- foreach (j=1:NIterations, .combine = rbind) %dopar% {
      holdout_fun(Compounds[i])
    }
    results$Iteration <- rep(1:NIterations, each = nrow(results)/NIterations)
    fwrite(results, file.path(Outpath, paste0(Compounds[i], ".txt")),
           quote = F, row.names = F, sep = "\t")
  }

}

#######################
## ORIGINAL ADJUSTED ##
#######################

# SCORE FUNCTION: In this case, param 1 is mean and param 2 is sd
normal_adj_score <- function(`Retention Index`, `Retention Index Ref`, Mean, SD) {
  return(exp(-((Mean - `Retention Index`)^2) / (2 * SD^2)))
}

# ESTIMATION FUNCTION: normal 
estimate_normal <- function(Dist) {
  tryCatch({MASS::fitdistr(Dist, "normal")}, error = function(e) {return(NULL)})
}

holdout_analysis(
  subset_data, normal_adj_score, "Normal Mean", "Normal SD", estimate_normal, 
  "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Subset_Holdout/NormAdjHoldout/"
)

########################
## NORMAL PROBABILITY ##
########################

# SCORE FUNCTION: In this case, param 1 is mean and param 2 is sd
NormScore <- function(`Retention Index`, `Retention Index Ref`, Mean, SD) {
  ScoreVal <- pnorm(`Retention Index`, Mean, SD)
  if (ScoreVal > 0.5) {ScoreVal <- 1 - ScoreVal}
  return(round(ScoreVal * 2, 8))
}

holdout_analysis(
  subset_data, NormScore, "Normal Mean", "Normal SD", estimate_normal, 
  "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Subset_Holdout/NormProbHoldout/"
)

########################
## GAMMA PROBABILITY ##
########################

# SCORE FUNCTION: In this case, param 1 is gamma shape, and param 2 is gamma rate
GammaScore <- function(`Retention Index`, `Retention Index Ref`, GammaShape, GammaRate) {
  ScoreVal <- pgamma(`Retention Index`, GammaShape, GammaRate)
  if (ScoreVal > 0.5) {ScoreVal <- 1 - ScoreVal}
  return(round(ScoreVal * 2, 8))
}

# ESTIMATION FUNCTION: gamma
estimate_gamma <- function(Dist) {
  tryCatch({MASS::fitdistr(Dist, "gamma")}, error = function(e) {return(NULL)})
}

holdout_analysis(
  subset_data, GammaScore, "Gamma Shape", "Gamma Rate", estimate_gamma, 
  "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Subset_Holdout/GammaProbHoldout/"
)

##########################
## LOGISTIC PROBABILITY ##
##########################

# SCORE FUNCTION: In this case, param 1 is logistic location, and param 2 is logistic scale
LogisticScore <- function(`Retention Index`, `Retention Index Ref`, LogisticLocation, LogisticScale) {
  ScoreVal <- plogis(`Retention Index`, LogisticLocation, LogisticScale)
  if (is.na(ScoreVal) || is.null(ScoreVal)) {
    return(exp(-((`Retention Index Ref` - `Retention Index`)^2) / (2 * 3^2)))
  }
  if (ScoreVal > 0.5) {ScoreVal <- 1 - ScoreVal}
  return(round(ScoreVal * 2, 8))
}

# ESTIMATION FUNCTION: logistic
estimate_logistic <- function(Dist) {
  tryCatch({MASS::fitdistr(Dist, "logistic")}, error = function(e) {return(NULL)})
}

holdout_analysis(
  subset_data, LogisticScore, "Logistic Location", "Logistic Scale", estimate_logistic, 
  "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Subset_Holdout/LogisticProbHoldout/"
)

############################
## LOG NORMAL PROBABILITY ##
############################

# SCORE FUNCTION: In this case, param 1 is log normal mean, and param 2 is the log normal sd
LogNormalScore <- function(`Retention Index`, `Retention Index Ref`, LogNormalMean, LogNormalSD) {
  ScoreVal <- plnorm(`Retention Index`, LogNormalMean, LogNormalSD)
  if (ScoreVal > 0.5) {ScoreVal <- 1 - ScoreVal}
  return(round(ScoreVal * 2, 8))
}

# ESTIMATION FUNCTION: log normal
estimate_lognormal <- function(Dist) {
  tryCatch({MASS::fitdistr(Dist, "lognormal")}, error = function(e) {return(NULL)})
}

holdout_analysis(
  subset_data, LogNormalScore, "Log Normal Mean", "Log Normal SD", estimate_lognormal, 
  "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Subset_Holdout/LogNormProbHoldout/"
)

