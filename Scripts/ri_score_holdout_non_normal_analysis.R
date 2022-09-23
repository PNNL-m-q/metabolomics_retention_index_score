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
#' 
#' @param InputFolder The folder with the CoreMS data 
#' 
#' @return A data.table
format_data <- function(InputFolder) {
  
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
  
  ############################
  ## 5. Get New Mean and SD ##
  ############################
  
  message("Calculating global statistics...")
  
  # Get the global mean and standard deviation
  Globals <- Pre_Holdout %>%
    filter(`Truth Annotation` == "True.Positive" & `Compound Name` %in% Compounds) %>%
    dplyr::select(c(`Compound Name`, `Retention Index`)) %>%
    group_by(`Compound Name`) %>%
    nest() %>%
    mutate(
      Mean = map(data, function(x) {mean(x$`Retention Index`)}) %>% unlist(),
      SD = map(data, function(x) {sd(x$`Retention Index`)}) %>% unlist(),
      Skew = map(data, function(x) {e1071::skewness(x$`Retention Index`)}) %>% unlist(),
      N = map(data, function(x) {length(x$`Retention Index`)}) %>% unlist()
    ) %>%
    dplyr::select(-data) %>%
    rbind(
      data.table(
        `Compound Name` = unique(Annotations$`Compound Name`)[unique(Annotations$`Compound Name`) %in% Compounds == FALSE],
        Mean = 0, SD = 3, Skew = 0
      )
    )
  
  # Append global mean and standard deviation, then calculate new retention index score
  Pre_Holdout <- Pre_Holdout %>%
    merge(Globals, by = "Compound Name") 
  
  return(Pre_Holdout)
  
}  

# Get rank statistics 
getRankStats <- function(DF) {
  
  DF %>% 
    group_by(`Sample Name`, `Peak Index`) %>% 
    mutate(
      `max(RI Score) - RI Score` = max(`Retention Index Score`) - `Retention Index Score`,
      `min(RI Score)` = min(`Retention Index Score`),
      `max(RI Score)` = max(`Retention Index Score`), 
      `sd(RI Score)` = sd(`Retention Index Score`)
    ) %>%
    ungroup()
  
}

##################################
## DEFINE TRELLISCOPE FUNCTIONS ##
##################################

SampleMetadata <- fread("~/Git_Repos/metabolomics_scoring_metrics/Data/Metadata/Sample_Metadata.csv")
SampleMetadata$`Sample Type`[SampleMetadata$`Sample Type` == "Standard Mixture"] <- "Standard"

sample_type_plot <- function(CompoundName, TruthAnnotation) {
  
  # Read in the simulation data
  SimData <- fread(file.path("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Logistic_Prob_Holdout/", paste0(CompoundName, ".txt")))
  
  if (nrow(SimData %>% filter(`Truth Annotation` == TruthAnnotation)) == 0) {
    return((ggplot() + theme_void() + ggtitle(paste("No", TruthAnnotation, "Detected"))))
  }
  
  # Get both for min and max ranges
  SummaryRankStats <- SimData %>%
    filter(`Truth Annotation` == TruthAnnotation) %>%
    dplyr::select(`Sample Name`, `Peak Index`, `Retention Index Rank`, `New RI Rank`) %>%
    group_by(`Sample Name`, `Peak Index`) %>%
    summarise(
      `Min Rank Change` = min(`Retention Index Rank` - `New RI Rank`, na.rm = T),
      `Median Rank Change` = median(`Retention Index Rank` - `New RI Rank`, na.rm = T),
      `Max Rank Change` = max(`Retention Index Rank` - `New RI Rank`, na.rm = T)
    ) %>% 
    unique() %>%
    merge(SampleMetadata[,1:2], by = "Sample Name") %>%
    arrange(`Min Rank Change`) %>%
    mutate(
      `Sample ID` = as.numeric(factor(paste(`Sample Name`, `Peak Index`), 
                                      levels = unique(paste(`Sample Name`, `Peak Index`))))
    )
  
  # Generate a specific color vector
  ColorVector <- c("CSF" = "forestgreen", 
                   "Fungi" = "steelblue", 
                   "Plasma" = "red",
                   "Soil Crust" = "brown", 
                   "Standard" = "black",
                   "Urine" = "darkorange")
  
  # Subset color vector 
  ColorVector <- ColorVector[names(ColorVector) %in% unique(SummaryRankStats$`Sample Type`)]
  
  # Set a min and max for the plot
  themin <-  min(SimData$`Retention Index Rank` - SimData$`New RI Rank`, na.rm = T)
  if (themin > 0) {themin <- -1}
  themax <- max(SimData$`Retention Index Rank` - SimData$`New RI Rank`, na.rm = T) + 5
  if (themax < 1) {themax <- 1}
  
  # Create the base ggplot 
  Plot <- ggplot(SummaryRankStats, aes(x = `Sample ID`, y = `Median Rank Change`, color = `Sample Type`)) +
    geom_point() + geom_segment(aes(x = `Sample ID`, xend = `Sample ID`, y = `Min Rank Change`, yend = `Max Rank Change`)) +
    geom_hline(yintercept = 0) + theme_bw() + ylab("Original Rank - New Rank") +
    scale_color_manual(values = ColorVector) + ylim(themin, themax) +
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(TruthAnnotation) + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  
  return(Plot)
  
}

sample_type_plot_TP <- function(x) {sample_type_plot(x, "True.Positive")}
sample_type_plot_TN <- function(x) {sample_type_plot(x, "True.Negative")}

make_cog_TN <- function(CompoundName, TruthAnnotation = "True.Negative") {
  
  # Read in simulation data
  SimData <- fread(file.path("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Logistic_Prob_Holdout/", paste0(CompoundName, ".txt"))) %>%
    filter(`Truth Annotation` == TruthAnnotation)
  
  # Create and return the cognostics 
  if (nrow(SimData) == 0) {
    cogs <- tibble(
      `Original TN Rank Median` = cog(0, "Original TN Rank Median"),
      `Local TN Rank Median` = cog(0, desc = "Local TN Rank Median"),
      `Original TN Top 1` = cog(0, "Original TN Top 1"),
      `Original TN Top 2` = cog(0, "Original TN Top 2"),
      `Original TN Top 5` = cog(0, "Original TN Top 5"),
      `Local TN Top 1` = cog(0, "Local TN Top 1"),
      `Local TN Top 2` = cog(0, "Local TN Top 2"),
      `Local TN Top 5` = cog(0, "Local TN Top 5"),
    )
  } else {
    
    # Define rank proportion function
    rankProp <- function(ranks, N) {
      round(length(ranks[ranks <= N]) / length(ranks), 4)
    }
    
    cogs <- tibble(
      `Original TN Rank Median` = cog(median(SimData$`Retention Index Rank`, na.rm = T), "Original TN Rank Median"),
      `Local TN Rank Median` = cog(median(SimData$`New RI Rank`, na.rm = T), desc = "Local TN Rank Median"),
      `Original TN Top 1` = cog(rankProp(SimData$`Retention Index Rank`, 1), "Original TN Top 1"),
      `Original TN Top 2` = cog(rankProp(SimData$`Retention Index Rank`, 2), "Original TN Top 2"),
      `Original TN Top 5` = cog(rankProp(SimData$`Retention Index Rank`, 5), "Original TN Top 5"),
      `Local TN Top 1` = cog(rankProp(SimData$`New RI Rank`, 1), "Local TN Top 1"),
      `Local TN Top 2` = cog(rankProp(SimData$`New RI Rank`, 2), "Local TN Top 2"),
      `Local TN Top 5` = cog(rankProp(SimData$`New RI Rank`, 5), "Local TN Top 5"),
    )
  }
  
  return(cogs)
}

make_cog_TP <- function(CompoundName, TruthAnnotation = "True.Positive") {
  
  # Read in simulation data
  SimData <- fread(file.path("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Logistic_Prob_Holdout/", paste0(CompoundName, ".txt"))) %>%
    filter(`Truth Annotation` == TruthAnnotation)
  
  # Create and return the cognostics 
  if (nrow(SimData) == 0) {
    cogs <- tibble(
      `Original TP Rank Median` = cog(0, "Original TP Rank Median"),
      `Local TP Rank Median` = cog(0, desc = "Local TP Rank Median"),
      `Original TP Top 1` = cog(0, "Original TP Top 1"),
      `Original TP Top 2` = cog(0, "Original TP Top 2"),
      `Original TP Top 5` = cog(0, "Original TP Top 5"),
      `Local TP Top 1` = cog(0, "Local TP Top 1"),
      `Local TP Top 2` = cog(0, "Local TP Top 2"),
      `Local TP Top 5` = cog(0, "Local TP Top 5"),
    )
  } else {
    
    # Define rank proportion function
    rankProp <- function(ranks, N) {
      round(length(ranks[ranks <= N]) / length(ranks), 4)
    }
    
    cogs <- tibble(
      `Original TP Rank Median` = cog(median(SimData$`Retention Index Rank`, na.rm = T), "Original TP Rank Median"),
      `Local TP Rank Median` = cog(median(SimData$`New RI Rank`), desc = "Local TP Rank Median"),
      `Original TP Top 1` = cog(rankProp(SimData$`Retention Index Rank`, 1), "Original TP Top 1"),
      `Original TP Top 2` = cog(rankProp(SimData$`Retention Index Rank`, 2), "Original TP Top 2"),
      `Original TP Top 5` = cog(rankProp(SimData$`Retention Index Rank`, 5), "Original TP Top 5"),
      `Local TP Top 1` = cog(rankProp(SimData$`New RI Rank`, 1), "Local TP Top 1"),
      `Local TP Top 2` = cog(rankProp(SimData$`New RI Rank`, 2), "Local TP Top 2"),
      `Local TP Top 5` = cog(rankProp(SimData$`New RI Rank`, 5), "Local TP Top 5"),
    )
  }
  
  return(cogs)
}

##############
## ANALYSES ##
##############

MQ_Holdout <- format_data("~/Git_Repos/metabolomics_scoring_metrics/Data/CoreMS_Identification/")
Parameters <- fread("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Parameter_Estimates.csv")
colnames(Parameters)[9] <- "Compound Name" 
Compounds <- MQ_Holdout[MQ_Holdout$Mean != 0, "Compound Name"] %>% unique()

# Standard normal score (using probability)
NormScore <- function(`Retention Index`, `Retention Index Ref`, Mean, SD) {
  
  # First, adjust the mean
  if (Mean == 0) {
    Mean <- `Retention Index Ref`
  } else {
    Mean <- `Retention Index Ref` + Mean
  }
  
  # Then determine which side of the probability distribution the value falls on,
  # and calculate the score.
  if (`Retention Index` < Mean) {
    return(pnorm(`Retention Index`, Mean, SD) * 2)
  } else if (`Retention Index` > Mean) {
    return((1 - pnorm(`Retention Index`, Mean, SD))*2)
  } else {return(1)}
  
}

###########
## GAMMA ##
###########

GammaScore <- function(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, GammaShape, GammaRate) {
  
  if (`Compound Name` %in% Compounds) {
    
    # Calculate gamma probability score if the compound is in our subset. Otherwise, 
    # calculate the normal probability score. 
    if (`Retention Index` < Mean) {
      return(pgamma(`Retention Index`, GammaShape, GammaRate)*2)
    } else if (`Retention Index` > Mean) {
      return((1 - pgamma(`Retention Index`, GammaShape, GammaRate))*2)
    } else {return(1)}
    
  } else {
    NormScore(`Retention Index`, `Retention Index Ref`, Mean, SD)
  }
  
} 

# Calculate the new RI score before holdout 
Gamma <- MQ_Holdout %>%
  merge(Parameters[,c("Compound Name", "Gamma Shape", "Gamma Rate")], by = "Compound Name", all.x = T) %>%
  mutate(
    `New RI Score` = pmap(list(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, `Gamma Shape`, `Gamma Rate`), GammaScore) %>% unlist()
  )
  
# Define a function to calculate new ranks for a metabolite 
holdout_test_gamma <- function(metabolite) {
  
  # Get the frequency 
  Freq <- Gamma[Gamma$`Compound Name` == metabolite, "N"][1] %>% unlist()
  
  # Pick a random number 1-10 to be the group that is held out
  group <- sample(1:10, 1)
  
  # Grab the local mean and standard deviation following hold out
  Dist <- Gamma %>%
    filter(`Compound Name` == metabolite & `Truth Annotation` == "True.Positive") %>%
    mutate(Group = sample((1:Freq %% 10) + 1)) %>%
    filter(Group != group) %>%
    dplyr::select(`Retention Index`) %>%
    unlist()
  LocalGammaParams <- tryCatch({fitdistr(Dist, "gamma")}, error = function(e) {return(NULL)})
  
  if (is.null(LocalGammaParams)) {
    
    FinalGamma <- Gamma
    FinalGamma[,c("Gamma Shape", "Gamma Rate", "New RI Score", "New RI Rank")] <- NA
    
  } else {
  
    # Pull parameters
    newShape <- LocalGammaParams$estimate[1]
    newRate <- LocalGammaParams$estimate[2]
    
    FinalGamma <- Gamma %>%
      mutate(`New RI Score` = pmap(
        list(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, `New RI Score`),
        function(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, `New RI Score`) {
          ifelse(`Compound Name` == metabolite, 
                 GammaScore(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, newShape, newRate), 
                 `New RI Score`)
        }) %>% unlist()
      ) %>%
      group_by(`Sample Name`, `Peak Index`) %>%
      mutate(`New RI Rank` = rank(-`New RI Score`)) %>%
      filter(`Truth Annotation` != "") %>%
      mutate(
        `Gamma Shape` = newShape,
        `Gamma Rate` = newRate
      ) %>%
      ungroup() %>%
      filter(`Compound Name` == metabolite)
  }
    
  # Return final holdout
  return(FinalGamma)
    
}

# Use parallel computing 
registerDoParallel(detectCores())
for (i in 1:length(Compounds)) {
  message(paste("On compound:", Compounds[i]))
  results <- foreach (j=1:50, .combine = rbind) %dopar% {
   holdout_test_gamma(Compounds[i])
  }
  results$Iteration <- rep(1:50, each = nrow(results)/50)
  fwrite(results, file.path("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Gamma_Prob_Holdout/", paste0(Compounds[i], ".txt")),
         quote = F, row.names = F, sep = "\t")
}

# Trelliscopes
data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TP),
    cog = map_cog(Compound, make_cog_TP)
  ) %>%
  trelliscope(path = "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/5-Retention-Index-Distribution-Adjustment/Gamma_Holdout", 
              name = "True Positive", nrow = 1, ncol = 1, thumb = T)
data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TN),
    cog = map_cog(Compound, make_cog_TN)
  ) %>%
  trelliscope(path = "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/5-Retention-Index-Distribution-Adjustment/Gamma_Holdout", 
              name = "True Negative", nrow = 1, ncol = 1, thumb = T)

##############
## LOGISTIC ##
##############

LogisticScore <- function(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, LogisticLocation, LogisticScale) {
  
  if (`Compound Name` %in% Compounds) {
    
    # Calculate logistic probability score if the compound is in our subset. Otherwise, 
    # calculate the normal probability score. 
    # Get parameter estimates 
    if (`Retention Index` < Mean) {
      return(plogis(`Retention Index`, LogisticLocation, LogisticScale)*2)
    } else if (`Retention Index` > Mean) {
      return((1 - plogis(`Retention Index`, LogisticLocation, LogisticScale))*2)
    } else {return(1)}
    
  } else {
    NormScore(`Retention Index`, `Retention Index Ref`, Mean, SD)
  }
  
} 

# Calculate the new RI score before holdout 
Logistic <- MQ_Holdout %>%
  merge(Parameters[,c("Compound Name", "Logistic Location", "Logistic Scale")], by = "Compound Name", all.x = T) %>%
  mutate(
    `New RI Score` = pmap(list(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, `Logistic Location`, `Logistic Scale`), LogisticScore) %>% unlist()
  )

# Define a function to calculate new ranks for a metabolite 
holdout_test_logistic <- function(metabolite) {
  
  # Get the frequency 
  Freq <- Logistic[Logistic$`Compound Name` == metabolite, "N"][1] %>% unlist()
  
  # Pick a random number 1-10 to be the group that is held out
  group <- sample(1:10, 1)
  
  # Grab the local mean and standard deviation following hold out
  Dist <- Logistic %>%
    filter(`Compound Name` == metabolite & `Truth Annotation` == "True.Positive") %>%
    mutate(Group = sample((1:Freq %% 10) + 1)) %>%
    filter(Group != group) %>%
    dplyr::select(`Retention Index`) %>%
    unlist()
  LocalLogisParams <- tryCatch({fitdistr(Dist, "logistic")}, error = function(e) {return(NULL)})
  
  if (is.null(LocalLogisParams)) {
    
    FinalLogis <- Logistic %>% filter(`Truth Annotation` %in% c("True.Positive", "True.Negative") &
                                      `Compound Name` == metabolite)
    FinalLogis$`Logistic Location` <- FinalLogis$`Logistic Scale` <- FinalLogis$`New RI Score` <- FinalLogis$`New RI Rank` <- NA
    
  } else {
    
    # Pull parameters
    newLocation <- LocalLogisParams$estimate[1]
    newScale <- LocalLogisParams$estimate[2]
    
    FinalLogis <- Logistic %>%
      mutate(`New RI Score` = pmap(
        list(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, `New RI Score`),
        function(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, `New RI Score`) {
          ifelse(`Compound Name` == metabolite, 
                 LogisticScore(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, newLocation, newScale), 
                 `New RI Score`)
        }) %>% unlist()
      ) %>%
      group_by(`Sample Name`, `Peak Index`) %>%
      mutate(`New RI Rank` = rank(-`New RI Score`)) %>%
      filter(`Truth Annotation` != "") %>%
      mutate(
        `Logistic Location` = newLocation,
        `Logistic Scale` = newScale
      ) %>%
      ungroup() %>%
      filter(`Compound Name` == metabolite)
  }
  
  # Return final holdout
  return(FinalLogis)
  
}

# Use parallel computing 
registerDoParallel(detectCores())
for (i in 1:length(Compounds)) {
  message(paste("On compound:", Compounds[i]))
  results <- foreach (j=1:50, .combine = rbind) %dopar% {
    holdout_test_logistic(Compounds[i])
  }
  results$Iteration <- rep(1:50, each = nrow(results)/50)
  fwrite(results, file.path("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Logistic_Prob_Holdout/", paste0(Compounds[i], ".txt")),
         quote = F, row.names = F, sep = "\t")
}

# Trelliscopes
data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TP)
  ) %>%
  trelliscope(path = "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/4-Retention-Index-Distribution-Adjustment-Rank/Logistic_Holdout", 
              name = "True Positive", nrow = 1, ncol = 1, thumb = T)
data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TN)
  ) %>%
  trelliscope(path = "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/4-Retention-Index-Distribution-Adjustment-Rank/Logistic_Holdout", 
              name = "True Negative", nrow = 1, ncol = 1, thumb = T)

################
## LOG NORMAL ##
################

LNormScore <- function(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, MeanLog, SDLog) {
  
  if (`Compound Name` %in% Compounds) {
    
    # Calculate gamma probability score if the compound is in our subset. Otherwise, 
    # calculate the normal probability score. 
    if (`Retention Index` < Mean) {
      return(plnorm(`Retention Index`, MeanLog, SDLog)*2)
    } else if (`Retention Index` > Mean) {
      return((1 - plnorm(`Retention Index`, MeanLog, SDLog))*2)
    } else {return(1)}
    
  } else {
    NormScore(`Retention Index`, `Retention Index Ref`, Mean, SD)
  }
  
} 

# Calculate the new RI score before holdout 
LogNormal <- MQ_Holdout %>%
  merge(Parameters[,c("Compound Name", "Log Normal Mean", "Log Normal SD")], by = "Compound Name", all.x = T) %>%
  mutate(
    `New RI Score` = pmap(list(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, `Log Normal Mean`, `Log Normal SD`), LNormScore) %>% unlist()
  )

# Define a function to calculate new ranks for a metabolite 
holdout_test_lnorm <- function(metabolite) {
  
  # Get the frequency 
  Freq <- LogNormal[LogNormal$`Compound Name` == metabolite, "N"][1] %>% unlist()
  
  # Pick a random number 1-10 to be the group that is held out
  group <- sample(1:10, 1)
  
  # Grab the local mean and standard deviation following hold out
  Dist <- LogNormal %>%
    filter(`Compound Name` == metabolite & `Truth Annotation` == "True.Positive") %>%
    mutate(Group = sample((1:Freq %% 10) + 1)) %>%
    filter(Group != group) %>%
    dplyr::select(`Retention Index`) %>%
    unlist()
  LocalLognormParams <- tryCatch({fitdistr(Dist, "lognormal")}, error = function(e) {return(NULL)})
  
  if (is.null(LocalLognormParams)) {
    
    FinalLNorm <- Logistic %>% filter(`Truth Annotation` %in% c("True.Positive", "True.Negative") &
                                        `Compound Name` == metabolite)
    FinalLNorm$`Log Normal Mean` <- FinalLNorm$`Log Normal Scale` <- FinalLNorm$`New RI Score` <- FinalLNorm$`New RI Rank` <- NA
    
  } else {
    
    # Pull parameters
    newLNormMean <- LocalLognormParams$estimate[1]
    newLNormSD <- LocalLognormParams$estimate[2]
    
    FinalLNorm <- LogNormal %>%
      mutate(`New RI Score` = pmap(
        list(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, `New RI Score`),
        function(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, `New RI Score`) {
          ifelse(`Compound Name` == metabolite, 
                 LNormScore(`Retention Index`, `Retention Index Ref`, Mean, SD, `Compound Name`, newLNormMean, newLNormSD), 
                 `New RI Score`)
        }) %>% unlist()
      ) %>%
      group_by(`Sample Name`, `Peak Index`) %>%
      mutate(`New RI Rank` = rank(-`New RI Score`)) %>%
      filter(`Truth Annotation` != "") %>%
      mutate(
        `Log Normal Mean` = newLNormMean,
        `Log Normal SD` = newLNormSD
      ) %>%
      ungroup() %>%
      filter(`Compound Name` == metabolite)
  }
  
  # Return final holdout
  return(FinalLNorm)
  
}

# Use parallel computing 
registerDoParallel(detectCores())
for (i in 1:length(Compounds)) {
  message(paste("On compound:", Compounds[i]))
  results <- foreach (j=1:50, .combine = rbind) %dopar% {
    holdout_test_lnorm(Compounds[i])
  }
  results$Iteration <- rep(1:50, each = nrow(results)/50)
  fwrite(results, file.path("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/LogNormal_Prob_Holdout/", paste0(Compounds[i], ".txt")),
         quote = F, row.names = F, sep = "\t")
}

# Trelliscopes
data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TP)
  ) %>%
  trelliscope(path = "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/4-Retention-Index-Distribution-Adjustment-Rank/Lognormal_Holdout", 
              name = "True Positive", nrow = 1, ncol = 1, thumb = T)
data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TN)
  ) %>%
  trelliscope(path = "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/4-Retention-Index-Distribution-Adjustment-Rank/Lognormal_Holdout", 
              name = "True Negative", nrow = 1, ncol = 1, thumb = T)
