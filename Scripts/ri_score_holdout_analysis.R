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
      select(`Sample Name`, `Peak Index`, `Retention Index Score`, `Truth Annotation`,
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
    select(`Compound Name`) %>%
    unlist() %>% 
    table(dnn = "Compound") %>% 
    data.frame() %>% 
    filter(Freq >= 30 & Compound != "[PNNLMET0040] Impurity 001 [12.148]") 
  
  Compounds <- Counts %>% select(Compound) %>% unlist() %>% as.character()
  
  # Get TN Count
  TN_Counts <- Pre_Holdout %>%
    filter(`Truth Annotation` == "True.Negative") %>%
    select(`Compound Name`) %>%
    unlist() %>% 
    table(dnn = "Compound") %>% 
    data.frame() %>% 
    filter(Compound != "[PNNLMET0040] Impurity 001 [12.148]" & Compound %in% Compounds) 
  TN_Compounds <- TN_Counts %>% select(Compound) %>% unlist() %>% as.character()
  
  ############################
  ## 5. Get New Mean and SD ##
  ############################
  
  message("Calculating global statistics...")
  
  # Get the global mean and standard deviation
  Globals <- Pre_Holdout %>%
    filter(`Truth Annotation` == "True.Positive" & `Compound Name` %in% Compounds) %>%
    select(c(`Compound Name`, `Retention Index`)) %>%
    group_by(`Compound Name`) %>%
    nest() %>%
    mutate(
      Mean = map(data, function(x) {mean(x$`Retention Index`)}) %>% unlist(),
      SD = map(data, function(x) {sd(x$`Retention Index`)}) %>% unlist(),
      Skew = map(data, function(x) {e1071::skewness(x$`Retention Index`)}) %>% unlist(),
      N = map(data, function(x) {length(x$`Retention Index`)}) %>% unlist()
    ) %>%
    select(-data) %>%
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

#' Perform a 10% holdout analysis 
#' 
#' @param Data The formatted data from the format_data function. Required.
#' @param ScoreFunction A function that calculates the score based on columns in data. Required.
#' @param Outpath Folder to pass files to. Required.
#' @param NIterations A number of iterations. Default is 50.
#' @param Center Whether the mean should be centered by retention index reference. Default is TRUE. 
#'
#' @return A data.table with rank changes for true positives and true negatives of 
#'    metabolites within the n>=30 subset. Each true negative must be from a bin
#'    with a true positive identification, which is calculated in format_data. 
ten_percent_holdout_analysis <- function(Data,
                                         ScoreFunction,
                                         Outpath, 
                                         NIterations = 50,
                                         Center = TRUE) {
  
  # Adjust the means 
  if (Center) {
    Data$Mean <- ifelse(Data$Mean != 0, Data$Mean - Data$`Retention Index Ref`, 0)
  }

  # Calculate the new RI score before holdout 
  Data <- Data %>%
    mutate(
      `New RI Score` = pmap(list(`Retention Index`, `Retention Index Ref`, Mean, SD), ScoreFunction) %>% unlist()
    )
  
  # Define a function to calculate new ranks for a metabolite 
  holdout_test <- function(metabolite) {
    
    # Get the frequency 
    Freq <- Data[Data$`Compound Name` == metabolite, "N"][1] %>% unlist()
    
    # Pick a random number 1-10 to be the group that is held out
    group <- sample(1:10, 1)
    
    # Grab the local mean and standard deviation following hold out
    LocalMeanSD <- Data %>%
      filter(`Compound Name` == metabolite & `Truth Annotation` == "True.Positive") %>%
      mutate(Group = sample((1:Freq %% 10) + 1)) %>%
      filter(Group != group) %>%
      group_by(`Compound Name`) 
    
    if (Center) {
      LocalMeanSD <- LocalMeanSD %>%
        summarise(Mean = mean(`Retention Index` - `Retention Index Ref`), 
                  SD = sd(`Retention Index` - `Retention Index Ref`)
      )
    } else {
      LocalMeanSD <- LocalMeanSD %>%
        summarise(Mean = mean(`Retention Index`), 
                  SD = sd(`Retention Index`)
        )
    }
    
    # Adjust the score for this metabolite based on the 10% holdout. Adjust this score for all
    # annotation cases. Filter down to just what we need here - True Positives, True Negatives, 
    # Compound Name, Sample Name, Peak Index, Truth Annotation, New RI Score, New RI Rank, Mean and SD.
    # Subset down to the metabolite of interest.
    newMean <- LocalMeanSD$Mean %>% unlist()
    newSD <- LocalMeanSD$SD %>% unlist()
    
    Final <- Data %>%
      mutate(`New RI Score` = pmap(
        list(`Compound Name`, `Retention Index`, `Retention Index Ref`, `New RI Score`),
        function(`Compound Name`, `Retention Index`, `Retention Index Ref`, `New RI Score`) {
          ifelse(`Compound Name` == metabolite, 
                 ScoreFunction(`Retention Index`, `Retention Index Ref`, newMean, newSD), 
                 `New RI Score`)
        }) %>% unlist()
      ) %>%
      group_by(`Sample Name`, `Peak Index`) %>%
      mutate(`New RI Rank` = rank(-`New RI Score`)) %>%
      filter(`Truth Annotation` != "") %>%
      mutate(
        Mean = LocalMeanSD$Mean,
        SD = LocalMeanSD$SD
      ) %>%
      ungroup() %>%
      filter(`Compound Name` == metabolite)
    
    # Return final holdout
    return(Final)
    
  }
  
  # Iterate through each metabolite and save results
  Compounds <- unique(Data$`Compound Name`[Data$Mean != 0])
  
  # Use parallel computing 
  registerDoParallel(detectCores())
  for (i in 1:length(Compounds)) {
    message(paste("On compound:", Compounds[i]))
    results <- foreach (j=1:NIterations, .combine = rbind) %dopar% {
      holdout_test(Compounds[i])
    }
    results$Iteration <- rep(1:NIterations, each = nrow(results)/NIterations)
    fwrite(results, file.path(Outpath, paste0(Compounds[i], ".txt")),
           quote = F, row.names = F, sep = "\t")
  }
  
}

##################################
## DEFINE TRELLISCOPE FUNCTIONS ##
##################################

SampleMetadata <- fread("~/Git_Repos/metabolomics_scoring_metrics/Data/Metadata/Sample_Metadata.csv")
SampleMetadata$`Sample Type`[SampleMetadata$`Sample Type` == "Standard Mixture"] <- "Standard"

sample_type_plot <- function(CompoundName, TruthAnnotation) {
  
  # Read in the simulation data
  SimData <- fread(file.path("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Normal_Prob_Holdout/", paste0(CompoundName, ".txt")))
  
  if (nrow(SimData %>% filter(`Truth Annotation` == TruthAnnotation)) == 0) {
    return((ggplot() + theme_void() + ggtitle(paste("No", TruthAnnotation, "Detected"))))
  }
  
  # Get both for min and max ranges
  SummaryRankStats <- SimData %>%
    filter(`Truth Annotation` == TruthAnnotation) %>%
    select(`Sample Name`, `Peak Index`, `Retention Index Rank`, `New RI Rank`) %>%
    group_by(`Sample Name`, `Peak Index`) %>%
    summarise(
      `Min Rank Change` = min(`Retention Index Rank` - `New RI Rank`),
      `Median Rank Change` = median(`Retention Index Rank` - `New RI Rank`),
      `Max Rank Change` = max(`Retention Index Rank` - `New RI Rank`)
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
  themin <-  min(SimData$`Retention Index Rank` - SimData$`New RI Rank`)
  if (themin > 0) {themin <- -1}
  themax <- max(SimData$`Retention Index Rank` - SimData$`New RI Rank`) + 5
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
  SimData <- fread(file.path("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Normal_Prob_Holdout/", paste0(CompoundName, ".txt"))) %>%
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
      `Original TN Rank Median` = cog(median(SimData$`Retention Index Rank`), "Original TN Rank Median"),
      `Local TN Rank Median` = cog(median(SimData$`New RI Rank`), desc = "Local TN Rank Median"),
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
  SimData <- fread(file.path("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Normal_Prob_Holdout/", paste0(CompoundName, ".txt"))) %>%
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
      `Original TP Rank Median` = cog(median(SimData$`Retention Index Rank`), "Original TP Rank Median"),
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

## Analysis for 3-Retention-Index-Normal-Adjustment

# First, pull the data 
MQ_Holdout <- format_data("~/Git_Repos/metabolomics_scoring_metrics/Data/CoreMS_Identification/")

# Write the new scoring function
normal_adj_score <- function(`Retention Index`, `Retention Index Ref`, Mean, SD) {
  RID <- `Retention Index` - `Retention Index Ref`
  return(exp(-((RID - Mean)^2) / (2 * SD^2)))
}

# Conduct the 10% holdout analysis
ten_percent_holdout_analysis(MQ_Holdout, normal_adj_score, "~/Desktop/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Normal_Adj_Holdout/")

# Get compounds in subset
Compounds <- unique(MQ_Holdout$`Compound Name`[MQ_Holdout$Mean != 0])

# Make the trellicope display
data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TP),
    cog = map_cog(Compound, make_cog_TP)
  ) %>%
  trelliscope(path = "~/Desktop/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/3-Retention-Index-Normal-Adjustment/Trelli1", 
              name = "True Positive", nrow = 1, ncol = 1, thumb = T)

data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TN),
    cog = map_cog(Compound, make_cog_TN)
  ) %>%
  trelliscope(path = "~/Desktop/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/3-Retention-Index-Normal-Adjustment/Trelli1", 
              name = "True Negative", nrow = 1, ncol = 1, thumb = T)

## Analysis for 4-Retention-Index-Probability-Test

# First, pull the data 
MQ_Holdout <- format_data("~/Git_Repos/metabolomics_scoring_metrics/Data/CoreMS_Identification/")

# Write the new scoring function
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

ten_percent_holdout_analysis(MQ_Holdout, NormScore, "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Normal_Prob_Holdout/")

# Get compounds in subset
Compounds <- unique(MQ_Holdout$`Compound Name`[MQ_Holdout$Mean != 0])

# Make the trellicope display
data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TP),
    cog = map_cog(Compound, make_cog_TP)
  ) %>%
  trelliscope(path = "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/4-Retention-Index-Probability-Test/Trelli1", 
              name = "True Positive", nrow = 1, ncol = 1, thumb = T)
data.frame(
  "Compound" = Compounds
) %>%
  mutate(
    panel = map_plot(Compound, sample_type_plot_TN),
    cog = map_cog(Compound, make_cog_TN)
  ) %>%
  trelliscope(path = "~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/Markdowns/4-Retention-Index-Probability-Test/Trelli1/", 
              name = "True Negative", nrow = 1, ncol = 1, thumb = T)

# Since the scores look similar, let's do a test case
NonProb <- fread("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Normal_Adj_Holdout/D-galactose [major].txt")
Prob <- fread("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Normal_Prob_Holdout/D-galactose [major].txt")

merge(
  NonProb %>% filter(`Truth Annotation` == "True.Positive" & Iteration == 25) %>% select(c("Sample Name", "New RI Score")),
  Prob %>% filter(`Truth Annotation` == "True.Positive" & Iteration == 25) %>% select(c("Sample Name", "New RI Score")),
  by = "Sample Name"
) %>%
  ggplot(aes(x = `New RI Score.x`, y = `New RI Score.y`)) + geom_point() +
  xlab("Normal Kernel Adjusted Score") + ylab("Normal Probability Adjusted Score")
