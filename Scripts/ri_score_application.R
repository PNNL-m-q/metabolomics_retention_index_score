set.seed(450)

library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

#############################
## Pull all necessary data ##
#############################

# First, let's pull the estimated parameters 
Parameters <- fread("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Parameter_Estimates.csv")

# Let's pull the list of compounds 
Compounds <- Parameters$Compound

# Get the corems files
CoreMS_Files <- list.files("~/Git_Repos/metabolomics_scoring_metrics/Data/CoreMS_Identification/", full.names = T)

# Pull all data - use the format data modules from the holdout scripts
MQ_All <- do.call(rbind, lapply(CoreMS_Files, function(x) {
  
  # Read file
  Data <- fread(x) %>%
    dplyr::select(`Sample Name`, `Compound Name`, `Peak Index`, 
                  `Retention Index`, `Retention Index Ref`, `Truth Annotation`,
                  `Retention Index Score`)
  
  return(Data)
  
}))

# Fix truth annotations
MQ_All$`Truth Annotation` <- gsub(".", " ", MQ_All$`Truth Annotation`, fixed = T)
MQ_All$`Truth Annotation`[MQ_All$`Truth Annotation` == ""] <- "Unknown"
MQ_All$`Truth Annotation` <- factor(MQ_All$`Truth Annotation`, levels = c("True Positive", "True Negative", "Unknown"))

# Since we only adjust for cases where there is a true positive, let's filter 
# to just examples where the true positive is in our subset 
MQ_Sub <- MQ_All %>%
  filter(`Sample Name` != "") %>%
  group_by(`Sample Name`, `Peak Index`) %>%
  filter("True Positive" %in% `Truth Annotation` && 
         all(`Compound Name`[which(`Truth Annotation` == "True Positive")] %in% Compounds))

# Get global statistics
Globals <- MQ_Sub %>%
  filter(`Truth Annotation` == "True Positive") %>%
  ungroup() %>%
  group_by(`Compound Name`) %>%
  summarise(
    Mean = mean(`Retention Index`),
    SD = sd(`Retention Index`)
  )

# Add global statistics to MQ_Sub data 
MQ_Sub <- MQ_Sub %>%
  merge(Globals, all.x = T, by = "Compound Name") %>%
  mutate(
    Mean = ifelse(is.na(Mean), `Retention Index Ref`, Mean),
    SD = ifelse(is.na(SD), 3, SD)
  )

################################
## Add Retention Index Scores ##
################################

adj_ri <- function(RI, Mean, SD) {
  return(exp(-1*((RI - Mean)^2) / (2 * SD^2)))
}

NormScore <- function(RI, Mean, SD) {
  if (RI < Mean) {
    return(pnorm(RI, Mean, SD) * 2)
  } else if (RI > Mean) {
    return((1 - pnorm(RI, Mean, SD))*2)
  } else {return(1)}
}

GammaScore <- function(RI, Mean, SD, CN) {
  
  if (CN %in% Compounds) {
    
    GammaShape <- Parameters[Parameters$Compound == CN, "Gamma Shape"] %>% unlist()
    GammaRate <- Parameters[Parameters$Compound == CN, "Gamma Rate"] %>% unlist()
    if (RI < Mean) {
      return(pgamma(RI, GammaShape, GammaRate)*2)
    } else if (RI > Mean) {
      return((1 - pgamma(RI, GammaShape, GammaRate))*2)
    } else {return(1)}
    
  } else {
    NormScore(RI, Mean, SD)
  }
  
} 

LogisticScore <- function(RI, Mean, SD, CN) {
  
  if (CN %in% Compounds) {
    
    LogisticLocation <- Parameters[Parameters$Compound == CN, "Logistic Location"] %>% unlist()
    LogisticScale <- Parameters[Parameters$Compound == CN, "Logistic Scale"] %>% unlist()
    if (RI < Mean) {
      val <- plogis(RI, LogisticLocation, LogisticScale)*2
    } else if (RI > Mean) {
      val <- (1 - plogis(RI, LogisticLocation, LogisticScale))*2
    } else {val <- 1}
    if (!is.na(val) && val > 1) {val <- 1}
    
    return(val)
    
  } else {
    NormScore(RI, Mean, SD)
  }
  
} 

LogNormalScore <- function(RI, Mean, SD, CN) {
  
  if (CN %in% Compounds) {
    
    LogNormalMean <- Parameters[Parameters$Compound == CN, "Log Normal Mean"] %>% unlist()
    LogNormalSD <- Parameters[Parameters$Compound == CN, "Log Normal SD"] %>% unlist()
    if (RI < Mean) {
      val <- plnorm(RI, LogNormalMean, LogNormalSD)*2
    } else if (RI > Mean) {
      val <- (1 - plnorm(RI, LogNormalMean, LogNormalSD))*2
    } else {val <- 1}
    if (!is.na(val) && val > 1) {val <- 1}
    
    return(val)
    
  } else {
    NormScore(RI, Mean, SD)
  }
  
} 


# Add scores
MQ_Sub <- MQ_Sub %>%
  mutate(
    `Normal Adj` = pmap(list(`Retention Index`, Mean, SD), adj_ri) %>% unlist(),
    Normal = pmap(list(`Retention Index`, Mean, SD), NormScore) %>% unlist(),
    Gamma = pmap(list(`Retention Index`, Mean, SD, `Compound Name`), GammaScore) %>% unlist()
  )

MQ_Sub <- MQ_Sub %>%
  mutate(
    Logistic = pmap(list(`Retention Index`, Mean, SD, `Compound Name`), LogisticScore) %>% unlist(),
    `Log Normal` = pmap(list(`Retention Index`, Mean, SD, `Compound Name`), LogNormalScore) %>% unlist()
  )

#################
## Make ggplot ##
#################

SampleMetadata <- fread("~/Git_Repos/metabolomics_scoring_metrics/Data/Metadata/Sample_Metadata.csv")

MQ_Score <- MQ_Sub %>%
  merge(SampleMetadata[,1:2], by = "Sample Name") %>%
  mutate(`Sample Type` = gsub("Standard Mixture", "Standard", `Sample Type`)) %>%
  dplyr::select(`Truth Annotation`, `Retention Index Score`, `Normal Adj`, Normal,
                Gamma, Logistic, `Log Normal`, `Sample Type`) %>%
  rename(Original = `Retention Index Score`, `Original Adjusted` = `Normal Adj`) %>%
  pivot_longer(c(Original, `Original Adjusted`, Normal, Gamma, Logistic, `Log Normal`)) %>%
  rename(Score = value, `Scoring Method` = name) %>%
  mutate(
    `Scoring Method` = factor(`Scoring Method`, 
       levels = c("Original", "Original Adjusted", "Normal", "Gamma", "Logistic", "Log Normal")),
    Score = round(Score, 8)
  ) 


ggplot(MQ_Score, 
       aes(x = `Scoring Method`, y = Score, fill = `Truth Annotation`)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") + facet_wrap(.~`Sample Type`)  

ggplot(MQ_Score %>% filter(`Truth Annotation` == "True Positive"), 
       aes(x = `Scoring Method`, y = Score)) + 
  geom_boxplot(fill = "forestgreen") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") + facet_wrap(.~`Sample Type`) + 
  ggtitle("True Positives")


almost_log <- function(val) {
  if (is.na(val)) {return(NA)} else if (val == 0) {return(0)} else {return(log(val))}
}

ggplot(MQ_Score %>% filter(`Truth Annotation` == "True Negative"), 
       aes(x = `Scoring Method`, y = lapply(Score, almost_log) %>% unlist())) + 
  geom_boxplot(fill = "firebrick") + theme_bw() + ylab("Log Score") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") + facet_wrap(.~`Sample Type`)  + 
  ggtitle("True Negatives")

ggplot(MQ_Score %>% filter(`Truth Annotation` == "Unknown"), 
       aes(x = `Scoring Method`, y = lapply(Score, almost_log) %>% unlist())) + 
  geom_boxplot(fill = "steelblue") + theme_bw() + ylab("Log10 Score") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") + facet_wrap(.~`Sample Type`)  + 
  ggtitle("Distribution of Unknowns")


#############################
## Write tables of results ##
#############################

MQ_Score %>% 
  filter(`Scoring Method` == "Original") %>% 
  dplyr::select(`Truth Annotation`, Score) %>%
  fwrite("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Applied_Score/Original.csv", quote = F, row.names = F)

MQ_Score %>% 
  filter(`Scoring Method` == "Original Adjusted") %>% 
  dplyr::select(`Truth Annotation`, Score) %>%
  fwrite("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Applied_Score/Original_Adjusted.csv", quote = F, row.names = F)

MQ_Score %>% 
  filter(`Scoring Method` == "Normal") %>% 
  dplyr::select(`Truth Annotation`, Score) %>%
  fwrite("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Applied_Score/Normal.csv", quote = F, row.names = F)

MQ_Score %>% 
  filter(`Scoring Method` == "Gamma") %>% 
  dplyr::select(`Truth Annotation`, Score) %>%
  fwrite("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Applied_Score/Gamma.csv", quote = F, row.names = F)

MQ_Score %>% 
  filter(`Scoring Method` == "Logistic") %>% 
  dplyr::select(`Truth Annotation`, Score) %>%
  fwrite("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Applied_Score/Logistic.csv", quote = F, row.names = F)

MQ_Score %>% 
  filter(`Scoring Method` == "Log Normal") %>% 
  dplyr::select(`Truth Annotation`, Score) %>%
  fwrite("~/Git_Repos/metabolomics_scoring_metrics/Retention_Index/RI_Specific_Data/Applied_Score/Log_Normal.csv", quote = F, row.names = F)

