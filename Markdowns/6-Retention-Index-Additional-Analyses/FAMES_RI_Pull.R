CoreMS_Files <- list.files("/Volumes/mqCC207720/3_CSV_Files/Reprocess_Sum_Abundance_02262021/", full.names = T)

# Read in annotations
Annotations <- fread("~/Downloads/AnnoTP.txt")

# Subset to N > 30
MetaboliteCount <- table(Annotations$`Compound Name`, dnn = c("Compound.Name")) %>%
  data.frame() 
Anno30 <- Annotations[Annotations$`Compound Name` %in% 
                        MetaboliteCount[MetaboliteCount$Freq >= 30, "Compound.Name"],] %>%
  filter(`Compound Name` != "[PNNLMET0040] Impurity 001 [12.148]")
Compounds <- unique(Anno30$`Compound Name`)


# Run the for loop 
for (file in CoreMS_Files) {
  
  message(file)
  
  basename <- strsplit(file, "/", fixed = T) %>% unlist() %>% tail(1) %>% gsub(pattern = ".csv", replacement = "_TP_FAMES.txt", fixed = T)
  
  # Read data
  data <- fread(file)
  
  # Subset to any relevant information
  AnnoSub <- Anno30 %>% filter(Anno30$`Sample Name` %in% unique(data$`Sample name`))
  
  if (nrow(AnnoSub) != 0) {
    
    # Pull only FAMEs information
    FAMES <- data %>% 
      
      # Only keep the FAME file and any files in the true positive subset
      filter(grepl("FAME", `Sample name`) & grepl("C8|C9|C10|C12|C14|C16|C18|C20|C22|C24|C26|C28", `Compound Name`) &
               `Retention index` == `Retention index Ref`) %>%
      
      # Select only columns of interest
      dplyr::select(`Sample name`, `Peak Index`, `Retention Time`, `Retention Time Ref`, 
                    `Retention index`, `Retention index Ref`, `Compound Name`) %>%
      
      # Rename 
      rename(`Sample Name` = `Sample name`, `Retention Index` = `Retention index`, `Retention Index Ref` = `Retention index Ref`)
    
    # Sub AnnoSub for merging
    rbind(
      FAMES, 
      AnnoSub %>% 
        dplyr::select("Sample Name", "Peak Index", "Retention Time", "Retention Time Ref", 
                      "Retention Index", "Retention Index Ref", "Compound Name")
    ) %>%
      fwrite(file.path("~/Downloads/FAMES_RI_Pull", basename), quote = F, row.names = F, sep = "\t")
    
  }
  
}





