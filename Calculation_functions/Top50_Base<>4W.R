##########################################################
# Top 50% of TIL-filtered baseline compared to 4W sample #
##########################################################

# Calculating amount of the 4W repertoire taken up by the baseline repertoire after filtered by the TIL Infusion Product clones

# @Param patient: specific patient code
# @Param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @Param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes

Top50_BASEto4W <- function(patient, sampcohort, chain, clnefrc){
    # Loading specific locus nad sampcohort 4W, TIL, and baseline samples for patient
    Load_data(patient, sampcohort, chain, clnefrc)
    # Filtering the base data by the clones in the TIL data
    Base_filtered_data <- Base_data[Base_data$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]
    # Fills in the FW_clonefraction column with the clone fraction of the matching clonotype in the 4W sample
    FW_clonefraction <- NA
    Base_filtered_data <- cbind(Base_filtered_data, FW_clonefraction)
    Base_filtered_data$FW_clonefraction <- FW_data$cloneFraction[match(Base_filtered_data$aaSeqCDR3, FW_data$aaSeqCDR3)]
    # Converts all NA values to zeroes
    Total_Base_data[is.na(Total_Base_data)] <- 0
    # Takes the top 50% of the baseline sample
    Top_50_clnefrc <- Top50_Base_data <- data.frame()
    f <- 1
    Top_50_clnefrc <- 0
    for(i in Total_Base_data$cloneFraction){
        if(Top_50_clnefrc < 0.5){
            Top_50_clnefrc <- Top_50_clnefrc + i
            Top50_Base_data <- rbind(Top50_Base_data, Total_Base_data[f,])
            f <- f + 1
        }
    }
    FW_rel_clnefrc <- sum(Top50_Base_data$FW_clonefraction)
    print(FW_rel_clnefrc)
    print(Top_50_clnefrc)
}
