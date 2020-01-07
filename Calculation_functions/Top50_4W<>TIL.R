###############################################################
# Top 50% clones in 4W sample compared to TIL Infusion Produt #
###############################################################

# Calculating the amount of the TIL Infusion Product taken up by the clones in the top 50% of the 4W repertoire

# @Param patient: specific patient code
# @Param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @Param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes

Top50_4WtoTIL <- function(patient, sampcohort, chain, clnefrc){
    # Load data including 4W data, TIL data, and Baseline data
    Load_data(patient, sampcohort, chain, clnefrc)
    # Creating a dataframe which contains TIL data and the associated clone fraction in the 4W sample
        Total_FW_data <- FW_data
        TIL_clonefraction <- NA
        Total_FW_data <- cbind(Total_FW_data, TIL_clonefraction)
    # Takes the top 50% (by clone fraction) of the Total_TIL_data dataset
    Top_50_clnefrc <- 0
    f <- 1
    for(i in Total_FW_data$cloneFraction){
        if(Top_50_clnefrc < 0.5){
            Top_50_clnefrc <- Top_50_clnefrc + i
            Top50_FW_data <- rbind(Top50_FW_data, Total_FW_data[f,])
            f <- f + 1
        }
    }
    # Fills the FW_clonefraction row with the clone fraction from the 4W clonotype that matched the specific clonotype in the TIL Infusion Product
    Total_FW_data$TIL_clonefraction <- TIL_data$cloneFraction[match(Total_FW_data$aaSeqCDR3, TIL_data$aaSeqCDR3)]
    Total_FW_data[is.na(Total_FW_data)] <- 0
    # Finds and prints the sum of the FW_clonefraction and cloneFraction
    TIL_rel_clnefrc <- sum(Top50_FW_data$TIL_clonefraction)
    print(TIL_rel_clnefrc)
    print(Top_50_clnefrc)
}
