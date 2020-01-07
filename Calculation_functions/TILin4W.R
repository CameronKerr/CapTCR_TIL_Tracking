##########################################
# Amount of TIL clones in 4W Calculation #
##########################################

# Calculating the amount of the 4W repertoire taken up by TIL-related clones

# @Param patient: specific patient code
# @Param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @Param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes

TIL_calc <- function(patient, sampcohort, chain, clnefrc){
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
    options(scipen = 999)
    #Loads patient, sampcohort, chain specific data
    Load_data(patient, sampcohort, chain, clnefrc)
    # Calculating the response of TIL-clones in the 4W sample
    FW_data <- FW_data[FW_data$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]
    response <- sum(FW_data$cloneFraction)
    print(response)
}
