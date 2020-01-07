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

############################################################
# Top 50% clones in TIL Infusion compared to 4 week sample #
############################################################

# Calculating the amount

# @Param patient: specific patient code
# @Param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @Param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes

Top50_TILto4W <- function(patient, sampcohort, chain, clnefrc){
    # Load data including 4W data, TIL data, and Baseline data
    Load_data(patient, sampcohort, chain, clnefrc)
    # Creating a dataframe which contains TIL data and the associated clone fraction in the 4W sample
        Total_TIL_data <- TIL_data
        FW_clonefraction <- NA
        Total_TIL_data <- cbind(Total_TIL_data, FW_clonefraction)
    # Takes the top 50% (by clone fraction) of the Total_TIL_data dataset
    Top_50_clnefrc <- 0
    f <- 1
    for(i in Total_TIL_data$cloneFraction){
        if(Top_50_clnefrc < 0.5){
            Top_50_clnefrc <- Top_50_clnefrc + i
            Top50_TIL_data <- rbind(Top50_TIL_data, Total_TIL_data[f,])
            f <- f + 1
        }
    }
    # Fills the FW_clonefraction row with the clone fraction from the 4W clonotype that matched the specific clonotype in the TIL Infusion Product
    Total_TIL_data$FW_clonefraction <- FW_data$cloneFraction[match(Total_TIL_data$aaSeqCDR3, FW_data$aaSeqCDR3)]
    Total_TIL_data[is.na(Total_TIL_data)] <- 0
    # Finds and prints the sum of the FW_clonefraction and cloneFraction
    FW_rel_clnefrc <- sum(Top50_TIL_data$FW_clonefraction)
    print(FW_rel_clnefrc)
    print(Top_50_clnefrc)
}

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

###################################
# Baseline TIL-clonotype counting #
###################################

# Calculating the amount of TIL-related clones in the baseline repertoire

# @Param patient: specific patient code
# @Param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @Param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes

Base_TILclone <- function(patient, sampcohort, chain, clnefrc){
    # Loading specific locus nad sampcohort 4W, TIL, and baseline samples for patient
    Load_data(patient, sampcohort, chain, clnefrc)
    # Filtering the base data by the clones in the TIL data
    Base_filtered_data <- Base_data[Base_data$aaSeqCDR3 %in% TIL_data$aaSeqCDR3,]
    # Take all the clonotypes in the baseline sample above the set clone fraction
    Base_filtered_data <- Base_filtered_data[Base_filtered_data$cloneFraction > clnefrc,]
    # Find the number of TIL-clonotypes in the baseline sample
    Base_TILclone_count <- nrow(Base_filtered_data)
    print(Base_TILclone_count)
}
