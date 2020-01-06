################
# Loading Data #
################

# Loading specific locus and sample cohort of patient
# @Param patient: specific patient code
# @Param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @Param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
Load_data <- function(patient, sampcohort, chain, clnefrc){
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
    #Compile a big file with patient's mixcr files loaded in from the TRA loci
    i <- 1
    for (f in samporder){
      mixcrfle <- read.table(paste("/PATH/", sampcohort, "/CLONES_", chain, patient, "/", f,".txt", sep = ""),
                     header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE,
                     na.strings = c("", "NA"))
      if(i == 1){
        compldfle <- mixcrfle[!duplicated(mixcrfle$aaSeqCDR3),]
        compldfle <- cbind(cloneno = row.names(compldfle),
                         filename = f,
                         compldfle)
        i <- i + 1
       }
      else{
        compldfle1 <- mixcrfle[!duplicated(mixcrfle$aaSeqCDR3),]
        compldfle1 <- cbind(cloneno = row.names(compldfle1),
                          filename = f,
                          compldfle1)
        compldfle <- rbind(compldfle, compldfle1)
        rm(compldfle1)
       }
    }
    
    # Subset df
    CDR3_fraction <- compldfle[, c("filename","aaSeqCDR3","cloneFraction", "cloneCount")]
    
    ## append the empty clonotypes after here.
    CDR3_fraction <<- CDR3_fraction[!(CDR3_fraction$cloneFraction == 0),]
    
    # Number of samples
    mysamples <- unique(CDR3_fraction$filename)
    
    # Creating custom dataframes for the TIL Infusion product, 4W sample, and baseline sample
    TIL_data <<- CDR3_fraction[which(CDR3_fraction$filename=='TIL Infusion Product'),]
    TIL_data <<- TIL_data[TIL_data$cloneFraction > clnefrc,]
    FW_data <<- dplyr::filter(CDR3_fraction, grepl('4W', filename))
    Base_data <<- CDR3_fraction[which(CDR3_fraction$filename==samporder[1]),]
   }
