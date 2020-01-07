############################
# Clone Tracking Function #
###########################

require(ggalluvial)

# Tracking immune cell clonotypes present in TIL Infusion product
# @Param patient: specific patient code
# @Param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @Param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes

clonetrack.fx <- function(patient, sampcohort, chain, clnefrc, tracking){
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
    
    # Loading specific locus nad sampcohort 4W, TIL, and baseline samples for patient
    Load_data(patient, sampcohort, chain, clnefrc)
    
    # Pulls gDNA TIL Infusion product if the sample does not have an infusion
    if(length(TIL_data$aaSeqCDR3) == 0){
        samporder <- c(samporder[1], "TIL Infusion Product", samporder[2])
        TIL_data <- data.frame()
         TIL_data <- read.table(paste("/PATH/gDNA/CLONES_", chain, patient, "/","TIL Infusion Product.txt", sep = ""),
                   header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
        TIL_data <- TIL_data[!duplicated(TIL_data$aaSeqCDR3),]
        TIL_data <- cbind(cloneno = row.names(TIL_data),
                     filename = 'TIL Infusion Product',
                   TIL_data)
        TIL_data <- TIL_data[, c("filename","aaSeqCDR3","cloneFraction", "cloneCount")]
        
        # Subset to include only clonotypes with more than specified clonal fraction
        TIL_data <- TIL_data[TIL_data$cloneFraction > clnefrc,]
        # Append the empty clonotypes after here.
        TIL_data <- TIL_data[!(TIL_data$cloneFraction == 0),]
        # Create a dataframe for all the TIL Infusion Product clonotypes
        #TIL_data <- CDR3_fraction[which(CDR3_fraction$filename=='TIL Infusion Product'),]
        # Generating dataframe with added gDNA TIL Infusion product and updated y-axis order for plots
        CDR3_fraction <- rbind(Base_data, TIL_data, FW_data)
        CDR3_fraction$filename <- factor(CDR3_added_fraction$filename, levels = c(samporder))
    }
    
    #Assign colors to TIL Infusion clonotypes
    TIL_clns <- TIL_data$aaSeqCDR3
    nonTIL_clns <- CDR3_fraction$aaSeqCDR3
    message("Total number of TIL Infusion clonotypes: ")
    print(length(TIL_clns))
    TIL_data <- na.omit(TIL_data)
    
    ## TIL filtering
        if(length(TIL_clns) > 150){
            message("Total number of TIL Infusion clonotypes > 50 ")
            TIL_data_ordered <- unique(TIL_data$aaSeqCDR3[order(TIL_data$cloneCount, decreasing = TRUE)])
            mycolors <- distinctColorPalette(150)
            mycolors <- c(mycolors, rep("white", length(TIL_data_ordered)-150),
                         rep("white",length(nonTIL_clns)))
            names(mycolors) <- c(TIL_data_ordered, nonTIL_clns)
        }else{
            mycolors <- distinctColorPalette(length(TIL_clns))
            mycolors <- c(mycolors, rep("white", length(nonTIL_clns)))
            names(mycolors) <- c(TIL_clns, nonTIL_clns)
        }
    message("these are what we color: ")
    print(mycolors[mycolors != "white"])
    #CDR3_fraction <- rbind(CDR3_fraction, TIL_data)
    levels(CDR3_fraction$filename) <- c(samporder)
    
    p <-  ggplot(CDR3_fraction, aes(x = filename,
                              y = cloneFraction,
                              fill = aaSeqCDR3,
                              stratum = aaSeqCDR3,
                              alluvium = aaSeqCDR3))
        myp <- p + geom_alluvium(decreasing = FALSE) +
        geom_stratum(decreasing = FALSE, stat = "alluvium") +
        scale_fill_manual(breaks = names(mycolors[mycolors != "white"]),
                      values = mycolors) +
        theme(axis.title.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.position = "blank",
          plot.margin = unit(c(0.2,0,0,0),"cm")) +
        ggtitle(paste("Clone Tracking Plot -", patient, chain, clnefrc, tracking, sep=" "))
    png(file = paste(patient, sampcohort, chain, clnefrc, tracking, ".png", sep="_"),
        width = 4000,
        #height = 1500,
        res=300)
   print(myp)
   dev.off()
}
