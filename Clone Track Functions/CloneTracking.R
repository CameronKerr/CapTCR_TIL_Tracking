############################
# Clone Tracking Function #
###########################

require(ggalluvial)
library(randomcoloR)

# Tracking immune cell clonotypes present in the reference sample
# @param patient: specific patient code 
# @param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param clnefrc: cut-off from 0 to 1 to track and plot only a subset of clonotypes
# @param tracking: Name of the file which contains the desired tracking clonotypes
# @param src: Path to MiXCR files. Must have the form of a folder with multiple sub-folders of each patient with the specific sample cohort starting with "CLONES_"


clonetrack.fx <- function(patient, sampcohort, chain, clnefrc, tracking, src){
    samporder <- eval(as.name(paste(patient, sampcohort, sep="")))
    #Compile a big file with patient's mixcr files loaded in from the TRA loci
    #Files must contain the initial source but are seperated based on the sample cohort. Inside the specific sample cohort folder, a folder starting with 'CLONES_' should contain the mixcr output files
    
    i <- 1
    for (f in samporder){
    mixcrfle <- read.table(paste(src, sampcohort, "/CLONES_", chain, patient, "/", f,".txt", sep = ""), 
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
    
    # Subset to include only clonotypes with more than specified clonal fraction    
    CDR3_fraction <- CDR3_fraction[CDR3_fraction$cloneFraction > clnefrc,] 
    
    ## append the empty clonotypes after here. 
    CDR3_fraction <- CDR3_fraction[!(CDR3_fraction$cloneFraction == 0),]
    
    # Number of samples
    mysamples <- unique(CDR3_fraction$filename)
    
    # Getting the CDR3 sequences from the reference
    Ref_data <- data.frame()
    nonRef_data <- data.frame()
    f <- 1
    for(i in CDR3_fraction$filename){
        if(i == Tracking){
        Ref_data <- rbind(Ref_data, CDR3_fraction[f,])
        }
        else{
        nonRef_data <- rbind(nonRef_data, CDR3_fraction[f,])
        }
        f <- f + 1
    }
    
    
    # Pulls gDNA reference if the sample does not have an infusion
    if(length(ref_data$aaSeqCDR3) == 0){
        ref_data <- data.frame()
        nonref_data <- data.frame()
         ref_fraction <- read.table(paste(src, chain, patient, "/",tracking , sep = ""), 
                   header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
        ref_fraction <- ref_fraction[!duplicated(ref_fraction$aaSeqCDR3),]
        ref_fraction <- cbind(cloneno = row.names(ref_fraction), 
                     filename = paste(reference, sep="), 
                   ref_fraction)
        ref_fraction <- ref_fraction[, c("filename","aaSeqCDR3","cloneFraction", "cloneCount")]
        # Subset to include only clonotypes with more than specified clonal fraction    
        ref_fraction <- ref_fraction[ref_fraction$cloneFraction > clnefrc,] 
        ## append the empty clonotypes after here. 
        ref_fraction <- ref_fraction[!(ref_fraction$cloneFraction == 0),]
        f <- 1
        for(i in ref_fraction$filename){
            if(i == tracking){
                ref_data <- rbind(ref_data, ref_fraction[f,])
            }
            else{
                nonref_data <- rbind(nonref_data, ref_fraction[f,])
               }
            f <- f + 1
        }
        # Generating dataframe with added gDNA TIL Infusion product and updated y-axis order for plots
        nonref_data <- CDR3_fraction
        first <- 0
        for(i in CDR3_fraction$filename){
            if(i == CDR3_fraction$filename[1]){
            first <- first + 1
            }
        }
        last <- first
        for(i in CDR3_fraction$filename){
            if(i == CDR3_fraction$filename[first+1]){
            last <- last + 1
            }
        }
        CDR3_fraction_first <- CDR3_fraction[1:first,]
        CDR3_fraction_ordered <- rbind(CDR3_fraction_first, ref_data)
        CDR3_fraction_last <- CDR3_fraction[(first + 1):last,]
        CDR3_fraction <- rbind(CDR3_fraction_ordered, CDR3_fraction_last) 
        CDR3_fraction$filename <- factor(CDR3_fraction$filename, levels = c(paste(CDR3_fraction$filename[first]), tracking, paste(CDR3_fraction$filename[length(CDR3_fraction$filename)])))
    }
   
   #Assign colors to Reference clonotypes
    Ref_clns <- Ref_data$aaSeqCDR3
    nonRef_clns <- nonRef_data$aaSeqCDR3
    message("Total number of Reference clonotypes: ")     
    print(length(Ref_clns))
    Ref_data <- na.omit(Ref_data)
    ## Reference filtering 
        if(length(Ref_clns) > 150){
            message("Total number of Reference clonotypes > 50 ")
            Ref_data_ordered <- unique(Ref_data$aaSeqCDR3[order(Ref_data$cloneCount, decreasing = TRUE)])
            mycolors <- distinctColorPalette(150)
            
            mycolors <- c(mycolors, rep("white", length(Ref_data_ordered)-150),
                         rep("white",length(nonRef_clns)))
            names(mycolors) <- c(Ref_data_ordered, nonRef_clns)
        }else{
            mycolors <- distinctColorPalette(length(Ref_clns))
            mycolors <- c(mycolors, rep("white", length(nonRef_clns)))
            names(mycolors) <- c(Ref_clns, nonRef_clns)
            
        }
    message("these are what we color: ")  
    print(mycolors[mycolors != "white"]) 

#Create alluvial based on reference CDR3 data
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
    myp
    png(file = paste(patient, sampcohort, chain, clnefrc, tracking, ".png", sep="_"),
        width = 4000,
        height = 1500,
        res=300)
   print(myp)
   dev.off()
}
