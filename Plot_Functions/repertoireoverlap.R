######################
# Repertoire Overlap #
######################

## The 'repLap' function takes each sample and computes the number of overlapping CDR3 sequences between samples

# @Param patient: specific patient code
# @Param sampcohort: desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG

library(immunarch)

repLap <- function(patient, sampcohort, chain){
    pngfilename <- paste(patient, sampcohort, "repOverlap.png", sep="")
    i <- eval(as.name(paste(patient, sampcohort, sep="")))
    
    # Takes in the patient name and loads the data using the 'repLoad' function from immunarch. #
    # Immunarch's 'repOverlap' function is the used to compute the overlapping CDR3 sequences.  #
    # between each of the samples in each loci.
    
    TLML_data <- repLoad(paste("/Users/cameronk/OneDrive - UHN/TLML/TLML_Clones/", sampcohort, "/CLONES_", chain, patient,"/", sep=""))
    TLML_ov <- repOverlap(TLML_data, .method="morisita", .verbose=F)
    TLML_myp <- vis(TLML_ov)
    TLML_myp <- TLML_myp + labs(y = "Patient TIL Infusion", x = "Patient TIL Infusion") +
        theme(plot.title = element_text(colour = "black", hjust = 0.5, face = "italic", size = 8), axis.text = element_text(colour = "black", size = 10),
                axis.title = element_text(colour="black", hjust = 0.5, size = 8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(angle = 90, hjust = 1, size=8),) +
        labs(title = "Alpha Locus", hjust = 0.05)
    TLML_myp <- TLML_myp + scale_x_discrete(limits = c(i))
    TLML_myp <<- TLML_myp + scale_y_discrete(limits = c(i))
    
    png(file = pngfilename,
        #width = 2280,
        #height = 1220,
        units = "px")
    print(TLML_myp)
    dev.off()
}
