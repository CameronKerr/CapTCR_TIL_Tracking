######################
# Repertoire Overlap #
######################

## The 'repLap' function takes each sample and computes the number of overlapping CDR3 sequences between samples. ##                                                                             ##
# @param patient: specific patient code
# @param sampcohort: desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param src: Path to MiXCR-output files

library(immunarch)
repLap <- function(patient, sampcohort, src){
    pngfilename <- paste(patient, sampcohort, "repOverlap.png", sep="")
    i <- eval(as.name(paste(patient, sampcohort, sep="")))
    # Takes in the patient name and loads the data using the 'repLoad' function from immunarch. #
    # Immunarch's 'repOverlap' function is the used to compute the overlapping CDR3 sequences.  #
    # between each of the samples in each loci.                                             
    TRA_TILdata <- repLoad(paste(src, sampcohort, "/CLONES_TRA", patient,"/", sep=""))
    TRB_TILdata <- repLoad(paste(src, sampcohort, "/CLONES_TRB", patient, "/", sep=""))
    TRA_TIL_ov <- repOverlap(TRA_TILdata, .method="public", .verbose=F)
    TRA_myp <- vis(TRA_TIL_ov)
    TRA_myp <- TRA_myp + labs(y = "Patient TIL Infusion", x = "Patient TIL Infusion") + 
        theme(plot.title = element_text(colour = "black", hjust = 0.5, face = "italic", size = 15), axis.text = element_text(colour = "black", size = 10),
                axis.title = element_text(colour="black", hjust = 0.5, size = 15),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90, hjust = 1, size=13),) +
        labs(title = "Alpha Locus", hjust = 0.05)
    TRA_myp <- TRA_myp + scale_x_discrete(limits = c(i))
    TRA_myp <- TRA_myp + scale_y_discrete(limits = c(i))
    TRB_TIL_ov <- repOverlap(TRB_TILdata, .method="public", .verbose=F)
    TRB_myp <- vis(TRB_TIL_ov)
    TRB_myp <- TRB_myp + labs(y = "Patient TIL Infusion", x = "Patient TIL Infusion") + 
        theme(plot.title = element_text(colour = "black", hjust = 0.5, face = "italic", size = 15), axis.text = element_text(colour = "black", size = 10),
                axis.title = element_text(colour="black", hjust = 0.5, size = 15),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90, hjust = 1, size=13),) +
        labs(title = "Beta Locus", hjust = 0.05)
    TRB_myp <- TRB_myp + scale_x_discrete(limits = c(i))
    TRB_myp <- TRB_myp + scale_y_discrete(limits = c(i))
    png(file = pngfilename,
        width = 1280,
        height = 620,
        units = "px")
    print(grid.arrange(TRA_myp, TRB_myp, nrow=1))
    dev.off()
}
