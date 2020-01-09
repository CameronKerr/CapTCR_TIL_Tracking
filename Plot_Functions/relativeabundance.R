###################################
### Relative Abundance Function ###
###################################

# Calculating and plotting clonal relative abundance
# @Param Bhatt patient: specific patient code
# @Param Bhatt sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param Bhatt chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG

Rel_Abundance <- function(patient, sampcohort, chain){
    pngfilename <- paste(patient, sampcohort, chain, "relative_abundance.png", sep="")
    TIL_data <- repLoad(paste("/PATH/", sampcohort, "/CLONES_", chain, patient,"/", sep=""))
    TIL_relab <- repClonality(TIL_data, .method="homeo",
                             .clone.types=c(Small=0.0005, Medium=0.005, Large=0.05, Hyperexpanded=0.5))
    relab_plot <- vis(TIL_relab)
    relab_plot <- relab_plot + theme(plot.title = element_text(colour = "black", hjust = 0.5, face = "italic", size = 20), axis.text = element_text(colour = "black", size = 10),
                axis.title = element_text(colour="black", hjust = 0.5, size = 20),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(angle = 90, hjust = 1, size=15),) +
        labs(title = paste("Relative Abundance", patient, sampcohort, sep=" "), hjust = 0.05)
    relab_plot <- relab_plot + scale_x_discrete(limits = c(eval(as.name(paste(patient,sampcohort,sep="")))))
    png(file = pngfilename,
        width = 2280,
        height = 1220,
        units = "px")
    print(relab_plot)
    dev.off()
}
