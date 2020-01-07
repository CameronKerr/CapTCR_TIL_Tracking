###################
# Diversity Plots #
###################

## The 'samplediv' function finds the inverse simpson diversity index of the alpha locus samples, ##
## from the gDNA and cDNA samples and compares them in a line plot.                               ##

# @Param patient: specific patient code
# @Param sampcohort: Desired sample cohort, could be gDNA, cDNA, or cfDNA
# @Param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG

library(immunarch)

sampleDiv <- function(patient, sampcohort, chain){
    pngfilename <- paste(patient, "sampleDiv.png")
    order <- eval(as.name(paste(patient, sampcohort, sep="")))
    # Alpha Locus Data Load using the 'repLoad' function from immunarch. #
    gDNA_Data = repLoad(paste("PATH/gDNA/CLONES_TRA", patient, "/", sep=""))
    cDNA_Data = repLoad(paste("/PATH/cDNA/CLONES_TRA", patient, "/", sep=""))
    gDNA_Data$data <- gDNA_Data
    cDNA_Data$data <- cDNA_Data
    
    # Takes the data generated above and runs it through immunarch's 'repDiversity' function with the #
    # inverse simpson setting.                                                                        #
    gDNA_div = repDiversity(gDNA_Data$data, .method = "inv.simp")
    cDNA_div = repDiversity(cDNA_Data$data, .method = "inv.simp")
    
    # Takes the diversity of the sample cohorts and generates a dataframe whith the inverse simpson  #
    # diversity is shown for both cohorts. If the value is NA, it will put an NA value in its place. #
    div <- data.frame(gDNA_div$Sample, NA, NA)
    colnames(div) <- c("Sample", "gDNA", "cDNA")
    for (i in gDNA_div$Sample) {
        div[match(i, div$Sample),]["gDNA"] = gDNA_div[match(i, gDNA_div$Sample),]["Value"]
    }
    for (i in cDNA_div$Sample) {
        div[match(i, div$Sample),]["cDNA"] = cDNA_div[match(i, cDNA_div$Sample),]["Value"]
    }
    target <- factor(c(paste(order)))
    div$Sample <- factor(div$Sample, levels = c(paste(order)))
    div <- div[match(target, div$Sample),]
    #Alpha Locus Plot
    f <- 1
    iterations <- sum(!is.na(div$cDNA))
    output <- matrix(ncol=1, nrow=iterations)
    for (i in which(!is.na(div$cDNA))){
        assign(paste("a", i, sep = ""), i)
        output[f] <- i
        f <- f + 1
    }
    
    # Takes the dataframes made in the last step and creates the arranged plots from them. It first creates the variable   #
    # 'iterations' which is the amount of times there is not an NA value in the cDNA column. It uses the variable in order #
    # to create a line plot where a line can be drawn through values with NA.                                              #
    iterations <- sum(!is.na(div$cDNA))
    div_plot <- ggplot(div, aes(x=div$Sample)) +
        theme_light() +
        theme(plot.title = element_text(colour = "black", face = "italic", hjust = 0.5, size = 15)) +
        geom_line(aes(y = div$cDNA, color="RNA", group=1), size = 1) +
        geom_point(aes(y = div$cDNA,), color = "orange2", size=3) +
        geom_line(aes(y = div$gDNA, color="DNA", group=1), size = 1) +
        geom_point(aes(y = div$gDNA), color = "orangered4", size=3) +
        scale_color_discrete(name = "Sample Type", labels=c("RNA", "RNA")) +
        scale_color_manual(values = c("orangered4", "orange2")) +
        labs(color = 'Sample Type') +
        geom_segment(aes(x = div$Sample[(output[1])], xend = div$Sample[(output[2])], y = div$cDNA[(output[1])], yend = div$cDNA[(output[2])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div$Sample[(output[2])], xend = div$Sample[(output[3])], y = div$cDNA[(output[2])], yend = div$cDNA[(output[3])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div$Sample[(output[3])], xend = div$Sample[(output[4])], y = div$cDNA[(output[3])], yend = div$cDNA[(output[4])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div$Sample[(output[4])], xend = div$Sample[(output[5])], y = div$cDNA[(output[4])], yend = div$cDNA[(output[5])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div$Sample[(output[5])], xend = div$Sample[(output[6])], y = div$cDNA[(output[5])], yend = div$cDNA[(output[6])]), linetype = "solid", color = "orange2", size = 1) +
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black")) +
        labs(y = paste("Inverse Simpson", "Diversity"),
               x = "none") +
        ggtitle(paste("Diversity Tracking", patient, sampcohort, chain, sep=" ")) +
        theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 15, angle = 90, hjust = 1, color = "black"),
            legend.position= "right", legend.direction = "vertical",
            text = element_text(size = 15))
    div_plot <- div_plot + scale_x_discrete(name = "", limits=c(paste(order)))
    # Creates and exports png file of the arranged alpha and beta loci #
    png(file = pngfilename,
        width = 1280,
        height = 620,
        units = "px")
    print(div_plot)
    dev.off()
 }
