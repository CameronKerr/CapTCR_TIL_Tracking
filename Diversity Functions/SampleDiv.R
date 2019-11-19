###################
# Diversity Plots #
###################

## The 'samplediv' function finds the inverse simpson diversity index of the alpha locus samples, ##
## from the gDNA and cDNA samples and compares them in a line plot.                               ##
# @param patient: specific patient code
# @param src: Path to MiXCR files. Must have the form of a folder with multiple sub-folders of each patient with the specific sample cohort starting with "CLONES_"

library(immunarch)
sampleDiv <- function(patient, src){
    pngfilename <- paste(patient, "sampleDiv.png")
    # Alpha Locus Data Load using the 'repLoad' function from immunarch. #
    alpha_gDNA_Data = repLoad(paste(src, "CLONES_TRA", patient, "/", sep=""))
    alpha_cDNA_Data = repLoad(paste(src, "CLONES_TRA", patient, "/", sep=""))
    alpha_gDNA_Data$data <- alpha_gDNA_Data
    alpha_cDNA_Data$data <- alpha_cDNA_Data
    # Beta Locus Data Load using the 'repLoad' function from immunarch. #
    beta_gDNA_Data = repLoad(paste(src, "CLONES_TRB", patient, "/", sep=""))
    beta_cDNA_Data = repLoad(paste(src, "CLONES_TRB", patient, "/", sep=""))
    beta_gDNA_Data$data <- beta_gDNA_Data
    beta_cDNA_Data$data <- beta_cDNA_Data
    
    
    # Takes the data generated above and runs it through immunarch's 'repDiversity' function with the #
    # inverse simpson setting.                                                                        #
    gDNA_div_alpha = repDiversity(alpha_gDNA_Data$data, .method = "inv.simp")
    cDNA_div_alpha = repDiversity(alpha_cDNA_Data$data, .method = "inv.simp")
    gDNA_div_beta = repDiversity(beta_gDNA_Data$data, .method = "inv.simp")
    cDNA_div_beta = repDiversity(beta_cDNA_Data$data, .method = "inv.simp")
    
    # Takes the diversity of the sample cohorts and generates a dataframe whith the inverse simpson  #
    # diversity is shown for both cohorts. If the value is NA, it will put an NA value in its place. #
    # Alpha Locus #
    div_alpha <- data.frame(gDNA_div_alpha$Sample, NA, NA)
    colnames(div_alpha) <- c("Sample", "gDNA", "cDNA")
    for (i in gDNA_div_alpha$Sample) {
        div_alpha[match(i, div_alpha$Sample),]["gDNA"] = gDNA_div_alpha[match(i, gDNA_div_alpha$Sample),]["Value"]
    }
    for (i in cDNA_div_alpha$Sample) {
        div_alpha[match(i, div_alpha$Sample),]["cDNA"] = cDNA_div_alpha[match(i, cDNA_div_alpha$Sample),]["Value"]
    }
    target <- factor(c(paste(order)))
    div_alpha$Sample <- factor(div_alpha$Sample, levels = c(paste(order)))
    div_alpha <- div_alpha[match(target, div_alpha$Sample),]
    # Beta Locus #
    div_beta <- data.frame(gDNA_div_beta$Sample, NA, NA)
    colnames(div_beta) <- c("Sample", "gDNA", "cDNA")
    for (i in gDNA_div_beta$Sample) {
        div_beta[match(i, div_beta$Sample),]["gDNA"] = gDNA_div_beta[match(i, gDNA_div_beta$Sample),]["Value"]
    }
    for (i in cDNA_div_beta$Sample) {
        div_beta[match(i, div_beta$Sample),]["cDNA"] = cDNA_div_beta[match(i, cDNA_div_beta$Sample),]["Value"]
    }
    target <- factor(c(paste(order)))
    div_beta$Sample <- factor(div_beta$Sample, levels = c(paste(order)))
    div_beta <- div_beta[match(target, div_beta$Sample),]
    
    #Alpha Locus Plot
    f <- 1
    iterations <- sum(!is.na(div_alpha$cDNA))
    output <- matrix(ncol=1, nrow=iterations)
    for (i in which(!is.na(div_alpha$cDNA))){
        assign(paste("a", i, sep = ""), i)
        output[f] <- i
        f <- f + 1
    }
    
    # Takes the dataframes made in the last step and creates the arranged plots from them. It first creates the variable   #
    # 'iterations' which is the amount of times there is not an NA value in the cDNA column. It uses the variable in order #
    # to create a line plot where a line can be drawn through values with NA.                                              #
   
   # Alpha Locus Code #
    iterations <- sum(!is.na(div_alpha$cDNA))
    div_plot_alpha <- ggplot(div_alpha, aes(x=div_alpha$Sample)) +
        theme_light() +
        theme(plot.title = element_text(colour = "black", face = "italic", hjust = 0.5, size = 15)) +
        geom_line(aes(y = div_alpha$cDNA, color="RNA", group=1), size = 1) +
        geom_point(aes(y = div_alpha$cDNA,), color = "orange2", size=3) +
        geom_line(aes(y = div_alpha$gDNA, color="DNA", group=1), size = 1) +
        geom_point(aes(y = div_alpha$gDNA), color = "orangered4", size=3) +
        scale_color_discrete(name = "Sample Type", labels=c("RNA", "RNA")) +
        scale_color_manual(values = c("orangered4", "orange2")) +
        labs(color = 'Sample Type') +
        geom_segment(aes(x = div_alpha$Sample[(output[1])], xend = div_alpha$Sample[(output[2])], y = div_alpha$cDNA[(output[1])], yend = div_alpha$cDNA[(output[2])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div_alpha$Sample[(output[2])], xend = div_alpha$Sample[(output[3])], y = div_alpha$cDNA[(output[2])], yend = div_alpha$cDNA[(output[3])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div_alpha$Sample[(output[3])], xend = div_alpha$Sample[(output[4])], y = div_alpha$cDNA[(output[3])], yend = div_alpha$cDNA[(output[4])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div_alpha$Sample[(output[4])], xend = div_alpha$Sample[(output[5])], y = div_alpha$cDNA[(output[4])], yend = div_alpha$cDNA[(output[5])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div_alpha$Sample[(output[5])], xend = div_alpha$Sample[(output[6])], y = div_alpha$cDNA[(output[5])], yend = div_alpha$cDNA[(output[6])]), linetype = "solid", color = "orange2", size = 1) +
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black")) +
        labs(y = paste("Inverse Simpson", "Diversity"),
               x = "none") +
        ggtitle("Alpha Locus") +
        theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 15, angle = 90, hjust = 1, color = "black"),
            legend.position= "right", legend.direction = "vertical",
            text = element_text(size = 15))

    # Beta Locus Code #
    f <- 1
    iterations <- sum(!is.na(div_beta$cDNA))
    output <- matrix(ncol=1, nrow=iterations)
    for (i in which(!is.na(div_beta$cDNA))){
      assign(paste("a", i, sep = ""), i)
      output[f] <- i
      f <- f + 1
    }
   
   div_plot_beta <- ggplot(div_beta, aes(x=div_beta$Sample)) +
        theme_light() +
        theme(plot.title = element_text(colour = "black", face = "italic", hjust = 0.5, size = 15)) +
        geom_line(aes(y = div_beta$cDNA, color="RNA", group=1), size = 1) +
        geom_point(aes(y = div_beta$cDNA,), color = "orange2", size=3) +
        geom_line(aes(y = div_beta$gDNA, color="DNA", group=1), size = 1) +
        geom_point(aes(y = div_beta$gDNA), color = "orangered4", size=3) +
        scale_color_discrete(name = "Sample Type", labels=c("RNA", "RNA")) +
        scale_color_manual(values = c("orangered4", "orange2")) +
        labs(color = 'Sample Type') +
        geom_segment(aes(x = div_beta$Sample[(output[1])], xend = div_beta$Sample[(output[2])], y = div_beta$cDNA[(output[1])], yend = div_beta$cDNA[(output[2])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div_beta$Sample[(output[2])], xend = div_beta$Sample[(output[3])], y = div_beta$cDNA[(output[2])], yend = div_beta$cDNA[(output[3])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div_beta$Sample[(output[3])], xend = div_beta$Sample[(output[4])], y = div_beta$cDNA[(output[3])], yend = div_beta$cDNA[(output[4])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div_beta$Sample[(output[4])], xend = div_beta$Sample[(output[5])], y = div_beta$cDNA[(output[4])], yend = div_beta$cDNA[(output[5])]), linetype = "solid", color = "orange2", size = 1) +
        geom_segment(aes(x = div_beta$Sample[(output[5])], xend = div_beta$Sample[(output[6])], y = div_beta$cDNA[(output[5])], yend = div_beta$cDNA[(output[6])]), linetype = "solid", color = "orange2", size = 1) +
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black")) +
        labs(y = paste("Inverse Simpson", "Diversity"),
               x = "none") +
        ggtitle("Beta Locus") +
        theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 15, angle = 90, hjust = 1, color = "black"),
            legend.position= "right", legend.direction = "vertical",
            text = element_text(size = 15))

   
   # Creates and exports png file of the arranged alpha and beta loci #
    png(file = pngfilename,
        width = 1280,
        height = 620,
        units = "px")
    print(grid.arrange(div_plot_alpha, div_plot_beta, nrow=1))
    dev.off()
 } 
