#########################
# Quality Control Plots #
#########################

# Takes in the tabular log dataset and creates a quality control plot which displays multiple data points including total sequencing reads, aligned reads, and more in a scatter plot format. 

## Setting up envioronment ##
# Loading required libraries ggplot and gridExtra. #
library(ggplot2)
library(gridExtra)

# Stops the plots from showing values in scientific notation. #
options(scipen = 999)

## Aligning Plots Function ##
# Creates function called 'align_plots1' which takes the input of multiple plots and y labels. The function then binds the multiple plots together using the build in 'ggplotGrob' function from ggplot

# @Param alignstatsfile: Log file from MiXCR output which displays all the align stats
# @Param assemblestatsfile: Log file from MiXCR output which displays all the assemble stats
# @Param samplelist: The order of samples desired for the resulting plot
# @Param plotname: The desired name for the plot
# @Param inpath: The path towards the alignstatsfile and assemblestatsfile
# @Param inpath: The desired path to put the plot

mixcrQC.fx <- function(alignstatsfile, assemblestatsfile, samplelist,
                       plotname, inpath, plotpath){
    align_plots1 <- function (...) {
    pl <- list(...)
    stopifnot(do.call(all, lapply(pl, inherits, "gg")))
    gl <- lapply(pl, ggplotGrob)
    bind2 <- function(x, y) gtable:::rbind_gtable(x, y, "first")
    combined <- Reduce(bind2, gl[-1], gl[[1]])
    wl <- lapply(gl, "[[", "widths")
    combined$widths <- do.call(grid::unit.pmax, wl)
    grid::grid.newpage()
    grid::grid.draw(combined)
  }
  
    ## Plot Theme ##
    # Creates the theme 'mytheme' which features no panel or coloured background. #
    mytheme <- theme(axis.title.y = element_text(size = 25),
                       axis.title.x = element_blank(),
                       axis.line = element_line(color = "black"),
                       axis.text = element_text(size = 22),
                       axis.text.x = element_text(angle = 90, hjust = 1, size = 15)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "transparent",colour = NA),
              legend.key = element_rect(fill = "white", colour = "white"),
              plot.margin = unit(c(0,0,0,0),"cm"))
     ## Reading Log Files from MiXCR ##
    # Reads the outputed align and assemble stats from the MiXCR software. #
      alignstats <- read.csv(file = paste(inpath,alignstatsfile, sep =""),
                             sep = ",", header = T)
      assemblestats <- read.csv(file = paste(inpath, assemblestatsfile, sep =""),
                                sep = ",", header = T)
                                
    ## Sample List Matching. #
    # Matches the sample names to the previously generated SampleID column in the align and assemble
    # files. NOTE: The SampleID column was added manually in order to ensure proper x axis ordering
    # and size of the labels.
      for(i in samplelist){
        alignstats$samplename[grepl(i, alignstats$SampleID)] <- i
      }
      alignstats <- alignstats[!is.na(alignstats$samplename),]
      for(j in samplelist){
        assemblestats$samplename[grepl(j, assemblestats$SampleID)] <- j
      }
      assemblestats <- assemblestats[!is.na(assemblestats$samplename),]
      
    # Ensures that the samplename is in the same order as the samplelist. #
      alignstats$samplename <- factor(alignstats$samplename, levels = samplelist)
      assemblestats$samplename <- factor(assemblestats$samplename, levels = samplelist)
      
    ## Plot Creation from Quality Control Data ##
    #'myplot_totalseq' creates the plot that uses the data 'Total.sequencing.reads' from the align stats file.                                                                                  
      myplot_totalseq <- ggplot(aes(x = samplename, y = Total.sequencing.reads),
                                data = alignstats) +
        geom_point(size = 5) +
        mytheme +
        theme(axis.text.x = element_blank()) +
        # This is to avoid printing min, max and median on y-axis
        scale_y_continuous(limits = c(0, max(alignstats$Total.sequencing.reads)))
                                                                               
      myplot_alignedseq <- ggplot(aes(x = samplename, y = Successfully.aligned.reads),
                                  data = alignstats) +
        geom_point(size = 5) +
        mytheme +
        theme(axis.text.x = element_blank())
        
      myplot_averagereads <- ggplot(aes(x = samplename, y = Average.number.of.reads.per.clonotype),
                                    data = assemblestats) +
        geom_point(size = 5) +
        mytheme +
        theme(axis.text.x = element_blank())
                                                                        #
      myplot_readsused <- ggplot(aes(x = samplename, y = Reads.used.in.clonotypes.before.clustering..percent.of.total),
                                 data = assemblestats) +
        geom_point(size = 5) +
        mytheme +
        theme(axis.text.x = element_blank())

      myplot_totalclon <- ggplot(aes(x = samplename, y = Final.clonotype.count),
                                 data = assemblestats) +
        geom_point(size = 5) +
        mytheme
        
    # Exports all plots as pdfs combined using the 'align_plots1' function created earlier 
      png(file = paste(plotpath, plotname, sep = ""),
          width = 12000,
          height = 7000,
          res = 300)
      align_plots1(myplot_totalseq + ylab("Total reads"),
                   myplot_alignedseq + ylab("Total aligned reads"),
                   myplot_readsused + ylab("Reads used \nin clonotypes"),
                   myplot_averagereads + ylab("Average \nreads/clonotype"),
                   myplot_totalclon + ylab("Total clonotypes"))
      dev.off()
    }
    
## Reading the align and assemble stats files ##
alignstats <- read.csv(file = "/PATH/align_stats.csv",
                      sep = ",", header = T)
assemblestats <- read.csv(file = "/PATH/assemble_stats.csv",
                      sep = ",", header = T)
myfilenames <- alignstats$SampleID

## Seperating Sample Cohorts ##
# The goal is to create a seperate pdf for each sample cohort (gDNA, cDNA, and cfDNA). This grabs the sample cohorts from the myfilenames list and stores them as a variable for each sample cohort.  
mycdnas <- myfilenames[grepl("cDNA", myfilenames)]
mygenomic <- myfilenames[grepl("gDNA", myfilenames)]
mycfdna <- myfilenames[grepl("cfDNA", myfilenames)]
    
    
## Final Creation of Plots ##
# Uses the mixcrQC.fx function with the align stats, assemble stats, sample cohort-specific file names #
# as the input. Outputs pdf files containing the stacked plots.                                        #
mixcrQC.fx("align_stats.csv", "assemble_stats.csv", mycdnas,
          "cDNA_QC.png","/PATH/", "/PATH/")
mixcrQC.fx("align_stats.csv", "assemble_stats.csv", mygenomic,
          "gDNA_QC.png","/PATH", "/PATH/")
