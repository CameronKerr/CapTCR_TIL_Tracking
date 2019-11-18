############################
# TIL Clonality Comparison #
############################

## Takes the clone tracking function that takes the top 10 clones in each of the files and tracks it over time. ##
## This specific program tracks the TIL Infusion products over time.                                            ##
# @param sampcohort: desired sample cohort, could be gDNA, cDNA, or cfDNA
# @param src: the path towards the files. Each folder inside should then start with 'CLONES_' followed by the chain and the patient ID
# @param chain: Desired chain to analyze, could be TRA, TRB, TRD, TRG
# @param samplist: List of the names of samples that are to be compared

TIL_clonality <- function(sampcohort, src, chain, samplist){
    # Clontrackfx function takes in the file paths, file list, and names and is able to create custom dataframes #
    # that are used later in order to create the clone tracking plot.                                            #
    clonckfx <- function(inpath, file_list, name){
      n_rows <- sapply(file_list, function(x)
        nrow(read.table(paste(inpath, x, sep = ""),
                        sep = "\t", header = T)))
      totalcdr3 <- unlist(sapply(file_list, function(x) read.table(paste(inpath, x, sep = ""),
                                                                   sep = "\t", header = T)[,33] ))
      rownames(clonTRAck) <- totalcdr3
      colnames(clonTRAck) <- file_list
      for(f in file_list){
        mixcrfle <- read.table(paste(inpath,f, sep = ""),
                               header = T, sep = "\t",
                               stringsAsFactors = F, check.names = F)
        for(i in 1:nrow(clonTRAck)){
          clonTRAck[i,colnames(clonTRAck)[match(f,colnames(clonTRAck))]] <-
            mixcrfle$cloneCount[match(rownames(clonTRAck)[i],
                                      mixcrfle$aaSeqCDR3)]
        }
      }
      clonTRAck <- clonTRAck[!duplicated(rownames(clonTRAck)),]
      assign(paste(name),clonTRAck, envir = .GlobalEnv)
    }

    # Creates a list of files present in the certain paths. Then reads them into the 'mixcrfle' and 'compldfle' #
    # variables to use later. Also creates the variables for taking the top 10 clones in the files and coloring #
    # them using the colors in the 'top.colors' variable.                                                       #
       
     for (b in samplist){
            dir <- paste(src, sampcohort, "/CLONES_", chain , b, "/", sep="")
            file_list <- list.files(path = dir, pattern="")

            readlist = list()
            i <- 1
            for(f in file_list){
              mixcrfle <- read.table(paste(dir,
                                       f,
                                       sep = ""),
                                 header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE,
                                 na.strings = c("", "NA"), nrows=10)
              f <- substr(f, 0, nchar(f))
              readlist[[i]] <- mixcrfle$cloneCount
              names(readlist)[i] <- f
              i <- i + 1
            }

            inpath <- dir

            flelst <- list.files(dir, pattern = paste("", "", sep = ""))
            flelst <- file_list
            i <- 1
            for (f in flelst){
                  mixcrfle <- read.table(paste(inpath, f, sep = ""),
                                     header = TRUE, sep = "\t",
                                     stringsAsFactors = FALSE,
                                     na.strings = c("", "NA"), nrows=10)
              if(nrow(mixcrfle) < 1){next()}
              if(i == 1){
                compldfle <- mixcrfle
                compldfle <- cbind(cloneno = row.names(compldfle), filename = f, compldfle)
                i <- i + 1
              }
              else{
                compldfle1 <- mixcrfle
                compldfle1 <- cbind(cloneno = row.names(compldfle1), filename = f, compldfle1)
                compldfle <- rbind(compldfle, compldfle1)
                rm(compldfle1)
              }
            }
            compldfle$filename <- as.character(compldfle$filename)
            compldfle$cloneno <- as.character(compldfle$cloneno)
            compldfle$filename <- substr(compldfle$filename, 0, nchar(compldfle$filename)-4)
            compldfle$cloneno <- formatC(compldfle$cloneno, width=nchar(max(as.numeric(compldfle$cloneno))), flag="0")

            compldfle <- compldfle[!is.na(compldfle$filename),]
            top.colors=c("#8f0108", "#bd0026", "#f21d1d", "#f05120", "#ed5a2d", 
                         "#e09a6e", "#e3bfa8", "#e8e1dc", "#f7f4f2", "#faeeed")
            clonenocol = c(top.colors, rep("#ffffff",
                                           nlevels(as.factor(compldfle$cloneno))-10))

            mixcrfle$clonePercentage <- mixcrfle$cloneFraction*100 
            compldfle$clonePercentage <- compldfle$cloneFraction*100

            compldfle <- subset(compldfle, filename=='TIL Infusion Product')
            for (file in compldfle$filename){
                compldfle$filename <- paste (b)
            }

            assign(b), compldfle)
        }
    }
allfile <- rbind(paste(samplist)

    # Creates the plots for the alpha and beta loci, and applies the 'clonenocol' variable created in the previous #
    # step which assign the top 10 clones to 'top.colors' and the rest to white. A custom theme is used again.     #
    barpt <- ggplot(data = allfile, aes(y = clonePercentage, x = filename, fill=cloneno)) +
          theme_light() +
          theme(plot.title = element_text(colour = "black", hjust = 0.5, face = "italic", size = 13)) +
          geom_bar(colour = "#D3D3D3", stat="identity", width = 0.8) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black")) +
          labs(y = paste("Percentage of Reads", "Per Unique Clonotype"),
               x = "none") +
          scale_fill_manual(values = clonenocol) +
          ggtitle("Alpha Locus") +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_text(size = 15, angle = 90, hjust = 1, color = "black"),
                legend.position= "none", legend.direction = "horizontal",
                text = element_text(size = 15))
                
    png(paste("Sample_Clonality_Comparison_", sampcohort ,".png",sep="",
        width = 1280,
        height = 620,
        units = "px"))
    print(barpt)
    dev.off()
}
