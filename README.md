# CapTCR_TIL_Tracking

CapTCR TIL Tracking is a research project being run in the [Pugh Lab](http://pughlab.org/) at the Princess Margaret Cancer Research Centre. The project uses the innovative tool that the lab has developed called [CapTCR](https://www.ncbi.nlm.nih.gov/pubmed/30530777) to longitudinally track adoptive-cell transfer of Tumor-infiltrating lymphocytes (TILs) in sold-tumor cancer patients from peripheral blood samples. 

## Materials and Methods

Through this study, patients peripheral blood was taken and genomic DNA (gDNA), RNA or copy DNA (cDNA), and circuilating free DNA (cfDNA) was extracted. The sample cohorts were taken through the single or double-step capture process of CapTCR and then sequenced. Resulting files were then taken through the bioinformatics pipeline. The first step was to align and assemble the clones using [MiXCR 1.3.1](https://www.nature.com/articles/nmeth.3364) and then run them through quality control downsampling. 

After quality control, the samples were analyzed using a custom Bash/Python/R pipeline. 

### [Specific CDR3 Clone Tracking](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Clone%20Track%20Functions/CloneTracking.R) ###

To demonstrate accurate longitudinal tracking of specific clones, the clone tracking function was used to measure the relative size of each clone and track them over time. This function takes the aaCDR3 sequences from the reference (TIL Infusion Product) and track the abundance of them in the other samples and display the results as an alluvial. What this allows us to see is if an accurate response was seen in the patients and view certain characteristics of the infusion and baseline samples (richness, clonality, diversity, similarity) that may influence patient response. 

The clone tracking function will take in the patient code, specific sample cohort wanting to be analyzed (gDNA, cDNA, or cfDNA), chain/locus (TRA, TRB, TRD, TRG), and the clone fraction which will be used to specifying how much each clone has to take up in order to be incorporated into the analysis. The function will track the clones according to the specific reference. For this function, the reference is the clonotypes shown in the TIL infusion. The maximum clonotypes the function will track is 150. These 150 will each be assigned a specific color and the size will be tracked longitdinally throughout the samples based on the recurrence of its amino acid CDR3 sequence.

Here is an example of the clone tracking plot displaying the clonotypes >0.002 from the TRA chain of the gDNA samples in a specific patient:
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/CloneTrackingPlot.png)

### [Clonal Abundance Comparison](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Clone%20Track%20Functions/ClonalAbundanceTracking.R) ###

To demonstrate the relativity between patient's TIL Infusion clonality, the clonal abundance comparison function was used to view the top 10 clones in the TIL Infusion of each patient. This function takes in a list of samples that you want analyzed, the specific sample cohort (gDNA, cDNA, or cfDNA), chain/locus (TRA, TRB, TRD, TRG) and will create a stacked bar plot of the top 10 clones of those specific samples.

Here is an example of the clone tracking plot displaying the TIL Infusion products of all the patients in the research project for the TRA and TRB chains:
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/ClonalAbundanceInfusionComparison.JPG)

### [Sample Diversity Function](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Diversity%20Functions/SampleDiv.R) ###

To demonstrate the effect of diversity on patients 
