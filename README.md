# CapTCR TIL Tracking

CapTCR TIL Tracking is a research project being run in the [Pugh Lab](http://pughlab.org/) at the Princess Margaret Cancer Research Centre. The project uses the innovative tool that the lab has developed called [CapTCR](https://www.ncbi.nlm.nih.gov/pubmed/30530777) to longitudinally track adoptive-cell transfer of Tumor-infiltrating lymphocytes (TILs) in sold-tumor cancer patients from peripheral blood samples. 

## Materials and Methods

Through this study, patients peripheral blood was taken and genomic DNA (gDNA), RNA or copy DNA (cDNA), and circuilating free DNA (cfDNA) was extracted. The sample cohorts were taken through the single or double-step capture process of CapTCR and then sequenced. Resulting files were then taken through the bioinformatics pipeline. The first step was to align and assemble the clones using [MiXCR 1.3.1](https://www.nature.com/articles/nmeth.3364) and then run them through quality control downsampling. 

After quality control, the samples were analyzed using a custom Bash/Python/R pipeline. 

### [Specific CDR3 Clone Tracking](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Clone%20Track%20Functions/CloneTracking.R) ###

To demonstrate accurate longitudinal tracking of specific clones, the clone tracking function was used to measure the relative size of each clone and track them over time. This function takes the aaCDR3 sequences from the reference (TIL Infusion Product) and track the abundance of them in the other samples and display the results as an alluvial. What this allows us to see is if an accurate response was seen in the patients and view certain characteristics of the infusion and baseline samples (richness, clonality, diversity, similarity) that may influence patient response. 

The clone tracking function will take in the patient code, specific sample cohort wanting to be analyzed (gDNA, cDNA, or cfDNA), chain/locus (TRA, TRB, TRD, TRG), and the clone fraction which will be used to specifying how much each clone has to take up in order to be incorporated into the analysis. The function will track the clones according to the specific reference. For this function, the reference is the clonotypes shown in the TIL infusion. The maximum clonotypes the function will track is 150. These 150 will each be assigned a specific color and the size will be tracked longitdinally throughout the samples based on the recurrence of its amino acid CDR3 sequence.

Here is an example of the clone tracking plot displaying the clonotypes >0.002 from the TRA chain of the gDNA samples in a specific patient (Patient and sample names have been cut out due to patient confidentiality):
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/CloneTrackingPlot.png)

### [Clonal Abundance Comparison](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Clone%20Track%20Functions/ClonalAbundanceTracking.R) ###

To demonstrate the relativity between patient's TIL Infusion clonality, the clonal abundance comparison function was used to view the top 10 clones in the TIL Infusion of each patient. This function takes in a list of samples that you want analyzed, the specific sample cohort (gDNA, cDNA, or cfDNA), chain/locus (TRA, TRB, TRD, TRG) and will create a stacked bar plot of the top 10 clones of those specific samples.

Here is an example of the clone tracking plot displaying the TIL Infusion products of all the patients in the research project for the TRA and TRB chains (Patient and sample names have been cut out due to patient confidentiality):
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/ClonalAbundanceInfusionComparison.JPG)

### [Sample Diversity Function](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Diversity%20Functions/SampleDiv.R) ###

To demonstrate the effect of diversity on patients and to compare sample type and loci measure of diversity the sample diversity function was used to compute the inverse-simpson diversity for each sample and compare them over time. This allows us to see a patients progression as a view of their whole T-cell repertoire as a measure of diversity. The inverse-simpson diversity was computed using the [Immunarch 4.0 package](https://github.com/immunomind/immunarch). 

The index allows us to measure and quantify the average abundance of types in the sample datasets and find the effective number of types. The sampleDiv function will produce a line plot with both the diversity computed for the DNA and RNA samples.

Here is an example of the diversity plot for a patient's TRA chain DNA and RNA samples (Patient and sample names have been cut out due to patient confidentiality):
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/SampleDiversityPlot.png)

### [Repertoire Overlap Function](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Repertoire%20Overlap/RepOverlap.R) ###

To demonstrate how the dynamics of a persons T-cell repertoire effects the patients response, the repLap function allowed us to analyze how much of a persons immune response is overlapping between multiple samples. This analysis was done with a specific focus on looking at Baseline <> TIL Infusion Product overlap, Baseline <> Post-Infusion overlap and TIL Infusion Product <> Post-Infusion overlap.

The repLap function will take in a patient's sample for a specific sample cohort and will compute the overlap between the samples. It produces a heatmap with the number of overlapping clonotypes shown between the samples. We used the Immunarch 4.0 'repOverlap' function to sort through the samples and compute the overlap. 

Here is an example of a patient's repertoire overlap for their DNA samples TRA and TRB chains (Patient and sample names have been cut out due to patient confidentiality):
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/RepertoireOverlap.png)

## Results ##

### Optimal Steps for CapTCR Molecular Depletion ###

### Tracking Immune Response in Peripheral Blood Samples ###

### Effect of Clonality on Immune Response ###

### Effect of Baseline <> Post-TIL Infusion Product Repertoire Overlap on Immune Response ###

### Selected Clonal Expansion of TIL Infusion ###
