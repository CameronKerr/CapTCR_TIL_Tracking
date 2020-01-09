# CapTCR TIL Tracking

CapTCR TIL Tracking is a research project being run in the [Pugh Lab](http://pughlab.org/) at the Princess Margaret Cancer Research Centre. The project uses the innovative tool that the lab has developed called [CapTCR](https://www.ncbi.nlm.nih.gov/pubmed/30530777) to longitudinally track adoptive-cell transfer of Tumor-infiltrating lymphocytes (TILs) in sold-tumor cancer patients from peripheral blood samples. 

## Materials and Methods

Through this study, patients peripheral blood was taken and genomic DNA (gDNA), RNA or copy DNA (cDNA), and circuilating free DNA (cfDNA) was extracted. The sample cohorts were taken through the single or double-step capture process of CapTCR and then sequenced. Resulting files were then taken through the bioinformatics pipeline. The first step was to align and assemble the clones using [MiXCR 1.3.1](https://www.nature.com/articles/nmeth.3364) and then run them through quality control downsampling. 

After quality control, the samples were analyzed using a custom Bash/Python/R pipeline. 

### [Specific CDR3 Clone Tracking](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Plot_Functions/clonetracking.R) ###

To demonstrate accurate longitudinal tracking of specific clones, the clone tracking function was used to measure the relative size of each clone and track them over time. This function takes the aaCDR3 sequences from the reference (TIL Infusion Product) and track the abundance of them in the other samples and display the results as an alluvial. What this allows us to see is if an accurate response was seen in the patients and view certain characteristics of the infusion and baseline samples (richness, clonality, diversity, similarity) that may influence patient response. 

The clone tracking function will take in the patient code, specific sample cohort wanting to be analyzed (gDNA, cDNA, or cfDNA), chain/locus (TRA, TRB, TRD, TRG), and the clone fraction which will be used to specifying how much each clone has to take up in order to be incorporated into the analysis. The function will track the clones according to the specific reference. For this function, the reference is the clonotypes shown in the TIL infusion. The maximum clonotypes the function will track is 150. These 150 will each be assigned a specific color and the size will be tracked longitdinally throughout the samples based on the recurrence of its amino acid CDR3 sequence.

Here is an example of the clone tracking plot displaying the clonotypes >0.002 from the TRA chain of the gDNA samples in a specific patient (Patient and sample names have been cut out due to patient confidentiality):
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/CloneTrackingPlot.png)

### [Sample Diversity Function](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Plot_Functions/diversityevaluation.R) ###

To demonstrate the effect of diversity on patients and to compare sample type and loci measure of diversity the sample diversity function was used to compute the inverse-simpson diversity for each sample and compare them over time. This allows us to see a patients progression as a view of their whole T-cell repertoire as a measure of diversity. The inverse-simpson diversity was computed using the [Immunarch 4.0 package](https://github.com/immunomind/immunarch). 

The index allows us to measure and quantify the average abundance of types in the sample datasets and find the effective number of types. The sampleDiv function will produce a line plot with both the diversity computed for the DNA and RNA samples.

Here is an example of the diversity plot for a patient's TRA chain DNA and RNA samples (Patient and sample names have been cut out due to patient confidentiality):
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/SampleDiversityPlot.png)

### [Sample Clonality Comparison](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Plot_Functions/relativeabundance.R) ###

To demonstrate the effect of relative clonality on patient expansion and response, I used the clonal abundance measure from the Immunarch 4.0 package. This categorized wether the clones took up 0.005%, 0.05%, 0.5%, or 5% of the total sample. It then plotted these sections as a bar plot. 

Here is an example of the relative abundance for a patient's TRA chain DNA and RNA samples (Patient and sample names have been cut out due to patient confidentiality): 

![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/Alpha_Beta_relAb.png)

### [Repertoire Overlap Function](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Plot_Functions/repertoireoverlap.R) ###

To demonstrate how the dynamics of a persons T-cell repertoire effects the patients response, the repLap function allowed us to analyze how much of a persons immune response is overlapping between multiple samples. This analysis was done with a specific focus on looking at Baseline <> TIL Infusion Product overlap, Baseline <> Post-Infusion overlap and TIL Infusion Product <> Post-Infusion overlap.

The repLap function will take in a patient's sample for a specific sample cohort and will compute the overlap between the samples. It produces a heatmap with the number of overlapping clonotypes shown between the samples. We used the Immunarch 4.0 'repOverlap' function to sort through the samples and compute the overlap. 

Here is an example of a patient's repertoire overlap for their DNA samples TRA and TRB chains (Patient and sample names have been cut out due to patient confidentiality):
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/RepertoireOverlap.png)

### [TIL Clonality Comparison](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Plot_Functions/ClonalAbundanceTracking.R) ###

To demonstrate the relativity between patient's TIL Infusion clonality, the clonal abundance comparison function was used to view the top 10 clones in the TIL Infusion of each patient. This function takes in a list of samples that you want analyzed, the specific sample cohort (gDNA, cDNA, or cfDNA), chain/locus (TRA, TRB, TRD, TRG) and will create a stacked bar plot of the top 10 clones of those specific samples.

Here is an example of the clone tracking plot displaying the TIL Infusion products of all the patients in the research project for the TRA and TRB chains (Patient and sample names have been cut out due to patient confidentiality):
![image](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Images/ClonalAbundanceInfusionComparison.JPG)

### [Calculation Functions](https://github.com/CameronKerr/CapTCR_TIL_Tracking/blob/master/Calculation_functions/all.R) ###

To get an exact representation of many of the characteristics I was analyzing, I created multiple functions to calculated these numbers. A major example of these are the top 50% calculations. These functions take the top 50% of whatever desired sample (i.e TIL Infusion Product) and compares those clones to the resulting size of them in another sample (i.e 4 week sample). This was specifically useful to measure the level of selected expansion that occured in all patients. 

## Results ##

### Two-step and three-step CapTCR methods yeilded increased signal strength ###

Previously, the Pugh lab has described a hybrid capture method using a probe set targeting subregions of all annotated TCR locus V and J genes. The main challenge faced was the large number of reads representing fragments from unrearranged V and J genes when enriching. To demonstrate the most efficient order of molecular depletion in CapTCR, the lab experiment with an iterative two-step and three-step capture approach and analyzed the resulting on-target reads. The results displayed higher signal strengths with the iterative two-step and three-step capture methods.

### The level of expansion to the adoptive cell transfer of TILs can be tracked longitudinally using CapTCR with high accuracy ###

To demonstrate the effective monitoring of TIL ACT in peripheral blood samples, we tracked TIL-specific clonotypes over the patient timeline. The clone tracking function took the TIL-specific CDR3 sequences and filtered through the peripheral blood samples for the ones containing those sequences. Overall, we were able to see multiple levels of expansion, clonality, and diversity between patients. 

### The level of expansion to the adoptive cell transfer of TILs can be tracked longitudinally using CapTCR with high accuracy ###

Multiple studies of analyzing the T-cell repertoire after TIL infusion have shown that patients with high expansion had a high clonality post-infusion. Authors showed that despite infusions with polyclonality, a single immunodominant clone emerged post-infusion. To demonstrate the effect of the clonality of the baseline, TIL infusion, and post-infusion samples on expansion, we analyzed the diversity and clonal structures of the patients. In our study, a pattern was seen during the post-infusion samples. Despite clonality in baseline or TIL-infusions, a higher clonality correlated with a higher expansion. 

### Small clonotypes undergo selected-expansion of TIL-infusion clones post-infusion ###

In previous studies, selected-expansion occurs for the resulting TIL-dominant clones. Studies have shown in the past that for some patients, clones with a pre-infusion fraction of <0.001% demonstrated a high expansion. To demonstrate the selected-expansion of TIL-related clones we used the plots tracking TIL-specific clonotypes clone fraction to view specific clonal expansion. The results were suprising with the top percentage of clones in TIL Infusion Product ending up taking a small amount of clones 4 weeks post-infusion.


### Selected Clonal Expansion of TIL Infusion ###

Multiple in vivo dynamics hypotheses have been put forward in order to explain the selected-expansion of certain clonotypes. Authors in some studies have suggested that TIL culturing methods (antigen-pulsed dendritic cells in the presence of IL-2) result in naive T cell expansion in vitro in the TIL-infusion product but defer expansion in vivo throughout the peripheral blood samples. To demonstrate the correlation between baseline, TIL-infusion, and post-infusion overlap and patient survival, we analyzed the TCR repertoire overlap between patient samples. We observed that patients who had a higher level of ;arger than average TIL-related clones in their baseline would have a higher expansion post-infusion.
