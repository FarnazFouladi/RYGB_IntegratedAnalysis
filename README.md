# RYGB_IntegratedAnalysis
Steps toward making the figures in the manuscript:

1. set work directory at RYGB_IntegratedAnalysis

2. Running mixed linear models on each dataset to compare the gut microbiome at 
each time point to baseline.

 a. 16S datasets were run through both DADA2 and Kraken2:
 Afshar_DADA2_16S.R
 Afshar_kraken2_16S.R
 Assal_DADA2_16S.R
 Assal_kraken2_16S.R
 BS_DADA2_16S.R
 BS_kraken2_16S.R
 Ilhan_DADA2_16S.R
 Ilhan_kraken2_16S.R
 
 b. Metagenomics were run through Kraken2 and HUMAnN2:
 BS_kraken2_Metagenomics.R
 BS_pathway.R
 Palleja_kraken2_Metagenomics.R
 Palleja_pathway.R
 
3. Ordination plot (Figure 1A and B):
combineCountTables.R (We need to join count tables for ordination plot and machine learning)
pco.R

4. Scatter plots for 16S datasets at the genus level (Figure 1C and D, Supplemental Figure 2) and also the boxplots in Figure 1E and F

compareStudies_16S.R

5. Scatter plot for 16S datasets at the sequence variant level (Supplemental Figure 3)

compareStudies_SV.R

6. Scatter plot for the Supplemental Figure 4

SequenceVariantAnalysis.R

7. Diversity (Supplemental Figure 1)

Diversity.R

8. Heatmap in Figure 2

Heatmap.R

9. Scatter plots for Kraken2 and HUMAnN2 (Figure 3 and Supplemental 
Figures 5 and 6 )

compareStudies_Metagenomics.R
compareStudies_pathways.R

9. Machine leaning (Figure 4 and Supplemental Figure 7)

train_model.R (LASSO)
train_model_cluster.R (Random Forest takes long and I had to run it on the cluster
They are similar scripts!)

Figure_performance.R

10. Opportunistic pathogenes (Figure 5)

OpportunisticPathogens.R














 
 
 