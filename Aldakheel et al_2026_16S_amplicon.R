# This document contains all statistical analyses conducted for the manuscript.
# Due to the stochastic and iterative nature of certain analyses, some figure
# parameters may vary slightly upon reanalysis; however, the core results remain unchanged.
# For editing purposes, certain figure colors, shapes, and layout arrangements
# have been modified without affecting the underlying results.
# All data required to reproduce the analyses can be found here:
  #https://github.com/feixiang1209/PET-Cork-bacterial-consortia


library(ape)
library(taxize) 
library(seqinr)  
library(ShortRead) 
library(ggplot2); packageVersion("ggplot2")
library(dplyr)
library("dada2"); packageVersion("dada2") 
library(phyloseq)





fnFs <- sort(list.files(path, pattern="_F.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

#####Filter 
#Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

#####Inspect read quality profiles
#We start by visualizing the quality profiles of the forward reads:

plotQualityProfile(fnFs[1:20])
plotQualityProfile(fnRs[1:20])




#####Trimming

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 


write.table(out, paste(dir,"filtered_reads.csv",sep=""),row.names=T, sep=",")



#####Learn the Error Rates
#The DADA2 algorithm makes use of a parametric error model (err) and every amplicon dataset has a different set of error rates. The??learnErrors??method learns this error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. As in many machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

#Run and Wait
errF <- learnErrors(filtFs, multithread=TRUE)

#Run and Wait
errR <- learnErrors(filtRs, multithread=TRUE)



plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#####Dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names

names(derepFs) <- sample.names
names(derepRs) <- sample.names



#####Sample inference

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE) 

#####Merged paired reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


head(mergers[[1]])


#####Construct sequence table
head(mergers)
seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))


#####Remove chimeras

# THIS will give me the OTU table
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

write.table(seqtab.nochim, paste(dir,"seqtabOTU.csv",sep=""),row.names=T, sep=",")

dim(seqtab.nochim)


sum(seqtab.nochim)/sum(seqtab)

#####Track reads throuhg the pipeline 
#As a final check of our progress, we???ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track



write.table(track, paste(dir,"track_table.csv",sep=""),row.names=T, sep=",")

hist(log(colSums(seqtab.nochim))) 

saveRDS(seqtab.nochim, paste(dir,"seqtab.rds",sep=""))


############Assign taxonomy

seqtab <- readRDS(paste(dir,"seqtab.rds",sep=""))  
ncol(seqtab) 

# [1] 3516




taxseq<- read.fasta(file ="silva_nr99_v138.2_train_set.fa.gz", as.string = TRUE, strip.desc=TRUE)
fast1<-getAnnot(taxseq)##     head(taxseq)
head(fast1)  

taxass50<- assignTaxonomy(seqtab, "silva_nr99_v138.2_train_set.fa.gz", multithread = T, verbose = T,
                          taxLevels=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species"))  # minboot default is 50   MIDORI_LONGEST_1.1_COI_RDP.fasta     MIDORI_LONGEST_cut.fasta


taxa_sp <- addSpecies(taxass50, "silva_species_assignment_v138.1.fa.gz")


############################## put together
taxall50sp<-data.frame(taxa_sp)
mess<-taxall50sp %>% 
  dplyr::filter(complete.cases("Species")) %>% 
  filter(!duplicated(NA))

dim(taxall50sp)
rownames(taxall50sp)<-c(1:3516)
taxall50sp

write.table(taxall50sp, paste(dir,"taxa50sp.csv",sep=""),row.names=T, sep=",")

ASVtable_transpose<-t(seqtab) #transpose matrix

write.csv(ASVtable_transpose, file = "ASVtable_transpose.csv", na = "NA")  #export table

saveRDS(taxa_sp, paste(dir,"taxa_sp.rds",sep=""))

save.image("RImage_aftertax.RData")



######### END of DADA2  #############


##############################
#######   PHYLOSEQ   #########
##############################

### BEFORE START ######
#load the Dada2 environment (after TAX) with ALL the objects



### BEFORE START ######
#load the Dada2 environment (after TAX) with ALL the objects or load ONLY the "taxa_sp.rds" and "seqtab.rds" in order to create your Phyloseq object



### Load libraries
library(ape)
library(taxize) 
library(seqinr)  
library(ShortRead) 
library(ggplot2); packageVersion("ggplot2")
library(dplyr)
library("dada2"); packageVersion("dada2") 
library(phyloseq)

taxa_sp <- readRDS(paste(dir,"taxa_sp.rds",sep=""))

seqtab <- readRDS(paste(dir,"seqtab.rds",sep=""))

metadata <- read.csv("metadata.csv")

SampleID <- rownames(seqtab)

rownames(metadata) <- SampleID

physeq <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
                   tax_table(taxa_sp),
                   sample_data(metadata))



###### CLEANING DATA SET #######

#Removing Mitochondria 
physeq_NoM <- subset_taxa(physeq, (tax_table(physeq)[,"Family"]!="Mitochondria") | is.na(tax_table(physeq)[,"Family"]))


#Removing Chloro 
physeq_NoMC <- subset_taxa(physeq_NoM, (tax_table(physeq_NoM)[,"Order"]!="Chloroplast") | is.na(tax_table(physeq_NoM)[,"Order"]))



#Removing Euka 
physeq_NoMCE <- subset_taxa(physeq_NoMC, (tax_table(physeq_NoMC)[,"Kingdom"]!="Eukaryota") | is.na(tax_table(physeq)[,"Kingdom"]))

#Removing low abundance taxas (singleton)
physeq_1 = prune_taxa(taxa_sums(physeq_NoMCE) > 1, physeq_NoMCE) 


#### RARECURVE ANALYSIS ####

#Access to the ASV table for the specific phyloseq object
ASV_1 = otu_table(physeq_1)

##We convert the the ASV table into a matrix. It will warn you about the consecuences of the transformation
class(ASV_1) = "matrix"

rowSums(ASV_1) #Find the sample with smallest number of n reads

###Set up the line for the sample with the minimum number of reads
(raremax <- min(rowSums(ASV_1)))


####RUN the RARECURVE to compare the sequencing depth of each sample
rarecurve(ASV_1, step = 10, sample = raremax, col = "orange", cex = 0.5, label = F)


# Sample coverage index (Good’s coverage)  diagnostics
otu_mat <- as(otu_table(physeq_1), "matrix")


N  <- rowSums(otu_mat)
F1 <- apply(otu_mat, 1, function(x) sum(x == 1))
goods_C <- 1 - (F1 / N)

summary(goods_C)
sum(goods_C < 0.99)   # how many samples below 99% coverage
## 0 samples W/O enough coverage 

######ALPHA DIVERSITY ######

alpha_diversityestimate<-estimate_richness(physeq_1)

#Export the alpha diversity data. Now you can plot and run statistics on these values
write.csv(alpha_diversityestimate,"alpha_diversity.csv")


###REPORT CONTAMINANT BY SUBSTRATE ###
library(decontam); packageVersion("decontam")


physeq_ExpX1 = subset_samples(physeq_1, ExpX1=="Yes")
physeq_ExpX1 = prune_taxa(taxa_sums(physeq_ExpX1) > 0, physeq_ExpX1) 
#ExpX1-582
#ExpX2-1604
#ExpX3-1463
#ExpX4-299
#ExpX5-246



# Create a logical variable: TRUE if control, FALSE otherwise
sample_data(physeq_ExpX1)$is.neg <- sample_data(physeq_ExpX1)$Sample_or_control == "Control sample"

# Now run decontam with the logical column
contamdf.prev <- isContaminant(physeq_ExpX1, method = "prevalence", neg = "is.neg")

head(contamdf.prev)

table(contamdf.prev$contaminant)
#FALSE  TRUE 
#ExpX1-565     17 
#ExpX2-1601     3 
#ExpX3-1433     30 
#ExpX4-295     4 
#ExpX1-246      


which(contamdf.prev$contaminant)

#ExpX1[1]  14  16  19  24  25  37  49  51  57  71  76  82  85 117 137 223 255
#ExpX2[1] 37 200 238
#ExpX3[1]  14  38  45  56  60  62  65  78  82  88 123 170 179 181 189 194 213 214 250 273 277 328 337 339 356 360 367 400 461 507
#ExpX4[1] 88 104 108 109
#ExpX5 NO

physeq_ExpX1Dp <- prune_taxa(!contamdf.prev$contaminant, physeq_ExpX1)
#ExpX1-565
#ExpX2-1601
#ExpX3-1433
#ExpX4-295
#ExpX5-246

physeq_ExpX1Dp = prune_taxa(taxa_sums(physeq_ExpX1Dp) > 0, physeq_ExpX1Dp) 


asv_tab <- as.data.frame(otu_table(physeq_ExpX1Dp))
asv_tab <- as.data.frame(t(otu_table(physeq_ExpX1Dp)))
tax_tab <- as.data.frame(tax_table(physeq_ExpX1Dp))
asv_tax <- cbind(asv_tab, tax_tab[rownames(asv_tab), ])

write.csv(asv_tax, "ASV-TAX_ExpX1.csv")

write.csv(contamdf.prev, "contamASV_ExpX1.csv")

saveRDS(physeq_ExpX1Dp, file = "physeq_ExpX1Dp.rds")

### Repeat from lane 98 changing ExpX1 for X2, X3 and X4




### ONLY TAKE SAMPLES FROM EXPERIMENT 1, Transfer 1 to 5  ########
physeq_Exp1 = subset_samples(physeq_1, Exp1=="Yes")


#Clean dataset
physeq_Exp1 = prune_taxa(taxa_sums(physeq_Exp1) > 0, physeq_Exp1) 




######BETA DIVERSITY ######

set.seed(1)

physeq_Exp1_Group <- ordinate(
  physeq = physeq_Exp1, 
  method = "CAP",
  distance = "bray",
  formula = ~Substrate
)

# CAP plot
Cap_plot_Group <- plot_ordination(
  physeq = physeq_Exp1, 
  ordination = physeq_Exp1_Group , 
  color="Substrate",
  axes = c(1,2)
)

color20 <- c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6")


Cap_plot_Group + 
  geom_point(size=3)  + 
  geom_point(size=6, alpha = 0.7) + 
  geom_text(aes(label=SampleID),col="black", size=2.5) +
  
  stat_ellipse(aes(group =Substrate), level= 0.95, linetype=1, color="gray") +
  scale_color_manual(values = c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6"))




### Adonis ###

library(vegan)

physeq_Exp1_bray<-phyloseq::distance(physeq_Exp1,method="bray")

sampledf <- data.frame(sample_data(physeq_Exp1))

adonis2(physeq_Exp1_bray ~Transfer , data = sampledf)

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = physeq_Exp1_bray ~ Transfer, data = sampledf)
#Df SumOfSqs      R2      F Pr(>F)    
#Model     4   2.9823 0.15642 1.8542  0.001 ***
#Residual 40  16.0838 0.84358                  
#Total    44  19.0661 1.00000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##### RELATIVE ABUNDANCE AND BAR PLOTS ########
#Including NEG


physeq_Exp1N = subset_samples(physeq_1, Exp1N=="Yes")


#Clean dataset
physeq_Exp1N = prune_taxa(taxa_sums(physeq_Exp1N) > 0, physeq_Exp1N) 


#library(BiocManager)
#BiocManager::install("microbiome")
library(microbiome)


#install.packages("forcats")
library(forcats)
#It will manipulate the arrange and the legend of the bar plots

##Transfor to relative abundance
pseq.compositional <- microbiome::transform(physeq_Exp1N, "compositional")





#Setting up color palette
color20 <- c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6")


### PLOT FAMILY ##### 

# merge all taxa that are detected rare (change taxonomic level in each plot)
physeq_Exp1N_family <- aggregate_rare(pseq.compositional, level="Family", detection = 5/100, prevalence = 5/100)



#All % for the table
Plot_family <- physeq_Exp1N_family %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at Order level
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.03) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort data frame alphabetically by Order


ggplot(Plot_family, aes(x = Transfer, y = Abundance, fill =fct_reorder(Family, Abundance,.desc=TRUE ))) + 
  facet_grid(Substrate~.) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values =color20) + #### chnage the colors here. Puedo usar manual o paleta.
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance > 0.05 \n") +
  ggtitle("Family Composition") + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.2))




### PLOT GENUS ##### 

# merge all taxa that are detected rare (change taxonomic level in each plot)
physeq_Exp1_genus <- aggregate_rare(pseq.compositional, level="Genus", detection = 5/100, prevalence = 5/100)



#All % for the table
Plot_genus <- physeq_Exp1_genus %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Order level
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.03) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by Order



ggplot(Plot_genus, aes(x = SampleID, y = Abundance, fill =fct_reorder(Genus, Abundance,.desc=TRUE ))) + 
  facet_grid(Transfer~.) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values =color20) + #### chnage the colors here. Puedo usar manual o paleta.
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance > 0.03 \n") +
  ggtitle("Genus Composition") + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.2))





### ONLY TAKE SAMPLES FROM EXPERIMENT 2, Transfer 6 to 8  ########
physeq_Exp2 = subset_samples(physeq_1, Exp2=="Yes")


#Clean dataset
physeq_Exp2 = prune_taxa(taxa_sums(physeq_Exp2) > 0, physeq_Exp2) 




### PLOT GENUS ##### 


##Transfor to relative abundance
pseq.compositional <- microbiome::transform(physeq_Exp2, "compositional")


# merge all taxa that are detected rare (change taxonomic level in each plot)
physeq_Exp2_genus <- aggregate_rare(pseq.compositional, level="Genus", detection = 7/100, prevalence = 7/100)



#All % for the table
Plot_genus <- physeq_Exp2_genus %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Order level
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.03) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by Order



ggplot(Plot_genus, aes(x = SampleID, y = Abundance, fill =fct_reorder(Genus, Abundance,.desc=TRUE ))) + 
  facet_grid(Substrate~.) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values =color20) + #### chnage the colors here. Puedo usar manual o paleta.
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance > 0.03 \n") +
  ggtitle("Genus Composition") + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.2))

#theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.2))



## BETA DIVERSITY by substrate in Exp2 (T6 to T8) #####

### ONLY TAKE SAMPLES FROM CORK ########
physeq_Exp2c = subset_samples(physeq_Exp2, Substrate=="cork")


#Clean dataset
physeq_Exp2c = prune_taxa(taxa_sums(physeq_Exp2c) > 0, physeq_Exp2c) 

rownames(sample_data(physeq_Exp2c))

set.seed(1)

physeq_Exp2c_Group <- ordinate(
  physeq = physeq_Exp2c, 
  method = "CAP",
  distance = "bray",
  formula = ~Transfer
)

# CAP plot
Cap_plot_Group <- plot_ordination(
  physeq = physeq_Exp2c, 
  ordination = physeq_Exp2c_Group , 
  color="Transfer",
  axes = c(1,2)
)

color20 <- c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6")


Cap_plot_Group + 
  geom_point(size=3)  + 
  geom_point(size=6, alpha = 0.7) + 
  geom_text(aes(label=SampleID),col="black", size=2.5) +
  
  stat_ellipse(aes(group =Transfer), level= 0.95, linetype=1, color="gray") +
  scale_color_manual(values = c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6"))




### Adonis ###

library(vegan)

physeq_Exp2c_bray<-phyloseq::distance(physeq_Exp2c,method="bray")

sampledf <- data.frame(sample_data(physeq_Exp2c))

adonis2(physeq_Exp2c_bray ~Transfer , data = sampledf)

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = physeq_Exp2c_bray ~ Transfer, data = sampledf)
#Df SumOfSqs      R2     F Pr(>F)
#Model     2   0.4442 0.13445 1.165  0.292
#Residual 15   2.8598 0.86555             
#Total    17   3.3040 1.00000     



### ONLY TAKE SAMPLES FROM PET ########
physeq_Exp2p = subset_samples(physeq_Exp2, Substrate=="PET")


#Clean dataset
physeq_Exp2p = prune_taxa(taxa_sums(physeq_Exp2p) > 0, physeq_Exp2p) 

rownames(sample_data(physeq_Exp2p))


set.seed(1)

physeq_Exp2p_Group <- ordinate(
  physeq = physeq_Exp2p, 
  method = "CAP",
  distance = "bray",
  formula = ~Transfer
)

# CAP plot
Cap_plot_Group <- plot_ordination(
  physeq = physeq_Exp2p, 
  ordination = physeq_Exp2p_Group , 
  color="Transfer",
  axes = c(1,2)
)

color20 <- c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6")


Cap_plot_Group + 
  geom_point(size=3)  + 
  geom_point(size=6, alpha = 0.7) + 
  geom_text(aes(label=SampleID),col="black", size=2.5) +
  
  stat_ellipse(aes(group =Transfer), level= 0.95, linetype=1, color="gray") +
  scale_color_manual(values = c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6"))




### Adonis ###

library(vegan)

physeq_Exp2p_bray<-phyloseq::distance(physeq_Exp2p,method="bray")

sampledf <- data.frame(sample_data(physeq_Exp2p))

adonis2(physeq_Exp2p_bray ~Transfer , data = sampledf)

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = physeq_Exp2p_bray ~ Transfer, data = sampledf)
#Df SumOfSqs      R2      F Pr(>F)  
#Model     2   1.1563 0.27002 2.7742  0.022 *
#  Residual 15   3.1261 0.72998                
#Total    17   4.2824 1.00000                
---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
  
  
  
  
  ## BETA DIVERSITY by substrate in Exp2 (T6 to T8) #####

### ONLY TAKE SAMPLES FROM PET+CORK ########
physeq_Exp2pc = subset_samples(physeq_Exp2, Substrate=="PET+cork")


#Clean dataset
physeq_Exp2pc = prune_taxa(taxa_sums(physeq_Exp2pc) > 0, physeq_Exp2pc) 

rownames(sample_data(physeq_Exp2pc))

set.seed(1)

physeq_Exp2pc_Group <- ordinate(
  physeq = physeq_Exp2pc, 
  method = "CAP",
  distance = "bray",
  formula = ~Transfer
)

# CAP plot
Cap_plot_Group <- plot_ordination(
  physeq = physeq_Exp2pc, 
  ordination = physeq_Exp2pc_Group , 
  color="Transfer",
  axes = c(1,2)
)

color20 <- c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6")


Cap_plot_Group + 
  geom_point(size=3)  + 
  geom_point(size=6, alpha = 0.7) + 
  geom_text(aes(label=SampleID),col="black", size=2.5) +
  
  stat_ellipse(aes(group =Transfer), level= 0.95, linetype=1, color="gray") +
  scale_color_manual(values = c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6"))




### Adonis ###

library(vegan)

physeq_Exp2pc_bray<-phyloseq::distance(physeq_Exp2pc,method="bray")

sampledf <- data.frame(sample_data(physeq_Exp2pc))

adonis2(physeq_Exp2pc_bray ~Transfer , data = sampledf)

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = physeq_Exp2c_bray ~ Transfer, data = sampledf)
#Df SumOfSqs      R2     F Pr(>F)
#Model     2   1.7648 0.53895 6.4293  0.001 ***
#  Residual 11   1.5097 0.46105                  
#Total    13   3.2745 1.00000    



### ONLY TAKE SAMPLES FROM PET ########
physeq_Exp2p = subset_samples(physeq_Exp2, Substrate=="PET")


#Clean dataset
physeq_Exp2p = prune_taxa(taxa_sums(physeq_Exp2p) > 0, physeq_Exp2p) 

rownames(sample_data(physeq_Exp2p))


set.seed(1)

physeq_Exp2p_Group <- ordinate(
  physeq = physeq_Exp2p, 
  method = "CAP",
  distance = "bray",
  formula = ~Transfer
)

# CAP plot
Cap_plot_Group <- plot_ordination(
  physeq = physeq_Exp2p, 
  ordination = physeq_Exp2p_Group , 
  color="Transfer",
  axes = c(1,2)
)

color20 <- c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6")


Cap_plot_Group + 
  geom_point(size=3)  + 
  geom_point(size=6, alpha = 0.7) + 
  geom_text(aes(label=SampleID),col="black", size=2.5) +
  
  stat_ellipse(aes(group =Transfer), level= 0.95, linetype=1, color="gray") +
  scale_color_manual(values = c("#3166CC", "#CB4728", "#F19E38", "#449431", "#8D2D94", "#4397C1", "#CC5077", "#71A22F", "#A93A34", "#3D6291", "#8E4995", "#52A899", "#AAAA39", "#5F35C4", "#D87A2D", "#801F13", "#5D1B64", "#4F9066", "#5B73A2", "#3B3EA6"))




### Adonis ###

library(vegan)

physeq_Exp2p_bray<-phyloseq::distance(physeq_Exp2p,method="bray")

sampledf <- data.frame(sample_data(physeq_Exp2p))

adonis2(physeq_Exp2p_bray ~Transfer , data = sampledf)

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = physeq_Exp2p_bray ~ Transfer, data = sampledf)
#Df SumOfSqs      R2      F Pr(>F)  
#Model     2   1.1563 0.27002 2.7742  0.028 *
#  Residual 15   3.1261 0.72998                
#Total    17   4.2824 1.00000                
---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
  
  
  ####### EXP T #########
physeq_ExpT = subset_samples(physeq_1, ExpT=="Yes")
physeq_ExpT = prune_taxa(taxa_sums(physeq_ExpT) > 0, physeq_ExpT) 





physeq_ExpT <- microbiome::transform(physeq_ExpT, "compositional")

library(microbiome)

physeq_TAXgenus <- aggregate_taxa(physeq_ExpT, 'Genus')

#Taxonomy table at specific level
TAXgenus <- tax_table(physeq_TAXgenus)

write.table(TAXgenus, paste(dir,"TAXgenus.csv",sep=""),row.names=T, sep=",")

## Extract ASV table from phyloseq object  ####

pie<-as.matrix(physeq_TAXgenus@otu_table)
pie<-as.data.frame(pie)

ASV_genus<-pie 

write.csv(ASV_genus, file = "ASV_genus.csv", na = "NA")


###############################
#####Take an specific taxa (with ALL the ASVs) across de samples #####


######For Table 1 ##########


physeq.T1 <- subset_taxa(physeq_TAXgenus, Genus %in% c("Brevibacillus", "Bacillus", "Paenibacillus"))

ExpT1<-as.matrix(physeq.T1@otu_table)
ExpT1<-as.data.frame(ExpT1)

ASV_ExpT1<-ExpT1
ASV_ExpT1<-t(ExpT1)
write.csv(ASV_ExpT1, file = "ASV_ExpT1.csv", na = "NA")



# --- Extract data from the phyloseq object ---
otu_df  <- as.data.frame(otu_table(physeq.T1))
tax_df  <- as.data.frame(tax_table(physeq.T1))
meta_df <- as.data.frame(sample_data(physeq.T1))

# --- Ensure OTUs are rows (transpose if needed) ---
if (!taxa_are_rows(physeq.T1)) {
  otu_df <- t(otu_df)
}

# --- Combine OTU table with taxonomy by taxa (rows) ---
otu_tax <- cbind(otu_df, tax_df)

# --- Transpose so samples become rows, matching metadata ---
otu_tax_t <- as.data.frame(t(otu_tax))

# --- Merge metadata with abundance+taxonomy table ---
# Make sure sample names match between metadata and OTU table
merged_df <- cbind(meta_df, otu_tax_t[rownames(meta_df), ])

# --- Export the merged table as a single CSV ---
write.csv(merged_df, file = "phyloseq_T1_merged.csv", row.names = TRUE)




######For Table 2 ##########


physeq.T2 <- subset_taxa(physeq_TAXgenus, Genus %in% c("Amycolatopsis", "Nocardiopsis", "Streptomyces"))

ExpT2<-as.matrix(physeq.T2@otu_table)
ExpT2<-as.data.frame(ExpT2)

ASV_ExpT2<-ExpT2
ASV_ExpT2<-t(ExpT2)
write.csv(ASV_ExpT2, file = "ASV_ExpT2.csv", na = "NA")



# --- Extract data from the phyloseq object ---
otu_df  <- as.data.frame(otu_table(physeq.T2))
tax_df  <- as.data.frame(tax_table(physeq.T2))
meta_df <- as.data.frame(sample_data(physeq.T2))

# --- Ensure OTUs are rows (transpose if needed) ---
if (!taxa_are_rows(physeq.T2)) {
  otu_df <- t(otu_df)
}

# --- Combine OTU table with taxonomy by taxa (rows) ---
otu_tax <- cbind(otu_df, tax_df)

# --- Transpose so samples become rows, matching metadata ---
otu_tax_t <- as.data.frame(t(otu_tax))

# --- Merge metadata with abundance+taxonomy table ---
# Make sure sample names match between metadata and OTU table
merged_df <- cbind(meta_df, otu_tax_t[rownames(meta_df), ])

# --- Export the merged table as a single CSV ---
write.csv(merged_df, file = "phyloseq_T2_merged.csv", row.names = TRUE)



######For Table 3 ##########


physeq.T3 <- subset_taxa(physeq_TAXgenus, Genus %in% c("Nitratireductor", "Pseudaminobacter", "Ramlibacter"))

ExpT3<-as.matrix(physeq.T3@otu_table)
ExpT3<-as.data.frame(ExpT3)

ASV_ExpT3<-ExpT3
ASV_ExpT3<-t(ExpT3)
write.csv(ASV_ExpT3, file = "ASV_ExpT3.csv", na = "NA")



# --- Extract data from the phyloseq object ---
otu_df  <- as.data.frame(otu_table(physeq.T3))
tax_df  <- as.data.frame(tax_table(physeq.T3))
meta_df <- as.data.frame(sample_data(physeq.T3))

# --- Ensure OTUs are rows (transpose if needed) ---
if (!taxa_are_rows(physeq.T3)) {
  otu_df <- t(otu_df)
}

# --- Combine OTU table with taxonomy by taxa (rows) ---
otu_tax <- cbind(otu_df, tax_df)

# --- Transpose so samples become rows, matching metadata ---
otu_tax_t <- as.data.frame(t(otu_tax))

# --- Merge metadata with abundance+taxonomy table ---
# Make sure sample names match between metadata and OTU table
merged_df <- cbind(meta_df, otu_tax_t[rownames(meta_df), ])

# --- Export the merged table as a single CSV ---
write.csv(merged_df, file = "phyloseq_T3_merged.csv", row.names = TRUE)



#############################################
########### END OF THE SCRIPT ##############
#############################################





