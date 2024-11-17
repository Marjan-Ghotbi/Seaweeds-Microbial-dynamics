

# PERMANOVA for seaweed and their proxy biofilms

# Load required packages 
library(qiime2R)
library(phyloseq)
library(vegan)
library(permute)
library(lattice)
library(pairwiseAdonis)
library(microbiomeutilities)
library(microbiome)
library(vegan)


# Read phyloseq object of prokaryotis community 
#(85 samples)
prok.phyloseq<- readRDS("prokaryotes-phyloseq.rds" )
prok.phyloseq

#phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 12685 taxa and 85 samples ]:
#sample_data() Sample Data:        [ 85 samples by 18 sample variables ]:
#tax_table()   Taxonomy Table:     [ 12685 taxa by 7 taxonomic ranks ]:
#phy_tree()    Phylogenetic Tree:  [ 12685 tips and 12651 internal nodes ]:

#We do analysis on three different subsets 

#Tab.S5-A. PERMANOVA to compare a pair set of seaweeds and their PBs after removal of control to have a balanced design(68 samples)



physeq.pair <- subset_samples ( prok.phyloseq, Isolation_source != "PC-Control" )
physeq<-physeq.pair
physeq


#Tab.S5-B. PERMANOVA to compare all PBs and control filters (49 samples of PC)
physeqP <- subset_samples ( prok.phyloseq, Substrate == "Polycarbonate" )
physeq<-physeqP
physeq

#Tab.S5-C. PERMANOVA to compare seaweeds biofilm composition 
physeqS <- subset_samples ( prok.phyloseq, Substrate == "Seaweed" )
physeq<-physeqS
physeq


# to check the minimum depth of each subset
sample_depths <- sample_sums(physeq)
min_depth <- min(sample_depths)
print(min_depth)


#Set the seed and rarefy at 2000 for each of the sebsets 
physeq_rare<-rarefy_even_depth(physeq,rngseed = 02032021,sample.size = 2000)
physeq_rare

# we can check the metadata and samples in each category to make sure filtration worked 
sample_data(physeq_rare)
table(meta(physeq_rare)$Isolation_source, useNA = "always")


#PC-F. serratus     PC-F. vesiculosus PC-G. vermiculophylla        Sw-F. serratus     Sw-F. vesiculosus 
#      12                    11                    12                    12                    12 
#Sw-G. vermiculophylla                  <NA> 
#     12                                 0 


library(vegan)
physeq_rel = transform_sample_counts(physeq_rare, function(x)100* x/sum(x) )

bray_dist = phyloseq::distance(physeq_rel, method="bray", weighted=T)
ordination = phyloseq::ordinate(physeq_rel, method="PCoA", distance=bray_dist)

#We used Bray-Curtis dissimilarity matrices but in case Unifrac is required 
wunifrac_dist = phyloseq::distance(physeq_rel, method="unifrac", weighted=T)
ordination = phyloseq::ordinate(physeq_rel, method="PCoA", distance=wunifrac_dist)




# Introduce these factors
sample_data(physeq)


Panel<-sample_data(physeq_rel)$Panel
bottle<-sample_data(physeq_rel)$bottle
Substrate<-sample_data(physeq_rel)$Substrate
Treatment<-sample_data(physeq_rel)$Treatment
Timepoint<-sample_data(physeq_rel)$Timepoint
Isolation_source <-sample_data(physeq_rel)$Isolation_source

#Either intriduce combined_strata concept here or we can use Panel:bottle in the formula for strata 
# Panel: to account for non-independence among observations from bottles installed on the same panel 
# bottle: to account for non-independence between filters and seaweeds from the same bottle over two experiments

metadata<-meta(physeq_rel)
dim(metadata)


metadata$combined_strata <- interaction (Panel, bottle, drop = TRUE)
metadata$combined_strata

levels(metadata$combined_strata)
#pairs 18
#PBs and control   7

ad<-adonis2(bray_dist ~ Substrate* Treatment * Timepoint, data=meta,  strata = combined_strata , permutations = 9999, method="bray")
ad

#or

ad<-adonis2(bray_dist ~ Substrate* Treatment * Timepoint, strata = Panel:bottle, permutations = 9999, method="bray")
ad




#Pairs 
pairwise.adonis2 (bray_dist ~ Treatment , factors = Treatment, data=metadata,    strata =  18 ,    perm = 9999,    p.adjust= "Holm",  method = "bray")

#PBs and Control 
# strata =  7


############################## 
# Read phyloseq object of microalgae community 
#82 samples

micro.phyloseq<- readRDS("Microalgae-phyloseq.rds" )
micro.phyloseq

#phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 12685 taxa and 85 samples ]:
#sample_data() Sample Data:        [ 85 samples by 18 sample variables ]:
#tax_table()   Taxonomy Table:     [ 12685 taxa by 7 taxonomic ranks ]:
#phy_tree()    Phylogenetic Tree:  [ 12685 tips and 12651 internal nodes ]:

#We do analysis on three different subsets 

#Tab.S5-A. PERMANOVA to compare a pair set of seaweeds and their PBs after removal of control to have a balanced design(68 samples)

physeq.pair <- subset_samples ( micro.phylose, Isolation_source != "PC-Control" )
physeq<-physeq.pair
physeq


#Tab.S5-B. PERMANOVA to compare all PBs and control filters (49 samples of PC)
physeqP <- subset_samples ( micro.phylose, Substrate == "Polycarbonate" )
physeq<-physeqP
physeq

#Tab.S5-C. PERMANOVA to compare seaweeds biofilm composition 
physeqS <- subset_samples ( micro.phylose, Substrate == "Seaweed" )
physeq<-physeqS
physeq


# to check the minimum depth of each subset
sample_depths <- sample_sums(physeq)
min_depth <- min(sample_depths)
print(min_depth)


#Set the seed and rarefy at 2000 for each of the sebsets 
physeq_rare<-rarefy_even_depth(physeq,rngseed = 02032021,sample.size = 2000)
physeq_rare

# we can check the metadata and samples in each category to make sure filtration worked 
sample_data(physeq_rare)
table(meta(physeq_rare)$Isolation_source, useNA = "always")


#PC-F. serratus     PC-F. vesiculosus PC-G. vermiculophylla        Sw-F. serratus     Sw-F. vesiculosus 
#      12                    11                    12                    12                    12 
#Sw-G. vermiculophylla                  <NA> 
#     12                                 0 


library(vegan)
physeq_rel = transform_sample_counts(physeq_rare, function(x)100* x/sum(x) )

bray_dist = phyloseq::distance(physeq_rel, method="bray", weighted=T)
ordination = phyloseq::ordinate(physeq_rel, method="PCoA", distance=bray_dist)

#We used Bray-Curtis dissimilarity matrices but in case Unifrac is required 
wunifrac_dist = phyloseq::distance(physeq_rel, method="unifrac", weighted=T)
ordination = phyloseq::ordinate(physeq_rel, method="PCoA", distance=wunifrac_dist)


# Introduce these factors
sample_data(physeq)


Panel<-sample_data(physeq_rel)$Panel
bottle<-sample_data(physeq_rel)$bottle
Substrate<-sample_data(physeq_rel)$Substrate
Treatment<-sample_data(physeq_rel)$Treatment
Timepoint<-sample_data(physeq_rel)$Timepoint
Isolation_source <-sample_data(physeq_rel)$Isolation_source

#Either intriduce combined_strata concept here or we can use Panel:bottle in the formula for strata 
# Panel: to account for non-independence among observations from bottles installed on the same panel 
# bottle: to account for non-independence between filters and seaweeds from the same bottle over two experiments
metadata<-meta(physeq_rel)
dim(metadata)


metadata$combined_strata <- interaction (Panel, bottle, drop = TRUE)
metadata$combined_strata

levels(metadata$combined_strata)
#pairs 18
#PBs and control   7



ad<-adonis2(bray_dist ~ Substrate* Treatment * Timepoint, strata = combined_strata , permutations = 9999, method="bray")
ad

#or

ad<-adonis2(bray_dist ~ Substrate* Treatment * Timepoint, strata = Panel:bottle, permutations = 9999, method="bray")
ad


pairwise.adonis2 (bray_dist ~ Treatment , factors = Treatment, data=metadata   strata =  18 ,    perm = 9999,    p.adjust= "Holm",  method = "bray")
