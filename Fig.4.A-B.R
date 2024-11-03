

#Seaweeds Biofilm  PCoA for prokaryotic and microalgae community 
#Fig.4- A & B  

#Prokaryotic community 


library(qiime2R)
library(phyloseq)
library(vegan)
library(permute)
library(lattice)
library(pairwiseAdonis)
library("metagMisc")
library(microbiomeutilities)
library(speedyseq)
library(microbiome)
library(eulerr)
library(microbiome)
library(microbiomeutilities)
library(vegan)


# if prok= euk no of samples =83



#prokaryotes (85 samples)
# First we make a phyloseq object

#85 samples 

#How we constructed phyloseq object after decontamination and reclassification of Unassigned ASVs
#The phyloseq object is saved and provided 

prok.physeq<-qza_to_phyloseq( features="table_SW-PC-no-Mit_Chl-Removed_Unassigned.Decontam_Ultimate_Silva_IQ3return_31Dec2023_gc.norm.qza",
                         tree="rooted-tree.qza",
                         "reformat.taxonomy_after_IQ3return_Final-ultimate.qza",
                         metadata = "metadata-T1-T2-SWPC-nowater-for-Adonis.txt")
prok.physeq


saveRDS(physeq, file = "prok.physeq.rds")
physeq <- readRDS("prok.physeq.rds")

#phyloseq-class experiment-level object
#otu_table()   OTU Table:          [ 12685 taxa and 85 samples ]:
#sample_data() Sample Data:        [ 85 samples by 18 sample variables ]:
#tax_table()   Taxonomy Table:     [ 12685 taxa by 7 taxonomic ranks ]:
#phy_tree()    Phylogenetic Tree:  [ 12685 tips and 12651 internal nodes ]:
#taxa are rows



sample_depths <- sample_sums(physeq)
min_depth <- min(sample_depths)
print(min_depth)

#[1] 2074.04 min sample depth for prokaryotic data 


physeq_rare<-rarefy_even_depth(physeq,rngseed = 02032021,sample.size = 2000)

physeq_rare

# we can check the associated metadata 
sample_data(physeq_rare)

table(meta(physeq_rare)$Isolation_source, useNA = "always")

#PC-Control          PC-F. serratus     PC-F. vesiculosus        PC-G. vermiculophylla     Sw-F. serratus     Sw-F. vesiculosus     Sw-G. vermiculophylla               
# 14                     12                     11                   12                       12                    12                      12



library(vegan)
physeq_rel = transform_sample_counts(physeq_rare, function(x)100* x/sum(x) )

bray_dist = phyloseq::distance(physeq_rel, method="bray", weighted=T)
ordination = phyloseq::ordinate(physeq_rel, method="PCoA", distance=bray_dist)
Bray_distance<-as.matrix(bray_dist)


library(ggplot2)
beta
beta<-phyloseq::plot_ordination(physeq_rel, ordination , color= "Treatment", shape= "Substrate") + scale_colour_manual(values=c("G. vermiculophylla" = "red","G. vermiculophylla" = " red", "F. serratus"="blue" ,"F. serratus"="blue",  "F. vesiculosus"="#FF66CC","F. vesiculosus"="#FF66CC","Control" ="#00CCCC","Ambient water"="green"))+ 
  scale_shape_manual(values=c("Ambient water"=15, "Polycarbonate"= 16, "Seaweed"=17))+theme_bw()+
  scale_alpha_manual(values = c("7-days"=0.4,"14-days"=1.5),name="Exposure time")


bbeta<- beta + geom_point(aes( color=Treatment,  shape= Substrate,   alpha= Exposure_time ),size=5)+theme_bw()   #facet_grid(cols = vars(Treatment), scales="free_y", switch = 'y')
bbbeta<-bbeta +labs(title="PCoA based on Bray Curtis distance matrix", x ="PCoA 1 (19.5 %)", y = " PCoA 2 (11.1 %)")
bbeta
ooo=bbbeta+theme(strip.text.x = element_text(family="Times New Roman", size = 12, color="black",face = "bold"))
v=ooo+  theme(strip.text.y=element_text(family="Times New Roman", size =9, color="black",face = "bold"))
R5=v+theme(strip.background = element_rect(colour = "black", fill = "gray98"))+ ggtitle("")
R5

#we can add title to ggtitle("") if needed 

ggsave("Prok-PCoA-Bray-6-7.pdf", height=6, width=7, device="pdf")



############################## 

#Microalgae community
#We first make phyloseq object

microalgae.physeq<-qza_to_phyloseq( features="Silva_chloroplast_SW-PC.table-NoMacrohost-50-52-58.qza",
                         tree="rooted-tree.qza",
                         "classification.qza",
                         metadata = "metadata-Microalgae-min50-58-52-nowater.txt")
microalgae.physeq


saveRDS(physeq, file = "microalgae.physeq.rds")
physeq <- readRDS("microalgae.physeq.rds")


sample_depths <- sample_sums(physeq)
min_depth <- min(sample_depths)
print(min_depth)

#min depth : 103


physeq_rare<-rarefy_even_depth(physeq,rngseed = 02032020,sample.size = 100)

physeq_rel = transform_sample_counts(physeq_rare, function(x)100* x/sum(x) )

bray_dist = phyloseq::distance(physeq_rel, method="bray", weighted=T)
bray_dist


#Important consideration: Here we use Bray-curtis matrices, but for Unifrac it is important to have exact right tree of the samples we have. In case some samples are deleted for decobtam or rarefaction
#now the old tree is not a right match and  a new tree is needed


physeq_rel = transform_sample_counts(physeq, function(x)100* x/sum(x) )

bray_dist = phyloseq::distance(physeq_rel, method="bray", weighted=T)
ordination = phyloseq::ordinate(physeq_rel, method="PCoA", distance=bray_dist)



beta<-phyloseq::plot_ordination(physeq_rel, ordination , color= "Treatment", shape= "Substrate") + scale_colour_manual(values=c("G. vermiculophylla" = "red","G. vermiculophylla" = " red", "F. serratus"="blue" ,"F. serratus"="blue",  "F. vesiculosus"="#FF66CC","F. vesiculosus"="#FF66CC","Control" ="#00CCCC","Ambient water"="green"))+ 
  scale_shape_manual(values=c("Ambient water"=15, "Polycarbonate"= 16, "Seaweed"=17))+theme_bw()+
  scale_alpha_manual(values = c("7-days"=0.4,"14-days"=1.5),name="Exposure time")
beta


bbeta<- beta + geom_point(aes( color=Treatment,  shape= Substrate,   alpha= Exposure_time ),size=5)+theme_bw()   #facet_grid(cols = vars(Treatment), scales="free_y", switch = 'y')
bbbeta<-bbeta +labs(title="PCoA based on Bray Curtis distance matrix", x ="PCoA 1 (28.7 %)", y = " PCoA 2 (15.6 %)")
bbeta

ooo=bbbeta+theme(strip.text.x = element_text(family="Times New Roman", size = 12, color="black",face = "bold"))
v=ooo+  theme(strip.text.y=element_text(family="Times New Roman", size =9, color="black",face = "bold"))
R5=v+theme(strip.background = element_rect(colour = "black", fill = "gray98"))+ ggtitle("Microalgae_WUni")
R5


ggsave("Microalgae-PCoA-Bray-6-7.pdf", height=6, width=7, device="pdf")





