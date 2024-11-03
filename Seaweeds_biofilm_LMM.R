# Load necessary libraries

library(ggplot2)
library(lme4)
library(multcomp)
library(dplyr)
library(lmerTest)
library(lme4)
library(emmeans)
library(qiime2R)

meta.shannon <-read_q2metadata("Tab.S1.Shannon_Environment.txt")

dim(meta.shannon)
set.seed(123)

duplicated(names(meta.shannon))

#to clean for duplicated names
#Keep data that you need considering NAs in the data at each step

meta.shannon2<-meta.shannon[,-(24:27)]
dim(meta.shannon2)

# We remove control filters, and later ambient water and compare paired samples of seaweeds and their PBs to have a balance design 

Paired <- meta.shannon2 %>% filter (Isolation_source!= "PC-Control")
dim(Paired)
# 73, 21

# We omit  ambient water samples and any other sample that does not have Shannon index 
Paired_na<- na.delete(Paired)
dim(Paired_na)
#68 21

# We can check diversity values and run a shapiro test to check if the data is normally distributed or not.
#We can alwys condider scaling for normalization
Paired_na$prokaryote_shannon
shapiro_test(Paired_na$prokaryote_shannon)

Paired_na$microalgae_shannon
shapiro_test(Paired_na$microalgae_shannon)


# We can check the number of samples in each group for both experiments (Timepoints)
table(meta(Paired_na)$Isolation_source, useNA = "always")

#PC-F. serratus     PC-F. vesiculosus PC-G. vermiculophylla        Sw-F. serratus     Sw-F. vesiculosus   Sw-G. vermiculophylla   
#12                    11                    12                    11                    11                 11



#We need to introduce all factors

Isolation_source<- as.factor(Paired_na$Isolation_source)  
Timepoint<- as.factor(Paired_na$Timepoint)  
Substrate <- as.factor(Paired_na$Substrate) 
Treatment<- as.factor(Paired_na$Treatment) 
Panel<- as.factor(Paired_na$Panel) 
bottle<- as.factor(Paired_na$bottle) 


#Prokaryotic community LMM analysis 
#we scale the shanoon indices

scale_Paird <- Paired_na %>% mutate_at(c("prokaryote_shannon"), ~(scale(.) %>% as.vector)) # it only scales prok_shannon col
scale_Paird 

alphab_Paired<-as.data.frame(scale_Paird)
alphab_Paired


# Now we run the LMM

lmeModel1 = lmer( prokaryote_shannon ~   Substrate * Treatment * Timepoint  + (1|Panel:bottle), data = alphab_Paired)
lmeModel1

lmeModel1
summary(lmeModel1)
Anova(lmeModel1)

pairs(emmeans(lmeModel1, specs = "Treatment", by ="Timepoint"), adjust = "Holm")


#######################
# We only select for non-living substrate (Polycarbonate filters) , and  PBs and controls 

nonliving<- meta.shannon2 %>% filter( Substrate == "Polycarbonate")
dim(nonliving)


table(meta(nonliving)$Isolation_source, useNA = "always")

scale_nonliving <- nonliving %>% mutate_at(c("prokaryote_shannon"), ~(scale(.) %>% as.vector)) # it only scales prok_shannon col
scale_nonliving

alphab_nonliving<-as.data.frame(scale_nonliving)
alphab_nonliving

lmeModel2 = lmer(prokaryote_shannon ~   Treatment*Timepoint  + (1|Panel:bottle), data = alphab_nonliving)
lmeModel2

summary(lmeModel2)
Anova(lmeModel2)

pairs(emmeans(lmeModel2, specs = "Timepoint"), adjust = "Holm")
pairs(emmeans(lmeModel2, specs = "Treatment", by ="Timepoint"), adjust = "Holm")


########################

# We do the same analysis for microalgae 
#paired samples

scale_Paird <- Paired_na %>% mutate_at(c("microalgae_shannon"), ~(scale(.) %>% as.vector)) # it only scales prok_shannon col
scale_Paird

alphab_Paired<-as.data.frame(scale_Paird)
alphab_Paired


# Now we run the LMM

lmeModel1 = lmer( microalgae_shannon ~   Substrate * Treatment * Timepoint  + (1|Panel:bottle), data = alphab_Paired)
lmeModel1

summary(lmeModel1)
Anova(lmeModel1)

pairs(emmeans(lmeModel1, specs = "Treatment", by ="Timepoint"), adjust = "Holm")


#############################
#microalgae on non-living substrate

scale_nonliving <- nonliving %>% mutate_at(c("microalgae_shannon"), ~(scale(.) %>% as.vector)) # it only scales prok_shannon col
scale_nonliving

alphab_nonliving<-as.data.frame(scale_nonliving)
alphab_nonliving

lmeModel1 = lmer(microalgae_shannon ~   Treatment*Timepoint  + (1|Panel:bottle), data = alphab_nonliving)
lmeModel1

summary(lmeModel1)
Anova(lmeModel1)

pairs(emmeans(lmeModel1, specs = "Timepoint"), adjust = "Holm")
pairs(emmeans(lmeModel1, specs = "Treatment", by ="Timepoint"), adjust = "Holm")


########################


# Analysis of microalgae enumeration

meta.shannon <-read_q2metadata("Tab.S1.Shannon_Environment.txt")


dim(meta.shannon)
set.seed(123)

#to clean for duplicated names

duplicated(names(meta.shannon))
meta.shannon2<-meta.shannon[,-(26)]

meta.shannon2$Plastid_count
shapiro_test(meta.shannon2$Plastid_count)



#We delet NAs and will keep plastid count just on PBs and control filters
Plastid<- na.delete(meta.shannon2)
dim(Plastid)



scale_Plastid <- Plastid %>% mutate_at(c("Plastid_count"), ~(scale(.) %>% as.vector)) # it only scales prok_shannon col
scale_Plastid


lmeModel1 = lmer(Plastid_count ~ Treatment * Timepoint  + (1|Panel:bottle), data = scale_Plastid)
lmeModel1
summary(lmeModel1)
Anova(lmeModel1)

pairs(emmeans(lmeModel1, specs = "Timepoint" ), adjust = "Holm")
pairs(emmeans(lmeModel1, specs = "Treatment", by ="Timepoint"), adjust = "Holm")

