
system.time( load("C:/Users/ernie/OneDrive/Desktop/Current Projects/RDE/RDE_16S_table.bin") )

system.time( load("C:/Users/ernie/OneDrive/Desktop/Current Projects/RDE/RDE_16S_taxa.bin") )


library(phyloseq)

###We can construct a simple sample data.frame from the information encoded in the filenames. Usually this step would instead involve reading the sample data in from a file.

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
samdf <- read.csv("C:/Users/ernie/OneDrive/Desktop/Current Projects/RDE/16S_meta_RDE.csv")
rownames(samdf) <- samples.out


###We now construct a phyloseq object directly from the dada2 outputs.

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#Change ASV names to something more manageable than the full DNA sequence

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


#Check sequence depth. Hoping for at least 10,000 sequences per sample. 

rowSums(otu_table(ps))

min(rowSums(otu_table(ps)))

ps <- rarefy_even_depth(ps,rngseed=TRUE)

ps

rowSums(otu_table(ps))

###some initial processing of the data for multivariate stats using vegan

library(vegan)

ps <- subset_samples(ps, id !="RDE24")

ps2 <- subset_samples(ps, Diversity != "Neg")

ps2 <- subset_samples(ps2, Diversity != "Pos")

a <- data.frame(otu_table(ps))

a2 <- data.frame(otu_table(ps2))

a <- data.frame(decostand(a,method="hellinger"))

a2 <- data.frame(decostand(a2,method="hellinger"))

b <- vegdist(a, method="bray")

b2 <- vegdist(a2, method="bray")

c <- data.frame(sample_data(ps))

c2 <- data.frame(sample_data(ps2))


#Alpha diversity analyses

bac_alpha <- estimate_richness(ps, measures=c("Shannon","Observed"))

bac_alpha <- cbind(bac_alpha, c) 

bac_alpha2 <-  bac_alpha[bac_alpha$Diversity!="Pos",]

bac_alpha2 <- bac_alpha2[bac_alpha2$Diversity!="Neg",]

glm1 = lm(Shannon~ LandUse*Diversity*Stress,data=bac_alpha2)

shapiro.test(resid(glm1))

library(car)

Anova(glm1)

library(emmeans)

em1 <- emmeans(glm1, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

aggregate(Observed~Diversity,data=bac_alpha,FUN=mean)

aggregate(Shannon~Diversity,data=bac_alpha,FUN=mean)

aggregate(Shannon~LandUse,data=bac_alpha,FUN=mean)


#Plot 16S Shannon diversity for all treatments

sum1 <- aggregate(Shannon~ LandUse+Diversity+Stress, data=bac_alpha, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum1 <- do.call(data.frame, sum1)
sum1

sum1$se <- sum1$Shannon.sd / sqrt(sum1$Shannon.n)
head(sum1)

sum1$Diversity <- factor(sum1$Diversity, levels = c("Pos","D0","D1","D2","Neg"))


#Now for some plotting

library(ggplot2)
library(RColorBrewer)

brewer.pal(n=5,"Set1")


sum1$Stress <- factor(sum1$Stress, levels =c("Control","Antibiotic"))

line_plot_shan <- ggplot(sum1, aes(x=Diversity, y=Shannon.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=Shannon.mean-se, ymax=Shannon.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_manual(values=c("#984EA3","#E41A1C", "#377EB8", "#4DAF4A","#FF7F00")) +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('16S Shannon Diversity')) +
  facet_wrap(~ LandUse)
plot(line_plot_shan)

#Plot 16S Richness diversity for all treatments

sum1 <- aggregate(Observed~ LandUse+Diversity+Stress, data=bac_alpha, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum1 <- do.call(data.frame, sum1)
sum1

sum1$se <- sum1$Observed.sd / sqrt(sum1$Observed.n)
head(sum1)

sum1$Diversity <- factor(sum1$Diversity, levels = c("Pos","D0","D1","D2","Neg"))


sum1$Stress <- factor(sum1$Stress, levels =c("Control","Antibiotic"))


#Now for some plotting

library(ggplot2)
library(RColorBrewer)

brewer.pal(n=5,"Set1")


line_plot_shan <- ggplot(sum1, aes(x=Diversity, y=Observed.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=Observed.mean-se, ymax=Observed.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_manual(values=c("#984EA3","#E41A1C", "#377EB8", "#4DAF4A","#FF7F00")) +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('16S ASV Richness')) +
  facet_wrap(~ LandUse)
plot(line_plot_shan)


#Plot 16S diversity excluding diversity controls

sum2 <- aggregate(Shannon~ LandUse+Diversity+Stress, data=bac_alpha2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum2 <- do.call(data.frame, sum2)
sum2

sum2$se <- sum2$Shannon.sd / sqrt(sum2$Shannon.n)
head(sum2)


sum2$Stress <- factor(sum2$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

line_plot_shan2 <- ggplot(sum2, aes(x=Diversity, y=Shannon.mean, color=Diversity, shape=Stress)) + 
  geom_point(data=sum2,position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=5) +
  geom_errorbar(data=sum2,position=position_dodge(width=0.9),aes(ymin=Shannon.mean-se, ymax=Shannon.mean+se), width=.5, size=.75) +
  annotate("text", 2.5, 5, label="Land~Use: ~italic(P) < 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", 2.5, 4.75, label="Diversity: ~italic(P) < 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", 2.5, 4.5, label="Stress: ~italic(P) == 0.597", parse=TRUE, size=4.5, fontface="bold")+
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('16S Shannon Diversity')) +
  facet_wrap(~ LandUse)
plot(line_plot_shan2)



#Multivariate analyses

library(vegan)

set.seed(101)
adonis2(b2~LandUse*Diversity*Stress, data=c2,permutations=999)


#Dispersion test

Diversity <- c2$Diversity

disp <- betadisper(b2,Diversity,type="centroid")
anova(disp)


#PCoA

e <- cmdscale(b2, k=2,eig=TRUE)

#The code below calculates the % of variation in composition accounted for by the first two PCoA axes

var1 <- round(e$eig[1]/sum(e$eig)*100,1)
var2 <- round(e$eig[2]/sum(e$eig)*100,1)

#All of the code below is taking the PCOA result and converting it into a data frame that we can plot

data.scores <- as.data.frame(e$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$id <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$LandUse <- c2$LandUse 
data.scores$Diversity <- c2$Diversity
data.scores$Stress <-c2$Stress
head(data.scores) 

data.scores$Stress <- factor(data.scores$Stress, levels=c("Control","Antibiotic"))

library(ggplot2)

ord_1 <- ggplot(data.scores, aes(x=V1, y=V2,shape=Stress,color=Diversity,alpha=LandUse)) + 
   geom_point(data=data.scores, aes(x=V1,y=V2,shape=Stress, color=Diversity),size=5) + # add the point markers
  annotate("text", 0, .05, label="Land~Use: ~italic(P) < 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", 0, 0, label="Diversity: ~italic(P) < 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", 0, -.05, label="Stress: ~italic(P) == 0.258", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", 0, -.1, label="`Land Use x Diversity`: ~italic(P) < 0.001", parse=TRUE, size=4.5, fontface="bold")+
  scale_alpha_discrete(range = c(0.35, 1))+
  scale_color_brewer(palette="Set1") +
  labs(alpha='Land Use') +
  xlab("PCoA 1 (32.2%)") +
  ylab("PCoA 2 ( 13.6%)") +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color='Diversity') +
  labs(shape='Stress') +
  theme(text=element_text(size=15)) +
  ggtitle("16S ASVs") +
  theme(legend.text=element_text(size=15)) 
ord_1

library(cowplot)

jpeg(filename="Figure1.jpeg", bg="transparent", res=500, units = "in", height=6, width=16) 

f4 <- plot_grid(line_plot_shan2, ord_1, ncol = 2, align="hv", label_size=25,labels=c('A', 'B'))

f4

dev.off()


#Include diversity controls


e2 <- cmdscale(b, k=2,eig=TRUE)
e2 <- metaMDS(b, k=2, trymax=1000)


#All of the code below is taking the PCOA result and converting it into a data frame that we can plot

data.scores2 <- as.data.frame(e2$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores2$id <- rownames(data.scores2)  # create a column of site names, from the rownames of data.scores
data.scores2$LandUse <- c$LandUse 
data.scores2$Diversity <- c$Diversity
data.scores2$Stress <- c$Stress
head(data.scores2) 

data.scores2$Diversity <- factor(data.scores2$Diversity, levels = c("Pos","D0","D1","D2","Neg"))

data.scores2$Stress <- factor(data.scores2$Stress, levels=c("Control","Antibiotic"))


ord_2 <- ggplot(data.scores2, aes(x=MDS1, y=MDS2,shape=Stress,color=Diversity,alpha=LandUse)) + 
  geom_point(data=data.scores2, aes(x=MDS1,y=MDS2,shape=Stress, color=Diversity),size=4) + # add the point markers
   scale_alpha_discrete(range = c(0.35, 1))+
  scale_color_manual(values=c("#984EA3","#E41A1C", "#377EB8", "#4DAF4A","#FF7F00")) +
  labs(alpha='Land Use') +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color='Diversity') +
  labs(shape='Stress') +
  theme(text=element_text(size=15)) +
  ggtitle("16S ASVs") +
  theme(legend.text=element_text(size=15)) 
ord_2


#Bar plots of bacterial phyla

#Simple bar plots of relative abundance of bacterial ph

rank_names(ps)
table(tax_table(ps)[, "Phylum"], exclude = NULL)

#bacterial phyla stuff

rk <- transform_sample_counts(ps, function(x) x/sum(x))

# agglomerate taxa
glom_rk <- tax_glom(rk, taxrank = 'Phylum')

# create dataframe from phyloseq object
dat_rk <- psmelt(glom_rk)
head(dat_rk)

dat_rk0 <- psmelt(glom_rk)

# convert phylum to a character vector 
dat_rk$Phylum <- as.character(dat_rk$Phylum)

dat_rk0$Phylum <- as.character(dat_rk0$Phylum)

library(plyr)
# group dataframe by phylum, calculate mean rel. abundance
means <- ddply(dat_rk, ~Phylum, function(x) c(mean=mean(x$Abundance)))
means

# find Phyla whose mean rel. abund. is less than 1%
remainder <- means[means$mean <= 0.01,]$Phylum

# change their name to "Other"
dat_rk[dat_rk$Phylum %in% remainder,]$Phylum <- 'Other'

#Get the data ready for analysis and processing


rk_phyla <- aggregate(Abundance~id+Phylum+LandUse+Diversity+Stress, dat_rk, FUN=sum)

rk_phyla0 <- aggregate(Abundance~id+Phylum+LandUse+Diversity+Stress, dat_rk0, FUN=sum)


library(reshape)

rk_phyla0.1 <- cast(rk_phyla0, id+LandUse+Diversity+Stress ~ Phylum, value="Abundance")

rk_phyla1 <- cast(rk_phyla, id+LandUse+Diversity+Stress ~ Phylum, value="Abundance")


#Calculate bacerial r:K index. Can be useful for looking at ecological responses to treatments

rk_phyla1$rk <- (rk_phyla1$Bacteroidota+rk_phyla1$Proteobacteria)/(rk_phyla1$Acidobacteriota+rk_phyla1$Actinobacteriota)

rk_phyla2 <- rk_phyla1[rk_phyla1$Diversity!="Pos",]

rk_phyla2 <- rk_phyla2[rk_phyla2$Diversity!="Neg",]



p1 <- lm(rk~ LandUse*Diversity*Stress, data=rk_phyla2)
shapiro.test(resid(p1))
p1.1<-glm(rk~LandUse*Diversity*Stress,data=rk_phyla2,family=Gamma(link=log))

aggregate(rk ~ Diversity+Stress,data=rk_phyla2,FUN=mean)

AIC(p1,p1.1)

library(car)

Anova(p1.1)


em1 <- emmeans(p1, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")


sum_rk <- aggregate(rk~ LandUse+Diversity+Stress, data=rk_phyla2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_rk <- do.call(data.frame, sum_rk)
sum_rk

sum_rk$se <- sum_rk$rk.sd / sqrt(sum_rk$rk.n)
head(sum_rk)

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

brewer.pal(n=4,"Set1")

sum_rk$Stress <- factor(sum_rk$Stress, levels=c("Control","Antibiotic"))

line_plot_rk <- ggplot(sum_rk, aes(x=Diversity, y=rk.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=rk.mean-se, ymax=rk.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('16S r:K')) +
  facet_wrap(~ LandUse)
plot(line_plot_rk)



p1 <- lm(Proteobacteria~ LandUse*Diversity*Stress, data=rk_phyla2)
shapiro.test(resid(p1))
p1.1<-glm(Proteobacteria~LandUse*Diversity*Stress,data=rk_phyla2,family=Gamma(link=log))

aggregate(Proteobacteria ~Diversity+Stress, data=rk_phyla2, FUN=mean)

Anova(p1.1)

em1 <- emmeans(p1.1, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")


p1 <- lm(Actinobacteriota~ LandUse*Diversity*Stress, data=rk_phyla2)
shapiro.test(resid(p1))
p1.1<-glm(Actinobacteriota~LandUse*Diversity*Stress,data=rk_phyla2,family=Gamma(link=log))

aggregate(Actinobacteriota ~Stress+Diversity, data=rk_phyla2, FUN=mean)


Anova(p1.1)

em1 <- emmeans(p1.1, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")


p1 <- lm(Bacteroidota~ LandUse*Diversity*Stress, data=rk_phyla2)
shapiro.test(resid(p1))
p1.1<-glm(Bacteroidota+.001~LandUse*Diversity*Stress,data=rk_phyla2,family=Gamma(link=log))

Anova(p1.1)

p1 <- lm(Firmicutes~ LandUse*Diversity*Stress, data=rk_phyla2)
shapiro.test(resid(p1))
p1.1<-glm(Firmicutes+.001~LandUse*Diversity*Stress,data=rk_phyla2,family=Gamma(link=log))

aggregate(Firmicutes ~Stress+Diversity+LandUse, data=rk_phyla2, FUN=mean)


Anova(p1.1)

em1 <- emmeans(p1.1, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")


#Now plotting

#This is across the two sites


jpeg(filename="bac_bar0.jpeg", bg="transparent", res=500, units = "in", height=5, width=5) 


bac_bar <-ggplot(rk_phyla, aes(fill=Phylum, y=Abundance, x=LandUse)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Relative Abundance") +
  scale_fill_brewer(palette="Set3") +
  theme_classic() +
  theme(legend.box.margin=margin(-10,-10,-10,-10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Bacterial Phyla") +
  xlab("") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) 
bac_bar

dev.off()

#This is across the diversity treatments


bac_bar <-ggplot(rk_phyla, aes(fill=Phylum, y=Abundance, x=Diversity)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Relative Abundance") +
  scale_fill_brewer(palette="Set3") +
  theme_classic() +
  theme(legend.box.margin=margin(-10,-10,-10,-10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Bacterial Phyla") +
  xlab("Diversity") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  facet_wrap(~ LandUse)
bac_bar

#This is the diversity treatments exclusing diversity controls and Stress soils

rk_phyla0.1 <-  rk_phyla[rk_phyla$Diversity!="Neg",]

rk_phyla0.1 <-  rk_phyla0.1[rk_phyla0.1$Diversity!="Pos",]

rk_phyla0.1 <-  rk_phyla0.1[rk_phyla0.1$Stress!="Antibiotic",]

jpeg(filename="bac_bar.jpeg", bg="transparent", res=500, units = "in", height=6, width=6) 

bac_bar <-ggplot(rk_phyla0.1, aes(fill=Phylum, y=Abundance, x=Diversity)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Relative Abundance") +
  scale_fill_brewer(palette="Set3") +
  theme_classic() +
  theme(legend.box.margin=margin(-10,-10,-10,-10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Bacterial Phyla") +
  xlab("Diversity") +
  theme(axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  facet_wrap(~ LandUse)
bac_bar

dev.off()

#This is controls and Stress within the cultivated site

rk_phyla2 <- rk_phyla[rk_phyla$LandUse=="Cultivated",]

rk_phyla2 <- rk_phyla2[rk_phyla2$Diversity!="Pos",]

rk_phyla2 <- rk_phyla2[rk_phyla2$Diversity!="Neg",]

rk_phyla2$Stress <- factor(rk_phyla2$Stress,levels=c("Control","Antibiotic"))

jpeg(filename="bac_bar1.jpeg", bg="transparent", res=500, units = "in", height=6, width=9) 


bac_bar2 <-ggplot(rk_phyla2, aes(fill=Phylum, y=Abundance, x=Stress)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Relative Abundance") +
  scale_fill_brewer(palette="Set3") +
  theme_classic() +
  theme(legend.box.margin=margin(-10,-10,-10,-10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Cultivated") +
  xlab("") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=18)) +
  facet_wrap(~ Diversity)
bac_bar2

dev.off()

#This is controls and Stress within the prairie site

rk_phyla3 <- rk_phyla[rk_phyla$LandUse=="Prairie",]

rk_phyla3 <- rk_phyla3[rk_phyla3$Diversity!="Pos",]

rk_phyla3 <- rk_phyla3[rk_phyla3$Diversity!="Neg",]


rk_phyla3$Stress <- factor(rk_phyla3$Stress, levels=c("Control","Antibiotic"))

jpeg(filename="bac_bar2.jpeg", bg="transparent", res=500, units = "in", height=6, width=9) 



bac_bar3 <-ggplot(rk_phyla3, aes(fill=Phylum, y=Abundance, x=Stress)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous("Relative Abundance") +
  scale_fill_brewer(palette="Set3") +
  theme_classic() +
  #theme(legend.box.margin=margin(-10,-10,-10,-10))+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Prairie") +
  xlab("") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=18)) +
  facet_wrap(~ Diversity)
bac_bar3

dev.off()

library(cowplot)

jpeg(filename="Figure2.jpeg", bg="transparent", res=500, units = "in", height=12, width=9) 

f4 <- plot_grid(bac_bar2, bac_bar3, ncol = 1, align="hv", label_size=25,labels=c('A', 'B'))

f4

dev.off()


#16S copy number stuff (picrust2)

library(biomformat)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(car)
library(lme4)
library(reshape)
library(emmeans)
library(cowplot)

Bac_otus2 <- "C:/Users/ernie/OneDrive/Desktop/Current Projects/RDE/16S-table-RDE.biom"

x1 <- read_biom(Bac_otus2)

x2 <- as.matrix(biom_data(x1))

head(x2)

colSums(x2)

copies2 <- read.csv("C:/Users/ernie/OneDrive/Desktop/Current Projects/RDE/picrust_RDE_filtered.csv")

x3.1 <- x2[rownames(x2) %in% copies2$sequence,]

colSums(x3.1)

min(colSums(x3.1))

x4.1 <- t(x3.1)

x5.1 <- rrarefy(x4.1, 14747)

rowSums(x5.1)

library(funrar)

x6.1 <- make_relative(x5.1)

rowSums(x6.1)

copies2 <- copies2[order(copies2$sequence),]

x7.1 <- t(x6.1)

x8.1 <- cbind(rownames(x7.1), data.frame(x7.1, row.names=NULL,check.names=FALSE))

colnames(x8.1)[1] <- c("sequence")

x8.1 <- x8.1[order(x8.1$sequence),]

x9.1 = sapply(x8.1[,c(2:81)], '*', copies2$X16S_rRNA_Count) 

colSums(x9.1)

table.1 <- data.frame(t(x9.1))

table.1$copynumber <- rowSums(table.1)

table2.1 <- cbind(rownames(table.1), data.frame(table.1, row.names=NULL,check.names=FALSE))

colnames(table2.1)[1] <- c("id")

table2.1 <- table2.1[order(table2.1$id),]

meta <- read.csv("C:/Users/ernie/OneDrive/Desktop/Current Projects/RDE/16S_meta_RDE.csv")

meta <- meta[order(meta$id),]

table2.1 <- cbind(table2.1, meta)

table3.1 <- table2.1[,c(1,5843:5849)]

table4.1 <- table3.1[table3.1$Diversity!="Pos",]

table4.1 <- table4.1[table4.1$Diversity!="Neg",]

p1 <- lm(copynumber~ LandUse*Diversity*Stress, data=table4.1)
shapiro.test(resid(p1))
p1.1<-glm(copynumber~LandUse*Diversity*Stress,data=table4.1,family=Gamma(link=log))

AIC(p1,p1.1)

Anova(p1.1)

library(emmeans)

em1 <- emmeans(p1.1, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_operons <- aggregate(copynumber~ LandUse+Diversity+Stress, data=table4.1, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_operons <- do.call(data.frame, sum_operons)
sum_operons

sum_operons$se <- sum_operons$copynumber.sd / sqrt(sum_operons$copynumber.n)
head(sum_operons)

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

brewer.pal(n=4,"Set1")

sum_operons$Stress <- factor(sum_operons$Stress, levels=c("Control","Antibiotic"))

plot_operons <- ggplot(sum_operons, aes(x=Diversity, y=copynumber.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=copynumber.mean-se, ymax=copynumber.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('Average 16S operon number')) +
  facet_wrap(~ LandUse)
plot(plot_operons)


#basal Resp

s <- read.csv("C:/Users/ernie/OneDrive/Desktop/Current Projects/RDE/RDE_summary.csv")

s2 <-  s[s$Diversity!="Pos",]

s2 <- s2[s2$Diversity!="Neg",]

m1 <- lm(Cmin~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m1))

library(car)

qqPlot(resid(m1))

m2 <- glm(Cmin~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m1,m2)

Anova(m2)

aggregate(Cmin~LandUse+Diversity,data=s2,FUN=mean)

aggregate(Cmin~LandUse,data=s2,FUN=mean)

aggregate(Cmin~Stress, data=s2,FUN=mean)

em1 <- emmeans(m2, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_cmin <- aggregate(Cmin ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_cmin <- do.call(data.frame, sum_cmin)
sum_cmin

sum_cmin$se <- sum_cmin$Cmin.sd / sqrt(sum_cmin$Cmin.n)
head(sum_cmin)


sum_cmin$Stress <- factor(sum_cmin$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

plotCmin <- ggplot(sum_cmin, aes(x=Diversity, y=Cmin.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=5) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=Cmin.mean-se, ymax=Cmin.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position = "none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  scale_y_continuous(bquote(~'C Mineralization'~~ (mu*'g'~'CO'[2]~'-C'~~g^-1~~h^-1)))+
  facet_wrap(~ LandUse)
plot(plotCmin)

#SIR

m1 <- lm(SIR~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m1))

library(car)

qqPlot(resid(m1))

m2 <- glm(SIR~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m1,m2)

Anova(m2)

aggregate(SIR~Stress+Diversity+LandUse,data=s2,FUN=mean)

aggregate(SIR~Diversity+LandUse,data=s2,FUN=mean)

em1 <- emmeans(m2, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_sir <- aggregate(SIR ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_sir <- do.call(data.frame, sum_sir)
sum_sir

sum_sir$se <- sum_sir$SIR.sd / sqrt(sum_sir$SIR.n)
head(sum_sir)

sum_sir$Stress <- factor(sum_sir$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

plotSIR <- ggplot(sum_sir, aes(x=Diversity, y=SIR.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=5) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=SIR.mean-se, ymax=SIR.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position = "none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  scale_y_continuous(bquote(~'SIR'~~ (mu*'g'~'CO'[2]~'-C'~~gdw^-1~~h^-1)))+
  facet_wrap(~ LandUse)
plot(plotSIR)

#qCO2

m2 <- lm(qCO2~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(qCO2~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m3)

aggregate(qCO2~Stress+LandUse+Diversity,data=s2,FUN=mean)

aggregate(qCO2~Diversity,data=s2,FUN=mean)

em1 <- emmeans(m3, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_q <- aggregate(qCO2 ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_q <- do.call(data.frame, sum_q)
sum_q

sum_q$se <- sum_q$qCO2.sd / sqrt(sum_q$qCO2.n)
head(sum_q)

sum_q$Stress <- factor(sum_q$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)


plotQCO2 <- ggplot(sum_q, aes(x=Diversity, y=qCO2.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=5) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=qCO2.mean-se, ymax=qCO2.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position = "none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  scale_y_continuous(bquote(~'qCO'[2]~~ (mu*'g'~'CO'[2]~'-C'~~'mg Mic C'^-1~~h^-1))) +
  facet_wrap(~ LandUse)
plot(plotQCO2)


#Total enzyme activity

m2 <- lm(totalH~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(totalH~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m3)

aggregate(totalH~LandUse,data=s2,FUN=mean)

aggregate(totalH~Stress+Diversity+LandUse,data=s2,FUN=mean)

em1 <- emmeans(m2, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_enz <- aggregate(totalH ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_enz <- do.call(data.frame, sum_enz)
sum_enz

sum_enz$se <- sum_enz$totalH.sd / sqrt(sum_enz$totalH.n)
head(sum_enz)

sum_enz$Stress <- factor(sum_enz$Stress,levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)


plotENZ <- ggplot(sum_enz, aes(x=Diversity, y=totalH.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=5) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=totalH.mean-se, ymax=totalH.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.position = "none") +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(expression("Total enzyme activity" ~~("nmol" ~"gdw"^-1* ~"h"^-1))) +
  facet_wrap(~ LandUse)
plot(plotENZ)

#N mineralization

m2 <- lm(Nmin~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(Nmin~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m2)

aggregate(Nmin~LandUse,data=s2,FUN=mean)

aggregate(Nmin~Stress+Diversity+LandUse,data=s2,FUN=mean)

em1 <- emmeans(m2, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_nmin <- aggregate(Nmin ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_nmin <- do.call(data.frame, sum_nmin)
sum_nmin

sum_nmin$se <- sum_nmin$Nmin.sd / sqrt(sum_nmin$Nmin.n)
head(sum_nmin)

sum_nmin$Stress <- factor(sum_nmin$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)


plotNMIN <- ggplot(sum_nmin, aes(x=Diversity, y=Nmin.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=5) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=Nmin.mean-se, ymax=Nmin.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(legend.position = "none") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  scale_y_continuous(bquote(~'N mineralization'~~ (mu*'g N'~~gdw^-1~~d^-1))) +
  facet_wrap(~ LandUse)
plot(plotNMIN)


#Nitrification

m2 <- lm(Nit~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

Anova(m2)

aggregate(Nit~LandUse,data=s2,FUN=mean)

aggregate(Nit~Stress+Diversity+LandUse,data=s2,FUN=mean)

em1 <- emmeans(m2, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_nit <- aggregate(Nit ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_nit <- do.call(data.frame, sum_nit)
sum_nit

sum_nit$se <- sum_nit$Nit.sd / sqrt(sum_nit$Nit.n)
head(sum_nit)

sum_nit$Stress <- factor(sum_nit$Stress,levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)


plotNIT <- ggplot(sum_nit, aes(x=Diversity, y=Nit.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=5) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=Nit.mean-se, ymax=Nit.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(legend.position = "none") +
  scale_y_continuous(bquote(~'Nitrification'~~ (mu*'g N'~~gdw^-1~~d^-1))) +
  facet_wrap(~ LandUse)
plot(plotNIT)


library(cowplot)

jpeg(filename="Figure3.jpeg", bg="transparent", res=600, units = "in", height=12, width=18) 

f4 <- plot_grid(plotCmin,plotSIR,plotQCO2,plotENZ,plotNMIN,plotNIT, ncol = 3, align="hv", label_size=30,labels=c('A', 'B','C','D','E','F'))

f4

dev.off()

#Extract legend

jpeg(filename="legend.jpeg", bg="transparent", res=600, units = "in", height=6, width=6) 


plotNIT <- ggplot(sum_nit, aes(x=Diversity, y=Nit.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=4) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=Nit.mean-se, ymax=Nit.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  scale_y_continuous(bquote(~'Nitrification'~~ (mu*'g N'~~gdw^-1~~d^-1))) +
  facet_wrap(~ LandUse)
plot(plotNIT)

dev.off()



#Multifunctionality

s3 <- s2[,c(1:5,16:20,25)]

library(vegan)

s4 <- data.frame(scale(s3[,c(6:11)]))

s5 <- cbind(s3[,1:5],s4)

s5$multi <- rowMeans(s5[,c(6:11)])

s2$multi <- s5$multi

a<- vegdist(s5[,c(6:11)], method="euclidean", na.rm=TRUE)

set.seed(101)
adonis2(a~LandUse*Diversity*Stress,data=s5,permutations=999)


#Now visualize the profiles using principle components analysis (PCA)

s6 <- na.omit(s5)

b <- princomp(s6[,c(6:11)], cor=FALSE,scores=TRUE)

summary(b)

b$loadings

data.scores3 <- as.data.frame(b$scores)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores3$LandUse <- s6$LandUse
data.scores3$Diversity <- s6$Diversity
data.scores3$Stress <- s6$Stress

multi_ord4 <- aggregate(Comp.1 ~ LandUse+Diversity+Stress, data.scores3,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord4

multi_ord5 <- aggregate(Comp.2 ~ LandUse+Diversity+Stress, data.scores3,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord5

multi_ord4 <- do.call(data.frame, multi_ord4)
multi_ord4

multi_ord5 <- do.call(data.frame, multi_ord5)
multi_ord5

multi_ord6 <- cbind(multi_ord4, multi_ord5[,c(4:6)])
multi_ord6

multi_ord6$Comp1se <- multi_ord6$Comp.1.sd / sqrt(multi_ord6$Comp.1.n)
multi_ord6

multi_ord6$Comp2se <- multi_ord6$Comp.2.sd / sqrt(multi_ord6$Comp.2.n)
multi_ord6

library(ggplot2)


multi_ord6$Stress <- factor(multi_ord6$Stress, levels=c("Control","Antibiotic"))

ord_multi <- ggplot(multi_ord6, aes(x=Comp.1.mean, y=Comp.2.mean, color=Diversity,shape=Stress,alpha=LandUse)) + 
  geom_errorbar(aes(ymin=multi_ord6$Comp.2.mean-multi_ord6$Comp2se, ymax=multi_ord6$Comp.2.mean+multi_ord6$Comp2se), width=0.1, size=.5, position=position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=multi_ord6$Comp.1.mean-multi_ord6$Comp1se, xmax=multi_ord6$Comp.1.mean+multi_ord6$Comp1se), height=0.1, size=.5) +
  geom_point(data=multi_ord6,aes(x=Comp.1.mean,y=Comp.2.mean, color=Diversity, shape=Stress),size=5) + # add the point markers
  xlab("PCA 1 (45%)" ) +
  ylab("PCA 2 (25%)") +
  annotate("text", -.5, -1, label="`Land Use:` ~italic(P) < 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.5, -1.4, label="`Diversity:` ~italic(P) < 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.5, -1.8, label="`Stress:` ~italic(P) == 0.005", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.5, -2.2, label="`Land Use x Diversity:` ~italic(P) < 0.001", parse=TRUE, size=4.5, fontface="bold")+
  annotate("text", -.5, -2.6, label="`Land Use x Stress:` ~italic(P) == 0.012", parse=TRUE, size=4.5, fontface="bold")+
  scale_color_brewer(palette = "Set1") +
  scale_alpha_discrete(range = c(0.35, 1))+
  labs(alpha='Land Use') +
  theme(legend.text=element_text(size=14)) +
  theme_classic() +
  theme(axis.text=element_text(size=18)) +
  ggtitle("Multifunctionality") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=18)) 
ord_multi


#Include controls

s3.1 <- s[,c(1:5,16:20,25)]

library(vegan)


s4.1 <- data.frame(scale(s3.1[,c(6:11)]))

s5.1 <- cbind(s3.1[,1:5],s4.1)

s5.1$multi <- rowMeans(s5.1[,c(6:11)])

s$multi <- s5.1$multi


s6.1 <- na.omit(s5.1)

b <- princomp(s6.1[,c(6:11)], cor=FALSE,scores=TRUE)

summary(b)

b$loadings

data.scores4 <- as.data.frame(b$scores)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores4$LandUse <- s6.1$LandUse
data.scores4$Diversity <- s6.1$Diversity
data.scores4$Stress <- s6.1$Stress

multi_ord4 <- aggregate(Comp.1 ~ LandUse+Diversity+Stress, data.scores4,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord4

multi_ord5 <- aggregate(Comp.2 ~ LandUse+Diversity+Stress, data.scores4,FUN= function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x)))
multi_ord5

multi_ord4 <- do.call(data.frame, multi_ord4)
multi_ord4

multi_ord5 <- do.call(data.frame, multi_ord5)
multi_ord5

multi_ord6 <- cbind(multi_ord4, multi_ord5[,c(4:6)])
multi_ord6

multi_ord6$Comp1se <- multi_ord6$Comp.1.sd / sqrt(multi_ord6$Comp.1.n)
multi_ord6

multi_ord6$Comp2se <- multi_ord6$Comp.2.sd / sqrt(multi_ord6$Comp.2.n)
multi_ord6

library(ggplot2)

ord_multi <- ggplot(multi_ord6, aes(x=Comp.1.mean, y=Comp.2.mean, color=Diversity,shape=Stress)) + 
  geom_errorbar(aes(ymin=multi_ord6$Comp.2.mean-multi_ord6$Comp2se, ymax=multi_ord6$Comp.2.mean+multi_ord6$Comp2se), width=0.2, size=.5, position=position_dodge(0.9)) +
  geom_errorbarh(aes(xmin=multi_ord6$Comp.1.mean-multi_ord6$Comp1se, xmax=multi_ord6$Comp.1.mean+multi_ord6$Comp1se), height=0.2, size=.5) +
  geom_point(data=multi_ord6,aes(x=Comp.1.mean,y=Comp.2.mean, color=Diversity, shape=Stress),size=4) + # add the point markers
  xlab("PCA 1 (39%)" ) +
  ylab("PCA 2 (26%)") +
  #annotate("text", 2.5, 1.5, label="`Land Use:` ~italic(P) == 0.666", parse=TRUE, size=4.5, fontface="bold")+
  #annotate("text", 2.5, 1.25, label="`Watershed Pair:` ~italic(P) == 0.001", parse=TRUE, size=4.5, fontface="bold")+
  #annotate("text", 2.5, 1, label="`Land Use x Pair:` ~italic(P) == 0.267", parse=TRUE, size=4.5, fontface="bold")+
  scale_color_brewer(palette = "Set1") +
  #labs(color="Land Use", shape="Watershed Pair")+
  theme(legend.text=element_text(size=14)) +
  theme_classic() +
  theme(axis.text=element_text(size=16)) +
  ggtitle("Multifunctionality") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=14)) +
  facet_wrap(~ LandUse)
ord_multi


###Random Forest models

#Incorporate microbial variables

s7 <- s[order(s$id),]

b_alpha <- bac_alpha

b_alpha2 <- b_alpha[,c(1,2,3)]

b_alpha3 <- rbind(b_alpha2, "RDE24" = c(NA, NA, "RDE24"))

b_alpha3 <- b_alpha3[order(b_alpha3$id),]

b_pcoa <- data.scores2

b_pcoa2 <- b_pcoa[,c(1:3)]

b_pcoa3 <- rbind(b_pcoa2, "RDE24" = c(NA, NA, "RDE24"))

b_pcoa3 <- b_pcoa3[order(b_pcoa3$id),]

b_tax <- rk_phyla1[,c(1,14)]

b_tax2 <- rbind(b_tax, "RDE24" = c("RDE24", NA))

b_tax3 <- b_tax2[order(b_tax2$id),]

s7$bac_shannon <- b_alpha3$Shannon

s7$bac_richness <- b_alpha3$Observed

s7$pcoa1 <- b_pcoa3$V1

s7$pcoa2 <- b_pcoa3$V2

s7$rk <- b_tax3$rk

s7$copynumber <- table3.1$copynumber

s7.1 <- s7[,c(6:34)]

s7.1 <- data.frame(sapply(s7.1,as.numeric))

library(caret)

pre <- preProcess(s7.1, method='knnImpute')
pre

# Use the imputation model to predict the values of missing data points
library(RANN)  # required for knnImpute

train <- predict(pre, s7.1) 


cor.test(train$bac_shannon,train$Cmin)

plot(y=train$Cmin,x=train$bac_shannon)


jpeg(filename="bac_scatter.jpeg", bg="transparent", res=500, units = "in", height=6, width=8) 


cor1 <- ggplot(train, aes(x=bac_shannon, y=Cmin)) + 
  geom_smooth(method='lm', colour="black",linetype="dashed") +
  geom_point(size=4) +   
  #scale_color_brewer(palette="Set1", name="Land Use") +
  annotate("text", label="r == 0.52~~~~italic(P) < 0.001", x=2, y=.6, size=6, parse=TRUE) +
  scale_y_continuous(bquote(~'Respiration'~~ (mu*'g'~'CO'[2]~'-C'~~gdw^-1~~h^-1)))+
  labs(x="16S Shannon Diversity")+
  theme_classic() +
  theme(text=element_text(size=18)) +
  theme(axis.text=element_text(size=18)) +
  theme(legend.position=c("none")) +
  theme(axis.line.x=element_line(colour="black", size=.5)) + 
  theme(axis.line.y=element_line(colour="black", size=.5)) 
cor1

dev.off()


#C mineralization random forest

resp <- train[,c(11,1:10,24,26:29,18)]

library(randomForest)

set.seed(101)

r_tune <- tuneRF(train[,c(1:10,24,26:29,18)],train$Cmin,ntreeTry=10000,stepFactor=1.5, improve=0.05,
                 trace=TRUE, plot=TRUE, doBest=FALSE)

print(r_tune)


set.seed(131)
resp.rf <- randomForest(Cmin ~ ., data=resp,mtry=4,ntree=10000,
                        importance=TRUE)
print(resp.rf)

## Show "importance" of variables: higher value mean more important:
round(importance(resp.rf), 2)

varImpPlot(resp.rf,type=1)

## Show "importance" of variables: higher value mean more important:
rimp <- data.frame(round(importance(resp.rf,type=1), 2))

library(data.table)
setDT(rimp, keep.rownames = TRUE)

rimp2 <- rimp[order(rimp$X.IncMSE,decreasing = TRUE),]

rimp3 <- rimp2[c(1:10),]

rimp3$var <- c("microbial","abiotic","abiotic","microbial","abiotic","abiotic","abiotic","microbial","microbial","microbial")

rimp3$rn[rimp3$rn == "pcoa1"] <- "16S PCoA 1"
rimp3$rn[rimp3$rn == "BG"] <- "BG activity"
rimp3$rn[rimp3$rn == "pcoa2"] <- "16S PCoA 2"
rimp3$rn[rimp3$rn == "tet_scaled"] <- "tet ARG abundance"
rimp3$rn[rimp3$rn == "rk"] <- "16S r:K"
rimp3$rn[rimp3$rn == "CN"] <- "Extractable C:N"
rimp3$rn[rimp3$rn == "TDN"] <- "Total Extractable N"
rimp3$rn[rimp3$rn == "bac_shannon"] <- "16S Shannon"
rimp3$rn[rimp3$rn == "DOC"] <- "Extractable DOC"

rimp3$rn <- factor(rimp3$rn, levels = rimp3$rn[order(rimp3$X.IncMSE)])


resp_imp <- ggplot(rimp3, aes(x=X.IncMSE, y=rn, fill=var)) + 
  geom_segment(aes(yend=rn, xend=0),linetype=2)+
  geom_point(aes(fill=var), shape=21,size=5) +
  ggtitle(bquote('C mineralization:'~R^2~'= 0.56'))+
  theme_classic() +
  ylab("")+
  scale_fill_manual(values=c("white", "grey50")) +
  theme(axis.title=element_text(size=18)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Variable Importance (% MSE Increase)") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(legend.position=c(.8,.15))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))
#scale_x_continuous(limits = c(0, 15))
plot(resp_imp)


###SIR random forest

sir <- train[,c(12,1:10,24,26:29)]

set.seed(101)

r_tune <- tuneRF(train[,c(1:10,24,26:29)],train$SIR,ntreeTry=10000,stepFactor=1.5, improve=0.05,
                 trace=TRUE, plot=TRUE, doBest=FALSE)

print(r_tune)

library(randomForest)

set.seed(131)
sir.rf <- randomForest(SIR ~ ., data=sir,mtry=4,ntree=10000,
                       importance=TRUE)
print(sir.rf)

## Show "importance" of variables: higher value mean more important:
round(importance(sir.rf), 2)

varImpPlot(sir.rf,type=1)

## Show "importance" of variables: higher value mean more important:
sirimp <- data.frame(round(importance(sir.rf,type=1), 2))

library(data.table)
setDT(sirimp, keep.rownames = TRUE)

sirimp2 <- sirimp[order(sirimp$X.IncMSE,decreasing = TRUE),]

sirimp3 <- sirimp2[c(1:10),]

sirimp3$var <- c("microbial","abiotic","abiotic","abiotic","microbial","microbial","microbial","microbial","abiotic","microbial")

sirimp3$rn[sirimp3$rn == "pcoa1"] <- "16S PCoA 1"
sirimp3$rn[sirimp3$rn == "pcoa2"] <- "16S PCoA 2"
sirimp3$rn[sirimp3$rn == "bac_scaled"] <- "16S abundance"
sirimp3$rn[sirimp3$rn == "MBC"] <- "Microbial Biomass"
sirimp3$rn[sirimp3$rn == "bac_shannon"] <- "16S Shannon"
sirimp3$rn[sirimp3$rn == "pcoa2"] <- "16S PCoA 2"
sirimp3$rn[sirimp3$rn == "CN"] <- "Extractable C:N"
sirimp3$rn[sirimp3$rn == "TDN"] <- "Total Extractable N"
sirimp3$rn[sirimp3$rn == "rk"] <- "16S r:K"
sirimp3$rn[sirimp3$rn == "DNA"] <- "Microbial Biomass"
sirimp3$rn[sirimp3$rn == "DOC"] <- "Extractable DOC"


sirimp3$rn <- factor(sirimp3$rn, levels = sirimp3$rn[order(sirimp3$X.IncMSE)])

sir_imp <- ggplot(sirimp3, aes(x=X.IncMSE, y=rn, fill=var)) + 
  geom_segment(aes(yend=rn, xend=0),linetype=2)+
  geom_point(aes(fill=var), shape=21,size=5) +
  ggtitle(bquote('SIR:'~R^2~'= 0.43'))+
  theme_classic() +
  ylab("")+
  scale_fill_manual(values=c("white", "grey50")) +
  theme(axis.title=element_text(size=18)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Variable Importance (% MSE Increase)") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(legend.position=c(.8,.15))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))
#scale_x_continuous(limits = c(0, 12))
plot(sir_imp)


#qCO2 random forest

qco2 <- train[,c(13,1:10,18,24:29)]

set.seed(101)

r_tune <- tuneRF(train[,c(1:10,18,24:29)],train$qCO2,ntreeTry=10000,stepFactor=1.5, improve=0.05,
                 trace=TRUE, plot=TRUE, doBest=FALSE)

print(r_tune)

library(randomForest)

set.seed(131)
q.rf <- randomForest(qCO2 ~ ., data=qco2,ntree=10000,mtry=4,
                     importance=TRUE)
print(q.rf)

## Show "importance" of variables: higher value mean more important:
round(importance(q.rf), 2)

varImpPlot(q.rf,type=1)

## Show "importance" of variables: higher value mean more important:
qimp <- data.frame(round(importance(q.rf,type=1), 2))

library(data.table)
setDT(qimp, keep.rownames = TRUE)

qimp2 <- qimp[order(qimp$X.IncMSE,decreasing = TRUE),]

qimp3 <- qimp2[c(1:10),]

qimp3$var <- c("microbial","microbial","abiotic","microbial","microbial","microbial","microbial","abiotic","microbial","abiotic")

qimp3$rn[qimp3$rn == "pcoa1"] <- "16S PCoA 1"
qimp3$rn[qimp3$rn == "bac_scaled"] <- "16S abundance"
qimp3$rn[qimp3$rn == "MBC"] <- "Microbial Biomass"
qimp3$rn[qimp3$rn == "bac_shannon"] <- "16S Shannon"
qimp3$rn[qimp3$rn == "pcoa2"] <- "16S PCoA 2"
qimp3$rn[qimp3$rn == "CN"] <- "Extractable C:N"
qimp3$rn[qimp3$rn == "TDN"] <- "Total Extractable N"
qimp3$rn[qimp3$rn == "rk"] <- "16S r:K"
qimp3$rn[qimp3$rn == "copynumber"] <- "16S operon number"
qimp3$rn[qimp3$rn == "DOC"] <- "Extractable DOC"
qimp3$rn[qimp3$rn == "bac_richness"] <- "16S Richness"


qimp3$rn <- factor(qimp3$rn, levels = qimp3$rn[order(qimp3$X.IncMSE)])

q_imp <- ggplot(qimp3, aes(x=X.IncMSE, y=rn, fill=var)) + 
  geom_segment(aes(yend=rn, xend=0),linetype=2)+
  geom_point(aes(fill=var), shape=21,size=5) +
  ggtitle(bquote(qCO[2]~':'~R^2~'= 0.41'))+
  theme_classic() +
  ylab("")+
  scale_fill_manual(values=c("white", "grey50")) +
  theme(axis.title=element_text(size=18)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Variable Importance (% MSE Increase)") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(legend.position=c(.8,.15))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))
#scale_x_continuous(limits = c(0, 12))
plot(q_imp)

summary(lm(train$qCO2~train$bac_shannon))

plot(train$bac_shannon,train$qCO2)



#Enzyme activity random forest


enz <- train[,c(20,1:5,7:10,24,26:29)]

set.seed(101)

r_tune <- tuneRF(train[,c(1:5,7:10,24,26:29)],train$totalH,ntreeTry=10000,stepFactor=1.5, improve=0.05,
                 trace=TRUE, plot=TRUE, doBest=FALSE)

print(r_tune)

library(randomForest)

set.seed(131)
enz.rf <- randomForest(totalH ~ ., data=enz,mtry=6,ntree=10000,
                       importance=TRUE)
print(enz.rf)

## Show "importance" of variables: higher value mean more important:
round(importance(enz.rf), 2)

varImpPlot(enz.rf,type=1)

enzimp <- data.frame(round(importance(enz.rf,type=1), 2))

library(data.table)
setDT(enzimp, keep.rownames = TRUE)

enzimp2 <- enzimp[order(enzimp$X.IncMSE,decreasing = TRUE),]

enzimp3 <- enzimp2[c(1:10),]

enzimp3$var <- c("microbial","microbial","microbial","abiotic","abiotic","abiotic","microbial","microbial","microbial","abiotic")

enzimp3$rn[enzimp3$rn == "pcoa1"] <- "16S PCoA 1"
enzimp3$rn[enzimp3$rn == "bac_scaled"] <- "16S abundance"
enzimp3$rn[enzimp3$rn == "DNA"] <- "Microbial Biomass"
enzimp3$rn[enzimp3$rn == "bac_shannon"] <- "16S Shannon"
enzimp3$rn[enzimp3$rn == "pcoa2"] <- "16S PCoA 2"
enzimp3$rn[enzimp3$rn == "CN"] <- "Extractable C:N"
enzimp3$rn[enzimp3$rn == "rk"] <- "16S r:K"
enzimp3$rn[enzimp3$rn == "TDN"] <- "Total Extractable N"
enzimp3$rn[enzimp3$rn == "DOC"] <- "Extractable DOC"
enzimp3$rn[enzimp3$rn == "tet.16S"] <- "tet:16S abundance"


enzimp3$rn <- factor(enzimp3$rn, levels = enzimp3$rn[order(enzimp3$X.IncMSE)])


enz_imp <- ggplot(enzimp3, aes(x=X.IncMSE, y=rn, fill=var)) + 
  geom_segment(aes(yend=rn, xend=0),linetype=2)+
  geom_point(aes(fill=var), shape=21,size=5) +
  ggtitle(bquote('Enzyme Activity:'~R^2~'= 0.77'))+
  theme_classic() +
  ylab("")+
  scale_fill_manual(values=c("white", "grey50")) +
  theme(axis.title=element_text(size=18)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Variable Importance (% MSE Increase)") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(legend.position=c(.8,.15))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))
plot(enz_imp)

#N mineralization random forest

min <- train[,c(14,1,3,7:10,24,26:29,17,19)]

set.seed(101)

r_tune <- tuneRF(train[,c(1,3,7:10,24,26:29,17,19)],train$Nmin,ntreeTry=10000,stepFactor=1.5, improve=0.05,
                 trace=TRUE, plot=TRUE, doBest=FALSE)

print(r_tune)

library(randomForest)

set.seed(131)
n.rf <- randomForest(Nmin ~ ., data=min,mtry=6,ntree=10000,
                     importance=TRUE)
print(n.rf)

## Show "importance" of variables: higher value mean more important:
round(importance(n.rf), 2)

varImpPlot(n.rf,type=1)

## Show "importance" of variables: higher value mean more important:
nminimp <- data.frame(round(importance(n.rf,type=1), 2))

library(data.table)
setDT(nminimp, keep.rownames = TRUE)

nminimp2 <- nminimp[order(nminimp$X.IncMSE,decreasing = TRUE),]

nminimp3 <- nminimp2[c(1:10),]

nminimp3$var <- c("microbial","microbial","microbial","microbial","abiotic","microbial","microbial","microbial","microbial","abiotic")

nminimp3$rn[nminimp3$rn == "pcoa1"] <- "16S PCoA 1"
nminimp3$rn[nminimp3$rn == "bac_scaled"] <- "16S abundance"
nminimp3$rn[nminimp3$rn == "DNA"] <- "Microbial Biomass"
nminimp3$rn[nminimp3$rn == "bac_shannon"] <- "16S Shannon"
nminimp3$rn[nminimp3$rn == "pcoa2"] <- "16S PCoA 2"
nminimp3$rn[nminimp3$rn == "CN"] <- "Extractable C:N"
nminimp3$rn[nminimp3$rn == "NAG"] <- "NAG activity"
nminimp3$rn[nminimp3$rn == "rk"] <- "16S r:K"
nminimp3$rn[nminimp3$rn == "copynumber"] <- "16S operon number"
nminimp3$rn[nminimp3$rn == "LAP"] <- "LAP activity"
nminimp3$rn[nminimp3$rn == "DOC"] <- "Extractable DOC"

nminimp3$rn <- factor(nminimp3$rn, levels = nminimp3$rn[order(nminimp3$X.IncMSE)])


nmin_imp <- ggplot(nminimp3, aes(x=X.IncMSE, y=rn, fill=var)) + 
  geom_segment(aes(yend=rn, xend=0),linetype=2)+
  geom_point(aes(fill=var), shape=21,size=5) +
  ggtitle(bquote('N Mineralization:'~R^2~'= 0.91'))+
  theme_classic() +
  ylab("")+
  scale_fill_manual(values=c("white", "grey50")) +
  theme(axis.title=element_text(size=18)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Variable Importance (% MSE Increase)") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(legend.position=c(.8,.15))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))
plot(nmin_imp)

dev.off()


#nitrification random forest

nit <- train[,c(15,1,3,7:9,22,21,24,26:29)]

library(randomForest)

set.seed(101)

r_tune <- tuneRF(train[,c(1,3,7:9,21,22,24,26:29)],train$Nit,ntreeTry=10000,stepFactor=1.5, improve=0.05,
                 trace=TRUE, plot=TRUE, doBest=FALSE)

print(r_tune)


set.seed(131)
nit.rf <- randomForest(Nit ~ ., data=nit,mtry=6,ntree=10000,
                       importance=TRUE)
print(nit.rf)

## Show "importance" of variables: higher value mean more important:
nit.imp <- data.frame(round(importance(nit.rf), 2))

varImpPlot(nit.rf)

## Show "importance" of variables: higher value mean more important:
nitimp <- data.frame(round(importance(nit.rf,type=1), 2))

library(data.table)
setDT(nitimp, keep.rownames = TRUE)

nitimp2 <- nitimp[order(nitimp$X.IncMSE,decreasing = TRUE),]

nitimp3 <- nitimp2[c(1:10),]

nitimp3$var <- c("microbial","microbial","abiotic","microbial","microbial","microbial","microbial","microbial","abiotic","microbial")

nitimp3$rn[nitimp3$rn == "pcoa1"] <- "16S PCoA 1"
nitimp3$rn[nitimp3$rn == "bac_scaled"] <- "16S abundance"
nitimp3$rn[nitimp3$rn == "DNA"] <- "Microbial Biomass"
nitimp3$rn[nitimp3$rn == "bac_shannon"] <- "16S Shannon"
nitimp3$rn[nitimp3$rn == "tet_scaled"] <- "tet ARG abundance"
nitimp3$rn[nitimp3$rn == "CN"] <- "Extractable C:N"
nitimp3$rn[nitimp3$rn == "AOB_scaled"] <- "AOB abundance"
nitimp3$rn[nitimp3$rn == "rk"] <- "16S r:K"
nitimp3$rn[nitimp3$rn == "copynumber"] <- "16S operon number"
nitimp3$rn[nitimp3$rn == "AOA_scaled"] <- "AOA abundance"
nitimp3$rn[nitimp3$rn == "DOC"] <- "Extractable DOC"
nitimp3$rn[nitimp3$rn == "pcoa2"] <- "16S PCoA 2"

nitimp3$rn <- factor(nitimp3$rn, levels = nitimp3$rn[order(nitimp3$X.IncMSE)])


nit_imp <- ggplot(nitimp3, aes(x=X.IncMSE, y=rn, fill=var)) + 
  geom_segment(aes(yend=rn, xend=0),linetype=2)+
  geom_point(aes(fill=var), shape=21,size=5) +
  ggtitle(bquote('Nitrification:'~R^2~'= 0.96'))+
  theme_classic() +
  ylab("")+
  scale_fill_manual(values=c("white", "grey50")) +
  theme(axis.title=element_text(size=18)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Variable Importance (% MSE Increase)") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(legend.position=c(.8,.15))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))
plot(nit_imp)


#Multifunctionality random forest


multi <- train[,c(23,1,3,7:9,24,26:29)]

set.seed(131)

r_tune <- tuneRF(train[,c(1,3,7:9,24,26:29)],train$multi,ntreeTry=10000,stepFactor=1.5, improve=0.05,
                 trace=TRUE, plot=TRUE, doBest=FALSE)

print(r_tune)

library(randomForest)

set.seed(131)
m.rf <- randomForest(multi ~ ., data=multi,mtry=3,ntree=10000,
                     importance=TRUE)
print(m.rf)

## Show "importance" of variables: higher value mean more important:
round(importance(m.rf), 2)

varImpPlot(m.rf,type=1)

mimp <- data.frame(round(importance(m.rf,type=1), 2))

library(data.table)
setDT(mimp, keep.rownames = TRUE)

mimp2 <- mimp[order(mimp$X.IncMSE,decreasing = TRUE),]

mimp3 <- mimp2[c(1:9),]

mimp3$var <- c("microbial","microbial","microbial","abiotic","microbial","abiotic","microbial","microbial","microbial")

mimp3$rn[mimp3$rn == "pcoa1"] <- "16S PCoA 1"
mimp3$rn[mimp3$rn == "bac_scaled"] <- "16S abundance"
mimp3$rn[mimp3$rn == "MBC"] <- "Microbial Biomass"
mimp3$rn[mimp3$rn == "bac_shannon"] <- "16S Shannon"
mimp3$rn[mimp3$rn == "pcoa2"] <- "16S PCoA 2"
mimp3$rn[mimp3$rn == "CN"] <- "Extractable C:N"
mimp3$rn[mimp3$rn == "TDN"] <- "Total Extractable N"
mimp3$rn[mimp3$rn == "rk"] <- "16S r:K"
mimp3$rn[mimp3$rn == "DNA"] <- "Microbial Biomass"
mimp3$rn[mimp3$rn == "DOC"] <- "Extractable DOC"
mimp3$rn[mimp3$rn == "tet_scaled"] <- "tet ARG abundance"
mimp3$rn[mimp3$rn == "copynumber"] <- "16S operon number"

mimp3$rn <- factor(mimp3$rn, levels = mimp3$rn[order(mimp3$X.IncMSE)])

m_imp <- ggplot(mimp3, aes(x=X.IncMSE, y=rn, fill=var)) + 
  geom_segment(aes(yend=rn, xend=0),linetype=2)+
  geom_point(aes(fill=var), shape=21,size=5) +
  ggtitle(bquote('Multifunctionality:'~R^2~'= 0.83'))+
  theme_classic() +
  ylab("")+
  scale_fill_manual(values=c("white", "grey50")) +
  theme(axis.title=element_text(size=18)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Variable Importance (% MSE Increase)") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  theme(legend.position=c(.8,.15))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))
#scale_x_continuous(limits = c(0, 12))
plot(m_imp)


summary(lm(train$multi~train$bac_shannon))

plot(train$bac_shannon,train$multi)



jpeg(filename="Figure4.jpeg", bg="transparent", res=600, units = "in", height=18, width=14) 

f4 <- plot_grid(resp_imp,sir_imp,q_imp,enz_imp,nmin_imp,nit_imp, ncol = 2, align="hv", label_size=35,labels=c('A', 'B','C','D','E','F'))

f4

dev.off()



jpeg(filename="Figure5.jpeg", bg="transparent", res=500, units = "in", height=6, width=14) 

f5 <- plot_grid(ord_multi, m_imp, ncol = 2, align="hv", label_size=25,labels=c('A', 'B'))

f5

dev.off()


#Supplementary stuff


#AOA abundance


m2 <- lm(AOA_scaled~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(AOA_scaled~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m3)

aggregate(AOA_scaled~Stress+Diversity, data=s2,FUN=mean)

em1 <- emmeans(m3, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_aoa <- aggregate(AOA_scaled ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_aoa <- do.call(data.frame, sum_aoa)
sum_aoa

sum_aoa$se <- sum_aoa$AOA_scaled.sd / sqrt(sum_aoa$AOA_scaled.n)
head(sum_aoa)


sum_aoa$Stress <- factor(sum_aoa$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

aoa<- ggplot(sum_aoa, aes(x=Diversity, y=AOA_scaled.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=AOA_scaled.mean-se, ymax=AOA_scaled.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab("AOA abundance") +
  facet_wrap(~ LandUse)
plot(aoa)

#AOB abundance


m2 <- lm(AOB_scaled~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(AOB_scaled~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m3)

aggregate(AOB_scaled~Stress+Diversity, data=s2,FUN=mean)

em1 <- emmeans(m3, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_AOB <- aggregate(AOB_scaled ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_AOB <- do.call(data.frame, sum_AOB)
sum_AOB

sum_AOB$se <- sum_AOB$AOB_scaled.sd / sqrt(sum_AOB$AOB_scaled.n)
head(sum_AOB)


sum_AOB$Stress <- factor(sum_AOB$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

AOB<- ggplot(sum_AOB, aes(x=Diversity, y=AOB_scaled.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=AOB_scaled.mean-se, ymax=AOB_scaled.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=20)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab("AOB abundance") +
  facet_wrap(~ LandUse)
plot(AOB)


#tetW copies

dna <- read.csv("C:/Users/ernie/OneDrive/Desktop/Current Projects/RDE/RDE_dna.csv")

dna2 <-  dna[dna$Diversity!="Pos",]

dna2 <- dna2[dna2$Diversity!="Neg",]

library(lme4)
m2 <- lm(tetW_scaled~Stress*LandUse*Diversity,data=dna2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(tetW_scaled+.01~Stress*LandUse*Diversity,data=dna2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m2)

em1 <- emmeans(m2, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_tetW <- aggregate(tetW_scaled ~ LandUse+Diversity+Stress, data=dna2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_tetW <- do.call(data.frame, sum_tetW)
sum_tetW

sum_tetW$se <- sum_tetW$tetW_scaled.sd / sqrt(sum_tetW$tetW_scaled.n)
head(sum)

sum_tetW$Stress <- factor(sum_tetW$Stress,levels=c("Control","Antibiotics"))


#Now for some plotting

library(ggplot2)
library(RColorBrewer)

line_plotCmin <- ggplot(sum_tetW, aes(x=Diversity, y=tetW_scaled.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=tetW_scaled.mean-se, ymax=tetW_scaled.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('tetW gene copies/g soil')) +
  facet_wrap(~ LandUse)
plot(line_plotCmin)


#tetM copies


library(lme4)
m2 <- lm(tetM_scaled~Stress*LandUse*Diversity,data=dna2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(tetM_scaled+.01~Stress*LandUse*Diversity,data=dna2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m2)

em1 <- emmeans(m2, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_tetM <- aggregate(tetM_scaled ~ LandUse+Diversity+Stress, data=dna2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_tetM <- do.call(data.frame, sum_tetM)
sum_tetM

sum_tetM$se <- sum_tetM$tetM_scaled.sd / sqrt(sum_tetM$tetM_scaled.n)
head(sum)

sum_tetM$Stress <- factor(sum_tetM$Stress,levels=c("Control","Antibiotics"))


#Now for some plotting

library(ggplot2)
library(RColorBrewer)

line_plotCmin <- ggplot(sum_tetM, aes(x=Diversity, y=tetM_scaled.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=tetM_scaled.mean-se, ymax=tetM_scaled.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('tetM gene copies/g soil')) +
  facet_wrap(~ LandUse)
plot(line_plotCmin)

#MBC

m2 <- lm(MBC~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(MBC~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m3)

em1 <- emmeans(m3, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_mbc <- aggregate(MBC ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_mbc <- do.call(data.frame, sum_mbc)
sum_mbc

sum_mbc$se <- sum_mbc$MBC.sd / sqrt(sum_mbc$MBC.n)
head(sum_mbc)

sum_mbc$Stress <- factor(sum_mbc$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

line_plotCmin <- ggplot(sum_mbc, aes(x=Diversity, y=MBC.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=MBC.mean-se, ymax=MBC.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('Microbial Biomass C')) +
  facet_wrap(~ LandUse)
plot(line_plotCmin)

#DNA


m2 <- lm(DNA~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

Anova(m2)

em1 <- emmeans(m2, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_dna <- aggregate(DNA ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_dna <- do.call(data.frame, sum_dna)
sum_dna

sum_dna$se <- sum_dna$DNA.sd / sqrt(sum_dna$DNA.n)
head(sum)

sum_dna$Stress <- factor(sum_dna$Stress,levels=c("Control","Antibiotic"))


#Now for some plotting

library(ggplot2)
library(RColorBrewer)

line_plotCmin <- ggplot(sum_dna, aes(x=Diversity, y=DNA.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=DNA.mean-se, ymax=DNA.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('DNA concentration (ng/g soil)')) +
  facet_wrap(~ LandUse)
plot(line_plotCmin)

#AP

m2 <- lm(AP~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(AP~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m3)

em1 <- emmeans(m3, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_ap <- aggregate(AP ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_ap <- do.call(data.frame, sum_ap)
sum_ap

sum_ap$se <- sum_ap$AP.sd / sqrt(sum_ap$AP.n)
head(sum_ap)

sum_ap$Stress <- factor(sum_ap$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

line_plotCmin <- ggplot(sum_ap, aes(x=Diversity, y=AP.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=AP.mean-se, ymax=AP.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('AP activity')) +
  facet_wrap(~ LandUse)
plot(line_plotCmin)

#BG

m2 <- lm(BG~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(BG~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m3)

em1 <- emmeans(m3, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_bg <- aggregate(BG ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_bg <- do.call(data.frame, sum_bg)
sum_bg

sum_bg$se <- sum_bg$BG.sd / sqrt(sum_bg$BG.n)
head(sum_bg)


sum_bg$Stress <- factor(sum_bg$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

line_plotCmin <- ggplot(sum_bg, aes(x=Diversity, y=BG.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=BG.mean-se, ymax=BG.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('BG activity')) +
  facet_wrap(~ LandUse)
plot(line_plotCmin)

#NAG

m2 <- lm(NAG~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(NAG~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m3)

em1 <- emmeans(m3, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_nag <- aggregate(NAG ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_nag <- do.call(data.frame, sum_nag)
sum_nag

sum_nag$se <- sum_nag$NAG.sd / sqrt(sum_nag$NAG.n)
head(sum_nag)

sum_nag$Stress <- factor(sum_nag$Stress, levels=c("Control","Antibiotic"))


#Now for some plotting

library(ggplot2)
library(RColorBrewer)

line_plotCmin <- ggplot(sum_nag, aes(x=Diversity, y=NAG.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=NAG.mean-se, ymax=NAG.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('NAG activity')) +
  facet_wrap(~ LandUse)
plot(line_plotCmin)


#LAP

m2 <- lm(LAP~Stress*LandUse*Diversity,data=s2)

shapiro.test(resid(m2))

library(car)

qqPlot(resid(m2))

m3 <- glm(LAP~Stress*LandUse*Diversity,data=s2,family=Gamma(link=log))

AIC(m2,m3)

Anova(m3)

em1 <- emmeans(m2, ~ Stress|LandUse*Diversity)
contrast(em1, method = "pairwise")

sum_lap <- aggregate(LAP ~ LandUse+Diversity+Stress, data=s2, FUN =function(x) c(mean=mean(x), sd=sd(x),n=length(x)))

sum_lap <- do.call(data.frame, sum_lap)
sum_lap

sum_lap$se <- sum_lap$LAP.sd / sqrt(sum_lap$LAP.n)
head(sum_lap)

sum_lap$Stress <- factor(sum_lap$Stress, levels=c("Control","Antibiotic"))

#Now for some plotting

library(ggplot2)
library(RColorBrewer)

line_plotCmin <- ggplot(sum_lap, aes(x=Diversity, y=LAP.mean, color=Diversity, shape=Stress)) + 
  geom_point(position=position_dodge(width=0.9),aes(shape=Stress,color=Diversity), size=3) +
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=LAP.mean-se, ymax=LAP.mean+se), width=.5, size=.75) +
  labs(list(x ="")) +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(axis.title=element_text(size=20)) +
  theme(text=element_text(size=15)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Diversity Treatment") +
  theme(axis.line.x=element_line(colour="black", size=1)) + 
  theme(axis.line.y=element_line(colour="black", size=1)) + 
  theme(legend.title=element_blank())+
  ylab(bquote('LAP activity')) +
  facet_wrap(~ LandUse)
plot(line_plotCmin)
