#########################################################
#### normalization tests SNL-plankton CSIA
#########################################################

# make the PCA df-normalized, remove proportion plant, keeping sources
# PCA.df above has sources in, removed later in 'PCA.df.trim'
# make ID column to run the normalization

PCA.df.normAA<-PCA.df

#make ID column to run the normalization
PCA.df.normAA$ID<-1:nrow(PCA.df.normAA)
for(i in 1:length(PCA.df.normAA$ID)){
  PCA.df.normAA$Ile.n[i] <- (PCA.df.normAA$Ile[i]-mean(as.numeric(PCA.df.normAA[i,3:7])))
  PCA.df.normAA$Leu.n[i] <- (PCA.df.normAA$Leu[i]-mean(as.numeric(PCA.df.normAA[i,3:7])))
  PCA.df.normAA$Phe.n[i] <- (PCA.df.normAA$Phe[i]-mean(as.numeric(PCA.df.normAA[i,3:7])))
  PCA.df.normAA$Thr.n[i] <- (PCA.df.normAA$Thr[i]-mean(as.numeric(PCA.df.normAA[i,3:7])))
  PCA.df.normAA$Val.n[i] <- (PCA.df.normAA$Val[i]-mean(as.numeric(PCA.df.normAA[i,3:7])))
}

# original data for LDA: (1) raw producers 'prod' df
head(prod)
head(PCA.df.normAA)

LDA.df.normAA<-PCA.df.normAA %>% 
  dplyr::select(Lake, phy.group, Ile.n, Leu.n, Phe.n, Thr.n, Val.n)

#subset for producers and zoops
# all sources
LDA.prod.norm<-LDA.df.normAA[(LDA.df.normAA$Lake=="Source"),]
LDA.prod.norm<-LDA.prod.norm %>% dplyr::select(-Lake) # drop lake
LDA.prod.norm<-droplevels(LDA.prod.norm) # remove dropped levels

# all zoops
LDA.zoop.norm<-LDA.df.normAA[!(LDA.df.normAA$Lake=="Source"),]
LDA.zoop.norm<-LDA.zoop.norm %>% dplyr::select(-Lake) # drop lake
LDA.zoop.norm<-droplevels(LDA.zoop.norm) # remove dropped levels


#### Run the LDA
#Run an LDA with a jackknifing model fit to look at error rate
#'CV = TRUE' makes the LDA run a jackknifed - leave one out - model fit
All.lda.n <- lda(phy.group ~ Ile.n + Leu.n + Phe.n + Thr.n + Val.n, data = LDA.prod.norm, CV = TRUE)

#Create a table which compares the classification from the LDA model to the actual classification
All.reclass.n <- table(LDA.prod.norm$phy.group, All.lda.n$class)

#Total percent of samples correctly classified is the sum of the diagonal of this table
sum(diag(prop.table(All.reclass.n)))

#Percent of each producer group correctly reclassified
diag(prop.table(All.reclass.n, 1))

#Create a training LDA function from the library data
#Note - you can't use the 'All.lda' object above because the 'CV = TRUE' command was used to create it, and for some reason this won't work with the predict() function
All.train.n <- lda(phy.group~ Ile.n + Leu.n + Phe.n + Thr.n + Val.n, data = LDA.prod.norm)
All.train.n

# coefficients of LD1
lda.info.n <- as.data.frame(All.train.n$scaling)

#Create a data frame with these LDA coordinates
AlldataPredict.n <- data.frame(Group = LDA.prod.norm$phy.group, SampleID = prod$SampleID, predict(All.train.n)$x)

####### 
#Classify zooplankton
zoop.predict.n <- predict(object = All.train.n, newdata = LDA.zoop.norm)
zoop.predict.data.n <- data.frame(SampleID = zoops$SampleID, Type = zoops$Type, Lake = zoops$Lake, zoop.predict.n)


##############
# combine the zooplankton and source LDA
##############

####### Producers
# AlldataPredict.mod original columns as - "Group", "SampleID", "LD1"
# rearrange
AlldataPredict.mod.n<- AlldataPredict.n %>%
  dplyr::select("SampleID", "Group", "LD1")

# rename columns
colnames(AlldataPredict.mod.n)<-c("SampleID", "orig.class", "LD1")

# make new columns to allow merge
AlldataPredict.mod.n$Lake<-"Source"
AlldataPredict.mod.n$LD.class<-"Source"

# rearrange
AlldataPredict.mod.n<- AlldataPredict.mod.n %>%
  dplyr::select("SampleID", "orig.class", "Lake", "LD.class", "LD1")


############# Zooplankton
# grab columns we want
zoop.predict.data.NMX.n<- zoop.predict.data.n %>%
  dplyr::select("SampleID", "Type", "Lake", "class", "LD1")

# rename factors
zoop.predict.data.NMX.n <- zoop.predict.data.NMX.n %>% 
  dplyr::rename("orig.class" = "Type")

zoop.predict.data.NMX.n <- zoop.predict.data.NMX.n %>% 
  dplyr::rename("LD.class" = "class")


####### bind the 2 dfs
LDA.df.n<-rbind(zoop.predict.data.NMX.n, AlldataPredict.mod.n)

# make new level for "phy.group" based on taxonomy
LDA.df.n$phy.group<- ifelse(LDA.df.n$orig.class == "calanoid", "copepoda",
                          ifelse(LDA.df.n$orig.class == "calanoid_cyclopoid", "copepoda",
                                 ifelse(LDA.df.n$orig.class == "Ceriodaphnia", "cladocera",
                                        ifelse(LDA.df.n$orig.class == "Daphnia", "cladocera",
                                               ifelse(LDA.df.n$orig.class == "Holopedium", "cladocera",
                                                      ifelse(LDA.df.n$orig.class == "Ceriodaphnia", "cladocera",
                                                             ifelse(LDA.df.n$orig.class == ">350um", ">350um",
                                                                    ifelse(LDA.df.n$orig.class == "Plants", "Plants",
                                                                           "POM"))))))))

# make new level for "Elevation" (in m), "0" for ordering, but actually is "NA"
LDA.df.n$Elevation<- ifelse(LDA.df.n$Lake == "Lukens", "2513",
                          ifelse(LDA.df.n$Lake == "MayPond", "2714",
                                 ifelse(LDA.df.n$Lake == "Sunrise2", "2830",
                                        ifelse(LDA.df.n$Lake == "Blue", "3054",
                                               ifelse(LDA.df.n$Lake == "Greenstone", "3091",
                                                      ifelse(LDA.df.n$Lake == "Cascade", "3140",
                                                             ifelse(LDA.df.n$Lake == "EasternBrook", "3155",
                                                                    ifelse(LDA.df.n$Lake == "Granite2", "3178",
                                                                           "0"))))))))


# rearrange
LDA.df.n<- LDA.df.n %>%
  dplyr::select("SampleID", "Lake", "Elevation", "orig.class", "phy.group", "LD.class", "LD1")

write.csv(LDA.df.n, "output/ProducerNORM_zoop_LDA_score.csv", row.names=FALSE)

#########################
#### make the plots
##### can load in from output '("output/Producer_zoop_LDA_score.csv")'

LDA.df.n$Elevation<-as.numeric(LDA.df.n$Elevation)

# reorder lake by relative elevation
LDA.df.n$Lake<-reorder(LDA.df.n$Lake, LDA.df.n$Elevation)

# reorder factors 
LDA.df.n$orig.class<-factor(LDA.df.n$orig.class, levels=c(">350um", "calanoid", 
                                                      "calanoid_cyclopoid", "Ceriodaphnia",
                                                      "Daphnia", "Holopedium", "POM", "Plants"))

phy.lda.col.n<-c("gray65", "firebrick1","firebrick1", "goldenrod1", "goldenrod1","goldenrod1", "dodgerblue","springgreen4")
phy.lda.pts.n<-c(21, 21,24, 21,24,22, 21,21)

###
LDA.spp.lake.n<-ggplot(LDA.df.n, aes(x=LD1, y=Lake, 
                                 fill = orig.class, shape= orig.class))+
  geom_point(size=3)+
  geom_vline(xintercept=0.4, linetype="dashed", color = "gray60")+
  annotate(geom="text", label="Plant-Supported", x=-3.5, y=9.5, color="springgreen4", size=3) +
  annotate(geom="text", label="POM-Supported", x=2.5, y=9.5, color="dodgerblue", size=3) +
  scale_fill_manual(values=phy.lda.col)+
  scale_shape_manual(values=phy.lda.pts)+
  ylab("Site: low-to-high elevation") +
  xlab("LD1") +
  theme_classic()

LDA.spp.lake.n
dev.copy(pdf, "figures/Fig.S3.LDA.spp.lakeNORM.pdf", width = 8, height = 6)
dev.off()

####################
#LDA box plot by phy.group
phy.5.colors.n<-c("dodgerblue", "gray85", "orange","coral", "springgreen4")

LDA.df.n$phy.group<-factor(LDA.df.n$phy.group, levels=c("POM", ">350um", "cladocera", 
                                                    "copepoda",
                                                    "Plants"))

LDA.phy.boxplot.n<-ggplot(LDA.df.n, aes(x=LD1, y=phy.group, 
                                    fill = phy.group))+
  geom_boxplot(alpha=0.5)+
  geom_jitter(aes(fill=phy.group), colour="black",pch=21, size=2) +
  geom_vline(xintercept=0.4, linetype="dashed", color = "gray60")+
  annotate(geom="text", label="Plant-Supported", x=-3.5, y=5.5, color="springgreen4", size=3) +
  annotate(geom="text", label="POM-Supported", x=2.5, y=5.5, color="dodgerblue", size=3) +
  scale_fill_manual(values=phy.5.colors)+
  scale_color_manual(values=phy.5.colors)+
  ylab("Consumer or Source") +
  xlab("LD1") +
  theme_classic()

LDA.phy.boxplot.n
dev.copy(pdf, "figures/Figx.LDA.phy.boxplotNORM.pdf", width = 6, height = 6)
dev.off()


############################################################################
############################################################################
# Larsen 2013 data for microbes, plants, fungi, Daphnia
############################################################################
############################################################################
lars.df<-read.csv("data/Larsen data/Larsen_EAAs.csv")


#make ID column to run the normalization
lars.df$SampleID<-1:nrow(lars.df)
for(i in 1:length(lars.df$Sample.type)){
  lars.df$Ile.n[i] <- (lars.df$Ile[i]-mean(as.numeric(lars.df[i,4:8])))
  lars.df$Leu.n[i] <- (lars.df$Leu[i]-mean(as.numeric(lars.df[i,4:8])))
  lars.df$Phe.n[i] <- (lars.df$Phe[i]-mean(as.numeric(lars.df[i,4:8])))
  lars.df$Thr.n[i] <- (lars.df$Thr[i]-mean(as.numeric(lars.df[i,4:8])))
  lars.df$Val.n[i] <- (lars.df$Val[i]-mean(as.numeric(lars.df[i,4:8])))
}

lars.df.norm<- lars.df %>% 
  dplyr::select(Data.origin, Sample.type, Ile.n, Leu.n, Phe.n, Thr.n, Val.n)

# rename column
lars.df.norm<- lars.df.norm %>% 
  dplyr::rename("phy.group" = "Sample.type")

#remove Seston abnd Macroalgae
lars.df.norm<-lars.df.norm[(!lars.df.norm$phy.group=="Seston" & !lars.df.norm$phy.group=="Macroalgae"),]

# now contains mean normalized values for: Bacteria, Fungi, Microalgae, Soil, Daphnia
# seaprate by producers and zoops

lars.norm.prod<-lars.df.norm[!(lars.df.norm$phy.group=="Daphnia"),]
lars.norm.Daph<-lars.df.norm[(lars.df.norm$phy.group=="Daphnia"),]

# pull data together and merge
head(LDA.prod.norm)
LDA.prod.norm$Data.origin<-"Wall.Besser"
LDA.zoop.norm$Data.origin<-"Wall.Besser"

# rearrange
LDA.prod.norm<- LDA.prod.norm %>% 
  dplyr::select(Data.origin, phy.group, Ile.n, Leu.n, Phe.n, Thr.n, Val.n)

LDA.zoop.norm<- LDA.zoop.norm %>% 
  dplyr::select(Data.origin, phy.group, Ile.n, Leu.n, Phe.n, Thr.n, Val.n)

# ready to merge as producers and zoops as separate dfs
Pooled.prod.EAA.n<-rbind(LDA.prod.norm, lars.norm.prod)
Pooled.zoop.EAA.n<-rbind(LDA.zoop.norm, lars.norm.Daph)

# merge all of it in one df
All.norm.pooled<-rbind(Pooled.prod.EAA.n,Pooled.zoop.EAA.n)


#### #### #### #### #### #### #### #### #### 
#### lests check out the PCA first
# the response variables
PCA.pool.norm<- All.norm.pooled %>%
  dplyr::select(Ile.n, Leu.n, Phe.n, Thr.n, Val.n)

# the factors 
PCA.pool.norm.fac<- All.norm.pooled %>%
  dplyr::select(Data.origin, phy.group)

# run the PCA on scaled and centered data
PC.pooled.norm<- prcomp(PCA.pool.norm, center = TRUE, scale= TRUE) 

PC.pooled.norm.summ<-summary(PC.pooled.norm)
#plot(PC.plank, type="lines", main="PC.area eigenvalues")



###### plot for PCA by Lake
#LDA box plot by phy.group
phy.10.colors<-c("chartreuse4", "dodgerblue", "gold2", "pink2", "cyan3", "lightsalmon2","gray65", "darkorange",  "firebrick2", "darkgoldenrod3")

## PC1 and PC2
PCA.pooled.norm <- ggbiplot(PC.pooled.norm, choices = 1:2, obs.scale = 1, var.scale = 1, 
                          groups=PCA.pool.norm.fac$phy.group,
                          ellipse = TRUE, circle = FALSE, alpha=0, ellipse.prob=0.70) + #alpha = zero makes points clear
  geom_point(aes(colour=PCA.pool.norm.fac$phy.group, shape=PCA.pool.norm.fac$Data.origin), size = 2) +
  scale_x_continuous(breaks=pretty_breaks(n=5))+
  scale_color_manual(values=phy.10.colors)+
  #scale_shape_manual(values=c(16,1,22,2,3))+
  theme(axis.ticks.length=unit(-0.25, "cm"), axis.text.y=element_text(margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")), axis.text.x=element_text(margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))) +
  ggtitle("ESS-norm")+
  theme_classic()+
  theme(legend.text=element_text(size=10), 
        panel.background = element_rect(colour = "black", linewidth=1),
        element_blank(), aspect.ratio=0.8, axis.ticks.length=unit(-0.25, "cm"),
        axis.text.y=element_text(margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
        axis.text.x=element_text(margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")))

### export it
pdf(file= "figures/Figx.PCA.norm.POOL.pdf", height=7, width=8)
plot_grid(PCA.pooled.norm)
dev.off()
