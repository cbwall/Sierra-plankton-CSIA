############################################################################
############################################################################
# Pooled food sources: algae, plants, microbes, fungi
############################################################################
############################################################################

# data from present study + microbes from Larsen et al (2009) and Besser et al (2025)

bact.source.df<-read.csv("data/producers lit/SNL_wMicrobes.csv")

make.numeric <- c("Ile", "Leu", "Phe", "Thr", "Val")
bact.source.df[make.numeric] <- lapply(bact.source.df[make.numeric], as.numeric)

#make ID column to run the normalization
bact.source.df$SampleID<-1:nrow(bact.source.df)
for(i in 1:length(bact.source.df$SampleID)){
  bact.source.df$Ile.n[i] <- (bact.source.df$Ile[i]-mean(as.numeric(bact.source.df[i,10:14])))
  bact.source.df$Leu.n[i] <- (bact.source.df$Leu[i]-mean(as.numeric(bact.source.df[i,10:14])))
  bact.source.df$Phe.n[i] <- (bact.source.df$Phe[i]-mean(as.numeric(bact.source.df[i,10:14])))
  bact.source.df$Thr.n[i] <- (bact.source.df$Thr[i]-mean(as.numeric(bact.source.df[i,10:14])))
  bact.source.df$Val.n[i] <- (bact.source.df$Val[i]-mean(as.numeric(bact.source.df[i,10:14])))
}

bact.source.df.norm<- bact.source.df %>% 
  dplyr::select(Source, Lake, Sample.names, Group, group4, Ile.n, Leu.n, Phe.n, Thr.n, Val.n)

bact.source.df.meas<- bact.source.df %>% 
  dplyr::select(Source, Lake, Sample.names, Group, group4, Ile, Leu, Phe, Thr, Val)


##############
# LDA
bact.source.df.meas<- bact.source.df.meas %>% 
  dplyr::mutate(Group = fct_recode(Group, "Benthic algae" = "Algae", "Terrestrial plants" = "Terrestrial Plants"))

bact.source.df.meas$Group<-factor(bact.source.df.meas$Group,
                                   levels=c("POM", "Benthic algae", "Terrestrial plants", "Bacteria", "Zooplankton"))

#bact.source.df.meas$pool.sourc<- ifelse(bact.source.df.meas$Group == "Terrestrial plants", "Terrestrial plants",
#                                        ifelse(bact.source.df.meas$Group == "Bacteria", "Bacteria",
#                                               ifelse(bact.source.df.meas$Group == "Zooplankton", "Zooplankton",
#                                               "POM-Alg")))
                                        
prod.bact <- bact.source.df.meas[!(bact.source.df.meas$Group == "Zooplankton"),] # remove zooplankton
prod.bact$Group<-droplevels(prod.bact$Group) # drop Zoop level

#Run an LDA with a jackknifing model fit to look at error rate
#'CV = TRUE' makes the LDA run a jackknifed - leave one out - model fit
All.lda.bact <- lda(Group ~ Ile + Leu + Phe + Thr + Val, data = prod.bact, CV = TRUE)

#Create a table which compares the classification from the LDA model to the actual classification
All.reclass.bact <- table(prod.bact$Group, All.lda.bact$class)

#Total percent of samples correctly classified is the sum of the diagonal of this table: 85%
sum(diag(prop.table(All.reclass.bact)))

mean(All.lda.bact$class == prod.bact$Group) # mean accuracy

#Percent of each producer Group correctly reclassified: 75% POM, 78% benthic algae, 80% plants, 91% bacteria
diag(prop.table(All.reclass.bact, 1))

#Create a training LDA function from the library data
#Note - you can't use the 'All.lda' object above because the 'CV = TRUE' command was used to create it, and for some reason this won't work with the predict() function

All.train.bact <- lda(Group~ Ile + Leu + Phe + Thr + Val, data = prod.bact)
All.train.bact

#plot LD1 and LD2 to show the sources
Source.LDA.bact <- ggord(All.train.bact, prod.bact$Group, arrow=0,
                    txt = NULL, veclsz = 0, vectyp="blank",# remove vectors and arrows
                    cols = c("dodgerblue","springgreen4", "lightsalmon3", "goldenrod"), grp_title = "Producers",
                    ellipse_pro = 0.95, xlim = c(-8,6), ylim = c(-5,5)) + theme_classic() + theme(legend.position = 'top')


#Create a data frame with these LDA coordinates
AllZoopPredict.bact <- data.frame(Source=prod.bact$Source, 
                                  Sample.names = prod.bact$Sample.names,
                                  Lake = prod.bact$Lake,
                                  Group = prod.bact$Group, 
                                  group4 = prod.bact$group4,
                                  class="Producer", 
                                  predict(All.train.bact)$x)

############### Now clasiify Zooplankton
Zoop.bact <- bact.source.df.meas[(bact.source.df.meas$Group == "Zooplankton"),]
Zoop.bact$Group<-droplevels(Zoop.bact$Group)

#Classify Zoop
Zoop.predict <- predict(object = All.train.bact, newdata = Zoop.bact)
Zoop.predict.data <- data.frame(Source = Zoop.bact$Source, 
                                Sample.names = Zoop.bact$Sample.names,
                                Lake = Zoop.bact$Lake,
                                Group = Zoop.bact$Group, 
                                group4 = Zoop.bact$group4, 
                                class= Zoop.predict$class, 
                                Zoop.predict$x)

Zoop.All.reclass <- table(Zoop.bact$Group, Zoop.predict.data$class) 
# Zooplankton =  23 POM, 4 benthic algae, 0 terrestrial plants, 1 as bacteria

####### plot consumers
# Producers == 'AllZoopPredict.bact'
# Zoop == `Zoop.predict.data`

####### bind the producer and zooplankton data frames
Bact.zoop.LDA.df<-rbind(AllZoopPredict.bact, Zoop.predict.data)
Bact.zoop.LDA.df$Group<-factor(Bact.zoop.LDA.df$Group,
                           levels=c("POM", "Benthic algae", "Terrestrial plants", "Bacteria", "Zooplankton"))

# colors for 5 groups
lda.5.colors.bact<-c("dodgerblue","springgreen4", "lightsalmon3", "goldenrod", "coral")

# version 1 with pooled plankton 
LDA.bact.Group <- ggplot(Bact.zoop.LDA.df, aes(x = LD1, y = LD2))+
  geom_point(size = 3.5, color="black",
             aes(fill = Group, color = Group, shape = Group), alpha = 0.9) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 0.5, aes(lty = Group, color = Group))+
  scale_color_manual(values = lda.5.colors.bact)+
  scale_fill_manual(values = lda.5.colors.bact)+
  scale_shape_manual(values = c(24, 24, 25, 23, 21))+
  scale_linetype_manual(values = c(2,2,2,2, 0))+
  xlab("LD1 (60.0%)") + ylab("LD2 (26.0%)") +
  theme_classic() + guides(lty = "none") 

LDA.bact.Group
dev.copy(pdf, "figures/Figx.LDA.EAA.bact_biplot.pdf", height = 6, width = 6)
dev.off()

write.csv(Bact.zoop.LDA.df, "output/Bact.zoop.LDA.df.csv")
