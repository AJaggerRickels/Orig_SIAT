require(profileR);require(ggplot2);require(Hmisc);require(reshape2);require(dplyr)
require(plotrix);require(lm.beta);require(corrplot);require(expss);require(gtsummary)
require(mediation);require(multilevel);require(bda);require(gvlma);require(lmSupport)
require(ggpubr);require(cowplot);require(gridExtra);require(car);require(aod)
require(reghelper);require(ggprism);require(pBrackets);require(grid);library(rstatix);
require(crosstable);require(MatchIt)
require(xlsx);require(optmatch);require(Rglpk)


library(lsr)


library(readxl)
set.seed(1) 
#ballab laptop
#setwd('/Users/ballab/Library/CloudStorage/OneDrive-Personal/Active/Research/RP010_SIAT_tfMRI_JaggerRickels/Databases/CurrentData')
#desktop mac
setwd('/Users/audreyanarickels/Library/CloudStorage/OneDrive-Personal/Active/Research/RP010_SIAT_tfMRI_JaggerRickels/Databases/CurrentData/')
#personal computer
#setwd('/Users/acjag/OneDrive/Active/Research/RP010_SIAT_tfMRI_JaggerRickels/Databases/CurrentData')

TPall = read.csv('FinalDataFile.csv', header = TRUE)
BrainData_clust = read.csv('Cluster_coef_CongF_p1neg7_clust38mean_tbl.csv',header = TRUE)
BrainData1_wm = read.csv('WholeBrain_coef_CongF_p1neg7_clust38mean_tbl.csv',header = TRUE)
BrainData2_wm = read.csv('WholeBrain_coef_CongBYwordF_p1neg3_clust38mean_tbl.csv',header = TRUE)

TPall <- TPall %>% mutate_all(~ifelse(is.nan(.), NA, .))
#TPall = subset(TPall, !is.na(TPall$BG_AGE))
TP = TPall[-c(1,23,29,34,38),]
TP_removed = TPall[c(1,23,29,34,38),]


#TPall2 = subset(TPall, TPall$Neurologic_Flag!=1&TPall$Cognitive_Flag!=1&TPall$Psychiatric_Flag!=1&TPall$ModSevTBI_Flag!=1&TPall$MSVT_Fail!=1 )
#TPall2$SA = as.numeric(TPall2$BSS_Q20!=0)
#TPall2$SI = as.numeric(TPall2$BSS_Total>0)
#TPall2$SI2 = as.numeric((TPall2$SA==0)&(TPall2$SI==1))
indx = cbind(grep("^Subj$", colnames(TPall)),grep("^Dscore$", colnames(TPall)),grep("^DeathWords_Me_mean$", colnames(TPall)),
             grep("^DeathWords_Other_mean$", colnames(TPall)),grep("^LifeWords_Me_mean$", colnames(TPall)),grep("^LifeWords_Other_mean$", colnames(TPall)),
             grep("^MeWords_Death_mean$", colnames(TPall)),grep("^MeWords_Life_mean$", colnames(TPall)), grep("^OtherWords_Death_mean$", colnames(TPall)),
             grep("^OtherWords_Life_mean$", colnames(TPall)), grep("^DeathWords_Me_SD$", colnames(TPall)),grep("^DeathWords_Other_SD$", colnames(TPall)),grep("^LifeWords_Me_SD$", colnames(TPall)),grep("^LifeWords_Other_SD$", colnames(TPall)),
             grep("^MeWords_Death_SD$", colnames(TPall)),grep("^MeWords_Life_SD$", colnames(TPall)), grep("^OtherWords_Death_SD$", colnames(TPall)),
             grep("^OtherWords_Life_SD$", colnames(TPall)),grep("^CAPS5_C$", colnames(TPall)),grep("^BSS_Q1$", colnames(TPall)),grep("^BSS_Q2$", colnames(TPall)),grep("^BSS_Q3$", colnames(TPall)),grep("^BSS_Q4$", colnames(TPall)),
             grep("^BSS_Q5$", colnames(TPall)),grep("^BSS_Q20$", colnames(TPall)),grep("^BSS_Total$", colnames(TPall)),grep("^DASS_Anxiety_TOT$", colnames(TPall)),
             grep("^PTSD_C_CAPS5$", colnames(TPall)),grep("^Current_MoodDx$", colnames(TPall)),grep("^Current_AnxietyDx$", colnames(TPall)),grep("^Current_SubstanceUseDx$", colnames(TPall)),
             grep("^DASS_Depression_TOT$", colnames(TPall)),grep("^PAIN_AVG$", colnames(TPall)),grep("^PSQI_GLOBAL$", colnames(TPall)),grep("^BG_AGE$", colnames(TPall)),
             grep("^BG_GEND$", colnames(TPall)),grep("^BG_ETHN_Hispanic$", colnames(TPall)),grep("^BG_RACE_White$", colnames(TPall)),grep("^BG_RACE_Black$", colnames(TPall)),
             grep("^BG_RACE_AmIndian$", colnames(TPall)),grep("^BG_RACE_Asian$", colnames(TPall)),grep("^BG_RACE_PacificIsland$", colnames(TPall)),grep("^BG_RACE_Other$", colnames(TPall)))
trial1 = TP[,c(indx)]
trial1$BSS_group = trial1$BSS_Total>2|trial1$BSS_Q20>0 #trial1$BSS_Q20>0 | trial1$BSS_Total>0;
trial1$DeathWordsIncC = trial1$DeathWords_Me_mean-trial1$DeathWords_Other_mean
trial1$LifeWordsIncC = trial1$LifeWords_Other_mean-trial1$LifeWords_Me_mean
trial1$MeWordsIncC = trial1$MeWords_Death_mean-trial1$MeWords_Life_mean
trial1$OtherWordsIncC = trial1$OtherWords_Life_mean-trial1$OtherWords_Death_mean

trial1$DeathWordsIncCV = (TP$DeathWords_Me_sd/TP$DeathWords_Me_mean)
trial1$DeathWordsCogCV = (TP$DeathWords_Other_sd/TP$DeathWords_Other_mean)
trial1$DeathWordsCV = (trial1$DeathWordsIncCV+trial1$DeathWordsCogCV)/2
trial1$LifeWordsIncCV = (TP$LifeWords_Other_sd/TP$LifeWords_Other_mean)
trial1$LifeWordsCogCV = (TP$LifeWords_Me_sd/TP$LifeWords_Me_mean)
trial1$LifeWordsCV = (trial1$LifeWordsIncCV+trial1$LifeWordsCogCV)/2
  
trial1$MeWordsIncCV = (TP$MeWords_Death_sd/TP$MeWords_Death_mean)
trial1$MeWordsCogCV = (TP$MeWords_Life_sd/TP$MeWords_Life_mean)
trial1$MeWordsCV = (trial1$MeWordsIncCV+trial1$MeWordsCogCV)/2
trial1$OtherWordsIncCV = (TP$OtherWords_Life_sd/TP$OtherWords_Life_mean)
trial1$OtherWordsCogCV = (TP$OtherWords_Death_sd/TP$OtherWords_Death_mean)
trial1$OtherWordsCV = (trial1$OtherWordsIncCV+trial1$OtherWordsCogCV)/2
trial1_noSTB = subset(trial1, BSS_group == 0)
trial1_STB = subset(trial1, BSS_group == 1)

#BrainData1_wm_noSTB = subset(BrainData1_wm, trial1$BSS_group == 0)
#BrainData1_wm_STB = subset(BrainData1_wm, trial1$BSS_group == 1)

#BrainData_clust_noSTB = subset(BrainData_clust, trial1$BSS_group == 0)
#BrainData_clust_STB = subset(BrainData_clust, trial1$BSS_group == 1)

### Pre and Post demographics 
tb = trial1 %>%
  tbl_summary(
    by = BSS_group,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{p}%"),
    digits = list(all_continuous() ~ 2,
                  all_categorical() ~ 2),
    type = list(c(BSS_Q1,BSS_Q2,BSS_Q3,BSS_Q4,BSS_Q5,BSS_Q20,BSS_Total) ~ "continuous"),
    #missing = 'yes'
    
  )%>%
  add_overall()%>% 
  add_p()

tb

hist(trial1$Dscore); range(trial1$Dscore); mean(trial1$Dscore)

####Figure for looking at interaction between congruency and word effects based on 
####reaction time.
indx = cbind(grep("^DeathWordsIncC$", colnames(trial1)),grep("^LifeWordsIncC$", colnames(trial1)),grep("^MeWordsIncC$", colnames(trial1)),grep("^OtherWordsIncC$", colnames(trial1)))
bardata = trial1[,c(indx)]
bardata2 = cbind(trial1$Subj,bardata); names(bardata2)[1]="Subj"
melt_bardata = melt(bardata)
melt_bardata2 = melt(bardata2,id="Subj")
STerr = c((sd(trial1$DeathWordsIncC/sqrt(length(trial1$DeathWordsIncC)))),(sd(trial1$LifeWordsIncC/sqrt(length(trial1$LifeWordsIncC)))),(sd(trial1$MeWordsIncC/sqrt(length(trial1$MeWordsIncC)))),(sd(trial1$OtherWordsIncC/sqrt(length(trial1$OtherWordsIncC)))))
means = c(mean(trial1$DeathWordsIncC),mean(trial1$LifeWordsIncC),mean(trial1$MeWordsIncC),mean(trial1$OtherWordsIncC))
ymin = means-STerr
ymax = means+STerr
labs =  c("DeathWords","LifeWords","MeWords","OtherWords")
bardatas = as.data.frame(cbind(labs,means,STerr))
bardatas$means = as.numeric(bardatas$means)
bardatas$STerr = as.numeric(bardatas$STerr)

p = ggplot(bardatas, aes(x=labs, y = means, fill=labs)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=means-STerr,ymax=means+STerr), width=.2,
                position=position_dodge(.9))+
  xlab("Word Type")+ylab("Incongruent-Congruent RTs (seconds)")+geom_hline(yintercept = 0)+ 
  theme(panel.background = element_blank(), panel.grid.major.y = element_line(color = "gray",size =.5),legend.position = "none")+
  scale_x_discrete(labels=c("DeathWords"="Death","LifeWords"="Life","MeWords"="Me","OtherWords"="Other"))
  
p
################################################################################

#this section makes a figure showing the brain activation related to an interaciton
##between congruency and word (wordBYcongruency effect)
################################################################################
Incon = rep("Incongruent", each=37)
Cong = rep("Congruent", each=37)
Fig_I = as.data.frame(cbind(Incon,BrainData2_wm$DW_I,BrainData2_wm$LW_I,
                             BrainData2_wm$MW_I,BrainData2_wm$NW_I))
colnames(Fig_I) = c("I_C","Death","Life","Me","NotMe")
Fig_I_m = melt(Fig_I,id = "I_C")

Fig_C = as.data.frame(cbind(Cong,BrainData2_wm$DW_C,BrainData2_wm$LW_C,
                            BrainData2_wm$MW_C,BrainData2_wm$NW_C))
colnames(Fig_C) = c("I_C","Death","Life","Me","NotMe")
Fig_C_m = melt(Fig_C,id = "I_C")


Fig_IC = rbind(Fig_I_m,Fig_C_m)
colnames(Fig_IC) = c("I_C","Word","Behavior")

Fig_IC$Behavior = as.numeric(Fig_IC$Behavior)


FigIC_plot = ggplot(Fig_IC, aes(x=Word, y=Behavior, fill=I_C, colour = I_C)) +
  geom_bar(stat="summary", position=position_dodge(.93))+ 
  geom_errorbar(stat='summary',position=position_dodge(.93), width=.2, colour="black")+
  theme_classic()
  
FigIC_plot

Fig_IC_2 = as.data.frame(cbind((BrainData2_wm$DW_I-BrainData2_wm$DW_C),
                               (BrainData2_wm$LW_I-BrainData2_wm$LW_C),
                               (BrainData2_wm$MW_I-BrainData2_wm$MW_C),
                               (BrainData2_wm$NW_I-BrainData2_wm$NW_C)))
colnames(Fig_IC_2) = c("DeathWords","LifeWords","MeWords","NotMeWords")
Fig_IC_2_m = melt(Fig_IC_2)
colnames(Fig_IC_2_m) = c("WordCategory","IncongruentCongruent")
MeanSDs = cbind(rbind(mean(Fig_IC_2$DeathWords),mean(Fig_IC_2$LifeWords),
           mean(Fig_IC_2$MeWords),mean(Fig_IC_2$NotMeWords)),
          rbind(std.error(Fig_IC_2$DeathWords),std.error(Fig_IC_2$LifeWords),
                std.error(Fig_IC_2$MeWords),std.error(Fig_IC_2$NotMeWords)));

FigIC2_plot = ggplot(Fig_IC_2_m, aes(x=WordCategory, y=IncongruentCongruent, fill=WordCategory, colour = WordCategory)) +
  geom_bar(stat="summary", position=position_dodge(.93))+ 
  geom_errorbar(stat='summary',position=position_dodge(.93), width=.2, colour="black")+
  theme_classic()
FigIC2_plot
################################################################################
wordlength = rbind(cbind(3,9,6,4),cbind(7,7,4,5),cbind(8,4,2,4),cbind(7,6,1,5),cbind(8,5,4,4))
meanWL = cbind(mean(wordlength[,1]),mean(wordlength[,2]),mean(wordlength[,3]),mean(wordlength[,4]))
##repeated measures anova for the behavioral interaction between congruency
##and word
indx = cbind(grep("^DeathWords_Me_mean$", colnames(trial1)),grep("^LifeWords_Other_mean$", colnames(trial1)),
             grep("^MeWords_Death_mean$", colnames(trial1)),grep("^OtherWords_Life_mean$", colnames(trial1)),
             grep("^DeathWords_Other_mean$", colnames(trial1)),grep("^LifeWords_Me_mean$", colnames(trial1)),
             grep("^MeWords_Life_mean$", colnames(trial1)),grep("^OtherWords_Death_mean$", colnames(trial1)))
temp1 =  trial1[,c(indx)]
temp2 = temp1[,1:4];temp2 = cbind(c(rep('Incongruent',37)),temp2)
names(temp2)[1] = "Congruency"; names(temp2)[2] = "Death";names(temp2)[3] = "Life";names(temp2)[4] = "Me";names(temp2)[5] = "NotMe"
temp2 = cbind(trial1$Subj,temp2);names(temp2)[1] = "Subj"
temp3 = temp1[,5:8];temp3 = cbind(c(rep('Congruent',37)),temp3)
names(temp3)[1] = "Congruency";names(temp3)[2] = "Death"; names(temp3)[3] = "Life"; names(temp3)[4] = "Me"; names(temp3)[5] = "NotMe"
temp3 = cbind(trial1$Subj,temp3); names(temp3)[1] = "Subj"
temp4 = rbind(temp2,temp3); temp5 = melt(temp4, id = c('Subj',"Congruency"))
test = c(rep(meanWL[1],37),rep(meanWL[1],37),rep(meanWL[2],37),rep(meanWL[2],37),
         rep(meanWL[3],37),rep(meanWL[3],37),rep(meanWL[4],37),rep(meanWL[4],37))
temp5_5 = cbind(temp5,test)

RManova_data = temp5; names(RManova_data)[3] = "Word"; names(RManova_data)[4] = "RT"
res.aov = anova_test(data=RManova_data,dv=RT,wid=Subj,within=c(Congruency,Word))
get_anova_table(res.aov)
pwc = melt_bardata2 %>% pairwise_t_test(value ~ variable, paired = TRUE); pwc

RManova_data2 = temp5_5; names(RManova_data2)[3] = "Word"; names(RManova_data2)[4] = "RT"
res.aov = aov(RT~Congruency+Word+test+Error(Subj),data=RManova_data2)
get_anova_table(res.aov)

################################################################################

##anova to look at interaction between congruency and word type in the brain, 
##this is lookign at a single brain region in teh right visual cortex
RManova_data = cbind(temp5$Subj,Fig_IC); names(RManova_data)[1] = "Subj"; names(RManova_data)[4] = "RT"
res.aov = anova_test(data=RManova_data,dv=RT,wid=Subj,within=c(I_C,Word))
get_anova_table(res.aov)
################################################################################

STerr = c((sd(trial1$DeathWords_Me_mean/sqrt(length(trial1$DeathWords_Me_mean)))),
          (sd(trial1$DeathWords_Other_mean/sqrt(length(trial1$DeathWords_Other_mean)))),
          (sd(trial1$LifeWords_Other_mean/sqrt(length(trial1$LifeWords_Other_mean)))),
          (sd(trial1$LifeWords_Me_mean/sqrt(length(trial1$LifeWords_Me_mean)))),
          (sd(trial1$MeWords_Death_mean/sqrt(length(trial1$MeWords_Death_mean)))),
          (sd(trial1$MeWords_Life_mean/sqrt(length(trial1$MeWords_Life_mean)))),
          (sd(trial1$OtherWords_Life_mean/sqrt(length(trial1$OtherWords_Life_mean)))),
          (sd(trial1$OtherWords_Death_mean/sqrt(length(trial1$OtherWords_Death_mean)))))
means = c(mean(trial1$DeathWords_Me_mean),mean(trial1$DeathWords_Other_mean),
          mean(trial1$LifeWords_Other_mean),mean(trial1$LifeWords_Me_mean),
          mean(trial1$MeWords_Death_mean),mean(trial1$MeWords_Life_mean),
          mean(trial1$OtherWords_Life_mean),mean(trial1$OtherWords_Death_mean))
Wordlabs =  c("DeathWords","DeathWords","LifeWords","LifeWords",
          "MeWords","MeWords","OtherWords","OtherWords")
IvC = c("I","C","I","C",
        "I","C","I","C")

bardatas2 = as.data.frame(cbind(Wordlabs,IvC,means,STerr))
bardatas2$means = as.numeric(bardatas2$means)
bardatas2$STerr = as.numeric(bardatas2$STerr)

p<- ggplot(bardatas2, aes(x=IvC, y=means, group=Wordlabs, color=Wordlabs)) + 
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=means-STerr,ymax=means+STerr), width=.2)+
  xlab("")+ylab("Reaction Times (seconds)")+
  theme(panel.background = element_blank(), panel.grid.major.y = element_line(color = "gray",size =.5))+
  scale_x_discrete(labels=c("I"="Incongruent","C"="Congruent"))+
  scale_color_discrete(name = "Word Type",labels=c("Death","Life","Me","Other"))+
  ggtitle("Reaction Time")
p

test = as.data.frame(trial1$Dscore)
p=ggplot(test,aes(x=trial1$Dscore))+geom_histogram()+xlab("D Score")+
  ylab("Count")+ggtitle("Histogram of D Scores")
p

mean(test$`trial1$Dscore`)
sd(test$`trial1$Dscore`)


###scatter plot looking at people with pos-dscores and their clinical data.

p1 =  ggplot(data = trial1, aes(x=Subj , y = Dscore))+ # creates the graph
  geom_point(color = "gray") + # creates main scatter plot with gray points
  geom_point(data = trial1_STB, color = "blue") # colors A.J. Hinch's point blue
p2 =  ggplot(data = trial1, aes(x=Subj , y = CAPS5_C))+ # creates the graph
  geom_point(color = "gray") + # creates main scatter plot with gray points
  geom_point(data = trial1_STB, color = "blue") # colors A.J. Hinch's point blue
p3 =  ggplot(data = trial1, aes(x=Subj , y = DASS_Depression_TOT))+ # creates the graph
  geom_point(color = "gray") + # creates main scatter plot with gray points
  geom_point(data = trial1_STB, color = "blue") # colors A.J. Hinch's point blue
p4 =  ggplot(data = trial1, aes(x=Subj , y = DASS_Anxiety_TOT))+ # creates the graph
  geom_point(color = "gray") + # creates main scatter plot with gray points
  geom_point(data = trial1_STB, color = "blue") # colors A.J. Hinch's point blue
p5 =  ggplot(data = trial1, aes(x=Subj , y = PAIN_AVG))+ # creates the graph
  geom_point(color = "gray") + # creates main scatter plot with gray points
  geom_point(data = trial1_STB, color = "blue") # colors A.J. Hinch's point blue
p6 =  ggplot(data = trial1, aes(x=Subj , y = PSQI_GLOBAL))+ # creates the graph
  geom_point(color = "gray") + # creates main scatter plot with gray points
  geom_point(data = trial1_STB, color = "blue") # colors A.J. Hinch's point blue

ps = ggarrange(p1,p2,p3,p4,p5,p6,ncol = 3,nrow = 2)

#correlate brain and behavior and clinical variables 
Dscore_BrainWM = cor.test(trial1$Dscore,BrainData1_wm$CongEff,method = "pearson")
Dscore_BrainWM

Dscore_BrainWM2 = cor.test(trial1$LifeWordsIncC,BrainData2_wm$CongEff,method = "pearson")
Dscore_BrainWM2

Dscore_BrainWM = cor.test(trial1$Dscore,BrainData_clust$CongEff_1,method = "pearson")
Dscore_BrainWM

Dscore_BrainWM = cor.test(trial1$Dscore,BrainData_clust$CongEff_1,method = "pearson")
Dscore_BrainWM

Dscore_BrainClust = rcorr(cbind(trial1$Dscore,
                    BrainData_clust$CongEff_1,BrainData_clust$CongEff_2,
                    BrainData_clust$CongEff_3,BrainData_clust$CongEff_4,BrainData_clust$CongEff_5,
                    BrainData_clust$CongEff_6,BrainData_clust$CongEff_7,BrainData_clust$CongEff_8,
                    BrainData_clust$CongEff_9,BrainData_clust$CongEff_10,BrainData_clust$CongEff_11,
                    BrainData_clust$CongEff_12,BrainData_clust$CongEff_13,BrainData_clust$CongEff_14,
                    BrainData_clust$CongEff_15,BrainData_clust$CongEff_16,BrainData_clust$CongEff_17,
                    BrainData_clust$CongEff_18,BrainData_clust$CongEff_19,BrainData_clust$CongEff_20,
                    BrainData_clust$CongEff_21,BrainData_clust$CongEff_22),type = "pearson")
Dscore_BrainClust
Dscore_BrainClust_r = Dscore_BrainClust$r[2:23,1];
Dscore_BrainClust_p = Dscore_BrainClust$P[2:23,1];
#Dscore_BrainClust_pIx = Dscore_BrainClust_p<(.05/22);
Dscore_BrainClust_pIx = Dscore_BrainClust_p<0.05;
Dscore_BrainClust_rp = Dscore_BrainClust_pIx*Dscore_BrainClust_r;
Dscore_BrainClust_rp = Dscore_BrainClust_rp


test = p.adjust(Dscore_BrainClust_p, method = "BH")
Dscore_BrainClust_pIx = test<0.05;
Dscore_BrainClust_rp = Dscore_BrainClust_pIx*Dscore_BrainClust_r;


Dep_Dscore = t.test(subset(trial1$Dscore,trial1$Current_MoodDx==1),subset(trial1$Dscore,trial1$Current_MoodDx==0), var.equal = FALSE)
Anx_Dscore= t.test(subset(trial1$Dscore,trial1$Current_AnxietyDx==1),subset(trial1$Dscore,trial1$Current_AnxietyDx==0), var.equal = FALSE); 
Sub_Dscore = t.test(subset(trial1$Dscore,trial1$Current_SubstanceUseDx==1),subset(trial1$Dscore,trial1$Current_SubstanceUseDx==0),var.equal = FALSE); 
PTSD_Dscore = t.test(subset(trial1$Dscore,trial1$PTSD_C_CAPS5==1),subset(trial1$Dscore,trial1$PTSD_C_CAPS5==0),var.equal = FALSE); 

Dep_BrainWM = t.test(subset(BrainData1_wm$CongEff,trial1$Current_MoodDx==1),subset(BrainData1_wm$CongEff,trial1$Current_MoodDx==0), var.equal = FALSE)
Anx_BrainWM= t.test(subset(BrainData1_wm$CongEff,trial1$Current_AnxietyDx==1),subset(BrainData1_wm$CongEff,trial1$Current_AnxietyDx==0), var.equal = FALSE); 
Sub_BrainWM = t.test(subset(BrainData1_wm$CongEff,trial1$Current_SubstanceUseDx==1),subset(BrainData1_wm$CongEff,trial1$Current_SubstanceUseDx==0),var.equal = FALSE); 
PTSD_BrainWM = t.test(subset(BrainData1_wm$CongEff,trial1$PTSD_C_CAPS5==1),subset(BrainData1_wm$CongEff,trial1$PTSD_C_CAPS5==0),var.equal = FALSE); 



Dscore1_BrainWM_lm = lm(formula = BrainData1_wm$CongEff ~ Dscore , data = trial1)
Dscore2_BrainWM_lm = lm(formula = BrainData1_wm$CongEff ~ Dscore+PTSD_C_CAPS5 , data = trial1)
Dscore3_BrainWM_lm = lm(formula = BrainData1_wm$CongEff ~ Dscore+PTSD_C_CAPS5+Current_MoodDx , data = trial1)
Dscore4_BrainWM_lm = lm(formula = BrainData1_wm$CongEff ~ Dscore+PTSD_C_CAPS5+Current_MoodDx+Current_AnxietyDx , data = trial1)
Dscore5_BrainWM_lm = lm(formula = BrainData1_wm$CongEff ~ Dscore+PTSD_C_CAPS5+Current_MoodDx+Current_AnxietyDx+Current_SubstanceUseDx, data = trial1)


#relationship between brain and clinical variables no STB only group
Dscore_BrainWM_noSTB = cor.test(trial1_noSTB$Dscore,BrainData1_wm_noSTB$V1,method = "pearson")
Dep_BrainWM_noSTB = cor.test(trial1_noSTB$DASS_Depression_TOT, BrainData1_wm_noSTB$V1,method = "pearson")
Anx_BrainWM_noSTB=cor.test(trial1_noSTB$DASS_Anxiety_TOT, BrainData1_wm_noSTB$V1,method = "pearson")
Pain_BrainWM_noSTB=cor.test(trial1_noSTB$PAIN_AVG, BrainData1_wm_noSTB$V1,method = "pearson")
Sleep_BrainWM_noSTB=cor.test(trial1_noSTB$PSQI_GLOBAL, BrainData1_wm_noSTB$V1,method = "pearson")
PTSD5_BrainWM_noSTB = cor.test(trial1_noSTB$CAPS5_C, BrainData1_wm_noSTB$V1,method = "pearson")

#relationship between brain and clinical variables in STB group only
Dscore_BrainWM_STB = cor.test(trial1_STB$Dscore,BrainData1_wm_STB$V1,method = "pearson")
Dep_BrainWM_STB = cor.test(trial1_STB$DASS_Depression_TOT, BrainData1_wm_STB$V1,method = "pearson")
Anx_BrainWM_STB=cor.test(trial1_STB$DASS_Anxiety_TOT, BrainData1_wm_STB$V1,method = "pearson")
Pain_BrainWM_STB=cor.test(trial1_STB$PAIN_AVG, BrainData1_wm_STB$V1,method = "pearson")
Sleep_BrainWM_STB=cor.test(trial1_STB$PSQI_GLOBAL, BrainData1_wm_STB$V1,method = "pearson")
PTSD5_BrainWM_STB = cor.test(trial1_STB$CAPS5_C, BrainData1_wm_STB$V1,method = "pearson")






#correlate brain and behavior and clinical variables - This is done for each cluster of main effect
test =rcorr(cbind(TP$PAIN_AVG,BrainData$CongEff_1,BrainData$CongEff_2,BrainData$CongEff_3,BrainData$CongEff_4,
                  BrainData$CongEff_5,BrainData$CongEff_6,BrainData$CongEff_7,BrainData$CongEff_8,
                  BrainData$CongEff_9,BrainData$CongEff_10,BrainData$CongEff_11,BrainData$CongEff_12,
                  BrainData$CongEff_13),type = "pearson")
test =rcorr(cbind(TP$CAPS_C,BrainData$CongEff_1,BrainData$CongEff_2,BrainData$CongEff_3,BrainData$CongEff_4,
                  BrainData$CongEff_5,BrainData$CongEff_6,BrainData$CongEff_7,BrainData$CongEff_8,
                  BrainData$CongEff_9,BrainData$CongEff_10,BrainData$CongEff_11,BrainData$CongEff_12,
                  BrainData$CongEff_13),type = "pearson")


#relationship between D-Score and clinical variables
Dep_Dscore = cor.test(TP$DASS_Depression_TOT, TP$Dscore,method = "pearson")
Anx_Dscore=cor.test(TP$DASS_Anxiety_TOT, TP$Dscore,method = "pearson")
Pain_Dscore=cor.test(TP$PAIN_AVG, TP$Dscore,method = "pearson")
Sleep_Dscore=cor.test(TP$PSQI_GLOBAL, TP$Dscore,method = "pearson")
PTSD4_Dscore=cor.test(TP$CAPS_C, TP$Dscore,method = "pearson")
PTSD5_Dscore = cor.test(TP$CAPS5_C, TP$Dscore,method = "pearson")

#relationship between D-Score and clinical variables no STB only group
Dep_Dscore_noSTB = cor.test(trial1_noSTB$DASS_Depression_TOT, trial1_noSTB$Dscore,method = "pearson")
Anx_Dscore_noSTB=cor.test(trial1_noSTB$DASS_Anxiety_TOT, trial1_noSTB$Dscore,method = "pearson")
Pain_Dscore_noSTB=cor.test(trial1_noSTB$PAIN_AVG, trial1_noSTB$Dscore,method = "pearson")
Sleep_Dscore_noSTB=cor.test(trial1_noSTB$PSQI_GLOBAL, trial1_noSTB$Dscore,method = "pearson")
PTSD5_Dscore_noSTB = cor.test(trial1_noSTB$CAPS5_C, trial1_noSTB$Dscore,method = "pearson")

#relationship between D-Score and clinical variables in STB group only
Dep_Dscore_STB = cor.test(trial1_STB$DASS_Depression_TOT, trial1_STB$Dscore,method = "pearson")
Anx_Dscore_STB=cor.test(trial1_STB$DASS_Anxiety_TOT, trial1_STB$Dscore,method = "pearson")
Pain_Dscore_STB=cor.test(trial1_STB$PAIN_AVG, trial1_STB$Dscore,method = "pearson")
Sleep_Dscore_STB=cor.test(trial1_STB$PSQI_GLOBAL, trial1_STB$Dscore,method = "pearson")
PTSD5_Dscore_STB = cor.test(trial1_STB$CAPS5_C, trial1_STB$Dscore,method = "pearson")

Dscore_corr = rcorr(cbind(trial1$Dscore,trial1$DeathWordsIncC,trial1$LifeWordsIncC,trial1$MeWordsIncC,trial1$OtherWordsIncC),type = "pearson")
Dscore_corr_nosTB = rcorr(cbind(trial1_noSTB$Dscore,trial1_noSTB$DeathWordsIncC,trial1_noSTB$LifeWordsIncC,trial1_noSTB$MeWordsIncC,trial1_noSTB$OtherWordsIncC),type = "pearson")
Dscore_corr_STB = rcorr(cbind(trial1_STB$Dscore,trial1_STB$DeathWordsIncC,trial1_STB$LifeWordsIncC,trial1_STB$MeWordsIncC,trial1_STB$OtherWordsIncC),type = "pearson")

Dscore_corr = rcorr(cbind(abs(trial1$Dscore),trial1$DeathWordsIncC,trial1$LifeWordsIncC,trial1$MeWordsIncC,trial1$OtherWordsIncC),type = "pearson")
Dscore_corr_nosTB = rcorr(cbind(abs(trial1_noSTB$Dscore),trial1_noSTB$DeathWordsIncC,trial1_noSTB$LifeWordsIncC,trial1_noSTB$MeWordsIncC,trial1_noSTB$OtherWordsIncC),type = "pearson")
Dscore_corr_STB = rcorr(cbind(abs(trial1_STB$Dscore),trial1_STB$DeathWordsIncC,trial1_STB$LifeWordsIncC,trial1_STB$MeWordsIncC,trial1_STB$OtherWordsIncC),type = "pearson")


lm_dscore = lm(abs(Dscore)~DeathWordsIncC+LifeWordsIncC+MeWordsIncC+OtherWordsIncC+PTSD_C_CAPS5, data = trial1)
lm_dscore = lm(abs(Dscore)~DeathWords_Me_mean+DeathWords_Other_mean+
                 LifeWords_Other_mean+LifeWords_Me_mean+
                 MeWords_Death_mean+MeWords_Life_mean+
                 OtherWords_Life_mean+OtherWords_Death_mean, data = trial1)


lm_dscore = lm(, data = trial1)

