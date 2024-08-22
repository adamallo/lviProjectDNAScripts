#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggbeeswarm)
  library(cowplot)
  library(ggsignif)
  library(ggpubr)
  library(car)
  library(emmeans)
  library(glmmTMB)
})

#FUNCTIONS
myPairwiseWilcox=function(indt,var="propAltered"){
  outdt=data.table(patient=as.character(NULL),stage=as.character(NULL))
  for (patient in indt[,unique(patient)]){
    for(stage in indt[stage!="NORMAL",unique(stage)]){
      outdt=rbind(outdt,data.table(patient=patient,stage=stage))  
    }
  }
  setkey(indt,patient,stage)
  setkey(outdt,patient,stage)
  outdt=indt[outdt]
  outtable=data.table(t(combn(indt[stage!="NORMAL",unique(stage)],2,simplify=T)))
  outtable[,`:=`(pvalue=sapply(combn(indt[stage!="NORMAL",unique(stage)],2,simplify=F),FUN=function(x){wilcox.test(as.formula(paste0(var,"~stage")),data=outdt[stage==x[1] | stage==x[2],],paired=T,na.action="na.pass")$p.value}))]
  setnames(outtable,old=c("V1","V2"),new=c("Stage1","Stage2"))
  outtable[,`:=`(p.adj=p.adjust(pvalue,method="fdr"))]
  return(data.frame(outtable))
}

getPvalFromCoeffTable <- function(x,var="LVI_subtype",p.val.regex="Pr(.*)"){
  x[which(grepl(var,rownames(x))),which(grepl(p.val.regex,colnames(x)))]
}

#IO Config
##########
configFile <- paste(sep="/",Sys.getenv("lviProjectDNAScripts"),"configFile")
if(!file.exists(configFile)){
    stop("Configuration file configFile not found. Edit configFile.in to adapt it to your system, save it as configFile, and export the lviProjectDNAScripts environment variable before running this script")
}
source(configFile)
inDir <- acnaDir
dataDir <- acnaDir
summaryFile <- paste(sep="/",dataDir,summaryFileName)

#Files
pattern <- "^(.*)[:.:]([^.]*)[:.:]best[:.:]ACNAs[:.:]c(.*)[:.:]p(.*)[:.:]tsv$" ##[::] is R's way of scape something, similar to \ or \\ in most regex
date <- Sys.Date()

#IO config
##########
#Colors
stage.col <- c('NORMAL' = "#FDBF6F",
              "DCIS" = "#A6CEE3",
              "IBC" = "#B2DF8A",
              "LVI" = "#FB9A99",
              "MET" ="#CAB2D6")

type.col <- c("Proliferative" = "#D1DA60FF","Immune dense" = "#47A1DCFF")

cnaType.col <- c("Gain"="#377eb8","Loss"="#e41a1c")

baseHeight <- 6

#Labels
cohort.labs <- c("Disc" = "Discovery", "Val" = "Validation")

ggSaveWidth <- 5.2
ggSaveHeight <- 3.3
ggsignifStepIncrease <- 0.075
ommnibusTextSize <- 5
ggSaveOmmnibusTextSize <- 3.5
ggSaveGGsignifTextSize <- 3.5
ggSaveGGsignifSize <- 0.3
ggSaveLegendScaleFactor <- 1
jitterWidth <- 0.2
ylimScaling = 1.2
theTheme <- theme_bw()

#Data modification config
#########################
adjustCell <- 1E-6

#Samples discarded due to problems
#These normals clearly do not represent a normal sample, with plenty alterations very similar to the Met in the same patient
badSamples <- data.table(patient=c("344","740"),sample=c("NORMAL_2","NORMAL_1"))

if(file.exists(summaryFile)) {
  alteredData <- readRDS(summaryFile)
} else {
  #IO of previous results with absolute and relative copy number calls using QDNAseq + Rascal
  fileList <- dir(inDir, pattern, full.names = TRUE)
  allCalls <- rbindlist(lapply(fileList,function(x,pattern){theData=fread(x);theData[,`:=`(patient=gsub(basename(x),pattern=pattern,replacement="\\1"),
                                                                                           sample=gsub(basename(x),pattern=pattern,replacement="\\2"),
                                                                                           cellularity=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\3")),
                                                                                           ploidy=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\4")))]},pattern=pattern))
  
  setkey(badSamples,patient,sample)
  setkey(allCalls,patient,sample)
  
  allCalls <- allCalls[!badSamples,]
  
  #Altered calls
  ##############
  allCalls[,`:=`(log2_ratio_adjusted=log2(ifelse(absolute_copy_number<0,0.1,absolute_copy_number)/ploidy))]
  
  #We consider as alterations all segments of estimated absolute copy number that is different to the patient's ploidy, both rounded to integers
  alteredWithInfo <- allCalls[round(ploidy)!= round(absolute_copy_number),.(chrom=chromosome,
                                                                            start,
                                                                            end,
                                                                            relative_copy_number,
                                                                            sample_id=sample,
                                                                            Stage=gsub(sample,pattern="_.*$",replacement=""),
                                                                            patient=gsub(patient,pattern="PATIENT_",replacement = "P"),
                                                                            cellularity,
                                                                            ploidy)]

  #NOTE: The most important findings are robust to calling strategy. The exception is that the beta mixed model did not find significant differences between lvi_subtypes when using thresholds; the regular mixed model does, like either when using the current calling strategy.
  
  #Alternatively, using thresholds
  #library(rascal)
  #relative_to_absolute_copy_number(1.06,8,0.2)
  #relative_to_absolute_copy_number(0.94,8,0.2)
  #These thresholds correspond to a gain or a loss of 1 copy at ploidy 8 and 0.2 cellularity, which are the max ploidy and min cell we are considering.
  #minNormalRelative <- 0.94
  #maxNormalRelative <- 1.06
  #alteredWithInfo <- allCalls[relative_copy_number>=maxNormalRelative | relative_copy_number<=minNormalRelative,.(chrom=chromosome,start,end,relative_copy_number,sample_id=sample,Stage=gsub(sample,pattern="_.*$",replacement=""),patient=gsub(patient,pattern="PATIENT_",replacement = "P"))]
  
  alteredWithInfo[,`:=`(StageF=factor(Stage,levels=names(stage.col)))]
  
  #Proportion of altered genome
  alteredData <- alteredWithInfo[,.(altered=as.numeric(sum(end-start)),
                                    number=.N/2,
                                    stage=first(StageF),
                                    cellularity=first(cellularity),
                                    gainL=as.numeric(sum(ifelse(relative_copy_number>1,end-start,0))),
                                    lossL=as.numeric(sum(ifelse(relative_copy_number<1,end-start,0))),
                                    ploidy=first(ploidy)),by=.(patient,sample_id)]
  
  setkey(alteredData,patient,sample_id)
  totalGenome <- allCalls[,.(chrom=chromosome,start,end,relative_copy_number,sample_id=sample,Stage=gsub(sample,pattern="_.*$",replacement=""),patient=gsub(patient,pattern="PATIENT_",replacement = "P"))][,.(genomeLength=sum(end-start)),by=.(patient,sample_id)]
  setkey(totalGenome,patient,sample_id)
  alteredData[totalGenome,`:=`(propAltered=altered/genomeLength,propGain=gainL/lossL)]
  
  #Info integration
  sampleData <- fread(originalSampleDataFile)
  sampleData <- sampleData[,lapply(.SD,first),by=Sample_ID]#WARNING: there were repeated labels that I am removing here
  patientData <- sampleData[,.(cohort = unique(Cohort), LVI_subtype = unique(LVI_subtype)),by=.(patient = as.character(Case))]
  patientData[LVI_subtype == "NA",`:=`(LVI_subtype = NA)][,`:=`(LVI_subtype = factor(LVI_subtype))]
  setkey(patientData,patient)
  alteredData <- patientData[alteredData,]
  
  saveRDS(alteredData,file = summaryFile)
}

##################################
##Differences between Stages
##################################

stagesData <- alteredData[,.(patient,
                             cohort,
                             LVI_subtype,
                             stage,
                             propAltered,
                             propGain,
                             number,
                             cellularity=ifelse(cellularity==1,1-adjustCell,cellularity),
                             ploidy)][stage!="NORMAL",]

stagesWithNormalData <- alteredData[,.(patient,
                                  cohort,
                                  LVI_subtype,
                                  stage,
                                  propAltered = propAltered*100,
                                  propGain,
                                  number,
                                  cellularity=ifelse(cellularity==1,1-adjustCell,cellularity),
                                  ploidy)]

#Proportion of genome altered
#############################
#Beta mixed effects
stagesPropAlteredMEB <- glmmTMB(propAltered ~ stage + (1 | patient),
                                         family = beta_family(),
                                         data = stagesData)

Anova(stagesPropAlteredMEB, type = "III")
emmeans(stagesPropAlteredMEB, pairwise ~ stage, type = "response")

plot(residuals(stagesPropAlteredMEB, type = "pearson")) #They look good. Scattered around 0, not very far (>abs(3)) and evenly spread

#Mixed effects
#Less appropriate since propAltered are proportions.
# stagesPropAlteredMixedEffects <- aov(data=stagesData,propAltered~stage+Error(patient))
# summary(stagesPropAlteredMixedEffects)
# ggqqplot(stagesData,"propAltered") + facet_grid(~stage) #Prob good enough
# shapiro.test(stagesPropAlteredMixedEffects$patient$residuals) #Shapiro agrees, but I am unsure if just using the fixed residuals is enough
# shapiro.test(stagesPropAlteredMixedEffects$Within$residuals) #shapiro is definetely not ok for the random effects residuals
# leveneTest(propAltered ~ stage,stagesData) #levene ok

#WARNING: No significant differences so I am not adding any labels to the plot
propPlot <- ggplot(data=stagesWithNormalData,aes(x=stage,y=propAltered,color=stage))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(aes(color=stage))+
  scale_color_manual(values=stage.col,guide="none") +
  scale_x_discrete(name="") +
  scale_y_continuous(name='PGA (% genome length)') +
  theme_cowplot()

propPlotB <- ggplot(data=stagesWithNormalData,aes(x=stage,y=propAltered,fill=stage))+
  geom_boxplot(outlier.shape = NA)+
  #geom_quasirandom(aes(color=stage))+
  geom_jitter(aes(shape = cohort, group = stage), position=position_jitter(jitterWidth)) +
  scale_fill_manual(values=stage.col,guide="none") +
  scale_x_discrete(name="") +
  #scale_y_continuous(name="Proportion of the genome altered") +
  scale_y_continuous(name='PGA (% genome length)') +
  scale_shape_discrete(name="Cohort",labels=cohort.labs) +
  #theme_cowplot() +
  theme_bw()

save_plot(filename=paste(sep="/",plotDir,paste0(date,"_PGA_by_Stage_ALL_DF",".pdf")),plot = propPlot, base_height = baseHeight, create.dir = T)
ggsave(filename=paste(sep="/",plotDir,paste0(date,"_PGA_by_Stage_ALL.pdf")), plot = propPlotB, width = ggSaveWidth, height = ggSaveHeight)   


#Number of CNAs
###############
#Generalized mixed model with a poisson (or negative binomial) distribution (negative binomial if there is overdispersion, variance>mean)
stagesNAlterationsMEP <- glmmTMB(number ~ stage + (1 | patient),
                                family = poisson,
                                data = stagesData)
stagesNAlterationsMENB <- glmmTMB(number ~ stage + (1 | patient),
                                 family = nbinom2,
                                 data = stagesData)
overdispersion <- sum(resid(stagesNAlterationsMEP, type = "pearson")^2) / df.residual(stagesNAlterationsMEP)

if(overdispersion>1){stagesNAlterationsBestME <- stagesNAlterationsMENB}else{stagesNAlterationsBestME <- stagesNAlterationsMEP}

#If there was any doubt, the AIC of the negative binomial model is a much better fit
AIC(stagesNAlterationsMEP,stagesNAlterationsMENB)
plot(residuals(stagesNAlterationsBestME , type = "pearson")) #Scattered around 0, evenly spread, but some outliers >abs(3) 
Anova(stagesNAlterationsBestME, type = "III")
emmeans(stagesNAlterationsBestME, pairwise ~ stage, type = "response")
stagesNAlterationsBestMEeMeansResults <- as.data.frame(summary(emmeans(stagesNAlterationsBestME, pairwise ~ stage, type = "response"))$contrasts)

#Mixed effects
#Less appropriate since the number of alterations is a count variable
# stagesNIDMixedEffects  <- aov(data = stagesData,sqrt(number)~stage+Error(patient))
# summary(stagesNIDMixedEffects )
# emmeans(stagesNIDMixedEffects , pairwise ~ stage, type = "response")
# ggqqplot(stagesData,"number") + facet_grid(~stage) #needs sqrt
# shapiro.test(stagesNIDMixedEffects $patient$residuals) #shapiro is not ok, but the qqplots don't look very bad
# leveneTest(number ~ stage,stagesData)

#Significant differences, so I am adding labels to the plot
stagesNAlterationsBestMECoeffTable <- Anova(stagesNAlterationsBestME, type = "III")
stagesNAlterationsBestME.pval <- getPvalFromCoeffTable(stagesNAlterationsBestMECoeffTable,var = "stage")

#WARNING: ggsignif HARDCODED
maxY=max(stagesWithNormalData$number)*ylimScaling

numPlot <- ggplot(data=stagesWithNormalData,aes(x=stage,y=number,color=stage))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom()+
  geom_signif(color="black",
              comparisons = list(c("DCIS","IBC"),c("DCIS","MET")),
              annotations = c(sprintf("%0.1g",stagesNAlterationsBestMEeMeansResults$p.value[which(stagesNAlterationsBestMEeMeansResults$contrast == "DCIS / IBC")]),
                              sprintf("%0.1g",stagesNAlterationsBestMEeMeansResults$p.value[which(stagesNAlterationsBestMEeMeansResults$contrast == "DCIS / MET")])),,
              step_increase = ggsignifStepIncrease) +
  annotate("text",label=sprintf("Ommnibus Stage p = %.1g",stagesNAlterationsBestME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
  scale_color_manual(values=stage.col,guide="none") +
  scale_x_discrete(name="")+
  scale_y_log10(name="NGA",limits=c(NA,maxY)) +
  theme_cowplot()

numPlotB <- ggplot(data=stagesWithNormalData,aes(x=stage,y=number,fill=stage))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = cohort, group = stage), position=position_jitter(jitterWidth))+
  geom_signif(color="black",
              comparisons = list(c("DCIS","IBC"),c("DCIS","MET")),
              annotations = c(sprintf("%0.1g",stagesNAlterationsBestMEeMeansResults$p.value[which(stagesNAlterationsBestMEeMeansResults$contrast == "DCIS / IBC")]),
                              sprintf("%0.1g",stagesNAlterationsBestMEeMeansResults$p.value[which(stagesNAlterationsBestMEeMeansResults$contrast == "DCIS / MET")])),
              step_increase = ggsignifStepIncrease,
              textsize = ggSaveGGsignifTextSize,
              size = ggSaveGGsignifSize) +
  annotate("text",label=sprintf("Ommnibus Stage p = %.1g",stagesNAlterationsBestME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ggSaveOmmnibusTextSize) +
  scale_fill_manual(values=stage.col,guide="none") +
  scale_x_discrete(name="")+
  scale_y_continuous(name="NGA",limits=c(NA,maxY)) +
  scale_shape_discrete(name="Cohort",labels=cohort.labs) +
  theme_bw()

save_plot(filename=paste(sep="/",plotDir,paste0(date,"_NGA_by_Stage_All_DF",".pdf")),plot=numPlot,base_height = baseHeight)
ggsave(filename=paste(sep="/",plotDir,paste0(date,"_NGA_by_Stage_ALL.pdf")), plot = numPlotB, width = ggSaveWidth, height = ggSaveHeight) 

#Gain/Loss ratio
################
stagesGLRatioData <-  stagesData[,.SD][,`:=`(logPropGain=log(propGain))][!is.na(logPropGain) & is.finite(logPropGain),]
stagesGLRatioWithNormalData <-  stagesWithNormalData[,.SD][,`:=`(logPropGain=log(propGain))][!is.na(logPropGain) & is.finite(logPropGain),]


#Mixed effects
stagesGLRatioMixedEffects <- aov(data=stagesGLRatioData,logPropGain~stage+Error(patient))
summary(stagesGLRatioMixedEffects)


#Cellularity
############
#Beta mixed effects
stagesCellularityMEB <- glmmTMB(cellularity ~ stage + (1 | patient),
                                family = beta_family(),
                                data = stagesData)

Anova(stagesCellularityMEB, type = "III")
emmeans(stagesCellularityMEB, pairwise ~ stage, type = "response")
stagesCellularityMEBeMeansResults <- as.data.frame(summary(emmeans(stagesCellularityMEB, pairwise ~ stage, type = "response"))$contrasts)
plot(residuals(stagesCellularityMEB, type = "pearson")) #They look good

#Mixed effects
#Less appropriate since cellularity are proportions.
# stagesCellularityMixedEffects <- aov(data=stagesData,cellularity~stage+Error(patient))
# summary(stagesCellularityMixedEffects)
# ggqqplot(stagesData,"cellularity") + facet_grid(~stage) #Pretty bad
# shapiro.test(stagesCellularityMixedEffects$patient$residuals) #Shapiro okish, but I am unsure if just using the fixed residuals is enough
# shapiro.test(stagesCellularityMixedEffects$Within$residuals) #shapiro is definetely not ok for the random effects residuals
# leveneTest(cellularity ~ stage,stagesData) #levene ok

#Significant differences, so I am adding labels to the plot
stagesCellularityMEBCoeffTable <- Anova(stagesCellularityMEB, type = "III")
stagesCellularityMEB.pval <- getPvalFromCoeffTable(stagesCellularityMEBCoeffTable,var = "stage")

#WARNING: ggsignif HARDCODED
maxY=max(stagesWithNormalData$cellularity)*ylimScaling

cellularityPlot <- ggplot(data=stagesWithNormalData,aes(x=stage,y=cellularity,color=stage))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(aes(color=stage))+
  geom_signif(color="black",
              comparisons = list(c("LVI","MET"),c("IBC","LVI")),
              annotations = c(sprintf("%0.1g",stagesCellularityMEBeMeansResults$p.value[which(stagesCellularityMEBeMeansResults$contrast == "LVI / MET")]),
                              sprintf("%0.1g",stagesCellularityMEBeMeansResults$p.value[which(stagesCellularityMEBeMeansResults$contrast == "IBC / LVI")])),
              step_increase = ggsignifStepIncrease) +
  annotate("text",label=sprintf("Ommnibus Stage p = %.1g",stagesCellularityMEB.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
  scale_color_manual(values=stage.col,guide="none") +
  scale_x_discrete(name="") +
  scale_y_continuous(name='Cellularity (% cancer cells)',limits=c(NA,maxY)) +
  theme_cowplot()

cellularityPlotB <- ggplot(data=stagesWithNormalData,aes(x=stage,y=cellularity,fill=stage))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(shape = cohort, group = stage),position=position_jitter(jitterWidth)) +
  geom_signif(color="black",
              comparisons = list(c("LVI","MET"),c("IBC","LVI")),
              annotations = c(sprintf("%0.1g",stagesCellularityMEBeMeansResults$p.value[which(stagesCellularityMEBeMeansResults$contrast == "LVI / MET")]),
                              sprintf("%0.1g",stagesCellularityMEBeMeansResults$p.value[which(stagesCellularityMEBeMeansResults$contrast == "IBC / LVI")])),
              step_increase = ggsignifStepIncrease,
              textsize = ggSaveGGsignifTextSize,
              size = ggSaveGGsignifSize) +
  annotate("text",label=sprintf("Ommnibus Stage p = %.1g",stagesCellularityMEB.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ggSaveOmmnibusTextSize) +
  scale_fill_manual(values=stage.col,guide="none") +
  scale_x_discrete(name="") +
  #scale_y_continuous(name="Proportion of the genome altered") +
  scale_y_continuous(name='Cellularity (% cancer cells)',limits=c(NA,maxY)) +
  scale_shape_discrete(name="Cohort",labels=cohort.labs) +
  #theme_cowplot() +
  theme_bw() #+

save_plot(filename=paste(sep="/",plotDir,paste0(date,"_Cellularity_by_Stage_ALL_DF",".pdf")),plot = cellularityPlot, base_height = baseHeight)
ggsave(filename=paste(sep="/",plotDir,paste0(date,"_Cellularity_by_Stage_ALL.pdf")), plot = cellularityPlotB, width = ggSaveWidth, height = ggSaveHeight)   

#Results from this Section
##########################
#No differences between stages in proportion of genome altered or ratio between gains and losses. Differences in number of alterations between DCIS - IBC and DCIS - MET, but not robust to using patient means (see below)
#Differences in cellularity between LVI and IBC and LVI and MET. Biologically/technically expected for METs (and maybe IBC). LVI microdisection ensures good cellularity.

#We can use the by-patient data (means) to mitigate the dependency issue instead of using Mixed effects
#Doing this, we lose all signal
##############################################################################################################
# alteredByPatientThresholdData=alteredData[,.(propAltered=mean(propAltered),number=mean(number)),by=.(patient,stage)]
# kruskal.test(propAltered~stage,data=alteredByPatientThresholdData[stage!="NORMAL",])
# kruskal.test(sqrt(number)~stage,data=alteredByPatientThresholdData[stage!="NORMAL",])
# 
# stagesMeanNAlterationsMEP <- glmmTMB(number ~ stage,
#                                      family = poisson,
#                                      data = alteredByPatientThresholdData[stage!="NORMAL",])
# summary(stagesMeanNAlterationsMEP)
# 
# stagesMeanNAlterationsMENB <- glmmTMB(number ~ stage,
#                                       family = nbinom2,
#                                       data = alteredByPatientThresholdData[stage!="NORMAL",])
# summary(stagesMeanNAlterationsMENB)
# 
# overdispersion <- sum(resid(stagesMeanNAlterationsMEP, type = "pearson")^2) / df.residual(stagesMeanNAlterationsMEP)
# 
# if(overdispersion>1){stagesMeanNAlterationsBestME <- stagesMeanNAlterationsMENB}else{stagesMeanNAlterationsBestME <- stagesMeanNAlterationsMEP}
# 
# #If there was any doubt, the AIC of the negative binomial model is a much better fit
# AIC(stagesMeanNAlterationsMEP,stagesMeanNAlterationsMENB)
# 
# plot(residuals(stagesMeanNAlterationsBestME , type = "pearson")) #Scattered around 0, evenly spread, but some outliers >abs(3) 
# 
# Anova(stagesMeanNAlterationsBestME, type = "III")
# emmeans(stagesMeanNAlterationsBestME, pairwise ~ stage, type = "response")
# 
# #Non-parametric version needs the design not to have missing data, so I will have to remove patients that do not have all stages
# #The friedman.test function is a mess, so I need to make sure the data is a data.frame and all the factors are re-coded
# dataFriedman=as.data.frame(alteredByPatientThresholdData[stage!="NORMAL",`:=`(patient=as.factor(patient))][patient %in% alteredByPatientThresholdData[stage!="NORMAL"][,.N,by=patient][N==4,]$patient,][stage!="NORMAL",.(patient,stage,propAltered)])
# dataFriedman$patient=factor(dataFriedman$patient)
# dataFriedman$stage=factor(dataFriedman$stage)
# friedman.test(propAltered~stage|patient,data=dataFriedman)
# 
# dataFriedman=as.data.frame(alteredByPatientThresholdData[stage!="NORMAL",`:=`(patient=as.factor(patient))][patient %in% alteredByPatientThresholdData[stage!="NORMAL"][,.N,by=patient][N==4,]$patient,][stage!="NORMAL",.(patient,stage,number)])
# dataFriedman$patient=factor(dataFriedman$patient)
# dataFriedman$stage=factor(dataFriedman$stage)
# friedman.test(number~stage|patient,data=dataFriedman)
# 
# #We can do non-parametric pairwise comparisons and this way we'd have more data for each comparison
# myPairwiseWilcox(alteredByPatientThresholdData)
# myPairwiseWilcox(alteredByPatientThresholdData,var = "number")
# 
# 
# thresholdPropPlot=ggplot(data=alteredData,aes(x=stage,y=propAltered))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_quasirandom(aes(color=stage))+
#   scale_color_manual(values=stage.col,guide="none") +
#   scale_x_discrete(name="Stage") +
#   scale_y_continuous(name="Proportion of the genome altered") +
#   theme_cowplot() +
#   theme(axis.title=element_text(face="bold"))
# 
# save_plot(filename=paste(sep="/",plotDir,paste0("thresholdPropPlot_DF",".pdf")),plot=thresholdPropPlot,base_height = 6)
# 
# thresholdPropPlot
# 
# thresholdBurdenPlot=ggplot(data=alteredData,aes(x=stage,y=number))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_quasirandom(aes(color=stage))+
#   scale_color_manual(values=stage.col,guide="none") +
#   scale_x_discrete(name="Stage") +
#   scale_y_continuous(name="Number of copy number alterations") +
#   theme_cowplot() +
#   theme(axis.title=element_text(face="bold"))
# 
# save_plot(filename=paste(sep="/",plotDir,paste0("thresholdBurdenPlot_DF",".pdf")),plot=thresholdBurdenPlot,base_height = 6)
# 
# thresholdBurdenPlot

#Stats for results section of the manuscript
#Means using all samples
alteredData[,.(meanAltered=mean(propAltered),sdAltered=sd(propAltered),meanN=mean(number),sdN=sd(number)),by=.(stage)]
#Means of means per patient
alteredData[,.(LVI_subtype=unique(LVI_subtype),propAltered=mean(propAltered),propGain=mean(propGain),number=mean(number),propGain=mean(propGain)),by=.(patient,stage)][,.(meanAltered=mean(propAltered),sdAltered=sd(propAltered),meanN=mean(number),sdN=sd(number)),by=.(stage)]


alteredData[!is.na(LVI_subtype),.(meanAltered=mean(propAltered),sdAltered=sd(propAltered),meanN=mean(number),sdN=sd(number)),by=.(stage,LVI_subtype)]

##################################
##Differences between LVI subtypes
##################################

subtypeData <- alteredData[,.(patient,
                              cohort,
                              LVI_subtype,
                              stage,
                              propAltered,
                              propGain,
                              number,
                              cellularity=ifelse(cellularity==1,1-adjustCell,cellularity),
                              ploidy)][!is.na(LVI_subtype),][stage!="NORMAL",]

subtypeWithNormalData <- alteredData[,.(patient,
                              cohort,
                              LVI_subtype,
                              stage,
                              propAltered = propAltered*100,
                              propGain,
                              number,
                              cellularity=ifelse(cellularity==1,1-adjustCell,cellularity),
                              ploidy)][!is.na(LVI_subtype),]

#Proportion of genome altered
#############################

#Beta mixed effects --> LVI subtype significantly different propAltered
propAlteredLVISubtypesBME <- glmmTMB(propAltered ~ LVI_subtype + stage + (1 | patient),
                                         family = beta_family(),
                                         data = subtypeData)

summary(propAlteredLVISubtypesBME)
Anova(propAlteredLVISubtypesBME, type = "III")
plot(residuals(propAlteredLVISubtypesBME, type = "pearson")) #They look good. Scattered around 0, not very far (>abs(3)) and evenly spread

#Mixed effects --> LVI subtype significantly different propAltered
#Less appropriate since propAltered are proportions.
# propAlteredLVISubtypesME <- aov(data=subtypeData ,propAltered~LVI_subtype+stage+Error(patient))
# ggqqplot(subtypeData,"propAltered") + facet_grid(LVI_subtype~stage) #Prob good enough
# shapiro.test(propAlteredLVISubtypesME$patient$residuals) #shapiro ok
# leveneTest(propAltered ~ stage*LVI_subtype,subtypeData) #Levene ok
# summary(propAlteredLVISubtypesME)
# propAlteredLVISubtypesMEeMeansResults <-as.data.frame(summary(emmeans(propAlteredLVISubtypesME,pairwise~stage))$contrasts)
# propAlteredLVISubtypesMEeMeansResults

#Plot
propAlteredLVISubtypesBMECoeffTable <- summary(propAlteredLVISubtypesBME)[["coefficients"]]$cond
propAlteredLVISubtypesBME.pval <- getPvalFromCoeffTable(propAlteredLVISubtypesBMECoeffTable)
#propAlteredLVISubtypesME.pval <- getPvalFromCoeffTable(summary(propAlteredLVISubtypesME)[[1]][[1]]) #From the regular mixed effects model

propPlotLVISubtype <- ggplot(data=subtypeWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                             aes(x=stage,y=propAltered,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_quasirandom(dodge.width = 0.85)+
  scale_color_manual(name = "LVI subtype", values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="PGA (% genome length)") +
  annotate("text",label=sprintf("LVI subtype p = %.1g",propAlteredLVISubtypesBME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
  theme_cowplot()

propPlotLVISubtypeB <- ggplot(data=subtypeWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                              aes(x=stage,y=propAltered,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_jitter(position=position_jitterdodge(jitter.width = jitterWidth,dodge.width = 0.85))+
  scale_color_manual(name = "LVI subtype",values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="PGA (% genome length)") +
  annotate("text",label=sprintf("LVI subtype p = %.1g",propAlteredLVISubtypesBME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ggSaveOmmnibusTextSize) +
  theme_bw() +
  theme(
    legend.text = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor),
    legend.title = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor * 1.25),
    legend.key.size = unit(as.numeric(theTheme$legend.key.size) * ggSaveLegendScaleFactor, 'lines')
    )

save_plot(filename=paste(sep="/",plotDir,paste0(date,"_PGA_by_Subtype_Stage_All_DF",".pdf")),plot = propPlotLVISubtype,base_height = baseHeight)
ggsave(filename=paste(sep="/",plotDir,paste0(date,"_PGA_by_Subtype_Stage_All.pdf")), plot = propPlotLVISubtypeB, width = ggSaveWidth, height = ggSaveHeight) 

#Number of CNAs --> LVI subtype significantly different number of CNAs
###############

#Generalized mixed model with a poisson (or negative binomial) distribution (negative binomial if there is overdispersion, variance>mean)
numberLVISubtypeMEP <- glmmTMB(number ~ LVI_subtype + stage + (1 | patient),
                                 family = poisson,
                                 data = subtypeData)
numberLVISubtypeMENB <- glmmTMB(number ~ LVI_subtype + stage + (1 | patient),
                                  family = nbinom2,
                                  data = subtypeData)
overdispersion <- sum(resid(numberLVISubtypeMEP, type = "pearson")^2) / df.residual(numberLVISubtypeMEP)

if(overdispersion>1){numberLVISubtypeBestME <- numberLVISubtypeMENB}else{numberLVISubtypeBestME <- numberLVISubtypeMEP}

#If there was any doubt, the AIC of the negative binomial model is a much better fit
AIC(numberLVISubtypeMEP,numberLVISubtypeMENB)
plot(residuals(numberLVISubtypeBestME , type = "pearson")) #Scattered around 0, evenly spread, but some outliers >abs(3) 
Anova(numberLVISubtypeBestME, type = "III")
#emmeans(numberLVISubtypeBestME, pairwise ~ stage, type = "response") #I think we should not include this in the plot to reduce complexity. We do have the differences between stages without considering LVI subtype.
#emmeans(numberLVISubtypeBestME, pairwise ~ stage + LVI_subtype, type = "response")

#Mixed effects
#Less appropriate since the number of alterations is a count variable
# numberLVISubtypesME <- aov(data=subtypeData, sqrt(number)~LVI_subtype+stage+Error(patient))
# ggqqplot(subtypeData,"sqrt(number)") + facet_grid(LVI_subtype~stage) #Better with the sqrt, but still some outliers
# shapiro.test(numberLVISubtypesME$patient$residuals) #shapiro ok
# leveneTest(sqrt(number) ~ stage*LVI_subtype,subtypeData) #Levene ok
# summary(numberLVISubtypesME)

#Plot
numberLVISubtypesBestMECoeffTable <- summary(numberLVISubtypeBestME)[["coefficients"]]$cond
numberLVISubtypesBestME.pval <- getPvalFromCoeffTable(numberLVISubtypesBestMECoeffTable)
#numberLVISubtypesME.pval <- getPvalFromCoeffTable(summary(numberLVISubtypesME)[[1]][[1]]) #From the regular mixed effects model

numberPlotLVISubtype <- ggplot(data=subtypeWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                             aes(x=stage,y=number,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_quasirandom(dodge.width = 0.85)+
  scale_color_manual(name = "LVI subtype", values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="NGA") +
  annotate("text",label=sprintf("LVI subtype p = %.1g",numberLVISubtypesBestME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
  theme_cowplot()

numberPlotLVISubtypeB <- ggplot(data=subtypeWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                              aes(x=stage,y=number,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_jitter(position=position_jitterdodge(jitter.width = jitterWidth,dodge.width = 0.85))+
  scale_color_manual(name = "LVI subtype",values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="NGA") +
  annotate("text",label=sprintf("LVI subtype p = %.1g",numberLVISubtypesBestME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ggSaveOmmnibusTextSize) +
  theme_bw() +
  theme(
    legend.text = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor),
    legend.title = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor * 1.25),
    legend.key.size = unit(as.numeric(theTheme$legend.key.size) * ggSaveLegendScaleFactor, 'lines')
  )

save_plot(filename=paste(sep="/",plotDir,paste0(date,"_NGA_by_Subtype_Stage_All_DF",".pdf")),plot = numberPlotLVISubtype,base_height = baseHeight)
ggsave(filename=paste(sep="/",plotDir,paste0(date,"_NGA_by_Subtype_Stage_All.pdf")), plot = numberPlotLVISubtypeB, width = ggSaveWidth, height = ggSaveHeight) 


#Gain/Loss ratio
################
subtypeGLRatioData <-  subtypeData[,.SD][,`:=`(logPropGain=log(propGain))][!is.na(logPropGain) & is.finite(logPropGain),]
subtypeGLRatioWithNormalData <-  subtypeWithNormalData[,.SD][,`:=`(logPropGain=log(propGain))][!is.na(logPropGain) & is.finite(logPropGain),]


#Mixed effects
logPropGainLVISubtypesME <- aov(data=subtypeGLRatioData ,logPropGain~LVI_subtype+stage+Error(patient))
ggqqplot(subtypeGLRatioData,"logPropGain") + facet_grid(LVI_subtype~stage) #Prob good enough
shapiro.test(logPropGainLVISubtypesME$patient$residuals) #shapiro ok
shapiro.test(logPropGainLVISubtypesME$Within$residuals) #shapiro not k
leveneTest(logPropGain ~ stage*LVI_subtype,subtypeGLRatioData) #Levene ok
summary(logPropGainLVISubtypesME)

#Plot
logPropGainLVISubtypesME.pval <- getPvalFromCoeffTable(summary(logPropGainLVISubtypesME)[[1]][[1]]) #From the regular mixed effects model

propGLPlotLVISubtype <- ggplot(data=subtypeGLRatioWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                             aes(x=stage,y=propGain,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_quasirandom(dodge.width = 0.85)+
  scale_color_manual(name = "LVI subtype", values=type.col) +
  scale_x_discrete(name="") +
  scale_y_log10(name="Gain/Loss ratio") +
  annotate("text",label=sprintf("LVI subtype p = %.1g",logPropGainLVISubtypesME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
  theme_cowplot()

ylimScaling = 1.1
maxY=max(subtypeGLRatioWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))]$logPropGain)*ylimScaling

propGLPlotLVISubtypeB <- ggplot(data=subtypeGLRatioWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                              aes(x=stage,y=logPropGain,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_jitter(position=position_jitterdodge(jitter.width = jitterWidth,dodge.width = 0.85))+
  scale_color_manual(name = "LVI subtype",values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Gain/Loss ratio",limits=c(NA,maxY)) +
  annotate("text",label=sprintf("LVI subtype p = %.1g",logPropGainLVISubtypesME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ggSaveOmmnibusTextSize) +
  theme_bw() +
  theme(
    legend.text = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor),
    legend.title = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor * 1.25),
    legend.key.size = unit(as.numeric(theTheme$legend.key.size) * ggSaveLegendScaleFactor, 'lines')
  )
propGLPlotLVISubtypeB 
save_plot(filename=paste(sep="/",plotDir,paste0(date,"_GLRatio_by_Subtype_Stage_All_DF",".pdf")),plot = propGLPlotLVISubtype,base_height = baseHeight)
ggsave(filename=paste(sep="/",plotDir,paste0(date,"_GLRatio_by_Subtype_Stage_All.pdf")), plot = propGLPlotLVISubtypeB, width = ggSaveWidth, height = ggSaveHeight) 


#Cellularity --> No differences between subtypes across stages, nor for LVI
############

#Beta mixed effects --> LVI subtype significantly different cellularity
cellularityLVISubtypesBME <- glmmTMB(cellularity ~ LVI_subtype + stage + (1 | patient),
                                     family = beta_family(),
                                     data = subtypeData)

summary(cellularityLVISubtypesBME)
Anova(cellularityLVISubtypesBME, type = "III")
plot(residuals(cellularityLVISubtypesBME, type = "pearson")) #They look good. Scattered around 0, not very far (>abs(3)) and evenly spread

#Or for LVI
cellularityLVISubtypesLVIOnlyBME <- glmmTMB(cellularity ~ LVI_subtype + (1 | patient),
                                     family = beta_family(),
                                     data = subtypeData[stage=="LVI",])
summary(cellularityLVISubtypesLVIOnlyBME)


####################################################
##Differences between LVI subtypes: by-patient means
####################################################
subtypeByPatientData <- alteredData[!is.na(LVI_subtype),][stage!="NORMAL",][,.(LVI_subtype=unique(LVI_subtype),
                                                                               propAltered=mean(propAltered),
                                                                               propGain=mean(propGain),
                                                                               number=mean(number)),by=.(patient,stage)]
subtypeByPatientWithNormalData <- alteredData[!is.na(LVI_subtype),][,.(LVI_subtype=unique(LVI_subtype),
                                                                       propAltered=mean(propAltered)*100,
                                                                       propGain=mean(propGain),
                                                                       number=mean(number)),by=.(patient,stage)]

#Proportion of genome altered --> LVI subtype significantly different propAltered
#############################

#Beta mixed effects
propAlteredLVISubtypesByPatientBME <- glmmTMB(propAltered ~ LVI_subtype + stage,
                                     family = beta_family(),
                                     data = subtypeByPatientData)

summary(propAlteredLVISubtypesByPatientBME)
Anova(propAlteredLVISubtypesByPatientBME, type = "III")
plot(residuals(propAlteredLVISubtypesByPatientBME, type = "pearson")) #They do not look great but not awful either. They are scattered around 0, not very far (>abs(3)), but not that evenly spread

#Mixed effects
#Less appropriate since propAltered are proportions.
# propAlteredLVISubtypesByPatientME <- aov(data=subtypeByPatientData ,propAltered~LVI_subtype+stage)
# ggqqplot(subtypeByPatientData,"propAltered") + facet_grid(LVI_subtype~stage) #Prob good enough
# shapiro.test(propAlteredLVISubtypesByPatientME$residuals) #shapiro ok for very little
# leveneTest(propAltered ~ stage*LVI_subtype,subtypeByPatientData) #Levene ok
# summary(propAlteredLVISubtypesByPatientME)

#Plot
propAlteredLVISubtypesByPatientBMECoeffTable <- summary(propAlteredLVISubtypesByPatientBME)[["coefficients"]]$cond
propAlteredLVISubtypesByPatientBME.pval <- getPvalFromCoeffTable(propAlteredLVISubtypesByPatientBMECoeffTable)
#propAlteredLVISubtypesByPatientME.pval <- getPvalFromCoeffTable(summary(propAlteredLVISubtypesByPatientME)[[1]]) #From the regular mixed effects model

propByPatientPlotLVISubtype <- ggplot(data=subtypeByPatientWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                             aes(x=stage,y=propAltered,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_quasirandom(dodge.width = 0.85)+
  scale_color_manual(name = "LVI subtype", values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="PGA (% genome length)") +
  annotate("text",label=sprintf("LVI subtype p = %.1g",propAlteredLVISubtypesByPatientBME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
  theme_cowplot()

propByPatientPlotLVISubtypeB <- ggplot(data=subtypeByPatientWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                              aes(x=stage,y=propAltered,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_jitter(position=position_jitterdodge(jitter.width = jitterWidth,dodge.width = 0.85))+
  scale_color_manual(name = "LVI subtype",values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="PGA (% genome length)") +
  annotate("text",label=sprintf("LVI subtype p = %.1g",propAlteredLVISubtypesByPatientBME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ggSaveOmmnibusTextSize) +
  theme_bw() +
  theme(
    legend.text = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor),
    legend.title = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor * 1.25),
    legend.key.size = unit(as.numeric(theTheme$legend.key.size) * ggSaveLegendScaleFactor, 'lines')
  )

save_plot(filename=paste(sep="/",plotDir,paste0(date,"_PGA_by_Subtype_Stage_ByPatient_DF",".pdf")),plot = propByPatientPlotLVISubtype,base_height = baseHeight)
ggsave(filename=paste(sep="/",plotDir,paste0(date,"_PGA_by_Subtype_Stage_ByPatient.pdf")), plot = propByPatientPlotLVISubtypeB, width = ggSaveWidth, height = ggSaveHeight) 

#Number of CNAs --> LVI subtype significantly different number of CNAs
###############

#Generalized mixed model with a poisson (or negative binomial) distribution (negative binomial if there is overdispersion, variance>mean)
numberLVISubtypeByPatientMEP <- glmmTMB(number ~ LVI_subtype + stage,
                               family = poisson,
                               data = subtypeByPatientData)
numberLVISubtypeByPatientMENB <- glmmTMB(number ~ LVI_subtype + stage,
                                family = nbinom2,
                                data = subtypeByPatientData)
overdispersion <- sum(resid(numberLVISubtypeByPatientMEP, type = "pearson")^2) / df.residual(numberLVISubtypeByPatientMEP)

if(overdispersion>1){numberLVISubtypeByPatientBestME <- numberLVISubtypeByPatientMENB}else{numberLVISubtypeByPatientBestME <- numberLVISubtypeByPatientMEP}

#If there was any doubt, the AIC of the negative binomial model is a much better fit
AIC(numberLVISubtypeByPatientMEP,numberLVISubtypeByPatientMENB)
plot(residuals(numberLVISubtypeByPatientBestME , type = "pearson")) #Scattered around 0, evenly spread, but some outliers >abs(3) 
Anova(numberLVISubtypeByPatientBestME, type = "III")

#Mixed effects
#Less appropriate since the number of alterations is a count variable
# numberLVISubtypesByPatientME <- aov(data=subtypeByPatientData, sqrt(number)~LVI_subtype+stage)
# ggqqplot(subtypeByPatientData,"sqrt(number)") + facet_grid(LVI_subtype~stage) #MET outliers, otherwise probably fine
# shapiro.test(numberLVISubtypesByPatientME$residuals) #shapiro very not ok
# leveneTest(sqrt(number) ~ stage*LVI_subtype,subtypeByPatientData) #Levene ok
# summary(numberLVISubtypesByPatientME)

#Plot
numberLVISubtypeByPatientBestMECoeffTable <- summary(numberLVISubtypeByPatientBestME)[["coefficients"]]$cond
numberLVISubtypeByPatientBestME.pval <- getPvalFromCoeffTable(numberLVISubtypeByPatientBestMECoeffTable)
#numberLVISubtypesME.pval <- getPvalFromCoeffTable(summary(numberLVISubtypesME)[[1]][[1]]) #From the regular mixed effects model

ylimScaling = 1.1
maxY=max(subtypeByPatientWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))]$number)*ylimScaling

numberByPatientPlotLVISubtype <- ggplot(data=subtypeByPatientWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                               aes(x=stage,y=number,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_quasirandom(dodge.width = 0.85)+
  scale_color_manual(name = "LVI subtype", values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="NGA") +
  annotate("text",label=sprintf("LVI subtype p = %.1g",numberLVISubtypeByPatientBestME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
  theme_cowplot()

numberByPatientPlotLVISubtypeB <- ggplot(data=subtypeByPatientWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                                aes(x=stage,y=number,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_jitter(position=position_jitterdodge(jitter.width = jitterWidth,dodge.width = 0.85))+
  scale_color_manual(name = "LVI subtype",values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="NGA",limits=c(NA,maxY)) +
  annotate("text",label=sprintf("LVI subtype p = %.1g",numberLVISubtypeByPatientBestME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ggSaveOmmnibusTextSize) +
  theme_bw() +
  theme(
    legend.text = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor),
    legend.title = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor * 1.25),
    legend.key.size = unit(as.numeric(theTheme$legend.key.size) * ggSaveLegendScaleFactor, 'lines')
  )

save_plot(filename=paste(sep="/",plotDir,paste0(date,"_NGA_by_Subtype_Stage_ByPatient_DF",".pdf")),plot = numberByPatientPlotLVISubtype,base_height = baseHeight)
ggsave(filename=paste(sep="/",plotDir,paste0(date,"_NGA_by_Subtype_Stage_ByPatient.pdf")), plot = numberByPatientPlotLVISubtypeB, width = ggSaveWidth, height = ggSaveHeight) 

#Gain/Loss ratio
#############################
subtypeGLRatioData <-  subtypeByPatientData [,.SD][,`:=`(logPropGain=log(propGain))][!is.na(logPropGain) & is.finite(logPropGain),]
subtypeGLRatioData <-  subtypeByPatientWithNormalData [,.SD][,`:=`(logPropGain=log(propGain))][!is.na(logPropGain) & is.finite(logPropGain),]

#Mixed effects
logPropGainLVISubtypesByPatientME <- aov(data=subtypeGLRatioData ,logPropGain~LVI_subtype+stage)
ggqqplot(subtypeGLRatioData,"logPropGain") + facet_grid(LVI_subtype~stage) #Prob good enough
shapiro.test(logPropGainLVISubtypesByPatientME$residuals) #shapiro NOT ok
ggqqplot(logPropGainLVISubtypesByPatientME$residuals) #but the residuals' qqplot is not that bad
leveneTest(logPropGain ~ stage*LVI_subtype,subtypeGLRatioData) #Levene ok
summary(logPropGainLVISubtypesByPatientME)

#Plot
logPropGainLVISubtypesByPatientME.pval <- getPvalFromCoeffTable(summary(logPropGainLVISubtypesByPatientME)[[1]])

ylimScaling = 1.1
maxY=max(subtypeGLRatioWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))]$logPropGain)*ylimScaling

propGLByPatientPlotLVISubtype <- ggplot(data=subtypeGLRatioWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                                      aes(x=stage,y=logPropGain,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_quasirandom(dodge.width = 0.85)+
  scale_color_manual(name = "LVI subtype", values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Gain/Loss ratio") +
  annotate("text",label=sprintf("LVI subtype p = %.1g",logPropGainLVISubtypesByPatientME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
  theme_cowplot()

propGLByPatientPlotLVISubtypeB <- ggplot(data=subtypeGLRatioWithNormalData[,.SD][,`:=`(LVI_subtype=factor(LVI_subtype,levels=c("ID","Prol"),labels=c("Immune dense","Proliferative")))],
                                       aes(x=stage,y=logPropGain,color=LVI_subtype))+
  geom_boxplot(position="dodge2")+
  geom_jitter(position=position_jitterdodge(jitter.width = jitterWidth,dodge.width = 0.85))+
  scale_color_manual(name = "LVI subtype",values=type.col) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Gain/Loss ratio",limits=c(NA,maxY)) +
  annotate("text",label=sprintf("LVI subtype p = %.1g",logPropGainLVISubtypesByPatientME.pval), y= Inf, x =Inf,vjust=1, hjust=1, size = ggSaveOmmnibusTextSize) +
  theme_bw() +
  theme(
    legend.text = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor),
    legend.title = element_text(size = theTheme$text$size * ggSaveLegendScaleFactor * 1.25),
    legend.key.size = unit(as.numeric(theTheme$legend.key.size) * ggSaveLegendScaleFactor, 'lines')
  )

save_plot(filename=paste(sep="/",plotDir,paste0(date,"_GLRatio_by_Subtype_Stage_ByPatient_DF",".pdf")),plot = propGLByPatientPlotLVISubtype,base_height = baseHeight)
ggsave(filename=paste(sep="/",plotDir,paste0(date,"_GLRatio_by_Subtype_Stage_ByPatient.pdf")), plot = propGLByPatientPlotLVISubtypeB, width = ggSaveWidth, height = ggSaveHeight) 

##############################################

#Results from this Section
#Proliferative patients have an increase in proportion of genome altered, number of alterations, and the gain/ratio ratio.This finding is robust to the mixed model used for the two first (gaussian vs. beta for proportions and negative binomial for counts) and using patient means.
#For ratios, there isn't a better alternative than logRatios with a gaussian that I know of 

#Stratifying the comparison across Stage by subtype:
####################################################

#In this section, we see that we find differences in number of CNAs between DCIS and MET in IDs, nothing in Prols. 
#I feel the results of this section don't add much/anything, and thus are commented out, but they are checked and updated
# 
# ###
# #ID
# ###
# propAlteredIDMEData <- alteredData[,.(patient,cohort,LVI_subtype,stage,propAltered,number)][LVI_subtype=="ID",][stage!="NORMAL",]
# 
# #Proportion of Genome Altered
# #Mixed effects --> No differences
# stagesPropAlteredIDMixedEffects <- aov(data=propAlteredIDMEData ,propAltered~stage+Error(patient))
# summary(stagesPropAlteredIDMixedEffects) #Results are equivalent with or without SQRT
# 
# #Beta mixed effects --> No differences
# #WARNING: this is more appropriate since propAltered are proportions.
# stagesPropAlteredIDMEB <- glmmTMB(propAltered ~ stage + (1 | patient),
#                                 family = beta_family(),
#                                 data = propAlteredIDMEData)
# Anova(stagesPropAlteredIDMEB, type = "III")
# 
# #Number of alterations --> Very close to differences
# #Mixed effects
# stagesNIDMixedEffects <- aov(data=propAlteredIDMEData ,sqrt(number)~stage+Error(patient))
# summary(stagesNIDMixedEffects) #Very close to being significant
# emmeans(stagesNIDMixedEffects , pairwise ~ stage, type = "response") #Between DCIS and MET
# #But the model fit is not great
# ggqqplot(stagesData,"sqrt(number)") + facet_grid(~stage) #Met pretty far from normal
# shapiro.test(stagesNIDMixedEffects$patient$residuals) #shapiro is ok
# shapiro.test(stagesNIDMixedEffects$Within$residuals) #shapiro is NOT OK
# leveneTest(number ~ stage,stagesData)
# 
# #Generalized mixed model with a poisson (or negative binomial) distribution (negative binomial if there is overdispersion, variance>mean)
# # --> Very significant differences between DCIS and MET
# #WARNING: this is more appropriate since these are counts
# stagesNIDAlterationsMEP <- glmmTMB(number ~ stage + (1 | patient),
#                                  family = poisson,
#                                  data = propAlteredIDMEData)
# stagesNIDAlterationsMENB <- glmmTMB(number ~ stage + (1 | patient),
#                                   family = nbinom2,
#                                   data = propAlteredIDMEData)
# overdispersion <- sum(resid(stagesNIDAlterationsMEP, type = "pearson")^2) / df.residual(stagesNIDAlterationsMEP)
# if(overdispersion>1){stagesNIDAlterationsBestME <- stagesNIDAlterationsMENB}else{stagesNIDAlterationsBestME <- stagesNIDAlterationsMEP}
# 
# #If there was any doubt, the AIC of the negative binomial model is a much better fit
# AIC(stagesNIDAlterationsMEP,stagesNIDAlterationsMENB)
# 
# plot(residuals(stagesNIDAlterationsBestME , type = "pearson")) #Scattered around 0, evenly spread, but some outliers >abs(3) 
# Anova(stagesNIDAlterationsBestME, type = "III")
# emmeans(stagesNIDAlterationsBestME, pairwise ~ stage, type = "response")
# ###########################################################
# 
# ###
# #Prol
# ###
# propAlteredProlMEData <- alteredData[,.(patient,cohort,LVI_subtype,stage,propAltered,number)][LVI_subtype=="Prol",][stage!="NORMAL",]
# 
# #Proportion of Genome Altered --> No differences
# #Mixed effects --> No differences
# stagesPropAlteredProlMixedEffects <- aov(data=propAlteredProlMEData ,propAltered~stage+Error(patient))
# summary(stagesPropAlteredProlMixedEffects) #Results are equivalent with or without SQRT
# 
# #Beta mixed effects --> No differences
# #WARNING: this is more appropriate since propAltered are proportions.
# stagesPropAlteredProlMEB <- glmmTMB(propAltered ~ stage + (1 | patient),
#                                   family = beta_family(),
#                                   data = propAlteredProlMEData)
# Anova(stagesPropAlteredProlMEB, type = "III")
# 
# #Number of alterations --> No differences
# #Mixed effects
# stagesNProlMixedEffects <- aov(data=propAlteredProlMEData ,sqrt(number)~stage+Error(patient))
# 
# #Generalized mixed model with a poisson (or negative binomial) distribution (negative binomial if there is overdispersion, variance>mean)
# #WARNING: this is more appropriate since these are counts
# stagesNProlAlterationsMEP <- glmmTMB(number ~ stage + (1 | patient),
#                                    family = poisson,
#                                    data = propAlteredProlMEData)
# stagesNProlAlterationsMENB <- glmmTMB(number ~ stage + (1 | patient),
#                                     family = nbinom2,
#                                     data = propAlteredProlMEData)
# overdispersion <- sum(resid(stagesNProlAlterationsMEP, type = "pearson")^2) / df.residual(stagesNProlAlterationsMEP)
# if(overdispersion>1){stagesNProlAlterationsBestME <- stagesNProlAlterationsMENB}else{stagesNProlAlterationsBestME <- stagesNProlAlterationsMEP}
# 
# #If there was any doubt, the AIC of the negative binomial model is a much better fit
# AIC(stagesNProlAlterationsMEP,stagesNProlAlterationsMENB)
# Anova(stagesNProlAlterationsBestME, type = "III")
# emmeans(stagesNProlAlterationsBestME, pairwise ~ stage, type = "response")

