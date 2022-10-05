library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(ggsignif)

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

#IO Config
##########

#Directories. This needs to be filled up to replicate the experiment
inDir="" #outDir of CNACalling.R
plotDir="/Users/diego/Documents/Ciencia/Postdoc/projects/dcis/lpwgs_study/lvi/plots/"

#Files
pattern="^(.*)[:.:]([^.]*)[:.:]best[:.:]ACNAs[:.:]c(.*)[:.:]p(.*)[:.:]tsv$" ##[::] is R's way of scape something, similar to \ or \\ in most regex

#Colors
stage.col = c('NORMAL' = "#FDBF6F",
              "DCIS" = "#A6CEE3",
              "IBC" = "#B2DF8A",
              "LVI" = "#FB9A99",
              "MET" ="#CAB2D6")

type.col=c("Proliferative" = "#D1DA60FF","Immune dense" = "#47A1DCFF")

cnaType.col=c("Gain"="#377eb8","Loss"="#e41a1c")

#Samples discarded due to problems
#Specifically, Met_1 and LVI_1 intensity/count plots look very noisy without clear alteration patterns
#Normal_2 does clearly not represent a normal sample, with plenty alterations very similar to the Met in the same patient
badSamples=data.table(patient=c("PATIENT_03","PATIENT_05","PATIENT_14"),sample=c("Met_1","LVI_1","NORMAL_2"))

#IO of previous results with absolute and relative copy number calls using QDNAseq + Rascal
fileList = dir(inDir, pattern, full.names = TRUE)
allCalls=rbindlist(lapply(fileList,function(x,pattern){theData=fread(x);theData[,`:=`(patient=gsub(basename(x),pattern=pattern,replacement="\\1"),
                                                                                      sample=gsub(basename(x),pattern=pattern,replacement="\\2"),
                                                                                      cellularity=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\3")),
                                                                                      ploidy=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\4")))]},pattern=pattern))

setkey(badSamples,patient,sample)
setkey(allCalls,patient,sample)

allCalls=allCalls[!badSamples,]

#Altered calls
##############
#library(rascal)
#relative_to_absolute_copy_number(1.06,8,0.2)
#relative_to_absolute_copy_number(0.94,8,0.2)
#These thresholds correspond to a gain or a loss of 1 copy at ploidy 8 and 0.2 cellularity, which are the max ploidy and min cell we are considering.
minNormalRelative=0.94
maxNormalRelative=1.06

alteredThresholdWithInfo=allCalls[relative_copy_number>=maxNormalRelative | relative_copy_number<=minNormalRelative,.(chrom=chromosome,start,end,relative_copy_number,sample_id=sample,Stage=gsub(sample,pattern="_.*$",replacement=""),patient=gsub(patient,pattern="PATIENT_",replacement = "P"))]
alteredThresholdWithInfo[Stage=="Met",`:=`(Stage="MET")]
alteredThresholdWithInfo[,`:=`(StageF=factor(Stage,levels=names(stage.col)))]

#Proportion of altered genome
alteredThresholdData=alteredThresholdWithInfo[,.(altered=sum(end-start),number=.N/2,stage=first(StageF)),by=.(patient,sample_id)]
setkey(alteredThresholdData,patient,sample_id)
totalGenome=allCalls[,.(chrom=chromosome,start,end,relative_copy_number,sample_id=sample,Stage=gsub(sample,pattern="_.*$",replacement=""),patient=gsub(patient,pattern="PATIENT_",replacement = "P"))][,.(genomeLength=sum(end-start)),by=.(patient,sample_id)]
setkey(totalGenome,patient,sample_id)
alteredThresholdData[totalGenome,`:=`(propAltered=altered/genomeLength)]

#Test between stages
#####################

#This testing approach is not correct since we are not taking into account the fact that most samples are not independent.
#kruskal.test(propAltered~stage,data=alteredThresholdData[stage!="NORMAL",])
#kruskal.test(number~stage,data=alteredThresholdData[stage!="NORMAL",])

#However, we can use the by-patient data (means) to mitigate the dependency issue
#################################################################################
alteredByPatientThresholdData=alteredThresholdData[,.(propAltered=mean(propAltered),number=mean(number)),by=.(patient,stage)]
kruskal.test(propAltered~stage,data=alteredByPatientThresholdData[stage!="NORMAL",])
kruskal.test(number~stage,data=alteredByPatientThresholdData[stage!="NORMAL",])

#Non-parametric version needs the design not to have missing data, so I will have to remove patients that do not have all stages
#The friedman.test function is a mess, so I need to make sure the data is a data.frame and all the factors are re-coded
dataFriedman=as.data.frame(alteredByPatientThresholdData[stage!="NORMAL",`:=`(patient=as.factor(patient))][patient %in% alteredByPatientThresholdData[stage!="NORMAL"][,.N,by=patient][N==4,]$patient,][stage!="NORMAL",.(patient,stage,propAltered)])
dataFriedman$patient=factor(dataFriedman$patient)
dataFriedman$stage=factor(dataFriedman$stage)
friedman.test(propAltered~stage|patient,data=dataFriedman)

dataFriedman=as.data.frame(alteredByPatientThresholdData[stage!="NORMAL",`:=`(patient=as.factor(patient))][patient %in% alteredByPatientThresholdData[stage!="NORMAL"][,.N,by=patient][N==4,]$patient,][stage!="NORMAL",.(patient,stage,number)])
dataFriedman$patient=factor(dataFriedman$patient)
dataFriedman$stage=factor(dataFriedman$stage)
friedman.test(number~stage|patient,data=dataFriedman)

#We can do non-parametric pairwise comparisons and this way we'd have more data for each comparison
myPairwiseWilcox(alteredByPatientThresholdData)
myPairwiseWilcox(alteredByPatientThresholdData,var = "number")


thresholdPropPlot=ggplot(data=alteredThresholdData,aes(x=stage,y=propAltered))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(aes(color=stage))+
  scale_color_manual(values=stage.col,guide="none") +
  scale_x_discrete(name="Stage") +
  scale_y_continuous(name="Proportion of the genome altered") +
  theme_cowplot() +
  theme(axis.title=element_text(face="bold"))

save_plot(filename=paste(sep="/",plotDir,paste0("thresholdPropPlot",".pdf")),plot=thresholdPropPlot,base_height = 6)

thresholdPropPlot

thresholdBurdenPlot=ggplot(data=alteredThresholdData,aes(x=stage,y=number))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(aes(color=stage))+
  scale_color_manual(values=stage.col,guide="none") +
  scale_x_discrete(name="Stage") +
  scale_y_continuous(name="Number of copy number alterations") +
  theme_cowplot() +
  theme(axis.title=element_text(face="bold"))

save_plot(filename=paste(sep="/",plotDir,paste0("thresholdBurdenPlot",".pdf")),plot=thresholdBurdenPlot,base_height = 6)

thresholdBurdenPlot

#Stats for results section of the manuscript
alteredThresholdData[,.(meanAltered=mean(propAltered),sdAltered=sd(propAltered),meanN=mean(number),sdN=sd(number)),by=.(stage)]