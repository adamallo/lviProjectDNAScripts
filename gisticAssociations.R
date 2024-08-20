library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)

#Config
setDTthreads(1) #I am doing coarse-grained parallelism

#IO Config
configFile <- paste(sep="/",Sys.getenv("lviProjectDNAScripts"),"configFile")
if(!file.exists(configFile)){
    stop("Configuration file configFile not found. Edit configFile.in to adapt it to your system, save it as configFile, and export the lviProjectDNAScripts environment variable before running this script")
}
source(configFile)
varFile <- paste(sep="/",acnaDir,"variance.tsv")
gisticFile <- gisticOutputFile

#IO
sampleData <- fread(sampleDataFile)
varData <- fread(varFile)
sampleData <- merge(sampleData,varData,by=c("patient","label"))
setkey(sampleData,fullLabel)
gisticData <- fread(gisticFile)
gisticData[,(colnames(gisticData)[ncol(gisticData)]):=NULL] #Remove last empty column

#Alterations to discard
#WARNING: I do not have RNAseq data so this is not being done here, but Belen previously found 3 alterations associated with the sequencing depth. Here, we are just removing them manually
invalidAlterations <- data.table(Type=c("gain","loss","loss"),Descriptor=c("1p35.2","3p14.1","11q13.1"))
invalidAlterations[,`:=`(CombinedDescriptor=paste(sep="_",Type,Descriptor))]

#Sample selection for when we want to use only one sample per patient and stage
selectedSamples <- sampleData[order(observedVariance),lapply(.SD,first),by=.(patient,Stage)]
setkey(selectedSamples,fullLabel)

#Data prep
gisticDataLong <- melt(gisticData[`Amplitude Threshold` != "Actual Copy Change Given",],id.vars=c(1:9),variable.name = "fullLabel",value.name = "Altered")
gisticDataLongBinary <- gisticDataLong[,`:=`(Altered=ifelse(Altered!=0,1,0), Type=ifelse(grepl(pattern = "Amplification", x = `Unique Name`),"gain","loss"))]
gisticDataLongBinary[,`:=`(CombinedDescriptor=paste(sep="_",Type,Descriptor))]
#Filtering invalid loci
gisticDataLongBinary <- gisticDataLongBinary[!CombinedDescriptor%in%invalidAlterations$CombinedDescriptor,]
#Merging sample and gistic data
finalGisticData <- merge(gisticDataLongBinary,sampleData,by="fullLabel")
finalGisticDataSelected <- merge(selectedSamples[,.(fullLabel)],finalGisticData,by = "fullLabel")
finalGisticDataUnionPerStage <- finalGisticData[,.(Altered=as.numeric(any(Altered==1)),LVI_subtype=unique(LVI_subtype)),by=.(CombinedDescriptor,patient,Stage)]
finalGisticDataUnionPerPatient <- finalGisticData[,.(Altered=as.numeric(any(Altered==1)),LVI_subtype=unique(LVI_subtype)),by=.(CombinedDescriptor,patient)]

#Stages

#All comparisons
StageCombinations <- combn(unique(finalGisticDataUnionPerStage$Stage),2)

getResultsStage <- function(thisData,stageFrom,stageTo){
  #Calculating associations
  results <- rbindlist(mclapply(mc.cores = mc.cores, unique(thisData$CombinedDescriptor),function(thisDescriptor){
    testStage <- fisher.test(table(thisData[CombinedDescriptor == thisDescriptor & Stage %in% c(stageFrom,stageTo),.(Altered,Stage)]))
    propsStage <- thisData[CombinedDescriptor == thisDescriptor & Stage %in% c(stageFrom,stageTo),.N, by =.(Altered,Stage)][,`:=`(Total=sum(N)),by=Stage][Altered==1,.(propAltered=N/Total),by=Stage]
    resultsStage <- data.table(CombinedDescriptor = thisDescriptor,p.value = testStage$p.value,props = propsStage$propAltered,Stage = propsStage$Stage, stageFrom=stageFrom,stageTo=stageTo) 
    return(resultsStage)
  }))
  results[,`:=`(adjusted.p.value=p.adjust(p.value,method="fdr"))]
  return(results)}

resultsStageAllCombs <- rbindlist(lapply(1:ncol(StageCombinations),function(i,thisData){
  getResultsStage(thisData,StageCombinations[1,i],StageCombinations[2,i])
},thisData = finalGisticDataUnionPerStage))

resultsStageAllCombs[p.value<0.05,]

resultsStageAllCombs1Sample <- rbindlist(lapply(1:ncol(StageCombinations),function(i,thisData){
  getResultsStage(thisData,StageCombinations[1,i],StageCombinations[2,i])
},thisData = finalGisticDataSelected))

resultsStageAllCombs1Sample[p.value<0.05,]


#Subtypes
getResultsSubtype <- function(thisData){
  #Calculating associations
  results <- rbindlist(mclapply(mc.cores = mc.cores, unique(thisData$CombinedDescriptor),function(thisDescriptor){
    testSubtype <- fisher.test(table(thisData[CombinedDescriptor == thisDescriptor,.(Altered,LVI_subtype)]))
    propsSubtype <- thisData[CombinedDescriptor == thisDescriptor,.N, by =.(Altered,LVI_subtype)][,`:=`(Total=sum(N)),by=LVI_subtype][Altered==1,.(propAltered=N/Total),by=LVI_subtype]
    resultsSubtype <- data.table(CombinedDescriptor = thisDescriptor,p.value = testSubtype$p.value,props = propsSubtype$propAltered,LVI_subtype=propsSubtype$LVI_subtype) 
    return(resultsSubtype)
  }))
  results[,`:=`(adjusted.p.value=p.adjust(p.value,method="fdr"))]
  return(results)}

resultsAcrossStagesUnionPerPatient <- getResultsSubtype(finalGisticDataUnionPerPatient)

resultsByStageUnion <- rbindlist(lapply(unique(finalGisticDataUnionPerStage$Stage),function(thisStage,theseData){
  theseResults <- getResultsSubtype(theseData[Stage==thisStage,])
  theseResults[,`:=`(Stage=thisStage)]
  return(theseResults)},theseData=finalGisticDataUnionPerStage))
resultsByStage1Sample <- rbindlist(lapply(unique(finalGisticDataSelected$Stage),function(thisStage,theseData){
  theseResults <- getResultsSubtype(theseData[Stage==thisStage,])
  theseResults[,`:=`(Stage=thisStage)]
  return(theseResults)},theseData=finalGisticDataSelected))

#I am removing this strategy because we still have multiple samples per patient with different numbers of samples, which may be an issue
#resultsAcrossStagesUnionPerStage <- getResultsSubtype(finalGisticDataUnionPerStage)
#resultsAcrossStages1Sample <- getResultsSubtype(finalGisticDataSelected)
#resultsAcrossStagesUnion[adjusted.p.value<0.05,]
#resultsAcrossStages1Sample[adjusted.p.value<0.05,]

resultsAcrossStagesUnionPerPatient[(!is.na(LVI_subtype)) & adjusted.p.value<0.05,]
resultsByStageUnion[(!is.na(LVI_subtype)) & adjusted.p.value<0.05,]
#resultsByStage1Sample[(!is.na(LVI_subtype)) & adjusted.p.value<0.05,] #I think the union is better than just getting one sample for this

resultsAcrossStagesUnionPerPatient[(!is.na(LVI_subtype)) & p.value<0.05,]
resultsByStageUnion[(!is.na(LVI_subtype)) & p.value<0.05,]
#resultsByStage1Sample[(!is.na(LVI_subtype)) & p.value<0.05,]  #I think the union is better than just getting one sample for this
