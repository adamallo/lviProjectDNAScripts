library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

#IO Config
##########
configFile <- paste(sep="/",Sys.getenv("lviProjectDNAScripts"),"configFile")
if(!file.exists(configFile)){
    stop("Configuration file configFile not found. Edit configFile.in to adapt it to your system, save it as configFile, and export the lviProjectDNAScripts environment variable before running this script")
}
source(configFile)
inDir <- acnaDir
varFile <- paste(sep="/",acnaDir,"variance.tsv")

if(gisticLVIOnly==T){
    gisticFile <- gisticLVIOutputFile
} else {
    gisticFile <- gisticOutputFile
}

#Files
pattern <- "^(.*)[:.:]([^.]*)[:.:]best[:.:]ACNAs[:.:]c(.*)[:.:]p(.*)[:.:]tsv$" ##[::] is R's way of escaping something, similar to \ or \\ in most regex

#Colors
stage.col <- c('NORMAL' = "#FDBF6F",
               "DCIS" = "#A6CEE3",
               "IBC" = "#B2DF8A",
               "LVI" = "#FB9A99",
               "MET" ="#CAB2D6")

type.col <- c("Proliferative" = "#D1DA60FF","Immune dense" = "#47A1DCFF")

cnaType.col <- c("Gain"="#377eb8","Loss"="#e41a1c")

#Samples discarded due to problems
#These normals clearly do not represent a normal sample, with plenty alterations very similar to the Met in the same patient
badSamples <- data.table(patient=c("344","740"),sample=c("NORMAL_2","NORMAL_1"))

#IO of previous results with absolute and relative copy number calls using QDNAseq + Rascal
fileList <- dir(inDir, pattern, full.names = TRUE)
allCalls <- rbindlist(lapply(fileList,function(x,pattern){theData=fread(x);theData[,`:=`(patient=gsub(basename(x),pattern=pattern,replacement="\\1"),
                                                                                         sample=gsub(basename(x),pattern=pattern,replacement="\\2"),
                                                                                         cellularity=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\3")),
                                                                                         ploidy=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\4")))]},pattern=pattern))

setkey(badSamples,patient,sample)
setkey(allCalls,patient,sample)

allCalls <- allCalls[!badSamples,]

allCalls[,`:=`(log2_ratio_adjusted=log2(ifelse(absolute_copy_number<0,0.1,absolute_copy_number)/ploidy))]
selectedCalls <- allCalls[round(ploidy)!= round(absolute_copy_number) & (log2_ratio_adjusted>0.1 | log2_ratio_adjusted<0.1),] #Equivalent to gistic's threshold strategy, if I remove the threshold, equivalent to genetics.R 

grCnas <- GRanges(
  seqnames = paste0("chr",selectedCalls$chromosome),
  ranges = IRanges(start = selectedCalls$start, end = selectedCalls$end),
  patient = selectedCalls$patient,
  sample = selectedCalls$sample,
  log2_ratio_adjusted = selectedCalls$log2_ratio_adjusted
)


#Prepare reference data
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
refGenes <- genes(txdb,single.strand.genes.only=FALSE)

#Filter reference data using GISTIC results (i.e., only consider genes that are commonly altered)
######################################################

#Parsing GISTIC results
gisticData <- fread(gisticFile)
gisticData[,(colnames(gisticData)[ncol(gisticData)]):=NULL] #Remove last empty column
gisticData[,`:=`(Type=ifelse(grepl(pattern = "Amplification", x = `Unique Name`),"gain","loss"))]
gisticData[,`:=`(CombinedDescriptor=paste(sep="_",Type,Descriptor))]

#Filtering
#I do not have the sequencing data so I can't do it here, but Belen previously found 3 alterations associated with the sequencing depth, so I am just removing them here
invalidAlterations <- data.table(Type=c("gain","loss","loss"),Descriptor=c("1p35.2","3p14.1","11q13.1"))
invalidAlterations[,`:=`(CombinedDescriptor=paste(sep="_",Type,Descriptor))]

cleanGisticData <- gisticData[`Amplitude Threshold` != "Actual Copy Change Given",][!CombinedDescriptor%in%invalidAlterations$CombinedDescriptor,]

grGistic <- GRanges(
  seqnames = gsub("([^:]*):.*","\\1",cleanGisticData$`Wide Peak Limits`),
  ranges = IRanges(start = as.numeric(gsub(".*:([0-9]*)-.*","\\1",cleanGisticData$`Wide Peak Limits`)), 
                   end = as.numeric(gsub(".*:[0-9]*-([0-9]*)[:(:].*","\\1",cleanGisticData$`Wide Peak Limits`))),
)

refGenesToConsider <- refGenes[queryHits(findOverlaps(refGenes,grGistic))]

#Getting the genes that are affected by the CNAs
overlaps <- findOverlaps(grCnas, refGenesToConsider)
overlappingGenes <- refGenes[subjectHits(overlaps)]

#Collecting and adding gene symbols (names)
geneIDs<- mcols(overlappingGenes)$gene_id
geneSymbols <- mapIds(org.Hs.eg.db, keys = geneIDs, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
mcols(overlappingGenes)$geneSymbol <- geneSymbols

#Combining all the information in long format
allGenes <- data.table(
  patient = mcols(grCnas)$patient[queryHits(overlaps)],
  sample = mcols(grCnas)$sample[queryHits(overlaps)],
  chromosome = as.vector(seqnames(grCnas)[queryHits(overlaps)]),
  start = start(grCnas)[queryHits(overlaps)],
  end = end(grCnas)[queryHits(overlaps)],
  geneID = mcols(overlappingGenes)$gene_id,
  geneSymbol = geneSymbols,
  log2_ratio_adjusted = mcols(grCnas)$log2_ratio_adjusted[queryHits(overlaps)]
)

#Merge with patient data
sampleData <- fread(sampleDataFile)
varData <- fread(varFile)
sampleData <- merge(sampleData,varData,by=c("patient","label"),all=T)
sampleData[,`:=`(patient=as.character(patient))]
setnames(sampleData,old=c("label"),new=c("sample"))
setkey(sampleData,patient,sample)
allGenesWithPatientData <- merge(allGenes,sampleData)
rm(allGenes)

#compareCluster configuration from Belen
getCK = function(x, 
                 fun = "enrichGO",
                 ont = "BP",
                 pAdjustMethod = "BH", 
                 pvalueCutoff = 0.01, 
                 qvalueCutoff = 0.05,
                 minGSSize = 5,
                 maxGSSize = 500,
                 readable = T){
  
  compareCluster(geneCluster = x, 
                 fun = fun,
                 OrgDb= org.Hs.eg.db,
                 ont = ont,
                 pAdjustMethod = pAdjustMethod,
                 pvalueCutoff  = pvalueCutoff,
                 qvalueCutoff  = qvalueCutoff,
                 minGSSize = minGSSize,
                 maxGSSize = maxGSSize,
                 readable = readable)
}


#LVI vs. IBC in Discovery cohort
#First we select the samples we will use so that we do not have multiple samples per patient and stage
theseSamples <- sampleData[order(observedVariance),lapply(.SD,first),by=.(patient,Stage)][,.SD,keyby=.(patient,sample)] #Using variance
#theseSamples <- sampleData[sample(1:.N),lapply(.SD,first),by=.(patient,Stage)] #At random

lviIbcDiscGenes <- list(lviG=allGenesWithPatientData[theseSamples,][Cohort == "Disc" & Stage == "LVI" & log2_ratio_adjusted>0,geneID],
                        lviL=allGenesWithPatientData[theseSamples,][Cohort == "Disc" & Stage == "LVI" & log2_ratio_adjusted<0,geneID],
                        ibcG=allGenesWithPatientData[theseSamples,][Cohort == "Disc" & Stage == "IBC" & log2_ratio_adjusted>0,geneID],
                        ibcL=allGenesWithPatientData[theseSamples,][Cohort == "Disc" & Stage == "IBC" & log2_ratio_adjusted<0,geneID])
ck_LviIbcDisc <- getCK(lviIbcDiscGenes)
if(!is.null(ck_LviIbcDisc)){
  ck_LviIbcDisc.simpl.0.35 = simplify(ck_LviIbcDisc,cutoff=0.35, by="p.adjust", select_fun=min) ##adjust cutoff to reduce more or less the terms.
  dotplot(ck_LviIbcDisc, showCategory = 40)
}

#LVI Across subtypes in Discovery cohort
lviSubtypesDiscGenes <- list(idG=allGenesWithPatientData[theseSamples,][Cohort == "Disc" & Stage == "LVI" & LVI_subtype=="ID" & log2_ratio_adjusted>0,geneID],
                             prolG=allGenesWithPatientData[theseSamples,][Cohort == "Disc" & Stage == "LVI" & LVI_subtype=="Prol" & log2_ratio_adjusted>0,geneID],
                             idL=allGenesWithPatientData[theseSamples,][Cohort == "Disc" & Stage == "LVI" & LVI_subtype=="ID" & log2_ratio_adjusted<0,geneID],
                             prolL=allGenesWithPatientData[theseSamples,][Cohort == "Disc" & Stage == "LVI" & LVI_subtype=="Prol" & log2_ratio_adjusted<0,geneID])
ck_lviSubtypesDisc <- getCK(lviSubtypesDiscGenes)
if(!is.null(ck_lviSubtypesDisc)){
  ck_lviSubtypesDisc.simpl.0.35 <- simplify(ck_lviSubtypesDisc,cutoff=0.35, by="p.adjust", select_fun=min) ##adjust cutoff to reduce more or less the terms.
  dotplot(ck_lviSubtypesDisc.simpl.0.35, showCategory = 40)
}

#LVI Across subtypes
lviSubtypesGenes <- list(idG=allGenesWithPatientData[theseSamples,][Stage == "LVI" & LVI_subtype=="ID" & log2_ratio_adjusted>0,geneID],
                             prolG=allGenesWithPatientData[theseSamples,][Stage == "LVI" & LVI_subtype=="Prol" & log2_ratio_adjusted>0,geneID],
                             idL=allGenesWithPatientData[theseSamples,][Stage == "LVI" & LVI_subtype=="ID" & log2_ratio_adjusted<0,geneID],
                             prolL=allGenesWithPatientData[theseSamples,][Stage == "LVI" & LVI_subtype=="Prol" & log2_ratio_adjusted<0,geneID])
ck_lviSubtypes <- getCK(lviSubtypesGenes)
if(!is.null(ck_lviSubtypes)){
  ck_lviSubtypes.simpl.0.35 <- simplify(ck_lviSubtypes,cutoff=0.35, by="p.adjust", select_fun=min) ##adjust cutoff to reduce more or less the terms.
  dotplot(ck_lviSubtypes.simpl.0.35, showCategory = 40)
}
