#!/usr/bin/env Rscript

suppressMessages({
  library(QDNAseq)
  library(DNAcopy)
  library(data.table)
  library(aCGH)
  library(CGHcall)
  library(cowplot)
  library(tibble)
  library(dplyr)
  library(rascal)
})

#FUNCTIONS
##########

##Inaccessible QDNAseq functions
###############################
#QDNAseq function
sdDiffTrim <- function(x, ..., trim=0.001, scale=TRUE) {
  require(matrixStats)
  if (scale)
    x <- x / mean(x, na.rm=TRUE)
  sdDiff(x, ..., trim=trim)
}

#QDNAseq function
getNoise=function(thisdata)
{
  condition=thisdata@featureData$use
  counts=thisdata@assayData$counts[condition,,drop=FALSE]
  fit=thisdata@assayData$fit[condition,,drop=FALSE]
  
  counts=counts/thisdata@featureData$bases[condition]*100L
  counts[thisdata@featureData$bases[condition]==0]=0L
  
  signal=counts/fit
  signal[fit<=0]=0
  
  noise=apply(signal,2,sdDiffTrim, na.rm=TRUE)
  return(noise)
}

##Other functions
getRascalCopyNumbers=function(x,sample=NULL) {
  
  if (length(sampleNames(x)) != 1){
    if(is.null(sample)){
      stop("The current implementation of this package does not allow to use multiple samples at the same time. Call this function indicating a single sample using the sample argument")
    } else {
      x = x[,sample]
    }
  }
  
  sample=sampleNames(x)
  
  if (any(class(x) == "QDNAseqAneuploidyCopyNumbers")){
    copy_number_values = aneuploidy(x)
    segmented_values=segmentedAneuploidy(x)
  } else if (any(class(x) == "QDNAseqCopyNumbers")){
    copy_number_values =copynumber(x)
    segmented_values=segmented(x)
  } else {
    stop("ERROR: x is not a QDNAseqAneuploidyCopyNumbers or QDNAseqCopyNumbers object")
  }
  relative_copy_numbers=Biobase::fData(x) %>%
    rownames_to_column(var = "id") %>%
    as_tibble() %>%
    select(id, chromosome, start, end) %>%
    mutate_at(vars(start, end), as.integer) %>%
    mutate(chromosome = factor(chromosome, levels = unique(chromosome))) %>%
    mutate(sample = sample) %>%
    mutate(copy_number = copy_number_values) %>%
    mutate(segmented = segmented_values) %>%
    select(sample, chromosome, start, end, copy_number, segmented)
}

getRascalSegments=function(x,sample=NULL){
  if(any(class(x)=="tbl")){
    segments <- copy_number_segments(x) %>% mutate(log2ratio = log2(copy_number))
  } else {
    stop("ERROR: getRascalSegments expects the output of getRascalCopyNumbers as input")
  }
  segments
}

getCellularityPloidySolutions=function(x,min_ploidy,max_ploidy,ploidy_step,min_cellularity,max_cellularity,cellularity_step,distance_function,sample=NULL) {
  
  if (length(sampleNames(x)) != 1){
    if(is.null(sample)){
      stop("The current implementation of this package does not allow to use multiple samples at the same time. Call this function indicating a single sample using the sample argument")
    } else {
      x = x[,sample]
    }
  }
  
  sample=sampleNames(x)
  
  if (any(class(x) == "QDNAseqAneuploidyCopyNumbers")){
    copy_number_values = aneuploidy(x)
    segmented_values=segmentedAneuploidy(x)
  } else if (any(class(x) == "QDNAseqCopyNumbers")){
    copy_number_values =copynumber(x)
    segmented_values=segmented(x)
  } else {
    stop("ERROR: x is not a QDNAseqAneuploidyCopyNumbers or QDNAseqCopyNumbers object")
  }
  relative_copy_numbers=Biobase::fData(x) %>%
    rownames_to_column(var = "id") %>%
    as_tibble() %>%
    select(id, chromosome, start, end) %>%
    mutate_at(vars(start, end), as.integer) %>%
    mutate(chromosome = factor(chromosome, levels = unique(chromosome))) %>%
    mutate(sample = sample) %>%
    mutate(copy_number = copy_number_values) %>%
    mutate(segmented = segmented_values) %>%
    select(sample, chromosome, start, end, copy_number, segmented)
  
  segments <- copy_number_segments(relative_copy_numbers) %>% mutate(log2ratio = log2(copy_number))
  
  solutions <- find_best_fit_solutions(
    segments$copy_number, segments$weight,
    min_ploidy = min_ploidy, max_ploidy = max_ploidy, ploidy_step = 0.01,
    min_cellularity = min_cellularity, max_cellularity = max_cellularity, cellularity_step = 0.01,
    distance_function = "MAD")
  
}

selectSolution=function(x,strategy=c("best","highcell")) {
  if(strategy=="highcell") {
    selectedRow=x %>% arrange(desc(cellularity),distance) %>% slice(1)
  } else {
    selectedRow=x %>% arrange(distance,desc(cellularity)) %>% slice(1)
  }
  selectedRow
}

writeACNA=function(x,bestSolution,outputFile) {	
  ACNA_data<- x %>% select(
    chromosome,
    start,
    end,
    bins=bin_count,
    log2_ratio=log2ratio,
    relative_copy_number=copy_number
  ) %>%
    mutate(absolute_copy_number = relative_to_absolute_copy_number(relative_copy_number, bestSolution$ploidy, bestSolution$cellularity)) %>%
    mutate(across(c(log2_ratio, relative_copy_number, absolute_copy_number), round, digits = 3))
  write.table(ACNA_data,file=gsub("(.[^.]*$)",paste0(".c",round(bestSolution$cellularity,digits=3),".p",round(bestSolution$ploidy,digits=3),"\\1"),outputFile),row.names=F,col.names=T,quote=F)
}

writeACNAPlot=function(x,bestSolution,
                       outputFile,
                       transform=c("relative","log2"),
                       min_ploidy=1.25,
                       max_ploidy=8,
                       min_cellularity=0.2,
                       max_cellularity=1,
                       max_absolute_copy_number=16,
                       min_absolute_copy_number=0,
                       max_bins_display=10000,
                       bin_color="black",
                       segment_color="red",
                       absolute_copy_number_step_color="blue") {	
  transform=match.arg(transform)
  segmentsPlot=getRascalSegments(x) %>% mutate(relative_copy_number=copy_number)
  copyNumbersPlot=x %>% mutate(log2ratio = log2(copy_number)) %>% mutate(relative_copy_number=copy_number)
  copyNumberSteps= tibble(absolute_copy_number = min_absolute_copy_number:max_absolute_copy_number) %>%
    mutate(relative_copy_number = absolute_to_relative_copy_number(absolute_copy_number, bestSolution$ploidy, bestSolution$cellularity)) %>%
    mutate(log2ratio = log2(relative_copy_number)) %>%
    mutate(copy_number = relative_copy_number)
  
  ylabel=NULL
  if(transform=="log2"){
    ylabel=expression(log[2]~ratio)
    segmentsPlot %>% mutate(copy_number=log2ratio)
    copyNumbersPlot %>% mutate(copy_number=log2ratio)
    copyNumberSteps %>% mutate(copy_number=log2ratio)
  } else {
    ylabel="Relative copy number"
  }
  
  ylims=1.1 * quantile(copyNumbersPlot$copy_number[is.finite(copyNumbersPlot$copy_number)], c(0.001, 0.999), na.rm = TRUE) #ylimits calculated as they do in their app with my own modification to remove infinite values. This was generating problems with -Infs
  
  this_genome_plot=genome_copy_number_plot(
    copy_number=copyNumbersPlot,
    segments=segmentsPlot,
    copy_number_steps=copyNumberSteps,
    max_points_to_display = max_bins_display,
    min_copy_number = ylims[1], max_copy_number = ylims[2],
    point_colour = bin_color,
    segment_colour = segment_color,
    copy_number_step_colour = absolute_copy_number_step_color,
    xlabel = "Chromosome", ylabel = ylabel)
  
  save_plot(plot = this_genome_plot,filename = outputFile,base_height = 6)
}
##END FUNCTIONS
configFile <- paste(sep="/",Sys.getenv("lviProjectDNAScripts"),"configFile")
if(!file.exists(configFile)){
    stop("Configuration file configFile not found. Edit configFile.in to adapt it to your system, save it as configFile, and export the lviProjectDNAScripts environment variable before running this script")
}
source(configFile)

outDir <- acnaDir
dir.create(outDir)

##Segmentation parameters
###MergeLevels
pv.thres=1e-10

##CNA calling (Rascal) parameters
min_ploidy <- 1.25
max_ploidy <- 5.5
min_cellularity <- 0.2
max_cellularity <- 1
max_absolute_copy_number <- 12
min_absolute_copy_number <- 0
max_bins_display <- 100000 ##Bins are subsampled to this number for plotting
bin_color <- "black"
segment_color <- "red"
absolute_copy_number_step_color <- "blue"

if(exists("randomSeed") & is.numeric(randomSeed)){
	set.seed(randomSeed)
}

#Data IO and manipulation:
sampleData <- fread(originalSampleDataFile)
setnames(sampleData,old=c("Sample_ID","Case"),new=c("sampleID","patient"))

sampleData <- sampleData[,lapply(.SD,first),by=sampleID]#WARNING: there were repeated labels that I am removing here

sampleData[,`:=`(label=paste(sep="_",Stage,1:.N)),by=.(patient,Stage)]
sampleData[,`:=`(fullLabel=paste(sep=".",patient,label))]
write.table(sampleData,
            file=sampleDataFile,
            sep="\t",
            quote=F,
            row.names = F)

##Iters contains the information of the samples for which we will run the CNA calling algorithm
sampleData[,`:=`(file=paste0(inDir,"/",sampleID,"_50_readCounts.RData"))]
iters <- sampleData[,.(file,patient,label)]

##Main loop.
#There is fine-grained parallelism inside some of the calling routines, so for now I am executing this loop in series. 
#This is not the best strategy and the analysis would most probably be faster with coarse-grain parallelism at this level.
future::plan("multisession", workers=mc_cores)
allSolutions=NULL
varTable=data.table(patient=as.character(NULL),label=as.character(NULL),observedVariance=as.numeric(NULL),expectedVariance=as.numeric(NULL))

for (iter in 1:nrow(iters)){
  thisFile=iters[iter,file]
  thisLabel=iters[iter,label]
  thisPatient=iters[iter,patient]
  
  #Data parsing
  readCounts=readRDS(thisFile)
  sampleNames(readCounts)=paste(sep="_",thisPatient,thisLabel)
  
  #Pre-processing
  readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE,mappability = 50) 
  readCountsFiltered <- estimateCorrection(readCountsFiltered,maxIter = 100)
  
  #Variance estimation (observed and expected)
  varTable=rbind(varTable,data.table(patient=thisPatient,label=thisLabel,observedVariance=getNoise(readCountsFiltered)^2,expectedVariance=QDNAseq:::expectedVariance(readCountsFiltered)))
  
  #Estimate corrected copy numbers
  copyNumbers <- correctBins(readCountsFiltered,maxIter=100,method = "median") ##This will calculate the fit if needed and otherwise use the one pre-calculated by the user
  smoothRelativeCopyNumbers <- smoothOutlierBins(normalizeBins(copyNumbers))  #normalizeBins divides by the median by default, generating relative copy numbers (as desired)
  
  #Segment copy numbers. 
  segmentedCopyNumbers=normalizeSegmentedBins(segmentBins(smoothRelativeCopyNumbers))
  #plot(segmentedCopyNumbers)
  
  #Mergelevels
  mergedSegmentationLevels=QDNAseq:::log2adhoc(aCGH::mergeLevels(QDNAseq:::log2adhoc(QDNAseq:::copynumber(segmentedCopyNumbers)),QDNAseq:::log2adhoc(QDNAseq:::segmented(segmentedCopyNumbers)),pv.thres=pv.thres)$vecMerged,inv=T)
  segmented(segmentedCopyNumbers)=mergedSegmentationLevels
  #plot(segmentedCopyNumbers)
  
  
  #Optimization of cellularity and ploidy using Rascal
  #Get solutions
  solutions=getCellularityPloidySolutions(segmentedCopyNumbers,min_ploidy = min_ploidy, max_ploidy = max_ploidy, ploidy_step = 0.01, min_cellularity = min_cellularity, max_cellularity = max_cellularity, cellularity_step = 0.01, distance_function = "MAD")
  
  if (nrow(solutions) == 0) {
    warning(paste0("The sample ",thisLabel," of patient ",thisPatient," did not reach any valid solutions"))
    next
  }
  #Write them
  write.table(solutions%>%mutate_at(vars(ploidy, cellularity), round, digits = 2) %>% mutate_at(vars(distance), round, digits = 3),file=paste(sep="/",outDir,paste(sep=".",thisPatient,thisLabel,"solutions","tsv")),sep="\t",row.names=F,col.names=T,quote = F)
  
  #Loop through two strategies
  #Get the data in Rascal format to re-use inside the loop
  copyNumbers=getRascalCopyNumbers(segmentedCopyNumbers)
  segments=getRascalSegments(copyNumbers)
  strategies=c("best","highcell")
  
  for(strategy in strategies) {
    bestSolution=selectSolution(solutions,strategy=strategy)
    
    ##Add the best solution and strategy to a table
    allSolutions=bind_rows(allSolutions,bestSolution %>% mutate(strategy=strategy) %>% mutate(sample=thisLabel) %>% mutate(patient=thisPatient))
    
    #These functions re-do quite a lot of stuff and thus are not very efficient but I am abstracting most stuff to make the code more readable. Modify if needed.
    #Write absolute copy numbers
    writeACNA(segments,bestSolution,outputFile=paste(sep="/",outDir,paste(sep=".",thisPatient,thisLabel,strategy,"ACNAs","tsv")))
    
    #Write plot
    writeACNAPlot(copyNumbers,bestSolution,outputFile=paste(sep="/",outDir,paste(sep=".",thisPatient,thisLabel,strategy,"ACNAs","png")),max_bins_display = 100000)
  }
}

##Print the best-solutions table
write.table(allSolutions%>%mutate_at(vars(ploidy, cellularity), round, digits = 2) %>% mutate_at(vars(distance), round, digits = 3),file=paste(sep="/",outDir,paste(sep=".","allSolutions","tsv")),row.names=F,col.names=T,quote = F)
write.table(varTable,file=paste(sep="/",outDir,paste(sep=".","variance","tsv")),row.names=F,col.names=T,quote = F)
