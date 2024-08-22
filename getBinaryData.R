#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
  library(phangorn)
})

#Config IO
configFile <- paste(sep="/",Sys.getenv("lviProjectDNAScripts"),"configFile")
if(!file.exists(configFile)){
    stop("Configuration file configFile not found. Edit configFile.in to adapt it to your system, save it as configFile, and export the lviProjectDNAScripts environment variable before running this script")
}
source(configFile)
inDir <- acnaDir
outDir <- breakpointDir

#Config
pattern <- "^(.*)[:.:]([^.]*)[:.:]best[:.:]ACNAs[:.:]c(.*)[:.:]p(.*)[:.:]tsv$" ##[::] is R's way of scape something, similar to \ or \\ in most regex


#Parsing data
fileList = dir(inDir, pattern, full.names = TRUE)

#Calls
allCalls=rbindlist(lapply(fileList,function(x,pattern){theData=fread(x);theData[,`:=`(patient=gsub(basename(x),pattern=pattern,replacement="\\1"),
                                                                                      sample=gsub(basename(x),pattern=pattern,replacement="\\2"),
                                                                                      cellularity=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\3")),
                                                                                      ploidy=as.numeric(gsub(basename(x),pattern=pattern,replacement="\\4")))]},pattern=pattern))

badSamples <- data.table(patient=c("344","740"),sample=c("NORMAL_2","NORMAL_1"))
setkey(badSamples,patient,sample)
setkey(allCalls,patient,sample)

allCalls=allCalls[!badSamples,]

callsWithInfo=allCalls[,.(chrom=chromosome,start,end,originalCN=absolute_copy_number,sample_id=sample,Stage=gsub(sample,pattern="_.*$",replacement=""),patient=gsub(patient,pattern="PATIENT_",replacement = "P"),aPloidy=round(as.numeric(ploidy)))]

#Rounding
callsWithInfo[,`:=`(cn=round(originalCN))]
callsWithInfo[cn<0,`:=`(cn=0)]

#Adding them to the cn and cellularity/ploidy data
finalUnsyncData=callsWithInfo[,.(patient,sample_id,chrom,start,end,cn)]

#All present breakpoints for all samples
longOnes=melt(finalUnsyncData, id.vars= c("patient","sample_id","chrom"), measure.vars=c("start","end"),value.name="pos")
longOnes[,`:=`(variable=NULL,state=1)]

setkey(longOnes,chrom,pos,sample_id)

#Preparing output directory
dir.create(outDir,recursive = T)

#Now, per patient, we need to fill-up the 0s (existing breakpoints that are not present)
for (thisPatient in longOnes[,unique(patient)]){
  #print(patient)
  ##Breakpoints in this patient, we start with all 0s and collapsing shared ones
  allBreakPointsThisPatient=longOnes[patient==thisPatient,.(state=0),by=.(chrom,pos)]
  
  ##Breakpoints in this patient, replicated in all samples, we start all with 0s here too
  allBreakPointsThisPatientBySample=rbindlist(
    lapply(
      longOnes[patient==thisPatient,unique(sample_id)],
      FUN=function(x){allBreakPointsThisPatient[,.(sample_id=x,chrom,pos,state)]}
    )
  )
  setkey(allBreakPointsThisPatientBySample,chrom,pos,sample_id)
  
  ##Substitute 0s for 1s in the corresponding positions
  allBreakPointsThisPatientBySample[longOnes[patient==thisPatient],`:=`(state=1)]
  
  #Just a quick check
  #hist(allBreakPointsThisPatientBySample[,.(sum(state)),by=.(chrom,pos)]$V1)
  
  ##Long to wide
  finalDataPseudoPhylip=dcast(allBreakPointsThisPatientBySample,chrom + pos ~ sample_id,value.var = "state")
  finalDataPseudoPhylip[,`:=`(chrom=NULL,pos=NULL)]
  
  ##Using phangorn to get the data
  finalPhyDatData=as.phyDat(finalDataPseudoPhylip,type="USER",levels=c(0,1))
  
  ##Write the data in phylip
  ##write.phyDat does not do much, basically uses ape's write.dna. Ape's write.dna does not do anything smart about the phylip names, only sets all at the max length of the names and adds 1 character.
  ##This is wrong, since real Phylip's names should be 10 characters long.
  ##I am fixing it here
  newNames=sprintf("%-9s",names(finalPhyDatData))
  if(length(unique(newNames))!=length(unique(names(finalPhyDatData)))) stop(paste0("ERROR: some labels for patient ",thisPatient," are not compatible with the Phylip format"))
  
  names(finalPhyDatData)=newNames
  
  write.phyDat(finalPhyDatData,paste(sep="/",outDir,paste0(thisPatient,".phy")),format="phylip")
  
}
