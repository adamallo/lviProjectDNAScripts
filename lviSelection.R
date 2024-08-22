suppressMessages({
  library(data.table)
})

configFile <- paste(sep="/",Sys.getenv("lviProjectDNAScripts"),"configFile")
if(!file.exists(configFile)){
    stop("Configuration file configFile not found. Edit configFile.in to adapt it to your system, save it as configFile, and export the lviProjectDNAScripts environment variable before running this script")
}
source(configFile)


inputFile <- gisticInputFile
outputFile <- gisticLVIInputFile
varFile <- paste(sep="/",acnaDir,"variance.tsv")

theseData <- fread(inputFile)
theseVariances <- fread(varFile)
theseVariances[,`:=`(Sample = paste(sep=".",patient,label),Stage = gsub(pattern = "_[0-9]*$",replacement = "",x = label))]

obsVarSelection <- theseVariances[Stage=="LVI",][order(observedVariance),lapply(.SD,first),by=.(patient)]
#relVarSelection <- theseVariances[Stage=="LVI",][,`:=`(relVariance = observedVariance/expectedVariance)][order(observedVariance/expectedVariance),lapply(.SD,first),by=.(patient)]

write.table(file = outputFile, x = theseData[Sample %in% obsVarSelection$Sample,], quote = F, sep = "\t", row.names = F)
