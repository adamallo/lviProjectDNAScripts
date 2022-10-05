library(data.table)
library(ape)
library(phangorn)
library(ggplot2)
library(cowplot)
library(tidytree)
library(ggtree)

#Config: This needs to be filled up to replicate the experiment
inDir="" # outDir of runPhylip.sh
outDir=""
outputDirPlot=""


pattern="^([^_]*).phy$"

stage.col = c('NORMAL' = "#FDBF6F",
              "DCIS" = "#A6CEE3",
              "IBC" = "#B2DF8A",
              "LVI" = "#FB9A99",
              "MET" ="#CAB2D6")

#Functions

## My own plot tree function for easy modification
#' @param tree tree to print
#' @param bs bootstrap support
#' @param bspos location of support (one of "node" or "branch")
#' @param title plot title
#' @param digits ndigits of the support
plotTreeBS = function(tree,bs,bspos=c("node","branch"),title,digits=2,...) {
  plot(tree,...)
  if(bspos=="node") {
    nodelabels(round(bs,digits = digits))
  } else {
    drawSupportOnEdges(round(bs,digits = digits))
  }
  title(main=title)
}

## My own save tree function for easy modification
#' @param file file where to save the tree plot
#' @param tree tree to save
#' @param bs bootstrap support
#' @param bspos location of support (one of "node" or "branch")
#' @param title plot title
#' @param digits ndigits of the support
saveTreePlotBS = function(file,tree,bs,bspos=c("node","branch"),title,digits=2,...) {
  pdf(file = file)
  plot(tree,...)
  if(bspos=="node") {
    nodelabels(round(bs,digits = digits))
  } else {
    drawSupportOnEdges(round(bs,digits = digits))
  }
  title(main=title)
  dev.off()
}

## My own save tree function for easy modification
#' @param file file where to save the tree plot
#' @param tree tree to save
#' @param bs bootstrap support
#' @param xMargin number of breakpoints (x axis) to pad the right-side of the trees to avoid problems with labels
saveBsAsTreePlot = function(file,tree,bs,xMargin=25,base_height=5,...) {

  #Converting tree data to ggtree format
  treeForPlot=tree%>%
    as_tibble()%>%
    mutate(label=gsub("_"," ",label))%>%
    mutate(label=gsub("Met","MET",label))%>%
    mutate(stage=as.factor(gsub(" .*$","",label)))
  suppressMessages(expr = {treeForPlot=full_join(treeForPlot,tibble(node=1:Nnode(tree) + Ntip(tree), bootstrap = round(bs*100,digits=0)))})
  
  #Pre-calculating plot information
  manualXlim=max(node.depth.edgelength(tree))+xMargin
  
  #Tree plot
  outTree=ggtree(as.treedata(treeForPlot))+
    geom_tiplab(aes(fill=stage),geom = "label",label.size=NA)+
    geom_nodelab(aes(label=round(bootstrap, 2)), hjust=1.5, vjust=-0.5, size=4) +
    scale_fill_manual(values=stage.col,guide='none') +
    theme_tree2(plot.caption=element_text(size=12))+
    labs(caption="Breakpoints") +
    xlim(0,manualXlim)
  
  #Generating ancestral state estimates
  #ML ancestral tree reconstruction of locations using a symmetrical Mk model. 
  #We use the marginal ancestral states as empirical Bayesian posterior probabilities. Marginal = F means marginal estimation (weird, we all know, see http://blog.phytools.org/2015/05/about-how-acemarginaltrue-does-not.html)
  aceResult=ace(treeForPlot$stage[!is.na(treeForPlot$stage)],tree,type="discrete",marginal = F,ip=0.01,use.expm=T,use.eigen = F)


  #If the best combination of parameters I found fails, we try some alternative packages and initial values for the ML estimation. Not the most elegant, but nothing really wrong with it.
  if(is.nan(aceResult$se)){
    aceResult=ace(treeForPlot$stage[!is.na(treeForPlot$stage)],tree,type="discrete",marginal = F) #Default
  }
    if(is.nan(aceResult$se)){
    aceResult=ace(treeForPlot$stage[!is.na(treeForPlot$stage)],tree,type="discrete",marginal = F,use.expm=T,use.eigen = F) #Same method as first attempt but with default ip=0.1
  }
  if(is.nan(aceResult$se)){
    aceResult=ace(treeForPlot$stage[!is.na(treeForPlot$stage)],tree,type="discrete",marginal = F,ip=0.001,use.expm=T,use.eigen = F) #Same method as first attempt but with smaller ip
  }
  if(is.nan(aceResult$se)){
    aceResult=ace(treeForPlot$stage[!is.na(treeForPlot$stage)],tree,type="discrete",marginal = F,use.expm=F,use.eigen = F) #Pure Ape method that is not expected to work very well
  }
  if(is.nan(aceResult$se)){
    aceResult=ace(treeForPlot$stage[!is.na(treeForPlot$stage)],tree,type="discrete",marginal = F,ip=0.01,use.expm=F,use.eigen = F)  #Pure Ape method that is not expected to work very well with different ip
  }

  if(!is.nan(aceResult$se)){
    asData=as.data.frame(aceResult$lik.anc)
    asData$node = rownames(asData)
    #Generating the barplots with the ancestral state information
    suppressWarnings(expr = {asPlots=nodebar(asData, cols=1:(ncol(asData)-1),position="dodge",color=stage.col)})
    #Adding internal bar plots
    outTreeWithAs=inset(outTree,asPlots,x="node",width=0.075,height=0.1)
  }else{
    warning(paste0("Problems with the Ancestral State Estimation for file",file,".Generating the plot without ancestral states."))
    outTreeWithAs=outTree
  }

  #Saving the results
  save_plot(filename = file,plot=outTreeWithAs,base_height = base_height)
}

#MAIN

fileList = dir(inDir, pattern, full.names = TRUE)

for (thisFile in fileList) {
  #message(paste0("File: ",thisFile))
  patient=gsub(basename(thisFile),pattern=pattern,replacement="\\1")
  #btrees=read.tree(paste(sep="/",inDir,paste0(patient,"_bootstrap_dollop.trees")))
  btrees=read.tree(paste(sep="/",inDir,paste0(patient,"_bootstrap_boostrapDollopFixed.trees")))
  tree=read.tree(paste(sep="/",inDir,paste0(patient,"_dollopFull.trees")))
  data=read.phyDat(file=thisFile,type="USER",levels=c(0,1))
  
  if(class(tree)=="phylo"){tree=c(tree)}
  
  for (treei in 1:length(tree)) {
    thisTree=tree[[treei]]
    normalTips=grep("NORMAL*",thisTree$tip.label) #Tip.label number of the normal samples
    
    if(length(normalTips)==0) {
      #Unrooted
      print("We cannot root this tree because we do not have normal samples")
      bs=prop.clades(thisTree,btrees,rooted = F)/length(btrees)
      
      #WARNING: ACCTRAN IS ESTIMATING BRANCH LENGHTS USING REGULAR PARSIMONY, BUT THE TREE HAS BEEN ESTIMATED USING DOLLO PARSIMONY
      finalTreeBL=acctran(thisTree,data)
      saveBsAsTreePlot(file=paste(sep="/",outputDirPlot,paste0(patient,"_dolloAndAncestral_","t",treei,".pdf")),finalTreeBL,bs)
    } else {
      if (length(normalTips)==1) {
        #Rooted with 1 normal sample, using its name
        outgroup=thisTree$tip.label[normalTips]
        rootedTree=root(thisTree,outgroup = outgroup, resolve.root = T)
        rootedBS=root(btrees,outgroup = outgroup, resolve.root = T)
      } else {
        #Rooted with the MRCA of all the normal samples, using the node number
        outgroupTree=getMRCA(phy = thisTree,tip = thisTree$tip.label[normalTips])
        if (outgroupTree==length(thisTree$tip.label)+1) {
          #We can't re-root because the root is already valid (most probably, the Normal's ancestor is the current root). See https://www.mail-archive.com/r-sig-phylo@r-project.org/msg03805.html
          rootedTree=thisTree
        } else {
          rootedTree=root(thisTree,node = outgroupTree, resolve.root = T)
        }
        outgroupBS=sapply(btrees,FUN=function(x){getMRCA(phy = x,tip = thisTree$tip.label[normalTips])})
        rootedBS=lapply(1:length(btrees),FUN=function(x){
          if (outgroupBS[x]!=length(btrees[[x]]$tip.label)+1){
            root(phy = btrees[[x]],node = outgroupBS[x])
          } else {
            btrees[[x]]
          }
        }
        )
        class(rootedBS)="multiPhylo"
      }
      
      bs=prop.clades(rootedTree,rootedBS,rooted = T)/length(rootedBS)
      bs[1]=NA ##I don't want to show the BS of the root since that is imposed by us
      
      #WARNING: ACCTRAN IS ESTIMATING BRANCH LENGHTS USING REGULAR PARSIMONY, BUT THE TREE HAS BEEN ESTIMATED USING DOLLO PARSIMONY
      finalTreeBL=acctran(rootedTree,data)
      saveBsAsTreePlot(file=paste(sep="/",outputDirPlot,paste0(patient,"_dolloAndAncestral_","t",treei,".pdf")),finalTreeBL,bs)
    }
  }
}