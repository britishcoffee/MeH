# Loading packages
foo <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE , quietly = TRUE) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE , quietly = TRUE)
    }
  }
}
foo( c("argparser" , "roperators", "dplyr", "foreach","doParallel") )

#  Load function
MeH.t=function(vector,conditions,compare) {
  ind1<-which(conditions == compare[1])+3 # +3 for chrom,bin and strand columns
  ind2<-which(conditions == compare[2])+3
  #l=length(vector)
  vector=as.data.frame(vector)
  mean2=mean(as.numeric(vector[ind2]),na.rm=TRUE)
  mean1=mean(as.numeric(vector[ind1]),na.rm=TRUE)
  diff=mean2-mean1
  if(sd(vector[ind1])<1e-5 && sd(vector[ind2])<1e-5)
    return(data.frame(chrom=vector[1],pos=vector[2],delta=diff,pvalue=NaN,mean2=mean2,mean1=mean1,strand=vector[3]))
  else {
    out=t.test(vector[ind1],vector[ind2])
    return(data.frame(chrom=vector[1],pos=vector[2],delta=out$est[2]-out$est[1],pvalue=as.numeric(out$p.value),mean2=out$est[2],mean1=out$est[1],strand=vector[3]))
  }
}
findgene = function(position) {
  pr = as.numeric(args$p)
  chr=as.character(position[1])
  #message(chr)
  BP=as.numeric(position[2])
  #message(BP)
  St=as.character(position[,3])
  Gene=geneloc$gene[which((geneloc$TSS<=BP)*(geneloc$TES>=BP)*(as.character(geneloc$chrom)==chr)*(as.character(geneloc$strand)==as.character(St))==1)][1]
  if (St=='f') {
    promoter=geneloc$gene[which((geneloc$TSS-pr<=BP)*(geneloc$TSS+pr>=BP)*(as.character(geneloc$chrom)==chr)*(geneloc$strand=="f")==1)][1]
  }
  if (St=='r') {
    promoter=geneloc$gene[which((geneloc$TES-pr<=BP)*(geneloc$TES+pr>=BP)*(as.character(geneloc$chrom)==chr)*(geneloc$strand=="r")==1)][1]
  }
  return(list(chrom=chr,bin=BP,Gene=Gene,Promoter=promoter,strand=St))
}

# Create a parser
cat("[",format(Sys.time(), "%X"),"]","Start proccessing","\n")

p <- arg_parser("Find DHR")

# Add command line arguments
p <- add_argument(p,"-m","Meh", help="input Meh resultes csv file", type="character")
p <- add_argument(p,"-s","conditions", help="Define conditions of all samples", type="character")
p <- add_argument(p,"-g","genelist", help="The gene list")
p <- add_argument(p,"-o", "output", help="output gene list result csv file", default="DHR")
p <- add_argument(p,"-c", "core", help="the number of core for analysis", default="4")
p <- add_argument(p,"-p", "promoter", help="the region of promoter", default="1000", type="numeric")
p <- add_argument(p,"-adjp", "adjust pvalue", help="Select differential heterogeneous regions based on adjust pvalue", default="0.05", type="numeric")
p <- add_argument(p,"-pvalue", "pvalue", help="Select differential heterogeneous regions based on pvalue", default="0.05", type="numeric")
p <- add_argument(p,"-delta", "delta", help="Select differential heterogeneous regions based on delta", default="1.4", type="numeric")


# Parse the command line arguments
args <- parse_args(p)

CG <- read.csv(args$m,header=TRUE)

#### Remove rows with missing data
CG=CG[which(apply(CG,1,function(x) sum(is.na(x)))==0),]

#### Define conditions of all samples; i.e., A and B for 2 conditions, each with two replicates, samples 1 and 2 are replicates of A and samples 3 and 4 are replicates for B.This is for comparisons to be carried out later on
conditions <- unlist(strsplit(args$s, ","))
us <- unique(conditions)

#### Calculate t-statistics and p-values for all bins between user specified conditions
registerDoParallel(cores=args$c)
# Compare condition B with A
Comp1<-data.frame(foreach(i = 1:dim(CG)[1],.combine = rbind) %dopar%
                    MeH.t(CG[i,],conditions=conditions,c(us)))
Comp1$padj=p.adjust(Comp1$pvalue)

#### Select differential heterogeneous regions based on user specified conditions; i.e., p-value of 0.05 and delta of 1 (positive or negative)

Comp1$DHR <- (Comp1$padj<args$adjp)*(abs(Comp1$delta)>args$delta)
Comp1$DHR <- (Comp1$pvalue<args$pvalue)*(abs(Comp1$delta)>args$delta)
Comp1$DHR.up <- (Comp1$pvalue<args$pvalue)*(Comp1$delta>args$delta)
Comp1$DHR.down <- (Comp1$pvalue<args$pvalue)*(Comp1$delta<args$delta)

#### DHG analysis if bed file is given as .txt with each row representing a gene and consists of gene name, chromosome number, TSS, TES and strand as 'f' (forward) or 'r' (reverse)
geneloc<-read.table(args$g,header=TRUE)
colnames(geneloc)<-c("gene","chrom","strand","TSS","TES")
geneloc$strand[as.character(geneloc$strand)=="+"]<-"f"
geneloc$strand[as.character(geneloc$strand)=="-"]<-"r"
geneloc$strand<-as.character(geneloc$strand)
geneloc$gene<-as.character(geneloc$gene)

genelist<-as.data.frame(foreach(i = 1:dim(Comp1)[1],.combine = rbind) %dopar% findgene(Comp1[i,c("chrom","bin","strand")]))
genelist$strand[as.character(genelist$strand)=="1"]<-"f"
genelist$strand[as.character(genelist$strand)=="2"]<-"r"
genelist <- as.data.frame(genelist)
Comp1 <- as.data.frame(Comp1)
Result_whole<-merge(Comp1,genelist,by=c("chrom","bin","strand"))
Result_whole$Gene <- unlist(Result_whole$Gene)
Result_whole$Promoter <- unlist(Result_whole$Promoter)


#### Get the up/down regulted DHG gene/promoter lists ####
DHG_Genebodys_up<-unique(unlist(genelist[which(Comp1$DHR.up==1),"Gene"])[!is.na(unlist(genelist[which(Comp1$DHR.up==1),"Gene"]))])
DHG_Genebodys_down<-unique(unlist(genelist[which(Comp1$DHR.down==1),"Gene"])[!is.na(unlist(genelist[which(Comp1$DHR.down==1),"Gene"]))])
DHG_Promoter_up<-unique(unlist(genelist[which(Comp1$DHR.up==1),"Promoter"])[!is.na(unlist(genelist[which(Comp1$DHR.up==1),"Promoter"]))])
DHG_Promoter_down<-unique(unlist(genelist[which(Comp1$DHR.down==1),"Promoter"])[!is.na(unlist(genelist[which(Comp1$DHR.down==1),"Promoter"]))])


infile1 <- as.character(gsub(" ", "", paste(args$o,"_MeH_Result",".csv")))
infile2 <- as.character(gsub(" ", "", paste(args$o,"_DHR_Result",".csv")))

Result_whole <- Result_whole[,-c(9:11)]
write.csv(f,infile1,row.names=F)

result <- file(infile2)
writeLines(paste("DHG Genebodys up: ",paste(DHG_Genebodys_up,collapse= ', ')), result)
close(result)
write(paste("DHG Genebodys down: ",paste(DHG_Genebodys_down,collapse= ', ')),infile2,append=TRUE)
write(paste("DHG Promoter up: ", paste(DHG_Promoter_up,collapse= ', ')),infile2,append=TRUE)
write(paste("DHG Promoter down: ",paste(DHG_Promoter_down,collapse= ', ')),infile2,append=TRUE)


cat("[",format(Sys.time(), "%X"),"]","Done!","\n")

