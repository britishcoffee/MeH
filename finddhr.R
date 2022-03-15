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
    return(data.frame(chrom=vector[1],pos=vector[2],delta=diff,pvalue=NaN,mean2=mean2,mean1=mean1))
  else {
    out=t.test(vector[ind1],vector[ind2])
    return(data.frame(chrom=vector[1],pos=vector[2],delta=out$est[2]-out$est[1],pvalue=as.numeric(out$p.value),mean2=out$est[2],mean1=out$est[1]))
  }
}
findgene = function(position) {
  chr=as.character(position[1])
  #message(chr)
  BP=as.numeric(position[2])
  #message(BP)
  St=as.character(position[3])
  Gene=geneloc$gene[which((geneloc$TSS<=BP)*(geneloc$TES>=BP)*(as.character(geneloc$chrom)==chr)*(as.character(geneloc$strand)==as.character(St))==1)][1]
  if (St=='f') {
    promoter=geneloc$gene[which((geneloc$TSS-1000<=BP)*(geneloc$TSS+1000>=BP)*(as.character(geneloc$chrom)==chr)*(geneloc$strand=="f")==1)][1]
  }
  if (St=='r') {
    promoter=geneloc$gene[which((geneloc$TES-1000<=BP)*(geneloc$TES+1000>=BP)*(as.character(geneloc$chrom)==chr)*(geneloc$strand=="r")==1)][1]
  }
  return(list(chrom=chr,bin=BP,Gene=Gene,Promoter=promoter,strand=St))
}

# Create a parser
p <- arg_parser("Find DHR")

# Add command line arguments
p <- add_argument(p,"-m","--Meh", help="input Meh resultes csv file", type="character")
p <- add_argument(p,"-s","--conditions", help="Define conditions of all samples", type="character")
p <- add_argument(p,"-r","--rep", help="The replication of samples", type="character")
p <- add_argument(p,"-o", "--output", help="output gene list result csv file", default="DHR")
p <- add_argument(p,"-c", "--core", help="the number of core for analysis", default="4")

# Parse the command line arguments
args <- parse_args(p)

CG <- read.csv(args$Meh,header=TRUE)

#### Remove rows with missing data
CG=CG[which(apply(CG,1,function(x) sum(is.na(x)))==0),]

#### Define conditions of all samples; i.e., A and B for 2 conditions, each with two replicates, samples 1 and 2 are replicates of A and samples 3 and 4 are replicates for B. This is for comparisons to be carried out later on
conditions <- rep(c(args$conditions),args$rep)

#### Calculate t-statistics and p-values for all bins between user specified conditions
registerDoParallel(cores=args$core)

# Compare condition B with A
Comp1<-data.frame(foreach(i = 1:dim(CG)[1],.combine = rbind) %dopar% 
                    MeH.t(CG[i,],conditions=conditions,rgs$conditions))
Comp1$padj=p.adjust(Comp1$pvalue)

#### Select differential heterogeneous regions based on user specified conditions; i.e., p-value of 0.05 and delta of 1 (positive or negative)
Comp1$DHR <- (Comp1$padj<0.05)*(Comp1$pvalue<0.05)*(abs(Comp1$delta)>1.4)

#### DHG analysis if bed file is given as .txt with each row representing a gene and consists of gene name, chromosome number, TSS, TES and strand as 'f' (forward) or 'r' (reverse)
geneloc<-read.table('genelist.txt',header=TRUE)
colnames(geneloc)<-c("gene","chrom","strand","TSS","TES")
geneloc$strand[as.character(geneloc$strand)=="+"]<-"f"
geneloc$strand[as.character(geneloc$strand)=="-"]<-"r"

genelist<-as.data.frame(foreach(i = 1:dim(Comp1)[1],.combine = rbind) %dopar% findgene(Comp1[i,c("chrom","bin","strand")]))

infile <- as.character(gsub(" ", "", paste(args$output,".csv")))

write.table(genelist,infile,col.names=T,row.names=F)


