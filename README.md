## System requirements

* python 2.7 +
* pandas package 0.21 +
* pysam package 0.16.0.1 +

### Can be fulfilled by running
```js
pip install MeHscr
```
## Methylation heterogeneity profiling

### 1. Download genome_scr.py
```js
wget https://raw.githubusercontent.com/britishcoffee/MeHscr/main/MeHscr.py
```
### 2. Open a folder named "MeHdata" under the same directory
```js
mkdir MeHdata
```
### 3. Place .bam and .bam.bai files of all samples you wish to obtain methylation heterogeneity profiles into folder MeHdata/
```js
scp [directory_to_bamfiles_of_all_samples].bam* ./MeHdata
# or within MeHdata/
ln -s [directory_to_bamfiles_of_all_samples].bam* ./
```
### 4. Also place .fa and .fa.fai of the reference genome into the folder
```js
scp [directory_to_reference_genome].fa* ./MeHdata
# or within MeHdata/
ln -s [directory_to_reference_genome].fa* ./
```
### 5. Run the program genome_scr.py (see examples below)

#### Examples

```ruby
# 'CG' only with window size of 4 cytosines and 4 cores parallel processing (default score is 
# pairwise-similarity-based method, default distance between methylation patterns is Hamming distance)
    python MeHscr.py -w 4 -c 4 --CG
```

```ruby
# 'CG', 'CHG' and 'CHH' with window size of 4 cytosines, weighted degree kernel for pairwise distances 
# between methylation patterns and 8 cores parallel processing
    python MeHscr.py -w 4 -c 8 --CG --CHG --CHH -d 2
```

### 6. Download DHR.R for subsequent analysis

#### Load required packages and functions
```R
install.packages("roperators")
library(roperators)
install.packages("dplyr")
library(dplyr)
install.packages("foreach")
library(foreach)

MeH.t = function(vector,conditions,compare) {
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

```
#### Load files for analysis by first setting the work directory to where your files are located
```R
setwd("~/MeHdata")
CG <- read.table('CG_Results.csv',header=TRUE,sep=",")
CHG <- read.table('CHG_Results.csv',header=TRUE,sep=",")
CHH <- read.table('CHH_Results.csv',header=TRUE,sep=",")
```

<img src="https://github.com/britishcoffee/Methylationhet/blob/main/READMEimages/image1.png?raw=true" width="600">

#### Define conditions of all samples; i.e., A and B for 2 conditions, each with two replicates, samples 1 and 2 are replicates of A and samples 3 and 4 are replicates for B. This is for comparisons to be carried out later on

```R
conditions <- c("A","A","B","B")
```

#### Calculate t-statistics and p-values for all bins between user specified conditions; An example is for A vs B here
```R
library(doParallel)
registerDoParallel(cores=4)
# Compare condition B with A
Comp1<-data.frame(foreach(i = 1:dim(CG)[1],.combine = rbind) %dopar% 
                      MeH.t(CG[i,],conditions=conditions,c("A","B")))
Comp1$padj=p.adjust(Comp1$pvalue)
```
#### Select differential heterogeneous regions based on user specified conditions; i.e., p-value of 0.05 and delta of 1.4 (positive or negative)
```R
Comp1$DHR <- (Comp1$padj<0.05)*(Comp1$pvalue<0.05)*(abs(Comp1$delta)>1.4)
```

<img src="https://github.com/britishcoffee/Methylationhet/blob/main/READMEimages/image6.png?raw=true" width="450">

#### DHG analysis if bed file is given as .txt with each row representing a gene and consists of gene name, chromosome, TSS, TES and strand as 'f' (forward) or 'r' (reverse)

```R
geneloc<-read.table('genelist.txt',header=TRUE)
colnames(geneloc)<-c("gene","chrom","strand","TSS","TES")
geneloc$strand[as.character(geneloc$strand)=="+"]<-"f"
geneloc$strand[as.character(geneloc$strand)=="-"]<-"r"
```
<img src="https://github.com/britishcoffee/Methylationhet/blob/main/READMEimages/image7.png?raw=true" width="300">

```R
genelist<-foreach(i = 1:dim(Comp1)[1],.combine = rbind) %dopar% findgene(Comp1[i,c("chrom","bin","strand")]) 
```

