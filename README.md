

## Methylation heterogeneity profiling

### 1. Download genome_scr.py

### 2. Open a folder named "MeHdata" under the same directory

### 3. Place .bam and .bam.bai files of all samples you wish to obtain methylation heterogeneity profiles into folder MeHdata/

### 4. Also place .fa and .fa.fai of the reference genome into the folder

### 5. Run the program genome_scr.py by using one of the following commands

#### 'CG' only with window size of 4 cytosines and 4 cores parallel processing

```js
    python genome_scr.py -w 4 -c 4 --CG
```
#### 'CG', 'CHG' and 'CHH' with window size of 4 cytosines, weighted degree kernel for pairwise distances between methylation patterns (default is Hamming distance if unspecified) and 8 cores parallel processing

```js
    python genome_scr.py -w 4 -c 8 --CG --CHG --CHH -d 2
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

MeH.t=function(vector,conditions,compare) {
  ind1<-which(conditions == compare[1])+3 # +2 for chrom,bin and strand columns
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
```
#### Load files for analysis by first setting the work directory to where your files are located
```R
setwd("~/MeHdata")
Dest <- read.table('CG_Results.csv',header=TRUE,sep=",")
Dest <- read.table('CHG_Results.csv',header=TRUE,sep=",")
Dest <- read.table('CHH_Results.csv',header=TRUE,sep=",")
```

![alt text](https://github.com/britishcoffee/Methylationhet/blob/main/image1.png?raw=true)

#### Remove rows with no data
```R
Dest=Dest[which(apply(Dest,1,function(x) sum(is.na(x)))==0),]
```

![alt text](https://github.com/britishcoffee/Methylationhet/blob/main/image2.png?raw=true)

#### Construct bins of 400bp (default, can be changed to any even number) and add a column of bin position to data
```R
bin_size=400
Dest$bin<-((Dest$pos-1) %/% bin_size)*bin_size+(bin_size/2)
```

![alt text](https://github.com/britishcoffee/Methylationhet/blob/main/image3.png?raw=true)

#### Obtain names of samples (will be identical to names of the bam files provided if unchanged)
```R
samples=colnames(Dest)[which(!colnames(Dest) %in% c("chrom","pos","strand","bin"))]
```

![alt text](https://github.com/britishcoffee/Methylationhet/blob/main/image4.png?raw=true)


#### Obtain results (average methylation heterogeneity for the bins) by taking averages of methylation heterogeneity for windows within the same bins
```R
# assign("[samplename]",data.frame(Dest %>% group_by(chrom,bin,strand) %>% summarise(mean(samplename))))
assign("C0031test1234",data.frame(Dest %>% group_by(chrom,bin,strand) %>% summarise(mean(C0031test1234))))
assign("C0033test1234",data.frame(Dest %>% group_by(chrom,bin,strand) %>% summarise(mean(C0033test1234))))
assign("C0035test1234",data.frame(Dest %>% group_by(chrom,bin,strand) %>% summarise(mean(C0035test1234))))
assign("C0037test1234",data.frame(Dest %>% group_by(chrom,bin,strand) %>% summarise(mean(C0037test1234))))
```

#### Merging results from different samples into a matrix called "new"
```R
new=c()
for (s in samples) {
    data=get(s)
    colnames(data)[4]=s
    if (!is.null(new)) new = merge(new,data,by=c("chrom","bin","strand"))
    else new=data
}
```

![alt text](https://github.com/britishcoffee/Methylationhet/blob/main/image5.png?raw=true)

#### Define conditions of all samples; i.e., A and B for 2 conditions, each with two replicates, samples 1 and 2 are replicates of A and samples 3 and 4 are replicates for B. This is for comparisons to be carried out later on

```R
conditions <- c("A","A","B","B")
```
#### Calculate t-statistics and p-values for all bins between user specified conditions
```R
library(doParallel)
registerDoParallel(cores=4)
# Compare condition B with A
Comp1<-data.frame(foreach(i = 1:dim(new)[1],.combine = rbind) %dopar% 
                      MeH.t(new[i,],conditions=conditions,c("A","B")))
```
#### Select differential heterogeneous regions based on user specified conditions; i.e., p-value of 0.05 and delta of 1 (positive or negative)
```R
Comp1$DHR <- (Comp1$pvalue<0.05)*(abs(Comp1$delta)>1)
```

![alt text](https://github.com/britishcoffee/Methylationhet/blob/main/image6.png?raw=true)

#### DHG analysis if bed file is given as .txt with each row representing a gene and consists of gene name, chromosome number, TSS, TES and strand as 'f' (forward) or 'r' (reverse)

```R
geneloc<-read.table('../hg19plusgenes.txt',header=FALSE)
```
![alt text](https://github.com/britishcoffee/Methylationhet/blob/main/image7.png?raw=true)
```R
genelist<-foreach(i = 1:dim(Comp1)[1],.combine = rbind) %dopar% findgene(Comp1[i,c(1,2,3)]) # same for all scores
```

