

<img src="https://github.com/britishcoffee/Methylationhet/blob/main/READMEimages/MeHscr.png?raw=true" width="300">

# MeH :sheep:

:mega:Genomewide methylation heterogeneity and differential heterogeneity analysis.


### Publication

Estimating methylation heterogeneity in bisulfite sequencing data by mathematical modelling

## Pipeline

*** figure here***

### Documentation

MeH users guide is available as a [PDF file](./Manual.pdf), containing the detail of each step. For questions please open an issue on [GitHub](https://github.com/britishcoffee/MeHscr/issues).

##  Table of Contents

   * [System requirements](#System requirements)
   * [Installation](#Installation)
   * [Methylation heterogeneity profiling](#Methylation heterogeneity profiling)
         * [Examples](#Examples)
   * [Subsequent analysis](#Subsequent analysis)
         * [Example](#Example)

## System requirements
* python 2.7 + 
* pandas package 0.24 +
* pysam package 0.16.0.1 +
* joblib package

## Installation

MeH can be installed for Linux, macOS, or Windows by either compiling  from source which has the advantage that it will be optimized to the specific system:

```bash
git clone https://github.com/britishcoffee/MeHscr.git
cd MeHscr
```
## Methylation heterogeneity profiling
Use the scrpit **MeHscr.py** to calculated the methylation heterogeneity.

**Input**

* Run all the files under folder "**MeHdata**", including:
  * .bam and .bam.bai files
  * .fa and .fa.fai of the reference genome 



**Useage** 

```ruby
$ python MeHscr.py -h
	
  usage: MeHscr.py [-h] [-w WINDOWSIZE] [-c CORES] [-m MEH] [-d DIST] [--CG]
                   [--CHG] [--CHH] [--opt] [--mlv]

  optional arguments:
    -h, --help            show this help message and exit
    -w WINDOWSIZE, --windowsize WINDOWSIZE
                          number of CGs
    -c CORES, --cores CORES
                          number of cores
    -m MEH, --MeH MEH     Methylation heterogeneity score 1:Abundance 2:PW
                          3:Phylogeny [Default: 2]
    -d DIST, --dist DIST  Distance between methylation patterns 1:Hamming 2:WDK [Default: 1]
    --CG                  Include genomic context CG
    --CHG                 Include genomic context CHG
    --CHH                 Include genomic context CHH
    --opt                 Outputs compositions of methylation patterns
    --mlv                 Outputs methylation levels

```

##### Examples

```ruby
# 'CG' only with window size of 4 cytosines and 4 cores parallel processing (default score is pairwise-similarity-based method, default distance between methylation patterns is Hamming distance)
python MeHscr.py -w 4 -c 4 --CG
# 'CG', 'CHG' and 'CHH' with window size of 4 cytosines, weighted degree kernel for pairwise distances between methylation patterns and 8 cores parallel processing
python MeHscr.py -w 4 -c 8 --CG --CHG --CHH -d 2
```

> The programme is running at folder "/MeHdata"

**Output**

* MeHscreening.log 

```
Sample AT31test has coverage 5240 for context CG out of data coverage 192834
Sample AT33test has coverage 5236 for context CG out of data coverage 193431
Sample AT35test has coverage 5203 for context CG out of data coverage 192548
Sample AT37test has coverage 5233 for context CG out of data coverage 192694
```

*  /MeHdata/sample.0.csv files for each sample

```bash
## CG_AT31test_0.csv in the example
chrom,pos,MeH,dis,strand
1,511,1.41421,139,f
1,791,2.7161,114,r
1,810,3.69631,102,r
1,840,4.11599,109,r
```

> Format desctiptions:
>
> (1) chromsome
> (2) position
> (3) Methlyation heterogeneity
> (4) distance  between methylation patterns
> (5) strand as 'f' for forward or 'r'  for reverse

*  /MeHdata/Results.csv files for summary results

```bash
## CG_Results.csv in the example
chrom,bin,strand,AT31test,AT33test,AT37test,AT35test
1,600,f,1.41421,4.42434,1.97092,2.219035
1,600,r,2.7161,2.59751,3.62414,2.79942
1,1000,r,3.90615,4.90306,6.5213,4.0907849999999994
1,2600,r,0.0,0.707105,0.0,0.0
```

> Format desctiptions:
>
> (1) chromsome
> (2) bin size
> (3) strand
> (4)-(6) Methlyation heterogeneity for each sample


## Subsequent analysis

Use the function of scrpit **DHR.R** to find differentailly heterogeneity regions.

##### Required packages

```R
# install.packages("roperators")
library(roperators)
# install.packages("dplyr")
library(dplyr)
# install.packages("foreach")
library(foreach)
# install.packages("doParallel")
library(doParallel)
```

##### Required Functions

```R
MeH.t=function(vector,conditions,compare) {
  ind1<-which(conditions == compare[1])+3 
  ind2<-which(conditions == compare[2])+3
  vector=as.data.frame(vector)
  mean2=mean(as.numeric(vector[ind2]),na.rm=TRUE)
  mean1=mean(as.numeric(vector[ind1]),na.rm=TRUE)
  diff=mean2-mean1
  if(sd(vector[ind1])<1e-5 && sd(vector[ind2])<1e-5) 
    return(data.frame(chrom=vector[1],pos=vector[2],strand=vector[3],delta=diff,pvalue=NaN,mean2=mean2,mean1=mean1))
  else {
    out=t.test(vector[ind1],vector[ind2])
    return(data.frame(chrom=vector[1],pos=vector[2],strand=vector[3],delta=out$est[2]-out$est[1],pvalue=as.numeric(out$p.value),mean2=out$est[2],mean1=out$est[1]))
  }
}

findgene = function(position) {
  chr=as.character(position[,1])
  #message(chr)
  BP=as.numeric(position[,2])
  #message(BP)
  St=as.character(position[,3])
  Gene=geneloc$gene[which((geneloc$TSS<=BP)*(geneloc$TES>=BP)*(as.character(geneloc$chrom)==chr)*(as.character(geneloc$strand)==as.character(St))==1)][1]
  #user can define theie own promoter region [default: 1000]
  if (St=='f') {
    promoter=geneloc$gene[which((geneloc$TSS-1000<=BP)*(geneloc$TSS+1000>=BP)*(as.character(geneloc$chrom)==chr)*(geneloc$strand=="f")==1)][1]
  }
  if (St=='r') {
    promoter=geneloc$gene[which((geneloc$TES-1000<=BP)*(geneloc$TES+1000>=BP)*(as.character(geneloc$chrom)==chr)*(geneloc$strand=="r")==1)][1]
  }
  return(list(chrom=chr,bin=BP,Gene=Gene,Promoter=promoter,strand=St))
}
```

##### Input

* Results.csv files for summary results
* genelist.txt

> genelist.txt can be modified based on gene.gff file consists of gene, chromosome, TSS, TES, and strand.

##### Example

1. Load files for analysis by first setting the work directory to where your files are located

```R
CG <- read.csv('MeHdata/CG_Results.csv',header=TRUE)
CG=CG[which(apply(CG,1,function(x) sum(is.na(x)))==0),]
```

```R
> head(CG)
  chrom  bin strand  AT31test  AT33test AT37test  AT35test
1     1  600      f 1.4142100 4.4243400  1.97092 2.2190350
2     1  600      r 2.7161000 2.5975100  3.62414 2.7994200
3     1 1000      r 3.9061500 4.9030600  6.52130 4.0907850
4     1 2600      r 0.0000000 0.7071050  0.00000 0.0000000
5     1 3800      f 0.3304952 0.2571291  0.00000 0.1844622
6     1 4200      f 0.0000000 0.0000000  0.00000 0.0000000
```

2. Define conditions of all samples

```R
# An example is for A vs B here
conditions <- c("A","A","B","B")
```

3. Calculate t-statistics and p-values for all bins between user specified conditions

```R
registerDoParallel(cores=4)
# Compare condition B with A
Comp1<-data.frame(foreach(i = 1:dim(CG)[1],.combine = rbind) %dopar% 
                      MeH.t(CG[i,],conditions=conditions,c("A","A","B","B")))
Comp1$padj=p.adjust(Comp1$pvalue)
stopImplicitCluster()
```

4. Select differential heterogeneous regions based on user specified conditions

```R
#  i.e., p-value of 0.05 and delta of 1.4 (positive or negative)
Comp1$DHR <- (Comp1$padj<0.05)*(abs(Comp1$delta)>1.4)
Comp1$DHR <- (Comp1$pvalue<0.05)*(abs(Comp1$delta)>1.4)
Comp1$DHR.up <- (Comp1$pvalue<0.05)*(Comp1$delta>1.4)
Comp1$DHR.down <- (Comp1$pvalue<0.05)*(Comp1$delta<(-1.4))
```

```R
> head(Comp1)
  chrom  bin strand delta pvalue     mean2     mean1 padj DHR DHR.up DHR.down
1     1  600      f     0      1 2.9192750 2.9192750    1   0      0        0
2     1  600      r     0      1 2.6568050 2.6568050    1   0      0        0
3     1 1000      r     0      1 4.4046050 4.4046050    1   0      0        0
4     1 2600      r     0      1 0.3535525 0.3535525    1   0      0        0
5     1 3800      f     0      1 0.2938122 0.2938122    1   0      0        0
6     1 4200      f     0    NaN 0.0000000 0.0000000  NaN  NA     NA       NA
```



5. DHG analysis if bed file is given as .txt with each row representing a gene and consists of gene name, chromosome, TSS, TES and strand

```R
geneloc <- read.table('MeHdata/genelist.txt',header=T)
colnames(geneloc) <- c("gene","chrom","TSS","TES","strand")
geneloc$strand<-as.character(geneloc$strand)
#geneloc$strand[as.character(geneloc$strand)=="+"] <- "f"
#geneloc$strand[as.character(geneloc$strand)=="-"] <- "r"
geneloc$gene<-as.character(geneloc$gene)
```
```R
> head(geneloc)
            gene chrom  TSS  TES strand
1           gene     1 3631 5899      f
2           mRNA     1 3631 5899      f
3        protein     1 3760 5630      f
4           exon     1 3631 3913      f
5 five_prime_UTR     1 3631 3759      f
6            CDS     1 3760 3913      f
```

6. DHG analysis if bed file is given as .txt with each row representing a gene and consists of gene name, chromosome, TSS, TES and strand

```R
genelist <- foreach(i = 1:dim(Comp1)[1],.combine = rbind) %dopar% findgene(Comp1[i,c("chrom","bin","strand")]) 
```

```R
> genelist
           chrom bin   Gene      Promoter     strand
result.1   "1"   600   "NA"      "NA"         "f"
result.2   "1"   600   "NA"      "NA"         "r"
result.3   "1"   1000  "NA"      "NA"         "r"
result.4   "1"   2600  "NA"      "NA"         "r"
result.5   "1"   3800  "NA"      "NA"         "f"
```



7. DHG analysis 

```R
DHG_Genebodys_up<-unique(unlist(genelist[which(Comp1$DHR.up==1),"Gene"])[!is.na(unlist(genelist[which(Comp1$DHR.up==1),"Gene"]))])
DHG_Genebodys_down<-unique(unlist(genelist[which(Comp1$DHR.down==1),"Gene"])[!is.na(unlist(genelist[which(Comp1$DHR.down==1),"Gene"]))])
DHG_Promoter_up<-unique(unlist(genelist[which(Comp1$DHR.up==1),"Promoter"])[!is.na(unlist(genelist[which(Comp1$DHR.up==1),"Promoter"]))])
DHG_Promoter_down<-unique(unlist(genelist[which(Comp1$DHR.down==1),"Promoter"])[!is.na(unlist(genelist[which(Comp1$DHR.down==1),"Promoter"]))])
```



##### Output

```R

```



## Contact

Sabrina -:email:  ytchang.sabrina@gmail.com

