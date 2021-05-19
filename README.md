

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
#### 'CG', 'CHG' and 'CHH' with window size of 4 cytosines, weighted degree kernel as for pairwise distances between methylation patterns and 8 cores parallel processing

```js
    python genome_scr.py -w 4 -c 8 --CG --CHG -d 2
```

### 6. Download DHR.R for subsequent analysis
