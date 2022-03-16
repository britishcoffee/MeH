# Tutorial

### Installation

#### Download from github

```sh
git clone https://github.com/britishcoffee/MeH.git
cd MeH
```



### Run example

#### Methylation heterogeneity profiling

The example data is in MeHdata folder.

##### input

* Methylomes: bam file and bam.bai
* Reference genome: fa and fai

```sh
python MeHscr.py -w 4 -c 4 --CG
```

##### output

* MeHscreening.log 
* /MeHdata/sample.0.csv files for each sample
  * CG_AT31test_0.csv
  * CG_AT33test_0.csv
  * CG_AT35test_0.csv
  * CG_AT37test_0.csv
* /MeHdata/CG_Results.csv files for summary results



#### Subsequent analysis

Transform Results.csv to .bedGraph for IGV visualisation

##### input

* MeHdata/CG_Results.csv files

```sh
Rscript tobed.R -m ./MeHdata/CG_Results.csv 
```

##### output

* Meh_result.bedGraph
  * PW_AT31test.bedGraph
  * PW_AT33test.bedGraph
  * PW_AT35test.bedGraph
  * PW_AT37test.bedGraph



finding DHR

##### input

* MeHdata/CG_Results.csv files
* genelist.txt for match the gene

```sh
Rscript finddhr.R -m ./MeHdata/CG_Results_test.csv -g ./MeHdata/genelist.txt -o ./MeHdata/CG -s W,W,D,D -p 1000 
```

##### output

* CG_DHR_Result.csv shows the list of DHR in down/up regulated genes/promoters
* CG_MeH_Result.csv shows the table with the differnece of MeH regions

