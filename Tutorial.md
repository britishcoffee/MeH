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

* CG_Results.csv

```sh
Rscript tobed.R -m /MeHdata/CG_Results.csv -o /MeHdata/Meh_result.bedGraph
```

##### output

* Meh_result.bedGraph



finding DHR

##### input

* Results.csv

```sh
Rscript finddhr.R -m /MeHdata/Results.csv -s A,B -r 2 -c 4
```

##### output

* DHR.csv

