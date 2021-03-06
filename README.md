# QIIME2 Metagenomics Analysis with Nextflow

## Description
This is a [Nextflow](www.nextflow.io) version of a workflow for a simple QIIME2 analysis. A more sophisticated version of this tool is used by the MU BAC internally that leverages our specific HPC and database resources. 

In support of open source and open science, we are making this version publicly available. Due to time and effort limitations we can't offer support on this version though so it is released **"as is"** and will only be updated as time allows from BAC team members. 

In short, ***you are mostly on your own with this, so enjoy and hack about, but good luck!*** 

## Use

#### Requirements
- Nextflow - tested with version 21.0+
- QIIME2 - tested with version (2021.8.0) 
- QIIME2 formatted classifiers - We built one ours using the naive-bayes model, the 16S classifier is based on the Silva (v.138) 99% reference database, targeted for PE250 reads. The ITS classifier is based on the Unite (v.10.05.2021) 99% reference database.

You will need to create a **Manifest File**, which is just a simple csv file, with your data using the format below:

```
sample-id,absolute-filepath,direction
```

As an example:
```
JK001,Data1_L001_R1_001.fastq.gz,forward
```

### Known Limitations

Additional metadata isn't currently supported. That is coming in future versions.

### Credit
Initial version written by [@KyleStiers](https://github.com/KyleStiers).  
Initial update and release by [@lcoghill](https://github.com/lcoghill).
