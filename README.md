# QIIME2 Metagenomics Analysis with Nextflow

## Description
This is a [Nextflow](www.nextflow.io) version of a workflow for a simple QIIME2 analysis. A more sophisticated version of this tool is used by the MU BAC internally that leverages our specific HPC and database resources. 
In support of open source and open sciencewe are making this version publicly available. Due to time and effort limitations we can't offer support on this version though so it is released **"as is"** and will 
only be updated as time allows from BAC team members. 

## Use

To run this you will need a modern version of Nextflow. *version 21.0+*

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
Initial version written by @kstiers.
Updated by @lcoghill.
