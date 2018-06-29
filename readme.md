# AgRenSeq\_GLM

## Description
AgRenSeq\_GLM is an extension of [AgRenSeq](https://github.com/steuernb/AgRenSeq), a pipeline to identify candidate resistance (_R_) genes in plants directly from a diversity panel. Please refer to the original pipeline for all the steps. Only the modifications made to the original pipeline are described below. These modifications pertain to the following steps of the original pipeline:
### 5: Create phenotype file
This step remains the same, except for two differences:

* By default, the input phenotype scores are assumed to be numeric values, as in the previous pipeline. However, it is no longer required that the scores be negative for susceptible and positive for resistant. The only requirement is that a higher score represents a higher resistance. 
* The phenotype scores can now be directly fed in the form of the Stackman's IT by setting the optional flag `--stackman` (refer to the usage and the parameters below). In this case, the input scores are first converted from the Stackman's IT to the AgRenSeq scores as specified in the original [AgRenSeq](https://github.com/steuernb/AgRenSeq) pipeline, before further analysis is done.

### 9: Generate association scores of _k_-mers and project onto the denovo assembly
This step has been re-implemented here in Python, and altered to incorporate:

* Correlation pre-filtering, to retain only those _k_-mers whose correlation to the phenotype is above a certain threshold. The correlation threshold is assumed to be 0.0 by default, but it can be set to any value between -1 and 1 by specifying the argument for the optional parameter `--correlationthreshold` (refer to the parameters below).
* Regression analysis to generate the measures of association which have been corrected for population structure by using the PCA dimensions of SNP markers matrix as covariates. The number of significant PCA dimensions used is set to 3 by default, but it can be set to any positive number (less than the number of SNP markers) by specifying the argument for the optional parameter `--pcadimensions` (refer to the parameters below).

Note that, unlike the original pipeline, the current version accepts the presence/absence matrix of _k_-mers only in the `gzipped` format.

Besides the above modifications to the steps described in the original pipeline, an additional step is added for sample size power study.

### 11: Random subsampling
A random subsample of the (usable) accessions is retained for the association analysis. The size of the subsample is specified as the argument for the optional parameter `--subsample` (refer to the parameters below).

For a given sample size, random subsampling - followed by association analysis - can be performed a number of times. The proportion of the runs where the gene(s) of interest are detected as the topmost candidate(s) gives an estimate of the power.


## Pre-requisites

### Python 3 and above

The code has been tested in Python 3.5.3. 

The following Python modules are required:

* `numpy` 
* `pandas` 
* `Biopython`: to parse assembly file
* `scikit-learn`: to compute PCA from SNP markers matrix
* `statsmodels`: for regression analysis
* `bitarray`


## Usage

```
python RunAssociation_GLM.py -i presenceMatrix.txt.gz -n nlr.txt -a assembly.fasta -p phenotype.txt -st -u usable.txt -s snp.txt -o output.txt
```

#### Parameters

Parameter (long version)| Parameter (short version) | Argument | Description
--- | --- | --- | ---
--inputmatrix | -i | presenceMatrix.txt.gz | Mandatory. The path to file containing the gzipped version of presence/absence matrix of k-mers created in step 4 of the [AgRenSeq](https://github.com/steuernb/AgRenSeq) pipeline.
--nlr | -n | nlr.txt | Mandatory. The path to file containing the list of contigs associated with nlrs created in step 7 of the [AgRenSeq](https://github.com/steuernb/AgRenSeq) pipeline.
--assembly | -a | assembly.fasta | Mandatory. The path to assembly file of accession onto which k-mers are mapped for plotting, created in step 6 of the [AgRenSeq](https://github.com/steuernb/AgRenSeq) pipeline.
--phenotype | -p | phenotype.txt | Mandatory. The path to phenotype file.
--usable | -u | usable.txt | Optional. The path to file containing the list of usable accessions.
--stackman | -st | snp.txt | Mandatory. The path to file containing matrix of SNP markers to compute PCA and correct for population structure.
--subsample | -sub | integer | Optional. Run association analysis by taking a random subsample of _this size_ from the given accessions.
--permute | -per |  | Optional. Permute the phenotype scores.
--pcadimensions | -dim | integer | Default 3. The Number of significant PCA dimensions used as covariates for regression analysis.
--correlationthreshold  | -c | float | Default 0.0. Only those k-mers whose correlation with phenotype is greater than _this value_ are retained for  regression analysis
--output | -o | output.txt | Mandatory. The path to output file to store the association values corresponding to nlr contigs in the given assembly.


