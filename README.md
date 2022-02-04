# Code repository for the CAGI6 SickKids challenge

A challenge submission by the [DROPpers](https://github.com/gagneurlab/drop).

The team members are: [Julien Gagneur](https://github.com/gagneur), [Christian Mertes](https://github.com/c-mertes), [Ines Scheller](https://github.com/ischeller), [Nicholas H. Smith](https://github.com/nickhsmith), [Vicente A. Yépez](vyepez88).

This is the code repository used to generate the results for our two submissions for the [CAGI6 SickKids challenge](https://genomeinterpretation.org/cagi6-sickkids.html). We tackled the challenge by predicting the molecular events underlying disease from a patient's genome and transcriptome using variant annotation, aberrant gene expression events, and human phenotype ontology.

The code consists of 4 parts that are described below:

1.  Aberrant event detection in RNA-seq data using [DROP](https://github.com/gagneurlab/drop).
2.  Annotating and filtering variants
3.  Computing phenotypic similarity scores
4.  Prioritizing events using XGBoost

A detailed description of our full analysis can be found [here](docs/Methods.pdf).

## Aberrant event detection in RNA-seq

We used [DROP](https://github.com/gagneurlab/drop) with the default configuration to call aberrant events. To run the full pipeline, we suggest in a nutshell (i) to install DROP through [bioconda](https://anaconda.org/bioconda/drop), (ii) put all relevant data into `Data/project_data/raw/`, and (iii) create a sample annotation in `Data/project_data/sample_annotation.tsv`. You can then run the full DROP pipeline with

    snakemake -j 20

The main pipeline configuration can be found [here](/config.yaml).

## Variant annotation and filtering

As described in the method, we used [VEP](https://m.ensembl.org/info/docs/tools/vep/index.html) to annotate the variants. In short, we annotated all default information from VEP, allele frequencies through gnomAD, added CADD, SpliceAI, and EVE scores, as well as ClinVar and UTRannotator information. The respective configuration and scripts can be found [here](/vep_GRCh37.config) and [here](/Snakefile_vep_anno.smk). After adapting the config to your local infrastructure and a successful run of the DROP pipeline, you should be able to run it with snakemake as following:

    snakemake -j 20 --snakefile Snakefile_vep_anno.smk

## Phenotypic similarity scores

We computed the phenotypic similarity scores as described by [Kopajtich et al](https://www.medrxiv.org/content/10.1101/2021.03.09.21253187v2). A more detailed version can be found also in our [Methods section](/docs/Methods.pdf). The scripts to run it can be found [here](/src/r/hpo).

## Prioritizing events using XGBoost

For the final submission of the SickKids challenge, we used XGBoost to predict the disease-causing gene given the HPO terms, genetic information, as well as RNA-seq-based aberrant events of an individual. The code for our model can be found [here](/src/r/run_xgboost_causal_variants.R) and [here](/src/r/xgsboost_cv.R). The model can be trained as soon as the RNA-seq outliers are called, the variants are annotated, filtered, and preprocessed, and the phenotypic similarity scores are calculated.

## Disclaimer

This code was put together for the CAGI6 SickKids challenge and is not production-ready. This repository is meant to be complementary to our method description and to help others to get started. If there is any question about the model/code please create a new [issue](https://github.com/gagneurlab/cagi6_sickkids/issues).
