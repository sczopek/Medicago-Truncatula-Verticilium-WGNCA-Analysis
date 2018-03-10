# Medicago Truncatula Verticillium WGCNA Analysis

This software uses performs WGCNA analysis on Medicago Truncatula sequence data.  WGCNA (Weighted Genome Coexpression Network Analysis) divides the RNA sequences into clusters that share a common count expression profile.  This cluster's count expression profile is vector, and this vector's coordinates are the average RNA counts for each experiment in the trial.  (Each cluster is characterized by a unique expression profile.)  The goal is to identify the clusters associated with the plant's fungal immune response providing a list of immunity candidate genes.

Since WGCNA analysis is computationally intensive the output files are provided, so there isn't a need to run this software on a cluster or on the cloud.  An overview for each file and the sample data is given, followed by some project notes.



## R Files

### LoadAndCluster.R

Loading, screening, and pre-clustering the data is this file's responsibilty.  After screening out the low quality RNA reads, the file prepares the data for processing.

### NetworkConstruction.R

The screened and preprocessed data is fed into the main WGNCA algorithm.  This algorithm groups RNA segments into clusters, each cluster is charactarized by a common expression profile across differing test conditions and replicant trials.

### CompareClinicalTraits.R

Each cluster is matched with a specific clinical treatment.  Identifying the clusters whose genes turn on and off in response to the viral infection is the goal of the analysis.

## R Data

### PlantTrait-01-datainput.RData

Output file containing the screened and preprocessed data.

### Plant-02-networkConstruction-auto.Rdata

Output file containing the associative gene networks that cluster the individual RNA sequences into groups.  Each group has a common expression profile across the different clinical conditions.

### analysisDifferentialExpression.Rproj

Output file that graphically explores the network data.

## Graphical Output

### HeatmapDifferentialExpression.pdf

A strong red or strong green color indicates a strong correlation between a cluster's average expression profile and a speicific clinical contrast.  Clinical contrast is a formal term that might mean comparing the different responses of fungal immune and the fungal susceptible plants after being infected with the Verticillium fungus.  The differentially expressed genes associated with this clinical contrast may be associated with fungal immunity.  This heatmap help identifies which cluster is associated with what clinical contrast.

## Input Files
### MtrEILNormalized.txt

The normalized RNA Illumina sequence count data.

### differentiallyExpressedContrast.xlsx

The Clinical contrasts used for cluster association.


          
## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.  All files can be run locally except for NetworkConstruction.R, which must be instead run on the cloud or a computing cluster.

### Prerequisites

R Statistical Software Suite

The output files are included since a computing cluster or the cloud might not be accessable.  Six 8-GB nodes are recommended for this job, should you choose to run this on the cloud.

### Installing

A step by step series of examples that tell you have to get a development env running

Download and save all files to a local directory.

In R run Run LoadAndCluster.R and CompareClinicalTraits.R.  If these programs are run in R-studio line enabling the multi-threading capability might throw an error preventing execution.

Have fun exploring the graphs.


## Versioning

v 0.2

## Authors

* **Scott Czopek** - *Initial work* - 1/22/17 - [WGCNA Analysis](https://github.com/sczopek/Medicago-Truncatula-Verticilium-WGNCA-Analysis)

## License

This project is free to copy and distribute.

## Acknowledgments

* I would like to thank my advisor Sergey Nuzdhin for giving me the chance to learn about genetics and assigning me this project.
* I would like to thank NIH for the funding to work on this project.
