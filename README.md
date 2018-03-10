# Medicago Truncatula Verticillium WGCNA Analysis

This software uses performs WGCNA analysis on Medicago Truncatula sequence data.  Since WGCNA analysis is computationally intensive the output files are provided, so there isn't a need to run this software on a cluster or on the cloud.  An overview for each file and the sample data is given, followed by some project notes.



## R Files

### LoadAndCluster.R

Loading, screening, and pre-clustering the data is this file's responsibilty.  After screening out the low quality RNA reads, the file prepares the data for processing.

### NetworkConstruction.R

The screened and preprocessed data is fed into the main WGNCA algorithm.  This algorithm groups RNA segments into clusters, each cluster is charactarized by a common expression profile across differing test conditions and replicant trials.

### CompareClinicalTraits.R

Each cluster is matched with a specific clinical trait and treatment.  Identifying the clusters whose genes turn on and off in response to the viral infection is the goal of the analysis.

## R Data

### PlantTrait-01-datainput.RData

Output file containing the screened and preprocessed data.

### Plant-02-networkConstruction-auto.Rdata

Output file containing the associative gene networks that cluster the individual RNA sequences into groups.  Each group has a common expression profile across the different clinical conditions.

### analysisDifferentialExpression.Rproj

Output file that graphically explores the network data.

## Graphical Output

### HeatmapDifferentialExpression.pdf

A strong red or strong green color indicates a strong correlation between a cluster's average expression profile and a speicific clinical treament.  Clinical treatments is a formal term that might mean comparing the fungal immune and the fungal susceptible plants different responses to infection.  The differentially expressed genes associated with this clinical trait might be associated with fungal immunity.  This heatmap help identifies which cluster is associated with what clinical trait.

## Input Files
### MtrEILNormalized.txt

The normalized RNA counts.

### differentiallyExpressedContrast.xlsx

The Clinical contrasts used for cluster association.


### Summary: 
This module implements the gibbs sampler algorithm, which is used to find common motifs in DNA sequences.  This probablistic search algorithm runs in polynomial time, which is an improvement over the brute force algorithm's exponential running time.
          
          
## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

R Statistical Software Suite

Ability to send the job to a computing cluster or the cloud.  Eight 4-GB nodes are recommended for this job.  If this isn't possible, the output files are included.

### Installing

A step by step series of examples that tell you have to get a development env running

Download and save all files to a local directory.

Run R-package.


## Running the tests

This project was tested function by function as part of a CourseEra project using their web testing service.

The algorithm locates and prints each binding site.  You can compare the algorithm's answer to the real answer, by opening the solution file.  In that file the real motif binding site appears in capital letters.


## Versioning

v 0.1

## Authors

* **Scott Czopek** - *Initial work* - 1/22/17 - [WGCNA Analysis](https://github.com/sczopek/Medicago-Truncatula-Verticilium-WGNCA-Analysis)

## License

This project is free to copy and distribute.

## Acknowledgments

* I would like to thank my advisor Sergey Nuzdhin for giving me the chance to learn about genetics.
* Thanks to CourseEra I know more about motif finding.
