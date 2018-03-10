##############
##############
##############

# Load and Cluster Normalized Sequence Data

##############
##############
##############

##  Modified by Scott Czopek
##  Originally written by Steven Horvath


# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# uncomment when run on parallel processor
# function not supported in RStudio
#enableWGCNAThreads()


#Read in the plant immune response data set
plantImmuneResponseData = read.table("MtrEILNormalized.txt");
plantImmuneResponseData = t(plantImmuneResponseData);
#Select all genes for differential expression WGCNA analysis
plantImmuneResponseData=plantImmuneResponseData[c(1:28),]
datExpr0 = as.data.frame(plantImmuneResponseData);



gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


sampleTree = flashClust(dist(datExpr0), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)






# plantData = datExpr
# geneNames<-names(plantData)
# for (name in geneNames){
#   plantData[[name]]<-as.numeric(plantData[[name]])
# }
# datExpr = plantData

# Plot a line to show the cut
abline(h = 50000, col = "red");
# Determine cluster under the line
# to keep all the samples changed cutheight=15 to cutheight=50000
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)






library(gdata)
traitData=read.xls("differentiallyExpressedContrast.xlsx")
rownames(traitData)=traitData[,1]
traitData=traitData[,-1]


names(traitData)
traitNames<-names(traitData)
for (name in traitNames){
  traitData[[name]]<-as.numeric(traitData[[name]])
}
traitData[is.na(traitData)] <- 0
datTraits = traitData
collectGarbage();




# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


save(datExpr, datTraits, file = "PlantTrait-01-dataInput.RData")