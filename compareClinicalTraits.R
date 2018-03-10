##############
##############
##############

# Clinical Trait Module Comparison

##############
##############
##############

##  Modified by Scott Czopek
##  Originally written by Steven Horvath




# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "PlantTrait-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Plant-02-networkConstruction-auto.RData");
lnames


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



# Early Vaa noVaa Comparison Resistant
# Define variable VaaE containing the VaaMiddle column of datTrait
VaaE = as.data.frame(datTraits$VaaMiddle);
names(VaaE) = "VaaMiddle";
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, VaaE, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(VaaE), sep="");
names(GSPvalue) = paste("p.GS.", names(VaaE), sep="");



module = "lightgreen"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for VaaE",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)




############################################
############################################
############################################
###    Write list for AgriGO SEA analysis
###    This is a list of genes
############################################

# (Export Genes By Module)
container = matrix(,length(MEs),length(datExpr))
for (modColor in 1:length(MEs))
{
    # create storage container for module's gene names here
    # or outside of loop
    writeIndex = 1
    names(container)[modColor] = names(MEs)[modColor]
    for(gene in 1:length(moduleColors))
    {
      
      if(names(MEs)[modColor]    == paste("ME", moduleColors[gene], sep = ""))
      {
        #reduce isoform name to gene name
        geneName = gsub("\\..*","",names(datExpr[gene]))
        # put gene from module in storage container
        container[modColor,writeIndex] = geneName
        writeIndex=writeIndex+1
      }
    }
    #
}

dir.create(file.path(getwd(), "ModuleGeneList"), showWarnings = FALSE)
for (modColor in 1:length(MEs))
{
  write.table(container[modColor,], 
              paste("ModuleGeneList/",names(container[modColor]),"Genes_DiffExpr.txt", sep = ""), 
              sep="\n", row.names=FALSE, col.names=FALSE,
              quote=FALSE,
              na = "") 
  
}








