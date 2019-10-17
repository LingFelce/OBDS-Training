# Network Analysis 17th October 2019
============================================
## Biological Networks - basic components
  Nodes (circles): gene/protein/metabolite/ontology
  Edge: line connecting 2 nodes
  Activation/inhibition (direct/indirect)
## Biological Networks - basic features/network topology
  Degree: number of connections of a node
  Shortest paths or distance
  Hubs: high degree nodes
  Scale free networks: small number of neighbours, small number of hubs. Biological networks unlikely to be scale free?
  A single gene knocked out in a network is unlikely to lead to a lethal effect due to redundancy; however knocking out a hub
  gene could lead to lethality eg P53 in cancer is a hub.
  Clusters/subcommunities: group of nodes more internally connected than they are with rest of network
  Centralities: estimate of how important node/edge is for connectivity
  Size/density of network; motifs/cliques/clusters/sub-networks
## Types/sources of biological networks:
  DNA-protein (transcriptional regulatory networks, methylation)
  RNA-RNA (miRNA regulatory)
  RNA-Protein (splicing regulatory)
  Gene-Gene (co-expression)/Protein-Protein (coexpression, colocalisation, Gene Ontology etc)
 
 What's the point of building/analysing biological networks?
  Functional features, meta-analysis.

Ingenuity Pathway Analysis - not free, but they do a lot of literature search to back up pathways

## WGCNA: weighted gene correlation network analysis
  construct gene co-expression network
  identify modules
  relate modules to external information
  study module relationships
  find key drivers in interesting molecules
    
# Remember in R Studio to clear cache before Knit to HTML!
==========================================================
# Code below is from R script - some plots don't appear nicely in .Rmd
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials
# Tutorial 1, parts 1, 2b, 3, 4 (didn't actually download packages though) and 5

# Tutorial for WGCNA package for R - 17th October 2019

# Network analysis of liver expression data in female mice (microarray)

######################################################################

# Data input and cleaning

# load WGCNA package
library(WGCNA)

# The following setting is important, do not omit.

options(stringsAsFactors = FALSE)

#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv")

# Take a quick look at what is in the data set:
dim(femData)
# 135 samples, 3600 genes
names(femData)

# Remove auxiliary data and transpose expression data for further analysis
datExpr0 = as.data.frame(t(femData[, -c(1:8)]))
names(datExpr0) = femData$substanceBXH
rownames(datExpr0) = names(femData)[-c(1:8)]
dim(datExpr0)
# 135 rows and 3600 columns

# Check for genes and samples with too many missing values (if not already done QC)
# remove genes that aren't expressed in at least 50% of samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# if last statement is TRUE, all genes have passed filtering.
# if not, remove genes and samples from data
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

# cluster samples to see if any obvious outliers (not same as clustering genes - later)
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# looks like one sample is an outlier, F2_221 - set a height and cut branch at that height (y-axis is height)
# Plot a line to show the cut
abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# 1 sample above line (cluster 0), other 134 samples below line (cluster 1)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ] #keep rows
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)
# 134 rows (samples), 3600 columns (genes)
# datExpr contains expression data ready for network analysis

# read in clinical data and match to samples
traitData = read.csv("ClinicalTraits.csv")
dim(traitData)
# 361 rows and 38 columns
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)] # remove columns Note and comments
allTraits = allTraits[, c(2, 11:36) ]  #keeping Mice column and columns 11 - 36
dim(allTraits)
# 361 rows and 27 columns
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr) #rownames are sample names
traitRows = match(femaleSamples, allTraits$Mice) #match samples names with Mice names in traits file
datTraits = allTraits[traitRows, -1] #create new table, remove first column
rownames(datTraits) = allTraits[traitRows, 1] #set sample names as rownames

collectGarbage() #frees up memory space

# Visualise how clinical traits relate to sample clustering

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
# will show dendrogram on top, underneath is heatmap with traits as rows

# save expression and trait data
save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")

#################################
# Automatic network construction and module detection

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames

# Step by step construction of gene network and identification of modules

# Constructing weighted gene network needs choice of soft thresholding power B to which co-expression similarity is raised to calculate adjacency.
# Set threshold for co-expression  (used to calculate adjacency) - hard or soft threshold
# User chooses set of candidate powers, and function returns network indices that should be inspected

# Code below generates mean connectivity graph showing scale-free topology fit index curve
# Choose power which is lowest power where curve flattens out 

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
# can also specify networkType - unsigned (default), signed, signed hybrid
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power - connectivity of genes
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# choose a power where you have few nodes with connectivity rather than lots of nodes with lots of connectivity - more informative
# dependent on assuming scale free topology

# Calculate adjacencies, using power 6 - therefore weighted
softPower = 6
adjacency = adjacency(datExpr, power = softPower)

# Minimise effect of noise, Turn adjacency into topological overlap matrix, calculate corresponding dissimilarity
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Use hierarchical clustering to produce dendrogram of genes
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
# each branch corresponds to a gene. Branches that are gathered close together are densely connected, highly co-expressed genes
# module identification = identification of individual branches ("branch cutting")

# Branch cutting
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut. Default method is hybrid
# deepSplit is sensitivity of cluster splitting. If higher (max 4 for hybrid) then more and smaller clusters produced.
# distM is distance matrix, using dissimilarity
# if pamRespectsDendro = TRUE, then objects and small clusters will only be assigned to clusters that belong to same branch that originally assigned to
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
# table returned 22 modules (0 = unassigned genes) module 1 has largest number of genes, module 22 smallest number.
# plot module assignment under gene dendrogram
# Convert numeric labels into colors - assign different colour to each module
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Merging modules with similar expression profiles
# quantify co-expression similarity of entire modules, calculate eigengenes (1st principal component) and cluster on correlation
# Calculate eigengenes - use branch clustering from above
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Choose cut heigh of 0.25, which corresponds to correlation of 0.75 to merge similar modules
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
# cutHeight = maximum dissimilarity that qualifies modules for merging
# nodes that are below red line are merged with corresponding node above red line (4 modules merged)
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# Plot dendrogram again, with original and merged colours underneath - 23 modules to 19 modules
sizeGrWindow(12, 9) # not needed if using R Studio!
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")

#################################

# Relating modules to external information and identifying important genes

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-stepByStep.RData")
lnames

# Quantifying module trait associations
# Identify modules associated with measured clinical traits.
# As already have summary profile (eigengene) for each module, correlate eigengenes with external traits and look for significant associations
# Define numbers of genes and samples
nGenes = ncol(datExpr) # genes are columns
nSamples = nrow(datExpr) # samples are rows
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) # reorder given eigen vectors so that similar ones (measured by correlation) are next to each other - reorders columns
moduleTraitCor = cor(MEs, datTraits, use = "p") # correlation between eigengenes and traits
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) # calculate p value for correlation (from above) and samples 

# Colour code each association by correlation value
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
# scale bar shows correlation - red +ve correlation/association, white 0 no correlation, green -ve correlation/assocation
# number shows correlation value, number in brackets is p value
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
# Certain traits are strongly correlated with modules (red or green) - concentrate on weight

# Gene relationship to trait and important modules - gene significance and module membership
# Quantify association of individual genes with trait (weight) by defining Gene Significance as correlation between gene and trait
# For each module, define quantitative Module Membership as correlation of module eigengene and gene expression profile
# Allows quantification of gene similarity on array to every module

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")) # correlation of gene expression and module eigengene
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) # p value of module membership and samples

names(geneModuleMembership) = paste("MM", modNames, sep="") # start each column name with MM
names(MMPvalue) = paste("p.MM", modNames, sep="") # start each column name with p.MM

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p")) # correlation between weight and gene expression
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)) # p value of gene significance and samples

names(geneTraitSignificance) = paste("GS.", names(weight), sep="") # start each column with GS.
names(GSPvalue) = paste("p.GS.", names(weight), sep="") # start each column with p.GS.

# Intramodular analysis - identifying genes with high GS and MM
# pick magenta as has highest positive correlation with weight (different from tutorial)
module = "magenta"
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# genes highly associated with trait are important central elements of module associated with trait

# Summary output of network analysis results
# Have found modules associated with certain traits. Merge information with gene annotation (currently probe ID names)
names(datExpr) #probe names (column headings)
names(datExpr)[moduleColors=="magenta"] # probe IDs specific to magenta module

# Probe annotation file from gene expression array manufacturer (Entrez IDs)
annot = read.csv(file = "GeneAnnotation.csv")
dim(annot) # 23388 rows 34 columns
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH) #substanceBXH is probe name
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create data frame for all probes with info:
# probe ID, gene symbol, Locus Link ID (Entrez code), module colour,
# gene significance for weight, module membership and p values in all modules
# modules ordered by significance for weight
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder, ]

# write gene to file
write.csv(geneInfo, file = "geneInfo.csv") # information about all probes but only related to weight trait

#################################
# Interfacing network analysis with other data such as functional annotation and gene ontology

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-stepByStep.RData")
lnames

# Output gene lists from interesting modules to analyse externally
# Read in the probe annotation
annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file 
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
intModules = c("magenta", "purple", "green")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

# Uses separate Bioconductor packages GO.db, AnnotationDBI and org.Mm.eg.db (for mouse)
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10)
tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)
keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab

################################
# Network visualisation

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6)
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
# heatmap can depict adjacencies or topological overlaps, light colours show low overlap, darker colours show overlap
# gene dendrograms and module colours also plotted along heatmap
# can only be used if network calculated using single block approach (step by step)
# generating heatmap takes a long time, can restrict number of genes to speed up plotting

nSelect = 400 # 400 genes - how selected? Top genes?
# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

# Visualising network of eigengenes - look at relationship between modules
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g)
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# replot dendrogram and heatmap so separate (not on top of each other)
# Plot the dendrogram
sizeGrWindow(6,6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

# Eigengene dendrogram and heatmap identify groups of correlated eigengenes - meta-modules
# modules that are next to each other are highly related - mutal correlations are stronger than correlations with weight
# meta-modules are defined as tight clusters of modules (correlation of eigengenes of at least 0.5)
