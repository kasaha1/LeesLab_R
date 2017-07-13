## To install Packages-------------
instPak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
  install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#------------- Packages ----
packages <- c("ggplot2", "dplyr", "data.table", "DGCA","matrixStats")
instPak (packages) 
#-----------------------------


# data load
data(darmanis)
data(design_mat)
darmanis
design_mat

# n_oligo_samples = 38; n_neuron_samples = 120 
# cell_type = c(rep("oligodendrocyte", n_oligo_samples), rep("neuron", n_neuron_samples))
# design_mat = makeDesign(cell_type) # make design_mat

# filtering
darmanis_mean_filtered = filterGenes(darmanis,filterTypes = c("central","dispersion"), filterCentralType = "median", filterCentralPercentile = 0.3,filterDispersionType = "cv",filterDispersionPercentile = 0.3) # low average removal, ;low SD/average removal


# general
ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat,
                     compare = c("oligodendrocyte", "neuron"),
                     adjust = "perm", nPerm = 10)

# get cor
cor_res = getCors(inputMat = darmanis, design = design_mat)  # get cor
dcPairs_res = pairwiseDCor(cor_res, compare = c("oligodendrocyte", "neuron")) # get different cor

# plotting the expression of gene pairs in multiple conditions

# darmanis = darmanis[ , -which.max(darmanis["COX6A1", ])] #remove one outlier sample before visualization 
# design_mat = design_mat[-which.max(darmanis["COX6A1", ]), ]
plotCors(inputMat = darmanis, design = design_mat, compare = c("oligodendrocyte", "neuron"), geneA = "RTN4", geneB = "COX6A1")

# Adjusting the resulting p-values without permutation testing
dd_pairs = dcTopPairs(dcPairs_res, nPairs = 100, classify = TRUE, adjust = "BH")
dd_pairs = dcTopPairs(dcPairs_res, nPairs = 100, classify = TRUE, adjust = "fndr",plotFdr = F)


#permutation

ddcor_res_perm = ddcorAll(inputMat = darmanis, design = design_mat,compare = c("oligodendrocyte", "neuron"),adjust = "perm", heatmapPlot = FALSE, nPerm = 10, splitSet = "RTN4")

# Imputing gene expression measurements NA -->  k-nearest neighbors
library(impute, quietly = TRUE)
darmanis_na = darmanis
darmanis_na["RTN4", 1] = NA #add an NA value to demonstrate
ddcor_res = ddcorAll(inputMat = darmanis_na, design = design_mat,
                     compare = c("oligodendrocyte", "neuron"),
                     adjust = "none", nPerm = 0, impute = TRUE)

# Gene-trait differential correlation analysis
data(ages_darmanis)
rownames_darmanis = rownames(darmanis)
darmanis_with_traits = rbind(ages_darmanis, darmanis)
rownames(darmanis_with_traits) = c("age", rownames_darmanis)
ddcor_res = ddcorAll(inputMat = darmanis_with_traits, design = design_mat,
                     compare = c("oligodendrocyte", "neuron"),
                     adjust = "none", nPerm = 0, splitSet = "age")


# heatmap
library(gplots, quietly = TRUE)
darmanis_top =  filterGenes(darmanis, 
                            filterTypes = c("central", "dispersion"), filterCentralPercentile = 0.75, 
                            filterDispersionPercentile = 0.75)
ddcor_res = ddcorAll(inputMat = darmanis_top, design = design_mat,
                     compare = c("oligodendrocyte", "neuron"),
                     adjust = "none", heatmapPlot = TRUE, nPerm = 0, nPairs = "all")

# Calculating the average change in correlation for each gene with all others
library(matrixStats)
ddcor_res = ddcorAll(inputMat = darmanis_top, design = design_mat,
                     compare = c("oligodendrocyte", "neuron"),
                     adjust = "perm", heatmapPlot = FALSE, nPerm = 20, nPairs = "all",
                     getDCorAvg = TRUE, dCorAvgType = "gene_average", dCorAvgMethod = "median")
head(ddcor_res$avg_dcor)

# Calculating the average change in correlation for all genes with all genes
library(matrixStats)
ddcor_res = ddcorAll(inputMat = darmanis_top, design = design_mat,
                     compare = c("oligodendrocyte", "neuron"),
                     adjust = "perm", heatmapPlot = FALSE, nPerm = 20, nPairs = "all",
                     getDCorAvg = TRUE, dCorAvgType = "total_average", dCorAvgMethod = "median")
head(ddcor_res$avg_dcor)



# GENE ONTOLOGY
library(GOstats, quietly = TRUE)
library(HGNChelper, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
ddcor_res_APP = ddcorAll(inputMat = darmanis, design = design_mat,
                         compare = c("oligodendrocyte", "neuron"),
                         adjust = "none", heatmapPlot = FALSE, nPerm = 0, splitSet = "APP")
ddcorGO_res = ddcorGO(ddcor_res_APP, universe = rownames(darmanis), 
                      gene_ontology = "all", HGNC_clean = TRUE, HGNC_switch = TRUE, annotation = "org.Hs.eg.db", calculateVariance = TRUE)
str(ddcorGO_res)
k <- ddcorGO_res$enrichment_significant_gain_of_correlation_genes$BP 



# MEGENA
module_genes = list(
  mod1 = rownames(darmanis)[1:100],
  mod2 = rownames(darmanis)[90:190],
  mod3 = rownames(darmanis)[190:290],
  mod4 = rownames(darmanis)[330:340],
  mod5 = rownames(darmanis)[350:360],
  mod6 = rownames(darmanis)[400:405])
modules = stack(module_genes)
modules$ind = as.character(modules$ind)
str(modules)
# Module-based differential correlation
moduleDC_res = moduleDC(inputMat = darmanis, design = design_mat,
                        compare = c("oligodendrocyte", "neuron"), genes = modules$values,
                        labels = modules$ind, nPerm = 50, number_DC_genes = 3,
                        dCorAvgMethod = "median")
head(moduleDC_res)

mod1_genes = modules[modules$ind == "mod1", "values"]
darmanis_mod1 = darmanis[mod1_genes, ]
moduleDC_res = ddcorAll(inputMat = darmanis_mod1, design = design_mat,
                        compare = c("oligodendrocyte", "neuron"), nPerm = 50,
                        getDCorAvg = TRUE, dCorAvgType = "gene_average",
                        dCorAvgMethod = "median")
head(moduleDC_res[["avg_dcor"]])
# Module-based gene ontology (GO) enrichment
library(GOstats, quietly = TRUE)
library(HGNChelper, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
moduleGO_res = moduleGO(genes = modules$values, labels = modules$ind,
                        universe = rownames(darmanis), pval_GO_cutoff = 1)

moduleGO_df = extractModuleGO(moduleGO_res) # extract results
library(ggplot2, quietly = TRUE)
plotModuleGO(moduleGO_df, nTerms = 4, text_size = 8, coord_flip = TRUE) # plot

# Differential correlation module detection through integration with MEGENA
library(MEGENA, quietly = TRUE)
ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat,
                     compare = c("oligodendrocyte", "neuron"),
                     adjust = "none", heatmapPlot = FALSE, nPerm = 0, nPairs = "all")
megena_res = ddMEGENA(ddcor_res, adjusted = FALSE, evalCompactness = FALSE)
