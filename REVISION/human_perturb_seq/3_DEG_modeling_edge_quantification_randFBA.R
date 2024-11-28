# edge quantification analysis of the interactions between core functions 

# the four scripts for edge quantification are generally similar, with alterations for specific analysis. 
# we commented one (3_DEG_modeling_edge_direction_test_randFBA.R) in greater details and the remaining fours 
# were concisely commented with more focus on places that are different than 3_DEG_modeling_edge_direction_test_randFBA.R

# this script focus on the edge quantification with testing hypothesis of overrepresentation of a specific interacting edge
# among all edges of a core function perturbation, using the randomization of core function associations


library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(Seurat)
library(SeuratDisk)
library(pheatmap)

# load data 
# load the met gene classifications
classMat = read.csv('output/FBA_classification_matrix.csv',row.names = 1)
colnames(classMat)[colnames(classMat) == 'ECM'] = 'pro_modi' # for the reuse of old code, just use the old name 

# load DEG result
adDE_res = read.csv('input/source_data/processed_data/anderson-darling_p-values_BH-corrected.csv', row.names = 1)
# although not noted in the figshare page where the data is downloaded, the column number (# perturbations) is equal to 
# that of the K562_gwps (K562_gwps_normalized_bulk_01.h5seurat), indicating this is the result of the genome-wide perturb-seq dataset in K562 cells 

# we will use the DE result and skip the thresholding on pseudobulk data (because it is z-score anyways - not fold change)
# (we can use raw UMI count (like Count-per-Ten-Thousound)) but that is more noisy and has batch problems)
# Load the .h5ad file into a Seurat object for meta information and sign of changes 
seurat_object <- LoadH5Seurat("./input/source_data/processed_data/K562_gwps_normalized_bulk_01.h5seurat")

# Assuming your Seurat object is named seurat_object
# get the normalzied expression matrix
exp_mat <- seurat_object@assays$RNA@counts # same as scale.data so it is scaled data
colnames(exp_mat) = paste('X',colnames(exp_mat),sep = '')

# get the metadata for each perturbations 
cell_annotations <- seurat_object@meta.data
head(cell_annotations)  # View the first few rows of the cell annotations

# get metadata for each genes for future filtering 
gene_annotations <- seurat_object@assays$RNA@meta.features
head(gene_annotations)  # View the first few rows of the cell annotations

# subset metabolic genes and metabolic perturbations and visualize 
metGenes = read.csv('input/human_model/genes.tsv', sep = '\t', header = T)

# exclude polymerases for now 
pols = read.csv('output/DNA_RNA_pol.csv')
metGenes = metGenes[!(metGenes$genes %in% pols$pols), ]

# get the met perturbations
perturbed_gene = colnames(adDE_res)
perturbed_gene = unlist(lapply(strsplit(perturbed_gene, '_'),function(x){x[[4]]}))
# subset matrix 
padj_mat_met = adDE_res[rownames(adDE_res) %in% metGenes$genes, perturbed_gene %in% metGenes$genes]
# the pass-filter matrix is around 1000 x 1700

# focus on responsive ones (this threshold of 10 may be removed and can be 8 or arbituary - the authors didnt explicitly define a hard cutoff)
cell_annotations <- seurat_object@meta.data
cell_annotations$perturbed_gene = perturbed_gene[match(rownames(cell_annotations), str_remove(colnames(adDE_res),'^X'))]
cell_annotations$perturbation_id = paste('X',rownames(cell_annotations),sep = '')
# padj_mat_met_resp = padj_mat_met[, colnames(padj_mat_met) %in% cell_annotations$perturbation_id[cell_annotations$anderson_darling_counts > 10]]
padj_mat_met_resp = padj_mat_met[, colnames(padj_mat_met) %in% colnames(exp_mat)]

# the final matrix is 975 by 472

# the z-score matrix 
exp_mat_met_resp = as.matrix(exp_mat[rownames(padj_mat_met_resp), colnames(padj_mat_met_resp)])

# mask change that is the RNAi targeted gene 
for (i in 1:ncol(padj_mat_met_resp)){
  if (any(perturbed_gene[colnames(padj_mat_met_resp)[i] == colnames(adDE_res)] == rownames(padj_mat_met_resp))){
    padj_mat_met_resp[perturbed_gene[colnames(padj_mat_met_resp)[i] == colnames(adDE_res)],i] = 1 # mask to insignificant (1)
  }
}


# set cutoff
score_cutoff = 0.01 # because of the sensitivity issue in scRNA-seq, it seems we should use 


# visualize the classification of all genes that are targeted in human 1 responsive (at least two up or two down) perturbations
ct1 = colSums(padj_mat_met_resp < score_cutoff & exp_mat_met_resp > 0,na.rm = T) # start with padj = 0.05
ct2 = colSums(padj_mat_met_resp < score_cutoff & exp_mat_met_resp < 0,na.rm = T)
human1_resp = union(names(ct1)[ct1 > 1], names(ct2)[ct2 > 1])
tmp = classMat[rownames(classMat) %in% cell_annotations[match(human1_resp, cell_annotations$perturbation_id),'perturbed_gene'],]
colSums(tmp)


# PART I: VISUALIZE DEG IN EACH PERTURBATION
# We first focus on perturbations whose targeted gene is associated with unique core function.

# NOTE: to keep this modeling simple, the class selection by DE similarity is not performed in this step; the rewiring 
# model calculation is fully based on FBA classification and the multi-class genes were used literally (assuming it affects multiple objectives)
# and the unclassified conditions will be left out from the analysis. This also guarantees the same set of RNAi was analyzed in every randomization 

# we get the conditions to analyze (icel_responsive (at least 2 up or 2 down DEG) and classified)
total_condition_analyzed = c()
for (i in 1:length(human1_resp)){
  if (rowSums(classMat[cell_annotations$perturbed_gene[match(human1_resp[i], cell_annotations$perturbation_id)],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
    total_condition_analyzed = c(total_condition_analyzed, human1_resp[i])
  }
}



# quantify the edges 

# 01072024: this simple model was kept as a control point to compare with old results; but in publication,
# we used the FBA-fused model all the time

model = list()
modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
for (obj in modelObj){
  model[[obj]] = list()
}
model[['energy']] = c('energy','lipid','pro_modi','pro_syn','nucl_acid') 
model[['lipid']] = c('lipid')
model[['pro_modi']] = c('pro_modi') 
model[['pro_syn']] = 'pro_syn'
model[['nucl_acid']] = 'nucl_acid'

model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
model_explained$UP_yes = 0
model_explained$UP_no = 0
model_explained$DOWN_yes = 0
model_explained$DOWN_no = 0

model_explained_single = model_explained

# add the fitted model - simple total count of DEG
# fitted model considers all possible interactions
model_fit = list()
for (obj in modelObj){
  model_fit[[obj]] = list()
  model_fit[[obj]][['up']] = list('energy' = 0,
                                  'lipid' = 0,
                                  'pro_modi' = 0,
                                  'pro_syn' = 0,
                                  'nucl_acid' = 0)
  model_fit[[obj]][['down']] = list('energy' = 0,
                                  'lipid' = 0,
                                  'pro_modi' = 0,
                                  'pro_syn' = 0,
                                  'nucl_acid' = 0)
}
# weighted total count of DEG - we normalize the total DEG count of each condition to 1 (to avoid confounding by conditions will super high DE number)
# (when we calculate, we DO NOT seperate the up and down genes)
model_fit_wtd = list()
for (obj in modelObj){
  model_fit_wtd[[obj]] = list()
  model_fit_wtd[[obj]][['up']] = list('energy' = 0,
                                  'lipid' = 0,
                                  'pro_modi' = 0,
                                  'pro_syn' = 0,
                                  'nucl_acid' = 0)
  model_fit_wtd[[obj]][['down']] = list('energy' = 0,
                                    'lipid' = 0,
                                    'pro_modi' = 0,
                                    'pro_syn' = 0,
                                    'nucl_acid' = 0)
}
model_fit_single = list('raw_counts' = model_fit,
                        'normalized_counts' = model_fit_wtd)
model_fit_all = model_fit_single
for (condInd in 1:length(total_condition_analyzed)){
    myCond = total_condition_analyzed[condInd]
    DEGs_up = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] > 0)]
    DEGs_down = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] < 0)]
    
    # determine affected obj
    labeledObj = colnames(classMat)[classMat[cell_annotations$perturbed_gene[cell_annotations$perturbation_id == myCond],]==1]
    labeledObj = intersect(labeledObj, modelObj)
    
    if (length(labeledObj) > 0){
      affectedObj = c()
      for (i in 1:length(labeledObj)){
        affectedObj = c(affectedObj, model[[labeledObj[i]]])
      }
      
      # compare with DEG 
      subClassMat_total = classMat[rownames(classMat) %in% c(DEGs_up,DEGs_down), modelObj]
      subClassMat_total_wtd = subClassMat_total / rowSums(subClassMat_total)
      
      # up genes
      subClassMat = classMat[rownames(classMat) %in% DEGs_up, modelObj]
      # calculate the DEG explained by model
      model_explained$UP_yes[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) > 0)
      model_explained$UP_no[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) == 0)
      
      # fill in the observed model 
      if (length(labeledObj) == 1){
        # only consider single-category genes
        subClassMat_single = subClassMat[rowSums(subClassMat) == 1,]
        
        # calculate the DEG explained by model
        subClassMat_single2 = subClassMat # mask all multiple genes
        subClassMat_single2[rowSums(subClassMat_single2[,1:5]) > 1,1:5] = 0
        model_explained_single$UP_yes[condInd] = sum(rowSums(subClassMat_single2[,affectedObj,drop = F]) > 0)
        model_explained_single$UP_no[condInd] = sum(rowSums(subClassMat_single2[,affectedObj,drop = F]) == 0)
        
        mycounts = colSums(subClassMat_single)
        for (ii in 1:length(mycounts)){
          model_fit_single$raw_counts[[labeledObj]][['up']][[names(mycounts)[ii]]] = 
            model_fit_single$raw_counts[[labeledObj]][['up']][[names(mycounts)[ii]]] + mycounts[ii]
        }
        if (any(rowSums(subClassMat_total) == 1)){
          mycounts_norm = mycounts / sum(subClassMat_total[rowSums(subClassMat_total) == 1,])
        }else{ # no singles; my counts will be all zeros
          mycounts_norm = mycounts
        }
        mycounts_norm[is.na(mycounts_norm)] = 0
        for (ii in 1:length(mycounts_norm)){
          model_fit_single$normalized_counts[[labeledObj]][['up']][[names(mycounts_norm)[ii]]] = 
            model_fit_single$normalized_counts[[labeledObj]][['up']][[names(mycounts_norm)[ii]]] + mycounts_norm[ii]
        }
      }else{
        # calculate the DEG explained by model
        model_explained_single$UP_yes[condInd] = NA
        model_explained_single$UP_no[condInd] = NA
      }
      # assign the literal model (all obj counted independently)
      # let's start with simple literal counting: multi-obj RNAi were counted multiple times and 
      # multi-obj DE was counted as independent DE; if it is too noisy, we can weight the multi cases 
      # by equally divide their contribution (multi-obj RNAi was counted as 0.5+0.5, etc)
      subClassMat_wtd = subClassMat / rowSums(subClassMat)
      mycounts = colSums(subClassMat_wtd,na.rm = T)
      for (jj in 1:length(labeledObj)){
        for (ii in 1:length(mycounts)){
          model_fit_all$raw_counts[[labeledObj[jj]]][['up']][[names(mycounts)[ii]]] = 
            model_fit_all$raw_counts[[labeledObj[jj]]][['up']][[names(mycounts)[ii]]] + (mycounts[ii] /length(labeledObj))
        }
        mycounts_norm = mycounts / sum(subClassMat_total_wtd,na.rm = T)
        mycounts_norm[is.na(mycounts_norm)] = 0
        for (ii in 1:length(mycounts_norm)){
          model_fit_all$normalized_counts[[labeledObj[jj]]][['up']][[names(mycounts_norm)[ii]]] = 
            model_fit_all$normalized_counts[[labeledObj[jj]]][['up']][[names(mycounts_norm)[ii]]] + (mycounts_norm[ii]/length(labeledObj))
        }
      }
      
      
      
      # down genes
      subClassMat = classMat[rownames(classMat) %in% DEGs_down, modelObj]
      # calculate the DEG explained by model
      model_explained$DOWN_yes[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) > 0)
      model_explained$DOWN_no[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) == 0)
    
      # fill in the observed model 
      # check if it is a single model 
      if (length(labeledObj) == 1){
        # only consider single-category genes
        subClassMat_single = subClassMat[rowSums(subClassMat) == 1,]
        
        # calculate the DEG explained by model
        subClassMat_single2 = subClassMat # mask all multiple genes
        subClassMat_single2[rowSums(subClassMat_single2[,1:5]) > 1,1:5] = 0
        model_explained_single$DOWN_yes[condInd] = sum(rowSums(subClassMat_single2[,setdiff(modelObj, affectedObj),drop = F]) > 0)
        model_explained_single$DOWN_no[condInd] = sum(rowSums(subClassMat_single2[,setdiff(modelObj, affectedObj),drop = F]) == 0)
        
        mycounts = colSums(subClassMat_single)
        for (ii in 1:length(mycounts)){
          model_fit_single$raw_counts[[labeledObj]][['down']][[names(mycounts)[ii]]] = 
            model_fit_single$raw_counts[[labeledObj]][['down']][[names(mycounts)[ii]]] + mycounts[ii]
        }
        if (any(rowSums(subClassMat_total) == 1)){
          mycounts_norm = mycounts / sum(subClassMat_total[rowSums(subClassMat_total) == 1,])
        }else{ # no singles; my counts will be all zeros
          mycounts_norm = mycounts
        }
        mycounts_norm[is.na(mycounts_norm)] = 0
        for (ii in 1:length(mycounts_norm)){
          model_fit_single$normalized_counts[[labeledObj]][['down']][[names(mycounts_norm)[ii]]] = 
            model_fit_single$normalized_counts[[labeledObj]][['down']][[names(mycounts_norm)[ii]]] + mycounts_norm[ii]
        }
      }else{
        model_explained_single$DOWN_yes[condInd] = NA
        model_explained_single$DOWN_no[condInd] = NA
        
      }
      
      # assign the literal model (all obj counted independently)
      # let's start with simple literal counting: multi-obj RNAi were counted multiple times and 
      # multi-obj DE was counted as independent DE; if it is too noisy, we can weight the multi cases 
      # by equally divide their contribution (multi-obj RNAi was counted as 0.5+0.5, etc)
      subClassMat_wtd = subClassMat / rowSums(subClassMat)
      mycounts = colSums(subClassMat_wtd,na.rm = T)
      for (jj in 1:length(labeledObj)){
        for (ii in 1:length(mycounts)){
          model_fit_all$raw_counts[[labeledObj[jj]]][['down']][[names(mycounts)[ii]]] = 
            model_fit_all$raw_counts[[labeledObj[jj]]][['down']][[names(mycounts)[ii]]] + (mycounts[ii]/length(labeledObj))
        }
        mycounts_norm = mycounts / sum(subClassMat_total_wtd, na.rm = T)
        mycounts_norm[is.na(mycounts_norm)] = 0
        for (ii in 1:length(mycounts_norm)){
          model_fit_all$normalized_counts[[labeledObj[jj]]][['down']][[names(mycounts_norm)[ii]]] = 
            model_fit_all$normalized_counts[[labeledObj[jj]]][['down']][[names(mycounts_norm)[ii]]] + (mycounts_norm[ii]/length(labeledObj))
        }
      }
    }else{
      stop('there is a bug!') # we only analyzed classified conditions
      # up genes
      subClassMat = classMat[rownames(classMat) %in% DEGs_up, modelObj]
      # calculate the DEG explained by model
      model_explained$UP_yes[condInd] = 0
      model_explained$UP_no[condInd] = nrow(subClassMat)
      
      # down genes
      subClassMat = classMat[rownames(classMat) %in% DEGs_down, modelObj]
      # calculate the DEG explained by model
      model_explained$DOWN_yes[condInd] = 0
      model_explained$DOWN_no[condInd] = nrow(subClassMat)
      
    }
}
# remove NA in singles 
model_explained_single = model_explained_single[!rowAlls(is.na(model_explained_single[,2:5])),]



# the overall average explained rate (up and down together) is 
tmp = model_explained[,2:5]
tmp = tmp / rowSums(tmp)
rewire_rate = rowSums(tmp[,c(1,3)])
rewire_rate[is.na(rewire_rate)] = 0
obs_rate = mean(rewire_rate)
obs_total_DE_rate = sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5])
# the overall average explained rate (up and down separately) is 
# up
tmp = model_explained[,2:3]
tmp = tmp / rowSums(tmp)
obs_rate_up = mean(tmp$UP_yes, na.rm = T)
obs_total_DE_rate_up = sum(model_explained[,c(2)]) / sum(model_explained[,2:3])
# down
tmp = model_explained[,4:5]
tmp = tmp / rowSums(tmp)
obs_rate_down = mean(tmp$DOWN_yes, na.rm = T)
obs_total_DE_rate_down = sum(model_explained[,c(4)]) / sum(model_explained[,4:5])


# the overall average explained rate (up and down together) is (for single)
tmp = model_explained_single[,2:5]
tmp = tmp / rowSums(tmp)
rewire_rate_single = rowSums(tmp[,c(1,3)])
rewire_rate_single[is.na(rewire_rate_single)] = 0
obs_rate_single = mean(rewire_rate_single)
obs_total_DE_rate_single = sum(model_explained_single[,c(2,4)]) / sum(model_explained_single[,2:5])
# the overall average explained rate (up and down seperately) is 
# up
tmp = model_explained_single[,2:3]
tmp = tmp / rowSums(tmp)
obs_rate_up_single = mean(tmp$UP_yes, na.rm = T)
obs_total_DE_rate_up_single = sum(model_explained_single[,c(2)]) / sum(model_explained_single[,2:3])
# down
tmp = model_explained_single[,4:5]
tmp = tmp / rowSums(tmp)
obs_rate_down_single = mean(tmp$DOWN_yes, na.rm = T)
obs_total_DE_rate_down_single = sum(model_explained_single[,c(4)]) / sum(model_explained_single[,4:5])


# quick visualization of the edges by proportion
library(igraph)
# the make_g_tbl function is different from that in 3_DEG_modeling_edge_direction_test_randFBA.R
make_g_tbl <- function(g){
  g_tbl = data.frame(RNAi = rep(modelObj,5), DEG = c(rep(modelObj[1], 5),
                                                     rep(modelObj[2], 5),
                                                     rep(modelObj[3], 5),
                                                     rep(modelObj[4], 5),
                                                     rep(modelObj[5], 5)))
  g_tbl$type = 'up'
  g_tbl$proportion = 0
  tmp = g_tbl
  tmp$type = 'down'
  g_tbl = rbind(g_tbl, tmp)
  for (i in 1:nrow(g_tbl)){
    g_tbl$proportion[i] = g[[g_tbl$RNAi[i]]][[g_tbl$type[i]]][[g_tbl$DEG[i]]]
  }
  # normalize within each type of RNAi (the proportion of a specific interaction over all interactions related to a core function perturbation)
  # the proportion calculated in this way will give a test on the hypothesis of overrepresentation of a specific edge among all possible edges instead of just the direction
  for (i in 1:length(modelObj)){
    g_tbl$proportion[g_tbl$RNAi == modelObj[i]] =  g_tbl$proportion[g_tbl$RNAi == modelObj[i]] / sum(g_tbl$proportion[g_tbl$RNAi == modelObj[i]])
  }
  return(g_tbl)
}

# visualization
g = model_fit_all$raw_counts
g_tbl = make_g_tbl(g)
# Creating the graph object
graph <- graph_from_data_frame(g_tbl[,1:2], directed = TRUE)
# Setting edge width and type
E(graph)$width <- g_tbl$proportion
E(graph)$type <- g_tbl$type
#E(graph)$arrow <- ifelse(g_tbl$type == 'up', 1, -1)
# Drawing the graph
plot(graph,
     edge.arrow.size = E(graph)$width / 10,
     edge.arrow.width = E(graph)$width / max(E(graph)$width),
     edge.width = E(graph)$width * 20,
     edge.arrow.mode = 1,
     edge.color = ifelse(E(graph)$type == 'up', 'red', 'blue'),
     vertex.color = "gold", 
     vertex.size = 20,
     vertex.frame.color = "gray", 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 1, 
     main = "Directed Graph")

# save the observed edge frequencies
g_tbl_single_raw = make_g_tbl(model_fit_single$raw_counts)
g_tbl_single_norm = make_g_tbl(model_fit_single$normalized_counts)
g_tbl_full_raw = make_g_tbl(model_fit_all$raw_counts)
g_tbl_full_norm = make_g_tbl(model_fit_all$normalized_counts)



# randomization to assess the significance 
nRand = 10000
set.seed(1126)
rand_stat = list()
rand_stat$default_CR_full = list(rand_rate = c(), rand_rate2 = c(), rand_rate_up = c(), rand_rate2_up = c(), rand_rate_down = c(), rand_rate2_down = c())
rand_stat$default_CR_single = list(rand_rate = c(), rand_rate2 = c(), rand_rate_up = c(), rand_rate2_up = c(), rand_rate_down = c(), rand_rate2_down = c())
rand_stat$edge_full = list(raw_count = g_tbl[,1:3],
                           norm_count = g_tbl[,1:3])
rand_stat$edge_single = list(raw_count = g_tbl[,1:3],
                             norm_count = g_tbl[,1:3])
  

# set up the randomization of gene classification

# we keep the number of labels for each obj and keep the analyzed conditions assigned to at least one label
mustHasClass = unique(cell_annotations$perturbed_gene[cell_annotations$perturbation_id %in% total_condition_analyzed])
justRandom = setdiff(rownames(classMat), mustHasClass)
hasClassInd = as.numeric(which(rowSums(classMat[,modelObj])>0))

for (nn in 1:nRand){
  # generate the random classification matrix
  new_gene_names = rep(NA, nrow(classMat))
  new_gene_names[sample(hasClassInd, length(mustHasClass))] = mustHasClass
  new_gene_names[is.na(new_gene_names)] = justRandom[sample(length(justRandom))]
  
  classMat_rand = classMat
  rownames(classMat_rand) = new_gene_names
  
  model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
  model_explained$UP_yes = 0
  model_explained$UP_no = 0
  model_explained$DOWN_yes = 0
  model_explained$DOWN_no = 0
  
  model_explained_single = model_explained
  
  # add the fitted model - simple total count of DEG
  # fitted model considers all possible interactions
  model_fit = list()
  for (obj in modelObj){
    model_fit[[obj]] = list()
    model_fit[[obj]][['up']] = list('energy' = 0,
                                    'lipid' = 0,
                                    'pro_modi' = 0,
                                    'pro_syn' = 0,
                                    'nucl_acid' = 0)
    model_fit[[obj]][['down']] = list('energy' = 0,
                                      'lipid' = 0,
                                      'pro_modi' = 0,
                                      'pro_syn' = 0,
                                      'nucl_acid' = 0)
  }
  # weighted total count of DEG - we normalize the total DEG count of each condition to 1 (to avoid confounding by conditions will super high DE number)
  # (when we calculate, we DO NOT seperate the up and down genes)
  model_fit_wtd = list()
  for (obj in modelObj){
    model_fit_wtd[[obj]] = list()
    model_fit_wtd[[obj]][['up']] = list('energy' = 0,
                                        'lipid' = 0,
                                        'pro_modi' = 0,
                                        'pro_syn' = 0,
                                        'nucl_acid' = 0)
    model_fit_wtd[[obj]][['down']] = list('energy' = 0,
                                          'lipid' = 0,
                                          'pro_modi' = 0,
                                          'pro_syn' = 0,
                                          'nucl_acid' = 0)
  }
  model_fit_single = list('raw_counts' = model_fit,
                          'normalized_counts' = model_fit_wtd)
  model_fit_all = model_fit_single
  for (condInd in 1:length(total_condition_analyzed)){
    myCond = total_condition_analyzed[condInd]
    DEGs_up = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] > 0)]
    DEGs_down = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] < 0)]
    
    # determine compromised obj
    labeledObj = colnames(classMat_rand)[classMat_rand[cell_annotations$perturbed_gene[cell_annotations$perturbation_id == myCond],]==1]
    labeledObj = intersect(labeledObj, modelObj)
    
    if (length(labeledObj) > 0){
      affectedObj = c()
      for (i in 1:length(labeledObj)){
        affectedObj = c(affectedObj, model[[labeledObj[i]]])
      }
      
      # compare with DEG 
      subClassMat_total = classMat_rand[rownames(classMat_rand) %in% c(DEGs_up,DEGs_down), modelObj]
      subClassMat_total_wtd = subClassMat_total / rowSums(subClassMat_total)
      
      # up genes
      subClassMat = classMat_rand[rownames(classMat_rand) %in% DEGs_up, modelObj]
      # calculate the DEG explained by model
      model_explained$UP_yes[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) > 0)
      model_explained$UP_no[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) == 0)
      
      # fill in the observed model 
      # check if it is a single model 
      if (length(labeledObj) == 1){
        # only consider single-category genes
        subClassMat_single = subClassMat[rowSums(subClassMat) == 1,]
        
        # calculate the DEG explained by model
        subClassMat_single2 = subClassMat # mask all multiple genes
        subClassMat_single2[rowSums(subClassMat_single2[,1:5]) > 1,1:5] = 0
        model_explained_single$UP_yes[condInd] = sum(rowSums(subClassMat_single2[,affectedObj,drop = F]) > 0)
        model_explained_single$UP_no[condInd] = sum(rowSums(subClassMat_single2[,affectedObj,drop = F]) == 0)
        
        mycounts = colSums(subClassMat_single)
        for (ii in 1:length(mycounts)){
          model_fit_single$raw_counts[[labeledObj]][['up']][[names(mycounts)[ii]]] = 
            model_fit_single$raw_counts[[labeledObj]][['up']][[names(mycounts)[ii]]] + mycounts[ii]
        }
        if (any(rowSums(subClassMat_total) == 1)){
          mycounts_norm = mycounts / sum(subClassMat_total[rowSums(subClassMat_total) == 1,])
        }else{ # no singles; my counts will be all zeros
          mycounts_norm = mycounts
        }
        mycounts_norm[is.na(mycounts_norm)] = 0
        for (ii in 1:length(mycounts_norm)){
          model_fit_single$normalized_counts[[labeledObj]][['up']][[names(mycounts_norm)[ii]]] = 
            model_fit_single$normalized_counts[[labeledObj]][['up']][[names(mycounts_norm)[ii]]] + mycounts_norm[ii]
        }
      }else{
        # calculate the DEG explained by model
        model_explained_single$UP_yes[condInd] = NA
        model_explained_single$UP_no[condInd] = NA
      }
      # assign the literal model (all obj counted independently)
      # let's start with simple literal counting: multi-obj RNAi were counted multiple times and 
      # multi-obj DE was counted as independent DE; if it is too noisy, we can weight the multi cases 
      # by equally divide their contribution (multi-obj RNAi was counted as 0.5+0.5, etc)
      subClassMat_wtd = subClassMat / rowSums(subClassMat)
      mycounts = colSums(subClassMat_wtd,na.rm = T)
      for (jj in 1:length(labeledObj)){
        for (ii in 1:length(mycounts)){
          model_fit_all$raw_counts[[labeledObj[jj]]][['up']][[names(mycounts)[ii]]] = 
            model_fit_all$raw_counts[[labeledObj[jj]]][['up']][[names(mycounts)[ii]]] + (mycounts[ii] /length(labeledObj))
        }
        mycounts_norm = mycounts / sum(subClassMat_total_wtd,na.rm = T)
        mycounts_norm[is.na(mycounts_norm)] = 0
        for (ii in 1:length(mycounts_norm)){
          model_fit_all$normalized_counts[[labeledObj[jj]]][['up']][[names(mycounts_norm)[ii]]] = 
            model_fit_all$normalized_counts[[labeledObj[jj]]][['up']][[names(mycounts_norm)[ii]]] + (mycounts_norm[ii]/length(labeledObj))
        }
      }
      
      
      
      # down genes
      subClassMat = classMat_rand[rownames(classMat_rand) %in% DEGs_down, modelObj]
      # calculate the DEG explained by model
      model_explained$DOWN_yes[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) > 0)
      model_explained$DOWN_no[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) == 0)
      
      # fill in the observed model 
      # check if it is a single model 
      if (length(labeledObj) == 1){
        # only consider single-category genes
        subClassMat_single = subClassMat[rowSums(subClassMat) == 1,]
        
        # calculate the DEG explained by model
        subClassMat_single2 = subClassMat # mask all multiple genes
        subClassMat_single2[rowSums(subClassMat_single2[,1:5]) > 1,1:5] = 0
        model_explained_single$DOWN_yes[condInd] = sum(rowSums(subClassMat_single2[,setdiff(modelObj, affectedObj),drop = F]) > 0)
        model_explained_single$DOWN_no[condInd] = sum(rowSums(subClassMat_single2[,setdiff(modelObj, affectedObj),drop = F]) == 0)
        
        mycounts = colSums(subClassMat_single)
        for (ii in 1:length(mycounts)){
          model_fit_single$raw_counts[[labeledObj]][['down']][[names(mycounts)[ii]]] = 
            model_fit_single$raw_counts[[labeledObj]][['down']][[names(mycounts)[ii]]] + mycounts[ii]
        }
        if (any(rowSums(subClassMat_total) == 1)){
          mycounts_norm = mycounts / sum(subClassMat_total[rowSums(subClassMat_total) == 1,])
        }else{ # no singles; my counts will be all zeros
          mycounts_norm = mycounts
        }
        mycounts_norm[is.na(mycounts_norm)] = 0
        for (ii in 1:length(mycounts_norm)){
          model_fit_single$normalized_counts[[labeledObj]][['down']][[names(mycounts_norm)[ii]]] = 
            model_fit_single$normalized_counts[[labeledObj]][['down']][[names(mycounts_norm)[ii]]] + mycounts_norm[ii]
        }
      }else{
        model_explained_single$DOWN_yes[condInd] = NA
        model_explained_single$DOWN_no[condInd] = NA
        
      }
      
      # assign the literal model (all obj counted independently)
      # let's start with simple literal counting: multi-obj RNAi were counted multiple times and 
      # multi-obj DE was counted as independent DE; if it is too noisy, we can weight the multi cases 
      # by equally divide their contribution (multi-obj RNAi was counted as 0.5+0.5, etc)
      subClassMat_wtd = subClassMat / rowSums(subClassMat)
      mycounts = colSums(subClassMat_wtd,na.rm = T)
      for (jj in 1:length(labeledObj)){
        for (ii in 1:length(mycounts)){
          model_fit_all$raw_counts[[labeledObj[jj]]][['down']][[names(mycounts)[ii]]] = 
            model_fit_all$raw_counts[[labeledObj[jj]]][['down']][[names(mycounts)[ii]]] + (mycounts[ii]/length(labeledObj))
        }
        mycounts_norm = mycounts / sum(subClassMat_total_wtd, na.rm = T)
        mycounts_norm[is.na(mycounts_norm)] = 0
        for (ii in 1:length(mycounts_norm)){
          model_fit_all$normalized_counts[[labeledObj[jj]]][['down']][[names(mycounts_norm)[ii]]] = 
            model_fit_all$normalized_counts[[labeledObj[jj]]][['down']][[names(mycounts_norm)[ii]]] + (mycounts_norm[ii]/length(labeledObj))
        }
      }
    }else{stop('there is a bug')}
    
  }
  # remove NA in singles 
  model_explained_single = model_explained_single[!rowAlls(is.na(model_explained_single[,2:5])),]
  
  # total explained rate
  tmp = model_explained[,2:5]
  tmp = tmp / rowSums(tmp)
  rewire_rate_rand = rowSums(tmp[,c(1,3)])
  rewire_rate_rand[is.na(rewire_rate_rand)] = 0
  # the overall average explained rate (up and down together) is 
  rand_stat$default_CR_full$rand_rate = c(rand_stat$default_CR_full$rand_rate, mean(rewire_rate_rand))
  rand_stat$default_CR_full$rand_rate2 = c(rand_stat$default_CR_full$rand_rate2, sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5]))
  # up
  tmp = model_explained[,2:3]
  tmp = tmp / rowSums(tmp)
  rand_stat$default_CR_full$rand_rate_up = c(rand_stat$default_CR_full$rand_rate_up, mean(tmp$UP_yes, na.rm = T))
  rand_stat$default_CR_full$rand_rate2_up = c(rand_stat$default_CR_full$rand_rate2_up, sum(model_explained[,c(2)]) / sum(model_explained[,2:3]))
  # down
  tmp = model_explained[,4:5]
  tmp = tmp / rowSums(tmp)
  rand_stat$default_CR_full$rand_rate_down = c(rand_stat$default_CR_full$rand_rate_down, mean(tmp$DOWN_yes, na.rm = T))
  rand_stat$default_CR_full$rand_rate2_down = c(rand_stat$default_CR_full$rand_rate2_down, sum(model_explained[,c(4)]) / sum(model_explained[,4:5]))
  
  
  # total explained rate - single
  tmp = model_explained_single[,2:5]
  tmp = tmp / rowSums(tmp)
  rewire_rate_rand = rowSums(tmp[,c(1,3)])
  rewire_rate_rand[is.na(rewire_rate_rand)] = 0
  # the overall average explained rate (up and down together) is 
  rand_stat$default_CR_single$rand_rate = c(rand_stat$default_CR_single$rand_rate, mean(rewire_rate_rand))
  rand_stat$default_CR_single$rand_rate2 = c(rand_stat$default_CR_single$rand_rate2, sum(model_explained_single[,c(2,4)]) / sum(model_explained_single[,2:5]))
  # up
  tmp = model_explained_single[,2:3]
  tmp = tmp / rowSums(tmp)
  rand_stat$default_CR_single$rand_rate_up = c(rand_stat$default_CR_single$rand_rate_up, mean(tmp$UP_yes, na.rm = T))
  rand_stat$default_CR_single$rand_rate2_up = c(rand_stat$default_CR_single$rand_rate2_up, sum(model_explained_single[,c(2)]) / sum(model_explained_single[,2:3]))
  # down
  tmp = model_explained_single[,4:5]
  tmp = tmp / rowSums(tmp)
  rand_stat$default_CR_single$rand_rate_down = c(rand_stat$default_CR_single$rand_rate_down, mean(tmp$DOWN_yes, na.rm = T))
  rand_stat$default_CR_single$rand_rate2_down = c(rand_stat$default_CR_single$rand_rate2_down, sum(model_explained_single[,c(4)]) / sum(model_explained_single[,4:5]))
  
  
  # the edge quantification 
  rand_stat$edge_single$raw_count[paste('rand',nn,sep = '')] = make_g_tbl(model_fit_single$raw_counts)$proportion
  rand_stat$edge_single$norm_count[paste('rand',nn,sep = '')] = make_g_tbl(model_fit_single$normalized_counts)$proportion
  rand_stat$edge_full$raw_count[paste('rand',nn,sep = '')] = make_g_tbl(model_fit_all$raw_counts)$proportion
  rand_stat$edge_full$norm_count[paste('rand',nn,sep = '')] = make_g_tbl(model_fit_all$normalized_counts)$proportion
  
  
  print(nn)
  # list(list(mean(rewire_rate_rand), sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5])))
}

save(file = 'randomization_result_edge_quantification.Rdata',rand_stat)
hist(rand_stat$default_CR_full$rand_rate)
hist(rand_stat$default_CR_single$rand_rate)
pdf(paste('figures/full_randomization_results.pdf',sep = ''),width = 7,height = 6)
hist(rand_stat$default_CR_full$rand_rate, 
     main = paste('CR model - full data\n','p <',signif((1+sum(rand_stat$default_CR_full$rand_rate>=obs_rate))/(1+length(rand_stat$default_CR_full$rand_rate)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate)*0.9,max(c(rand_stat$default_CR_full$rand_rate,obs_rate))*1.1),
     xlab = 'condition-wise average percentage of DEG explained')
abline(v = obs_rate)

hist(rand_stat$default_CR_full$rand_rate2, main = paste('CR model - full data\n', 'p <',signif((1+sum(rand_stat$default_CR_full$rand_rate2>=obs_total_DE_rate))/(1+length(rand_stat$default_CR_full$rand_rate2)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate2)*0.9,max(c(rand_stat$default_CR_full$rand_rate2,obs_total_DE_rate))*1.1),
     xlab = 'percentage of total DEG explained')
abline(v = obs_total_DE_rate)

hist(rand_stat$default_CR_full$rand_rate_up, main = paste('CR model - full data\n','p <',signif((1+sum(rand_stat$default_CR_full$rand_rate_up>=obs_rate_up))/(1+length(rand_stat$default_CR_full$rand_rate_up)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate_up)*0.9,max(c(rand_stat$default_CR_full$rand_rate_up,obs_rate_up))*1.1),
     xlab = 'condition-wise average percentage of up DEG explained')
abline(v = obs_rate_up)

hist(rand_stat$default_CR_full$rand_rate2_up, main = paste('CR model - full data\n','p <',signif((1+sum(rand_stat$default_CR_full$rand_rate2_up>=obs_total_DE_rate_up))/(1+length(rand_stat$default_CR_full$rand_rate2_up)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate2_up)*0.9,max(c(rand_stat$default_CR_full$rand_rate2_up,obs_total_DE_rate_up))*1.1),
     xlab = 'percentage of total UP DEG explained')
abline(v = obs_total_DE_rate_up)

hist(rand_stat$default_CR_full$rand_rate_down, main = paste('CR model - full data\n', 'p <',signif((1+sum(rand_stat$default_CR_full$rand_rate_down>=obs_rate_down))/(1+length(rand_stat$default_CR_full$rand_rate_down)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate_down)*0.9,max(c(rand_stat$default_CR_full$rand_rate_down,obs_rate_down))*1.1),
     xlab = 'condition-wise average percentage of down DEG explained')
abline(v = obs_rate_down)

hist(rand_stat$default_CR_full$rand_rate2_down, main = paste('CR model - full data\n', 'p <',signif((1+sum(rand_stat$default_CR_full$rand_rate2_down>=obs_total_DE_rate_down))/(1+length(rand_stat$default_CR_full$rand_rate2_down)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate2_down)*0.9,max(c(rand_stat$default_CR_full$rand_rate2_down,obs_total_DE_rate_down))*1.1),
     xlab = 'percentage of total down DEG explained')
abline(v = obs_total_DE_rate_down)


hist(rand_stat$default_CR_single$rand_rate, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate>=obs_rate_single))/(1+length(rand_stat$default_CR_single$rand_rate)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate)*0.9,max(c(rand_stat$default_CR_single$rand_rate,obs_rate_single))*1.1),
     xlab = 'condition-wise average percentage of DEG explained')
abline(v = obs_rate_single)

hist(rand_stat$default_CR_single$rand_rate2, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate2>=obs_total_DE_rate_single))/(1+length(rand_stat$default_CR_single$rand_rate2)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate2)*0.9,max(c(rand_stat$default_CR_single$rand_rate2,obs_total_DE_rate_single))*1.1),
     xlab = 'percentage of total DEG explained')
abline(v = obs_total_DE_rate_single)

hist(rand_stat$default_CR_single$rand_rate_up, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate_up>=obs_rate_up_single))/(1+length(rand_stat$default_CR_single$rand_rate_up)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate_up)*0.9,max(c(rand_stat$default_CR_single$rand_rate_up,obs_rate_up_single))*1.1),
     xlab = 'condition-wise average percentage of up DEG explained')
abline(v = obs_rate_up_single)

hist(rand_stat$default_CR_single$rand_rate2_up, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate2_up>=obs_total_DE_rate_up_single))/(1+length(rand_stat$default_CR_single$rand_rate2_up)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate2_up)*0.9,max(c(rand_stat$default_CR_single$rand_rate2_up,obs_total_DE_rate_up_single))*1.1),
     xlab = 'percentage of total UP DEG explained')
abline(v = obs_total_DE_rate_up_single)

hist(rand_stat$default_CR_single$rand_rate_down, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate_down>=obs_rate_down_single))/(1+length(rand_stat$default_CR_single$rand_rate_down)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate_down)*0.9,max(c(rand_stat$default_CR_single$rand_rate_down,obs_rate_down_single))*1.1),
     xlab = 'condition-wise average percentage of down DEG explained')
abline(v = obs_rate_down_single)

hist(rand_stat$default_CR_single$rand_rate2_down, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate2_down>=obs_total_DE_rate_down_single))/(1+length(rand_stat$default_CR_single$rand_rate2_down)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate2_down)*0.9,max(c(rand_stat$default_CR_single$rand_rate2_down,obs_total_DE_rate_down_single))*1.1),
     xlab = 'percentage of total down DEG explained')
abline(v = obs_total_DE_rate_down_single)

# quantify each edge and plot significance 
sig_cutoff = 0.2

# note: in rare cases, one class of RNAi may be missing. for example, the protein syn only has 16 genes so it is 
# possible that after shuffling, all RNAi conditions gets other categories and no protein syn RNAi. In this case, 
# the protein syn RNAi edges will get NA and should be just ignored in the p-value calculation (to not bias anything
# as putting all proportion as zero is also problematic as we normalizing within RNAi or within a RNAi-DEG pair)

# single and raw
g_tbl_single_raw$p_value = NA
for (i in 1:nrow(g_tbl_single_raw)){
  g_tbl_single_raw$p_value[i] = (1+sum(rand_stat$edge_single$raw_count[i,4:ncol(rand_stat$edge_single$raw_count)] >= g_tbl_single_raw$proportion[i], na.rm = T)) / (sum(!is.na(rand_stat$edge_single$raw_count[i,4:ncol(rand_stat$edge_single$raw_count)]))+1)
}
filtered_g = g_tbl_single_raw[g_tbl_single_raw$p_value < sig_cutoff,]
colnames(filtered_g) = c('from','to','type','width','p_value')
graph <- graph_from_data_frame(filtered_g, directed = TRUE)
# add a label attribute
E(graph)$label <- paste("Width:", signif(E(graph)$width,2), "P-value:", signif(E(graph)$p_value,2))
# Drawing the graph
plot(graph,
     edge.arrow.size = E(graph)$width / max(E(graph)$width),
     edge.arrow.width = E(graph)$width / max(E(graph)$width),
     edge.label = E(graph)$label,
     edge.label.cex = 0.8,
     edge.width = E(graph)$width * 30,
     edge.arrow.mode = '->',
     edge.color = ifelse(E(graph)$type == 'up', 'red', 'blue'),
     vertex.color = "gold", 
     vertex.size = 20,
     vertex.frame.color = "gray", 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 1, 
     main = paste("single - raw - p <", sig_cutoff))


# single and normalized
g_tbl_single_norm$p_value = NA
for (i in 1:nrow(g_tbl_single_norm)){
  g_tbl_single_norm$p_value[i] = (1+sum(rand_stat$edge_single$norm_count[i,4:ncol(rand_stat$edge_single$norm_count)] >= g_tbl_single_norm$proportion[i], na.rm = T)) / (sum(!is.na(rand_stat$edge_single$norm_count[i,4:ncol(rand_stat$edge_single$norm_count)]))+1)
}
filtered_g = g_tbl_single_norm[g_tbl_single_norm$p_value < sig_cutoff,]
colnames(filtered_g) = c('from','to','type','width','p_value')
graph <- graph_from_data_frame(filtered_g, directed = TRUE)
# add a label attribute
E(graph)$label <- paste("Width:", signif(E(graph)$width,2), "P-value:", signif(E(graph)$p_value,2))
# Drawing the graph
plot(graph,
     edge.arrow.size = E(graph)$width / max(E(graph)$width),
     edge.arrow.width = E(graph)$width / max(E(graph)$width),
     edge.label = E(graph)$label,
     edge.label.cex = 0.8,
     edge.width = E(graph)$width * 30,
     edge.arrow.mode = '->',
     edge.color = ifelse(E(graph)$type == 'up', 'red', 'blue'),
     vertex.color = "gold", 
     vertex.size = 20,
     vertex.frame.color = "gray", 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 1, 
     main = paste("single - normalized - p <", sig_cutoff))

# full and raw
g_tbl_full_raw$p_value = NA
for (i in 1:nrow(g_tbl_full_raw)){
  g_tbl_full_raw$p_value[i] = (1+sum(rand_stat$edge_full$raw_count[i,4:ncol(rand_stat$edge_full$raw_count)] >= g_tbl_full_raw$proportion[i], na.rm = T)) / (sum(!is.na(rand_stat$edge_full$raw_count[i,4:ncol(rand_stat$edge_full$raw_count)]))+1)
}
filtered_g = g_tbl_full_raw[g_tbl_full_raw$p_value < sig_cutoff,]
colnames(filtered_g) = c('from','to','type','width','p_value')
graph <- graph_from_data_frame(filtered_g, directed = TRUE)
# add a label attribute
E(graph)$label <- paste("Width:", signif(E(graph)$width,2), "P-value:", signif(E(graph)$p_value,2))
# Drawing the graph
plot(graph,
     edge.arrow.size = E(graph)$width / max(E(graph)$width),
     edge.arrow.width = E(graph)$width / max(E(graph)$width),
     edge.label = E(graph)$label,
     edge.label.cex = 0.8,
     edge.width = E(graph)$width * 30,
     edge.arrow.mode = '->',
     edge.color = ifelse(E(graph)$type == 'up', 'red', 'blue'),
     vertex.color = "gold", 
     vertex.size = 20,
     vertex.frame.color = "gray", 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 1, 
     main = paste("full - raw - p <", sig_cutoff))


# full and normalized
g_tbl_full_norm$p_value = NA
for (i in 1:nrow(g_tbl_full_norm)){
  g_tbl_full_norm$p_value[i] = (1+sum(rand_stat$edge_full$norm_count[i,4:ncol(rand_stat$edge_full$norm_count)] >= g_tbl_full_norm$proportion[i], na.rm = T)) / (sum(!is.na(rand_stat$edge_full$norm_count[i,4:ncol(rand_stat$edge_full$norm_count)]))+1)
}
filtered_g = g_tbl_full_norm[g_tbl_full_norm$p_value < sig_cutoff,]
colnames(filtered_g) = c('from','to','type','width','p_value')
graph <- graph_from_data_frame(filtered_g, directed = TRUE)
# add a label attribute
E(graph)$label <- paste("Width:", signif(E(graph)$width,2), "P-value:", signif(E(graph)$p_value,2))
# Drawing the graph
plot(graph,
     edge.arrow.size = E(graph)$width / max(E(graph)$width),
     edge.arrow.width = E(graph)$width / max(E(graph)$width),
     edge.label = E(graph)$label,
     edge.label.cex = 0.8,
     edge.width = E(graph)$width * 30,
     edge.arrow.mode = '->',
     edge.color = ifelse(E(graph)$type == 'up', 'red', 'blue'),
     vertex.color = "gold", 
     vertex.size = 20,
     vertex.frame.color = "gray", 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 1, 
     main = paste("full - normalized - p <", sig_cutoff))

dev.off()


# finally, save the results
edge_quant = g_tbl_single_raw
colnames(edge_quant)[4] = 'prop_single_raw'
colnames(edge_quant)[5] = 'prop_p_single_raw'
edge_quant$prop_single_norm = g_tbl_single_norm$proportion
edge_quant$prop_full_raw = g_tbl_full_raw$proportion
edge_quant$prop_full_norm = g_tbl_full_norm$proportion
edge_quant = edge_quant[,c(1:4,6:ncol(edge_quant),5)]
edge_quant$prop_p_single_norm = g_tbl_single_norm$p_value
edge_quant$prop_p_full_raw = g_tbl_full_raw$p_value
edge_quant$prop_p_full_norm = g_tbl_full_norm$p_value
write.csv(edge_quant,'figures/edge_quantification.csv')





# # combine the two ways of edge assessment and make the final figures of edge quantification
# edge_info = read.csv('figures/edge_quantification.csv') # this is the result when considering the proportion of an interaction among all interactions of a core function perturbation
# 
# # visualization
# 
# # first to see how robust it is among all various combinations of setups
# edge_info = edge_info[order(edge_info$RNAi), ]
# 
# 
# # these settings include the following options: 
# # how to quantify significant: by overrepresentation of a direction of interaction or overrepresentation of an specific interaction 
# # what to use: normalized (per condition proportion) or raw (total proportion)
# # what data to use: single-core-function genes only (single) or all genes with core functions (full)
# # now visualize every combinations 
# return_graph <- function(filter_by, display_by, sig_cutoff){
#   # decide filter by
#   edge_info$best_p = pmin(edge_info[,filter_by[1]], edge_info[,filter_by[2]])
#   filtered_g = edge_info[edge_info$best_p < sig_cutoff,]
#   # decide display by
#   filtered_g = filtered_g[,c("RNAi", "DEG", "type", display_by,"best_p")]
#   colnames(filtered_g) = c('from','to','type','width','p_value')
#   graph <- graph_from_data_frame(filtered_g, directed = TRUE)
#   # add a label attribute
#   E(graph)$label <- paste("Width:", signif(E(graph)$width,2), "P-value:", signif(E(graph)$p_value,2))
#   # set the color based on -log10pvalue 
#   signed_pvalues = -log10(E(graph)$p_value) * ifelse(E(graph)$type == 'up', 1, -1)
#   mycolors = "white"
#   colorScale = max(abs(signed_pvalues))
#   
#   # assign positive values
#   colors <- colorRampPalette(c("white", "red"))
#   # Generate the color vector
#   mycolors[signed_pvalues>0] <- colors(100)[findInterval(signed_pvalues[signed_pvalues>0], 
#                                                          seq(0, colorScale, 
#                                                              length = 100))]
#   
#   # assign negative values
#   colors <- colorRampPalette(c("blue", "white"))
#   # Generate the color vector
#   mycolors[signed_pvalues<0] <- colors(100)[findInterval(signed_pvalues[signed_pvalues<0], 
#                                                          seq(-colorScale,0, 
#                                                              length = 100))]
#   
#   # Drawing the graph
#   l <- layout_in_circle(graph)
#   p = plot(graph, layout = l,
#            edge.arrow.size = E(graph)$width / max(E(graph)$width) * 1.5,
#            edge.arrow.width = E(graph)$width / max(E(graph)$width) * 1.5,
#            edge.label = E(graph)$label,
#            edge.label.cex = 0.8,
#            edge.width = E(graph)$width / max(E(graph)$width) * 15,
#            edge.arrow.mode = '->',
#            edge.color = mycolors,
#            vertex.color = "gold", 
#            vertex.size = 20,
#            vertex.frame.color = "gray", 
#            vertex.label.color = "black", 
#            vertex.label.cex = 0.8, 
#            vertex.label.dist = 1, 
#            edge.curved=0.5,
#            main = paste('filtered by:', paste(filter_by,collapse = ' '),'\n',
#                         'displayed by:',display_by,'\n',
#                         "edge proprotion or edge direction p <", sig_cutoff,'\n',
#                         'color scale is', colorScale))
#   return(p)
# }
# 
# 
# pdf('figures/fitted_CR_interaction_diagram.pdf',width = 10)
# # the abbreviations refer to 
# # filter by: whether or not to display an edge (interaction) is filtered by pvalues calculated in this setup
# # display by: the thickness of the edge is displayed with values calculated by this setup
# # prop_p: p-values calculated with the proportion of an edge among all edges for a core function perturbation 
# # full: all genes with core function associations 
# # single: only genes with single core function associations
# # CR_p: p-values calculated with the proportion of edge direction (up or down, compensation or repression) between two core functions 
# # norm: normalized proportion (first calculate proportion within each condition and then get average proportion)
# # prop: proportion of an edge among all edges for a core function perturbation 
# # raw: raw proportion (aggregate DEGs in all conditions and then calculate the proportion of one type of interaction)
# sig_cutoff = 0.1
# 
# # for a preliminary test in humans, we skip the randomization p-value based on testing edge direction (in worm study, we
# # visualized all edges that is either significant for )
# filter_by = c('prop_p_full_norm','prop_p_full_norm')
# display_by = 'prop_full_norm'
# return_graph(filter_by, display_by,sig_cutoff)
# 
# sig_cutoff = 0.1
# filter_by = c('prop_p_single_norm','prop_p_single_norm')
# display_by = 'prop_single_norm'
# return_graph(filter_by, display_by,sig_cutoff)
# 
# sig_cutoff = 0.1
# filter_by = c('prop_p_full_raw','prop_p_full_raw')
# display_by = 'prop_full_raw'
# return_graph(filter_by, display_by,sig_cutoff)
# 
# sig_cutoff = 0.1
# filter_by = c('prop_p_single_raw','prop_p_single_raw')
# display_by = 'prop_single_raw'
# return_graph(filter_by, display_by,sig_cutoff)
# 
# dev.off()