# edge quantification analysis of the interactions between core functions 

# this script focus on the edge quantification with testing hypothesis of overrepresentation of a regulatory direction (up vs. down) of an interacting edge
# the layout of the code is generally same as visualizing DEGs, which an additional block to quantify the edges 

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


# we get the conditions to analyze (icel_responsive (at least 2 up or 2 down DEG) and classified)
total_condition_analyzed = c()
for (i in 1:length(human1_resp)){
  if (rowSums(classMat[cell_annotations$perturbed_gene[match(human1_resp[i], cell_annotations$perturbation_id)],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
    total_condition_analyzed = c(total_condition_analyzed, human1_resp[i])
  }
}


# the following model variable is not important for the edge quantification; it is there for a historical reason
# we kept the original code for reproducibility purpose.
# start here -- 
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
# end here -- 


# the model_fit variable is the variable for quantifying different interacting edges
# for historical reasons, we call the interacting edges as a fitted model in the code

# add the fitted model - simple total count of DEG
# fitted model considers all possible interactions
# this is the raw proportion
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

# weighted total count of DEG - we normalize the total DEG count of each condition to 1 (to avoid confounding by conditions with super high DE number)
# (when we calculate, we DO NOT seperate the up and down genes)
# this is the normalized proportion 
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

# we also make the same quantification using only genes with unique core function association 
model_fit_single = list('raw_counts' = model_fit,
                        'normalized_counts' = model_fit_wtd)
model_fit_all = model_fit_single

# quantify the proportions for each types of up and down DEG core functions in each type of core function perturbations
for (condInd in 1:length(total_condition_analyzed)){
    # loop through each condition and get the up and down DEGs
    myCond = total_condition_analyzed[condInd]
    DEGs_up = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] > 0)]
    DEGs_down = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] < 0)]
    
    # determine affected core function 
    labeledObj = colnames(classMat)[classMat[cell_annotations$perturbed_gene[cell_annotations$perturbation_id == myCond],]==1]
    labeledObj = intersect(labeledObj, modelObj)
    
    # if the RNAi targeted gene is associated with any function, we start to quantify the core functions of its DEGs 
    if (length(labeledObj) > 0){
      
      
      affectedObj = c() # this part of code should be ignored - it is to quantify how many DEGs can be explained by "model" shown above. It is here for historical reason and we didnt delete because we need to keep track of the history 
      for (i in 1:length(labeledObj)){ # this part of code should be ignored - it is to quantify how many DEGs can be explained by "model" shown above. It is here for historical reason and we didnt delete because we need to keep track of the history 
        affectedObj = c(affectedObj, model[[labeledObj[i]]])
      }
      
      # compare with DEG 
      
      # get the core function matrix for all DEGs
      subClassMat_total = classMat[rownames(classMat) %in% c(DEGs_up,DEGs_down), modelObj]
      # for weighted counts, we normalize the matrix
      subClassMat_total_wtd = subClassMat_total / rowSums(subClassMat_total)
      
      
      
      # part A: let's first analyze the up regulated genes
      
      subClassMat = classMat[rownames(classMat) %in% DEGs_up, modelObj]
      
      # ignore the following two lines
      model_explained$UP_yes[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) > 0) # this part of code should be ignored - it is to quantify how many DEGs can be explained by "model" shown above. It is here for historical reason and we didnt delete because we need to keep track of the history 
      model_explained$UP_no[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) == 0) # this part of code should be ignored - it is to quantify how many DEGs can be explained by "model" shown above. It is here for historical reason and we didnt delete because we need to keep track of the history 

      # perform edge quantification

      # we first quantify using only single-core-function genes 
      # calculate for "model_fit_single"
      
      # check if it is a single-core-function association (of the RNAi targeted gene) 
      # (if it is a multiple-core-function gene, the count has to be weighted (see the "else" clause))
      if (length(labeledObj) == 1){
        
        # subset the classification matrix for single-function genes - this is used to calculate the quantification based on only single-function genes
        subClassMat_single = subClassMat[rowSums(subClassMat) == 1,]
        
        # ignore the following four lines - this is to quantify genes explained by "model" defined in the above using only genes with unique core function 
        subClassMat_single2 = subClassMat # mask all multiple genes # this part of code should be ignored - it is to quantify how many DEGs can be explained by "model" shown above. It is here for historical reason and we didnt delete because we need to keep track of the history 
        subClassMat_single2[rowSums(subClassMat_single2[,1:5]) > 1,1:5] = 0 # this part of code should be ignored - it is to quantify how many DEGs can be explained by "model" shown above. It is here for historical reason and we didnt delete because we need to keep track of the history 
        model_explained_single$UP_yes[condInd] = sum(rowSums(subClassMat_single2[,affectedObj,drop = F]) > 0) # this part of code should be ignored - it is to quantify how many DEGs can be explained by "model" shown above. It is here for historical reason and we didnt delete because we need to keep track of the history 
        model_explained_single$UP_no[condInd] = sum(rowSums(subClassMat_single2[,affectedObj,drop = F]) == 0) # this part of code should be ignored - it is to quantify how many DEGs can be explained by "model" shown above. It is here for historical reason and we didnt delete because we need to keep track of the history 
        
        # count the number of interactions (DEGs and weighted DEGs) for each core function 
       
        # first count by genes with single core function associations
        
        # raw counts (proportion calculated by total DEGs belonging to a type of interaction over total DEGs in all conditions)
        mycounts = colSums(subClassMat_single)
        for (ii in 1:length(mycounts)){
          model_fit_single$raw_counts[[labeledObj]][['up']][[names(mycounts)[ii]]] = 
            model_fit_single$raw_counts[[labeledObj]][['up']][[names(mycounts)[ii]]] + mycounts[ii] # calculate the total across conditions
        }
        
        # normalized counts (proportion calculated by DEG proportions within each conditions first, and then averaged over all conditions)
        if (any(rowSums(subClassMat_total) == 1)){ # has any DEG with single core function association
          mycounts_norm = mycounts / sum(subClassMat_total[rowSums(subClassMat_total) == 1,]) # calculate the proportion within this condition by dividing the total number of DEGs (only consider DEG with unique core function association)
        }else{ # otherwise mycounts will be all zeros - no need to devide anything
          mycounts_norm = mycounts
        }
        mycounts_norm[is.na(mycounts_norm)] = 0 # avoid nan (they were 0/0)
        for (ii in 1:length(mycounts_norm)){
          model_fit_single$normalized_counts[[labeledObj]][['up']][[names(mycounts_norm)[ii]]] = 
            model_fit_single$normalized_counts[[labeledObj]][['up']][[names(mycounts_norm)[ii]]] + mycounts_norm[ii] # calculate the total across conditions
        }
        # the total counts calculated here will be converted to proportion quantities when deviding total number of conditions 
        
      }else{ # otherwise excluded from analysis 
        model_explained_single$UP_yes[condInd] = NA
        model_explained_single$UP_no[condInd] = NA
      }
      
      # we next quantify using all genes with core function associations 
      # calculate for "model_fit_all"
      
      # it is more complex when it comes to an RNAi targeted gene with multiple core functions
      # we need to count all core functions independently, meaning that multi-function RNAi will be
      # counted multiple times and multi-function DE was counted as independent DE; However, we will
      # equally divide their contribution (multi-obj RNAi was counted as 0.5+0.5, etc).
      subClassMat_wtd = subClassMat / rowSums(subClassMat) # devide by number of core functions for each gene, so such weighted count will be used next
      
      # raw counts (proportion calculated by total DEGs belonging to a type of interaction over total DEGs in all conditions)
      mycounts = colSums(subClassMat_wtd,na.rm = T)
      for (jj in 1:length(labeledObj)){
        for (ii in 1:length(mycounts)){
          model_fit_all$raw_counts[[labeledObj[jj]]][['up']][[names(mycounts)[ii]]] = 
            model_fit_all$raw_counts[[labeledObj[jj]]][['up']][[names(mycounts)[ii]]] + (mycounts[ii] /length(labeledObj)) # againt, devide the contribution by the number of appearance
        }
        
        # normalized counts (proportion calculated by DEG proportions within each conditions first, and then averaged over all conditions)
        mycounts_norm = mycounts / sum(subClassMat_total_wtd,na.rm = T)
        mycounts_norm[is.na(mycounts_norm)] = 0
        for (ii in 1:length(mycounts_norm)){
          model_fit_all$normalized_counts[[labeledObj[jj]]][['up']][[names(mycounts_norm)[ii]]] = 
            model_fit_all$normalized_counts[[labeledObj[jj]]][['up']][[names(mycounts_norm)[ii]]] + (mycounts_norm[ii]/length(labeledObj))
        }
      }
      # now we are done with interactions involving up regulated genes
      
      
      # part B: let's then analyze the up regulated genes
      # basically repeat the above codes. we didnt put it into a for loop to code more decently..
      
      # subset down genes
      subClassMat = classMat[rownames(classMat) %in% DEGs_down, modelObj]
      
      # calculate the DEG explained by model
      model_explained$DOWN_yes[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) > 0) # ignore this
      model_explained$DOWN_no[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) == 0) # ignore this
    
      # quantify the edges
      # check if it is a single model 
      if (length(labeledObj) == 1){
        # only consider single-core-function genes
        subClassMat_single = subClassMat[rowSums(subClassMat) == 1,]
        
        # calculate the DEG explained by model - ignore this
        subClassMat_single2 = subClassMat # mask all multiple genes - # ignore this
        subClassMat_single2[rowSums(subClassMat_single2[,1:5]) > 1,1:5] = 0 # ignore this
        model_explained_single$DOWN_yes[condInd] = sum(rowSums(subClassMat_single2[,setdiff(modelObj, affectedObj),drop = F]) > 0) # ignore this
        model_explained_single$DOWN_no[condInd] = sum(rowSums(subClassMat_single2[,setdiff(modelObj, affectedObj),drop = F]) == 0) # ignore this
        
        # quantify using only genes with unique core function associations first
        # same as above, quantify for raw counts
        mycounts = colSums(subClassMat_single)
        for (ii in 1:length(mycounts)){
          model_fit_single$raw_counts[[labeledObj]][['down']][[names(mycounts)[ii]]] = 
            model_fit_single$raw_counts[[labeledObj]][['down']][[names(mycounts)[ii]]] + mycounts[ii]
        }
        # quantify for normalized counts
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
      
      # same as above, now do this for all genes with core function associations but apply a weighting scheme 
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
      # ignore codes below this..
      
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
# remove NA in singles (we kept them in calculation because they are used for manual inspection when developing the code)
model_explained_single = model_explained_single[!rowAlls(is.na(model_explained_single[,2:5])),]


# ignore the following codes - they are for calculating the proportion of DEGs explained by "model" defined in very beginning. Old and unused code
# start here -- 
# the overall average explained rate (up and down together) is 
tmp = model_explained[,2:5]
tmp = tmp / rowSums(tmp)
rewire_rate = rowSums(tmp[,c(1,3)])
rewire_rate[is.na(rewire_rate)] = 0
obs_rate = mean(rewire_rate)
obs_total_DE_rate = sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5])
# the overall average explained rate (up and down seperately) is 
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

# -- end here 



# analyze the edge quantification 

# quick visualization of the edges by proportion
library(igraph)
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
  # normalize within each type of RNAi-DEG pair (i.e., between two core functions)
  for (i in 1:length(modelObj)){
    for(j in 1:length(modelObj)){
      g_tbl$proportion[g_tbl$RNAi == modelObj[i] & g_tbl$DEG == modelObj[j]] =  g_tbl$proportion[g_tbl$RNAi == modelObj[i] & g_tbl$DEG == modelObj[j]] / sum(g_tbl$proportion[g_tbl$RNAi == modelObj[i] & g_tbl$DEG == modelObj[j]])
    }
  }
  return(g_tbl)
}

# visualization (in a graph)
g = model_fit_single$raw_counts
g_tbl = make_g_tbl(g)
# Creating the graph object
graph <- graph_from_data_frame(g_tbl[,1:2], directed = TRUE)
# Setting edge width and type
E(graph)$width <- g_tbl$proportion
E(graph)$type <- g_tbl$type

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



# finally, assess whether each interaction edge is statistically significant

# use randomization to assess the significance - this script uses randomization of the FBA classification of gene core functions
nRand = 10000
set.seed(1126)
# we have a ton to calculate
rand_stat = list()
# default_CR_full and default_CR_single can be ignored - they are for "model" defined above
rand_stat$default_CR_full = list(rand_rate = c(), rand_rate2 = c(), rand_rate_up = c(), rand_rate2_up = c(), rand_rate_down = c(), rand_rate2_down = c())
rand_stat$default_CR_single = list(rand_rate = c(), rand_rate2 = c(), rand_rate_up = c(), rand_rate2_up = c(), rand_rate_down = c(), rand_rate2_down = c())
# the edge frequencies of each interaction 
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
  
  # calculate explained gene numbers - ignore this 
  model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
  model_explained$UP_yes = 0
  model_explained$UP_no = 0
  model_explained$DOWN_yes = 0
  model_explained$DOWN_no = 0
  
  model_explained_single = model_explained # ignore
  
  # calculate all the numbers using randomized data 
  # the following codes are same as above
  
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
  # weighted total count of DEG - we normalize the total DEG count of each condition to 1 (to avoid confounding by conditions with super high DE number)
  # (when we calculate, we DO NOT separate the up and down genes)
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
      
      # quantify using only single-core-function genes
      if (length(labeledObj) == 1){
        # only consider single-category genes
        subClassMat_single = subClassMat[rowSums(subClassMat) == 1,]
        
        # calculate the DEG explained by model - ignore
        subClassMat_single2 = subClassMat # mask all multiple genes
        subClassMat_single2[rowSums(subClassMat_single2[,1:5]) > 1,1:5] = 0
        model_explained_single$UP_yes[condInd] = sum(rowSums(subClassMat_single2[,affectedObj,drop = F]) > 0)
        model_explained_single$UP_no[condInd] = sum(rowSums(subClassMat_single2[,affectedObj,drop = F]) == 0)
        
        # quantify edges
        mycounts = colSums(subClassMat_single)
        for (ii in 1:length(mycounts)){
          model_fit_single$raw_counts[[labeledObj]][['up']][[names(mycounts)[ii]]] = 
            model_fit_single$raw_counts[[labeledObj]][['up']][[names(mycounts)[ii]]] + mycounts[ii]
        }
        if (any(rowSums(subClassMat_total) == 1)){
          mycounts_norm = mycounts / sum(subClassMat_total[rowSums(subClassMat_total) == 1,])
        }else{ 
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
      
      
      # quantify using all genes with core function association
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
        }else{ 
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
      
      # all genes with core functions
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
  
  # now put all numbers 
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

# save results 
save(file = 'randomization_result_edge_direction.Rdata',rand_stat)
hist(rand_stat$default_CR_full$rand_rate) # ignore
hist(rand_stat$default_CR_single$rand_rate) # ignore
pdf(paste('figures/full_randomization_results_edge_direction.pdf',sep = ''),width = 7,height = 6)
# ignore this
hist(rand_stat$default_CR_full$rand_rate, 
     main = paste('CR model - full data\n','p <',signif((1+sum(rand_stat$default_CR_full$rand_rate>=obs_rate))/(1+length(rand_stat$default_CR_full$rand_rate)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate)*0.9,max(c(rand_stat$default_CR_full$rand_rate,obs_rate))*1.1),
     xlab = 'condition-wise average percentage of DEG explained')
abline(v = obs_rate)

# ignore this
hist(rand_stat$default_CR_full$rand_rate2, main = paste('CR model - full data\n', 'p <',signif((1+sum(rand_stat$default_CR_full$rand_rate2>=obs_total_DE_rate))/(1+length(rand_stat$default_CR_full$rand_rate2)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate2)*0.9,max(c(rand_stat$default_CR_full$rand_rate2,obs_total_DE_rate))*1.1),
     xlab = 'percentage of total DEG explained')
abline(v = obs_total_DE_rate)

# ignore this
hist(rand_stat$default_CR_full$rand_rate_up, main = paste('CR model - full data\n','p <',signif((1+sum(rand_stat$default_CR_full$rand_rate_up>=obs_rate_up))/(1+length(rand_stat$default_CR_full$rand_rate_up)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate_up)*0.9,max(c(rand_stat$default_CR_full$rand_rate_up,obs_rate_up))*1.1),
     xlab = 'condition-wise average percentage of up DEG explained')
abline(v = obs_rate_up)

# ignore this
hist(rand_stat$default_CR_full$rand_rate2_up, main = paste('CR model - full data\n','p <',signif((1+sum(rand_stat$default_CR_full$rand_rate2_up>=obs_total_DE_rate_up))/(1+length(rand_stat$default_CR_full$rand_rate2_up)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate2_up)*0.9,max(c(rand_stat$default_CR_full$rand_rate2_up,obs_total_DE_rate_up))*1.1),
     xlab = 'percentage of total UP DEG explained')
abline(v = obs_total_DE_rate_up)

# ignore this
hist(rand_stat$default_CR_full$rand_rate_down, main = paste('CR model - full data\n', 'p <',signif((1+sum(rand_stat$default_CR_full$rand_rate_down>=obs_rate_down))/(1+length(rand_stat$default_CR_full$rand_rate_down)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate_down)*0.9,max(c(rand_stat$default_CR_full$rand_rate_down,obs_rate_down))*1.1),
     xlab = 'condition-wise average percentage of down DEG explained')
abline(v = obs_rate_down)

# ignore this
hist(rand_stat$default_CR_full$rand_rate2_down, main = paste('CR model - full data\n', 'p <',signif((1+sum(rand_stat$default_CR_full$rand_rate2_down>=obs_total_DE_rate_down))/(1+length(rand_stat$default_CR_full$rand_rate2_down)),2)),
     xlim = c(min(rand_stat$default_CR_full$rand_rate2_down)*0.9,max(c(rand_stat$default_CR_full$rand_rate2_down,obs_total_DE_rate_down))*1.1),
     xlab = 'percentage of total down DEG explained')
abline(v = obs_total_DE_rate_down)

# ignore this
hist(rand_stat$default_CR_single$rand_rate, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate>=obs_rate_single))/(1+length(rand_stat$default_CR_single$rand_rate)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate)*0.9,max(c(rand_stat$default_CR_single$rand_rate,obs_rate_single))*1.1),
     xlab = 'condition-wise average percentage of DEG explained')
abline(v = obs_rate_single)

# ignore this
hist(rand_stat$default_CR_single$rand_rate2, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate2>=obs_total_DE_rate_single))/(1+length(rand_stat$default_CR_single$rand_rate2)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate2)*0.9,max(c(rand_stat$default_CR_single$rand_rate2,obs_total_DE_rate_single))*1.1),
     xlab = 'percentage of total DEG explained')
abline(v = obs_total_DE_rate_single)

# ignore this
hist(rand_stat$default_CR_single$rand_rate_up, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate_up>=obs_rate_up_single))/(1+length(rand_stat$default_CR_single$rand_rate_up)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate_up)*0.9,max(c(rand_stat$default_CR_single$rand_rate_up,obs_rate_up_single))*1.1),
     xlab = 'condition-wise average percentage of up DEG explained')
abline(v = obs_rate_up_single)

# ignore this
hist(rand_stat$default_CR_single$rand_rate2_up, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate2_up>=obs_total_DE_rate_up_single))/(1+length(rand_stat$default_CR_single$rand_rate2_up)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate2_up)*0.9,max(c(rand_stat$default_CR_single$rand_rate2_up,obs_total_DE_rate_up_single))*1.1),
     xlab = 'percentage of total UP DEG explained')
abline(v = obs_total_DE_rate_up_single)

# ignore this
hist(rand_stat$default_CR_single$rand_rate_down, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate_down>=obs_rate_down_single))/(1+length(rand_stat$default_CR_single$rand_rate_down)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate_down)*0.9,max(c(rand_stat$default_CR_single$rand_rate_down,obs_rate_down_single))*1.1),
     xlab = 'condition-wise average percentage of down DEG explained')
abline(v = obs_rate_down_single)

# ignore this
hist(rand_stat$default_CR_single$rand_rate2_down, main = paste('CR model - single data\n', 'p <',signif((1+sum(rand_stat$default_CR_single$rand_rate2_down>=obs_total_DE_rate_down_single))/(1+length(rand_stat$default_CR_single$rand_rate2_down)),2)),
     xlim = c(min(rand_stat$default_CR_single$rand_rate2_down)*0.9,max(c(rand_stat$default_CR_single$rand_rate2_down,obs_total_DE_rate_down_single))*1.1),
     xlab = 'percentage of total down DEG explained')
abline(v = obs_total_DE_rate_down_single)

# now plot the edge quantification results
# quantify each edge and plot significance 
sig_cutoff = 0.2 # we put an arbiturary pvalue cutoff as 0.2 here. 
# this is not the published version. The published version is plotted later with combined results

# note: in rare cases, one class of RNAi may be missing. for example, the protein syn only has 16 genes so it is 
# possible that after shuffling, all RNAi conditions gets other categories and no protein syn RNAi. In this case, 
# the protein syn RNAi edges will get NA and should be just ignored in the p-value calculation (to not bias anything
# as putting all proportion as zero is also problematic as we normalizing within RNAi or within a RNAi-DEG pair)

# calculate the pvalues
# for the case that use single-core-function only genes and use raw count to quantify
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
     edge.arrow.size = E(graph)$width / max(E(graph)$width) *1.5,
     edge.arrow.width = E(graph)$width / max(E(graph)$width)*1.5,
     edge.label = E(graph)$label,
     edge.label.cex = 0.8,
     edge.width = E(graph)$width * 10,
     edge.arrow.mode = '->',
     edge.color = ifelse(E(graph)$type == 'up', 'red', 'blue'),
     vertex.color = "gold", 
     vertex.size = 20,
     vertex.frame.color = "gray", 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 1, 
     main = paste("single - raw - p <", sig_cutoff))


# then for the case that uses single-core function genes and normalized counts (proportion within each condition and then average)
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
     edge.arrow.size = E(graph)$width / max(E(graph)$width)*1.5,
     edge.arrow.width = E(graph)$width / max(E(graph)$width)*1.5,
     edge.label = E(graph)$label,
     edge.label.cex = 0.8,
     edge.width = E(graph)$width * 10,
     edge.arrow.mode = '->',
     edge.color = ifelse(E(graph)$type == 'up', 'red', 'blue'),
     vertex.color = "gold", 
     vertex.size = 20,
     vertex.frame.color = "gray", 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 1, 
     main = paste("single - normalized - p <", sig_cutoff))

# now is using (full) all genes with core function associations and raw count (all conditions aggregated)
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
     edge.arrow.size = E(graph)$width / max(E(graph)$width)*1.5,
     edge.arrow.width = E(graph)$width / max(E(graph)$width)*1.5,
     edge.label = E(graph)$label,
     edge.label.cex = 0.8,
     edge.width = E(graph)$width * 10,
     edge.arrow.mode = '->',
     edge.color = ifelse(E(graph)$type == 'up', 'red', 'blue'),
     vertex.color = "gold", 
     vertex.size = 20,
     vertex.frame.color = "gray", 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8, 
     vertex.label.dist = 1, 
     main = paste("full - raw - p <", sig_cutoff))


# the last case is using full set of genes (all genes with core functions) and normalized count
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
     edge.arrow.size = E(graph)$width / max(E(graph)$width)*1.5,
     edge.arrow.width = E(graph)$width / max(E(graph)$width)*1.5,
     edge.label = E(graph)$label,
     edge.label.cex = 0.8,
     edge.width = E(graph)$width * 10,
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


# In this script, we calculate the p-value with relative frequency against the opposite direction instead of total frequency (
# calculate the proportion by up/(up+down) for each RNAi-DEG node pair)
# save results
edge_quant = g_tbl_single_raw
colnames(edge_quant)[4] = 'CR_prop_single_raw'
colnames(edge_quant)[5] = 'CR_p_single_raw'
edge_quant$CR_prop_single_norm = g_tbl_single_norm$proportion
edge_quant$CR_prop_full_raw = g_tbl_full_raw$proportion
edge_quant$CR_prop_full_norm = g_tbl_full_norm$proportion
edge_quant = edge_quant[,c(1:4,6:ncol(edge_quant),5)]
edge_quant$CR_p_single_norm = g_tbl_single_norm$p_value
edge_quant$CR_p_full_raw = g_tbl_full_raw$p_value
edge_quant$CR_p_full_norm = g_tbl_full_norm$p_value
write.csv(edge_quant,'figures/edge_direction_test.csv')


# alternatively, the p-value can be assessed with an alternative hypothesis of whether an interaction is enriched among all interactions in the perturbation of a core function
# this case is assessed in the file "3_DEG_modeling_edge_quantification_randFBA.R" and 

# combine the two ways of edge assessment and make the final figures of edge quantification
edge_prop = read.csv('figures/edge_quantification.csv') # this is the result when considering the proportion of an interaction among all interactions of a core function perturbation
edge_direct = read.csv('figures/edge_direction_test.csv') # this is the result when considering the proportion of an direction (up vs down) among the two direction of the interactions between a core function perturbation and DEG response related to another core function
edge_info = cbind(edge_prop, edge_direct)
# reorder to facilitate inspection 
edge_info = edge_info[,c(2,3,4,
                         5,17,6,18,7,19,8,20,
                         9,21, 10,22,11,23,12,24)]
write.csv(edge_info,'figures/edge_test_info.csv')

# visualization

# first to see how robust it is among all various combinations of setups
edge_info = edge_info[order(edge_info$RNAi), ]
library(pheatmap)
pheatmap(-log10(edge_info[,12:ncol(edge_info)]),breaks = seq(0,3,length.out = 100),
         cluster_rows = F, cluster_cols = F, labels_row = paste(edge_info$RNAi, '->',edge_info$DEG),annotation_row = edge_info[,"type", drop = F])

pheatmap(-log10(edge_info[,12:ncol(edge_info)]),breaks = seq(0,3,length.out = 100),
         cluster_rows = T, cluster_cols = T, labels_row = paste(edge_info$RNAi, '->',edge_info$DEG),annotation_row = edge_info[,"type", drop = F])


# these settings include the following options: 
# how to quantify significant: by overrepresentation of a direction of interaction or overrepresentation of an specific interaction 
# what to use: normalized (per condition proportion) or raw (total proportion)
# what data to use: single-core-function genes only (single) or all genes with core functions (full)
# now visualize every combinations 
return_graph <- function(filter_by, display_by, sig_cutoff){
  # decide filter by
  edge_info$best_p = pmin(edge_info[,filter_by[1]], edge_info[,filter_by[2]])
  filtered_g = edge_info[edge_info$best_p < sig_cutoff,]
  # decide display by
  filtered_g = filtered_g[,c("RNAi", "DEG", "type", display_by,"best_p")]
  colnames(filtered_g) = c('from','to','type','width','p_value')
  graph <- graph_from_data_frame(filtered_g, directed = TRUE)
  # add a label attribute
  E(graph)$label <- paste("Width:", signif(E(graph)$width,2), "P-value:", signif(E(graph)$p_value,2))
  # set the color based on -log10pvalue 
  signed_pvalues = -log10(E(graph)$p_value) * ifelse(E(graph)$type == 'up', 1, -1)
  mycolors = "white"
  colorScale = max(abs(signed_pvalues))
  
  # assign positive values
  colors <- colorRampPalette(c("white", "red"))
  # Generate the color vector
  mycolors[signed_pvalues>0] <- colors(100)[findInterval(signed_pvalues[signed_pvalues>0], 
                                           seq(0, colorScale, 
                                           length = 100))]
  
  # assign negative values
  colors <- colorRampPalette(c("blue", "white"))
  # Generate the color vector
  mycolors[signed_pvalues<0] <- colors(100)[findInterval(signed_pvalues[signed_pvalues<0], 
                                                                        seq(-colorScale,0, 
                                                                            length = 100))]
  
  # Drawing the graph
  l <- layout_in_circle(graph)
  p = plot(graph, layout = l,
       edge.arrow.size = E(graph)$width / max(E(graph)$width) * 1.5,
       edge.arrow.width = E(graph)$width / max(E(graph)$width) * 1.5,
       edge.label = E(graph)$label,
       edge.label.cex = 0.8,
       edge.width = E(graph)$width / max(E(graph)$width) * 15,
       edge.arrow.mode = '->',
       edge.color = mycolors,
       vertex.color = "gold", 
       vertex.size = 20,
       vertex.frame.color = "gray", 
       vertex.label.color = "black", 
       vertex.label.cex = 0.8, 
       vertex.label.dist = 1, 
       edge.curved=0.5,
       main = paste('filtered by:', paste(filter_by,collapse = ' '),'\n',
                    'displayed by:',display_by,'\n',
                    "edge proprotion or edge direction p <", sig_cutoff,'\n',
                    'color scale is', colorScale))
  return(p)
}


pdf('figures/fitted_CR_interaction_diagram.pdf',width = 10)
# the abbreviations refer to 
# filter by: whether or not to display an edge (interaction) is filtered by pvalues calculated in this setup
# display by: the thickness of the edge is displayed with values calculated by this setup
# prop_p: p-values calculated with the proportion of an edge among all edges for a core function perturbation 
# full: all genes with core function associations 
# single: only genes with single core function associations
# CR_p: p-values calculated with the proportion of edge direction (up or down, compensation or repression) between two core functions 
# norm: normalized proportion (first calculate proportion within each condition and then get average proportion)
# prop: proportion of an edge among all edges for a core function perturbation 
# raw: raw proportion (aggregate DEGs in all conditions and then calculate the proportion of one type of interaction)
sig_cutoff = 0.11
filter_by = c('prop_p_full_norm','CR_p_full_norm')
display_by = 'prop_full_norm'
return_graph(filter_by, display_by,sig_cutoff)

sig_cutoff = 0.1
filter_by = c('prop_p_single_norm','CR_p_single_norm')
display_by = 'prop_single_norm'
return_graph(filter_by, display_by,sig_cutoff)

sig_cutoff = 0.1
filter_by = c('prop_p_full_raw','CR_p_full_raw')
display_by = 'prop_full_raw'
return_graph(filter_by, display_by,sig_cutoff)

sig_cutoff = 0.1
filter_by = c('prop_p_single_raw','CR_p_single_raw')
display_by = 'prop_single_raw'
return_graph(filter_by, display_by,sig_cutoff)

dev.off()


# assess the robustness of interpretation given different combinations of setups
# save p-value and proportion heatmaps with selected edge labeled (to show the qualitative consistency)
pdf('figures/fitted_CR_parameter_robustness.pdf',width = 10)
pheatmap(-log10(edge_info[,12:ncol(edge_info)]),breaks = seq(0,3,length.out = 100),
         cluster_rows = T, cluster_cols = T, labels_row = paste(edge_info$RNAi, '->',edge_info$DEG),
         annotation_row = edge_info[,"type", drop = F])

pheatmap(edge_info[,c(4,6,8,10)],#breaks = seq(0,3,length.out = 100),
         cluster_rows = T, cluster_cols = T, labels_row = paste(edge_info$RNAi, '->',edge_info$DEG),
         annotation_row = edge_info[,"type", drop = F])

pheatmap(edge_info[,c(5,7,9,11)],#breaks = seq(0,3,length.out = 100),
         cluster_rows = T, cluster_cols = T, labels_row = paste(edge_info$RNAi, '->',edge_info$DEG),
         annotation_row = edge_info[,"type", drop = F])
dev.off()


# check the special case of ECM
edge_info[edge_info$RNAi=='pro_modi' & edge_info$DEG=='pro_modi' & edge_info$type=='up','CR_p_full_norm']
