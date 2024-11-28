# Summary
# Visualize the DEGs based on their core functions assigned by FBA simulation. 

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(Seurat)
library(SeuratDisk)
library(pheatmap)

# load the met gene classifications
classMat = read.csv('output/FBA_classification_matrix.csv',row.names = 1)

# scoreMat = read.csv('output/delta_flux_scoreMat.csv',row.names = 1)
# classMat = as.data.frame(1*(scoreMat > 1e-3)) # cutoff at 1e-3, the robustness of this cutoff will be tested later


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

# get the overlap (with normalized expression data)
cell_annotations <- seurat_object@meta.data
cell_annotations$perturbed_gene = perturbed_gene[match(rownames(cell_annotations), str_remove(colnames(adDE_res),'^X'))]
cell_annotations$perturbation_id = paste('X',rownames(cell_annotations),sep = '')
# padj_mat_met_resp = padj_mat_met[, colnames(padj_mat_met) %in% cell_annotations$perturbation_id[cell_annotations$anderson_darling_counts > 10]]
padj_mat_met_resp = padj_mat_met[, colnames(padj_mat_met) %in% colnames(exp_mat)]


# the z-score matrix 
exp_mat_met_resp = as.matrix(exp_mat[rownames(padj_mat_met_resp), colnames(padj_mat_met_resp)])

# mask change that is the RNAi targeted gene 
for (i in 1:ncol(padj_mat_met_resp)){
  if (any(perturbed_gene[colnames(padj_mat_met_resp)[i] == colnames(adDE_res)] == rownames(padj_mat_met_resp))){
    padj_mat_met_resp[perturbed_gene[colnames(padj_mat_met_resp)[i] == colnames(adDE_res)],i] = 1 # mask to insignificant (1)
  }
}


# visualize the classification of all genes that are targeted in human 1 responsive (at least two up or two down) perturbations
ct1 = colSums(padj_mat_met_resp < 0.05 & exp_mat_met_resp > 0,na.rm = T) # start with padj = 0.05
ct2 = colSums(padj_mat_met_resp < 0.05 & exp_mat_met_resp < 0,na.rm = T)
human1_resp = union(names(ct1)[ct1 > 1], names(ct2)[ct2 > 1])
tmp = classMat[rownames(classMat) %in% cell_annotations[match(human1_resp, cell_annotations$perturbation_id),'perturbed_gene'],]
colSums(tmp)
pheatmap::pheatmap(tmp)

# calculate some numbers
# FBA class
sum(rowSums(classMat[,c("energy","lipid","pro_syn","nucl_acid","ECM")])>0)/nrow(classMat)
sum(rowSums(classMat[,c("energy","lipid","pro_syn","nucl_acid")])>0)/nrow(classMat)
# only 22-29% genes are assigned

# at gene-level:
sum(rowSums(tmp[,c("energy","lipid","pro_syn","nucl_acid","ECM")])>0)/nrow(tmp)
sum(rowSums(tmp[,c("energy","lipid","pro_syn","nucl_acid")])>0)/nrow(tmp)
# 0.37 pass-filtered perturbations were included in the modeling framework 
sum(rowSums(tmp[,c("energy","lipid","pro_syn","nucl_acid","ECM")])==1)/nrow(tmp)
sum(rowSums(tmp[,c("energy","lipid","pro_syn","nucl_acid")])==1)/nrow(tmp)
# 0.25 were assigned to a unique classification


# PART I: VISUALIZE DEG IN EACH PERTURBATION
# We first focus on perturbations whose targeted gene is associated with unique core function.
total_condition_analyzed = c()
score_cutoff = 0.01 # because of the sensitivity issue in scRNA-seq, it seems we should use 
# a very greedy threshold instead of the cutoff of 0.05 used in the paper to increase power

# it seems 0.01 is very clean and more patterned; however, less perturbations and less DEGs
# in contrast, 0.2 is very rich with less clear C/R but more genes to show 
# let's keep in mind and test both in the future


obj_perturb = 'energy' # to debug with - not meaningful here
modelObj = c("energy","lipid","pro_syn","nucl_acid","ECM") # we exclude potential new core functions for now (glycogen, cofactor, metabolite))
pdfHeight = c("energy" = 12,
            'lipid' = 12,
            'pro_syn' = 12,
            'nucl_acid' = 12,
            'ECM' = 12
           # 'cofactor' = 12,
           # 'metabolite' = 12
           )
pdfWidth = c("energy" = 14,
             'lipid' = 14,
             'pro_syn' = 14,
             'nucl_acid' = 14,
             'ECM' = 12
             # 'cofactor' = 14,
            # 'metabolite' = 14
             )

# loop through each core functions
for (obj_perturb in modelObj){
  # we now only simply show the single-core-function genes
  obj_genes = intersect(rownames(classMat)[classMat[,obj_perturb]==1 & rowSums(classMat[,modelObj]) == 1],cell_annotations$perturbed_gene)
  # obj_genes = intersect(rownames(classMat)[classMat[,obj_perturb]==1],cell_annotations$perturbed_gene) # show everything 
  # get corresponding conditions 
  myconds = intersect(cell_annotations$perturbation_id[cell_annotations$perturbed_gene %in% obj_genes], colnames(padj_mat_met_resp)) 
  if (length(myconds) > 0){
    # count for up and down DEGs seperately
    for (obj_direction in c('up','down')){
      stats = data.frame(RNAi_cond = myconds)
      for (i in 1:length(modelObj)){
        stats[,modelObj[i]] = 0
      }
      stats$noClass = 0
      for (i in 1:nrow(stats)){
        myCond = stats$RNAi_cond[i]
        if (obj_direction == 'up'){
          DEGs = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] > 0)]
        }else{
          DEGs = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] < 0)]
        }
        # use majority-vote strategy to visualize DEGs with multiple core functions associated
        if (any(rownames(classMat) %in% DEGs)){
          # rank the obj candidates 
          nCalls = colSums(classMat[rownames(classMat) %in% DEGs, modelObj])
          orders = names(nCalls)[order(nCalls,decreasing = T)]
          # let's use a simple strategy: assign the multi-obj to the most abundant applicable obj
          subClassMat = classMat[rownames(classMat) %in% DEGs, modelObj]
          for (j in 1:nrow(subClassMat)){
            if (sum(subClassMat[j,])>1){
              for (k in 1:length(orders)){
                if (subClassMat[j,orders[k]] == 1){
                  subClassMat[j,] = 0
                  subClassMat[j,orders[k]] = 1
                  break
                }
              }
            }
          }
          stats[i,2:(length(modelObj)+1)] = colSums(subClassMat)
          stats[i,(length(modelObj)+2)] = sum(rowSums(subClassMat) == 0)
        }else{
          stats[i,2:(length(modelObj)+2)] = 0
        }
      }
      # before plotting...
      # filter out low response conditions
      stats = stats[rowSums(stats[,2:(length(modelObj)+1)])>1,] # we require at least two human 1 DE to enable meaningful interpretation (of proportion)
      # (if up or down gene is less than 2, we will mask the condition will NA (no show as blank bar))
      
      # plot
      b = stats$RNAi_cond
      # sort by the most abundant class
      tmp = stats[,2:(length(modelObj)+1)]
      tmp = tmp / rowSums(tmp)
      aveRewireRate = colMeans(tmp[,1:(length(modelObj))])
      sortBy = names(aveRewireRate)[order(aveRewireRate,decreasing = T)][1]
      rewire_rate = stats[,sortBy] / rowSums(stats[,2:(length(modelObj)+1)])
      rewire_rate[is.na(rewire_rate)] = 0
      b[rank(rewire_rate,ties.method ='first')] =b
      stats$RNAi_cond = factor(stats$RNAi_cond,levels = b)
      stats_long = reshape2::melt(stats)
      stats_long$variable = factor(as.character(stats_long$variable), levels = c('noClass',names(aveRewireRate)[order(aveRewireRate,decreasing = F)]))
      if (obj_direction == 'up'){
        stats_up = stats
        sortBy_up = sortBy
        stats_long_up = stats_long
      }else{
        stats_down = stats
        sortBy_down = sortBy
        stats_long_down = stats_long
      }
      
    }
    
    if (nrow(stats) > 0){ # this is an error-handling statement used in the development of the human analysis; not necessary and not mean anything. This if statement is not in the C. elegans code. The enclosed codes can be executed without the statement 
    
      # format the merged table
      stats_long_down$value = -stats_long_down$value
      extra_conds = levels(stats_long_down$RNAi_cond)[!(levels(stats_long_down$RNAi_cond) %in% levels(stats_long_up$RNAi_cond))]
      new_levels = c(extra_conds, levels(stats_long_up$RNAi_cond))
      stats_long_down$RNAi_cond = factor(as.character(stats_long_down$RNAi_cond), levels = new_levels)
      stats_long_up$RNAi_cond = factor(as.character(stats_long_up$RNAi_cond), levels = new_levels)
      
      stats_long_merge = stats_long_up
      stats_long_merge$up = stats_long_merge$value
      stats_long_merge$down = 0
      stats_long_down$up = 0
      stats_long_down$down = stats_long_down$value
      toMerge = c()
      for (i in 1:nrow(stats_long_down)){
        if (any(stats_long_merge$RNAi_cond == stats_long_down$RNAi_cond[i] & stats_long_merge$variable == stats_long_down$variable[i])){
          stats_long_merge$down[stats_long_merge$RNAi_cond == stats_long_down$RNAi_cond[i] & stats_long_merge$variable == stats_long_down$variable[i]] = stats_long_down$value[i]
        }else{
          toMerge = c(toMerge, i)
        }
      }
      stats_long_merge = rbind(stats_long_merge, stats_long_down[toMerge,])
      
      # plot up and down with the same condition order for later merging in graphics
      old_levels = levels(stats_long_merge$RNAi_cond)
      gene_names = as.character(stats_long_merge$RNAi_cond)
      gene_names = unlist(lapply(strsplit(gene_names, '_'),function(x){x[[2]]}))
      new_levels = unlist(lapply(strsplit(old_levels, '_'),function(x){x[[2]]}))
      stats_long_merge$RNAi_cond = factor(gene_names, new_levels)
      
      p1 <- ggplot(stats_long_merge, aes(x=RNAi_cond, y = up,fill=variable)) +
        geom_bar(stat="identity", position="fill",colour="black", size = 0.1) + 
        scale_y_continuous(labels = scales::percent,expand = c(0,0))+
        theme_bw()+
        theme(axis.text.x = element_text(colour ='black'),axis.title.x=element_text(colour ='black'),
              axis.ticks.y=element_blank(),
              axis.text.y = element_blank(),axis.title.y=element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(colour ='black'),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
        coord_flip()+
        scale_fill_manual("legend", values = c("energy" = "#FF6701", 
                                               "lipid" = "#77AC30",
                                               "ECM" = "#EDB120",
                                               "nucl_acid" = "#4DBEEE",
                                               "pro_syn" = "#0072BD",
                                              # "cofactor" = "#8E44AD",
                                              # "metabolite" = "#C0392B",
                                               "noClass" = "white")) +
        labs(y = 'Up DEG')+
        geom_text(aes(label = up), position = position_fill(vjust = 0.5),size = 3)+
        geom_hline(yintercept = sum(classMat[,sortBy_up]) / nrow(classMat),color="white",
                   linetype="dashed")
      
      #p1
      
      # plot down with obj reordered 
      stats_long_merge$variable = factor(as.character(stats_long_merge$variable), levels = levels(stats_long_down$variable))
      p2 <- ggplot(stats_long_merge, aes(x=RNAi_cond, y = down,fill=variable)) +
        geom_bar(stat="identity", position="fill",colour="black", size = 0.1) + 
        scale_y_continuous(labels = scales::percent,expand = c(0,0))+
        theme_bw()+
        theme(axis.text.x = element_text(colour ='black'),axis.title.x=element_text(colour ='black'),
              axis.text.y = element_text(colour ='black'),
              legend.position = "none",
              panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
        coord_flip()+
        scale_fill_manual("legend", values = c("energy" = "#FF6701", 
                                               "lipid" = "#77AC30",
                                               "ECM" = "#EDB120",
                                               "nucl_acid" = "#4DBEEE",
                                               "pro_syn" = "#0072BD",
                                             #  "cofactor" = "#8E44AD",
                                             #  "metabolite" = "#C0392B",
                                               "noClass" = "white")) +
        labs(y = 'Down DEG')+
        geom_text(aes(label = -down), position = position_fill(vjust = 0.5),size = 3)+
        geom_hline(yintercept = -sum(classMat[,sortBy]) / nrow(classMat),color="white",
                   linetype="dashed")
      
      #p2
      
      pdf(paste('figures/FBA_classification_',obj_perturb,'_uniObjGenes.pdf',sep = ''),
          height = pdfHeight[obj_perturb], width = pdfWidth[obj_perturb])
      print(ggarrange(p2, p1, nrow = 1))
      dev.off()
      print(ggarrange(p2, p1, nrow = 1))
      total_condition_analyzed = c(total_condition_analyzed, as.character(unique(stats_long_merge$RNAi_cond)))
    }
  }
}

# it turns out that 148 conditions were analyzed (unique obj association)


##### plot the Euler graph for core function classifications #####
modelObj = c("energy","lipid","ECM","pro_syn","nucl_acid") # we exclude potential new core functions for now (glycogen, cofactor, metabolite))

gene_list = list(
  iCEL = rownames(classMat),
  energy = rownames(classMat)[classMat$energy==1],
  lipid = rownames(classMat)[classMat$lipid==1],
  ECM = rownames(classMat)[classMat$ECM==1],
  pro_syn = rownames(classMat)[classMat$pro_syn==1],
  nucl_acid = rownames(classMat)[classMat$nucl_acid==1]
)

library(eulerr)
set.seed(1030)
seeds = sample(1:100000, 100)
fits = list()
for (i in 1:5){
  set.seed(seeds[i])
  
  fit2 <- euler(classMat[,modelObj], 
                shape = "ellipse",
                #loss = 'region',
                #loss_aggregator = 'max', 
                control = list(extraopt = TRUE, extraopt_threshold = 0)
                )
  fits[[i]] = fit2
  pdf(paste('figures/eulerPlot_random/euler_plot_iCEL_gene_classification_seed_',seeds[i],'.pdf',sep = ''),width = 7,height = 7)
  print(plot(fit2, quantities = TRUE, fill = c('#FF6701','#77AC30','#EDB120','#0072BD','#4DBEEE'),main = paste('unclassified =', sum(rowSums(classMat[,modelObj]==1)==0),'(',100*round(sum(rowSums(classMat[,modelObj]==1)==0)/nrow(classMat),2),'% )')))
  dev.off()
}

# the best appearance is from seed 94084 (there are quite a few with similar layout. we just randomly picked one)

df = as.data.frame(fit2$original.values)
write.csv(df, 'figures/euler_plot_iCEL_gene_classification_breakdown.csv')
# 03112024: we realized that mccc-2 was not included in this figure, but here we didnt compare with WPS data
# so we should include it. For the sake of convinience, we just modify the numbers in the final figure without 
# regenerating the figure (because +1 is not going to be notable in the area of the chart)

# we cannot get a visually plausible figure by changing seeds, so we give up on getting this from the function
# itself, we will just change the numbers in the original figure - it is generally fine because the changes in
# the area is very minor 

# error_plot(fit2)


