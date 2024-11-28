# sensitivity test of the cutoff parameter in CR model 

# the CR model analysis contains a single parameter, the cutoff of the core function score. By default, we used 0.001 because this one is empirically
# learned by recalling most genes in ETC RNAi cluster for energy function and recalling most genes in lipid RNAi cluster for lipid core function. 

# To test the sensitivity of this parameter, we titrated this parameter across all possible range (from minimal score of 0 to maximal score of 1)

# since we use greater than cutoff to define associations, we tested cutoffs as follows: (0, 1e-4,1e-3, 1e-2, 0.1, 0.99)
# 0.99 representing the test of maximal score of 1.

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(Seurat)
library(SeuratDisk)
library(pheatmap)

FBA_cutoffs = c(0, 1e-6, 1e-5, 1e-4,1e-3, 1e-2, 0.1, 0.5, 0.99)

# set cutoff
score_cutoff = 0.01 # because of the sensitivity issue in scRNA-seq, it seems we should use 


pvalues = c()
total_DEG_explained = c()

# load DEG result
adDE_res = read.csv('input/source_data/processed_data/anderson-darling_p-values_BH-corrected.csv', row.names = 1)
seurat_object <- LoadH5Seurat("./input/source_data/processed_data/K562_gwps_normalized_bulk_01.h5seurat")

for (FBA_cutoff in FBA_cutoffs){
  # repeat the same code as in the FBA-fused CR model analysis, just with different threhsold 
  scoreMat = read.csv('output/delta_flux_scoreMat.csv',row.names = 1)
  colnames(scoreMat)[colnames(scoreMat) == 'ECM'] = 'pro_modi' # for the reuse of old code, just use the old name 
  classMat = as.data.frame(1*(scoreMat > FBA_cutoff))
  
  
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
  
  
  # visualize the classification of all genes that are targeted in human 1 responsive (at least two up or two down) perturbations
  ct1 = colSums(padj_mat_met_resp < score_cutoff & exp_mat_met_resp > 0,na.rm = T) # start with padj = 0.05
  ct2 = colSums(padj_mat_met_resp < score_cutoff & exp_mat_met_resp < 0,na.rm = T)
  human1_resp = union(names(ct1)[ct1 > 1], names(ct2)[ct2 > 1])
  tmp = classMat[rownames(classMat) %in% cell_annotations[match(human1_resp, cell_annotations$perturbation_id),'perturbed_gene'],]
  colSums(tmp)
  
  # calculate the proportion explained
  total_condition_analyzed = c()
  for (i in 1:length(human1_resp)){
    if (rowSums(classMat[cell_annotations$perturbed_gene[match(human1_resp[i], cell_annotations$perturbation_id)],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
      total_condition_analyzed = c(total_condition_analyzed, human1_resp[i])
    }
  }
  
  
  modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
  
  model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
  model_explained$UP_yes = 0
  model_explained$UP_no = 0
  model_explained$DOWN_yes = 0
  model_explained$DOWN_no = 0
  
  for (condInd in 1:length(total_condition_analyzed)){
      myCond = total_condition_analyzed[condInd]
      DEGs_up = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] > 0)]
      DEGs_down = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] < 0)]
      
      # determine compromised obj
      labeledObj = colnames(classMat)[classMat[cell_annotations$perturbed_gene[cell_annotations$perturbation_id == myCond],]==1]
      affectedObj = intersect(labeledObj, modelObj)
      
      # compare with DEG 
      # up genes
      subClassMat = classMat[rownames(classMat) %in% DEGs_up, modelObj]
      # calculate the DEG explained by model
      model_explained$UP_yes[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) > 0)
      model_explained$UP_no[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) == 0)
      # down genes
      subClassMat = classMat[rownames(classMat) %in% DEGs_down, modelObj]
      # calculate the DEG explained by model
      model_explained$DOWN_yes[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) > 0)
      model_explained$DOWN_no[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) == 0)
    
  }
  
  # plot 
  b = model_explained$RNAi_cond
  # sort by the total explained rate
  tmp = model_explained[,2:5]
  tmp = tmp / rowSums(tmp)
  rewire_rate = rowSums(tmp[,c(1,3)])
  rewire_rate[is.na(rewire_rate)] = 0
  b[rank(rewire_rate,ties.method ='first')] =b
  model_explained$RNAi_cond = factor(model_explained$RNAi_cond,levels = b)
  
  model_explained_long_up = reshape2::melt(model_explained[,c("RNAi_cond","UP_yes","UP_no")])
  model_explained_long_down = reshape2::melt(model_explained[,c("RNAi_cond","DOWN_yes","DOWN_no" )])
  model_explained_long = model_explained_long_up
  model_explained_long$down = model_explained_long_down$value
  colnames(model_explained_long) = c("RNAi_cond", "variable", "up", "down") 
  model_explained_long$variable = as.character(model_explained_long$variable)
  model_explained_long$variable[model_explained_long$variable == 'UP_yes'] = 'yes'
  model_explained_long$variable[model_explained_long$variable == 'UP_no'] = 'no'
  model_explained_long$variable = factor(model_explained_long$variable, levels = c('no','yes'))
  # calculate the average
  up_ave1 = mean(model_explained$UP_yes/(model_explained$UP_yes + model_explained$UP_no),na.rm = T) # mean of the rate
  down_ave1 = mean(model_explained$DOWN_yes/(model_explained$DOWN_yes + model_explained$DOWN_no),na.rm = T) # mean of the rate
  up_ave2 = sum(model_explained_long$up[model_explained_long$variable=='yes'])/sum(model_explained_long$up) # total rate
  down_ave2 = sum(model_explained_long$down[model_explained_long$variable=='yes'])/sum(model_explained_long$down) # total rate 
  # the total DEG explained is less sensitive to conditions that have very few DEG, so we use this as the reference line
  
  # plot up and down with the same condition order for later merging in graphics
  p1 <- ggplot(model_explained_long, aes(x=RNAi_cond, y = up,fill=variable)) +
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
    scale_fill_manual("legend", values = c("yes" = "#FF6701", 
                                           "no" = 'grey')) +
    labs(y = 'Up DEG')+
    geom_text(aes(label = up), position = position_fill(vjust = 0.5),size = 1)+
    geom_hline(yintercept = up_ave2,color="white",
               linetype="dashed")
  
  #p1
  
  # plot down 
  p2 <- ggplot(model_explained_long, aes(x=RNAi_cond, y = -down,fill=variable)) +
    geom_bar(stat="identity", position="fill",colour="black", size = 0.1) + 
    scale_y_continuous(labels = scales::percent,expand = c(0,0))+
    theme_bw()+
    theme(axis.text.x = element_text(colour ='black'),axis.title.x=element_text(colour ='black'),
          axis.text.y = element_text(colour ='black', size = 2),
          legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    coord_flip()+
    scale_fill_manual("legend", values = c("yes" = "#FF6701", 
                                           "no" = 'grey')) +
    labs(y = 'Down DEG')+
    geom_text(aes(label = -down), position = position_fill(vjust = 0.5),size = 1)+
    geom_hline(yintercept = -down_ave2,color="white",
               linetype="dashed")
  
  #p2
  
  pdf(paste('figures/sensitivity_analysis/FBA_DEG_explained_distribution_FBA_fused_model_',FBA_cutoff,'.pdf',sep = ''),
      height = 7, width = 14)
  print(ggarrange(p2, p1, nrow = 1))
  dev.off()
  
  
  # the overall average explained rate (up and down together) is 
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
  
  # randomization to assess the significance 
  nRand = 10000
  set.seed(1126)
  rand_rate = c()
  rand_rate2 = c()
  rand_rate_up = c()
  rand_rate2_up = c()
  rand_rate_down = c()
  rand_rate2_down = c()
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
    
    for (condInd in 1:length(total_condition_analyzed)){
      myCond = total_condition_analyzed[condInd]
      DEGs_up = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] > 0)]
      DEGs_down = rownames(exp_mat_met_resp)[which(padj_mat_met_resp[, myCond] < score_cutoff & exp_mat_met_resp[, myCond] < 0)]
      
      # determine compromised obj
      labeledObj = colnames(classMat_rand)[classMat_rand[cell_annotations$perturbed_gene[cell_annotations$perturbation_id == myCond],]==1]
      affectedObj = intersect(labeledObj, modelObj)
      
      # compare with DEG 
      # up genes
      subClassMat = classMat_rand[rownames(classMat_rand) %in% DEGs_up, modelObj]
      # calculate the DEG explained by model
      model_explained$UP_yes[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) > 0)
      model_explained$UP_no[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) == 0)
      # down genes
      subClassMat = classMat_rand[rownames(classMat_rand) %in% DEGs_down, modelObj]
      # calculate the DEG explained by model
      model_explained$DOWN_yes[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) > 0)
      model_explained$DOWN_no[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) == 0)
      
    }
    
    # total explained rate
    tmp = model_explained[,2:5]
    tmp = tmp / rowSums(tmp)
    rewire_rate_rand = rowSums(tmp[,c(1,3)])
    rewire_rate_rand[is.na(rewire_rate_rand)] = 0
    
    # the overall average explained rate (up and down together) is 
    rand_rate = c(rand_rate, mean(rewire_rate_rand))
    rand_rate2 = c(rand_rate2, sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5]))
    
    # up
    tmp = model_explained[,2:3]
    tmp = tmp / rowSums(tmp)
    rand_rate_up = c(rand_rate_up, mean(tmp$UP_yes, na.rm = T))
    rand_rate2_up = c(rand_rate2_up, sum(model_explained[,c(2)]) / sum(model_explained[,2:3]))
    # down
    tmp = model_explained[,4:5]
    tmp = tmp / rowSums(tmp)
    rand_rate_down = c(rand_rate_down, mean(tmp$DOWN_yes, na.rm = T))
    rand_rate2_down = c(rand_rate2_down, sum(model_explained[,c(4)]) / sum(model_explained[,4:5]))
    
    print(nn)
    # list(list(mean(rewire_rate_rand), sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5])))
  }
  
  # save results 
  save(file = paste('output/sensitivity_analysis/randomization_result_FBA_fused_model_cutoff_',FBA_cutoff,'.Rdata',sep = ''),rand_rate, rand_rate2, rand_rate_up, rand_rate2_up, rand_rate_down, rand_rate2_down)
  hist(rand_rate)
  hist(rand_rate2)
  pdf(paste('figures/sensitivity_analysis/overall_DEG_explained_randomization_FBA_fused_model_',FBA_cutoff,'.pdf',sep = ''),width = 7,height = 6)
  hist(rand_rate, main = paste('p <',signif((1+sum(rand_rate>=obs_rate))/(1+length(rand_rate)),2)),
       xlim = c(min(rand_rate)*0.9,max(c(rand_rate,obs_rate))*1.1),
       xlab = 'condition-wise average percentage of DEG explained')
  abline(v = obs_rate)
  
  hist(rand_rate2, main = paste('p <',signif((1+sum(rand_rate2>=obs_total_DE_rate))/(1+length(rand_rate2)),2)),
       xlim = c(min(rand_rate2)*0.9,max(c(rand_rate2,obs_total_DE_rate))*1.1),
       xlab = 'percentage of total DEG explained')
  abline(v = obs_total_DE_rate)
  
  hist(rand_rate_up, main = paste('p <',signif((1+sum(rand_rate_up>=obs_rate_up))/(1+length(rand_rate_up)),2)),
       xlim = c(min(rand_rate_up)*0.9,max(c(rand_rate_up,obs_rate_up))*1.1),
       xlab = 'condition-wise average percentage of up DEG explained')
  abline(v = obs_rate_up)
  
  hist(rand_rate2_up, main = paste('p <',signif((1+sum(rand_rate2_up>=obs_total_DE_rate_up))/(1+length(rand_rate2_up)),2)),
       xlim = c(min(rand_rate2_up)*0.9,max(c(rand_rate2_up,obs_total_DE_rate_up))*1.1),
       xlab = 'percentage of total UP DEG explained')
  abline(v = obs_total_DE_rate_up)
  
  hist(rand_rate_down, main = paste('p <',signif((1+sum(rand_rate_down>=obs_rate_down))/(1+length(rand_rate_down)),2)),
       xlim = c(min(rand_rate_down)*0.9,max(c(rand_rate_down,obs_rate_down))*1.1),
       xlab = 'condition-wise average percentage of down DEG explained')
  abline(v = obs_rate_down)
  
  hist(rand_rate2_down, main = paste('p <',signif((1+sum(rand_rate2_down>=obs_total_DE_rate_down))/(1+length(rand_rate2_down)),2)),
       xlim = c(min(rand_rate2_down)*0.9,max(c(rand_rate2_down,obs_total_DE_rate_down))*1.1),
       xlab = 'percentage of total down DEG explained')
  abline(v = obs_total_DE_rate_down)
  
  dev.off()
  
  # the two key metric: proprotion of genes explained and the significance
  pvalues = c(pvalues, (1+sum(rand_rate2>=obs_total_DE_rate))/(1+length(rand_rate2)))
  total_DEG_explained = c(total_DEG_explained,obs_total_DE_rate)
}

# plot the result
save(FBA_cutoffs,pvalues,total_DEG_explained, file = 'output/sensitivity_analysis/data_to_plot.Rdata')
pdf('figures/parameter_sensitivity_to_objFBA_cutoff.pdf',width = 7,height = 6)
plot(-log10(FBA_cutoffs+1e-7), -log10(pvalues),type = 'b')
plot(-log10(FBA_cutoffs+1e-7), total_DEG_explained,type = 'b')
dev.off()
