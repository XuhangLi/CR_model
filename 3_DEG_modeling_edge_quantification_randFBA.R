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

# load data 
# load the met gene classifications
classMat = read.csv('output/FBA_classification_matrix.csv',row.names = 1)
iCELnames = read.csv('input/iCEL_IDtbl.csv')
rownames(classMat) = iCELnames$WormBase_Gene_ID[match(rownames(classMat),iCELnames$ICELgene)]
# fix minor bug - mrpl-44 is a misannotated gene in the model; in the model it should be mccc-2 but mislabeled as mrpl-44 (which is mito ribo gene), since ribo genes are excluded from the analysis, we also exclude mrpl-44
classMat = classMat[-which(rownames(classMat) == 'WBGene00008514'),]

# load the DE similarity 
DEsim = read.csv('input/cosineSimilarity_FC_denoised_stationery_metabolic.csv', row.names = 1)
# remove the mrpl-44 condition from the analysis (we remove both from DEG and conditions since this DEG is not icel either)
DEsim = DEsim[-which(rownames(DEsim) == 'x.mrpl_44_met3_lib1'),]
DEsim = DEsim[,-which(colnames(DEsim) == 'x.mrpl_44_met3_lib1')]

# load DEG result
inputTb=read.csv('input/DE_merged_clean_pcutoff_0.005_master_table_FDR2d0.1_FC1.5_FCtype_log2FoldChange_raw_ALL.csv')
# remove the mrpl-44 condition from the analysis (we remove both from DEG and conditions since this DEG is not icel either)
inputTb = inputTb[-which(inputTb$RNAi == 'x.mrpl_44'),]
inputTb = inputTb[inputTb$WBID != 'WBGene00008514',]

# no self
inputTb$RNAiID = paste(inputTb$RNAi, inputTb$batchID)
inputTb$RNAi_geneName = inputTb$RNAi
inputTb$RNAi_geneName = str_replace(inputTb$RNAi_geneName,'^x.','')
inputTb$RNAi_geneName = str_replace_all(inputTb$RNAi_geneName,'_','-')
inputTb$RNAi_geneName = str_replace(inputTb$RNAi_geneName,'-L[0-9]$','')
inputTb = inputTb[-which(inputTb$RNAi_geneName == inputTb$Gene_name),]
inputTb$condID = paste(inputTb$RNAi,inputTb$batchID,sep = '_')
if (length(which(inputTb$WBID=="NoHit"))>0){
  inputTb=inputTb[-which(inputTb$WBID=="NoHit"),]}
# filtering the low responsive and non-metablic ones
conditionInfo = read.csv('input/RNAi_condition_metaInfo.csv',row.names = 1)
# only analyze the iCEL responsive
inputTb_metResponsiove = inputTb[inputTb$RNAiID %in% conditionInfo$RNAiID[conditionInfo$isICEL & conditionInfo$isResponsive], ]
inputTb_metResponsiove=inputTb_metResponsiove[,c("WBID","RNAi","log2FoldChange_raw","condID")]

# exclude the nonclassic metabolic genes from the analysis 
# only keep iCEL DEG
inputTb_metResponsiove = inputTb_metResponsiove[inputTb_metResponsiove$WBID %in% rownames(classMat),]
# exclude UGT (both for RNAi and DEG)
ugtlist = read.csv('input/ugt_in_model.csv')
inputTb_metResponsiove = inputTb_metResponsiove[!(inputTb_metResponsiove$WBID %in% 
                                                    iCELnames$WormBase_Gene_ID[match(ugtlist$ugt_in_model, iCELnames$ICELgene)]),]
inputTb_metResponsiove = inputTb_metResponsiove[!(conditionInfo$RNAi_WBID[match(inputTb_metResponsiove$condID, str_replace(conditionInfo$RNAiID,' ','_'))] %in% 
                                                    iCELnames$WormBase_Gene_ID[match(ugtlist$ugt_in_model, iCELnames$ICELgene)]),]
# exclude vha
vhalist = read.csv('input/vha_in_model.csv')
inputTb_metResponsiove = inputTb_metResponsiove[!(inputTb_metResponsiove$WBID %in% 
                                                    iCELnames$WormBase_Gene_ID[match(vhalist$vha_in_model, iCELnames$ICELgene)]),]
inputTb_metResponsiove = inputTb_metResponsiove[!(conditionInfo$RNAi_WBID[match(inputTb_metResponsiove$condID, str_replace(conditionInfo$RNAiID,' ','_'))] %in% 
                                                    iCELnames$WormBase_Gene_ID[match(vhalist$vha_in_model, iCELnames$ICELgene)]),]
# exclude RNA pol and DNA pol (eq. ribosome that is not in the model)
pollist = read.csv('input/pol_in_model.csv')
inputTb_metResponsiove = inputTb_metResponsiove[!(inputTb_metResponsiove$WBID %in% 
                                                    iCELnames$WormBase_Gene_ID[match(pollist$pol_in_model, iCELnames$ICELgene)]),]
inputTb_metResponsiove = inputTb_metResponsiove[!(conditionInfo$RNAi_WBID[match(inputTb_metResponsiove$condID, str_replace(conditionInfo$RNAiID,' ','_'))] %in% 
                                                    iCELnames$WormBase_Gene_ID[match(pollist$pol_in_model, iCELnames$ICELgene)]),]


# show the classification of all genes that is iCEL responsive (at least two up or two down)
upTbl = inputTb_metResponsiove[inputTb_metResponsiove$log2FoldChange_raw > 0, ]
ct1 = table(upTbl$condID)
downTbl = inputTb_metResponsiove[inputTb_metResponsiove$log2FoldChange_raw < 0, ]
ct2 = table(downTbl$condID)
icel_resp = union(names(ct1)[ct1 > 1], names(ct2)[ct2 > 1])
tmp = classMat[rownames(classMat) %in% conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% icel_resp],]
colSums(tmp)


sum(rowSums(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')])>0)/nrow(tmp)
sum(rowSums(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')])==1)/nrow(tmp)
# total number of analyzable condition is 
n_total = sum(rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp, str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0 )
n_total/length(icel_resp)
# check for the coverage of the unclustered conditions
uniConds = icel_resp[rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp, str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0]
#RNAiclusters = read.csv('./../../2_DE/output/RNAi_groups.csv')
#sum(RNAiclusters$clusters[(str_replace(RNAiclusters$RNAiID,' ','_') %in% uniConds)] == -1)/sum(RNAiclusters$clusters==-1)
length(intersect(rownames(classMat)[rowSums(classMat[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0], inputTb_metResponsiove$WBID))/(length(unique(inputTb_metResponsiove$WBID)))


# we get the conditions to analyze (icel_responsive (at least 2 up or 2 down DEG) and classified)
total_condition_analyzed = c()
for (i in 1:length(icel_resp)){
  if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
    total_condition_analyzed = c(total_condition_analyzed, icel_resp[i])
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
    DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
    DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
    
    # determine affected obj
    labeledObj = colnames(classMat)[classMat[conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') == myCond],]==1]
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
mustHasClass = unique(conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% total_condition_analyzed])
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
    DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
    DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
    
    # determine compromised obj
    labeledObj = colnames(classMat_rand)[classMat_rand[conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') == myCond],]==1]
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


#############
# # some supplementary figures
# gene_list = list(
#   iCEL = rownames(classMat),
#   energy = rownames(classMat)[classMat$energy==1],
#   lipid = rownames(classMat)[classMat$lipid==1],
#   pro_modi = rownames(classMat)[classMat$pro_modi==1],
#   pro_syn = rownames(classMat)[classMat$pro_syn==1],
#   nucl_acid = rownames(classMat)[classMat$nucl_acid==1]
# )
# 
# library(eulerr)
# fit2 <- euler(classMat[,modelObj], 
#               shape = "ellipse",
#               #loss = 'region',
#               #loss_aggregator = 'max', 
#               control = list(extraopt = TRUE, extraopt_threshold = 0)
#               )
# pdf(paste('figures/euler_plot_iCEL_gene_classification.pdf',sep = ''),width = 7,height = 7)
# plot(fit2, quantities = TRUE, fill = c('#FF6701','#77AC30','#EDB120','#0072BD','#4DBEEE'),main = paste('unclassified =', sum(rowSums(classMat[,modelObj]==1)==0),'(',100*round(sum(rowSums(classMat[,modelObj]==1)==0)/nrow(classMat),2),'% )'))
# dev.off()
# 
# # error_plot(fit2)
# 
# 
# 
# 
# 
# # given any RNAi node to any DEG node has only one or less regulation pattern, the default CR model describes
# # one of all possible models; it explains 50% of DEG based on independent counting strategy with raw count and no multi-class scaling 
# # (it counts as true as long as if any catgeory combination fits into expectation therefore it is similar to raw count); in which case
# # using the major direction (up or down) of each possible interaction will give a best fit. However, this model is overfiting and less meaningful
# # we alternatively derive a best decriptive model based on manual inspection of interaction summary; we down-weight the multi-class
# # to gain the best interpretation instead of best fitting (because the explainability criteria in multi-obj is greedy) (but these two model might
# # be the same)
# 
# 
# plot(unlist(model_fit_single$raw_counts$energy),
# unlist(model_fit_single$normalized_counts$energy))
# 
# unlist(model_fit_single$raw_counts$lipid$up)
# unlist(model_fit_single$normalized_counts$lipid$up)
# 
# unlist(model_fit_single$raw_counts$pro_modi$up)
# unlist(model_fit_single$normalized_counts$pro_modi$up)
# 
# unlist(model_fit_single$raw_counts$pro_syn$up)
# unlist(model_fit_single$normalized_counts$pro_syn$up)
# 
# unlist(model_fit_single$raw_counts$nucl_acid$up)
# unlist(model_fit_single$normalized_counts$nucl_acid$up)
# 
# sum(unlist(model_fit_single$normalized_counts$pro_modi$up))
# 
# barplot(model_fit_single$raw_counts$energy)
# 
# # 
# # unlist(model_fit_single$raw_counts$energy$down)
# # unlist(model_fit_single$normalized_counts$energy$down)
# # 
# # unlist(model_fit_single$raw_counts$lipid$down)
# # unlist(model_fit_single$normalized_counts$lipid$down)
# # 
# # unlist(model_fit_single$raw_counts$pro_modi$down)
# # unlist(model_fit_single$normalized_counts$pro_modi$down)
# # 
# # unlist(model_fit_single$raw_counts$pro_syn$down)
# # unlist(model_fit_single$normalized_counts$pro_syn$down)
# # 
# # unlist(model_fit_single$raw_counts$nucl_acid$down)
# # unlist(model_fit_single$normalized_counts$nucl_acid$down)
# # 
# # 
# # plot(unlist(model_fit_single$raw_counts$energy$up), 
# #      unlist(model_fit_all$raw_counts$energy$up))
# # plot(unlist(model_fit_single$raw_counts$energy$down), 
# #      unlist(model_fit_all$raw_counts$energy$down))
# # 
# # plot(unlist(model_fit_single$raw_counts$lipid$up), 
# #      unlist(model_fit_all$raw_counts$lipid$up))
# # plot(unlist(model_fit_single$raw_counts$lipid$down), 
# #      unlist(model_fit_all$raw_counts$lipid$down))
# # 
# # plot(unlist(model_fit_single$raw_counts$pro_modi$up), 
# #      unlist(model_fit_all$raw_counts$pro_modi$up))
# # plot(unlist(model_fit_single$raw_counts$pro_modi$down), 
# #      unlist(model_fit_all$raw_counts$pro_modi$down))
# # 
# # plot(unlist(model_fit_single$raw_counts$nucl_acid$up), 
# #      unlist(model_fit_all$raw_counts$nucl_acid$up))
# # plot(unlist(model_fit_single$raw_counts$nucl_acid$down), 
# #      unlist(model_fit_all$raw_counts$nucl_acid$down))
# # 
# # plot(unlist(model_fit_single$raw_counts$pro_syn$up), 
# #      unlist(model_fit_all$raw_counts$pro_syn$up))
# # plot(unlist(model_fit_single$raw_counts$pro_syn$down), 
# #      unlist(model_fit_all$raw_counts$pro_syn$down))
# # 
# # 
# # 
# # 
# # 
# # unlist(model_fit_all$raw_counts$energy$up)
# # unlist(model_fit_all$normalized_counts$energy$up)
# # 
# # unlist(model_fit_all$raw_counts$lipid$up)
# # unlist(model_fit_all$normalized_counts$lipid$up)
# # 
# # unlist(model_fit_all$raw_counts$pro_modi$up)
# # unlist(model_fit_all$normalized_counts$pro_modi$up)
# # 
# # unlist(model_fit_all$raw_counts$pro_syn$up)
# # unlist(model_fit_all$normalized_counts$pro_syn$up)
# # 
# # unlist(model_fit_all$raw_counts$nucl_acid$up)
# # unlist(model_fit_all$normalized_counts$nucl_acid$up)
# # 
# # sum(unlist(model_fit_all$normalized_counts$pro_modi$up))
# # 
# # 
# # unlist(model_fit_all$raw_counts$energy$down)
# # unlist(model_fit_all$normalized_counts$energy$down)
# # 
# # unlist(model_fit_all$raw_counts$lipid$down)
# # unlist(model_fit_all$normalized_counts$lipid$down)
# # 
# # unlist(model_fit_all$raw_counts$pro_modi$down)
# # unlist(model_fit_all$normalized_counts$pro_modi$down)
# # 
# # unlist(model_fit_all$raw_counts$pro_syn$down)
# # unlist(model_fit_all$normalized_counts$pro_syn$down)
# # 
# # unlist(model_fit_all$raw_counts$nucl_acid$down)
# # unlist(model_fit_all$normalized_counts$nucl_acid$down)
