# 02062023
# note: these responsive conditions may be redefined with losen cutoff for meta-analysis 

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)

# load the met gene classifications
classMat = read.csv('output/FBA_classification_matrix.csv',row.names = 1)
iCELnames = read.csv('./../../input_data/otherTbls/iCEL_IDtbl.csv')
rownames(classMat) = iCELnames$WormBase_Gene_ID[match(rownames(classMat),iCELnames$ICELgene)]
# fix minor bug - mrpl-44 is a misannotated gene in the model; in the model it should be mccc-2 but mislabeled as mrpl-44 (which is mito ribo gene), since ribo genes are excluded from the analysis, we also exclude mrpl-44
classMat = classMat[-which(rownames(classMat) == 'WBGene00008514'),]

# load the DE similarity 
DEsim = read.csv('./../../2_DE/output/cosineSimilarity_FC_denoised_stationery_metabolic.csv', row.names = 1)
# remove the mrpl-44 condition from the analysis (we remove both from DEG and conditions since this DEG is not icel either)
DEsim = DEsim[-which(rownames(DEsim) == 'x.mrpl_44_met3_lib1'),]
DEsim = DEsim[,-which(colnames(DEsim) == 'x.mrpl_44_met3_lib1')]

# load DEG result
inputTb=read.csv('./../../2_DE/output/DE_merged_clean_pcutoff_0.005_master_table_FDR2d0.1_FC1.5_FCtype_log2FoldChange_raw_ALL.csv')
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
conditionInfo = read.csv('./../../2_DE/output/RNAi_condition_metaInfo.csv',row.names = 1)
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


# manually inspecting a few conditions
# myCond = 'x.nduf_2.2_met6_lib1' # use this to tune the program
# upDEGs = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
# downDEGs = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
# 
# tmp = classMat[rownames(classMat) %in% upDEGs, ]
# colSums(tmp)
# pheatmap::pheatmap(tmp, labels_row = iCELnames$ICELgene[match(rownames(tmp), iCELnames$WormBase_Gene_ID)])
# pdf('tmp.pdf')
# pheatmap::pheatmap(tmp, labels_row = iCELnames$ICELgene[match(rownames(tmp), iCELnames$WormBase_Gene_ID)],fontsize_row = 3)
# dev.off()
# 
# tmp = classMat[rownames(classMat) %in% downDEGs, ]
# colSums(tmp)
# pheatmap::pheatmap(tmp)
# 
# 
# tmp = classMat[rownames(classMat) %in% conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% unique(inputTb_metResponsiove$condID)], ]
# pheatmap::pheatmap(tmp)
# sum(rowSums(tmp) == 0)

# show the classification of all genes that is iCEL responsive (at least two up or two down)
upTbl = inputTb_metResponsiove[inputTb_metResponsiove$log2FoldChange_raw > 0, ]
ct1 = table(upTbl$condID)
downTbl = inputTb_metResponsiove[inputTb_metResponsiove$log2FoldChange_raw < 0, ]
ct2 = table(downTbl$condID)
icel_resp = union(names(ct1)[ct1 > 1], names(ct2)[ct2 > 1])
tmp = classMat[rownames(classMat) %in% conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% icel_resp],]
colSums(tmp)
#pheatmap::pheatmap(tmp, labels_row = iCELnames$ICELgene[match(rownames(tmp), iCELnames$WormBase_Gene_ID)])
#pheatmap::pheatmap(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')], labels_row = iCELnames$ICELgene[match(rownames(tmp), iCELnames$WormBase_Gene_ID)])
# pdf('tmp.pdf')
# pheatmap::pheatmap(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')], labels_row = iCELnames$ICELgene[match(rownames(tmp), iCELnames$WormBase_Gene_ID)],fontsize_row = 2)
# dev.off()

sum(rowSums(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')])>0)/nrow(tmp)
# 77% valid RNAi targets was included in the modeling framework 
sum(rowSums(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')])==1)/nrow(tmp)
# 54% was assigned to a unique classification
# total number of analyzable condition is 
n_total = sum(rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp, str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0 )
n_total/length(icel_resp)
# 78% valid RNAi conditions are analyzed in the modeling framework 
# check for the coverage of the unclustered conditions
uniConds = icel_resp[rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp, str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0]
RNAiclusters = read.csv('./../../2_DE/output/RNAi_groups.csv')
sum(RNAiclusters$clusters[(str_replace(RNAiclusters$RNAiID,' ','_') %in% uniConds)] == -1)/sum(RNAiclusters$clusters==-1)
# although the overall coverage of the RNAi conditions is at similar level (75% vs 78%), the FBA modeling is unbiased 
# so it covers 44% of unclustered conditions and it assesses all iCEL DEG instead of just selected coexpression clusters
# so this justifies it as a systems-level validation
length(intersect(rownames(classMat)[rowSums(classMat[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0], inputTb_metResponsiove$WBID))/(length(unique(inputTb_metResponsiove$WBID)))
# it covers 66% of DEG space 


# simulate the simple rewiring model and calculate the percetage explained by the model 
# assume the major energy drain is: protein syn, lipid syn, and glycan syn
# define the rewiring model 
# the principle is simplified into one sentence: activate compromised obj and repress all other objs
# we need to define the obj compromise for each RNAi 

# OF NOTE: to keep this modeling simple, the class selection by DE similarity is not performed in this step; the rewiring 
# model calculation is fully based on FBA classification and the multi-class genes were used literally (assuming it affects multiple objectives)
# and the unclassified conditions will be left out from the analysis. This also guarantees the same set of RNAi was analyzed in every randomization 
# we get the conditions to analyze (icel_responsive (at least 2 up or 2 down DEG) and classified)
total_condition_analyzed = c()
for (i in 1:length(icel_resp)){
  if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
    total_condition_analyzed = c(total_condition_analyzed, icel_resp[i])
  }
}

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
    
    # determine compromised obj
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


# ==> may be quantify by single and calculate pvalue for each edge and together for single only again
# ==> we can provide two quantification figure, either single or weighted multiple
# ==> integrate the calculation of overall explainability by only single obj


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
  # normalize within each type of RNAi-DEG pair 
  for (i in 1:length(modelObj)){
    for(j in 1:length(modelObj)){
      g_tbl$proportion[g_tbl$RNAi == modelObj[i] & g_tbl$DEG == modelObj[j]] =  g_tbl$proportion[g_tbl$RNAi == modelObj[i] & g_tbl$DEG == modelObj[j]] / sum(g_tbl$proportion[g_tbl$RNAi == modelObj[i] & g_tbl$DEG == modelObj[j]])
    }
  }
  return(g_tbl)
}

# visualization
g = model_fit_single$raw_counts
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




# try a randomization to assess the significance 
nRand = 10
set.seed(1126)
# we have a ton to randomize
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


randomizeFBA = FALSE


for (nn in 1:nRand){
  load(paste('./randomGRN/randNet1000_ES_batch_',nn,'.Rdata',sep = ''))
  for (zz in 1:1000){
    if(randomizeFBA){
      # generate the random classification matrix
      new_gene_names = rep(NA, nrow(classMat))
      new_gene_names[sample(hasClassInd, length(mustHasClass))] = mustHasClass
      new_gene_names[is.na(new_gene_names)] = justRandom[sample(length(justRandom))]
      
      classMat_rand = classMat
      rownames(classMat_rand) = new_gene_names
    }else{
      # we dont randomize the classification matrix, we do on the GRN
      classMat_rand = classMat
    }
    
    # generate random GRN
    inputTb_metResponsiove_rand = randNets[[zz]]
  
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
      DEGs_up = inputTb_metResponsiove_rand$WBID[inputTb_metResponsiove_rand$RNAiID == myCond & inputTb_metResponsiove_rand$log2FoldChange_raw > 0]
      DEGs_down = inputTb_metResponsiove_rand$WBID[inputTb_metResponsiove_rand$RNAiID == myCond & inputTb_metResponsiove_rand$log2FoldChange_raw < 0]
      
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
    rand_stat$edge_single$raw_count[paste('rand',nn,'_',zz,sep = '')] = make_g_tbl(model_fit_single$raw_counts)$proportion
    rand_stat$edge_single$norm_count[paste('rand',nn,'_',zz,sep = '')] = make_g_tbl(model_fit_single$normalized_counts)$proportion
    rand_stat$edge_full$raw_count[paste('rand',nn,'_',zz,sep = '')] = make_g_tbl(model_fit_all$raw_counts)$proportion
    rand_stat$edge_full$norm_count[paste('rand',nn,'_',zz,sep = '')] = make_g_tbl(model_fit_all$normalized_counts)$proportion
    
    
    print(paste(nn,'...',zz))
  }
  # list(list(mean(rewire_rate_rand), sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5])))
}

save(file = 'randomization_result_edge_direction_randGRN.Rdata',rand_stat)
hist(rand_stat$default_CR_full$rand_rate)
hist(rand_stat$default_CR_single$rand_rate)
pdf(paste('figures/full_randomization_results_edge_direction_randGRN.pdf',sep = ''),width = 7,height = 6)
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


# when calculate the enrichment, calculate the reletive frequency against the opposite direction instead of total frequency (
# calculate the proportion by up/(up+down) for each RNAi-DEG node pair)
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
write.csv(edge_quant,'figures/edge_direction_test_randGRN.csv')



# combine the two ways of edge assessment and find the final decision of edge quantification
edge_prop = read.csv('figures/edge_quantification_randGRN.csv')
edge_direct = read.csv('figures/edge_direction_test_randGRN.csv')
edge_info = cbind(edge_prop, edge_direct)
edge_info = edge_info[,c(2,3,4,
                         5,17,6,18,7,19,8,20,
                         9,21, 10,22,11,23,12,24)]
write.csv(edge_info,'figures/edge_test_info_randGRN.csv')

# visualization
edge_info = edge_info[order(edge_info$RNAi), ]
library(pheatmap)
pheatmap(-log10(edge_info[,12:ncol(edge_info)]),breaks = seq(0,3,length.out = 100),
         cluster_rows = F, cluster_cols = F, labels_row = paste(edge_info$RNAi, '->',edge_info$DEG),annotation_row = edge_info[,"type", drop = F])

pheatmap(-log10(edge_info[,12:ncol(edge_info)]),breaks = seq(0,3,length.out = 100),
         cluster_rows = T, cluster_cols = T, labels_row = paste(edge_info$RNAi, '->',edge_info$DEG),annotation_row = edge_info[,"type", drop = F])



# use best p ('OR') among the normalized calculation for both CR and prop seems most reasonable. 
# the hub strongly influence interpretation, so normalization makes sense. Single and full generally 
# agrees, so we just go with full to be consistent with the Main figure. But need to prepare evidence of 
# consistency (i.e, same figure by single data) and also compare what edge is missed if we use normalized only 
# as compared with simply best p value among all. 

# p-value is only to filter the edges. The length is still given by the normalized proportion (tentitively with full data)

# how to filter: by CR, by total, or by both
# what to use: norm or raw 
# what data to use: single or full
# assumed parameter: display by total proportion; p value 0.1; full for now
# key to decide: interpret by norm or by raw; how to combine the two tests (CR and total)
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


pdf('figures/fitted_CR_interaction_diagram_randGRN.pdf',width = 10)

sig_cutoff = 0.1
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



# save p-value and proportion heatmaps with selected edge labeled (to show the qualitative consistency)
pdf('figures/fitted_CR_parameter_robustness_randGRN.pdf',width = 10)
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



