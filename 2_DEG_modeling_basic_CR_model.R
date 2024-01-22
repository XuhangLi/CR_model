# Summary
# Visualize the DEGs based on their core functions assigned by FBA simulation. 

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)

# load the met gene classifications
classMat = read.csv('output/FBA_classification_matrix.csv',row.names = 1)
iCELnames = read.csv('./../../input_data/otherTbls/iCEL_IDtbl.csv')
rownames(classMat) = iCELnames$WormBase_Gene_ID[match(rownames(classMat),iCELnames$ICELgene)]

# mrpl-44 is a misannotated gene in the model; in the model it should be mccc-2 but mislabeled as mrpl-44 (which is mito ribo gene), since ribo genes are excluded from the analysis, we also exclude mrpl-44
classMat = classMat[-which(rownames(classMat) == 'WBGene00008514'),]

# load the DE similarity 
DEsim = read.csv('./../../2_DE/output/cosineSimilarity_FC_denoised_stationery_metabolic.csv', row.names = 1)
# also remove the mrpl-44 condition from the analysis (we remove both from DEG and conditions since this DEG is not iCEL either)
DEsim = DEsim[-which(rownames(DEsim) == 'x.mrpl_44_met3_lib1'),]
DEsim = DEsim[,-which(colnames(DEsim) == 'x.mrpl_44_met3_lib1')]

# load DEG result
inputTb=read.csv('./../../2_DE/output/DE_merged_clean_pcutoff_0.005_master_table_FDR2d0.1_FC1.5_FCtype_log2FoldChange_raw_ALL.csv')
# remove the mrpl-44 condition from the analysis (we remove both from DEG and conditions since this DEG is not iCEL either)
inputTb = inputTb[-which(inputTb$RNAi == 'x.mrpl_44'),]
inputTb = inputTb[inputTb$WBID != 'WBGene00008514',]

# remove the DEGs that is the RNAi targeted gene 
inputTb$RNAiID = paste(inputTb$RNAi, inputTb$batchID)
inputTb$RNAi_geneName = inputTb$RNAi
inputTb$RNAi_geneName = str_replace(inputTb$RNAi_geneName,'^x.','')
inputTb$RNAi_geneName = str_replace_all(inputTb$RNAi_geneName,'_','-')
inputTb$RNAi_geneName = str_replace(inputTb$RNAi_geneName,'-L[0-9]$','')
inputTb = inputTb[-which(inputTb$RNAi_geneName == inputTb$Gene_name),]
inputTb$condID = paste(inputTb$RNAi,inputTb$batchID,sep = '_')
if (length(which(inputTb$WBID=="NoHit"))>0){
  inputTb=inputTb[-which(inputTb$WBID=="NoHit"),]}

# filtering out nonresponsive and non-metabolic perturbations
conditionInfo = read.csv('./../../2_DE/output/RNAi_condition_metaInfo.csv',row.names = 1)
# only analyze the iCEL responsive
inputTb_metResponsiove = inputTb[inputTb$RNAiID %in% conditionInfo$RNAiID[conditionInfo$isICEL & conditionInfo$isResponsive], ]
inputTb_metResponsiove=inputTb_metResponsiove[,c("WBID","RNAi","log2FoldChange_raw","condID")]

# filtering of some irrelevant or ambiguously annotated

# only keep iCEL DEG
inputTb_metResponsiove = inputTb_metResponsiove[inputTb_metResponsiove$WBID %in% rownames(classMat),]

# exclude the those non-canonical metabolic genes from the analysis 
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
# exclude RNA pol and DNA pol (eq. to that ribosome is not in the model as well)
pollist = read.csv('input/pol_in_model.csv')
inputTb_metResponsiove = inputTb_metResponsiove[!(inputTb_metResponsiove$WBID %in% 
                                                    iCELnames$WormBase_Gene_ID[match(pollist$pol_in_model, iCELnames$ICELgene)]),]
inputTb_metResponsiove = inputTb_metResponsiove[!(conditionInfo$RNAi_WBID[match(inputTb_metResponsiove$condID, str_replace(conditionInfo$RNAiID,' ','_'))] %in% 
                                                    iCELnames$WormBase_Gene_ID[match(pollist$pol_in_model, iCELnames$ICELgene)]),]


# visualize the classification of all genes that are targeted in iCEL responsive (at least two up or two down) perturbations
upTbl = inputTb_metResponsiove[inputTb_metResponsiove$log2FoldChange_raw > 0, ]
ct1 = table(upTbl$condID)
downTbl = inputTb_metResponsiove[inputTb_metResponsiove$log2FoldChange_raw < 0, ]
ct2 = table(downTbl$condID)
icel_resp = union(names(ct1)[ct1 > 1], names(ct2)[ct2 > 1])
tmp = classMat[rownames(classMat) %in% conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% icel_resp],]
colSums(tmp)
pheatmap::pheatmap(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')], labels_row = iCELnames$ICELgene[match(rownames(tmp), iCELnames$WormBase_Gene_ID)])

# calculate some numbers
# at gene-level:
sum(rowSums(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')])>0)/nrow(tmp)
# 78.5% pass-filtered perturbations were included in the modeling framework 
sum(rowSums(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')])==1)/nrow(tmp)
# 55.5% were assigned to a unique classification
# at perturbation-level:
n_total = sum(rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp, str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0 )
n_total/length(icel_resp)
# 80% valid RNAi conditions are analyzed in the modeling framework 
# we assume FBA can analyze many more unclustered RNAi conditions: check for the coverage of the unclustered conditions
uniConds = icel_resp[rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp, str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0]
RNAiclusters = read.csv('./../../2_DE/output/RNAi_groups.csv')
sum(RNAiclusters$clusters[(str_replace(RNAiclusters$RNAiID,' ','_') %in% uniConds)] == -1)/sum(RNAiclusters$clusters==-1)
# so it covers 45% of unclustered conditions and it assesses all iCEL DEG instead of just selected coexpression clusters
# so this justifies it as a systems-level validation
length(intersect(rownames(classMat)[rowSums(classMat[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0], inputTb_metResponsiove$WBID))/(length(unique(inputTb_metResponsiove$WBID)))
# it covers 66% of DEG space 

# PART I: VISUALIZE DEG IN EACH PERTURBATION
# We first focus on perturbations whose targeted gene is associated with unique core function.
total_condition_analyzed = c()

obj_perturb = 'lipid' # to debug with - not meaningful here
modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
pdfHeight = c("energy" = 12,
            'lipid' = 7,
            'pro_modi' = 5,
            'pro_syn' = 3,
            'nucl_acid' = 2)
pdfWidth = c("energy" = 20,
             'lipid' = 14,
             'pro_modi' = 14,
             'pro_syn' = 14,
             'nucl_acid' = 14)

# loop through each core functions
for (obj_perturb in modelObj){
  # we now only simply show the single-core-function genes
  obj_genes = intersect(rownames(classMat)[classMat[,obj_perturb]==1 & rowSums(classMat[,modelObj]) == 1],conditionInfo$RNAi_WBID)
  # get corresponding conditions
  myconds = intersect(str_replace(conditionInfo$RNAiID[conditionInfo$RNAi_WBID %in% obj_genes],' ','_'), inputTb_metResponsiove$condID)
  # count for up and down DEGs seperately
  for (obj_direction in c('up','down')){
    stats = data.frame(RNAi_cond = myconds)
    stats$energy = 0
    stats$lipid = 0
    stats$pro_modi = 0
    stats$pro_syn = 0 
    stats$nucl_acid = 0 
    stats$noClass = 0
    for (i in 1:nrow(stats)){
      myCond = stats$RNAi_cond[i]
      if (obj_direction == 'up'){
        DEGs = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
      }else{
        DEGs = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
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
        stats[i,2:6] = colSums(subClassMat)
        stats[i,7] = sum(rowSums(subClassMat) == 0)
      }else{
        stats[i,2:7] = 0
      }
    }
    # before plotting...
    # filter out low response conditions
    stats = stats[rowSums(stats[,2:7])>1,] # we require at least two iCEL DE to enable meaningful interpretation (of proportion)
    # (if up or down gene is less than 2, we will mask the condition will NA (no show as blank bar))
    
    # plot
    b = stats$RNAi_cond
    # sort by the most abundant class
    tmp = stats[,2:7]
    tmp = tmp / rowSums(tmp)
    aveRewireRate = colMeans(tmp[,1:5])
    sortBy = names(aveRewireRate)[order(aveRewireRate,decreasing = T)][1]
    rewire_rate = stats[,sortBy] / rowSums(stats[,2:7])
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
                                           "pro_modi" = "#EDB120",
                                           "nucl_acid" = "#4DBEEE",
                                           "pro_syn" = "#0072BD",
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
                                           "pro_modi" = "#EDB120",
                                           "nucl_acid" = "#4DBEEE",
                                           "pro_syn" = "#0072BD",
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
  total_condition_analyzed = c(total_condition_analyzed, as.character(unique(stats_long_merge$RNAi_cond)))
}

# it turns out that 160 conditions were analyzed (unique obj association)

# Next, we visualize the genes with multiple core functions or the unclassified genes using imputation based on perturbation-perturbation similarity

# we use DE similarity to assign them into the most similar known category 
multi_obj_genes = intersect(rownames(classMat)[rowSums(classMat[,modelObj]) > 1],conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% icel_resp])
unclassified_genes = intersect(rownames(classMat)[rowSums(classMat[,modelObj]) == 0],conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% icel_resp])
imputedClass = data.frame(gene = c(multi_obj_genes, unclassified_genes), ori_class = c(rep('multi',length(multi_obj_genes)), rep('unclass',length(unclassified_genes))))
imputedClass$class = NA

# To assign the most likely labels, we first filtered out all conditions without a cloest condition of at least 0.2 cosine with labels.
# Next we looked at the top 5 most similar conditions with labels; the label was imputed if (a) all top 5 are the same label, or (2) 
#  the most abundant class (max sum of cosine) is very dominating the entire thing (the sum of cosine for a class is higher than the sum of cosine of any other one by two fold or more)
# or (3) the most abundant class (max sum of cosine) is much more similar than others (the max cosine for a class is higher than the max of any other one by 0.1 or more)
# of note, if the most similar class is different than the most abundant class, we dont assign any class because it is unclear (example is kynu-1)
for (i in 1:nrow(imputedClass)){
  # impute for genes one by one
  myGene = imputedClass$gene[i]
  simMat = DEsim[rownames(DEsim) %in% str_replace(conditionInfo$RNAiID[which(conditionInfo$isICEL)],' ','_') & 
                   !(rownames(DEsim) %in% str_replace(conditionInfo$RNAiID[which(conditionInfo$RNAi_WBID == myGene & conditionInfo$isICEL)],' ','_')), 
                 colnames(DEsim) %in% str_replace(conditionInfo$RNAiID[which(conditionInfo$RNAi_WBID == myGene & conditionInfo$isICEL)],' ','_'), drop = F]
  maxSim = rowMaxs(as.matrix(simMat))
  names(maxSim) = rownames(simMat)
  maxSim = as.data.frame(maxSim)
  # we use a hard cutoff of 0.2 (see the DE similarity clustering section for fitting distributions and get this cutoff (0.2 is rounded))
  maxSim = maxSim[maxSim$maxSim > 0.2,,drop=F]
  if (nrow(maxSim) > 0){
    maxSim$WBID = conditionInfo$RNAi_WBID[match(rownames(maxSim), str_replace(conditionInfo$RNAiID,' ','_'))]
    maxSim$class = NA
    for (j in 1:nrow(maxSim)){
      if (sum(classMat[maxSim$WBID[j],modelObj]) == 1){
        maxSim$class[j] = modelObj[classMat[maxSim$WBID[j],modelObj] == 1]
      }
    }
    # only look at the most significant ones (TOP 5 non-NA)
    maxSim = maxSim[!is.na(maxSim$class),]
    maxSim = maxSim[order(maxSim$ maxSim,decreasing = T),,drop = F]
    maxSim = maxSim[1:min(c(5, nrow(maxSim))), ,drop = F]
    
    allClass = setdiff(unique(maxSim$class),NA)
    sumCosine = c()
    for(j in 1:length(allClass)){
      sumCosine[j] = sum(maxSim$maxSim[maxSim$class == allClass[j]],na.rm = T)
    }
    bestClass = allClass[sumCosine == max(sumCosine)]
    # best class should be dominate 
    if (length(allClass)==0){
      # no similar RNAi that is labeled
      imputedClass$class[i] = 'no_similar_with_label'
    }else if (length(allClass)==1){
      imputedClass$class[i] = bestClass
    }else if(all(sumCosine[bestClass == allClass] > sumCosine[bestClass != allClass] * 2)){
      # if one class is very dominating the entire thing 
      imputedClass$class[i] = bestClass
    }else if(max(maxSim$maxSim[which(maxSim$class == bestClass)]) - max(maxSim$maxSim[which(maxSim$class != bestClass)]) > 0.1){
      # if one class is very similar but none of the others are 
      imputedClass$class[i] = bestClass
    }
    else {
      names(sumCosine) = allClass
      print(i)
      print(str_replace(conditionInfo$RNAiID[which(conditionInfo$RNAi_WBID == myGene & conditionInfo$isICEL)],' ','_'))
      print(sumCosine)
      print(maxSim[order(maxSim$maxSim,decreasing = T),])
      imputedClass$class[i] = 'no_major_class'
    }
    
  }else{
    # no similar RNAi 
    imputedClass$class[i] = 'no_similar'
  }
}

# check for conflicts with old class (the FBA assigned classes)
for (i in 1:nrow(imputedClass)){
  if (sum(classMat[imputedClass$gene[i],modelObj]) > 0 & imputedClass$class[i] %in% modelObj){ # is imputed and is originally assigned
    ori_class = modelObj[classMat[imputedClass$gene[i],modelObj] == 1]
    imp_class = imputedClass$class[i]
    if (length(intersect(imp_class, ori_class)) == 0){# imputed class is not overlapped with the FBA class 
      print(ori_class)
      print(imputedClass[i,])
      imputedClass$class[i] = 'conflict_with_FBA_class'
    }
  }
}

# look at how many genes for each type
sum(imputedClass$class == 'no_major_class',na.rm = T)
sum(imputedClass$class == 'no_similar',na.rm = T)
sum(imputedClass$class == 'no_similar_with_label',na.rm = T)
sum(imputedClass$class == 'conflict_with_FBA_class',na.rm = T)
sum(imputedClass$class %in% modelObj)
sum(imputedClass$class %in% modelObj) / nrow(imputedClass)

# how many perturbations for each type
sum(str_replace(conditionInfo$RNAiID[conditionInfo$RNAi_WBID %in% imputedClass$gene[imputedClass$class == 'no_major_class']],' ','_')
    %in% icel_resp)
sum(str_replace(conditionInfo$RNAiID[conditionInfo$RNAi_WBID %in% imputedClass$gene[imputedClass$class == 'no_similar']],' ','_')
    %in% icel_resp)
sum(str_replace(conditionInfo$RNAiID[conditionInfo$RNAi_WBID %in% imputedClass$gene[imputedClass$class == 'no_similar_with_label']],' ','_')
    %in% icel_resp)
sum(str_replace(conditionInfo$RNAiID[conditionInfo$RNAi_WBID %in% imputedClass$gene[imputedClass$class == 'conflict_with_FBA_class']],' ','_')
    %in% icel_resp)
sum(str_replace(conditionInfo$RNAiID[conditionInfo$RNAi_WBID %in% imputedClass$gene[imputedClass$class %in% modelObj]],' ','_')
    %in% icel_resp)

# It turns out that we extend our analysis to 70/119 unanalyzed conditions, leaving only 49 conditions not analyzed by the modeling 

# the total coverage of the modeling is 
RNAi_visualized = c()
singles = c()
multiple_imputed = c()
no_class_imputed = c()
for (i in 1:length(icel_resp)){
  if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) ==1 |# single biomass
      conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))] %in% imputedClass$gene[imputedClass$class %in% modelObj] ){
    RNAi_visualized = c(RNAi_visualized, icel_resp[i])
    if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) ==1 ){
      singles = c(singles, icel_resp[i])
    }
    
    if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 1){
      multiple_imputed = c(multiple_imputed, icel_resp[i])
    }
    
    if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) == 0){
      no_class_imputed = c(no_class_imputed, icel_resp[i])
    }
  }
}
length(RNAi_visualized)/length(icel_resp)
# 82% responsive RNAi conditions will be analyzed in the visualization, including 160+35=195 that originally has a core function
# and 35 that originally does not has a core function (total to 230)
# note this is a different number than the number of conditions analyzed in the next step (FBA-fused CR model) because that one uses
# all conditions with any core function associations under a simply cutoff (~80% of all conditions)

# visualize these imputed genes
types = c('multi','unclass')
modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
pdfHeight = c("energy" = 3.5,
              'lipid' = 3.5,
              'pro_modi' = 2.5,
              'pro_syn' = 1.5,
              'nucl_acid' = 1.5)
pdfWidth = c("energy" = 14,
             'lipid' = 14,
             'pro_modi' = 14,
             'pro_syn' = 14,
             'nucl_acid' = 14)

for (obj_perturb in modelObj){
  # same strategy as above - majority vote 
  pdf(paste('figures/FBA_classification_',obj_perturb,'_multiObj_and_unclassified_genes.pdf',sep = ''),
     height = pdfHeight[obj_perturb], width = pdfWidth[obj_perturb])
  for (type in types){
    obj_genes = imputedClass$gene[imputedClass$class == obj_perturb & imputedClass$ori_class == type]
    myconds = intersect(str_replace(conditionInfo$RNAiID[conditionInfo$RNAi_WBID %in% obj_genes],' ','_'), inputTb_metResponsiove$condID)
    if (length(myconds) > 0){
      for (obj_direction in c('up','down')){
        stats = data.frame(RNAi_cond = myconds)
        stats$energy = 0
        stats$lipid = 0
        stats$pro_modi = 0
        stats$pro_syn = 0 
        stats$nucl_acid = 0 
        stats$noClass = 0
        for (i in 1:nrow(stats)){
          myCond = stats$RNAi_cond[i]
          if (obj_direction == 'up'){
            DEGs = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
          }else{
            DEGs = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
          }
          
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
            stats[i,2:6] = colSums(subClassMat)
            stats[i,7] = sum(rowSums(subClassMat) == 0)
          }else{
            stats[i,2:7] = 0
          }
        }
        
        # filter out low response conditions
        stats = stats[rowSums(stats[,2:7])>1,] # we require at least two iCEL DE to enable meaningful interpretation (of proportion)
        
        # plot
        b = stats$RNAi_cond
        # sort by the most abundant class
        tmp = stats[,2:7]
        tmp = tmp / rowSums(tmp)
        aveRewireRate = colMeans(tmp[,1:5])
        sortBy = names(aveRewireRate)[order(aveRewireRate,decreasing = T)][1]
        rewire_rate = stats[,sortBy] / rowSums(stats[,2:7])
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
      
      # format the merged table
      stats_long_down$value = -stats_long_down$value
      extra_conds = levels(stats_long_down$RNAi_cond)[!(levels(stats_long_down$RNAi_cond) %in% levels(stats_long_up$RNAi_cond))]
      new_levels = c(extra_conds, levels(stats_long_up$RNAi_cond))
      stats_long_down$RNAi_cond = factor(as.character(stats_long_down$RNAi_cond), levels = new_levels)
      stats_long_up$RNAi_cond = factor(as.character(stats_long_up$RNAi_cond), levels = new_levels)
      
      stats_long_merge = stats_long_up
      stats_long_merge$up = stats_long_merge$value
      stats_long_merge$down = 0
      if (nrow(stats_long_down) > 0){
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
      }
      
      
      # plot up and down with the same condition order for later merging in graphics
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
                                               "pro_modi" = "#EDB120",
                                               "nucl_acid" = "#4DBEEE",
                                               "pro_syn" = "#0072BD",
                                               "noClass" = "white")) +
        labs(y = 'Up DEG')+
        geom_text(aes(label = up), position = position_fill(vjust = 0.5),size = 3)+
        geom_hline(yintercept = sum(classMat[,sortBy_up]) / nrow(classMat),color="white",
                   linetype="dashed")
      
      p1
      
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
                                               "pro_modi" = "#EDB120",
                                               "nucl_acid" = "#4DBEEE",
                                               "pro_syn" = "#0072BD",
                                               "noClass" = "white")) +
        labs(y = 'Down DEG')+
        geom_text(aes(label = -down), position = position_fill(vjust = 0.5),size = 3)+
        geom_hline(yintercept = -sum(classMat[,sortBy]) / nrow(classMat),color="white",
                   linetype="dashed")
      
      p2
      plot = ggarrange(p2, p1, nrow = 1)
      print(annotate_figure(plot, top = text_grob(type, 
                                                color = "black", face = "bold", size = 14))   )
      total_condition_analyzed = c(total_condition_analyzed, as.character(unique(stats_long_merge$RNAi_cond)))
    }
  }
  dev.off()
}


# Note: the following part is not included in the paper!

# simulate the simple rewiring model and calculate the percentage explained by the model 
# consider any associated core function as perturbed functions and others as unaffected 
# only exception is that energy perturbation will perturb all functions because energy underlies everything
total_condition_analyzed = c()
for (i in 1:length(icel_resp)){
  if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
    total_condition_analyzed = c(total_condition_analyzed, icel_resp[i])
  }
}

# define such as model
model = list()
for (obj in modelObj){
  model[[obj]] = list()
}
model[['energy']] = c('energy','lipid','pro_modi','pro_syn','nucl_acid') # energy RNAi perturbs all
model[['lipid']] = c('lipid') # other RNAi perturbs themself
model[['pro_modi']] = c('pro_modi') 
model[['pro_syn']] = 'pro_syn'
model[['nucl_acid']] = 'nucl_acid'

# calculate the proportion of DEG explained
model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
model_explained$UP_yes = 0
model_explained$UP_no = 0
model_explained$DOWN_yes = 0
model_explained$DOWN_no = 0

for (condInd in 1:length(total_condition_analyzed)){
    myCond = total_condition_analyzed[condInd]
    DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
    DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
    
    # determine compromised obj
    labeledObj = colnames(classMat)[classMat[conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') == myCond],]==1]
    affectedObj = c()
    for (i in 1:length(labeledObj)){
      affectedObj = c(affectedObj, model[[labeledObj[i]]])
    }
    
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

# sort by the total explained rate
b = model_explained$RNAi_cond
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

pdf(paste('figures/FBA_DEG_explained_distribution_simplest_model.pdf',sep = ''),
    height = 7, width = 14)
print(ggarrange(p2, p1, nrow = 1))
dev.off()

# this result is similar to the FBA-fused model - so we go with the simplest FBA fused model

##### the following codes are randomization tests of the above proportion, which was not included in the ms #####
# change these code to commented state

# # the overall average explained rate (up and down together) is 
# obs_rate = mean(rewire_rate)
# obs_total_DE_rate = sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5])
# 
# # the overall average explained rate (up and down seperately) is 
# # up
# tmp = model_explained[,2:3]
# tmp = tmp / rowSums(tmp)
# obs_rate_up = mean(tmp$UP_yes, na.rm = T)
# obs_total_DE_rate_up = sum(model_explained[,c(2)]) / sum(model_explained[,2:3])
# # down
# tmp = model_explained[,4:5]
# tmp = tmp / rowSums(tmp)
# obs_rate_down = mean(tmp$DOWN_yes, na.rm = T)
# obs_total_DE_rate_down = sum(model_explained[,c(4)]) / sum(model_explained[,4:5])
# 
# # try a randomization to assess the significance 
# nRand = 10000
# set.seed(1126)
# rand_rate = c()
# rand_rate2 = c()
# rand_rate_up = c()
# rand_rate2_up = c()
# rand_rate_down = c()
# rand_rate2_down = c()
# # we keep the number of labels for each obj and keep the analyzed conditions assigned to at least one label
# mustHasClass = unique(conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% total_condition_analyzed])
# justRandom = setdiff(rownames(classMat), mustHasClass)
# hasClassInd = as.numeric(which(rowSums(classMat[,modelObj])>0))
# 
# # library(doParallel)
# # myCluster <- makeCluster(10)
# # registerDoParallel(myCluster)
# 
# # system.time(
# # x <- foreach(nn=1:nRand, .combine='c') %do% {
# for (nn in 1:nRand){
#   # generate the random classification matrix
#   new_gene_names = rep(NA, nrow(classMat))
#   new_gene_names[sample(hasClassInd, length(mustHasClass))] = mustHasClass
#   new_gene_names[is.na(new_gene_names)] = justRandom[sample(length(justRandom))]
#   
#   classMat_rand = classMat
#   rownames(classMat_rand) = new_gene_names
#   
#   model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
#   model_explained$UP_yes = 0
#   model_explained$UP_no = 0
#   model_explained$DOWN_yes = 0
#   model_explained$DOWN_no = 0
#   
#   for (condInd in 1:length(total_condition_analyzed)){
#     myCond = total_condition_analyzed[condInd]
#     DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
#     DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
#     
#     # determine compromised obj
#     labeledObj = colnames(classMat_rand)[classMat_rand[conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') == myCond],]==1]
#     affectedObj = c()
#     for (i in 1:length(labeledObj)){
#       affectedObj = c(affectedObj, model[[labeledObj[i]]])
#     }
#     
#     # compare with DEG 
#     # up genes
#     subClassMat = classMat_rand[rownames(classMat_rand) %in% DEGs_up, modelObj]
#     # calculate the DEG explained by model
#     model_explained$UP_yes[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) > 0)
#     model_explained$UP_no[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) == 0)
#     # down genes
#     subClassMat = classMat_rand[rownames(classMat_rand) %in% DEGs_down, modelObj]
#     # calculate the DEG explained by model
#     model_explained$DOWN_yes[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) > 0)
#     model_explained$DOWN_no[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) == 0)
#     
#   }
#   
#   # total explained rate
#   tmp = model_explained[,2:5]
#   tmp = tmp / rowSums(tmp)
#   rewire_rate_rand = rowSums(tmp[,c(1,3)])
#   rewire_rate_rand[is.na(rewire_rate_rand)] = 0
#   
#   # the overall average explained rate (up and down together) is 
#   rand_rate = c(rand_rate, mean(rewire_rate_rand))
#   rand_rate2 = c(rand_rate2, sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5]))
#   
#   # up
#   tmp = model_explained[,2:3]
#   tmp = tmp / rowSums(tmp)
#   rand_rate_up = c(rand_rate_up, mean(tmp$UP_yes, na.rm = T))
#   rand_rate2_up = c(rand_rate2_up, sum(model_explained[,c(2)]) / sum(model_explained[,2:3]))
#   # down
#   tmp = model_explained[,4:5]
#   tmp = tmp / rowSums(tmp)
#   rand_rate_down = c(rand_rate_down, mean(tmp$DOWN_yes, na.rm = T))
#   rand_rate2_down = c(rand_rate2_down, sum(model_explained[,c(4)]) / sum(model_explained[,4:5]))
#   
#   print(nn)
#   # list(list(mean(rewire_rate_rand), sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5])))
# }
# 
# save(file = 'randomization_result_simplest_model.Rdata',rand_rate, rand_rate2, rand_rate_up, rand_rate2_up, rand_rate_down, rand_rate2_down)
# hist(rand_rate)
# hist(rand_rate2)
# pdf(paste('figures/overall_DEG_explained_randomization_simplest_model.pdf',sep = ''),width = 7,height = 6)
# hist(rand_rate, main = paste('p <',signif((1+sum(rand_rate>=obs_rate))/(1+length(rand_rate)),2)),
#      xlim = c(min(rand_rate)*0.9,max(c(rand_rate,obs_rate))*1.1),
#      xlab = 'condition-wise average percentage of DEG explained')
# abline(v = obs_rate)
# 
# hist(rand_rate2, main = paste('p <',signif((1+sum(rand_rate2>=obs_total_DE_rate))/(1+length(rand_rate2)),2)),
#      xlim = c(min(rand_rate2)*0.9,max(c(rand_rate2,obs_total_DE_rate))*1.1),
#      xlab = 'percentage of total DEG explained')
# abline(v = obs_total_DE_rate)
# 
# hist(rand_rate_up, main = paste('p <',signif((1+sum(rand_rate_up>=obs_rate_up))/(1+length(rand_rate_up)),2)),
#      xlim = c(min(rand_rate_up)*0.9,max(c(rand_rate_up,obs_rate_up))*1.1),
#      xlab = 'condition-wise average percentage of up DEG explained')
# abline(v = obs_rate_up)
# 
# hist(rand_rate2_up, main = paste('p <',signif((1+sum(rand_rate2_up>=obs_total_DE_rate_up))/(1+length(rand_rate2_up)),2)),
#      xlim = c(min(rand_rate2_up)*0.9,max(c(rand_rate2_up,obs_total_DE_rate_up))*1.1),
#      xlab = 'percentage of total UP DEG explained')
# abline(v = obs_total_DE_rate_up)
# 
# hist(rand_rate_down, main = paste('p <',signif((1+sum(rand_rate_down>=obs_rate_down))/(1+length(rand_rate_down)),2)),
#      xlim = c(min(rand_rate_down)*0.9,max(c(rand_rate_down,obs_rate_down))*1.1),
#      xlab = 'condition-wise average percentage of down DEG explained')
# abline(v = obs_rate_down)
# 
# hist(rand_rate2_down, main = paste('p <',signif((1+sum(rand_rate2_down>=obs_total_DE_rate_down))/(1+length(rand_rate2_down)),2)),
#      xlim = c(min(rand_rate2_down)*0.9,max(c(rand_rate2_down,obs_total_DE_rate_down))*1.1),
#      xlab = 'percentage of total down DEG explained')
# abline(v = obs_total_DE_rate_down)
# 
# dev.off()

##### plot the Euler graph for core function classifications #####

gene_list = list(
  iCEL = rownames(classMat),
  energy = rownames(classMat)[classMat$energy==1],
  lipid = rownames(classMat)[classMat$lipid==1],
  pro_modi = rownames(classMat)[classMat$pro_modi==1],
  pro_syn = rownames(classMat)[classMat$pro_syn==1],
  nucl_acid = rownames(classMat)[classMat$nucl_acid==1]
)

library(eulerr)
set.seed(1030)
fit2 <- euler(classMat[,modelObj], 
              shape = "ellipse",
              #loss = 'region',
              #loss_aggregator = 'max', 
              control = list(extraopt = TRUE, extraopt_threshold = 0)
              )
pdf(paste('figures/euler_plot_iCEL_gene_classification.pdf',sep = ''),width = 7,height = 7)
plot(fit2, quantities = TRUE, fill = c('#FF6701','#77AC30','#EDB120','#0072BD','#4DBEEE'),main = paste('unclassified =', sum(rowSums(classMat[,modelObj]==1)==0),'(',100*round(sum(rowSums(classMat[,modelObj]==1)==0)/nrow(classMat),2),'% )'))
dev.off()

df = as.data.frame(fit2$original.values)
write.csv(df, 'figures/euler_plot_iCEL_gene_classification_breakdown.csv')

# we cannot get a visually plausible figure by changing seeds, so we give up on getting this from the function
# itself, we will just change the numbers in the original figure - it is generally fine because the changes in
# the area is very minor 

# error_plot(fit2)


