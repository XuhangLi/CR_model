# supplementary analysis for "2_DEG_modeling_basic_CR_model".

# codes are generally the same except for only using genes with single core function association. 
# this script is minimally commented to avoid redundant comments with "2_DEG_modeling_basic_CR_model"

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)

# load the met gene classifications
classMat = read.csv('output/FBA_classification_matrix.csv',row.names = 1)
iCELnames = read.csv('./../../input_data/otherTbls/iCEL_IDtbl.csv')
rownames(classMat) = iCELnames$WormBase_Gene_ID[match(rownames(classMat),iCELnames$ICELgene)]
# 02152023 since mrpl-44 is an multi-obj gene in FBA, so removing it or not will not change the analysis of single-obj genes here
# so we didnt rerun the code with mrpl-44 removed

# exclude the multiple obj genes 
modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
classMat[rowSums(classMat[,modelObj]) > 1,] = 0 


# load DEG result
inputTb=read.csv('./../../2_DE/output/DE_merged_clean_pcutoff_0.005_master_table_FDR2d0.1_FC1.5_FCtype_log2FoldChange_raw_ALL.csv')
inputTb$RNAiID = paste(inputTb$RNAi, inputTb$batchID)
inputTb$RNAi_geneName = inputTb$RNAi
inputTb$RNAi_geneName = str_replace(inputTb$RNAi_geneName,'^x.','')
inputTb$RNAi_geneName = str_replace_all(inputTb$RNAi_geneName,'_','-')
inputTb$RNAi_geneName = str_replace(inputTb$RNAi_geneName,'-L[0-9]$','')
inputTb = inputTb[-which(inputTb$RNAi_geneName == inputTb$Gene_name),]
inputTb$condID = paste(inputTb$RNAi,inputTb$batchID,sep = '_')
if (length(which(inputTb$WBID=="NoHit"))>0){
  inputTb=inputTb[-which(inputTb$WBID=="NoHit"),]}
# filtering the low responsive and non-metabolic ones
conditionInfo = read.csv('./../../2_DE/output/RNAi_condition_metaInfo.csv',row.names = 1)
inputTb_metResponsiove = inputTb[inputTb$RNAiID %in% conditionInfo$RNAiID[conditionInfo$isICEL & conditionInfo$isResponsive], ]
inputTb_metResponsiove=inputTb_metResponsiove[,c("WBID","RNAi","log2FoldChange_raw","condID")]

# exclude the non-canonical metabolic genes from the analysis 
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
pheatmap::pheatmap(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')], labels_row = iCELnames$ICELgene[match(rownames(tmp), iCELnames$WormBase_Gene_ID)])

# check numbers 
sum(rowSums(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')])>0)/nrow(tmp)
sum(rowSums(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')])==1)/nrow(tmp)
# total number of analyzable condition is 
n_total = sum(rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp, str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0 )
n_total/length(icel_resp)
# check for the coverage of the unclustered conditions
uniConds = icel_resp[rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp, str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0]
RNAiclusters = read.csv('./../../2_DE/output/RNAi_groups.csv')
sum(RNAiclusters$clusters[(str_replace(RNAiclusters$RNAiID,' ','_') %in% uniConds)] == -1)/sum(RNAiclusters$clusters==-1)


# visualization 
total_condition_analyzed = c()

obj_perturb = 'lipid'
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

for (obj_perturb in modelObj){
  
  # ignore these if..else.. statement: it is from an old script but here we filtered the classification 
  # matrix to single-core-function-only so the if..else.. statement didnt do anything.
  # basically equivalent to "obj_genes = intersect(rownames(classMat)[classMat[,obj_perturb]==1],conditionInfo$RNAi_WBID)"
  if (obj_perturb == 'energy'){
    # energy has to be only energy 
    obj_genes = intersect(rownames(classMat)[classMat[,obj_perturb]==1 & rowSums(classMat[,modelObj]) <= 1],conditionInfo$RNAi_WBID)
  }else{
    # if it is energy + one biomass obj; it is considered as the biomass obj
    # if it is more than one biomass obj, it is considered as multi-obj
    obj_genes = intersect(rownames(classMat)[classMat[,obj_perturb]==1 & rowSums(classMat[,setdiff(modelObj,'energy')]) <= 1],conditionInfo$RNAi_WBID)
  }
  myconds = intersect(str_replace(conditionInfo$RNAiID[conditionInfo$RNAi_WBID %in% obj_genes],' ','_'), inputTb_metResponsiove$condID)
  
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
      
      # apply the majority vote 
      if (any(rownames(classMat) %in% DEGs)){
        # rank the obj candidates 
        nCalls = colSums(classMat[rownames(classMat) %in% DEGs, modelObj])
        orders = names(nCalls)[order(nCalls,decreasing = T)]
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
  
  pdf(paste('figures/special_excludeMultiObj_FBA_classification_',obj_perturb,'_uniObjGenes.pdf',sep = ''),
      height = pdfHeight[obj_perturb], width = pdfWidth[obj_perturb])
  print(ggarrange(p2, p1, nrow = 1))
  dev.off()
  total_condition_analyzed = c(total_condition_analyzed, as.character(unique(stats_long_merge$RNAi_cond)))
}


#### the following section was omitted bc it is not in the ms #####
# # simulate the simple rewiring model and calculate the percetage explained by the model 
# # assume the major energy drain is: protein syn, lipid syn, and glycan syn
# # define the rewiring model 
# # the principle is simplified into one sentence: activate compromised obj and repress all other objs
# # we need to define the obj compromise for each RNAi 
# model = list()
# for (obj in modelObj){
#   model[[obj]] = list()
# }
# model[['energy']] = c('energy','lipid','pro_modi','pro_syn') # we assume Nucl Acid is not dependent on energy
# model[['lipid']] = c('lipid')
# model[['pro_modi']] = c('lipid','pro_modi') # pro_modi RNAi causes ER dysfunction which in turn compromise lipid syn
# model[['pro_syn']] = 'pro_syn'
# model[['nucl_acid']] = 'nucl_acid'
# 
# model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
# model_explained$UP_yes = 0
# model_explained$UP_no = 0
# model_explained$DOWN_yes = 0
# model_explained$DOWN_no = 0
# 
# for (condInd in 1:length(total_condition_analyzed)){
#     myCond = total_condition_analyzed[condInd]
#     DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
#     DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
#     
#     # determine compromised obj
#     labeledObj = colnames(classMat)[classMat[conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') == myCond],]==1]
#     affectedObj = c()
#     for (i in 1:length(labeledObj)){
#       affectedObj = c(affectedObj, model[[labeledObj[i]]])
#     }
#     
#     # compare with DEG 
#     # up genes
#     subClassMat = classMat[rownames(classMat) %in% DEGs_up, modelObj]
#     # calculate the DEG explained by model
#     model_explained$UP_yes[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) > 0)
#     model_explained$UP_no[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) == 0)
#     # down genes
#     subClassMat = classMat[rownames(classMat) %in% DEGs_down, modelObj]
#     # calculate the DEG explained by model
#     model_explained$DOWN_yes[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) > 0)
#     model_explained$DOWN_no[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) == 0)
#   
# }
# 
# # plot 
# b = model_explained$RNAi_cond
# # sort by the total explained rate
# tmp = model_explained[,2:5]
# tmp = tmp / rowSums(tmp)
# rewire_rate = rowSums(tmp[,c(1,3)])
# rewire_rate[is.na(rewire_rate)] = 0
# b[rank(rewire_rate,ties.method ='first')] =b
# model_explained$RNAi_cond = factor(model_explained$RNAi_cond,levels = b)
# 
# model_explained_long_up = reshape2::melt(model_explained[,c("RNAi_cond","UP_yes","UP_no")])
# model_explained_long_down = reshape2::melt(model_explained[,c("RNAi_cond","DOWN_yes","DOWN_no" )])
# model_explained_long = model_explained_long_up
# model_explained_long$down = model_explained_long_down$value
# colnames(model_explained_long) = c("RNAi_cond", "variable", "up", "down") 
# model_explained_long$variable = as.character(model_explained_long$variable)
# model_explained_long$variable[model_explained_long$variable == 'UP_yes'] = 'yes'
# model_explained_long$variable[model_explained_long$variable == 'UP_no'] = 'no'
# model_explained_long$variable = factor(model_explained_long$variable, levels = c('no','yes'))
# # calculate the average
# up_ave1 = mean(model_explained$UP_yes/(model_explained$UP_yes + model_explained$UP_no),na.rm = T) # mean of the rate
# down_ave1 = mean(model_explained$DOWN_yes/(model_explained$DOWN_yes + model_explained$DOWN_no),na.rm = T) # mean of the rate
# up_ave2 = sum(model_explained_long$up[model_explained_long$variable=='yes'])/sum(model_explained_long$up) # total rate
# down_ave2 = sum(model_explained_long$down[model_explained_long$variable=='yes'])/sum(model_explained_long$down) # total rate 
# # the total DEG explained is less sensitive to conditions that have very few DEG, so we use this as the reference line
# 
# # plot up and down with the same condition order for later merging in graphics
# p1 <- ggplot(model_explained_long, aes(x=RNAi_cond, y = up,fill=variable)) +
#   geom_bar(stat="identity", position="fill",colour="black", size = 0.1) + 
#   scale_y_continuous(labels = scales::percent,expand = c(0,0))+
#   theme_bw()+
#   theme(axis.text.x = element_text(colour ='black'),axis.title.x=element_text(colour ='black'),
#         axis.ticks.y=element_blank(),
#         axis.text.y = element_blank(),axis.title.y=element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(colour ='black'),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
#   coord_flip()+
#   scale_fill_manual("legend", values = c("yes" = "#FF6701", 
#                                          "no" = 'grey')) +
#   labs(y = 'Up DEG')+
#   geom_text(aes(label = up), position = position_fill(vjust = 0.5),size = 1)+
#   geom_hline(yintercept = up_ave2,color="white",
#              linetype="dashed")
# 
# #p1
# 
# # plot down 
# p2 <- ggplot(model_explained_long, aes(x=RNAi_cond, y = -down,fill=variable)) +
#   geom_bar(stat="identity", position="fill",colour="black", size = 0.1) + 
#   scale_y_continuous(labels = scales::percent,expand = c(0,0))+
#   theme_bw()+
#   theme(axis.text.x = element_text(colour ='black'),axis.title.x=element_text(colour ='black'),
#         axis.text.y = element_text(colour ='black', size = 2),
#         legend.position = "none",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
#   coord_flip()+
#   scale_fill_manual("legend", values = c("yes" = "#FF6701", 
#                                          "no" = 'grey')) +
#   labs(y = 'Down DEG')+
#   geom_text(aes(label = -down), position = position_fill(vjust = 0.5),size = 1)+
#   geom_hline(yintercept = -down_ave2,color="white",
#              linetype="dashed")
# 
# #p2
# 
# # pdf(paste('figures/FBA_DEG_explained_distribution.pdf',sep = ''),
# #     height = 7, width = 14)
# print(ggarrange(p2, p1, nrow = 1))
# # dev.off()

# 
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
# nRand = 1000
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
# # save(file = 'randomization_result.Rdata',rand_rate, rand_rate2, rand_rate_up, rand_rate2_up, rand_rate_down, rand_rate2_down)
# hist(rand_rate)
# hist(rand_rate2)
# # pdf(paste('figures/overall_DEG_explained_randomization.pdf',sep = ''),width = 7,height = 6)
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

# dev.off()
