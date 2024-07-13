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
score_cutoffs = c(0, 1e-6, 1e-5, 1e-4,1e-3, 1e-2, 0.1, 0.5, 0.99)
pvalues = c()
total_DEG_explained = c()
for (score_cutoff in score_cutoffs){
  # repeat the same code as in the FBA-fused CR model analysis, just with different threhsold 
    
  # load the met gene classifications
  scoreMat = read.csv('output/delta_flux_scoreMat.csv',row.names = 1)
  classMat = as.data.frame(1*(scoreMat > score_cutoff))
  
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
  
  # remove DEGs that are RNAi targeted genes
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
  
  
  # some numbers
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
  #length(intersect(rownames(classMat)[rowSums(classMat[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0], inputTb_metResponsiove$WBID))/(length(unique(inputTb_metResponsiove$WBID)))

  
  # calculate the proportion explained
  total_condition_analyzed = c()
  for (i in 1:length(icel_resp)){
    if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
      total_condition_analyzed = c(total_condition_analyzed, icel_resp[i])
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
      DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
      DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
      
      # determine compromised obj
      labeledObj = colnames(classMat)[classMat[conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') == myCond],]==1]
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
  
  pdf(paste('figures/sensitivity_analysis/SENSITIVITY_ANALYSIS_FBA_DEG_explained_distribution_FBA_fused_model_',score_cutoff,'.pdf',sep = ''),
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
    
    for (condInd in 1:length(total_condition_analyzed)){
      myCond = total_condition_analyzed[condInd]
      DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
      DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
      
      # determine compromised obj
      labeledObj = colnames(classMat_rand)[classMat_rand[conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') == myCond],]==1]
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
  save(file = paste('output/sensitivity_analysis/SENSITIVITY_ANALYSIS_randomization_result_FBA_fused_model_cutoff_',score_cutoff,'.Rdata',sep = ''),rand_rate, rand_rate2, rand_rate_up, rand_rate2_up, rand_rate_down, rand_rate2_down)
  hist(rand_rate)
  hist(rand_rate2)
  pdf(paste('figures/sensitivity_analysis/SENSITIVITY_ANALYSIS_overall_DEG_explained_randomization_FBA_fused_model_',score_cutoff,'.pdf',sep = ''),width = 7,height = 6)
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
save(score_cutoffs,pvalues,total_DEG_explained, file = 'output/sensitivity_analysis/data_to_plot.Rdata')
pdf('figures/parameter_sensitivity_to_objScore_cutoff.pdf',width = 7,height = 6)
plot(-log10(score_cutoffs+1e-7), -log10(pvalues),type = 'b')
plot(-log10(score_cutoffs+1e-7), total_DEG_explained,type = 'b')
dev.off()
