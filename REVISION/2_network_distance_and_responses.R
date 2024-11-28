# Summary
# Visualize the DEGs based on their core functions assigned by FBA simulation. 

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)

# load icel genes
classMat = read.csv('./../CR_model_final_run/output/FBA_classification_matrix.csv',row.names = 1)
iCELnames = read.csv('./../../input_data/otherTbls/iCEL_IDtbl.csv')
iCELgenes = iCELnames$WormBase_Gene_ID[match(rownames(classMat),iCELnames$ICELgene)]

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
inputTb_metResponsiove = inputTb_metResponsiove[inputTb_metResponsiove$WBID %in% iCELgenes,]

# exclude the those non-canonical metabolic genes from the analysis 
# exclude UGT (both for RNAi and DEG)
ugtlist = read.csv('./../CR_model_final_run/input/ugt_in_model.csv')
inputTb_metResponsiove = inputTb_metResponsiove[!(inputTb_metResponsiove$WBID %in% 
                                                    iCELnames$WormBase_Gene_ID[match(ugtlist$ugt_in_model, iCELnames$ICELgene)]),]
inputTb_metResponsiove = inputTb_metResponsiove[!(conditionInfo$RNAi_WBID[match(inputTb_metResponsiove$condID, str_replace(conditionInfo$RNAiID,' ','_'))] %in% 
                                                    iCELnames$WormBase_Gene_ID[match(ugtlist$ugt_in_model, iCELnames$ICELgene)]),]
# exclude vha
vhalist = read.csv('./../CR_model_final_run/input/vha_in_model.csv')
inputTb_metResponsiove = inputTb_metResponsiove[!(inputTb_metResponsiove$WBID %in% 
                                                    iCELnames$WormBase_Gene_ID[match(vhalist$vha_in_model, iCELnames$ICELgene)]),]
inputTb_metResponsiove = inputTb_metResponsiove[!(conditionInfo$RNAi_WBID[match(inputTb_metResponsiove$condID, str_replace(conditionInfo$RNAiID,' ','_'))] %in% 
                                                    iCELnames$WormBase_Gene_ID[match(vhalist$vha_in_model, iCELnames$ICELgene)]),]
# exclude RNA pol and DNA pol (eq. to that ribosome is not in the model as well)
pollist = read.csv('./../CR_model_final_run/input/pol_in_model.csv')
inputTb_metResponsiove = inputTb_metResponsiove[!(inputTb_metResponsiove$WBID %in% 
                                                    iCELnames$WormBase_Gene_ID[match(pollist$pol_in_model, iCELnames$ICELgene)]),]
inputTb_metResponsiove = inputTb_metResponsiove[!(conditionInfo$RNAi_WBID[match(inputTb_metResponsiove$condID, str_replace(conditionInfo$RNAiID,' ','_'))] %in% 
                                                    iCELnames$WormBase_Gene_ID[match(pollist$pol_in_model, iCELnames$ICELgene)]),]



# load the network distance 
rxnGeneMat=read.csv("./../../4_networkAnalysis/input/rxnGeneMat_iCEL1314.csv",row.names = 1)
oriDist=read.csv("./../../4_networkAnalysis/input/distanceMatrix_plain.txt",row.names = 1,sep = '\t')
wtdDist=read.csv("./../../4_networkAnalysis/input/distanceMatrix_weighted.txt",row.names = 1, sep = '\t')
distMat_regular <- pmin(as.matrix(oriDist), t(as.matrix(oriDist)),na.rm = T)
distMat_weighted <- pmin(as.matrix(wtdDist), t(as.matrix(wtdDist)),na.rm = T) # some mystery NaN, ignore
# convert to the gene-gene distance (using minimal by association)
inputTb_metResponsiove$RNAi_WBID = conditionInfo$RNAi_WBID[match(inputTb_metResponsiove$condID, str_replace(conditionInfo$RNAiID,' ','_'))]
geneNames = intersect(colnames(rxnGeneMat), c(inputTb_metResponsiove$RNAi_WBID, inputTb_metResponsiove$WBID))
ggDistMat_regular = matrix(NA, nrow = length(geneNames), ncol = length(geneNames))
ggDistMat_weighted = matrix(NA, nrow = length(geneNames), ncol = length(geneNames))

for (i in 1:ncol(ggDistMat_regular)){
  for (j in 1:ncol(ggDistMat_regular)){
    gene1 = geneNames[i]
    gene2 = geneNames[j]
    rxns1 = c(paste(rownames(rxnGeneMat)[rxnGeneMat[,gene1] == 1],'f',sep = ''),
              paste(rownames(rxnGeneMat)[rxnGeneMat[,gene1] == 1],'r',sep = ''))
    rxns2 = c(paste(rownames(rxnGeneMat)[rxnGeneMat[,gene2] == 1],'f',sep = ''),
              paste(rownames(rxnGeneMat)[rxnGeneMat[,gene2] == 1],'r',sep = ''))
    
    ggDistMat_regular[i,j] = min(distMat_regular[rownames(distMat_regular) %in% rxns1, colnames(distMat_regular) %in% rxns2, drop = F])
    ggDistMat_weighted[i,j] = min(distMat_weighted[rownames(distMat_weighted) %in% rxns1, colnames(distMat_weighted) %in% rxns2, drop = F])
    # the missing values (i.e., Reactions with only side metabolites that are not in the distance matrix will get inf in the min() function)
  }
  print(i)
}
rownames(ggDistMat_regular) = geneNames
colnames(ggDistMat_regular) = geneNames
rownames(ggDistMat_weighted) = geneNames
colnames(ggDistMat_weighted) = geneNames

for (i in 1:nrow(inputTb_metResponsiove)){
  inputTb_metResponsiove$regular_dist[i] = ggDistMat_regular[inputTb_metResponsiove$WBID[i], inputTb_metResponsiove$RNAi_WBID[i]]
  inputTb_metResponsiove$weighted_dist[i] = ggDistMat_weighted[inputTb_metResponsiove$WBID[i], inputTb_metResponsiove$RNAi_WBID[i]]
  
}

pdf('figures/show_case_network_distance_for_responses.pdf', height = 5, width = 5)
h1 = hist(ggDistMat_regular,probability = T, col = rgb(0.5,0.5,0.5, alpha = 0.5), ylim = c(0,0.25),main = '',
          xlab = 'Metabolic network distance')
hist(inputTb_metResponsiove$regular_dist,probability = T, add = T, col = rgb(1,0,0, alpha = 0.5), breaks = h1$breaks)
hist(inputTb_metResponsiove$regular_dist[inputTb_metResponsiove$condID == 'x.hpo_8_met2_lib3'],probability = T, add = T, col = rgb(0,0,1, alpha = 0.5), breaks = h1$breaks)
legend(8, 0.25, c('all','iCEL GRN interactions','hpo-8 RNAi interactions'), 
       fill = c(rgb(0.5,0.5,0.5, alpha = 0.5),rgb(1,0,0, alpha = 0.5),rgb(0,0,1, alpha = 0.5)))
dev.off()

# hist(ggDistMat_weighted,probability = T, col = rgb(0.5,0.5,0.5, alpha = 0.5))
# hist(inputTb_metResponsiove$weighted_dist,probability = T, add = T, col = rgb(1,0,0, alpha = 0.5))



