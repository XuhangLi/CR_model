# 02062023
# note: these responsive conditions may be redefined with losen cutoff for meta-analysis 
# we also tested CR model with energy not repressing anything and only activates energy; the result is slightly less significant than the current CR model
# we didnt test energy repressed others since it is not logical

# in general, this 150% increase from baseline is quite representative. Even if we check for the individual edge randomization and 
# for these big edges (the 'C' edges), it is still similar level of significance (~1.5 fold baseline). So there is no magic small 
# baseline unless we go for super simplified model (such as only five 'C'). 

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
pheatmap::pheatmap(tmp[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')], labels_row = iCELnames$ICELgene[match(rownames(tmp), iCELnames$WormBase_Gene_ID)])
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
RNAiclusters = read.csv('./../../../../2_DE/output/RNAi_groups.csv')
sum(RNAiclusters$clusters[(str_replace(RNAiclusters$RNAiID,' ','_') %in% uniConds)] == -1)/sum(RNAiclusters$clusters==-1)
# although the overall coverage of the RNAi conditions is at similar level (75% vs 78%), the FBA modeling is unbiased 
# so it covers 44% of unclustered conditions and it assesses all iCEL DEG instead of just selected coexpression clusters
# so this justifies it as a systems-level validation
length(intersect(rownames(classMat)[rowSums(classMat[,c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0], inputTb_metResponsiove$WBID))/(length(unique(inputTb_metResponsiove$WBID)))
# it covers 66% of DEG space 


# try to make plots for the 54% and 23% seperately (we also include the 44% in the 23% analysis)
total_condition_analyzed = c()

# to visualize the model, we define the multi-obj that contains the target obj as the target obj, others 
# rank by their total counts of calls and defined as the highest count to lowest (winer takes all)
# also try simply use this winer takes all rule for all ranks, so the plot is not biased by hypothesis
# first make plots for the 50% unique classified genes
obj_perturb = 'lipid'
# obj_response = 'lipid' # only used
# obj_direction = 'up'
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



myCond = 'x.gspd_1_met7_lib2'
DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
classMat_up = classMat[rownames(classMat) %in% DEGs_up, modelObj]
classMat_down = classMat[rownames(classMat) %in% DEGs_down, modelObj]

# manually sort the rows
sortRowsManual <- function(classMat_up){
  ct = colSums(classMat_up)
  classMat_up = classMat_up[,names(ct)[order(ct,decreasing = T)]]
  rowOrder = c()
  for (j in 1:5){
    for (i in 5:j){
      if (j >1){
        if(i<5){
          rowOrder = c(rowOrder, rownames(classMat_up)[rowAlls(classMat_up[,colnames(classMat_up)[j:i],drop = F] == 1) & 
                                                         rowAlls(classMat_up[,colnames(classMat_up)[(i+1):5],drop = F] == 0) & 
                                                         rowAlls(classMat_up[,colnames(classMat_up)[1:(j-1)],drop = F] == 0)])
        }else{
          rowOrder = c(rowOrder, rownames(classMat_up)[rowAlls(classMat_up[,colnames(classMat_up)[j:i],drop = F] == 1) & 
                                                         rowAlls(classMat_up[,colnames(classMat_up)[1:(j-1)],drop = F] == 0)])
        }
        
      }else{
        if(i<5){
          rowOrder = c(rowOrder, rownames(classMat_up)[rowAlls(classMat_up[,colnames(classMat_up)[j:i],drop = F] == 1) & 
                                                         rowAlls(classMat_up[,colnames(classMat_up)[(i+1):5],drop = F] == 0)])
        }else{
          rowOrder = c(rowOrder, rownames(classMat_up)[rowAlls(classMat_up[,colnames(classMat_up)[j:i],drop = F] == 1)])
        }
        }
      
    }
  }
  rowOrder = c(rowOrder, rownames(classMat_up)[rowAlls(classMat_up[,colnames(classMat_up)[1:5],drop = F] == 0)])
  
  rowOrder_final = c()
  for (i in 5:0){
    tmp = rownames(classMat_up)[rowSums(classMat_up) == i]
    if (length(tmp) > 0){
      tmp = c(rowOrder[rowOrder %in% tmp], setdiff(tmp, rowOrder))
    }
    rowOrder_final = c(rowOrder_final, tmp)
  }
  
  classMat_up = classMat_up[rowOrder_final, ]
  return(classMat_up)
}

library(pheatmap)
dev.off()

pdf('figures/FBA_classification_example_gspd_1.pdf')
classMat_up = sortRowsManual(classMat_up)
classMat_up = classMat_up[,c("lipid", "pro_modi","nucl_acid", "energy",  "pro_syn")]
pheatmap(classMat_up,cluster_rows = F, cluster_cols = F,main = 'UP')

classMat_down = sortRowsManual(classMat_down)
classMat_down = classMat_down[c(1:3,7:17,4:5,18:20,6,21:27),]
pheatmap(classMat_down,cluster_rows = F, cluster_cols = F,main = 'DOWN')
dev.off()



    
    