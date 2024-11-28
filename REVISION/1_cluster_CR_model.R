# Summary
# Visualize the DEGs based on their core functions assigned by FBA simulation. 

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)

# load the met gene classifications
classMat = read.csv('./../CR_model_final_run/output/FBA_classification_matrix.csv',row.names = 1)
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
# try to organize the perturbations with clustering 
proportionMat = c()
# get all genes that can be analyzed 
total_condition_analyzed = c()
for (i in 1:length(icel_resp)){
  if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
    total_condition_analyzed = c(total_condition_analyzed, icel_resp[i])
  }
}

# visualize all 223 conditions
# get corresponding conditions
myconds = total_condition_analyzed
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

# recover the matrix for clustering 
library(reshape2)
# we use the up gene profiles to cluster to enhance the pattern of specific compensation
prop_up <- acast(stats_long_merge[,c('RNAi_cond', 'variable', 'up')], RNAi_cond ~ variable, value.var = "up")
prop_up = t(scale(t(prop_up), scale = rowSums(prop_up), center = F))
prop_up[is.na(prop_up)] = 0 # avoid errors

# this clustering is critical, a little tuning here
cl_method = 'centroid'
distMat = dist(prop_up)

# library(lsa)
# tmp = cosine(t(prop_up))
# tmp[is.na(tmp)] = 1
# distMat = as.dist(1-tmp)

cl = hclust(distMat, method = cl_method)
roworders = cl$labels[cl$order]
# plot(cl)
stats_long_merge$RNAi_cond = factor(as.character(stats_long_merge$RNAi_cond), levels = roworders)

# manually assign the core function order for better visual interpretation 
stats_long_merge$variable = factor(as.character(stats_long_merge$variable), 
                                   levels = c("noClass",
                                              "pro_syn",
                                              "nucl_acid",
                                              "pro_modi",
                                              "energy",
                                              "lipid"))


# make the good labels
stats_long_merge$labels = str_remove(stats_long_merge$RNAi_cond,'^x.')
stats_long_merge$labels = str_remove(stats_long_merge$labels, '_met[0-9]+_lib.$')
stats_long_merge$labels = str_replace(stats_long_merge$labels,'_','-')

# plot up and down with the same condition order for later merging in graphics
p1 <- ggplot(stats_long_merge, aes(x=RNAi_cond, y = up,fill=variable)) +
  geom_bar(stat="identity", position="fill",colour="black", size = 0.1) + 
  scale_y_continuous(labels = scales::percent,expand = c(0,0))+
  theme_bw()+
  theme(axis.text.x = element_text(colour ='black'),axis.title.x=element_text(colour ='black'),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),axis.title.y=element_blank(),
        legend.position = "none", # skip to align figures
        #legend.title = element_blank(),
        #legend.text = element_text(colour ='black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  coord_flip()+
  scale_fill_manual("legend", values = c("energy" = "#FF6701", 
                                         "lipid" = "#77AC30",
                                         "pro_modi" = "#EDB120",
                                         "nucl_acid" = "#4DBEEE",
                                         "pro_syn" = "#0072BD",
                                         "noClass" = "white")) +
  labs(y = 'Up DEG')+
  geom_text(aes(label = up), position = position_fill(vjust = 0.5),size = 1.5)

# p1

# plot down with obj reordered 
stats_long_merge$variable = factor(as.character(stats_long_merge$variable), 
                                   levels = c("noClass",
                                              "pro_syn",
                                              "nucl_acid",
                                              "pro_modi",
                                              "lipid",
                                              "energy"))
p2 <- ggplot(stats_long_merge, aes(x=RNAi_cond, y = down,fill=variable)) +
  geom_bar(stat="identity", position="fill",colour="black", size = 0.1) + 
  scale_y_continuous(labels = scales::percent,expand = c(0,0))+
  scale_x_discrete(labels = stats_long_merge$labels[match(roworders, stats_long_merge$RNAi_cond)]) + 
  theme_bw()+
  theme(axis.text.x = element_text(colour ='black'),axis.title.x=element_text(colour ='black'),
        axis.text.y = element_text(colour ='black', size = 3),axis.title.y=element_blank(),
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
  geom_text(aes(label = -down), position = position_fill(vjust = 0.5),size = 1.5)

# p2

# create the RNAi core function matrix 
library(ComplexHeatmap)
library(circlize)
genenames = conditionInfo$RNAi_WBID[match(roworders[length(roworders):1], str_replace(conditionInfo$RNAiID, ' ','_'))]
heatMat = classMat[genenames, c('pro_syn','lipid','energy','pro_modi','nucl_acid')]
heatMat[heatMat[,1]==1,1] = 1
heatMat[heatMat[,2]==1,2] = 2
heatMat[heatMat[,3]==1,3] = 3
heatMat[heatMat[,4]==1,4] = 4
heatMat[heatMat[,5]==1,5] = 5

col_fun <- colorRamp2(c(0, 1, 2, 3, 4, 5), c("white", "#0072BD","#77AC30","#FF6701","#EDB120","#4DBEEE"))

colnames(heatMat) = c("pro", "Lp", "Eg", "ecm", "NA")
heatmap <- Heatmap(heatMat,cluster_rows = F, cluster_columns = F, col = col_fun, show_row_names = F, show_column_names = T,show_heatmap_legend = F,
                   rect_gp = gpar(col = 'black', lwd = 0.5))

# Draw the heatmap to a grob
heatmap_grob <- grid.grabExpr(draw(heatmap))



# Create the dendrogram plot
library(ggdendro)
dendro_data <- dendro_data(as.dendrogram(cl), type = "rectangle")
dendro_data$segments$y <- max(dendro_data$segments$y) - dendro_data$segments$y
dendro_data$segments$yend <- max(dendro_data$segments$yend) - dendro_data$segments$yend

dendro_plot <- ggplot(dendro_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  theme_void() + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"))



library(cowplot)
pdf('figures/CR_model_barplots_clustered.pdf', height = 9, width = 16)
print(ggarrange(dendro_plot, p2, p1, ggdraw() + draw_grob(heatmap_grob), nrow = 1, widths = c(1,4,4,0.4)))
dev.off()


# it turns out that 160 conditions were analyzed (unique obj association)
