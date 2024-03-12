# testing of CR model based on FBA

# this analysis is very similar to the visualization, however, with a simpler setting. We directly use the core function scores from FBA to define 
# the affected and unaffected core objectives under a gene perturbation; and accordingly, we analyze how many DEGs for the same core function is up
# regulated and how many for others are down regulated. This is a direct test of CR model based on FBA simulation and WPS data. The only parameter in
# this analysis the the core function score cutoff, which we assume 0.001 here and will test its robustness in another script. We name this type of 
# CR model analysis as "fused with FBA" (internally here).

# notes: about pros and cons for FBA-fused CR model
# in this case, most energy genes will be classified to all five core functions, as they underlie other core functions in FBA. 
# This case is not good for visualization (too ambiguous - is a five-obj gene an energy, or lipid, or else?). It has the same problem
# for edge quantification (same thing - what is the edge for a five-obj gene? each edge will get equal quantification so it is not interpretable)
# Conversely, it is totally fine for testing CR model: in terms of CR model, there is no concept of gene class - it is only about what objective 
# is perturbed by the gene perturbation, and does it correlate which that of the DEGs. More uniquely classified genes are needed when we wanted to
# gain insight of core-function crosstalk and human interpretation (such as visualizations).


library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)

# first load data -  same code copied from visualization except for loading the core function scores 

# load the met gene classifications
scoreMat = read.csv('output/delta_flux_scoreMat.csv',row.names = 1)
classMat = as.data.frame(1*(scoreMat > 1e-3)) # cutoff at 1e-3, the robustness of this cutoff will be tested later

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

# remove RNAi target genes as DEGs
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
conditionInfo = read.csv('input/RNAi_condition_metaInfo.csv',row.names = 1)
# only analyze the iCEL responsive
inputTb_metResponsiove = inputTb[inputTb$RNAiID %in% conditionInfo$RNAiID[conditionInfo$isICEL & conditionInfo$isResponsive], ]
inputTb_metResponsiove=inputTb_metResponsiove[,c("WBID","RNAi","log2FoldChange_raw","condID")]

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


# check for some numbers
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


# compare the CR model expectations with real data

# NOTE: to keep this modeling simple, the class selection by DE similarity is not performed in this step; the rewiring 
# model calculation is fully based on FBA classification and the multi-class genes were used literally (assuming it affects multiple objectives)
# and the unclassified conditions will be left out from the analysis. This also guarantees the same set of RNAi was analyzed in every randomization 

# we get the conditions to analyze (icel_responsive (at least 2 up or 2 down DEG) and classified)
total_condition_analyzed = c()
for (i in 1:length(icel_resp)){
  if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
    total_condition_analyzed = c(total_condition_analyzed, icel_resp[i])
  }
}
# we have 223 conditions whose RNAi is associated with core functions so could be analyzed 

# analyze the five core functions
modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
model_explained$UP_yes = 0
model_explained$UP_no = 0
model_explained$DOWN_yes = 0
model_explained$DOWN_no = 0
# loop through each condition
for (condInd in 1:length(total_condition_analyzed)){
    myCond = total_condition_analyzed[condInd]
    # gets up and down DEGs
    DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
    DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
    
    # determine affected core functions
    labeledObj = colnames(classMat)[classMat[conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') == myCond],]==1]
    affectedObj = intersect(labeledObj, modelObj)
    
    # compare CR model expectations with the core functions of DEGs
    # get core functions of up genes
    subClassMat = classMat[rownames(classMat) %in% DEGs_up, modelObj]
    # calculate the DEG explained by model
    model_explained$UP_yes[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) > 0) # up and also affect the affected core function
    model_explained$UP_no[condInd] = sum(rowSums(subClassMat[,affectedObj,drop = F]) == 0) # others
    # same thing for down genes
    subClassMat = classMat[rownames(classMat) %in% DEGs_down, modelObj]
    # calculate the DEG explained by model
    model_explained$DOWN_yes[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) > 0) # down and affect other (unaffected) core functions
    model_explained$DOWN_no[condInd] = sum(rowSums(subClassMat[,setdiff(modelObj, affectedObj),drop = F]) == 0)
    # that's it - as simple as possible
}

# plot the result (in a bar plot of percentage explained)
b = model_explained$RNAi_cond
# sort by the total explained rate
tmp = model_explained[,2:5]
tmp = tmp / rowSums(tmp)
rewire_rate = rowSums(tmp[,c(1,3)])
rewire_rate[is.na(rewire_rate)] = 0
b[rank(rewire_rate,ties.method ='first')] =b
model_explained$RNAi_cond = factor(model_explained$RNAi_cond,levels = b)
# reformating the data frame for plotting with ggplot2
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

pdf(paste('figures/FBA_DEG_explained_distribution_FBA_fused_model.pdf',sep = ''),
    height = 7, width = 14)
print(ggarrange(p2, p1, nrow = 1))
dev.off()


# the overall average explained rate (up and down together) is 
obs_rate = mean(rewire_rate)
obs_total_DE_rate = sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5]) # 58%

# the overall average explained rate (up and down separately) is 
# up
tmp = model_explained[,2:3]
tmp = tmp / rowSums(tmp)
obs_rate_up = mean(tmp$UP_yes, na.rm = T)
obs_total_DE_rate_up = sum(model_explained[,c(2)]) / sum(model_explained[,2:3]) # 61%
# down
tmp = model_explained[,4:5]
tmp = tmp / rowSums(tmp)
obs_rate_down = mean(tmp$DOWN_yes, na.rm = T)
obs_total_DE_rate_down = sum(model_explained[,c(4)]) / sum(model_explained[,4:5]) # 50%



# next, perform randomization of core objective function associations to assess the significance 

# we do 10,000 randomizations
nRand = 10000
set.seed(1126)
rand_rate = c()
rand_rate2 = c()
rand_rate_up = c()
rand_rate2_up = c()
rand_rate_down = c()
rand_rate2_down = c()

# this is a good treatment to keep it a conservative randomization: 
# (1) we keep the number of labels for each obj (see below by row shuffling)
# (2) we keep the analyzed conditions assigned to at least one label - so all the WPS conditions will be analyzed in randomization
mustHasClass = unique(conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% total_condition_analyzed])
justRandom = setdiff(rownames(classMat), mustHasClass)
hasClassInd = as.numeric(which(rowSums(classMat[,modelObj])>0))

# loop through randomizations
for (nn in 1:nRand){
  
  # generate the random classification matrix
  # we shuffle the rowname (gene) labels to randomize the rows of the matrix 
  new_gene_names = rep(NA, nrow(classMat))
  new_gene_names[sample(hasClassInd, length(mustHasClass))] = mustHasClass
  new_gene_names[is.na(new_gene_names)] = justRandom[sample(length(justRandom))]
  
  classMat_rand = classMat
  rownames(classMat_rand) = new_gene_names
  
  # recalculate the number of DEGs explained by CR model
  model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
  model_explained$UP_yes = 0
  model_explained$UP_no = 0
  model_explained$DOWN_yes = 0
  model_explained$DOWN_no = 0
  
  for (condInd in 1:length(total_condition_analyzed)){
    myCond = total_condition_analyzed[condInd]
    DEGs_up = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw > 0]
    DEGs_down = inputTb_metResponsiove$WBID[inputTb_metResponsiove$condID == myCond & inputTb_metResponsiove$log2FoldChange_raw < 0]
    
    # determine affected obj
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

# save the randomization result and plot the data
save(file = 'randomization_result_FBA_fused_model.Rdata',rand_rate, rand_rate2, rand_rate_up, rand_rate2_up, rand_rate_down, rand_rate2_down)
hist(rand_rate)
hist(rand_rate2)
pdf(paste('figures/overall_DEG_explained_randomization_FBA_fused_model.pdf',sep = ''),width = 7,height = 6)
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

#############
# we give an visualization of how much the core function associations can be overlapped between genes if we only classify by a simple cutoff on the scores
gene_list = list(
  iCEL = rownames(classMat),
  energy = rownames(classMat)[classMat$energy==1],
  lipid = rownames(classMat)[classMat$lipid==1],
  pro_modi = rownames(classMat)[classMat$pro_modi==1],
  pro_syn = rownames(classMat)[classMat$pro_syn==1],
  nucl_acid = rownames(classMat)[classMat$nucl_acid==1]
)

library(eulerr)
set.seed(1126)
fit2 <- euler(classMat[,modelObj], 
              shape = "ellipse",
              #loss = 'region',
              #loss_aggregator = 'max', 
              control = list(extraopt = TRUE, extraopt_threshold = 0)
              )
# plot(fit2, quantities = TRUE, fill = c('#FF6701','#77AC30','#EDB120','#0072BD','#4DBEEE'),main = paste('unclassified =', sum(rowSums(classMat[,modelObj]==1)==0),'(',100*round(sum(rowSums(classMat[,modelObj]==1)==0)/nrow(classMat),2),'% )'))
# 
# 
# 
pdf(paste('figures/euler_plot_iCEL_gene_classification_only_by_thresholding.pdf',sep = ''),width = 7,height = 7)
plot(fit2, quantities = TRUE, fill = c('#FF6701','#77AC30','#EDB120','#0072BD','#4DBEEE'),main = paste('unclassified =', sum(rowSums(classMat[,modelObj]==1)==0),'(',100*round(sum(rowSums(classMat[,modelObj]==1)==0)/nrow(classMat),2),'% )'))
dev.off()

# this is not published - we didn't check if the numbers are correct (because some unfitted area will miss valid numbers)
