# testing FBA-fused CR model with GRN randomization

# in addition to randomizing core function associations, we can also randomize the GRN to test if it is significant. 

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)

# load data for testing the FBA-fused model
# load the met gene classifications
scoreMat = read.csv('output/delta_flux_scoreMat.csv',row.names = 1)
classMat = as.data.frame(1*(scoreMat > 1e-3))

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

# remove DEGs related to RNAi targeted genes
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
modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')

# we get the conditions to analyze (icel_responsive (at least 2 up or 2 down DEG) and classified)
total_condition_analyzed = c()
for (i in 1:length(icel_resp)){
  if (rowSums(classMat[conditionInfo$RNAi_WBID[match(icel_resp[i], str_replace(conditionInfo$RNAiID,' ','_'))],c('energy','lipid','pro_modi','pro_syn','nucl_acid')]) > 0){
    total_condition_analyzed = c(total_condition_analyzed, icel_resp[i])
  }
}

modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
# 01052024: update to use the FBA-fused version

# check CR-model explanation 
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

print(ggarrange(p2, p1, nrow = 1))

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


# run the GRN randomization to access significance 


# first, generate the random iCEL GRN
# To focus on the edge wiring specificity instead any other possible confounding factors, we keep the GRN to be 
# randomized within the sub-GRN in this specific analysis, which means iCEL genes and the analyzed conditions

# let's generate the random GRN first
library(igraph)
# define the GRN using igraph
oriNet = inputTb_metResponsiove[inputTb_metResponsiove$condID %in% total_condition_analyzed, c("condID","WBID")]
oriNet_FC = inputTb_metResponsiove[inputTb_metResponsiove$condID %in% total_condition_analyzed, ]
oriNet_FC$RNAi_geneName = conditionInfo$RNAi_geneName[match(oriNet_FC$condID, str_replace(conditionInfo$RNAiID,' ','_'))]
oriNet_FC$RNAi_WBID = conditionInfo$RNAi_WBID[match(oriNet_FC$condID, str_replace(conditionInfo$RNAiID,' ','_'))]
allgenes = unique(oriNet_FC$WBID)
colnames(oriNet) = c('from','to')
g_ori <- graph_from_data_frame(oriNet, directed=TRUE)
# randomize using edge swapping
set.seed(1030)
for (j in 1:10){ # generate final data later 
  sampleSet = list()
  for (i in 1:1000){
    # randomly rewiring (swap edges) the network for 50-100x of the number of edges in the original network
    randNet = rewire(g_ori,with = keeping_degseq(niter = sample(50:100,1) * nrow(oriNet)))

    # convert randNet back to a DE table for analysis 
    randNet = as_data_frame(randNet)
    colnames(randNet) = c('RNAiID','WBID')
    randNet$RNAi_geneName = oriNet_FC$RNAi_geneName[match(randNet$RNAiID, oriNet_FC$condID)]
    randNet$RNAi_WBID = oriNet_FC$RNAi_WBID[match(randNet$RNAiID, oriNet_FC$condID)]
    
    # fix the RNAi-target edges (they should not exist; if exits, rewire to something else again)
    nSelf = which(randNet$WBID == randNet$RNAi_WBID)
    while(length(nSelf)>0){
      for (k in 1:length(nSelf)){
        rewireTo = 0
        while (rewireTo == 0) {
          tryTo = sample(nrow(randNet),1)
          # 05252023: fix minor bug in the criteria to stringently avoid self targeting edges and avoid bias in the swapping (the original is already close to perfect)
          if (!(randNet$WBID[tryTo] %in% randNet$WBID[randNet$RNAiID == randNet$RNAiID[nSelf[k]]]) & # not existing edge
              !(randNet$WBID[nSelf[k]] %in% randNet$WBID[randNet$RNAiID == randNet$RNAiID[tryTo]]) & # not existing edge
              randNet$WBID[tryTo] != randNet$RNAi_WBID[nSelf[k]] & # not the RNAi target
              randNet$WBID[nSelf[k]] != randNet$RNAi_WBID[tryTo]
          ){
            rewireTo = tryTo
          }
        }
        # rewire
        tmp = randNet$WBID[nSelf[k]]
        randNet$WBID[nSelf[k]] = randNet$WBID[rewireTo]
        randNet$WBID[rewireTo] = tmp
      }
      nSelf = which(randNet$WBID == randNet$RNAi_WBID) # this guarantee the original still has no self
    }
    
    
    # finally assign the rewired FC (fold changes)
    # since the rewire function in igraph does not support attributes, we have to manually add back the FC information
    # as the network was randomly shuffled, it doesnt matter if FC attributes goes with the edge in rewiring; the key 
    # is that the distribution of FC for each gene should be maintained. so we randomly assign the original FC vector for 
    # each gene into the new random network (the FC distribution of each sample does not matter as we dont expect it to 
    # have any distribution)
    randNet$log2FoldChange_raw = NA
    for (z in 1:length(allgenes)){
      randNet$log2FoldChange_raw[randNet$WBID == allgenes[z]] = 
        oriNet_FC$log2FoldChange_raw[oriNet_FC$WBID == allgenes[z]][sample(sum(oriNet_FC$WBID == allgenes[z]))]
    }
    
    
    sampleSet[[i]] = randNet
    
    # these are codes for sanity check of the random network 
    # a = table(randNet$RNAiID)
    # plot(as.numeric(a[names(N_DE)]), as.numeric(N_DE))
    # all(as.numeric(a[names(N_DE)]) == as.numeric(N_DE))
    # a = table(randNet$WBID)
    # geneCounts = table(inputTb_Responsiove$WBID)
    # plot(as.numeric(a[names(geneCounts)]), as.numeric(geneCounts))
    # all(as.numeric(a[names(geneCounts)])== as.numeric(geneCounts))
    # randNet$pair = paste(randNet$WBID, randNet$RNAiID)
    # length(unique(randNet$pair)) == nrow(randNet)
    # oriPairs = paste(oriNet$to,oriNet$from)
    # length(intersect(oriPairs,randNet$pair)) / nrow(randNet)
    # hist(randNet$log2FoldChange_raw[randNet$WBID == 'WBGene00007215'])
    # hist(inputTb_Responsiove$log2FoldChange_raw[inputTb_Responsiove$WBID == 'WBGene00007215'])
    
    print(i)
  }
  randNets = sampleSet
  save(randNets,file = paste('randomGRN/randNet1000_ES_batch_',j,'.Rdata',sep = ''))
}



# now, let's calculable the CR model explanation using the random networks
nRand = 10
set.seed(1126)
rand_rate = c()
rand_rate2 = c()
rand_rate_up = c()
rand_rate2_up = c()
rand_rate_down = c()
rand_rate2_down = c()

randomizeFBA = FALSE # whether or not to randomize the gene classification simultaneously 

# this is for randomizeFBA = TRUE
# we keep the number of labels for each obj and keep the analyzed conditions assigned to at least one label
mustHasClass = unique(conditionInfo$RNAi_WBID[str_replace(conditionInfo$RNAiID,' ','_') %in% total_condition_analyzed])
justRandom = setdiff(rownames(classMat), mustHasClass)
hasClassInd = as.numeric(which(rowSums(classMat[,modelObj])>0))

for (nn in 1:nRand){
  load(paste('./randomGRN/randNet1000_ES_batch_',nn,'.Rdata',sep = ''))
  for (ii in 1:1000){
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
    inputTb_metResponsiove_rand = randNets[[ii]]

    model_explained = data.frame(RNAi_cond = total_condition_analyzed)  
    model_explained$UP_yes = 0
    model_explained$UP_no = 0
    model_explained$DOWN_yes = 0
    model_explained$DOWN_no = 0
    
    for (condInd in 1:length(total_condition_analyzed)){
      myCond = total_condition_analyzed[condInd]
      DEGs_up = inputTb_metResponsiove_rand$WBID[inputTb_metResponsiove_rand$RNAiID == myCond & inputTb_metResponsiove_rand$log2FoldChange_raw > 0]
      DEGs_down = inputTb_metResponsiove_rand$WBID[inputTb_metResponsiove_rand$RNAiID == myCond & inputTb_metResponsiove_rand$log2FoldChange_raw < 0]
      
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
    print(paste(nn,'...',ii))
  }
  # list(list(mean(rewire_rate_rand), sum(model_explained[,c(2,4)]) / sum(model_explained[,2:5])))
}

# save and plot the results 
if (randomizeFBA){
  save(file = 'randomization_result_FBA_fused_model_randGRN_randClassification.Rdata',rand_rate, rand_rate2, rand_rate_up, rand_rate2_up, rand_rate_down, rand_rate2_down)
}else{
  save(file = 'randomization_result_FBA_fused_model_randGRN.Rdata',rand_rate, rand_rate2, rand_rate_up, rand_rate2_up, rand_rate_down, rand_rate2_down)
}
hist(rand_rate)
hist(rand_rate2)
pdf(paste('figures/overall_DEG_explained_FBA_fused_randomization_randGRN_model.pdf',sep = ''),width = 7,height = 6)
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


# additional notes:

# when normalizing the GRN, the background is even higher but the variance is lower. This is because most iCEL DEGs are the 
# energy and lipid genes for up and energy genes for down. Therefore, it fits into the models for many conditions (especially
# for energy and lipid conditions), causing a generally high bg. The gain of explanation is only 10% but obviously significant 

# the FBA randomization has a lower baseline because it is a randomization at the network gene level (without a preference for 
# lipid or energy gene), so the permutation of lipid/energy bias in a way decreased the baseline.

# interestingly, we found dual randomization (not shown) is similar to the randomization of FBA labeling 



# 
# # ---- trouble shooting visualization --- 
# # exclude the multiple obj genes 
# classMat_bak = classMat
# modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
# classMat[rowSums(classMat[,modelObj]) > 1,] = 0 
# 
# # try to make plots for the 54% and 23% seperately (we also include the 44% in the 23% analysis)
# total_condition_analyzed = c()
# 
# # to visualize the model, we define the multi-obj that contains the target obj as the target obj, others 
# # rank by their total counts of calls and defined as the highest count to lowest (winer takes all)
# # also try simply use this winer takes all rule for all ranks, so the plot is not biased by hypothesis
# # first make plots for the 50% unique classified genes
# obj_perturb = 'lipid'
# # obj_response = 'lipid' # only used
# # obj_direction = 'up'
# modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
# pdfHeight = c("energy" = 12,
#               'lipid' = 7,
#               'pro_modi' = 5,
#               'pro_syn' = 3,
#               'nucl_acid' = 2)
# pdfWidth = c("energy" = 20,
#              'lipid' = 14,
#              'pro_modi' = 14,
#              'pro_syn' = 14,
#              'nucl_acid' = 14)
# 
# for (obj_perturb in modelObj){
#   if (obj_perturb == 'energy'){
#     # energy has to be only energy 
#     obj_genes = intersect(rownames(classMat)[classMat[,obj_perturb]==1 & rowSums(classMat[,modelObj]) <= 1],conditionInfo$RNAi_WBID)
#   }else{
#     # if it is energy + one biomass obj; it is considered as the biomass obj
#     # if it is more than one biomass obj, it is considered as multi-obj
#     obj_genes = intersect(rownames(classMat)[classMat[,obj_perturb]==1 & rowSums(classMat[,setdiff(modelObj,'energy')]) <= 1],conditionInfo$RNAi_WBID)
#   }
#   myconds = intersect(str_replace(conditionInfo$RNAiID[conditionInfo$RNAi_WBID %in% obj_genes],' ','_'), inputTb_metResponsiove_rand$RNAiID)
#   
#   for (obj_direction in c('up','down')){
#     stats = data.frame(RNAi_cond = myconds)
#     stats$energy = 0
#     stats$lipid = 0
#     stats$pro_modi = 0
#     stats$pro_syn = 0 
#     stats$nucl_acid = 0 
#     stats$noClass = 0
#     for (i in 1:nrow(stats)){
#       myCond = stats$RNAi_cond[i]
#       if (obj_direction == 'up'){
#         DEGs = inputTb_metResponsiove_rand$WBID[inputTb_metResponsiove_rand$RNAiID == myCond & inputTb_metResponsiove_rand$log2FoldChange_raw > 0]
#       }else{
#         DEGs = inputTb_metResponsiove_rand$WBID[inputTb_metResponsiove_rand$RNAiID == myCond & inputTb_metResponsiove_rand$log2FoldChange_raw < 0]
#       }
#       
#       if (any(rownames(classMat) %in% DEGs)){
#         # rank the obj candidates 
#         nCalls = colSums(classMat[rownames(classMat) %in% DEGs, modelObj])
#         orders = names(nCalls)[order(nCalls,decreasing = T)]
#         # let's use a simple strategy: assign the multi-obj to the most abundant applicable obj
#         subClassMat = classMat[rownames(classMat) %in% DEGs, modelObj]
#         for (j in 1:nrow(subClassMat)){
#           if (sum(subClassMat[j,])>1){
#             for (k in 1:length(orders)){
#               if (subClassMat[j,orders[k]] == 1){
#                 subClassMat[j,] = 0
#                 subClassMat[j,orders[k]] = 1
#                 break
#               }
#             }
#           }
#         }
#         stats[i,2:6] = colSums(subClassMat)
#         stats[i,7] = sum(rowSums(subClassMat) == 0)
#       }else{
#         stats[i,2:7] = 0
#       }
#     }
#     
#     # filter out low response conditions
#     stats = stats[rowSums(stats[,2:7])>1,] # we require at least two iCEL DE to enable meaningful interpretation (of proportion)
#     
#     # plot
#     b = stats$RNAi_cond
#     # sort by the most abundant class
#     tmp = stats[,2:7]
#     tmp = tmp / rowSums(tmp)
#     aveRewireRate = colMeans(tmp[,1:5])
#     sortBy = names(aveRewireRate)[order(aveRewireRate,decreasing = T)][1]
#     rewire_rate = stats[,sortBy] / rowSums(stats[,2:7])
#     rewire_rate[is.na(rewire_rate)] = 0
#     b[rank(rewire_rate,ties.method ='first')] =b
#     stats$RNAi_cond = factor(stats$RNAi_cond,levels = b)
#     stats_long = reshape2::melt(stats)
#     stats_long$variable = factor(as.character(stats_long$variable), levels = c('noClass',names(aveRewireRate)[order(aveRewireRate,decreasing = F)]))
#     if (obj_direction == 'up'){
#       stats_up = stats
#       sortBy_up = sortBy
#       stats_long_up = stats_long
#     }else{
#       stats_down = stats
#       sortBy_down = sortBy
#       stats_long_down = stats_long
#     }
#     
#   }
#   
#   # format the merged table
#   stats_long_down$value = -stats_long_down$value
#   extra_conds = levels(stats_long_down$RNAi_cond)[!(levels(stats_long_down$RNAi_cond) %in% levels(stats_long_up$RNAi_cond))]
#   new_levels = c(extra_conds, levels(stats_long_up$RNAi_cond))
#   stats_long_down$RNAi_cond = factor(as.character(stats_long_down$RNAi_cond), levels = new_levels)
#   stats_long_up$RNAi_cond = factor(as.character(stats_long_up$RNAi_cond), levels = new_levels)
#   
#   stats_long_merge = stats_long_up
#   stats_long_merge$up = stats_long_merge$value
#   stats_long_merge$down = 0
#   stats_long_down$up = 0
#   stats_long_down$down = stats_long_down$value
#   toMerge = c()
#   for (i in 1:nrow(stats_long_down)){
#     if (any(stats_long_merge$RNAi_cond == stats_long_down$RNAi_cond[i] & stats_long_merge$variable == stats_long_down$variable[i])){
#       stats_long_merge$down[stats_long_merge$RNAi_cond == stats_long_down$RNAi_cond[i] & stats_long_merge$variable == stats_long_down$variable[i]] = stats_long_down$value[i]
#     }else{
#       toMerge = c(toMerge, i)
#     }
#   }
#   stats_long_merge = rbind(stats_long_merge, stats_long_down[toMerge,])
#   # fill in missing values
#   # for (cond in levels(stats_long_merge$RNAi_cond)){
#   #   for (obj in levels(stats_long_merge$variable)){
#   #     if (!any(stats_long_merge$RNAi_cond == cond & stats_long_merge$variable == obj)){
#   #       stats_long_merge = rbind(stats_long_merge, c(cond, obj, 0,0,0))
#   #     }
#   #   }
#   # }
#   
#   # plot up and down with the same condition order for later merging in graphics
#   p1 <- ggplot(stats_long_merge, aes(x=RNAi_cond, y = up,fill=variable)) +
#     geom_bar(stat="identity", position="fill",colour="black", size = 0.1) + 
#     scale_y_continuous(labels = scales::percent,expand = c(0,0))+
#     theme_bw()+
#     theme(axis.text.x = element_text(colour ='black'),axis.title.x=element_text(colour ='black'),
#           axis.ticks.y=element_blank(),
#           axis.text.y = element_blank(),axis.title.y=element_blank(),
#           legend.title = element_blank(),
#           legend.text = element_text(colour ='black'),
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
#     coord_flip()+
#     scale_fill_manual("legend", values = c("energy" = "#FF6701", 
#                                            "lipid" = "#77AC30",
#                                            "pro_modi" = "#EDB120",
#                                            "nucl_acid" = "#4DBEEE",
#                                            "pro_syn" = "#0072BD",
#                                            "noClass" = "white")) +
#     labs(y = 'Up DEG')+
#     geom_text(aes(label = up), position = position_fill(vjust = 0.5),size = 3)+
#     geom_hline(yintercept = sum(classMat[,sortBy_up]) / nrow(classMat),color="white",
#                linetype="dashed")
#   
#   #p1
#   
#   # plot down with obj reordered 
#   stats_long_merge$variable = factor(as.character(stats_long_merge$variable), levels = levels(stats_long_down$variable))
#   p2 <- ggplot(stats_long_merge, aes(x=RNAi_cond, y = down,fill=variable)) +
#     geom_bar(stat="identity", position="fill",colour="black", size = 0.1) + 
#     scale_y_continuous(labels = scales::percent,expand = c(0,0))+
#     theme_bw()+
#     theme(axis.text.x = element_text(colour ='black'),axis.title.x=element_text(colour ='black'),
#           axis.text.y = element_text(colour ='black'),
#           legend.position = "none",
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
#     coord_flip()+
#     scale_fill_manual("legend", values = c("energy" = "#FF6701", 
#                                            "lipid" = "#77AC30",
#                                            "pro_modi" = "#EDB120",
#                                            "nucl_acid" = "#4DBEEE",
#                                            "pro_syn" = "#0072BD",
#                                            "noClass" = "white")) +
#     labs(y = 'Down DEG')+
#     geom_text(aes(label = -down), position = position_fill(vjust = 0.5),size = 3)+
#     geom_hline(yintercept = -sum(classMat[,sortBy]) / nrow(classMat),color="white",
#                linetype="dashed")
#   
#   #p2
#   
#   # pdf(paste('figures/FBA_classification_',obj_perturb,'_uniObjGenes.pdf',sep = ''),
#   #     height = pdfHeight[obj_perturb], width = pdfWidth[obj_perturb])
#   print(ggarrange(p2, p1, nrow = 1))
#   # dev.off()
#   total_condition_analyzed = c(total_condition_analyzed, as.character(unique(stats_long_merge$RNAi_cond)))
# }
# classMat = classMat_bak
