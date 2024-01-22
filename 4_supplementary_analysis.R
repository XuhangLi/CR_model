# (1) generate some numbers 
classMat = read.csv('output/FBA_classification_matrix.csv',row.names = 1)
iCELnames = read.csv('./../../input_data/otherTbls/iCEL_IDtbl.csv')
rownames(classMat) = iCELnames$WormBase_Gene_ID[match(rownames(classMat),iCELnames$ICELgene)]

# in this statistics, we dont do filtering so we keep mccc-2 mislabeled as is
modelObj = c('energy','lipid','pro_modi','pro_syn','nucl_acid')
classMat = classMat[,modelObj]

library(matrixStats)
# confirm 736 (+1 mccc-2/mrpl-44) assigned genes 
classMat = as.matrix(classMat)
sum(rowAnys(classMat))

# load WPS meta table
RNAiInfo = read.csv('./../../data_and_tables_to_publish/RNAi_condition_information.csv')
RNAiInfo$RNAi_WBID[which(RNAiInfo$RNAi_WBID == 'WBGene00008514/WBGene00303021')] = 'WBGene00008514'

# total tested genes is 548
sum(rownames(classMat)[rowAnys(classMat)] %in% RNAiInfo$RNAi_WBID)
# total responsive genes is 260
sum(rownames(classMat)[rowAnys(classMat)] %in% RNAiInfo$RNAi_WBID[RNAiInfo$isResponsive])
# total iCEL responsive is 361
sum(rownames(classMat) %in% RNAiInfo$RNAi_WBID[RNAiInfo$isResponsive])
# total tested iCEL is 891
sum(rownames(classMat) %in% RNAiInfo$RNAi_WBID)

# check enrichment 
phyper(260-1,361, 891-361, 548,lower.tail = F) # 5.723005e-08
phyper(260-1,548, 891-548, 361,lower.tail = F) # 5.723005e-08


# load DEG table
DEGinfo = read.csv('./../../data_and_tables_to_publish/DEG_clusters.csv')
DEtbl = read.csv('./../../data_and_tables_to_publish/DE_analysis_result_of_metabolic_WPS.csv')
DEtbl = DEtbl[DEtbl$is_in_mGRN,]
# we assume all genes in the genome is tested: if it never appears as DEG, it should be either never
# expressed or never differentially expressed.

# total unique DEG in iCEL is 977
sum(rownames(classMat) %in% DEtbl$DEG_WBID)
# total unique DEG genes assigned to objective is 613
sum(rownames(classMat)[rowAnys(classMat)] %in% DEtbl$DEG_WBID)
# total genes in iCEL is 1314
nrow(classMat)
# total objective-assigned genes in iCEL is 737
sum(rowAnys(classMat))

# # total DEG calls in iCEL is 14084
# sum(DEtbl$DEG_WBID %in% rownames(classMat))
# # total DEG calls assigned to objective is 8803
# sum(DEtbl$DEG_WBID %in% rownames(classMat)[rowAnys(classMat)])
# --> the proportion is similar 

# check enrichment 
phyper(613-1,977, 1314-977, 737,lower.tail = F) # 1.128663e-16


# the iCEL genes and conditions in mGRN without any filteirng 
DEtbl = read.csv('./../../data_and_tables_to_publish/DE_analysis_result_of_metabolic_WPS.csv')
DEtbl = DEtbl[DEtbl$is_in_mGRN,]
RNAiInfo = read.csv('./../../data_and_tables_to_publish/RNAi_condition_information.csv')
# only keep iCEL genes
classMat = read.csv('output/FBA_classification_matrix.csv',row.names = 1)
iCELnames = read.csv('./../../input_data/otherTbls/iCEL_IDtbl.csv')
rownames(classMat) = iCELnames$WormBase_Gene_ID[match(rownames(classMat),iCELnames$ICELgene)]
iCELgenes = rownames(classMat)
DEtbl = DEtbl[DEtbl$RNAi_WBID %in% iCELgenes & 
                DEtbl$DEG_WBID %in% iCELgenes, ]

# the total iCEL perturbation is 354
length(intersect(DEtbl$RNAiID, RNAiInfo$RNAiID[RNAiInfo$isICEL]))
# the total iCEL unique DEG is 977
length(intersect(DEtbl$DEG_WBID, iCELnames$WormBase_Gene_ID))



