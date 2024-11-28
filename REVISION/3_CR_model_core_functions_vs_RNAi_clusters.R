library(stringr)

# load the met gene classifications
classMat = read.csv('./../CR_model_final_run/output/FBA_classification_matrix.csv',row.names = 1)
iCELnames = read.csv('./../../input_data/otherTbls/iCEL_IDtbl.csv')
rownames(classMat) = iCELnames$WormBase_Gene_ID[match(rownames(classMat),iCELnames$ICELgene)]
# mrpl-44 is a misannotated gene in the model; in the model it should be mccc-2 but mislabeled as mrpl-44 (which is mito ribo gene), since ribo genes are excluded from the analysis, we also exclude mrpl-44
classMat = classMat[-which(rownames(classMat) == 'WBGene00008514'),]

# load the RNAi clusters
RNAigroups = read.csv('./../../2_DE/clustering_python/HDBSCAN_RNAi_UMAP_densMap_embedding.csv',row.names = 1)
RNAinfo = read.csv('./../../2_DE/output/RNAi_condition_metaInfo.csv',row.names = 1)
RNAinfo = RNAinfo[(RNAinfo$isICEL|RNAinfo$isMetabolic) & RNAinfo$isResponsive,]
RNAinfo$conID = str_replace(RNAinfo$RNAiID,' ','_')
RNAinfo$clusters = RNAigroups$clusters[match(RNAinfo$conID, RNAigroups$conditions)]
RNAinfo$clusters[is.na(RNAinfo$clusters)] = -1
RNAigroups = RNAinfo
RNAi_group_names = read.csv('./../../2_DE/output/RNAi group annotation FINAL_APR21_2023.csv')
RNAinfo$RNAi_cluster_name = RNAi_group_names$summary[match(RNAinfo$clusters, RNAi_group_names$groupname)]

# subset to overlap 
RNAinfo = RNAinfo[RNAinfo$RNAi_WBID %in% rownames(classMat),]
classMat = classMat[rownames(classMat) %in% RNAinfo$RNAi_WBID,]

# we will compare the classifications for the 360 responsive iCEL genes 

# compare the overlap 
mygroup_CR = list(
  energy = rownames(classMat)[classMat$energy==1],
  lipid = rownames(classMat)[classMat$lipid==1],
  pro_modi = rownames(classMat)[classMat$pro_modi==1],
  pro_syn = rownames(classMat)[classMat$pro_syn==1],
  nucl_acid = rownames(classMat)[classMat$nucl_acid==1]
)
# it will be aributrary to match the categories; so we stay simple - focus on the large lipid cluster 
# lets pick a few 

mygroup_RNAi = list(
  lipid = RNAinfo$RNAi_WBID[RNAinfo$RNAi_cluster_name %in% c('Lipid synthesis & vacuolar ATPase')],
  energy = RNAinfo$RNAi_WBID[RNAinfo$RNAi_cluster_name %in% c('ETC cluster 1')],
  pro_modi = RNAinfo$RNAi_WBID[RNAinfo$RNAi_cluster_name %in% c('N-glycan synthesis')],
  pro_syn = RNAinfo$RNAi_WBID[RNAinfo$RNAi_cluster_name %in% c('AA-tRNA synthesis')]
)



library("VennDiagram")

pdf('figures/CR_core_function_vs_RNA_clusters.pdf',width = 7,height = 7)
g_names1 = c('Lipid genes by FBA','Energy genes by FBA',
             'ECM genes by FBA','Protein synthesis genes by FBA')
g_names2 = c('Lipid synthesis & vacuolar ATPase','ETC cluster 1','N-glycan synthesis','AA-tRNA synthesis')

for (i in 1:length(mygroup_RNAi)){
  
  gene_list = list(
    lipid_FBA = mygroup_CR[[names(mygroup_RNAi)[i]]],
    lipid_RNAi = mygroup_RNAi[[names(mygroup_RNAi)[i]]]
  )

  # Create and plot the Venn diagram with customization
  
  venn.plot <- venn.diagram(
    x = gene_list,
    category.names = c(g_names1[i], g_names2[i]),
    filename = NULL, # this will plot directly to R's graphics device
    output = TRUE,
    # Customizations
    fill = c("#FF9999", "#9999FF"), # Colors for the circles
    alpha = 0.5, # Transparency
    cat.col = c("#FF6666", "#6666FF"), # Colors for the category names
    cat.fontface = "bold", # Font style for the category names
    cat.cex = 1.5, # Font size for the category names
    cat.pos = c(-20, 20), # Position of the category names
    cat.dist = c(0.05, 0.05), # Distance of the category names from the circles
    fontface = "bold", # Font style for the numbers
    cex = 1.5, # Font size for the numbers
    margin = 0.1, # Margin around the diagram
    main.fontface = "bold", # Font style for the title
    main.cex = 2, # Font size for the title
    main.col = "black" # Color for the title
  )
  
  # Plot the Venn diagram
  grid.newpage()
  grid.draw(venn.plot)
}

dev.off()
