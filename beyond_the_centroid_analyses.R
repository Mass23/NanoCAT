library(ggplot2)
library(tidyverse)
library(performance)
library(mgcv)
library(ggpubr)
library(vegan)
library(nlme)
library(propagate)
library(lmerTest)
library(dplyr)
library(cowplot)
library(emmeans)
setwd('~/Documents/MACE/NanoTax16S')

dir.create('stats')
dir.create('figures')

study_a_mocks = c("X10_Strain_Even_Mix_Genomic_Material_A","X10_Strain_Even_Mix_Genomic_Material_B","X10_Strain_Even_Mix_Genomic_Material_C","X10_Strain_Even_Mix_Genomic_Material_D") 
study_b_mocks = c("Full_Zymo")
study_e_mocks = c("Sample1","Sample10","Sample11","Sample12","Sample2","Sample3",
                  "Sample4","Sample5","Sample6","Sample7","Sample8","Sample9")
study_f_mocks = c("zym_zym_laboratory_sample")
study_g_mocks = c("MOCK_R9_1","MOCK_R9_2")

# Functions to extract the best threshold
find_alternat_thr = function(){
  all_res = data.frame()
  for (study_name in c('study_a', 'study_b', 'study_e', 'study_f', 'study_g')){
    otu_table = read.csv(paste0(c('all_res/otu_tab_', study_name, '.tsv'), collapse = ''), sep='\t')
    grg_table = read.csv(paste0(c('all_res/tax_greengenes_', study_name, '.tsv'), collapse = ''))
    slv_table = read.csv(paste0(c('all_res/tax_silva_', study_name, '.tsv'), collapse = ''))
    gdb_table = read.csv(paste0(c('all_res/tax_gtdb_', study_name, '.tsv'), collapse = ''))
    
    rownames(otu_table) = otu_table$X.OTU.ID
    otu_table$X.OTU.ID = NULL
    
    rownames(grg_table) = grg_table$Cluster
    grg_table$Cluster = NULL
    
    rownames(slv_table) = slv_table$Cluster
    slv_table$Cluster = NULL
    
    rownames(gdb_table) = gdb_table$Cluster
    gdb_table$Cluster = NULL
    
    for (thr in seq(1.01, 1.43, 0.01)){
      filtered_slv = filter_alternat(otu_table, slv_table, thr)
      filtered_grg = filter_alternat(otu_table, grg_table, thr)
      filtered_gdb = filter_alternat(otu_table, gdb_table, thr)
      
      bc_thr_slv = bc_thr_alternat(filtered_slv$otu_alternat, filtered_slv$tax_alternat, thr, study_name, 'Silva')
      bc_thr_grg = bc_thr_alternat(filtered_grg$otu_alternat, filtered_grg$tax_alternat, thr, study_name, 'Greengenes 2')
      bc_thr_gdb = bc_thr_alternat(filtered_gdb$otu_alternat, filtered_gdb$tax_alternat, thr, study_name, 'GTDB')
      
      all_res = rbind(all_res, bc_thr_slv)
      all_res = rbind(all_res, bc_thr_grg)
      all_res = rbind(all_res, bc_thr_gdb)
    }
  }
  return(all_res)
}

bc_thr_alternat = function(otu_table, tax_table, thr, study, db){
  mock_data = data.frame()
  
  if (study == 'study_a'){
    # Mock study A
    study_a_mocks = c("X10_Strain_Even_Mix_Genomic_Material_A","X10_Strain_Even_Mix_Genomic_Material_B","X10_Strain_Even_Mix_Genomic_Material_C","X10_Strain_Even_Mix_Genomic_Material_D")
    mock_genus_a_real = data.frame(row.names = c("Others","g__Bacillus","g__Bifidobacterium","g__Clostridium",
                                                 "g__Deinococcus","g__Enterococcus","g__Escherichia",
                                                 "g__Lactobacillus","g__Rhodobacter","g__Staphylococcus","g__Streptococcus"),
                                   real = c(0,rep(0.1,10)))
    
      mock_genus_a = summarise_alternat(otu_table, tax_table, study_a_mocks, db)
      beta_to_real = compare_genus_to_real(mock_genus_a, mock_genus_a_real)
      mock_data = rbind(mock_data, data.frame(Database = rep(db, 4),
                                                Study = rep('Study A', 4),
                                                Threshold = rep(thr, 4),
                                                Value = beta_to_real$bc_diss))}                       
  if (study == 'study_b'){
    # Mock study B
    study_b_mocks = c("Full_Zymo")
    mock_genus_b_real = data.frame(row.names = c("Others","g__Pseudomonas","g__Escherichia","g__Salmonella",
                                                 "g__Lactobacillus","g__Enterococcus","g__Staphylococcus",
                                                 "g__Listeria","g__Bacillus"),
                                   real = c(0,rep(1/8,8)))
    mock_genus_b = summarise_alternat(otu_table, tax_table, study_b_mocks, db)
    beta_to_real = compare_genus_to_real(mock_genus_b, mock_genus_b_real)
    mock_data = rbind(mock_data, data.frame(Database = rep(db, 1),
                                                Study = rep('Study B', 1),
                                                Threshold = rep(thr, 1),
                                                Value = beta_to_real$bc_diss))}   
  if (study == 'study_e'){
    # Mock study E
    study_e_mocks = c("Sample1","Sample10","Sample11","Sample12","Sample2","Sample3",
                      "Sample4","Sample5","Sample6","Sample7","Sample8","Sample9")
    mock_genus_e_real = data.frame(row.names = c("Others","g__Bacillus","g__Listeria",
                                                 "g__Staphylococcus","g__Enterococcus",
                                                 "g__Lactobacillus","g__Salmonella",
                                                 "g__Escherichia", "g__Pseudomonas"),
                                   real = c(0,0.174,0.141,0.155,0.099,0.184,0.104,0.101,0.042))
    
    mock_genus_e = summarise_alternat(otu_table, tax_table, study_e_mocks, db)
    beta_to_real = compare_genus_to_real(mock_genus_e, mock_genus_e_real)
    mock_data = rbind(mock_data, data.frame(Database = rep(db, 12),
                                              Study = rep('Study E', 12),
                                              Threshold = rep(thr, 12),
                                              Value = beta_to_real$bc_diss))}   
    
    if (study == 'study_f'){
    # Mock study F
    study_f_mocks = c("zym_zym_laboratory_sample")
    mock_genus_f_real = data.frame(row.names = c("Others","g__Pseudomonas","g__Escherichia","g__Salmonella",
                                                 "g__Lactobacillus","g__Enterococcus","g__Staphylococcus",
                                                 "g__Listeria","g__Bacillus"),
                                   real = c(0,rep(1/8,8)))
    
    mock_genus_f = summarise_alternat(otu_table, tax_table, study_f_mocks, db)
    beta_to_real = compare_genus_to_real(mock_genus_f, mock_genus_f_real)
    mock_data = rbind(mock_data, data.frame(Database = rep(db, 1),
                                              Study = rep('Study F', 1),
                                              Threshold = rep(thr, 1),
                                              Value = beta_to_real$bc_diss))}
    
    if (study == 'study_g'){
    # Mock study G
    study_g_mocks = c("MOCK_R9_1","MOCK_R9_2")
    mock_genus_g_real = data.frame(row.names = c("Others","g__Schaalia","g__Bifidobacterium","g__Cutibacterium",
                                                 "g__Phocaeicola","g__Porphyromonas","g__Deinococcus",
                                                 "g__Cereibacter","g__Neisseria","g__Acinetobacter",
                                                 "g__Pseudomonas","g__Escherichia","g__Helicobacter",
                                                 "g__Bacillus","g__Clostridium","g__Enterococcus",
                                                 "g__Lactobacillus","g__Staphylococcus","g__Streptococcus"),
                                   real = c(0,rep(1/20,16),2/20,2/20))
    
    mock_genus_g = summarise_alternat(otu_table, tax_table, study_g_mocks, db)
    beta_to_real = compare_genus_to_real(mock_genus_g, mock_genus_g_real)
    mock_data = rbind(mock_data, data.frame(Database = rep(db, 1),
                                              Study = rep('Study G', 1),
                                              Threshold = rep(thr, 1),
                                              Value = beta_to_real$bc_diss))}

  return(mock_data)
}

filter_alternat = function(otu_table, tax_table, alternat_thr){
  tax_table$Taxonomy.hyCAT = tax_table$Taxonomy.centroid
  mask_alternat = (tax_table$Confidence.most.confident > (tax_table$Confidence.centroid * alternat_thr)) & (tax_table$Taxonomy.most.confident != 'assigned')
  tax_table$Taxonomy.hyCAT[mask_alternat] = tax_table$Taxonomy.most.confident[mask_alternat]
  
  assigned_otus_alternat = rownames(tax_table)[tax_table$Taxonomy.hyCAT == 'assigned']
  bacteria_otus_alternat = rownames(tax_table)[startsWith(tax_table$Taxonomy.hyCAT, 'd__Bacteria')]
  tax_table_alternat = tax_table[grepl('d__Bacteria', tax_table$Taxonomy.hyCAT),]
  tax_table_alternat = tax_table_alternat[!(grepl('Mitochondria', tax_table_alternat$Taxonomy.hyCAT)),]
  tax_table_alternat = tax_table_alternat[!(grepl('Chloroplast', tax_table_alternat$Taxonomy.hyCAT)),]
  filtered_asvs_alternat = rownames(tax_table_alternat)
  
  otu_table_alternat_filtered = otu_table[rownames(otu_table) %in% filtered_asvs_alternat,]
  #otu_table_alternat_filtered[otu_table_alternat_filtered < 2] = 0
  otu_table_alternat_filtered = otu_table_alternat_filtered[rowSums(otu_table_alternat_filtered) > 1,]
  tax_table_alternat_filtered = tax_table_alternat[rownames(tax_table_alternat) %in% rownames(otu_table_alternat_filtered),]

  samples_to_remove = colnames(otu_table_alternat_filtered)[colSums(otu_table_alternat_filtered) < 400]

  otu_table_alternat_filtered = otu_table_alternat_filtered[,!(colnames(otu_table_alternat_filtered) %in% samples_to_remove)]
  tax_table_alternat_filtered = tax_table_alternat[rownames(tax_table_alternat) %in% rownames(otu_table_alternat_filtered),]
  
  return(list(tax_alternat=tax_table_alternat_filtered, otu_alternat=otu_table_alternat_filtered))
}

summarise_alternat <- function(otu_alternat, tax_alternat, cols_to_keep, db){
  out_tables = list()

  sep_db = ';'
  if (db == 'greengenes'){sep_db = '; '}
  
  tax_alternat = tax_alternat %>% separate(Taxonomy.hyCAT, into = c('domain','phylum','class','order','family','genus','species') , sep_db)
  otu_alternat$genus = map_chr(rownames(otu_alternat), function(x) tax_alternat$genus[rownames(tax_alternat) == x])
  otu_alternat$genus = gsub(' g__', 'g__', otu_alternat$genus)
  
  otu_alternat$genus[otu_alternat$genus == 'g__Limosilactobacillus'] = 'g__Lactobacillus' # Limosilactobacillus is part of Lactobacillus as of 2020
  otu_alternat$genus[otu_alternat$genus == 'g__Cereibacter'] = 'g__Rhodobacter' # Rhodobacter is also Cereibacter https://gtdb.ecogenomic.org/searches?s=al&q=g__Rhodobacter
  otu_alternat$genus[otu_alternat$genus == 'g__Clostridium_sensu_stricto_1'] = 'g__Clostridium'
  otu_alternat$genus[otu_alternat$genus == 'g__Escherichia-Shigella'] = 'g__Escherichia'
  otu_alternat$genus[otu_alternat$genus == 'g__Bacillus_A'] = 'g__Bacillus'
  otu_alternat$genus[otu_alternat$genus == 'g__Cereibacter_A'] = 'g__Rhodobacter'
  otu_alternat$genus[otu_alternat$genus == 'g__Bifidobacterium_388775'] = 'g__Bifidobacterium'
  otu_alternat$genus[otu_alternat$genus == 'g__Clostridium_T'] = 'g__Clostridium'
  otu_alternat$genus[otu_alternat$genus == 'g__Deinococcus_B'] = 'g__Deinococcus'
  otu_alternat$genus[otu_alternat$genus == 'g__Cereibacter_A_494292'] = 'g__Rhodobacter'
  otu_alternat$genus[otu_alternat$genus == 'g__Enterococcus_H_360604'] = 'g__Enterococcus'
  otu_alternat$genus[otu_alternat$genus == 'g__Bacillus_P_294101'] = 'g__Bacillus'
  otu_alternat$genus[otu_alternat$genus == 'g__Pseudomonas_B_650326'] = 'g__Pseudomonas'
  otu_alternat$genus[otu_alternat$genus == 'g__Listeria_A'] = 'g__Listeria'
    
  otu_alternat$genus[is.na(otu_alternat$genus)] = 'Others'
  otu_alternat$genus[otu_alternat$genus == 'g__'] = 'Others'
  otu_alternat$genus[otu_alternat$genus == 'g__uncultured'] = 'Others'
    
  otu_table_genus_alternat = otu_alternat %>% group_by(genus) %>% 
      summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% as.data.frame()
  otu_table_genus_alternat <- column_to_rownames(otu_table_genus_alternat, "genus")
    
  otu_table_genus_alternat = sweep(otu_table_genus_alternat,2,colSums(otu_table_genus_alternat),'/')
  otu_table_genus_alternat = subset(otu_table_genus_alternat, rowSums(otu_table_genus_alternat) > 0)

  return(otu_table_genus_alternat  %>% dplyr::select(cols_to_keep) %>% filter(rowSums(across(where(is.numeric)))>0))
}

filter_tables = function(otu_table, tax_table, alternat_thr){
  tax_table$Taxonomy.hyCAT = tax_table$Taxonomy.centroid
  tax_table$Confidence.hyCAT = tax_table$Confidence.centroid
  
  mask_alternat = (tax_table$Confidence.most.confident > (tax_table$Confidence.centroid * alternat_thr)) & (tax_table$Taxonomy.most.confident != 'assigned')
  tax_table$Taxonomy.hyCAT[mask_alternat] = tax_table$Taxonomy.most.confident[mask_alternat]
  tax_table$Confidence.hyCAT[mask_alternat] = tax_table$Confidence.most.confident[mask_alternat]
  
  otu_table[otu_table < 2] = 0 # removing singletons
  otu_table = otu_table[rowSums(otu_table) > 1,]
  tax_table = tax_table[rownames(tax_table) %in% rownames(otu_table),]
  
  samples_to_remove = colnames(otu_table)[colSums(otu_table) < 400] # removing samples with less than 400 filtered reads
  otu_table = otu_table[,!(colnames(otu_table) %in% samples_to_remove)]
  tax_table = tax_table[rownames(tax_table) %in% rownames(otu_table),]

  unassigned_otus_centroid = rownames(tax_table)[tax_table$Taxonomy.centroid == 'Unassigned']
  unassigned_otus_mostconf = rownames(tax_table)[tax_table$Taxonomy.most.confident == 'Unassigned']
  unassigned_otus_alternat = rownames(tax_table)[tax_table$Taxonomy.hyCAT == 'Unassigned']
  
  bacteria_otus_centroid = rownames(tax_table)[startsWith(tax_table$Taxonomy.centroid, 'd__Bacteria')]
  bacteria_otus_mostconf = rownames(tax_table)[startsWith(tax_table$Taxonomy.most.confident, 'd__Bacteria')]
  bacteria_otus_alternat = rownames(tax_table)[startsWith(tax_table$Taxonomy.hyCAT, 'd__Bacteria')]
  
  colsums_before = colSums(otu_table) 
  colsums_presence = colSums(otu_table > 0)
  
  prop_unassigned_centroid = colSums(otu_table[rownames(otu_table) %in% unassigned_otus_centroid,] > 0) / colsums_presence
  rela_unassigned_centroid = colSums(otu_table[rownames(otu_table) %in% unassigned_otus_centroid,]) / colsums_before
  prop_bacteria_centroid = colSums(otu_table[rownames(otu_table) %in% bacteria_otus_centroid,] > 0) / colsums_presence
  rela_bacteria_centroid = colSums(otu_table[rownames(otu_table) %in% bacteria_otus_centroid,]) / colsums_before
  
  prop_unassigned_mostconf = colSums(otu_table[rownames(otu_table) %in% unassigned_otus_mostconf,] > 0) / colsums_presence
  rela_unassigned_mostconf = colSums(otu_table[rownames(otu_table) %in% unassigned_otus_mostconf,]) / colsums_before
  prop_bacteria_mostconf = colSums(otu_table[rownames(otu_table) %in% bacteria_otus_mostconf,] > 0) / colsums_presence
  rela_bacteria_mostconf = colSums(otu_table[rownames(otu_table) %in% bacteria_otus_mostconf,]) / colsums_before
  
  prop_unassigned_alternat = colSums(otu_table[rownames(otu_table) %in% unassigned_otus_alternat,] > 0) / colsums_presence
  rela_unassigned_alternat = colSums(otu_table[rownames(otu_table) %in% unassigned_otus_alternat,]) / colsums_before
  prop_bacteria_alternat = colSums(otu_table[rownames(otu_table) %in% bacteria_otus_alternat,] > 0) / colsums_presence
  rela_bacteria_alternat = colSums(otu_table[rownames(otu_table) %in% bacteria_otus_alternat,]) / colsums_before

  summary_assigned = data.frame(Sample = rep(colnames(otu_table), 6),
                                  Metric = c(rep('Proportion of assigned', length(prop_unassigned_centroid)),
                                             rep('Proportion of assigned', length(prop_unassigned_mostconf)),
                                             rep('Proportion of assigned', length(prop_unassigned_alternat)),
                                             rep('Relative abundance of assigned', length(rela_unassigned_centroid)),
                                             rep('Relative abundance of assigned', length(rela_unassigned_mostconf)),
                                             rep('Relative abundance of assigned', length(rela_unassigned_alternat))),
                                  Method = c(rep('Centroid', length(prop_unassigned_centroid)),
                                             rep('mcCAT', length(prop_unassigned_mostconf)),
                                             rep('hyCAT', length(prop_unassigned_alternat)),
                                             rep('Centroid', length(rela_unassigned_centroid)),
                                             rep('mcCAT', length(rela_unassigned_mostconf)),
                                             rep('hyCAT', length(rela_unassigned_alternat))),
                                  Value = c(1-prop_unassigned_centroid, 1-prop_unassigned_mostconf, 1-prop_unassigned_alternat,
                                            1-rela_unassigned_centroid, 1-rela_unassigned_mostconf, 1-rela_unassigned_alternat))
  
  summary_bacteria = data.frame(Sample = rep(colnames(otu_table), 6),
                                Metric = c(rep('Proportion of Bacteria', length(prop_bacteria_centroid)),
                                             rep('Proportion of Bacteria', length(prop_bacteria_mostconf)),
                                             rep('Proportion of Bacteria', length(prop_bacteria_alternat)),
                                             rep('Relative abundance of Bacteria', length(rela_bacteria_centroid)),
                                             rep('Relative abundance of Bacteria', length(rela_bacteria_mostconf)),
                                             rep('Relative abundance of Bacteria', length(rela_bacteria_alternat))),
                                  Method = c(rep('Centroid', length(prop_bacteria_centroid)),
                                             rep('mcCAT', length(prop_bacteria_mostconf)),
                                             rep('hyCAT', length(prop_bacteria_alternat)),
                                             rep('Centroid', length(rela_bacteria_centroid)),
                                             rep('mcCAT', length(rela_bacteria_mostconf)),
                                             rep('hyCAT', length(rela_bacteria_alternat))),
                                  Value = c(prop_bacteria_centroid, prop_bacteria_mostconf, prop_bacteria_alternat,
                                            rela_bacteria_centroid, rela_bacteria_mostconf, rela_bacteria_alternat))
  
  return(list(tax_centroid=tax_table, tax_mostconf=tax_table, tax_alternat=tax_table, 
              otu_centroid=otu_table, otu_mostconf=otu_table, otu_alternat=otu_table, 
              assigned=summary_assigned, bacteria = summary_bacteria))
  }

load_study = function(study_name, thr_slv, thr_grg, thr_gdb){
  otu_table = read.csv(paste0(c('all_res/otu_tab_', study_name, '.tsv'), collapse = ''), sep='\t')
  grg_table = read.csv(paste0(c('all_res/tax_greengenes_', study_name, '.tsv'), collapse = ''))
  slv_table = read.csv(paste0(c('all_res/tax_silva_', study_name, '.tsv'), collapse = ''))
  gdb_table = read.csv(paste0(c('all_res/tax_gtdb_', study_name, '.tsv'), collapse = ''))
  
  rownames(otu_table) = otu_table$X.OTU.ID
  otu_table$X.OTU.ID = NULL
  
  rownames(grg_table) = grg_table$Cluster
  grg_table$Cluster = NULL
  
  rownames(slv_table) = slv_table$Cluster
  slv_table$Cluster = NULL
  
  rownames(gdb_table) = gdb_table$Cluster
  gdb_table$Cluster = NULL
  
  filtered_slv = filter_tables(otu_table, slv_table, thr_slv)
  filtered_grg = filter_tables(otu_table, grg_table, thr_grg)
  filtered_gdb = filter_tables(otu_table, gdb_table, thr_gdb)
  
  return(list(silva=filtered_slv, greengenes=filtered_grg, gtdb=filtered_gdb))
}

get_proportion_assigned_to_level <- function(dataset){
  levels = c('domain','phylum','class','order','family','genus','species')
  levels_assign = list()
  all_stats = data.frame()
  
  for (db in c('silva', 'greengenes', 'gtdb')){
    print(db)
    tax_mostconf = dataset[[db]]$tax_mostconf %>% 
      separate(Taxonomy.most.confident, into = c('domain','phylum','class','order','family','genus','species') ,';') %>% dplyr::select(-species)
    tax_centroid = dataset[[db]]$tax_centroid %>%
      separate(Taxonomy.centroid, into = c('domain','phylum','class','order','family','genus','species') ,';') %>% dplyr::select(-species)
    tax_alternat = dataset[[db]]$tax_alternat %>%
      separate(Taxonomy.hyCAT, into = c('domain','phylum','class','order','family','genus','species') ,';') %>% dplyr::select(-species)
    
    # put empty as NA
    tax_mostconf$domain[tax_mostconf$domain == 'd__'] = NA
    tax_centroid$domain[tax_centroid$domain == 'd__'] = NA
    tax_alternat$domain[tax_alternat$domain == 'd__'] = NA
    
    tax_mostconf$phylum[tax_mostconf$phylum == 'p__'] = NA
    tax_centroid$phylum[tax_centroid$phylum == 'p__'] = NA
    tax_alternat$phylum[tax_alternat$phylum == 'p__'] = NA
    
    tax_mostconf$class[tax_mostconf$class == 'c__'] = NA
    tax_centroid$class[tax_centroid$class == 'c__'] = NA
    tax_alternat$class[tax_alternat$class == 'c__'] = NA
    
    tax_mostconf$order[tax_mostconf$order == 'o__'] = NA
    tax_centroid$order[tax_centroid$order == 'o__'] = NA
    tax_alternat$order[tax_alternat$order == 'o__'] = NA
    
    tax_mostconf$family[tax_mostconf$family == 'f__'] = NA
    tax_centroid$family[tax_centroid$family == 'f__'] = NA
    tax_alternat$family[tax_alternat$family == 'f__'] = NA
    
    tax_mostconf$genus[tax_mostconf$genus == 'g__'] = NA
    tax_centroid$genus[tax_centroid$genus == 'g__'] = NA
    tax_alternat$genus[tax_alternat$genus == 'g__'] = NA
    
    tax_mostconf$genus[tax_mostconf$genus == 'g__uncultured'] = NA
    tax_centroid$genus[tax_centroid$genus == 'g__uncultured'] = NA
    tax_alternat$genus[tax_alternat$genus == 'g__uncultured'] = NA
    
    tab_centroid = as.data.frame(as.matrix(dataset[[db]]$otu_centroid))
    tab_mostconf = as.matrix(dataset[[db]]$otu_mostconf)
    tab_alternat = as.matrix(dataset[[db]]$otu_alternat)
    
    samples = intersect(colnames(tab_centroid), colnames(tab_mostconf))
    for (sample in samples){
      sample_centroid = tab_centroid[,sample]
      sample_centroid = sample_centroid / sum(sample_centroid)
      names(sample_centroid) = rownames(tab_centroid)
      
      sample_mostconf = tab_mostconf[,sample]
      sample_mostconf = sample_mostconf / sum(sample_mostconf)
      names(sample_mostconf) = rownames(tab_mostconf)

      sample_alternat = tab_alternat[,sample]
      sample_alternat = sample_alternat / sum(sample_alternat)
      names(sample_alternat) = rownames(tab_alternat)
      
      n_otu_centroid = sum(sample_centroid > 0)
      n_otu_mostconf = sum(sample_mostconf > 0)
      n_otu_alternat = sum(sample_alternat > 0)
      
      for (level in levels[1:6]){
        otus_centroid = rownames(tax_centroid)[!is.na(tax_centroid[,level])]
        otus_mostconf = rownames(tax_mostconf)[!is.na(tax_mostconf[,level])]
        otus_alternat = rownames(tax_alternat)[!is.na(tax_alternat[,level])]
        
        # Number of OTUs
        prop_otu_centroid = sum(sample_centroid[names(sample_centroid) %in% otus_centroid] > 0) / n_otu_centroid
        prop_otu_mostconf = sum(sample_mostconf[names(sample_mostconf) %in% otus_mostconf] > 0) / n_otu_mostconf
        prop_otu_alternat = sum(sample_alternat[names(sample_alternat) %in% otus_alternat] > 0) / n_otu_alternat
        
        # Relative abundance of the OTUs
        rela_otu_centroid = sum(sample_centroid[names(sample_centroid) %in% otus_centroid])
        rela_otu_mostconf = sum(sample_mostconf[names(sample_mostconf) %in% otus_mostconf])
        rela_otu_alternat = sum(sample_alternat[names(sample_alternat) %in% otus_alternat])
        
        if (db == 'silva'){db = 'Silva'}
        if (db == 'greengenes'){db = 'Greengenes 2'}
        if (db == 'gtdb'){db = 'GTDB'}
        
        all_stats = rbind(all_stats, data.frame(Database = db, Level = level, Sample = sample,
                                                Method = 'Centroid', Metric = 'Proportion',
                                                Value = prop_otu_centroid))
        all_stats = rbind(all_stats, data.frame(Database = db, Level = level, Sample = sample,
                                                Method = 'mcCAT', Metric = 'Proportion',
                                                Value = prop_otu_mostconf))
        all_stats = rbind(all_stats, data.frame(Database = db, Level = level, Sample = sample,
                                                Method = 'hyCAT', Metric = 'Proportion',
                                                Value = prop_otu_alternat))
        all_stats = rbind(all_stats, data.frame(Database = db, Level = level, Sample = sample,
                                                Method = 'mcCAT', Metric = 'Relative abundance',
                                                Value = rela_otu_mostconf))
        all_stats = rbind(all_stats, data.frame(Database = db, Level = level, Sample = sample,
                                                Method = 'Centroid', Metric = 'Relative abundance',
                                                Value = rela_otu_centroid))
        all_stats = rbind(all_stats, data.frame(Database = db, Level = level, Sample = sample, 
                                                Method = 'hyCAT', Metric = 'Relative abundance',
                                                Value = rela_otu_alternat))
      }
    }
  }
  
  all_stats = all_stats[!(all_stats$Sample %in% study_a_mocks),]
  all_stats = all_stats[!(all_stats$Sample %in% study_b_mocks),]
  all_stats = all_stats[!(all_stats$Sample %in% study_e_mocks),]
  all_stats = all_stats[!(all_stats$Sample %in% study_f_mocks),]
  all_stats = all_stats[!(all_stats$Sample %in% study_g_mocks),]
  
  return(all_stats)
}

assigned_otus = function(parsed_a, parsed_b, parsed_c, parsed_d, parsed_e, parsed_f, parsed_g){
  all_assigned_silva = rbind(parsed_a$silva$assigned, parsed_b$silva$assigned, 
                               parsed_c$silva$assigned, parsed_d$silva$assigned,
                               parsed_e$silva$assigned, parsed_f$silva$assigned,
                               parsed_g$silva$assigned)
  all_assigned_silva$Database = 'Silva'
  
  all_assigned_greengenes = rbind(parsed_a$greengenes$assigned, parsed_b$greengenes$assigned, 
                                    parsed_c$greengenes$assigned, parsed_d$greengenes$assigned,
                                    parsed_e$greengenes$assigned, parsed_f$greengenes$assigned,
                                    parsed_g$greengenes$assigned)
  all_assigned_greengenes$Database = 'Greengenes 2'

  all_assigned_gtdb = rbind(parsed_a$gtdb$assigned, parsed_b$gtdb$assigned, 
                              parsed_c$gtdb$assigned, parsed_d$gtdb$assigned,
                              parsed_e$gtdb$assigned, parsed_f$gtdb$assigned,
                              parsed_g$gtdb$assigned)
  all_assigned_gtdb$Database = 'GTDB'
  
  all_assigned = rbind(all_assigned_silva, all_assigned_gtdb, all_assigned_greengenes)
  
  all_assigned = all_assigned[!(all_assigned$Sample %in% study_a_mocks),]
  all_assigned = all_assigned[!(all_assigned$Sample %in% study_b_mocks),]
  all_assigned = all_assigned[!(all_assigned$Sample %in% study_e_mocks),]
  all_assigned = all_assigned[!(all_assigned$Sample %in% study_f_mocks),]
  all_assigned = all_assigned[!(all_assigned$Sample %in% study_g_mocks),]
  
  return(all_assigned)   
}

bacterial_otus = function(parsed_a, parsed_b, parsed_c, parsed_d, parsed_e, parsed_f, parsed_g){
  all_bacteria_silva = rbind(parsed_a$silva$bacteria, parsed_b$silva$bacteria, 
                             parsed_c$silva$bacteria, parsed_d$silva$bacteria,
                             parsed_e$silva$bacteria, parsed_f$silva$bacteria,
                             parsed_g$silva$bacteria)
  all_bacteria_silva$Database = 'Silva'
  
  
  all_bacteria_greengenes = rbind(parsed_a$greengenes$bacteria, parsed_b$greengenes$bacteria, 
                                  parsed_c$greengenes$bacteria, parsed_d$greengenes$bacteria,
                                  parsed_e$greengenes$bacteria, parsed_f$greengenes$bacteria,
                                  parsed_g$greengenes$bacteria)
  all_bacteria_greengenes$Database = 'Greengenes 2'
  
  
  all_bacteria_gtdb = rbind(parsed_a$gtdb$bacteria, parsed_b$gtdb$bacteria, 
                            parsed_c$gtdb$bacteria, parsed_d$gtdb$bacteria,
                            parsed_e$gtdb$bacteria, parsed_f$gtdb$bacteria,
                            parsed_g$gtdb$bacteria)
  all_bacteria_gtdb$Database = 'GTDB'
  
  all_bacteria = rbind(all_bacteria_silva, all_bacteria_greengenes, all_bacteria_gtdb)
  
  all_bacteria = all_bacteria[!(all_bacteria$Sample %in% study_a_mocks),]
  all_bacteria = all_bacteria[!(all_bacteria$Sample %in% study_b_mocks),]
  all_bacteria = all_bacteria[!(all_bacteria$Sample %in% study_e_mocks),]
  all_bacteria = all_bacteria[!(all_bacteria$Sample %in% study_f_mocks),]
  all_bacteria = all_bacteria[!(all_bacteria$Sample %in% study_g_mocks),]
  
  return(all_bacteria)  
}

summarise_genus <- function(parsed_data, cols_to_keep){
  out_tables = list()
  for (db in c('silva', 'gtdb', 'greengenes')){ 
    otu_centroid = parsed_data[[db]]$otu_centroid
    otu_mostconf = parsed_data[[db]]$otu_mostconf 
    otu_alternat = parsed_data[[db]]$otu_alternat
    
    sep_db = ';'
    if (db == 'greengenes'){sep_db = '; '}
    tax_centroid = parsed_data[[db]]$tax_centroid %>% separate(Taxonomy.centroid,  
                                                               into = c('domain','phylum','class','order','family','genus','species') , sep_db)
    tax_mostconf = parsed_data[[db]]$tax_mostconf %>% separate(Taxonomy.most.confident,  
                                                               into = c('domain','phylum','class','order','family','genus','species') , sep_db)
    tax_alternat = parsed_data[[db]]$tax_alternat %>% separate(Taxonomy.hyCAT,  
                                                               into = c('domain','phylum','class','order','family','genus','species') , sep_db)
    
    #tax_table_centroid = tax_table[grepl('d__Bacteria', tax_table$Taxonomy.centroid),]
    #tax_table_centroid = tax_table_centroid[!(grepl('Mitochondria', tax_table_centroid$Taxonomy.centroid)),]
    #tax_table_centroid = tax_table_centroid[!(grepl('Chloroplast', tax_table_centroid$Taxonomy.centroid)),]
    #filtered_asvs_centroid = rownames(tax_table_centroid)
    
    #tax_table_mostconf = tax_table[grepl('d__Bacteria', tax_table$Taxonomy.most.confident),]
    #tax_table_mostconf = tax_table_mostconf[!(grepl('Mitochondria', tax_table_mostconf$Taxonomy.most.confident)),]
    #tax_table_mostconf = tax_table_mostconf[!(grepl('Chloroplast', tax_table_mostconf$Taxonomy.most.confident)),]
    #filtered_asvs_mostconf = rownames(tax_table_mostconf)
    
    #tax_table_alternat = tax_table[grepl('d__Bacteria', tax_table$Taxonomy.hyCAT),]
    #tax_table_alternat = tax_table_alternat[!(grepl('Mitochondria', tax_table_alternat$Taxonomy.hyCAT)),]
    #tax_table_alternat = tax_table_alternat[!(grepl('Chloroplast', tax_table_alternat$Taxonomy.hyCAT)),]
    #filtered_asvs_alternat = rownames(tax_table_alternat)
    
    otu_centroid$genus = map_chr(rownames(otu_centroid), function(x) tax_centroid$genus[rownames(tax_centroid) == x])
    otu_mostconf$genus = map_chr(rownames(otu_mostconf), function(x) tax_mostconf$genus[rownames(tax_mostconf) == x])
    otu_alternat$genus = map_chr(rownames(otu_alternat), function(x) tax_alternat$genus[rownames(tax_alternat) == x])
    
    otu_centroid$genus[otu_centroid$genus == 'g__Limosilactobacillus'] = 'g__Lactobacillus' # Limosilactobacillus is part of Lactobacillus as of 2020
    otu_centroid$genus[otu_centroid$genus == 'g__Cereibacter'] = 'g__Rhodobacter' # Rhodobacter is also Cereibacter https://gtdb.ecogenomic.org/searches?s=al&q=g__Rhodobacter
    otu_centroid$genus[otu_centroid$genus == 'g__Clostridium_sensu_stricto_1'] = 'g__Clostridium'
    otu_centroid$genus[otu_centroid$genus == 'g__Escherichia-Shigella'] = 'g__Escherichia'
    otu_centroid$genus[otu_centroid$genus == 'g__Bacillus_A'] = 'g__Bacillus'
    otu_centroid$genus[otu_centroid$genus == 'g__Cereibacter_A'] = 'g__Rhodobacter'
    otu_centroid$genus[otu_centroid$genus == 'g__Bifidobacterium_388775'] = 'g__Bifidobacterium'
    otu_centroid$genus[otu_centroid$genus == 'g__Clostridium_T'] = 'g__Clostridium'
    otu_centroid$genus[otu_centroid$genus == 'g__Deinococcus_B'] = 'g__Deinococcus'
    otu_centroid$genus[otu_centroid$genus == 'g__Cereibacter_A_494292'] = 'g__Rhodobacter'
    otu_centroid$genus[otu_centroid$genus == 'g__Enterococcus_H_360604'] = 'g__Enterococcus'
    otu_centroid$genus[otu_centroid$genus == 'g__Bacillus_P_294101'] = 'g__Bacillus'
    otu_centroid$genus[otu_centroid$genus == 'g__Pseudomonas_B_650326'] = 'g__Pseudomonas'
    otu_centroid$genus[otu_centroid$genus == 'g__Listeria_A'] = 'g__Listeria'
    
    otu_mostconf$genus[otu_mostconf$genus == 'g__Limosilactobacillus'] = 'g__Lactobacillus' # Limosilactobacillus is part of Lactobacillus as of 2020
    otu_mostconf$genus[otu_mostconf$genus == 'g__Cereibacter'] = 'g__Rhodobacter' # Rhodobacter is also Cereibacter https://gtdb.ecogenomic.org/searches?s=al&q=g__Rhodobacter
    otu_mostconf$genus[otu_mostconf$genus == 'g__Clostridium_sensu_stricto_1'] = 'g__Clostridium'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Escherichia-Shigella'] = 'g__Escherichia'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Bacillus_A'] = 'g__Bacillus'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Cereibacter_A'] = 'g__Rhodobacter'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Bifidobacterium_388775'] = 'g__Bifidobacterium'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Clostridium_T'] = 'g__Clostridium'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Deinococcus_B'] = 'g__Deinococcus'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Cereibacter_A_494292'] = 'g__Rhodobacter'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Enterococcus_H_360604'] = 'g__Enterococcus'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Bacillus_P_294101'] = 'g__Bacillus'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Pseudomonas_B_650326'] = 'g__Pseudomonas'
    otu_mostconf$genus[otu_mostconf$genus == 'g__Listeria_A'] = 'g__Listeria'
    
    otu_alternat$genus[otu_alternat$genus == 'g__Limosilactobacillus'] = 'g__Lactobacillus' # Limosilactobacillus is part of Lactobacillus as of 2020
    otu_alternat$genus[otu_alternat$genus == 'g__Cereibacter'] = 'g__Rhodobacter' # Rhodobacter is also Cereibacter https://gtdb.ecogenomic.org/searches?s=al&q=g__Rhodobacter
    otu_alternat$genus[otu_alternat$genus == 'g__Clostridium_sensu_stricto_1'] = 'g__Clostridium'
    otu_alternat$genus[otu_alternat$genus == 'g__Escherichia-Shigella'] = 'g__Escherichia'
    otu_alternat$genus[otu_alternat$genus == 'g__Bacillus_A'] = 'g__Bacillus'
    otu_alternat$genus[otu_alternat$genus == 'g__Cereibacter_A'] = 'g__Rhodobacter'
    otu_alternat$genus[otu_alternat$genus == 'g__Bifidobacterium_388775'] = 'g__Bifidobacterium'
    otu_alternat$genus[otu_alternat$genus == 'g__Clostridium_T'] = 'g__Clostridium'
    otu_alternat$genus[otu_alternat$genus == 'g__Deinococcus_B'] = 'g__Deinococcus'
    otu_alternat$genus[otu_alternat$genus == 'g__Cereibacter_A_494292'] = 'g__Rhodobacter'
    otu_alternat$genus[otu_alternat$genus == 'g__Enterococcus_H_360604'] = 'g__Enterococcus'
    otu_alternat$genus[otu_alternat$genus == 'g__Bacillus_P_294101'] = 'g__Bacillus'
    otu_alternat$genus[otu_alternat$genus == 'g__Pseudomonas_B_650326'] = 'g__Pseudomonas'
    otu_alternat$genus[otu_alternat$genus == 'g__Listeria_A'] = 'g__Listeria'
  
    otu_centroid$genus[is.na(otu_centroid$genus)] = 'Others'
    otu_mostconf$genus[is.na(otu_mostconf$genus)] = 'Others'
    otu_alternat$genus[is.na(otu_alternat$genus)] = 'Others'
    
    otu_centroid$genus[otu_centroid$genus == 'g__'] = 'Others'
    otu_mostconf$genus[otu_mostconf$genus == 'g__'] = 'Others'
    otu_alternat$genus[otu_alternat$genus == 'g__'] = 'Others'
    
    otu_centroid$genus[otu_centroid$genus == 'g__uncultured'] = 'Others'
    otu_mostconf$genus[otu_mostconf$genus == 'g__uncultured'] = 'Others'
    otu_alternat$genus[otu_alternat$genus == 'g__uncultured'] = 'Others'
    
    otu_table_genus_centroid = otu_centroid %>% group_by(genus) %>% 
      summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% as.data.frame()
    otu_table_genus_mostconf = otu_mostconf %>% group_by(genus) %>% 
      summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% as.data.frame()
    otu_table_genus_alternat = otu_alternat %>% group_by(genus) %>% 
      summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% as.data.frame()
    
    otu_table_genus_centroid <- column_to_rownames(otu_table_genus_centroid, "genus")
    otu_table_genus_mostconf <- column_to_rownames(otu_table_genus_mostconf, "genus")
    otu_table_genus_alternat <- column_to_rownames(otu_table_genus_alternat, "genus")
    
    #otu_table_genus_centroid = sweep(otu_table_genus_centroid,2,colSums(otu_table_genus_centroid),'/')
    #otu_table_genus_mostconf = sweep(otu_table_genus_mostconf,2,colSums(otu_table_genus_mostconf),'/')
    #otu_table_genus_alternat = sweep(otu_table_genus_alternat,2,colSums(otu_table_genus_alternat),'/')
    
    otu_table_genus_centroid = subset(otu_table_genus_centroid, rowSums(otu_table_genus_centroid) > 0)
    otu_table_genus_mostconf = subset(otu_table_genus_mostconf, rowSums(otu_table_genus_mostconf) > 0)
    otu_table_genus_alternat = subset(otu_table_genus_alternat, rowSums(otu_table_genus_alternat) > 0)
    
    out_tables[[db]] = list(centroid=otu_table_genus_centroid %>% dplyr::select(cols_to_keep) , 
                            mostconf=otu_table_genus_mostconf %>% dplyr::select(cols_to_keep) ,
                            alternat=otu_table_genus_alternat %>% dplyr::select(cols_to_keep))}
  return(out_tables)
}

compare_genus_to_real <- function(table, expected){
  count_table = merge(table, expected, by = "row.names", all = TRUE) %>% column_to_rownames(var='Row.names')
  count_table[is.na(count_table)] = 0
  filtered_count_table = count_table[rowSums(count_table) > 0,]
  counts_all = colSums(filtered_count_table)[colnames(filtered_count_table) != 'real']
  counts_expected = colSums(filtered_count_table[filtered_count_table$real > 0,])[colnames(filtered_count_table) != 'real']
  counts_spurious = colSums(filtered_count_table[filtered_count_table$real == 0,])[colnames(filtered_count_table) != 'real']
  
  rela_table = sweep(table,2,colSums(table),'/')
  merged_table <- merge(rela_table, expected, by = "row.names", all = TRUE) %>% column_to_rownames(var='Row.names')
  merged_table[is.na(merged_table)] = 0
  filtered_table = merged_table[rowSums(merged_table) > 0,]

  bc_diss = as.matrix(vegdist(t(merged_table), method = 'bray'))['real',]
  bc_diss = bc_diss[names(bc_diss) != 'real']
  
  jc_diss = as.matrix(vegdist(t(merged_table), method = 'jaccard'))['real',]
  jc_diss = jc_diss[names(jc_diss) != 'real']
  sum_dev = colSums(abs(dplyr::select(merged_table, -real) - merged_table$real))
  
  no_others = filtered_table[rownames(filtered_table) != 'Others',]
  unexpected_taxa = no_others[!(rownames(no_others) %in% rownames(expected)),]
  expected_taxa = no_others[rownames(no_others) %in% rownames(expected),]
  n_unexpected_taxa = colSums(unexpected_taxa > 0)
  n_unexpected_taxa = n_unexpected_taxa[names(n_unexpected_taxa) != 'real']
  n_taxa_found = colSums(expected_taxa > 0)
  n_taxa_found = n_taxa_found[names(n_taxa_found) != 'real']
  n_taxa_expected = nrow(expected) - 1 # -1 for the Others category
  
  return(list(sum_dev = sum_dev, bc_diss = bc_diss, n_otu = n_taxa_found, jc_diss = jc_diss, prop_exp = n_taxa_found / n_taxa_expected, 
              num_unexp = n_unexpected_taxa, ex_rich = rep(n_taxa_expected, length(sum_dev)), counts_expected = counts_expected, counts_all = counts_all, counts_spurious=counts_spurious))}

get_mock_results <- function(parsed_a, parsed_b, parsed_c, parsed_d, parsed_e, parsed_f, parsed_g, CONCOMPRA_genus){
  mock_data = data.frame()
  
  # Mock study A
  mock_genus_a = summarise_genus(parsed_a, study_a_mocks)
  mock_genus_a_real = data.frame(row.names = c("Others","g__Bacillus","g__Bifidobacterium","g__Clostridium",
                                               "g__Deinococcus","g__Enterococcus","g__Escherichia",
                                               "g__Lactobacillus","g__Rhodobacter","g__Staphylococcus","g__Streptococcus"),
                                 real = c(0,rep(0.1,10)))
  
  for (db in c('silva', 'greengenes', 'gtdb')){
    print(db)
    beta_to_real_centroid = compare_genus_to_real(mock_genus_a[[db]]$centroid, mock_genus_a_real)
    beta_to_real_mostconf = compare_genus_to_real(mock_genus_a[[db]]$mostconf, mock_genus_a_real)
    beta_to_real_alternat = compare_genus_to_real(mock_genus_a[[db]]$alternat, mock_genus_a_real)
  
    subset_CONCOMPRA = CONCOMPRA_genus[[db]] %>% dplyr::select(all_of(study_a_mocks)) 
    beta_to_real_CONCOMPRA = compare_genus_to_real(subset_CONCOMPRA, mock_genus_a_real)
    datasets = list(beta_to_real_centroid, beta_to_real_mostconf, beta_to_real_alternat, beta_to_real_CONCOMPRA)
    methods = c('Centroid', 'mcCAT', 'hyCAT', 'CONCOMPRA')
    
    for (i in 1:4){
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$sum_dev),
                                            Database = rep(db, 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('Summed deviation', 4),
                                            Value = datasets[[i]]$sum_dev))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$bc_diss),
                                            Database = rep(db, 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('Bray-Curtis dissimilarity', 4),
                                            Value = datasets[[i]]$bc_diss))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$n_otu),
                                            Database = rep(db, 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('OTU number', 4),
                                            Value = datasets[[i]]$n_otu))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$jc_diss),
                                            Database = rep(db, 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('Jaccard dissimilarity', 4),
                                            Value = datasets[[i]]$jc_diss))
    mock_data = rbind(mock_data, data.frame(Sample = study_a_mocks,
                                            Database = rep(db, 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('Expected richness', 4),
                                            Value = datasets[[i]]$ex_rich))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$prop_exp),
                                            Database = rep(db, 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('Proportion expected', 4),
                                            Value = datasets[[i]]$prop_exp))                   
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$num_unexp),
                                            Database = rep(db , 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('Number unexpected', 4),
                                            Value = datasets[[i]]$num_unexp))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_expected),
                                            Database = rep(db , 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('Counts expected', 4),
                                            Value = datasets[[i]]$counts_expected))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                            Database = rep(db , 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('Counts all', 4),
                                            Value = datasets[[i]]$counts_all))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                            Database = rep(db , 4),
                                            Study = rep('Study A', 4),
                                            Method = rep(methods[i], 4),
                                            Metric = rep('Counts spurious', 4),
                                            Value = datasets[[i]]$counts_spurious))}}                       

  ###################################################################################################
  # Mock study B
  mock_genus_b = summarise_genus(parsed_b, study_b_mocks)
  mock_genus_b_real = data.frame(row.names = c("Others","g__Pseudomonas","g__Escherichia","g__Salmonella",
                                                   "g__Lactobacillus","g__Enterococcus","g__Staphylococcus",
                                                   "g__Listeria","g__Bacillus"),
                                     real = c(0,rep(1/8,8)))

  # silva B
  for (db in c('silva', 'greengenes', 'gtdb')){
  beta_to_real_silva_centroid = compare_genus_to_real(mock_genus_b[[db]]$centroid, mock_genus_b_real)
  beta_to_real_silva_mostconf = compare_genus_to_real(mock_genus_b[[db]]$mostconf, mock_genus_b_real)
  beta_to_real_silva_alternat = compare_genus_to_real(mock_genus_b[[db]]$alternat, mock_genus_b_real)
  subset_CONCOMPRA = CONCOMPRA_genus[[db]] %>% dplyr::select(any_of(study_b_mocks))
  beta_to_real_CONCOMPRA = compare_genus_to_real(subset_CONCOMPRA, mock_genus_b_real)
  datasets = list(beta_to_real_silva_centroid, beta_to_real_silva_mostconf, beta_to_real_silva_alternat, beta_to_real_CONCOMPRA)
  methods = c('Centroid', 'mcCAT', 'hyCAT', 'CONCOMPRA')

    for (i in 1:4){
  mock_data = rbind(mock_data, data.frame(Sample = study_b_mocks,
                                          Database = rep(db, 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('Summed deviation', 1),
                                          Value = datasets[[i]]$sum_dev))
  mock_data = rbind(mock_data, data.frame(Sample = study_b_mocks,
                                          Database = rep(db, 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('Bray-Curtis dissimilarity', 1),
                                          Value = datasets[[i]]$bc_diss))
  mock_data = rbind(mock_data, data.frame(Sample = study_b_mocks,
                                          Database = rep(db, 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('OTU number', 1),
                                          Value = datasets[[i]]$n_otu))
  mock_data = rbind(mock_data, data.frame(Sample = study_b_mocks,
                                          Database = rep(db, 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('Jaccard dissimilarity', 1),
                                          Value = datasets[[i]]$jc_diss))
  mock_data = rbind(mock_data, data.frame(Sample = study_b_mocks,
                                          Database = rep(db, 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('Expected richness', 1),
                                          Value = datasets[[i]]$ex_rich))
  mock_data = rbind(mock_data, data.frame(Sample = study_b_mocks,
                                          Database = rep(db, 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('Proportion expected', 1),
                                          Value = datasets[[i]]$prop_exp))                   
  mock_data = rbind(mock_data, data.frame(Sample = study_b_mocks,
                                          Database = rep(db , 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('Number unexpected', 1),
                                          Value = datasets[[i]]$num_unexp))
  mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_expected),
                                          Database = rep(db , 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('Counts expected', 1),
                                          Value = datasets[[i]]$counts_expected))
  mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                          Database = rep(db , 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('Counts all', 1),
                                          Value = datasets[[i]]$counts_all))
  mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                          Database = rep(db , 1),
                                          Study = rep('Study B', 1),
                                          Method = rep(methods[i], 1),
                                          Metric = rep('Counts spurious', 1),
                                          Value = datasets[[i]]$counts_spurious))}}   
  
  ###################################################################################################
  # Mock study E
  mock_genus_e = summarise_genus(parsed_e, study_e_mocks)
  mock_genus_e_real = data.frame(row.names = c("Others","g__Bacillus","g__Listeria",
                                                   "g__Staphylococcus","g__Enterococcus",
                                                   "g__Lactobacillus","g__Salmonella",
                                                   "g__Escherichia", "g__Pseudomonas"),
                                     real = c(0,0.174,0.141,0.155,0.099,0.184,0.104,0.101,0.042))
  
  for (db in c('silva', 'greengenes', 'gtdb')){
    beta_to_real_silva_centroid = compare_genus_to_real(mock_genus_e[[db]]$centroid, mock_genus_e_real)
    beta_to_real_silva_mostconf = compare_genus_to_real(mock_genus_e[[db]]$mostconf, mock_genus_e_real)
    beta_to_real_silva_alternat = compare_genus_to_real(mock_genus_e[[db]]$alternat, mock_genus_e_real)
    subset_CONCOMPRA = CONCOMPRA_genus[[db]] %>% dplyr::select(all_of(study_e_mocks))
    subset_CONCOMPRA = subset_CONCOMPRA[rowSums(subset_CONCOMPRA) > 0,]
    beta_to_real_CONCOMPRA = compare_genus_to_real(subset_CONCOMPRA, mock_genus_e_real)
    datasets = list(beta_to_real_silva_centroid, beta_to_real_silva_mostconf, beta_to_real_silva_alternat, beta_to_real_CONCOMPRA)
    methods = c('Centroid', 'mcCAT', 'hyCAT', 'CONCOMPRA')
    
    for (i in 1:4){
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$sum_dev),
                                            Database = rep(db, 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('Summed deviation', 12),
                                            Value = datasets[[i]]$sum_dev))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$sum_dev),
                                            Database = rep(db, 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('Proportion expected', 12),
                                            Value = datasets[[i]]$prop_exp)) 
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$bc_diss),
                                            Database = rep(db, 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('Bray-Curtis dissimilarity', 12),
                                            Value = datasets[[i]]$bc_diss))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$n_otu),
                                            Database = rep(db, 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('OTU number', 12),
                                            Value = datasets[[i]]$n_otu))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$jc_diss),
                                            Database = rep(db, 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('Jaccard dissimilarity', 12),
                                            Value = datasets[[i]]$jc_diss))
    mock_data = rbind(mock_data, data.frame(Sample = study_e_mocks,
                                            Database = rep(db, 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('Expected richness', 12),
                                            Value = datasets[[i]]$ex_rich))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$num_unexp),
                                            Database = rep(db , 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('Number unexpected', 12),
                                            Value = datasets[[i]]$num_unexp))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_expected),
                                            Database = rep(db , 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('Counts expected', 12),
                                            Value = datasets[[i]]$counts_expected))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                            Database = rep(db , 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('Counts all', 12),
                                            Value = datasets[[i]]$counts_all))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                            Database = rep(db , 12),
                                            Study = rep('Study E', 12),
                                            Method = rep(methods[i], 12),
                                            Metric = rep('Counts spurious', 12),
                                            Value = datasets[[i]]$counts_spurious))}}   
  
  ###################################################################################################
  # Mock study F
  mock_genus_f = summarise_genus(parsed_f, study_f_mocks)
  mock_genus_f_real = data.frame(row.names = c("Others","g__Pseudomonas","g__Escherichia","g__Salmonella",
                                                   "g__Lactobacillus","g__Enterococcus","g__Staphylococcus",
                                                   "g__Listeria","g__Bacillus"),
                                     real = c(0,rep(1/8,8)))
  
  for (db in c('silva', 'greengenes', 'gtdb')){
    beta_to_real_silva_centroid = compare_genus_to_real(mock_genus_f[[db]]$centroid, mock_genus_f_real)
    beta_to_real_silva_mostconf = compare_genus_to_real(mock_genus_f[[db]]$mostconf, mock_genus_f_real)
    beta_to_real_silva_alternat = compare_genus_to_real(mock_genus_f[[db]]$alternat, mock_genus_f_real)
    subset_CONCOMPRA = CONCOMPRA_genus[[db]] %>% dplyr::select(all_of(study_f_mocks)) 
    beta_to_real_CONCOMPRA = compare_genus_to_real(subset_CONCOMPRA, mock_genus_f_real)
    datasets = list(beta_to_real_silva_centroid, beta_to_real_silva_mostconf, beta_to_real_silva_alternat, beta_to_real_CONCOMPRA)
    methods = c('Centroid', 'mcCAT', 'hyCAT', 'CONCOMPRA')
    
    for (i in 1:4){
    mock_data = rbind(mock_data, data.frame(Sample = study_f_mocks,
                                            Database = rep(db, 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('Summed deviation', 1),
                                            Value = datasets[[i]]$sum_dev))
    mock_data = rbind(mock_data, data.frame(Sample = study_f_mocks,
                                            Database = rep(db, 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('Proportion expected', 1),
                                            Value = datasets[[i]]$prop_exp)) 
    mock_data = rbind(mock_data, data.frame(Sample = study_f_mocks,
                                            Database = rep(db, 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('Bray-Curtis dissimilarity', 1),
                                            Value = datasets[[i]]$bc_diss))
    mock_data = rbind(mock_data, data.frame(Sample = study_f_mocks,
                                            Database = rep(db, 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('OTU number', 1),
                                            Value = datasets[[i]]$n_otu))
    mock_data = rbind(mock_data, data.frame(Sample = study_f_mocks,
                                            Database = rep(db, 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('Jaccard dissimilarity', 1),
                                            Value = datasets[[i]]$jc_diss))
    mock_data = rbind(mock_data, data.frame(Sample = study_f_mocks,
                                            Database = rep(db, 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('Expected richness', 1),
                                            Value = datasets[[i]]$ex_rich))
    mock_data = rbind(mock_data, data.frame(Sample = study_f_mocks,
                                            Database = rep(db , 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('Number unexpected', 1),
                                            Value = datasets[[i]]$num_unexp))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_expected),
                                            Database = rep(db , 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('Counts expected', 1),
                                            Value = datasets[[i]]$counts_expected))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                            Database = rep(db , 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('Counts all', 1),
                                            Value = datasets[[i]]$counts_all))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                            Database = rep(db , 1),
                                            Study = rep('Study F', 1),
                                            Method = rep(methods[i], 1),
                                            Metric = rep('Counts spurious', 1),
                                            Value = datasets[[i]]$counts_spurious))}}
  
  ###################################################################################################
  # Mock study G
  mock_genus_g = summarise_genus(parsed_g, study_g_mocks)
  mock_genus_g_real = data.frame(row.names = c("Others","g__Schaalia","g__Bifidobacterium","g__Cutibacterium",
                                                        "g__Phocaeicola","g__Porphyromonas","g__Deinococcus",
                                                        "g__Cereibacter","g__Neisseria","g__Acinetobacter",
                                                        "g__Pseudomonas","g__Escherichia","g__Helicobacter",
                                                        "g__Bacillus","g__Clostridium","g__Enterococcus",
                                                        "g__Lactobacillus","g__Staphylococcus","g__Streptococcus"),
                                     real = c(0,rep(1/20,16),2/20,2/20))
  
  for (db in c('silva', 'greengenes', 'gtdb')){
    beta_to_real_silva_centroid = compare_genus_to_real(mock_genus_g[[db]]$centroid, mock_genus_g_real)
    beta_to_real_silva_mostconf = compare_genus_to_real(mock_genus_g[[db]]$mostconf, mock_genus_g_real)
    beta_to_real_silva_alternat = compare_genus_to_real(mock_genus_g[[db]]$alternat, mock_genus_g_real)
    subset_CONCOMPRA = CONCOMPRA_genus[[db]] %>% dplyr::select(all_of(study_g_mocks))
    subset_CONCOMPRA = subset_CONCOMPRA[rowSums(subset_CONCOMPRA) > 0,]
    beta_to_real_CONCOMPRA = compare_genus_to_real(subset_CONCOMPRA, mock_genus_g_real)
    datasets = list(beta_to_real_silva_centroid, beta_to_real_silva_mostconf, beta_to_real_silva_alternat, beta_to_real_CONCOMPRA)
    methods = c('Centroid', 'mcCAT', 'hyCAT', 'CONCOMPRA')
    
    for (i in 1:4){
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$sum_dev),
                                            Database = rep(db, 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('Summed deviation', 2),
                                            Value = datasets[[i]]$sum_dev))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$prop_exp),
                                            Database = rep(db, 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('Proportion expected', 2),
                                            Value = datasets[[i]]$prop_exp))     
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$bc_diss),
                                            Database = rep(db, 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('Bray-Curtis dissimilarity', 2),
                                            Value = datasets[[i]]$bc_diss))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$n_otu),
                                            Database = rep(db, 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('OTU number', 2),
                                            Value = datasets[[i]]$n_otu))
    mock_data = rbind(mock_data, data.frame(Sample = study_g_mocks,
                                            Database = rep(db, 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('Jaccard dissimilarity', 2),
                                            Value = datasets[[i]]$jc_diss))
    mock_data = rbind(mock_data, data.frame(Sample = study_g_mocks,
                                            Database = rep(db, 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('Expected richness', 2),
                                            Value = datasets[[i]]$ex_rich))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$num_unexp),
                                            Database = rep(db , 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('Number unexpected', 2),
                                            Value = datasets[[i]]$num_unexp))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_expected),
                                            Database = rep(db , 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('Counts expected', 2),
                                            Value = datasets[[i]]$counts_expected))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                            Database = rep(db, 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('Counts all', 2),
                                            Value = datasets[[i]]$counts_all))
    mock_data = rbind(mock_data, data.frame(Sample = names(datasets[[i]]$counts_all),
                                            Database = rep(db , 2),
                                            Study = rep('Study G', 2),
                                            Method = rep(methods[i], 2),
                                            Metric = rep('Counts spurious', 2),
                                            Value = datasets[[i]]$counts_spurious))}} 
  return(mock_data %>% filter(Sample != 'x'))}

load_CONCOMPRA = function(study_name){
  otu_table = read.csv('all_res/otu_tab_CONCOMPRA.tsv', sep=',')
  grg_table = read.csv('all_res/tax_greengenes_CONCOMPRA.tsv', sep='\t')
  slv_table = read.csv('all_res/tax_silva_CONCOMPRA.tsv', sep='\t')
  gdb_table = read.csv('all_res/tax_gtdb_CONCOMPRA.tsv', sep='\t')
  
  grg_table$Cluster = map_chr(grg_table$Cluster, function(x) strsplit(x, ';')[[1]][1])
  slv_table$Cluster = map_chr(slv_table$Cluster, function(x) strsplit(x, ';')[[1]][1])
  gdb_table$Cluster = map_chr(gdb_table$Cluster, function(x) strsplit(x, ';')[[1]][1])
  
  rownames(otu_table) = otu_table$X.OTU.ID
  otu_table$X.OTU.ID = NULL
  
  rownames(grg_table) = grg_table$Cluster
  grg_table$Cluster = NULL
  
  rownames(slv_table) = slv_table$Cluster
  slv_table$Cluster = NULL
  
  rownames(gdb_table) = gdb_table$Cluster
  gdb_table$Cluster = NULL
  
  grg_table = grg_table[grepl('d__Bacteria', grg_table$Taxonomy),]
  grg_table = grg_table[!(grepl('Mitochondria', grg_table$Taxonomy)),]
  grg_table = grg_table[!(grepl('Chloroplast', grg_table$Taxonomy)),]
  filtered_asvs_grg = rownames(grg_table)
  
  slv_table = slv_table[grepl('d__Bacteria', slv_table$Taxonomy),]
  slv_table = slv_table[!(grepl('Mitochondria', slv_table$Taxonomy)),]
  slv_table = slv_table[!(grepl('Chloroplast', slv_table$Taxonomy)),]
  filtered_asvs_slv = rownames(slv_table)
  
  gdb_table = gdb_table[grepl('d__Bacteria', gdb_table$Taxonomy),]
  gdb_table = gdb_table[!(grepl('Mitochondria', gdb_table$Taxonomy)),]
  gdb_table = gdb_table[!(grepl('Chloroplast', gdb_table$Taxonomy)),]
  filtered_asvs_gdb = rownames(gdb_table)
  
  otu_table_slv = otu_table[rownames(otu_table) %in% filtered_asvs_slv,]
  otu_table_slv_filtered = otu_table_slv[rowSums(otu_table_slv) > 1,]
  tax_table_slv_filtered = slv_table[rownames(slv_table) %in% rownames(otu_table_slv_filtered),]
  
  otu_table_grg = otu_table[rownames(otu_table) %in% filtered_asvs_grg,]
  otu_table_grg_filtered = otu_table_grg[rowSums(otu_table_grg) > 1,]
  tax_table_grg_filtered = grg_table[rownames(grg_table) %in% rownames(otu_table_grg_filtered),]
  
  otu_table_gdb = otu_table[rownames(otu_table) %in% filtered_asvs_gdb,]
  otu_table_gdb_filtered = otu_table_gdb[rowSums(otu_table_gdb) > 1,]
  tax_table_gdb_filtered = gdb_table[rownames(gdb_table) %in% rownames(otu_table_gdb_filtered),]
  
  filtered_slv = list(otu_table=otu_table_slv_filtered,tax_table=tax_table_slv_filtered)
  filtered_grg = list(otu_table=otu_table_grg_filtered,tax_table=tax_table_grg_filtered)
  filtered_gdb = list(otu_table=otu_table_gdb_filtered,tax_table=tax_table_gdb_filtered)
  
  return(list(silva=filtered_slv, greengenes=filtered_grg, gtdb=filtered_gdb))
}

summarise_CONCOMPRA <- function(parsed_data){
  out_tables = list()
  for (db in c('silva', 'gtdb', 'greengenes')){ 
    otu_table = parsed_data[[db]]$otu_table
    
    sep_db = ';'
    if (db == 'greengenes'){sep_db = '; '}
    tax_table = parsed_data[[db]]$tax_table %>% separate(Taxonomy,  
                                                         into = c('domain','phylum','class','order','family','genus','species') , sep_db)
    
    otu_table$genus = map_chr(rownames(otu_table), function(x) tax_table$genus[rownames(tax_table) == x])
    
    otu_table$genus[otu_table$genus == 'g__Limosilactobacillus'] = 'g__Lactobacillus' # Limosilactobacillus is part of Lactobacillus as of 2020
    otu_table$genus[otu_table$genus == 'g__Cereibacter'] = 'g__Rhodobacter' # Rhodobacter is also Cereibacter https://gtdb.ecogenomic.org/searches?s=al&q=g__Rhodobacter
    otu_table$genus[otu_table$genus == 'g__Clostridium_sensu_stricto_1'] = 'g__Clostridium'
    otu_table$genus[otu_table$genus == 'g__Escherichia-Shigella'] = 'g__Escherichia'
    otu_table$genus[otu_table$genus == 'g__Bacillus_A'] = 'g__Bacillus'
    otu_table$genus[otu_table$genus == 'g__Cereibacter_A'] = 'g__Rhodobacter'
    otu_table$genus[otu_table$genus == 'g__Bifidobacterium_388775'] = 'g__Bifidobacterium'
    otu_table$genus[otu_table$genus == 'g__Clostridium_T'] = 'g__Clostridium'
    otu_table$genus[otu_table$genus == 'g__Deinococcus_B'] = 'g__Deinococcus'
    otu_table$genus[otu_table$genus == 'g__Cereibacter_A_494292'] = 'g__Cereibacter'
    otu_table$genus[otu_table$genus == 'g__Enterococcus_H_360604'] = 'g__Enterococcus'
    otu_table$genus[otu_table$genus == 'g__Bacillus_P_294101'] = 'g__Bacillus'
    otu_table$genus[otu_table$genus == 'g__Pseudomonas_B_650326'] = 'g__Pseudomonas'
    otu_table$genus[otu_table$genus == 'g__Listeria_A'] = 'g__Listeria'
    
    otu_table$genus[is.na(otu_table$genus)] = 'Others'
    otu_table$genus[otu_table$genus == 'g__'] = 'Others'
    otu_table$genus[otu_table$genus == 'g__uncultured'] = 'Others'
    
    otu_table_genus = otu_table %>% group_by(genus) %>% 
      summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% as.data.frame()
    
    otu_table_genus <- column_to_rownames(otu_table_genus, "genus")
    #otu_table_genus = sweep(otu_table_genus,2,colSums(otu_table_genus),'/')
    otu_table_genus = subset(otu_table_genus, rowSums(otu_table_genus) > 0)
    
    colnames(otu_table_genus) = map_chr(colnames(otu_table_genus), function(x) strsplit(x, split = '_porechopped')[[1]][1])
    
    out_tables[[db]] = otu_table_genus}
  return(out_tables)
}

###############################################################################
# 0. First optimise the threshold value based on mock communities for the hyCAT 
# taxonomic assignment, then load data
alternat_thr_data = find_alternat_thr()

alternat_thr_data %>% group_by(Database, Threshold, Study) %>% 
  group_by(Database, Threshold) %>% summarise(mean=mean(Value), q25 = quantile(Value, probs = 0.25), q75 = quantile(Value, probs = 0.75)) %>%
  ggplot(aes(x=Threshold, colour=Database,y =mean)) + geom_smooth(method='loess') + geom_line() +
  geom_line(aes(x=Threshold, y = q25), linetype = 'dashed') + 
  geom_line(aes(x=Threshold, y = q75), linetype = 'dashed') + 
  geom_vline(xintercept = 1.15, colour='#1b9e77', linewidth=1.2) + geom_vline(xintercept = 1.2, colour='#7570b3', linewidth=1.2) + geom_vline(xintercept = 1.4, colour='#d95f02', linewidth=1.2) + 
  xlab('Coefficient for hyCAT approach') + ylab('Bray-Curtis dissimilarity') + theme_bw() + scale_color_brewer(palette = 'Dark2')
ggsave('figures/Figure_S1.pdf', width=9, height = 8) 

# Main script, first load data
parsed_a = load_study('study_a', 1.2, 1.15, 1.4)
parsed_b = load_study('study_b', 1.2, 1.15, 1.4)
parsed_c = load_study('study_c', 1.2, 1.15, 1.4)
parsed_d = load_study('study_d', 1.2, 1.15, 1.4)
parsed_e = load_study('study_e', 1.2, 1.15, 1.4)
parsed_f = load_study('study_f', 1.2, 1.15, 1.4)
parsed_g = load_study('study_g', 1.2, 1.15, 1.4)
reads_stats = read.csv('all_res/stats_reads.tsv', sep='\t')
reads_stats$Sample = map_chr(reads_stats$file, function(x) gsub('_chopped.fastq.gz', '', strsplit(x, split='/')[[1]][3]))
reads_stats$Sample = gsub('10_Strain_Even', 'X10_Strain_Even', reads_stats$Sample)
reads_stats$Sample = gsub('2A2_', 'X2A2_', reads_stats$Sample)
reads_stats$Sample = gsub('2AT_', 'X2AT_', reads_stats$Sample)
reads_stats$Sample = gsub('2B1_', 'X2B1_', reads_stats$Sample)
reads_stats$Sample = gsub('51_Kilo', 'X51_Kilo', reads_stats$Sample)
reads_stats$Sample = gsub('16S_', 'X16S_', reads_stats$Sample)

###############################################################################
# 1. assigned and Bacterial OTUs comparison
res_assigned = assigned_otus(parsed_a, parsed_b, parsed_c, parsed_d, parsed_e, parsed_f, parsed_g)
res_assigned$Method = factor(res_assigned$Method, levels = c('Centroid', 'mcCAT', 'hyCAT'))
res_bacterial = bacterial_otus(parsed_a, parsed_b, parsed_c, parsed_d, parsed_e, parsed_f, parsed_g)
res_bacterial$Method = factor(res_bacterial$Method, levels = c('Centroid', 'mcCAT', 'hyCAT'))

res_assigned$avg_qual = map_dbl(res_assigned$Sample, function(x) reads_stats$AvgQual[reads_stats$Sample == x])
res_bacterial$avg_qual = map_dbl(res_bacterial$Sample, function(x) reads_stats$AvgQual[reads_stats$Sample == x])
res_assigned$read_n = map_dbl(res_assigned$Sample, function(x) reads_stats$num_seqs[reads_stats$Sample == x])
res_bacterial$read_n = map_dbl(res_bacterial$Sample, function(x) reads_stats$num_seqs[reads_stats$Sample == x])

res_assigned$Study = map_chr(res_assigned$Sample, function(x) strsplit(reads_stats$file[reads_stats$Sample == x], split = '_results/')[[1]][1])
res_bacterial$Study = map_chr(res_bacterial$Sample, function(x) strsplit(reads_stats$file[reads_stats$Sample == x], split = '_results/')[[1]][1])

res_assigned$Database = factor(res_assigned$Database, levels = c('Silva', 'Greengenes 2', 'GTDB'))
res_bacterial$Database = factor(res_bacterial$Database, levels = c('Silva', 'Greengenes 2', 'GTDB'))
                                                              
a = res_assigned %>% filter(Metric == 'Proportion of assigned') %>% 
  ggplot(aes(x=Method, fill = Method, y=Value*100)) + geom_boxplot() + facet_grid(~Database) + geom_hline(yintercept = 100, linetype='dashed', colour='grey', linewidth=0.5) + 
  stat_compare_means(method = "kruskal.test", label.y = 53, size = 2) +  # global Kruskal–Wallis test
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods='holm',
                     comparisons = list(c("Centroid", "mcCAT"),c("Centroid", "hyCAT"),c("mcCAT", "hyCAT"))) +
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank()) + scale_y_continuous(limits = c(50, 116)) +
  xlab('') + ylab('Proportion of assigned OTUs [%]') + scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07"))
b = res_assigned %>% filter(Metric == 'Relative abundance of assigned') %>% 
  ggplot(aes(x=Method, fill=Method, y=Value*100)) + geom_boxplot() + facet_grid(~Database) + geom_hline(yintercept = 100, linetype='dashed', colour='grey', linewidth=0.5) + 
  stat_compare_means(method = "kruskal.test", label.y = 53, size = 2) +  # global Kruskal–Wallis test
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods='holm',
                     comparisons = list(c("Centroid", "mcCAT"),c("Centroid", "hyCAT"),c("mcCAT", "hyCAT"))) +
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank()) + scale_y_continuous(limits = c(50, 116)) +
  xlab('') + ylab('Proportion of assigned reads [%]') + scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07"))
c = res_bacterial %>% filter(Metric == 'Proportion of Bacteria') %>% 
  ggplot(aes(x=Method, fill=Method, y=Value*100)) + geom_boxplot() + facet_grid(~Database) + geom_hline(yintercept = 100, linetype='dashed', colour='grey', linewidth=0.5) + 
  stat_compare_means(method = "kruskal.test", label.y = 37, size = 2) +  # global Kruskal–Wallis test
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods='holm',
                     comparisons = list(c("Centroid", "mcCAT"),c("Centroid", "hyCAT"),c("mcCAT", "hyCAT"))) +
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank()) + scale_y_continuous(limits = c(35, 118)) +
  xlab('') + ylab('Proportion of bacterial OTUs [%]') + scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07"))
d = res_bacterial %>% filter(Metric == 'Relative abundance of Bacteria') %>% 
  ggplot(aes(x=Method, fill=Method, y=Value*100)) + geom_boxplot() + facet_grid(~Database) + geom_hline(yintercept = 100, linetype='dashed', colour='grey', linewidth=0.5) + 
  stat_compare_means(method = "kruskal.test", label.y = 37, size = 2) +  # global Kruskal–Wallis test
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods='holm',
                     comparisons = list(c("Centroid", "mcCAT"),c("Centroid", "hyCAT"),c("mcCAT", "hyCAT"))) +
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank()) + scale_y_continuous(limits = c(35, 118)) +
  xlab('') + ylab('Proportion of bacterial reads [%]') + scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07"))
ggarrange(a, b, c, d, nrow = 2, ncol = 2, labels = c('A', 'B', 'C', 'D'), common.legend = T, align = 'h')
ggsave('figures/Figure_2.pdf', width = 10, height = 7)

res_assigned %>% group_by(Metric, Method, Database) %>% summarise(median=median(Value), iqr=quantile(Value,probs=0.75)-quantile(Value,probs=0.25))
res_bacterial %>% group_by(Metric, Method, Database) %>% summarise(median=median(Value), iqr=quantile(Value,probs=0.75)-quantile(Value,probs=0.25))

###############################################################################
# 2. Assigned and Bacterial models with quality
get_regression_results <- function(metric, metric_axis, database, all_data,
                                   link = "logit", eps = 1e-6, level = 0.95){
  library(dplyr)
  library(ggplot2)
  library(betareg)
  library(tidyr)
  
  alpha <- 1 - level
  
  #---------- subset ----------
  df <- all_data %>%
    filter(Metric == metric, Database == database) %>%
    mutate(Method = droplevels(Method))
  
  if (!nrow(df)) stop("No data after filtering.")
  
  #---------- response in (0,1) ----------
  n <- nrow(df)
  df <- df %>%
    mutate(Value_beta = (Value * (n - 1) + 0.5) / n)  # Smithson-Verkuilen transform
  
  #---------- fit ----------
  fit <- betareg(Value_beta ~ 0 + Method:avg_qual, data = df, link = link)
  
  #---------- extract mean-model coefficients ----------
  beta_hat <- coef(fit, model = "mean")
  vc       <- vcov(fit)[names(beta_hat), names(beta_hat)]
  
  #---------- predictions with manual SE ----------
  newdat <- expand.grid(
    avg_qual = seq(min(df$avg_qual), max(df$avg_qual), length.out = 200),
    Method   = levels(df$Method)
  )
  X <- model.matrix(~ 0 + Method:avg_qual, data = newdat)
  eta <- as.vector(X %*% beta_hat)
  se_eta <- sqrt(diag(X %*% vc %*% t(X)))
  crit <- qnorm(1 - alpha/2)
  
  inv_link <- switch(link,
                     "logit"   = function(eta) 1 / (1 + exp(-eta)),
                     "cloglog" = function(eta) 1 - exp(-exp(eta)),
                     stop("Unsupported link")
  )
  
  newdat$pred  <- inv_link(eta)
  newdat$lower <- inv_link(eta - crit * se_eta)
  newdat$upper <- inv_link(eta + crit * se_eta)
  newdat$Database <- database
  
  #---------- line types ----------
  meth_levels <- levels(df$Method)
  line_types <- rep("solid", length(meth_levels))
  if (database %in% c("Silva","GTDB") && length(meth_levels) > 1) {
    line_types[length(meth_levels)] <- "dashed"
  }
  names(line_types) <- meth_levels
  
  #---------- plot ----------
  p <- ggplot() +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_point(data = df,
               aes(x = avg_qual, y = Value * 100, color = Method),
               alpha = 0.5,
               position = position_jitter(width = 0.2, height = 0)) +
    geom_line(data = newdat,
              aes(x = avg_qual, y = pred * 100, color = Method, linetype = Method),
              linewidth = 1, alpha = 0.9) +
    geom_ribbon(data = newdat,
                aes(x = avg_qual, ymin = lower * 100, ymax = upper * 100,
                    fill = Method, group = Method),
                alpha = 0.12, inherit.aes = FALSE) +
    facet_grid(~Database) +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "Sequencing quality [phred score]", y = metric_axis) +
    scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07")) +
    scale_colour_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07")) +
    scale_linetype_manual(values = line_types)
  
  #---------- return ----------
  list(
    plot = p,
    fit = fit,
    coef_table = summary(fit)$coefficients$mean,
    predictions = newdat
  )
}

# SILVA
silv_reg_propassigned = get_regression_results('Proportion of assigned', 'Proportion of assigned OTUs [%]', 'Silva', res_assigned)
gre2_reg_propassigned = get_regression_results('Proportion of assigned', 'Proportion of assigned OTUs [%]', 'Greengenes 2', res_assigned)
gtdb_reg_propassigned = get_regression_results('Proportion of assigned', 'Proportion of assigned OTUs [%]', 'GTDB', res_assigned)

leg_p = gre2_reg_propassigned$plot + theme(legend.position = "bottom")
legend <- ggpubr::get_legend(leg_p)
a = silv_reg_propassigned$plot + theme(legend.position = "none") + coord_cartesian(ylim = c(55,100))
b = gre2_reg_propassigned$plot + theme(legend.position = "none",
                                       axis.title.y = element_blank(),
                                       axis.text.y  = element_blank(),
                                       axis.ticks.y = element_blank()) + coord_cartesian(ylim = c(55,100))
c = gtdb_reg_propassigned$plot + theme(legend.position = "none",
                                       axis.title.y = element_blank(),
                                       axis.text.y  = element_blank(),
                                       axis.ticks.y = element_blank()) + coord_cartesian(ylim = c(55,100))
aligned <- align_plots(a, b, c, align = "v", axis = "l")
plots = plot_grid(aligned[[1]], aligned[[2]], aligned[[3]], ncol = 3, align = "hv", axis = "tblr", labels = c('A', 'B', 'C'))
ggarrange(plots, legend, ncol = 1, nrow = 2, heights = c(1,0.1))
ggsave('figures/Figure_3.pdf', width=9, height = 5)


silv_reg_relaassigned = get_regression_results('Relative abundance of assigned', 'Proportion of assigned reads [%]', 'Silva', res_assigned)
gre2_reg_relaassigned = get_regression_results('Relative abundance of assigned', 'Proportion of assigned reads [%]', 'Greengenes 2', res_assigned)
gtdb_reg_relaassigned = get_regression_results('Relative abundance of assigned', 'Proportion of assigned reads [%]', 'GTDB', res_assigned)

leg_p = gre2_reg_relaassigned$plot + theme(legend.position = "bottom")
legend <- ggpubr::get_legend(leg_p)
a = silv_reg_relaassigned$plot + theme(legend.position = "none") + coord_cartesian(ylim = c(55,100))
b = gre2_reg_relaassigned$plot + theme(legend.position = "none",
                                       axis.title.y = element_blank(),
                                       axis.text.y  = element_blank(),
                                       axis.ticks.y = element_blank()) + coord_cartesian(ylim = c(55,100))
c = gtdb_reg_relaassigned$plot + theme(legend.position = "none",
                                       axis.title.y = element_blank(),
                                       axis.text.y  = element_blank(),
                                       axis.ticks.y = element_blank()) + coord_cartesian(ylim = c(55,100))
aligned <- align_plots(a, b, c, align = "v", axis = "l")
plots = plot_grid(aligned[[1]], aligned[[2]], aligned[[3]], ncol = 3, align = "hv", axis = "tblr", labels = c('A', 'B', 'C'))
ggarrange(plots, legend, ncol = 1, nrow = 2, heights = c(1,0.1))
ggsave('figures/Figure_S2.pdf', width=9, height = 5)

sink("stats/Figure3_betareg_models_quality_propassignment.txt")
print('Silva - Prop assignment ~ Quality:')
silv_reg_propassigned$coef_table
print('Greengenes 2 - Prop assignment ~ Quality:')
gre2_reg_propassigned$coef_table
print('GTDB - Prop assignment ~ Quality:')
gtdb_reg_propassigned$coef_table
sink()

sink("stats/Figure_S2_betareg_models_quality_relaassignment.txt")
print('Silva - Rela assignment ~ Quality:')
silv_reg_relaassigned$mod
print('Greengenes 2 - Rela assignment ~ Quality:')
gre2_reg_relaassigned$mod
print('GTDB - Rela assignment ~ Quality:')
gtdb_reg_relaassigned$mod
sink()

# Bacterial (Proportion and Relative abundance) models
silv_reg_propbacterial = get_regression_results('Proportion of Bacteria', 'Proportion of Bacterial OTUs [%]', 'Silva', res_bacterial)
gre2_reg_propbacterial = get_regression_results('Proportion of Bacteria', 'Proportion of Bacterial OTUs [%]', 'Greengenes 2', res_bacterial)
gtdb_reg_propbacterial = get_regression_results('Proportion of Bacteria', 'Proportion of Bacterial OTUs [%]', 'GTDB', res_bacterial)

silv_reg_relabacterial = get_regression_results('Relative abundance of Bacteria', 'Proportion of Bacterial reads [%]', 'Silva', res_bacterial)
gre2_reg_relabacterial = get_regression_results('Relative abundance of Bacteria', 'Proportion of Bacterial reads [%]', 'Greengenes 2', res_bacterial)
gtdb_reg_relabacterial = get_regression_results('Relative abundance of Bacteria', 'Proportion of Bacterial reads [%]', 'GTDB', res_bacterial)

sink("stats/Table_S2_linear_models_quality_propbacterial.txt")
print('Silva - Prop bacterial ~ Quality:')
silv_reg_propbacterial$mod
print('Greengenes2 - Prop bacterial ~ Quality:')
gre2_reg_propbacterial$mod
print('GTDB - Prop bacterial ~ Quality:')
gtdb_reg_propbacterial$mod
sink()

sink("stats/Table_S2_linear_models_quality_relabacterial.txt")
print('Silva - Rela bacterial ~ Quality:')
silv_reg_relabacterial$mod
print('Greengenes 2 - Rela bacterial ~ Quality:')
gre2_reg_relabacterial$mod
print('GTDB - Rela bacterial ~ Quality:')
gtdb_reg_relabacterial$mod
sink()



###############################################################################
# 3. Taxonomic levels analysis
levels_a = get_proportion_assigned_to_level(parsed_a)
levels_a$Study = 'Study A'
levels_b = get_proportion_assigned_to_level(parsed_b)
levels_b$Study = 'Study B'
levels_c = get_proportion_assigned_to_level(parsed_c)
levels_c$Study = 'Study C'
levels_d = get_proportion_assigned_to_level(parsed_d)
levels_d$Study = 'Study D'
levels_f = get_proportion_assigned_to_level(parsed_f)
levels_f$Study = 'Study F'
levels_g = get_proportion_assigned_to_level(parsed_g)
levels_g$Study = 'Study G'

levels_all = rbind(levels_a, levels_b, levels_c, levels_d, levels_f, levels_g)
levels_all$Level[levels_all$Level == 'phylum'] = 'Phylum'
levels_all$Level[levels_all$Level == 'class'] = 'Class'
levels_all$Level[levels_all$Level == 'order'] = 'Order'
levels_all$Level[levels_all$Level == 'family'] = 'Family'
levels_all$Level[levels_all$Level == 'genus'] = 'Genus'
levels_all$Level <- factor(levels_all$Level, levels = c("Phylum", "Class", "Order", "Family", "Genus"))
levels_all$Method <- factor(levels_all$Method, levels = c("Centroid", "mcCAT", "hyCAT"))

levels_all$Database = factor(levels_all$Database, levels = c('Silva', 'Greengenes 2', 'GTDB'))

levels_all %>% filter(Metric == 'Relative abundance', Level != 'domain') %>% ggplot(aes(x=Method, fill=Method, y=Value*100)) + 
  geom_boxplot() + facet_grid(Database~Level, scales = 'free_y') + 
  stat_compare_means(method = "kruskal.test", label.y = -9, size = 2) +  # global Kruskal–Wallis test
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods='holm', size = 3, label.y = c(100, 110, 120),
                     comparisons = list(c("Centroid", "mcCAT"),c("Centroid", "hyCAT"),c("mcCAT", "hyCAT"))) +
  scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07")) + 
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank()) +
  xlab('') + ylab('Proportion of reads assigned to taxonomic level [%]') + 
  scale_y_continuous(limits = c(-10,127), breaks = c(0,25,50,75,100))
ggsave('figures/Figure_4.pdf', width = 10, height = 8)

levels_all %>% filter(Metric == 'Proportion', Level != 'domain') %>% ggplot(aes(x=Method, fill=Method, y=Value*100)) + 
  geom_boxplot() + facet_grid(Database~Level, scales = 'free') + 
  stat_compare_means(ref.group = 'Centroid',label = "p.signif", paired = T, label.y = 105) + scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods='holm', size = 3,
                     comparisons = list(c("Centroid", "mcCAT"),c("Centroid", "hyCAT"),c("mcCAT", "hyCAT"))) +
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank()) +
  xlab('') + ylab('Proportion of OTUs assigned to taxonomic level [%]') + 
  scale_y_continuous(limits = c(-10,135), breaks = c(0,25,50,75,100))
ggsave('figures/Figure_S3.pdf', width = 8, height = 5)


###############################################################################
# 4. Mock communities benchmark analysis
parsed_CONCOMPRA = load_CONCOMPRA()
genus_CONCOMPRA = summarise_CONCOMPRA(parsed_CONCOMPRA)

mock_data = get_mock_results(parsed_a, parsed_b, parsed_c, parsed_d, parsed_e, parsed_f, parsed_g, genus_CONCOMPRA)
mock_data$Method = factor(mock_data$Method, levels = c('Centroid', 'mcCAT', 'hyCAT', 'CONCOMPRA'))
mock_data$avg_qual = map_dbl(mock_data$Sample, function(x) reads_stats$AvgQual[reads_stats$Sample == x])
mock_data$read_n = map_dbl(mock_data$Sample, function(x) reads_stats$num_seqs[reads_stats$Sample == x])
mock_data$Database[mock_data$Database == 'silva'] = 'Silva'
mock_data$Database[mock_data$Database == 'greengenes'] = 'Greengenes 2'
mock_data$Database[mock_data$Database == 'gtdb'] = 'GTDB'

mock_data$Database = factor(mock_data$Database, levels = c('Silva', 'Greengenes 2', 'GTDB'))
mock_data = mock_data %>% remove_rownames() %>% distinct()

a = mock_data %>% filter(Metric == 'Bray-Curtis dissimilarity') %>%
  ggplot(aes(x=Method, y=1-Value, fill=Method)) + geom_boxplot() + 
  stat_compare_means(method = "kruskal.test", label.y = -0.09, size = 2) +  # global Kruskal–Wallis test
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.methods='holm', size = 3,
                     comparisons = list(c("Centroid", "mcCAT"),c("Centroid", "hyCAT"),c("Centroid",'CONCOMPRA'),
                                        c("mcCAT", "hyCAT"), c('mcCAT', "CONCOMPRA"), c('hyCAT', "CONCOMPRA"))) +
  facet_grid(~Database) + 
  theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank(), legend.position = 'top') +
  xlab('') + ylab('Bray-Curtis similarity') + scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07","#E7B800"))

b_dat <- mock_data %>%
  dplyr::filter(Metric == 'Counts expected') %>%
  dplyr::mutate(prop_expected = Value) 
b_dat$Metric = "Proportion expected reads"

b <-  ggplot(data = b_dat, mapping = aes(x = Method, y = prop_expected, fill = Method)) +
  geom_boxplot(width = 0.6, colour = "black", show.legend = FALSE) + facet_grid(~Database) +
  ggpubr::stat_compare_means(method = "kruskal.test", label.y =-1000, size = 2, paired = T) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", p.adjust.method='holm', paired = T, size=3,
                     comparisons = list(c("Centroid", "mcCAT"),c("Centroid", "hyCAT"),c("Centroid",'CONCOMPRA'),
                                        c("mcCAT", "hyCAT"), c('mcCAT', "CONCOMPRA"), c('hyCAT', "CONCOMPRA"))) +
  theme_linedraw() +
  ylab("Number of reads assigned to expected genera [#]") + xlab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07","#E7B800"))

ggarrange(a, b, ncol = 1, nrow = 2, labels = c('A', 'B'), common.legend = T, align = 'v')
ggsave('figures/Figure_5.pdf', width = 7, height = 8)

mock_data %>% filter(Metric %in% c('Counts expected', 'Bray-Curtis dissimilarity')) %>% 
  group_by(Method, Database, Metric) %>% summarise(median = median(Value),
                                                   iqr = quantile(Value,probs=0.75)-quantile(Value,probs=0.25))  %>% arrange(Metric) %>% print(n=30)


get_regression_results_qual <- function(metric, metric_axis, database, all_data,
                                   link = "logit", eps = 1e-6, level = 0.95){
  library(dplyr)
  library(ggplot2)
  library(betareg)
  library(tidyr)
  
  alpha <- 1 - level
  
  #---------- subset ----------
  df <- all_data %>%
    filter(Metric == metric, Database == database) %>%
    mutate(Method = droplevels(Method))
  
  if (!nrow(df)) stop("No data after filtering.")
  
  #---------- response in (0,1) ----------
  n <- nrow(df)
  df <- df %>%
    mutate(Value_beta = (Value * (n - 1) + 0.5) / n)  # Smithson-Verkuilen transform
  
  #---------- fit ----------
  fit <- betareg(Value_beta ~ 0 + Method:avg_qual, data = df, link = link)
  
  #---------- extract mean-model coefficients ----------
  beta_hat <- coef(fit, model = "mean")
  vc       <- vcov(fit)[names(beta_hat), names(beta_hat)]
  
  #---------- predictions with manual SE ----------
  newdat <- expand.grid(
    avg_qual = seq(min(df$avg_qual), max(df$avg_qual), length.out = 200),
    Method   = levels(df$Method)
  )
  X <- model.matrix(~ 0 + Method:avg_qual, data = newdat)
  eta <- as.vector(X %*% beta_hat)
  se_eta <- sqrt(diag(X %*% vc %*% t(X)))
  crit <- qnorm(1 - alpha/2)
  
  inv_link <- switch(link,
                     "logit"   = function(eta) 1 / (1 + exp(-eta)),
                     "cloglog" = function(eta) 1 - exp(-exp(eta)),
                     stop("Unsupported link")
  )
  
  newdat$pred  <- inv_link(eta)
  newdat$lower <- inv_link(eta - crit * se_eta)
  newdat$upper <- inv_link(eta + crit * se_eta)
  newdat$Database <- database
  
  #---------- line types ----------
  meth_levels <- levels(df$Method)
  line_types <- rep("solid", length(meth_levels))
  if (database %in% c("Silva","GTDB") && length(meth_levels) > 1) {
    line_types[length(meth_levels)] <- "dashed"
  }
  names(line_types) <- meth_levels
  
  #---------- plot ----------
  p <- ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_point(data = df,
               aes(x = avg_qual, y = Value , color = Method),
               alpha = 0.5,
               position = position_jitter(width = 0.2, height = 0)) +
    geom_line(data = newdat,
              aes(x = avg_qual, y = pred , color = Method),
              linewidth = 1, alpha = 0.9) +
    geom_ribbon(data = newdat,
                aes(x = avg_qual, ymin = lower , ymax = upper ,
                    fill = Method, group = Method),
                alpha = 0.12, inherit.aes = FALSE) +
    facet_grid(~Database) +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "Sequencing quality [phred score]", y = metric_axis) +
    scale_fill_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07","#E7B800")) +
    scale_colour_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07","#E7B800")) 
  
  #---------- return ----------
  list(
    plot = p,
    fit = fit,
    coef_table = summary(fit)$coefficients$mean,
    predictions = newdat
  )
}

mock_bs_data = mock_data %>% filter(Metric == 'Bray-Curtis dissimilarity') %>% mutate(Value=1-Value)
silv_reg_bc = get_regression_results_qual('Bray-Curtis dissimilarity', 'Bray-Curtis similarity', 'Silva', mock_bs_data)
silv_reg_bc$plot
ggsave('figures/Figure_S4.pdf', width = 6, height = 5)



###############################################################################
# 4. Compare mock communities and quality/read number for CAT and CONCOMPRA using Silva
expected_comp = expand.grid(Sample = unique(mock_data$Sample), Database = unique(mock_data$Database), Method = unique(mock_data$Method))
expected_comp$OTU_number = map_int(1:nrow(expected_comp), function(i) mock_data %>% filter(Metric == 'OTU number', Database == expected_comp$Database[i],
                                                                                          Sample == expected_comp$Sample[i], Method == expected_comp$Method[i]) %>% pull(Value))
expected_comp$bc_sim = map_dbl(1:nrow(expected_comp), function(i) 1 - mock_data %>% filter(Metric == 'Bray-Curtis dissimilarity', Database == expected_comp$Database[i],
                                                                                           Sample == expected_comp$Sample[i], Method == expected_comp$Method[i]) %>% pull(Value))
expected_comp$Expected_richness = map_int(1:nrow(expected_comp), function(i) mock_data %>% filter(Metric == 'Expected richness', Database == expected_comp$Database[i],
                                                                                           Sample == expected_comp$Sample[i], Method == expected_comp$Method[i]) %>% pull(Value))
expected_comp$Spurious_richness = map_int(1:nrow(expected_comp), function(i) mock_data %>% filter(Metric == 'Number unexpected', Database == expected_comp$Database[i],
                                                                                                  Sample == expected_comp$Sample[i], Method == expected_comp$Method[i]) %>% pull(Value))

expected_comp$Expected_count = map_int(1:nrow(expected_comp), function(i) mock_data %>% filter(Metric == 'Counts expected', Database == expected_comp$Database[i],
                                                                                                  Sample == expected_comp$Sample[i], Method == expected_comp$Method[i]) %>% pull(Value))
expected_comp$Spurious_count = map_int(1:nrow(expected_comp), function(i) mock_data %>% filter(Metric == 'Counts spurious', Database == expected_comp$Database[i],
                                                                                                  Sample == expected_comp$Sample[i], Method == expected_comp$Method[i]) %>% pull(Value))
expected_comp$Spurious_proportion = expected_comp$Spurious_count / (expected_comp$Spurious_count + expected_comp$Expected_count)

expected_comp$Read_number = map_int(1:nrow(expected_comp), function(i) mock_data %>% filter(Metric == 'Number unexpected', Database == expected_comp$Database[i],
                                                                                      Sample == expected_comp$Sample[i], Method == expected_comp$Method[i]) %>% pull(read_n))
expected_comp$avg_qual = map_dbl(1:nrow(expected_comp), function(i) mock_data %>% filter(Metric == 'Number unexpected', Database == expected_comp$Database[i],
                                                                                            Sample == expected_comp$Sample[i], Method == expected_comp$Method[i]) %>% pull(avg_qual))
expected_comp$Study = map_chr(1:nrow(expected_comp), function(i) mock_data %>% filter(Metric == 'Number unexpected', Database == expected_comp$Database[i],
                                                                                                   Sample == expected_comp$Sample[i], Method == expected_comp$Method[i]) %>% pull(Study))

expected_comp$sample_weights = map_dbl(expected_comp$Study, function(x) 1 / sum(expected_comp$Study == x))


library(repmod)
sink('stats/Figure_6_models.txt')
a = ggplot(expected_comp %>% filter(Database == 'Silva'), aes(x=avg_qual, colour = Method, y = bc_sim)) +
  geom_point() + geom_smooth(method = 'rlm', method.args = list(maxit = 100)) +
  theme_linedraw() +
  ylab("Bray-Curtis similarity") + xlab("Average sequencing quality [phred score]") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07","#E7B800")) + coord_cartesian(ylim=c(0, 1))
model <- rlm(bc_sim ~ Method*avg_qual, data = expected_comp %>% filter(Database == 'Silva'))
summary(model)
rob.pvals(model)

b = ggplot(expected_comp %>% filter(Database == 'Silva'), aes(x=Read_number, colour = Method, y = bc_sim)) +
  geom_point() + geom_smooth(method = 'rlm', method.args = list(maxit = 100)) +
  theme_linedraw() +
  ylab("Bray-Curtis similarity") + xlab("Number of raw reads [#]") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07","#E7B800")) + coord_cartesian(ylim=c(0, 1))
model <- rlm(bc_sim ~ Method*Read_number, data = expected_comp %>% filter(Database == 'Silva'))
summary(model)
rob.pvals(model)

c = ggplot(expected_comp %>% filter(Database == 'Silva'), aes(x=avg_qual, colour = Method, y = Expected_count)) +
  geom_point() + geom_smooth(method = 'rlm', method.args = list(maxit = 100)) +
  theme_linedraw() +
  ylab("Number of reads assigned to expected genera [#]") + xlab("Average sequencing quality [phred score]") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07","#E7B800")) + coord_cartesian(ylim=c(0, 16500))
model <- rlm(Expected_count ~ Method*avg_qual, data = expected_comp %>% filter(Database == 'Silva'))
summary(model)
rob.pvals(model)

d = ggplot(expected_comp %>% filter(Database == 'Silva'), aes(x=avg_qual, colour = Method, y = Spurious_proportion*100)) +
  geom_point() + geom_smooth(method = 'rlm', method.args = list(maxit = 100)) +
  theme_linedraw() +
  ylab("Proportion of reads assigned to spurious genera [%]") + xlab("Average sequencing quality [phred score]") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = c("#7E7E61", "#00AFBB", "#FC4E07","#E7B800")) + coord_cartesian(ylim=c(0, 100))
model <- rlm(Spurious_proportion*100 ~ Method*avg_qual, data = expected_comp %>% filter(Database == 'Silva'))
summary(model)
rob.pvals(model)

ggarrange(a,b,c,d, labels = c('A','B','C','D'), ncol = 2, nrow = 2, common.legend = T)
ggsave('figures/Figure_6.pdf', width = 9.2, height = 9.2)
sink()
