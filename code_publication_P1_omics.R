## In your working directory, there should be a 'data' folder with all the datasets
## 'PG_data' folder: all the pre-generated data
## 'PG_data/GSEA_result' folder: all the GSEA output data

## This is the first part which is mainly for omics analysis (CCLE)
## The second part (another R file) contains the code for machine learning and clinical prediction (survival & therapy response)
## The first part contains necessary code to generate the data for the second part

##------- loading the package -------
## gene name and ID annotation
library(org.Hs.eg.db)
library(biomaRt)  ## for gene name/id annotation
## Tidyverse and supporting packages for Tidyverse
library(dplyr); library(tidyr); library(tibble); library(stringr); library(readr)
library(ggplot2); library(ggpubr); library(ggrepel)
library(patchwork)  ## plot multiple plot together by direct adding
library(gridExtra)  ## multiple plots
library(scales)  ## axis scales
library(GGally); library(extrafont)
## Machine Learning
library(caret) ## ML assistant tool
## survival analysis
library(survival); library(survminer)


##------- setting the working directory ---------
setwd("~/Your_wd") ## change Your_wd into the designated working directory
## CODE CHECK
setwd("/Users/ofeklin/Library/CloudStorage/OneDrive-Technion/Cancer/CCLE/thesis_code_wd")

##------- defining function ---------
## cor.test_error2na_miss_observation(x, y, method) helps to escape the error from insufficient observations from the table, and interpret this error into returned NA
cor.test_error2na_miss_observation <- function (x, y, method){
  if (sum(!is.na(x)) <= 2 | sum(!is.na(y)) <= 2) {return(NA)}
  else {cor.test(x,y,method = method)}
}

## auto generate label location in ggplot
axis_location <- function(data, x, y, location = c("leftup", "rightup", "leftdown", "rightdown", "leftmid", "upmid", "rightmid", "downmid")){
  ## x, y are number; should be the column number of x and y data
  ## location are the location you want in the plot
  x_max <- max(data[,x])*0.9; x_min <- min(data[,x])*1.1; x_avg <- mean(data[,x]); y_max <- max(data[,y])*0.9; y_min <- min(data[,y])*1.1; y_avg <- mean(data[,y])
  if(location=="leftup"){
    hold <- c(x_min,y_max)
  } else if(location=="rightup"){
    hold <- c(x_max,y_max)
  } else if(location=="leftdown"){
    hold <- c(x_min,y_min)
  } else if(location=="rightdown"){
    hold <- c(x_max,y_min)
  } else if(location=="leftmid"){
    hold <- c(x_min,y_avg)
  } else if(location=="upmid"){
    hold <- c(x_avg,y_max)
  } else if(location=="rightmid"){
    hold <- c(x_max,y_avg)
  } else if(location=="downmid"){
    hold <- c(x_avg,y_min)
  } else {stop("Please input correct location.")}
  return(hold)
}

## PrankGSEA create a data table used for GSEA Prank analysis and (optional) output an .rnk.txt file
## which is composed of (gene) name column and parameter (e.g. q-value) column
## provide filename, if .rnk.txt output is desired
PrankGSEA <- function(data, name_colnumber, parameter_colnumber, cutoff=0.05, filename=NA) {
  ## data - the input of result data table
  ## name_colnumber - columne number where the gene name is located
  ## parameter_colnumber - columne number where the ranking parameter is located (e.g. sorted q-value)
  ## cutoff - cutoff value for parameter value (if alternative="pvalue", parameter <cutoff will be included; 
  ##                                            if alternative="estimate", abs(parameter)>cutoff will be included)
  ## alternative - "pvalue" for pvalue parameter; "estimate" for cor/rho/tau estimates
  ## filename - output file name without postfix
  data <- data[,c(name_colnumber,parameter_colnumber)]
  ## if no cutoff is required, set p_cutoff with any number > 1, and est_cutoff with any number < 0
  data <- data[data[,2]<cutoff & !is.na(data[,2]),]
  data <- data[order(data[,2]),]
  
  if (is.na(filename)) {  ## don't output file
    return(data)
  } else {  ## output file without header
    write.table(data, paste(filename,".rnk.txt",sep = ""), sep = "\t", row.names = F, col.names = F)
    return(data)
  }
}

## for preranked list with positive/negative correlation coefficient, sorted according to p-value
PrankGSEA.est <- function(data, name_colnumber, sort_colnumber, parameter_colnumber, p_cutoff = 0.05, est_cutoff=0.05, filename=NA) {
  data <- data[,c(name_colnumber,sort_colnumber,parameter_colnumber)]
  data <- data[order(data[,2]),]  ## sort according to p-value
  ## if no cutoff is required, set p_cutoff with any number > 1, and est_cutoff with any number < 0
  data <- data[data[,2]<p_cutoff & abs(data[,3])>est_cutoff & (!is.na(data[,2])&!is.na(data[,3])),]  ## apply cutoff
  data <- data[,c(1,3)]  ## only keep the name and parameter
  
  if (is.na(filename)) {  ## don't output file
    return(data)
  } else {  ## output file without header
    write.table(data, paste(filename,".rnk.txt",sep = ""), sep = "\t", row.names = F, col.names = F)
    return(data)
  }
}

## for DataTable creation
DT.generator <- function(data, meta, meta_col = c("DEPMAPID","PRIMARY_SITE","DOUBLING_TIME"), prof.name = "DOUBLING_TIME") {
  ## data is (d.f.) expression data. Has "DEPMAPID" column. Has gene on each column
  ## meta is (d.f.) sample metadata (usually named "DEPMAPID_DT"). Has column called "DEPMAPID". Has column called "DOUBLING_TIME"
  ## meta_col (index/name vector) has the index/name of column(s) from meta table
  ## prof.name (chr.) is used to specify the column containing the proliferation (DT or GR) data; make a change if it is not "DOUBLING_TIME"
  library(dplyr)
  dt <- meta[,meta_col] %>% right_join(data, by = "DEPMAPID") 
  dt <- dt[!is.na(dt[,prof.name]),]
  return(dt)
}

## used for dplyr::select() to select desired columns
col.selector <- function(order_df, col_order, order, custom = c(), decreasing = F, colnum_id = 1) {
  ## order_df is (d.f.) scoring data (e.g. summary_RNA_expression). Has id column (name or numeric id). Has score column (p value/q value).
  ## col_order is (string or number) the name/number of column whose scores are used to select the genes (usually one of the p values or q values)
  ## order (value series) is the range of value used to select the genes (e.g. for top 300, input 1:300)
  ## custom (string) is custom column name, and those column will be included at the head (e.g. "DEPMAPID", "DOUBLING_TIME")
  ## decreasing is parameter used in order(). colnum_id used to select common names to match data and meta
  ind_col <- c(custom, order_df[order(order_df[,col_order], decreasing = decreasing)[order], colnum_id])
  return(ind_col)
}


## survival analysis
## make sure to switch part of the code (silence and unsilence) when applying to TCGA or GENT2
SurvivalData.generator <- function(srvl, col_pDT, col_fustat, col_checkNA = c("DSS_STATUS","DSS_MONTHS","PFS_STATUS","PFS_MONTHS","OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS"), tail_quantile = F, 
                                   custom_fustat = F,custom_fustat_label = c("Title","Low","High")) {
  ## a version 2 function - can handle MKI67 and other biomarkers besides pDT
  ## the function helps to generate labels on pDT (high/low) ($pDT_Group) and fustat labels ($fustat)
  ## srvl is your input data with pDT, OS, DSS, DFS data
  ## col_pDT is the column number where pDT is located
  ## col_fustat is the column number where fustat/DSS_STATUS is located
  ## col_checkNA indicates in which column the NA or empty cases are going to be excluded - number or colname
  ## tail_quantile is the command to take head and tail quantile or not: F for taking median; T for taking quantile
  
  ## input check
  if(!all(c("OS_STATUS", "OS_MONTHS", "DSS_STATUS", "DSS_MONTHS", "DFS_STATUS", "DFS_MONTHS", "PFS_STATUS", "PFS_MONTHS") %in% colnames(srvl))) {
    stop("Make sure your input has: OS, DSS, PFS with time and status columns - in total 6 columns")
  }
  library(dplyr)
  library(survival)
  library(survminer)
  
  colnames(srvl)[col_pDT] <- "pDT"  ## for code consistency
  
  if(typeof(col_checkNA) == "double" | typeof(col_checkNA) == "integer") {
    col_checkNA <- c("DSS_STATUS","DSS_MONTHS","PFS_STATUS","PFS_MONTHS","OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS")[col_checkNA]
  } else if (typeof(col_checkNA) == "character") { } else {stop("Provide proper col_checkNA values!")}
  
  ## exclude the case with missing information
  for (CNA in col_checkNA) {
    srvl <- srvl %>% dplyr::filter(!is.na(!!as.symbol(CNA)) & nchar(!!as.symbol(CNA)) != 0)
  }
  # srvl <- srvl %>% dplyr::filter(!is.na(PFS_STATUS) & !is.na(PFS_MONTHS) & nchar(PFS_STATUS) != 0 & nchar(PFS_MONTHS) != 0) %>%
  #   dplyr::filter(!is.na(DSS_STATUS) & !is.na(DSS_MONTHS) & nchar(DSS_STATUS) != 0 & nchar(DSS_MONTHS) != 0)
  
  ## mark cases according to level of pDT
  srvl <- srvl %>% mutate(pDT_Group = ifelse(pDT >= median(srvl$pDT), 1, 0))
  if(tail_quantile) {  ## if taking quantile is TRUE
    srvl <- srvl[srvl$pDT <= quantile(srvl$pDT)["25%"] | srvl$pDT >= quantile(srvl$pDT)["75%"],]
  }
  if(!custom_fustat){
    srvl$pDT_Group <- factor(srvl$pDT_Group, levels = c(0,1), labels = c("Low","High"))
  } else if(custom_fustat) {
    srvl$pDT_Group <- factor(srvl$pDT_Group, levels = c(0,1), labels = custom_fustat_label[2:3]); colnames(srvl)[colnames(srvl)=="pDT_Group"]<-custom_fustat_label[1]
  }
  
  # ## for TCGA ##
  # srvl$fustat <- substr(srvl[,col_fustat],1,1) %>% as.numeric() #%>% as.factor()
  ## for GENT2 ##
  if (colnames(srvl)[col_fustat] != "DFS_STATUS") {
    srvl <- srvl %>% mutate(fustat = ifelse(!!as.symbol(colnames(srvl)[col_fustat]) == "Dead", 1, 0))
  } else if (colnames(srvl)[col_fustat] == "DFS_STATUS") {
    srvl <- srvl %>% mutate(fustat = ifelse(!!as.symbol(colnames(srvl)[col_fustat]) == "Recurrence", 1, 0))
  } else {stop("survival type can't be determined => fustat can't be created.")}
  
  cat("0 and 1 for low pDT and high pDT, respectively; fustat label is consistent with input data.\n")
  return(srvl)
}


## for fixing CCLE and TCGA gene identification compatiblity
## identify missing genes
missing.gene.compatible <- function (table1, table2, significant = "", cutoff = "") {
  ## table1 could be your meta_gene_GExp for CCLE data or your correlation result
  ## table2 could be your TCGA table
  ## table1 will be mapped with table2 according to Entrez ID first, then HUGO symbol
  ## make sure your table1 and table2 have: Hugo_Symbol, Entrez_Gene_id
  ## significant = c("pvalue_pearson", "pvalue_spearman"). Define your significance and scores should be available. If not provided, all the genes will be executed
  ## cutoff is your cutoff on signficant
  ## output returns "missing" includes all the observations/rows from table1 that is missing in table2 (can't be found in either name or ID)
  
  ## change name from different version of codes in the case of table1 <- summary_RNA_expression
  if ("Gene_Transcript" %in% colnames(table1)) {colnames(table1)[colnames(table1) == "Gene_Transcript"] <- "Hugo_Symbol"}
  
  ## check required columns from input tables
  if (!all(c("Hugo_Symbol","Entrez_Gene_Id") %in% colnames(table1))) {stop("Hugo or Entrez is missing in table1")
  } else if (!all(c("Hugo_Symbol","Entrez_Gene_Id") %in% colnames(table2))) {stop("Hugo or Entrez is missing in table2")}
  
  ##------- prepare ------
  colnames(table1)  ## make sure your table1 has: Hugo_Symbol, Entrez_Gene_id, Hugo_Entrez
  if (is.null(table1$Hugo_Entrez)) {table1$Hugo_Entrez <- paste(table1$Hugo_Symbol, table1$Entrez_Gene_Id, sep = "_")}  ## if not, use this to create combined gene name/id
  if (is.null(table2$Hugo_Entrez)) {table2$Hugo_Entrez <- paste(table2$Hugo_Symbol, table2$Entrez_Gene_Id, sep = "_")}  ## if not, use this to create combined gene name/id
  
  ## a rank according to Spearman p value
  if (significant != "") {table1$order_sign <- match(table1[,significant],sort(table1[,significant], na.last = T))}
  
  ##------- data compatibility ---------
  table(table1$Hugo_Entrez %in% table2$Hugo_Entrez)  ## name_ID can be found in table2
  table(table1$Entrez_Gene_Id %in% table2$Entrez_Gene_Id)  ## ID can be found in table2
  table(table1$Hugo_Symbol %in% table2$Hugo_Symbol)  ## name can be found in table2
  
  ## index indicating genes can be found by which way
  ind_gene <- data.frame(Entrez_Gene_Id = table1$Entrez_Gene_Id %in% table2$Entrez_Gene_Id, Hugo_Symbol = table1$Hugo_Symbol %in% table2$Hugo_Symbol)
  ## genes in correlation list are missing in table2 and can't be found by Entrez ID or Hugo symbol
  missing <- table1[sapply(1:nrow(ind_gene), function(X) !any(ind_gene$Entrez_Gene_Id[X],ind_gene$Hugo_Symbol[X])), ]
  # missing <- table1[!(table1$Entrez_Gene_Id %in% table2$Entrez_Gene_Id) ,]
  if (nrow(missing) == 0) {stop("No missing genes found.")}
  
  if (significant != "") {
    missing <- missing[order(missing$pvalue_spearman),]  ## recorder according to p value
    missing <- missing[missing$pvalue_spearman<0.05 & !is.na(missing$pvalue_spearman),]  ## only care about the significant and non-NA
  }
  
  ## should be expected that genes can't be found in table2 in either name or ID
  if (any(sapply(1:nrow(missing), function(X) missing$Hugo_Symbol[X] %in% table2$Hugo_Symbol | missing$Entrez_Gene_Id[X] %in% table2$Entrez_Gene_Id))) {stop("genes from missing list can't be mapped into table2")}
  
  return(missing)
}

gene.identif.mapping <- function (table1, table2, suffix.table1, suffix.table2) {
  ## same requirement as missing.gene.compatible
  ## in this case, all the genes from table1 should be found in table2 either by symbol or id
  ## output returns "gene_table1_table2" - the gene lists with common identifications
  
  ## change name from different version of codes in the case of table1 <- summary_RNA_expression
  if ("Gene_Transcript" %in% colnames(table1)) {colnames(table1)[colnames(table1) == "Gene_Transcript"] <- "Hugo_Symbol"}
  
  ## check required columns from input tables
  if (!all(c("Hugo_Symbol","Entrez_Gene_Id") %in% colnames(table1))) {stop("Hugo or Entrez is missing in table1")
  } else if (!all(c("Hugo_Symbol","Entrez_Gene_Id") %in% colnames(table2))) {stop("Hugo or Entrez is missing in table2")}
  
  ## common gene id
  gene_table1_table2 <- table1[,c("Hugo_Symbol","Entrez_Gene_Id","Hugo_Entrez")]
  table(gene_table1_table2$Entrez_Gene_Id %in% table2$Entrez_Gene_Id)  ## genes can be matched by ID
  gene_table1_table2 <- left_join(gene_table1_table2, table2[,c("Hugo_Symbol","Entrez_Gene_Id")], 
                                  by = c("Entrez_Gene_Id" = "Entrez_Gene_Id"), keep = T, suffix = c(".table1",".table2"))
  
  ind <- which(is.na(gene_table1_table2$Hugo_Symbol.table2))
  gene_table1_table2[ind,] <- left_join(gene_table1_table2[ind,c("Hugo_Symbol.table1","Entrez_Gene_Id.table1","Hugo_Entrez")], 
                                        table2[,c("Hugo_Symbol","Entrez_Gene_Id")],by = c("Hugo_Symbol.table1" = "Hugo_Symbol"), 
                                        keep = T, suffix = c(".table1",".table2"))
  
  ## double check to prevent bug
  if (!all(gene_table1_table2$Hugo_Entrez %in% table1$Hugo_Entrez)) {stop("genes in gene_table1_table2 are missing in table1")}
  
  if (exists("suffix.table1")) {colnames(gene_table1_table2) <- gsub(".table1", suffix.table1, colnames(gene_table1_table2))}
  if (exists("suffix.table2")) {colnames(gene_table1_table2) <- gsub(".table2", suffix.table2, colnames(gene_table1_table2))}
  
  return(gene_table1_table2)
}

## add.inform.col(origin_table, supplement, colnumber_key_origin, colnumber_key_supplement) 
##  adds the supplementary information to the original table, via matching keys, in column order
add.inform.col <- function (origin_table, supplement, colnumber_key_origin, colnumber_key_supplement) {
  ## the origin_table and supplement should include the column containing the keys to match
  ## origin_table should be complete
  ## supplement should only contain the desired column
  ## column number is related to respective table - either original table or supplement table
  col_number <- ncol(origin_table)
  supplement_exclude_keyColumn <- as.data.frame(supplement[,-colnumber_key_supplement])
  col_added_number <- ncol(supplement_exclude_keyColumn)
  
  i = 1
  while (i <= col_added_number) {origin_table <- data.frame(origin_table,NA); i<-i+1}
  colnames(origin_table)[(col_number+1):(col_number+col_added_number)] <- colnames(supplement_exclude_keyColumn)
  origin_table[,(col_number+1):(col_number+col_added_number)] <- supplement_exclude_keyColumn[match(origin_table[,colnumber_key_origin],supplement[,colnumber_key_supplement]),]
  return(origin_table)
}


##------- loading data ----------
## cancer gene annotation 
CancerGene <- read.csv("data/Cancer_Gene_GRCh38_COSMIC_v96.csv")  ## data loading
## clinical data
DEPMAPID_DT <- read.delim("data/data_clinical_sample.txt", skip = 4, header = 4)
## to translate the cell line ID or cell name into doubling time
DEPMAPID_DT <- DEPMAPID_DT[,c("SAMPLE_ID", "DEPMAPID","DOUBLING_TIME", "DOUBLING_TIME_FROM_VENDOR", "PRIMARY_SITE")]   ## with doubling time according to DEPMAPID

## mutation data - CCLE
mutations <- read.delim("data/data_mutations.txt",skip = 2)
## clinical data
RAW_clinical <- read.delim("data/data_clinical_sample.txt", skip = 4, header = 4)

## CODE CHECKING
# RAW_clinical <- read.delim("/Users/ofeklin/OneDrive - Technion/Cancer/CCLE/[Database]ccle_broad_2019/data_clinical_sample.txt", skip = 4, header = 4)
# mutations <- read.delim("/Users/ofeklin/Library/CloudStorage/OneDrive-Technion/Cancer/CCLE/[Database]ccle_broad_2019/data_mutations.txt",skip = 2) ## for code testing
summmary_mutation_2check <- read.csv("/Users/ofeklin/Library/CloudStorage/OneDrive-Technion/Cancer/CCLE/[Database]ccle_broad_2019/WorkingSpace/DT_VS_GeneticsAlternation/table/SNV_result_ttest_wilcox.csv")
summmary_cnv_2check <- read.csv("/Users/ofeklin/Library/CloudStorage/OneDrive-Technion/Cancer/CCLE/CNV/table/CNV_result.csv")
summmary_express_2check <- read.csv("/Users/ofeklin/Library/CloudStorage/OneDrive-Technion/Cancer/CCLE/RNA_expression/table/expression_result.csv")

## CNV - CCLE
CNV <- read.csv("data/CCLE_gene_cn.csv", check.names = F)

## RNA expression - CCLE
RNA_expression <- read.csv("data/CCLE_expression.csv", check.names = F)
rownames(RNA_expression) <- RNA_expression[,1]
colnames(RNA_expression)[1] <- "DEPMAPID"
## Gene Hugo symbols and Entrez IDs for CCLE data
meta_gene_GExp <- data.frame(colnames = colnames(RNA_expression)[-1]) %>% separate(colnames, into = c("Hugo_Symbol", "Entrez_Gene_Id"), sep = " \\(")
meta_gene_GExp$Entrez_Gene_Id <- as.numeric(sub("\\).*","", meta_gene_GExp$Entrez_Gene_Id))  ## purify the id
## change with proper names
colnames(RNA_expression) <- sub(" \\(", "_", colnames(RNA_expression))
colnames(RNA_expression) <- sub("\\)", "", colnames(RNA_expression))
colnames(RNA_expression)[-1] <- sapply(2:ncol(RNA_expression), function(X) sub(" \\(.*", "", colnames(RNA_expression)[X]))

## pre-generated RNA expression correlation list (CCLE), which can be generated by the correlation analysis code. It is needed for Machine Learning.
summary_RNA_expression <- read.csv("PG_data/summary_RNA_expression.csv") ## whole CL

## Drug Sensitivity
drugsensitive <- read.csv("data/primary-screen-replicate-collapsed-logfold-change.csv")
DEPMAPID_DT <- read.delim("data/data_clinical_sample.txt", skip = 4, header = 4)
## to translate the cell line ID or cell name into doubling time
DEPMAPID_DT <- DEPMAPID_DT[,c("SAMPLE_ID", "DEPMAPID","DOUBLING_TIME")]   ## with doubling time according to DEPMAPID
DEPMAPID_DT$LINEAGE <- sapply(DEPMAPID_DT[,1], function(X) sub(".*?_", "", X))  ## extract the lineage from the sample ID
## drug information
Drug_inform <- read.delim("data/repurposing_drugs_20200324.txt", skip = 9, header = 9)
Drug_inform_broadID <- read.delim("data/repurposing_samples_20200324.txt", skip = 9, header = 9)

## build the database manually whose source is from internet
Drug_inform_broadID_external <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
colnames(Drug_inform_broadID_external) <- colnames(Drug_inform_broadID)
Drug_inform_broadID_external[1,1:2] <- c("BRD-A03506276-001-01-5","XL888")
Drug_inform_broadID_external[2,1:2] <- c("BRD-A53952395-003-25-1","prilocaine")
Drug_inform_broadID_external[3,1:2] <- c("BRD-A55484088-050-02-5","BNTX")
Drug_inform_broadID_external[4,1:2] <- c("BRD-A64064900-001-02-3","sangivamycin")
Drug_inform_broadID_external[5,1:2] <- c("BRD-A67516570-001-02-8","lafutidine")
Drug_inform_broadID_external[6,1:2] <- c("BRD-K69280563-001-01-8","vinorelbine")

Drug_inform_external <- data.frame(NA,NA,NA,NA,NA,NA)
colnames(Drug_inform_external) <- colnames(Drug_inform)



## Machine Learning: pre-trained model
model_ExpDT_reg <- readRDS("PG_data/model_ExpDT_reg.rds")
model_ExpDT_gent2_reg <- readRDS("PG_data/model_ExpDT_gent2_reg.rds")
model_ExpDT_ctr_MC20161 <- readRDS("PG_data/model_ExpDT_ctr_MC20161.rds")
model_ExpDT_ctr_MC12399 <- readRDS("PG_data/model_ExpDT_ctr_MC12399.rds")
model_ExpDT_immuno_collection <- readRDS("PG_data/model_ExpDT_immuno_collection.rds")
model_ExpDT_targeteddrug_collection <- readRDS("PG_data/model_ExpDT_targeteddrug_collection.rds")


## TCGA
## TCGA - patient data
TCGA.PCA18.cbp_patient_RAW <- readRDS("data/TCGA.PCA18.cbp_patient_RAW.rds")
## TCGA - RNA expression data
## TCGA.PCA18.cbp_exp.rsem_RAW is a version whose gene name and ID issue have been resolved by the following code:
TCGA.PCA18.cbp_exp.rsem_RAW <- readRDS("data/TCGA.PCA18.cbp_exp.rsem_RAW.rds")
## TCGA - RNA expression data which is in a ready-to-be-used formate
TCGA.PCA18.cbp_exp.rsem <- readRDS("data/TCGA.PCA18.cbp_exp.rsem.rds")


## GENT2
## list of genes from expression data in GENT2: it is extracted from the SQL dump file
list_gene_P2 <- read.csv("data/list_gene_P2.csv")
## list of samples from expression data in GENT2: it is extracted from the SQL dump file
list_sample_P2 <- read.csv("data/list_sample_P2.csv")

## sample data
sql_survival <- read_csv("data/Subtype_Data.csv", col_names = F)
## label some key columns
colnames(sql_survival) <- c("Series_id","Sample_id","Primary_site","Disease","Subtype","Stage","Molecular_subgroup","Histological_type",
                            "OS_Month","OS_Status",
                            "DFS_Month","DFS_Status",
                            "DSS_Month","DSS_Status",
                            "PFS_Month","PFS_Status",
                            "TNM_Stage","X18","DukesStage","X20","X21","X22","X23","X24")


## GENT2 data with 50 genes
GeneSymbolP2_df_sub_df <- read.csv("data/GeneSymbolP2_df_sub_df.csv", row.names = 1, check.names = F)
## GENT2 data with KI-67
GeneSymbolP2_df_MKI67 <- readRDS("data/GeneSymbolP2_df_MKI67.rds")


## Therapy response: 
## Expression data with 4 cohorts (Hugo, Riaz, Van Allen, MGH) from Freeman et al. (2022) - before batch effect correction
## Freeman, S. S. et al. Combined tumor and immune signals from genomes or transcriptomes predict outcomes of checkpoint inhibition in melanoma. Cell Rep. Med. 3, 100500–100500 (2022).
CPB4_exp_pretr <- readxl::read_xlsx("data/Freeman_Supplementary_Table_4.xlsx", 
                                    sheet = (1+8))
colnames(CPB4_exp_pretr)[1]="Name"  ## fix the column name
## patient data downloaded from Freeman et al. (2022)
SKCM_NIR_patient_list_fullcohort <- readxl::read_xlsx("data/Freeman_Supplementary_Table_paper.xlsx", 
                                                      sheet = (1+12))
SKCM_NIR_patient_4cohort <- list(
  MGH = readxl::read_xlsx("data/Lee_Supplementary_Table_paper.xlsx", 
                          sheet = c(1+6), skip=1),
  Allen = readxl::read_xlsx("data/Lee_Supplementary_Table_paper.xlsx", 
                            sheet = c(1+7), skip=1),
  Riaz = readxl::read_xlsx("data/Lee_Supplementary_Table_paper.xlsx", 
                           sheet = c(1+8), skip=1),
  Hugo = readxl::read_xlsx("data/Lee_Supplementary_Table_paper.xlsx", 
                           sheet = c(1+9), skip=1)
)

## fixing the missing value due to merged cell from Excel, by auto filling the NA row from previous parent row
## pick one of the following:
current_table <- SKCM_NIR_patient_4cohort$MGH
current_table <- SKCM_NIR_patient_4cohort$Riaz
## silence and un-silence part of the code for each case
for(i in 1:nrow(current_table)){
  if(all(is.na(current_table[i,1:7]))){ ## search for NA row
    ## double check whether the patient id from current NA row is the same as the one from parent row (the one that is going to fill this NA row)
    ## for MGH
    if(regexpr(substr(current_table$Patient[i-1],1,(nchar(current_table$Patient[i-1])-1)), current_table$Sample_id[i])){
    ## for Riaz
    # if(current_table$Patient[i-1] == sub("_.*","",current_table$Sample_id[i])){
      current_table[i,1:7] <- current_table[(i-1),1:7]
    } else {stop("Attempting to fill NA row with incorrect previous row")}
  }
}
## in order to have a consistent data, the # of NA row in original table and # of duplicated row on patient id should be the same
## for MGH
table(is.na(SKCM_NIR_patient_4cohort$MGH$Patient))["TRUE"] == table(duplicated(current_table$Patient))["TRUE"]
SKCM_NIR_patient_4cohort$MGH <- current_table; rm(current_table)
## for Riaz
table(is.na(SKCM_NIR_patient_4cohort$Riaz$Patient))["TRUE"] == table(duplicated(current_table$Patient))["TRUE"]
SKCM_NIR_patient_4cohort$Riaz <- current_table; rm(current_table)
SKCM_NIR_patient_4cohort$Riaz$Sample_id <- gsub("_p","_P",SKCM_NIR_patient_4cohort$Riaz$Sample_id)  ## change all the "pre" and "post" into "Pre" and "Post": to be consistent with supplementary table



## Lee et al. (2021)
## Lee, J. S. et al. Synthetic lethality-mediated precision oncology via the tumor transcriptome. Cell 184, 2487-2502.e13 (2021).
## 11 of 13 cohorts (not including Kim and Snyder cohorts) for immunotherapy although not all of them are going to be analyzed
## expression data were transformed from matrix into data frame and gene names were annotated into row/column names based on the meta data
Immuno_13_opt <- readRDS("data/Lee_Immuno_13_opt.rds")
## 10 cohorts for targeted drug therapy although not all of them are going to be analyzed
## expression data were transformed from matrix into data frame and gene names were annotated into row/column names based on the meta data
Targeted_10 <- readRDS("data/Lee_Targeted_10.rds")

## CTR database
## expression data and clinical data
RNAseq_data <- readRDS("data/CTR_RNAseq_data.rds")
RNAseq_clinic <- readRDS("data/CTR_RNAseq_clinic.rds")
RNAseq_clinic_df <- readRDS("data/CTR_RNAseq_clinic_df.rds")  ## same as the last one but integrated into one data frame
## MC data including MC20161 and MC12399
MC_data <- readRDS("data/CTR_MC_data.rds")
MC_clinic_df <- readRDS("data/CTR_MC_clinic_df.rds")



##------- Doubling Time distribution across cancer tissues ----------
plot_table <- DEPMAPID_DT %>% filter(!is.na(DOUBLING_TIME), !is.na(PRIMARY_SITE), !(PRIMARY_SITE %in% c("Na")))
plot_table <- plot_table[,c("DOUBLING_TIME", "PRIMARY_SITE")]
plot_table <- plot_table %>% filter(DOUBLING_TIME <900)
plot_table$DOUBLING_TIME <- log10(plot_table$DOUBLING_TIME)
plot_table$PRIMARY_SITE[plot_table$PRIMARY_SITE == "Haematopoietic_And_Lymphoid_Tissue"] <- "Blood_And_Lymphoid"

cancer_order <- plot_table %>% group_by(PRIMARY_SITE) %>% summarise(median = median(DOUBLING_TIME)) %>%  ## calculate median according to labels
  arrange(median) %>% select(PRIMARY_SITE) ## sorted cancer types according to median 
cancer_label <- cancer_order$PRIMARY_SITE
names(cancer_label) <- cancer_label
cancer_label <- gsub("_"," ",cancer_label)

## distribution plot with violin plot
hold_plot <-
ggplot(plot_table, aes(PRIMARY_SITE, DOUBLING_TIME, fill = PRIMARY_SITE)) + 
  geom_violin(aes(PRIMARY_SITE, DOUBLING_TIME, fill = PRIMARY_SITE), width=3) + geom_boxplot(aes(PRIMARY_SITE, DOUBLING_TIME),width=0.1) + geom_point(shape = 1, size = 0.5, alpha = 0.8)+
  geom_signif(test = "wilcox.test" , map_signif_level = TRUE, comparisons = list(c("Soft_Tissue", "Bone")), family = "Times New Roman", textsize = 6,vjust = -0.1, tip_length=0.15, y_position=2.8) + ## signif. bar
  geom_hline(yintercept = median(plot_table$DOUBLING_TIME), linetype = 2, col = "red", alpha = 0.5) +
  scale_x_discrete(limits = cancer_order$PRIMARY_SITE, label = cancer_label) +  ## sorted cancer types according to median 
  # scale_x_discrete(label = temp) +  ## sorted cancer types according to median 
  scale_y_continuous(limits = c(1.25,3))+
  # geom_label_repel(min.segment.length = 0, max.overlaps = 5, family = "Times New Roman") +
  theme_classic()+
  theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5),
        # theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5,margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(size = 20, hjust = 0.5,margin=margin(t=0,b=10,l=0,r=0)),
        axis.text =  element_text(size = 20), axis.text.x = element_text(angle=45,hjust=1,size =20), #axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),
        legend.position = "none")  +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,1,0.5,0.5), "cm"))+ ## axis title 10 units from two sides; right plot margin=1cm
  labs(y= bquote('Doubling Time log' ['10'] ('hr')), x = NULL, 
       title = "Cancer Cell Line Doubling Time across Different Cancer Type")#,subtitle = "(Comparision Based on t-test)")

## This is Figure 1A ##
ggsave(filename = "plot/DT_lineage_violin.pdf",
       hold_plot,
       width = 15, height = 8)

###### Correlation Analysis - SNV ############

## create a list of mutation genes with DT, according to clinical sample
mutations[,"DT"] <- RAW_clinical[match(mutations$Broad_ID, RAW_clinical$DEPMAPID), "DOUBLING_TIME"]
mutations[,"Cell_Line"] <- RAW_clinical[match(mutations$Broad_ID, RAW_clinical$DEPMAPID), "SAMPLE_ID"]

## to list the genetic variant record details:
Variant_9 <- table(mutations$Variant_Classification)[c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site")]

## the recorded alternated genes
GeneList <- unique(mutations$Hugo_Symbol) ## with original order


## calculate mean DT according to each mutation genes (method according to each cell lines)

## reorganize cell line table
## template for t test of DT via loop for each cell line (CL)
CL <- RAW_clinical
CL$Mutation_Gene_Significance_pValue <- NA
CL <- CL[-which(is.na(CL$DEPMAPID)),] ## exclude the cell lines without available DEPMAPID
CL <- CL[-which(is.na(CL$DOUBLING_TIME)),] ## exclude the cell lines without available doubling time
## this Fn column is for storage of data to tell, whether the CL has the mutation (Strong/Weak), or just WT (i.e. whether the CL$DEPMAPID can be founded in the mutations table)
CL$FUNCTION1 <- NA 
## this column is for storage of data to tell, whether the mutated CL has the Strong_Mutation or Weak Mutation
CL$Mutation_Strength <- NA 
head(CL)

## reorganize mutation table
## sort the Mutations data set according to Hugo_symbol, so that the mutation data is list and catalogued according to each gene one by one
mutations_CL <- mutations[order(mutations$Hugo_Symbol),]
## omit the empty (NA) column in mutations_CL
mutations_CL <- mutations_CL[, -which(apply(mutations_CL,MARGIN = 2,FUN = function(X) all(is.na(X))))]
## interpret whether the mutation is strong (9 variants) or weak (others) mutation or not: Strong_Mutation marked as TRUE, Weak_Mutation marked as FALSE
mutations_CL$Mutation_Significance <- sapply(1:nrow(mutations_CL), FUN = function(X) mutations_CL$Variant_Classification[X] %in% names(Variant_9))

## storage of the CL scanning result for each gene[i] from CL$Mutation_Strength column
Gene_VS_CL <- data.frame(DEPMAPID = CL$DEPMAPID, DT = CL$DOUBLING_TIME, row.names = CL$SAMPLE_ID)
## store the p value for each cell line
# ## Comp_Mean_CL will be the storage of test result for each single mutated gene among all the cell line
# Comp_Mean_CL <- data.frame(c())
# head(Comp_Mean_CL)

## for each mutation gene [i] from GeneList (original order), classify each cell line among Strong_Mutation, Weak_Mutation and WT
for (i in 1:length(GeneList)) {

  ## for each single gene [i], build a sub data table, MutCollection, which contains all the strong&weak mutation in this specific gene [i]
  ## MutCollection is temporary and rewritten for each gene [i]
  MutCollection <- mutations_CL[mutations_CL$Hugo_Symbol == sort(GeneList)[i],]
  
  ## CHECK_POINT: if mutation in gene[i] doesn't exist among all the cell line at all, skip everything, marked it down, and jump for the next loop i
  if (any(sapply(1:nrow(MutCollection), FUN = function(X) MutCollection$Broad_ID[X] %in% CL$DEPMAPID))) {
    ## do nothing and continue
  } else {
    # Comp_Mean_CL[i,] <- "All CLs are WT"  ## no mutations are founded
    # rownames(Comp_Mean_CL)[i] <- sort(GeneList)[i]
    Gene_VS_CL[,i+2] <- NA
    colnames(Gene_VS_CL)[i+2] <- sort(GeneList)[i]
    next
  }
  
  ## FUNCTION1 in CL: all the cell lines has strong/weak mutation in Gene[i] marked as TRUE
  
  CL$FUNCTION1<-CL$DEPMAPID %in% MutCollection$Broad_ID
  CL$Mutation_Strength<-MutCollection[match(CL$DEPMAPID , MutCollection$Broad_ID),"Mutation_Significance"]
  ## recover the mutated CL marked as FALSE due to multiple mutation in single CL in single gene - i.e. for specific gene, there are multiple DEPMAPID appearing, including TRUE & FALSE in Mutation_Significance column.
  ## and FALSE is selected and marked to the CL due to the order of the MutCollection - the first one from the multiple observations.
  CL$Mutation_Strength[CL$DEPMAPID %in% MutCollection$Broad_ID[MutCollection$Mutation_Significance == TRUE]] <- TRUE  
  CL$Mutation_Strength[is.na(CL$Mutation_Strength)]<-"WT"
  
  ## store the scanning result
  Gene_VS_CL[,i+2] <- CL$Mutation_Strength
  colnames(Gene_VS_CL)[i+2] <- sort(GeneList)[i]
  
}

##------- t-test ---------
library(dplyr); library(tidyr)
Comp_Mean_CL_ttest <- data.frame(matrix(NA,1,6))
colnames(Comp_Mean_CL_ttest) <- c("p.value", "DT_Mutation_Mean", "DT_WT_Mean", "Genetic_Significance", "Mutated_Count", "WT_Count")

i=1
for(current_gene in colnames(Gene_VS_CL)[-1:-2]){
  current_working_table <- Gene_VS_CL[,c("DEPMAPID","DT",current_gene)]
  colnames(current_working_table)[3] <- "gene"
  ind <- which(current_working_table$gene == "TRUE")
  indd <- which(current_working_table$gene != "TRUE")
  Comp_Mean_CL_ttest[i,] <- NA
  
  if(length(ind)>=2 & length(indd)>=2) {
    result <- t.test(current_working_table$DT[ind], current_working_table$DT[indd])
    Comp_Mean_CL_ttest$p.value[i] <- result$p.value
    Comp_Mean_CL_ttest$DT_Mutation_Mean[i] <- mean(current_working_table$DT[ind])
    Comp_Mean_CL_ttest$DT_WT_Mean[i] <- mean(current_working_table$DT[indd])
    Comp_Mean_CL_ttest$Genetic_Significance[i] <- NA
    rm(result)
  } else if(length(ind)<2){
    Comp_Mean_CL_ttest[i,1:4] <- NA
    Comp_Mean_CL_ttest$Genetic_Significance[i] <- "solo value in Strong_Mutation"
    
  } else if(length(indd)<2){
    Comp_Mean_CL_ttest[i,1:4] <- NA
    Comp_Mean_CL_ttest$Genetic_Significance[i] <- "solo value or data missing in WT&Weak_Mutation"
    
  }
  
  Comp_Mean_CL_ttest$Mutated_Count[i] <- length(ind)
  Comp_Mean_CL_ttest$WT_Count[i] <- length(indd)
  rownames(Comp_Mean_CL_ttest)[i] <- current_gene
  rm(ind, indd, current_working_table)
  i=i+1
}

Comp_Mean_Valid_CL_ttest <- Comp_Mean_CL_ttest %>% filter(!(Genetic_Significance %in% c("solo value in Strong_Mutation", "solo value or data missing in WT&Weak_Mutation")))
Comp_Mean_Valid_CL_ttest$q_value <- p.adjust(Comp_Mean_Valid_CL_ttest$p.value, method = "BH")

##------- Wilcox test ---------
library(dplyr); library(tidyr)
Comp_Mean_CL_Wilcox <- data.frame(matrix(NA,1,6))
colnames(Comp_Mean_CL_Wilcox) <- c("p.value", "DT_Mutation_Mean", "DT_WT_Mean", "Genetic_Significance", "Mutated_Count", "WT_Count")

i=1
for(current_gene in colnames(Gene_VS_CL)[-1:-2]){
  current_working_table <- Gene_VS_CL[,c("DEPMAPID","DT",current_gene)]
  colnames(current_working_table)[3] <- "gene"
  ind <- which(current_working_table$gene == "TRUE")
  indd <- which(current_working_table$gene != "TRUE")
  Comp_Mean_CL_Wilcox[i,] <- NA
  
  if(length(ind)>=2 & length(indd)>=2) {
    result <- wilcox.test(current_working_table$DT[ind], current_working_table$DT[indd])
    Comp_Mean_CL_Wilcox$p.value[i] <- result$p.value
    Comp_Mean_CL_Wilcox$DT_Mutation_Mean[i] <- mean(current_working_table$DT[ind])
    Comp_Mean_CL_Wilcox$DT_WT_Mean[i] <- mean(current_working_table$DT[indd])
    Comp_Mean_CL_Wilcox$Genetic_Significance[i] <- NA
    rm(result)
  } else if(length(ind)<2){
    Comp_Mean_CL_Wilcox[i,1:4] <- NA
    Comp_Mean_CL_Wilcox$Genetic_Significance[i] <- "solo value in Strong_Mutation"
    
  } else if(length(indd)<2){
    Comp_Mean_CL_Wilcox[i,1:4] <- NA
    Comp_Mean_CL_Wilcox$Genetic_Significance[i] <- "solo value or data missing in WT&Weak_Mutation"
    
  }
  
  Comp_Mean_CL_Wilcox$Mutated_Count[i] <- length(ind)
  Comp_Mean_CL_Wilcox$WT_Count[i] <- length(indd)
  rownames(Comp_Mean_CL_Wilcox)[i] <- current_gene
  rm(ind, indd, current_working_table)
  i=i+1
}

Comp_Mean_Valid_CL_Wilcox <- Comp_Mean_CL_Wilcox %>% filter(!(Genetic_Significance %in% c("solo value in Strong_Mutation", "solo value or data missing in WT&Weak_Mutation")))
Comp_Mean_Valid_CL_Wilcox$q_value <- p.adjust(Comp_Mean_Valid_CL_Wilcox$p.value, method = "BH")
table(Comp_Mean_Valid_CL_Wilcox)

## a table contains both t.test and wilcox result
temp1 <- Comp_Mean_Valid_CL_ttest %>% rownames_to_column("Gene")
colnames(temp1)[colnames(temp1)=="p.value"] <- "p_value"
temp2 <- Comp_Mean_Valid_CL_Wilcox %>% rownames_to_column("Gene")
colnames(temp2)[colnames(temp2)=="p.value"] <- "p_value"
hold <- temp1[,colnames(temp1)!="Genetic_Significance"] %>% left_join(temp2[,c("Gene","p_value","q_value")], by="Gene", suffix=c("_t.test","_wilcox"))
hold <- hold %>% relocate(q_value_t.test, p_value_wilcox, q_value_wilcox, .after = "p_value_t.test")
SNV_result_ttest_wilcox <- hold; rm(hold, temp1, temp2)

## Cancer gene annotation
SNV_result_ttest_wilcox$CancerGene <- SNV_result_ttest_wilcox$Gene %in% CancerGene$Gene.Symbol
SNV_result_ttest_wilcox$CancerGene_H <- CancerGene$Hallmark[match(SNV_result_ttest_wilcox$Gene,CancerGene$Gene.Symbol)]

## This is Table S2.1 CCLE SNV correlation ##
SNV_result_ttest_wilcox

##------- plotting ------
## single gene plot
current_gene <- "ZNRF3"
plot_table <- Gene_VS_CL %>% select(DT, ZNRF3)
colnames(plot_table)[2] <- "group"
ind <- plot_table$group == "TRUE"
plot_table$group[ind] <- "Mutated"
plot_table$group[!ind] <- "WT"
plot_table$DT <- log10(plot_table$DT)
input_bquote <- bquote(.(current_gene) ~ ', t-test q-value' ~ .(sub("e.*","",signif(Comp_Mean_Valid_CL_ttest[current_gene,"q_value"],digits=4))%>%as.numeric()) 
                       %*% 10^.(sub(".*e","",signif(Comp_Mean_Valid_CL_ttest[current_gene,"q_value"],digits=4))%>%as.numeric()) )

hold_plot <-
  ggplot(plot_table, aes(group, DT)) + geom_violin(aes(fill = group)) + geom_boxplot(aes(fill = group),width=0.15, alpha = 0.8) + geom_point(shape = 1, alpha = 0.5, size =2)+
  geom_signif(test = "t.test" ,map_signif_level = TRUE, comparisons = list(c("Mutated", "WT")), textsize =8,vjust = 0.7,margin_top = 0.07) +
  theme_classic() + 
  theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5), 
        axis.text = element_text(size = 20), 
        legend.position = c(0.13,0.85)) +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = "Sample Group", y = bquote('Doubling Time log' ['10'] ('hr')), fill = "Group", title = input_bquote)

## This Figure 1B ##
ggsave(filename = "plot/Mutation_ZNRF3_violin.pdf",
       hold_plot,
       width = 6.5, height = 5)

## q-q plot: manual calculation
plot_data <- data.frame(Observed = -log10(sort(Comp_Mean_Valid_CL_ttest$p.value)), Expected = 1:nrow(Comp_Mean_Valid_CL_ttest))
plot_data$Expected <- -log10(plot_data$Expected/nrow(Comp_Mean_Valid_CL_ttest))

hold_plot<-
ggplot(plot_data, aes(Expected, Observed)) + geom_point(shape = 1,size=2,alpha=0.7) +
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5))+
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('-log'['10'] ('Expected p-value')), y = bquote('-log'['10'] ('Observed p-value')), 
       title = paste("Q-Q Plot of t-test p-values")) +
  geom_abline(intercept = 0, slope = 1,col='grey', linetype =2, linewidth = 1)

## This Figure 1C ##
ggsave(filename = "plot/Mutation_ttest_pvalue_qqplot_manual.pdf",
       hold_plot,
       width = 5, height = 5)


###### Correlation Analysis - CNV ############


##------- data manipulate -------
colnames(CNV)[-1] <- sub(" \\(.*", "", colnames(CNV[-1]))
colnames(CNV)[1] <- "DEPMAPID"

## add DT and lineage data
CNV$DOUBLING_TIME <- DEPMAPID_DT[match(CNV[,1],DEPMAPID_DT[,"DEPMAPID"]),"DOUBLING_TIME"]
CNV$LINEAGE <- DEPMAPID_DT[match(CNV[,1],DEPMAPID_DT[,"DEPMAPID"]),"PRIMARY_SITE"]
# which(colnames(CNV) == "DOUBLING_TIME" | colnames(CNV) == "LINEAGE")  ## DOUBLING_TIME and LINEAGE location
CNV <- CNV[,c(1, 25370, 25371, 2:25369)]  ## rearrange
CNV <- CNV[!is.na(CNV$DOUBLING_TIME),] ## remove all the rows without doubling time (NA)

##------- analysis -------
## store the result
result_pearson <- apply(CNV[,-1:-3], 2, function(X) cor.test(CNV[,2], X, method= "pearson"))
result_spearman <- apply(CNV[,-1:-3], 2, function(X) cor.test(CNV[,2], X, method= "spearman"))

## result output
summary_CNV <- data.frame(Gene = colnames(CNV)[-1:-3], 
                          qvalue_pearson = NA, qvalue_spearman = NA, pvalue_pearson = NA, pvalue_spearman = NA, estimate_pearson = NA, estimate_spearman = NA)
for (i in 1:length(result_pearson)) {
  summary_CNV[i,4] <- result_pearson[[i]]["p.value"]
  summary_CNV[i,6] <- result_pearson[[i]]["estimate"]
  summary_CNV[i,5] <- result_spearman[[i]]["p.value"]
  summary_CNV[i,7] <- result_spearman[[i]]["estimate"]
  
  if(summary_CNV[i,1] != colnames(CNV)[i+3]){warning(paste("check the data consistency of",colnames(CNV)[i+3]))}
}
rm(result_pearson, result_spearman)
## q value
summary_CNV[,2] <- p.adjust(summary_CNV[,4], method = "BH")
summary_CNV[,3] <- p.adjust(summary_CNV[,5], method = "BH")

## cancer gene annotation
summary_CNV$CancerGene <- summary_CNV$Gene %in% CancerGene$Gene.Symbol
summary_CNV <- left_join(summary_CNV, CancerGene[,c("Gene.Symbol","Hallmark")], by = c(Gene="Gene.Symbol"))
colnames(summary_CNV)[colnames(summary_CNV) == "Hallmark"] <- "CancerGene_H"
## eventually, we are more interested in Spearman correlation
summary_CNV <- summary_CNV[,c(-2,-4,-6)]

## This is Table S2.2 CCLE CNV correlation ##
summary_CNV

## CODE CHECK
table(summmary_cnv_2check$Gene == summary_CNV$Gene)
table(abs(summmary_cnv_2check$pvalue_spearman-summary_CNV$pvalue_spearman)<0.0000000000001)
table(abs(summmary_cnv_2check$pvalue_spearman-summary_CNV$pvalue_spearman)<0.00000000000000001)

##------- plotting --------
## single gene plot
current_gene <- "BRMS1L"; gene_name <- "BRMS1L"
plot_table <- CNV[,c("DEPMAPID", current_gene, "DOUBLING_TIME", "LINEAGE")] %>% 
  filter(!is.na(DOUBLING_TIME))
colnames(plot_table)[2] <- "gene"
colnames(plot_table)[colnames(plot_table)=="LINEAGE"] <- "PRIMARY_SITE"

## location of the label
ptext_loca <- axis_location(plot_table,3,2,"rightup"); ptext_loca[1] <- log10(ptext_loca[1])
## stastistical test result label
input_bquote <- paste("Spr. q-value = ", sub("e.*","",signif(summary_CNV$qvalue_spearman[summary_CNV$Gene==gene_name],digits=4)))

hold_plot <-
ggplot(plot_table, aes(log10(DOUBLING_TIME), gene)) + geom_point(color = "dodgerblue4",size=2,alpha=0.7) + #scale_fill_distiller(palette = "Blues")+
  geom_smooth(method = "lm", se=F, col='grey', linetype =2) +
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), legend.position = "none", plot.subtitle = element_text(size = 20, hjust = 0.5)) +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('Doubling Time log' ['10'] ('hr')), y = bquote('Copy Number log' ['2'] ('CN+1')), col = "Cancer Type",
       title = paste(gene_name,sep="")) +
  scale_x_continuous(limits = c(min(log10(plot_table$DOUBLING_TIME)),max(log10(plot_table$DOUBLING_TIME)*1.05))) +
  stat_regline_equation(label.x = ptext_loca[1]*0.83, label.y = ptext_loca[2]*1.05, size = 6,family = "Times New Roman")+ ## trendline equation
  stat_cor(aes(label = paste(after_stat(rr.label), sep = "~`,`~")), method = "spearman",
           label.x = ptext_loca[1]*0.84, label.y = ptext_loca[2], size = 6,family = "Times New Roman") + ## R^2
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.91, y=ptext_loca[2]*0.95, size = 6,
           label = input_bquote)+
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.91, y=ptext_loca[2]*0.9, size = 6,
           label = paste("Rho =",signif(summary_CNV$estimate_spearman[summary_CNV$Gene==gene_name],digits=4)))

## This Figure 1D ##
ggsave(filename = "plot/DT_CNV_SKAP1_dotplot.pdf",
       hold_plot,
       width = 7, height = 6)

## q-q plot: manual calculation
plot_data <- data.frame(Observed = -log10(sort(summary_CNV$pvalue_spearman)), Expected = 1:length(na.omit(summary_CNV$pvalue_spearman)))
plot_data$Expected <- -log10(plot_data$Expected/length(na.omit(summary_CNV$pvalue_spearman)))
hold_plot<-
ggplot(plot_data, aes(Expected, Observed)) + geom_point(shape = 1,size=2,alpha=0.7) +
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5))+
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('-log'['10'] ('Expected p-value')), y = bquote('-log'['10'] ('Observed p-value')),
       title = paste("Q-Q Plot of Spearman's p-values")) +
  geom_abline(intercept = 0, slope = 1,col='grey', linetype =2, linewidth = 1)

## This Figure 1E ##
ggsave(filename = "plot/CNV_Spearman_pvalue_qqplot_manual.pdf",
       hold_plot,
       width = 5, height = 5)


###### Correlation Analysis - RNA expression ############


rownames(RNA_expression) <- RNA_expression[,1]
colnames(RNA_expression)[1] <- "DEPMAPID"
## change with proper names
colnames(RNA_expression) <- sub(" \\(", "_", colnames(RNA_expression))
colnames(RNA_expression) <- sub("\\)", "", colnames(RNA_expression))
colnames(RNA_expression)[-1] <- sapply(2:ncol(RNA_expression), function(X) sub(" \\(.*", "", colnames(RNA_expression)[X]))

DataTable <- DT.generator(RNA_expression, DEPMAPID_DT, c("DEPMAPID","PRIMARY_SITE","DOUBLING_TIME")); current_DT_GR <- "DT"
colnames(DataTable)[2] <- "LINEAGE"  ## to be consistent with the old naming

##------- correlation analysis @ expression -------
## store the result
result_pearson <- apply(DataTable[,-1:-3], 2, function(X) cor.test(DataTable$DOUBLING_TIME, X, method= "pearson"))
result_spearman <- apply(DataTable[,-1:-3], 2, function(X) cor.test(DataTable$DOUBLING_TIME, X, method= "spearman"))

## result output
summary_DataTable <- data.frame(Gene_Transcript = sub("_.*","",colnames(DataTable)[-1:-3]), 
                                qvalue_pearson = NA, qvalue_spearman = NA, pvalue_pearson = NA, pvalue_spearman = NA, estimate_pearson = NA, estimate_spearman = NA,
                                Entrez_Gene_id = sub(".*_","",colnames(DataTable)[-1:-3]), Hugo_Entrez = colnames(DataTable)[-1:-3])

for (i in 1:length(result_pearson)) {
  summary_DataTable[i,4] <- result_pearson[[i]]["p.value"]
  summary_DataTable[i,6] <- result_pearson[[i]]["estimate"]
  summary_DataTable[i,5] <- result_spearman[[i]]["p.value"]
  summary_DataTable[i,7] <- result_spearman[[i]]["estimate"]
  
  if(summary_DataTable$Hugo_Entrez[i] != colnames(DataTable)[i+3]){warning(paste("check the data consistency of",colnames(DataTable)[i+3]))}
}; rm(result_pearson, result_spearman) ## job done

## q value
summary_DataTable[,2] <- p.adjust(summary_DataTable[,4], method = "BH")
summary_DataTable[,3] <- p.adjust(summary_DataTable[,5], method = "BH")
## order according to 
summary_DataTable <- summary_DataTable[order(summary_DataTable$qvalue_spearman),]  ## reorder the list according to q value
summary_DataTable$order_SPq <- match(summary_DataTable$pvalue_spearman,sort(summary_DataTable$pvalue_spearman, na.last = T))
## cancer gene annotation
summary_DataTable$CancerGene <- summary_DataTable$Gene_Transcript %in% CancerGene$Gene.Symbol
summary_DataTable <- left_join(summary_DataTable, CancerGene[,c("Gene.Symbol","Hallmark")], by = c(Gene_Transcript="Gene.Symbol"))
colnames(summary_DataTable)[colnames(summary_DataTable) == "Hallmark"] <- "CancerGene_H"

## This is Table S2.3 CCLE RNA expression correlation ##
## eventually, we are more interested in Spearman correlation
summary_RNA_expression <- summary_DataTable[,c(-2,-4,-6)]; rm(summary_DataTable)

## CODE CHECK
table(summary_RNA_expression$Gene_Transcript == summmary_express_2check$Gene_Transcript)
table(abs(summary_RNA_expression$pvalue_spearman-summmary_express_2check$pvalue_spearman)<0.00000001)
table(abs(summary_RNA_expression$pvalue_spearman-summmary_express_2check$pvalue_spearman)<0.0000000000000000001)

##------- plotting --------
## single gene plot
current_gene <- "NPM3_10360"; gene_name <- "NPM3"
plot_table <- RNA_expression[,c("DEPMAPID", current_gene)] %>% 
  left_join(DEPMAPID_DT[,c("DEPMAPID","DOUBLING_TIME","PRIMARY_SITE")]) %>%
  filter(!is.na(DOUBLING_TIME))

colnames(plot_table)[2] <- "gene"

ptext_loca <- axis_location(plot_table,3,2,"rightdown"); ptext_loca[1] <- log10(ptext_loca[1])

input_bquote <- bquote(paste('Spr. q-value ='~ .(sub("e.*","",signif(summary_RNA_expression$qvalue_spearman[summary_RNA_expression$Gene_Transcript==gene_name],digits=4))%>%as.numeric())
                             %*% 10^.(sub(".*e","",signif(summary_RNA_expression$qvalue_spearman[summary_RNA_expression$Gene_Transcript==gene_name],digits=4))%>%as.numeric())))


hold_plot <-
ggplot(plot_table, aes(log10(DOUBLING_TIME), gene)) + geom_point(color = "dodgerblue4",size=2,alpha=0.7) + #scale_fill_distiller(palette = "Blues")+
  geom_smooth(method = "lm", se=F, col='grey', linetype =2) +
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), legend.position = "none", plot.subtitle = element_text(size = 20, hjust = 0.5)) +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('Doubling Time log' ['10'] ('hr')), y = bquote('RNA Expression Level log' ['2'] ('TPM+1')), col = "Cancer Type",
       title = paste(gene_name)) +
  scale_x_continuous(limits = c(min(log10(plot_table$DOUBLING_TIME)),max(log10(plot_table$DOUBLING_TIME)*1.05))) +
  stat_regline_equation(label.x = ptext_loca[1]*0.82, label.y = ptext_loca[2]*5, size = 6,family = "Times New Roman")+ ## trendline equation
  stat_cor(aes(label = paste(after_stat(rr.label), sep = "~`,`~")), method = "spearman",
           label.x = ptext_loca[1]*0.83, label.y = ptext_loca[2]*4, size = 6,family = "Times New Roman") + ## R^2
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*3, size = 6,
           label = input_bquote)+  ## qvalue
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*2, size = 6,
           label = paste("Rho =",signif(summary_RNA_expression$estimate_spearman[summary_RNA_expression$Gene_Transcript==gene_name],digits=4)))  ## est.

## This Figure 1F ##
ggsave(filename = "plot/DT_Exp_NPM3_dotplot.pdf",
       hold_plot,
       width = 8, height = 6)

## q-q plot: manual calculation
plot_data <- data.frame(Observed = -log10(sort(summary_RNA_expression$pvalue_spearman)), Expected = 1:length(na.omit(summary_RNA_expression$pvalue_spearman)))
plot_data$Expected <- -log10(plot_data$Expected/length(na.omit(summary_RNA_expression$pvalue_spearman)))

hold_plot<-
ggplot(plot_data, aes(Expected, Observed)) + geom_point(shape = 1,size=2,alpha=0.7) +
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5))+
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('-log'['10'] ('Expected p-value')), y = bquote('-log'['10'] ('Observed p-value')), 
       title = paste("Q-Q Plot of Spearman's p-values")) +
  geom_abline(intercept = 0, slope = 1,col='grey', linetype =2, linewidth = 1)

## This Figure 1G ##
ggsave(filename = "plot/Exp_Spearman_pvalue_qqplot_manual.pdf",
       hold_plot,
       width = 5, height = 5)



##------- GSEA - file constructing and output reading --------
### constructing GSEA input file - preranked gene list
GSEA_summary_RNA_expression_estSpm <- PrankGSEA.est(summary_RNA_expression, 1,5,7, 
                                                    p_cutoff =  1, est_cutoff =  -2, "GSEA_summary_RNA_expression_estSpm")

## auto read GSEA output files
GSEAres_RNAexp <- list()  ## build empty working object
## set mother folder wd
folderDir <- c("PG_data/GSEA_result") ## set wd to mother folder
## set index for target sub-folder in mother folder - this will be looped in first loop (i)
ind_getfolder <- c(grep("est_Spm", list.files(folderDir))) ## select all folders marked with "est"
for (i in 1:length(ind_getfolder)) {
  ## set sub-list i
  GSEAres_RNAexp[[i]] <- list() ## create empty list @ i wrapped in GSEAres list
  ## set target sub-folder
  fileDir <- list.files(folderDir, full.names = T)[ind_getfolder[i]]  ## grab current folder wd
  ## set index for target files in target sub-folder
  ind_getfile <- grep(".tsv", list.files(fileDir))
  
  ## loop insider target sub-folder
  for (j in 1:length(ind_getfile)) {
    GSEAres_RNAexp[[i]][[j]] <- read.delim(list.files(fileDir, full.names = T)[ind_getfile[j]], check.names = F, skip = 1)[,c(-2,-3,-12)]
    if (as.logical(length(grep("GO",list.files(fileDir)[ind_getfile[j]])))) {  ## reorganize DF
      ## separate collection name from GS name
      GSEAres_RNAexp[[i]][[j]]$GS.Type <- sapply(1:nrow(GSEAres_RNAexp[[i]][[j]]), 
                                                 function(X) substr(GSEAres_RNAexp[[i]][[j]]$NAME[X],1,4))
      GSEAres_RNAexp[[i]][[j]]$NAME <- sapply(1:nrow(GSEAres_RNAexp[[i]][[j]]), 
                                              function(X) sub(".*?(\\_+)", "",GSEAres_RNAexp[[i]][[j]]$NAME[X]))  ## this is used to find the first pattern only
    }
    if (as.logical(length(grep("_H_",list.files(fileDir)[ind_getfile[j]])))) {  ## reorganize DF
      ## separate collection name from GS name
      GSEAres_RNAexp[[i]][[j]]$GS.Type <- sapply(1:nrow(GSEAres_RNAexp[[i]][[j]]), 
                                                 function(X) substr(GSEAres_RNAexp[[i]][[j]]$NAME[X],1,8))
      GSEAres_RNAexp[[i]][[j]]$NAME <- sapply(1:nrow(GSEAres_RNAexp[[i]][[j]]), 
                                              function(X) sub(".*?(\\_+)", "",GSEAres_RNAexp[[i]][[j]]$NAME[X]))  ## this is used to find the first pattern only
    }
  }
  names(GSEAres_RNAexp[[i]]) <- gsub(".tsv","",list.files(fileDir)[ind_getfile])  ## rename DF in the current sub-list i
}; names(GSEAres_RNAexp) <- sub("weight_RNAexpression_","",list.files(folderDir)[ind_getfolder]) ## rename each sub-list
rm(folderDir, ind_getfolder, ind_getfile, i, j)

## This is Table S4.1 GSEA on Hallmark GS and S4.2 GSEA on GO GS
GSEAres_RNAexp

##------- GSEA result analysis --------------
## make a proper name for GSs
GSEAname_hallmark <- data.frame(
  Upper_case = c("ADIPOGENESIS", "ALLOGRAFT_REJECTION", "ANDROGEN_RESPONSE", "ANGIOGENESIS", "APICAL_JUNCTION",
                 "APICAL_SURFACE", "APOPTOSIS", "BILE_ACID_METABOLISM", "CHOLESTEROL_HOMEOSTASIS", "COAGULATION",
                 "COMPLEMENT", "DNA_REPAIR", "E2F_TARGETS", "EPITHELIAL_MESENCHYMAL_TRANSITION", "ESTROGEN_RESPONSE_EARLY",
                 "ESTROGEN_RESPONSE_LATE", "FATTY_ACID_METABOLISM", "G2M_CHECKPOINT", "GLYCOLYSIS", "HEDGEHOG_SIGNALING",
                 "HEME_METABOLISM", "HYPOXIA", "IL2_STAT5_SIGNALING", "IL6_JAK_STAT3_SIGNALING", "INFLAMMATORY_RESPONSE",
                 "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "KRAS_SIGNALING_DN", "KRAS_SIGNALING_UP", "MITOTIC_SPINDLE",
                 "MTORC1_SIGNALING", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "MYOGENESIS", "NOTCH_SIGNALING",
                 "OXIDATIVE_PHOSPHORYLATION", "P53_PATHWAY", "PANCREAS_BETA_CELLS", "PEROXISOME", "PI3K_AKT_MTOR_SIGNALING",
                 "PROTEIN_SECRETION", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "SPERMATOGENESIS", "TGF_BETA_SIGNALING", "TNFA_SIGNALING_VIA_NFKB",
                 "UNFOLDED_PROTEIN_RESPONSE", "UV_RESPONSE_DN", "UV_RESPONSE_UP", "WNT_BETA_CATENIN_SIGNALING", "XENOBIOTIC_METABOLISM"),
  
  name = c("Adipogenesis", "Allograft rejection", "Androgen response", "Angiogenesis", "Apical junction",
           "Apical surface", "Apoptosis", "Bile acid metabolism", "Cholesterol homeostasis", "Coagulation",
           "Complement", "DNA repair", "E2F targets", "Epithelial mesenchymal transition", "Estrogen response early",
           "Estrogen response late", "Fatty acid metabolism", "G2M checkpoint", "Glycolysis", "Hedgehog signaling",
           "Heme metabolism", "Hypoxia", "IL2 STAT5 signaling", "IL6 JAK STAT3 signaling", "Inflammatory response",
           "Interferon alpha response", "Interferon gamma response", "KRAS signaling DN", "KRAS signaling UP", "Mitotic spindle",
           "MTORC1 signaling", "MYC targets V1", "MYC targets V2", "Myogenesis", "Notch signaling",
           "Oxidative phosphorylation", "P53 pathway", "Pancreas beta cells", "Peroxisome", "PI3K AKT MTOR signaling",
           "Protein secretion", "Reactive oxygen species pathway", "Spermatogenesis", "TGF beta signaling", "TNFA signaling via NFKB",
           "Unfolded protein response", "UV response DN", "UV response UP", "WNT beta catenin signaling", "Xenobiotic metabolism")
)

## for positive NES
plot_data <- GSEAres_RNAexp$est_Spm$estSpm_H_wgh_PosReport %>% filter(`FDR q-val` <0.1) %>% arrange(NES)
plot_data$NAME <- GSEAname_hallmark$name[match(plot_data$NAME, GSEAname_hallmark$Upper_case)]
## GS label
GS_order <- plot_data %>% arrange(-NES) %>% select(NAME)
GS_label <- GS_order$NAME
names(GS_label) <- GS_label
GS_label[1] <- "KRAS Signaling DN"
GS_label[-1] <- str_to_title(GS_label[-1])

plot1 <-
ggplot(plot_data, aes(x=NES, y=NAME, col = -log10(`FDR q-val`), size = SIZE)) +
  geom_point() + scale_y_discrete(limit = plot_data$NAME,label = GS_label) + 
  scale_color_gradient(low = "blue", high = "red",na.value = "red", name = bquote(-log [10] ('FDR q-value'))) + 
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5), 
                          legend.position = c(0.75,0.3)) +
  guides(size="none")+
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x="NES", y="Gene Set Name", legend =c("a","b"), size = "Size")

## for negative NES
plot_data <- GSEAres_RNAexp$est_Spm$estSpm_H_wgh_NegReport %>% filter(`FDR q-val` <0.1) %>% arrange(desc(`FDR q-val`), desc(NES))
plot_data$NAME <- GSEAname_hallmark$name[match(plot_data$NAME, GSEAname_hallmark$Upper_case)]
## GS label
GS_order <- plot_data %>% select(NAME)
GS_label <- GS_order$NAME
names(GS_label) <- GS_label
GS_label <- gsub("response","Response",GS_label)
GS_label <- gsub("targets","Targets",GS_label)
GS_label <- gsub("signaling","Signaling",GS_label)
GS_label <- gsub("checkpoint","Checkpoint",GS_label)
GS_label[-c(1,9,11,13,14,15,16)] <- str_to_title(GS_label[-c(1,9,11,13,14,15,16)])
plot2 <-
ggplot(plot_data, aes(x=NES, y=NAME, col = -log10(`FDR q-val`), size = SIZE)) +
  geom_point() + scale_y_discrete(limit = plot_data$NAME, label=GS_label) + 
  scale_color_gradient(low = "blue", high = "red",na.value = "red", name = bquote(-log [10] ('FDR q-value'))) + 
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5), 
                          legend.position = c(0.4,0.3)) +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x="NES", y="Gene Set Name", legend =c("a","b"), size = "Size")

## This is Figure 2 ##
ggarrange(plot2, plot1, labels = c("A", "B"))


##------- Drug Sensitivity Analysis (CCLE) ----------
## data manipulation and analysis
## data tailoring
## remove redundant part in colnames
colnames(drugsensitive) <- sapply(colnames(drugsensitive), function(X) gsub("\\::.*", "", X))
colnames(drugsensitive)[1] <- c("DEPMAPID")
drugsensitive <- drugsensitive[nchar(drugsensitive[,1]) != 21, ] ## cell lines marked with "FAILED_STR" are removed

DataTable <- DT.generator(drugsensitive, DEPMAPID_DT, c("DEPMAPID","LINEAGE","DOUBLING_TIME"), prof.name = "DOUBLING_TIME")

## statistical test
result_pearson <- apply(DataTable[,-1:-3], 2, function(X) cor.test(DataTable[,3], X, method= "pearson"))
result_spearman <- apply(DataTable[,-1:-3], 2, function(X) cor.test(DataTable[,3], X, method= "spearman"))

## result output
summary_DataTable <- data.frame(DRUG_ID = colnames(DataTable)[-1:-3], 
                                qvalue_pearson = NA, qvalue_spearman = NA, pvalue_pearson = NA, pvalue_spearman = NA, estimate_pearson = NA, estimate_spearman = NA)
## interpret the drug ID into common name
summary_DataTable <- add.inform.col(summary_DataTable, Drug_inform_broadID[,1:2],1,1)
colnames(summary_DataTable)[ncol(summary_DataTable)] <- "pert_iname"
## add drug information according to common name
summary_DataTable <- add.inform.col(summary_DataTable, Drug_inform,8,1)
## supplement the missing name
ind <- summary_DataTable$DRUG_ID %in% Drug_inform_broadID_external$broad_id
summary_DataTable$pert_iname[ind] <- Drug_inform_broadID_external[match(summary_DataTable$DRUG_ID[ind], Drug_inform_broadID_external$broad_id), "pert_iname"]

for (i in 1:length(result_pearson)) {
  summary_DataTable[i,4] <- result_pearson[[i]]["p.value"]
  summary_DataTable[i,6] <- result_pearson[[i]]["estimate"]
  summary_DataTable[i,5] <- result_spearman[[i]]["p.value"]
  summary_DataTable[i,7] <- result_spearman[[i]]["estimate"]
  
  if(summary_DataTable$DRUG_ID[i] != colnames(DataTable)[i+3]){warning(paste("check the data consistency of",colnames(DataTable)[i+3]))}
}; rm(result_pearson, result_spearman) ## job done

## q value
summary_DataTable[,2] <- p.adjust(summary_DataTable[,4], method = "BH")
summary_DataTable[,3] <- p.adjust(summary_DataTable[,5], method = "BH")

## correct the moa
summary_DataTable$moa[summary_DataTable$moa %in% "nuclear factor erythroid derived|like (NRF2) activator"] <- "nuclear factor (erythroid-derived 2)-like 2 (NRF2) activator"

## output
summary_drug <- summary_DataTable; rm(summary_DataTable)

## drug information 
## interpret the drug ID into common name
summary_drugsensitive <- add.inform.col(summary_drugsensitive, Drug_inform_broadID[,1:2],1,1)
colnames(summary_drugsensitive)[ncol(summary_drugsensitive)] <- "pert_iname"
## supplement the missing name
ind <- summary_drugsensitive$DRUG_ID %in% Drug_inform_broadID_external$broad_id
summary_drugsensitive$pert_iname[ind] <- Drug_inform_broadID_external[match(summary_drugsensitive$DRUG_ID[ind], Drug_inform_broadID_external$broad_id), "pert_iname"]
## add drug information according to common name
summary_drugsensitive <- add.inform.col(summary_drugsensitive, Drug_inform,6,1)



### Single drug plots
## pick one of the following: (the name, current_gene, is used in order to stay consistent with the previous section)
## VU0238429
current_gene <- summary_drug_sig$DRUG_ID[order(summary_drug_sig$qvalue_spearman)[1]]; gene_name <- summary_drug_sig[order(summary_drug_sig$qvalue_spearman)[1],c("pert_iname","moa")]
## captopril 
current_gene <- "BRD-K54529596-001-26-3"; gene_name <- summary_drug_sig[summary_drug_sig$DRUG_ID==current_gene,c("pert_iname","moa")]
## fanetizole
current_gene <- "BRD-K60798049-001-03-8"; gene_name <- summary_drug_sig[summary_drug_sig$DRUG_ID==current_gene,c("pert_iname","moa")]

plot_table <- drugsensitive[,c("DEPMAPID", current_gene)] %>% 
  left_join(DEPMAPID_DT[,c("DEPMAPID","DOUBLING_TIME","LINEAGE")]) %>%
  filter(!is.na(DOUBLING_TIME))

colnames(plot_table)[2] <- "gene"
plot_table <- plot_table %>% filter(!is.na(plot_table$gene))

ptext_loca <- axis_location(plot_table,3,2,"rightup"); ptext_loca[1] <- log10(ptext_loca[1])

input_bquote <- bquote(paste('Spr. q-value ='~ .(sub("e.*","",signif(summary_drug$qvalue_spearman[summary_drug$DRUG_ID==current_gene],digits=4)%>%scientific())%>%as.numeric())
                             %*% 10^.(sub(".*e","",signif(summary_drug$qvalue_spearman[summary_drug$DRUG_ID==current_gene],digits=4)%>%scientific())%>%as.numeric())))

## silence and un-silence the following code for specific plot
hold_plot <-
  ggplot(plot_table, aes(log10(DOUBLING_TIME), gene)) + geom_point(color = "dodgerblue4",size=2,alpha=0.7) + #scale_fill_distiller(palette = "Blues")+
  geom_smooth(method = "lm", se=F, col='grey', linetype =2) +
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), legend.position = "none", plot.subtitle = element_text(size = 18, hjust = 0.5)) +
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('Doubling Time log' ['10'] ('hr')), y = "Drug Effect", col = "Cancer Type",
       title = paste(gene_name[1]), subtitle = str_to_title(gene_name[2])) +
  scale_x_continuous(limits = c(min(log10(plot_table$DOUBLING_TIME)),max(log10(plot_table$DOUBLING_TIME)*1.05))) +
  
  ## VU0238429
  stat_regline_equation(label.x = ptext_loca[1]*0.82, label.y = ptext_loca[2]*4.7-12, size = 6,family = "Times New Roman")+ ## trendline equation
  stat_cor(aes(label = paste(after_stat(rr.label), sep = "~`,`~")), method = "spearman",
           label.x = ptext_loca[1]*0.83, label.y = ptext_loca[2]*4-12, size = 6,family = "Times New Roman") + ## R^2
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*3-12, size = 6,
           label = input_bquote)+  ## qvalue
  annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*2-12, size = 6,
           label = paste("Rho =",signif(summary_drug$estimate_spearman[summary_drug$DRUG_ID==current_gene],digits=4)))  ## est.

## Captopril
# labs(title = str_to_title(gene_name[1]), subtitle = str_to_title(gene_name[2])) +
# stat_regline_equation(label.x = ptext_loca[1]*0.82, label.y = ptext_loca[2]*4.2-6.5, size = 6,family = "Times New Roman")+ ## trendline equation
# stat_cor(aes(label = paste(after_stat(rr.label), sep = "~`,`~")), method = "spearman",
#          label.x = ptext_loca[1]*0.83, label.y = ptext_loca[2]*3.9-6.5, size = 6,family = "Times New Roman") + ## R^2
# annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*3.6-6.5, size = 6,
#          label = input_bquote)+  ## qvalue
# annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*3.3-6.5, size = 6,
#          label = paste("Rho =",signif(summary_drug$estimate_spearman[summary_drug$DRUG_ID==current_gene],digits=4)))  ## est.

## Fanetizole
# labs(title = str_to_title(gene_name[1]), subtitle = str_to_title(gene_name[2])) +
# stat_regline_equation(label.x = ptext_loca[1]*0.82, label.y = ptext_loca[2]*4.45-4, size = 6,family = "Times New Roman")+ ## trendline equation
# stat_cor(aes(label = paste(after_stat(rr.label), sep = "~`,`~")), method = "spearman",
#          label.x = ptext_loca[1]*0.83, label.y = ptext_loca[2]*4.2-4, size = 6,family = "Times New Roman") + ## R^2
# annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*3.9-4, size = 6,
#          label = input_bquote)+  ## qvalue
# annotate("text",family = "Times New Roman", x=ptext_loca[1]*0.9, y=ptext_loca[2]*3.6-4, size = 6,
#          label = paste("Rho =",signif(summary_drug$estimate_spearman[summary_drug$DRUG_ID==current_gene],digits=4)))  ## est.


## plot output - Figure 9A, C, D ##
## pick one of the following
ggsave(filename = "plot/DS_VU238429.pdf", ## This is Figure 9A
       # ggsave(filename = "plot/DS_Captopril.pdf", ## This is Figure 9C
       # ggsave(filename = "plot/DS_Fanetizole.pdf", ## This is Figure 9D
       hold_plot,
       width = 6, height = 5)



## QQ plot
## manual calculation
plot_data <- data.frame(Observed = -log10(sort(summary_drug$pvalue_spearman)), Expected = 1:length(na.omit(summary_drug$pvalue_spearman)))
plot_data$Expected <- -log10(plot_data$Expected/length(na.omit(summary_drug$pvalue_spearman)))

hold_plot<-
  ggplot(plot_data, aes(Expected, Observed)) + geom_point(shape = 1,size=2,alpha=0.7) +
  theme_classic() + theme(text = element_text(size = 20, family = "Times New Roman"), plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5))+
  theme(plot.title = element_text(margin=margin(t=10,b=10,l=0,r=0)), plot.subtitle = element_text(margin=margin(t=0,b=10,l=0,r=0)), ## title and subtitle 10 units from up and down
        axis.title.y = element_text(margin = margin(l=10,r=10,t=0,b=0)), axis.title.x = element_text(margin = margin(t=10,b=10,l=0,r=0)),plot.margin = unit(c(0,0.7,0.2,0.2), "cm"))+ ## axis title 10 units from two sides; right plot margin=0.7cm
  labs(x = bquote('-log'['10'] ('Expected p-value')), y = bquote('-log'['10'] ('Observed p-value')), 
       title = paste("Q-Q Plot of Spearman's p-values")) +
  geom_abline(intercept = 0, slope = 1,col='grey', linetype =2, linewidth = 1)

## This is Figure 9B ##
ggsave(filename = "plot/DS_Spearman_pvalue_qqplot_manual.pdf",
       hold_plot,
       width = 5, height = 5)






