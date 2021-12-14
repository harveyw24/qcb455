#################################### Setup #####################################

library(affy)
library(dplyr)
library(data.table)

# Main directory
MAIN_DIR = 'D:/Data'

# Datasets of CEL files
DATASETS = list('GSE11121',
                'GSE18864', 
                'GSE20711', 
                'GSE23593', 
                'GSE27120', 
                'GSE32646',
                'GSE36771', 
                'GSE42568', 
                'GSE50948', 
                'GSE5460', 
                'GSE11001', 
                'GSE87007', 
                'GSE88770', 
                'GSE7390', 
                'GSE78958', 
                'GSE45255', 
                'GSE61304', 
                'GSE63471', 
                'GSE21653', 
                'GSE26639', 
                'GSE17907',
                'GSE10810', 
                'GSE25066', 
                'GSE47109', 
                'GSE95700', 
                'GSE5327', 
                'GSE48390', 
                'GSE58984',
                'GSE103091', 
                'GSE45827', 
                'GSE65194', 
                'GSE1456A',
                'GSE102484')
# 'GSE1456B' Removed

#################################### Part 0 ####################################
############################# Load and Unzip Data ##############################
##################################### (HW) #####################################

# Unzip all CEL.gz files in data directory
library(R.utils)
for (dataset in DATASETS) {
    # Set working directory in each dataset subfolder
    celpath = sprintf('%s/%s', MAIN_DIR, dataset)
    setwd(celpath)

    files <- list.files(celpath)

    # unzip all .gz files in subfolder
    for (file in files) {
        gunzip(file, remove=FALSE)
    }
}

# Check to make sure no .gz files are left
for (dataset in DATASETS) {
    celpath = sprintf('%s/%s', MAIN_DIR, dataset)
    setwd(celpath)

    files <- list.files(celpath)

    # Print names of any files with non .gz file type
    for (file in files) {
        if (substring(file, nchar(file)-1, nchar(file)) == 'gz') {
            print(paste(celpath, file))
        }
    }
}



#################################### Part 1 ####################################
################################ Apply ReadAffy ################################
########################## Copied from provided code ###########################


# Iterate through all dataset subfolders
for (dataset in DATASETS){
    celpath = sprintf('%s/%s', MAIN_DIR, dataset)
    setwd(celpath)
    print(celpath)
    
    # Change the cel data to Affy to use in RMA (ReadAffy)
    cancer_affy <- affy::ReadAffy(celfile.path = celpath)
    
    # Save the file
    save(cancer_affy, file=sprintf('%s/%s/%s_Affy.RData', MAIN_DIR, dataset, dataset))
}



################################### Part 2 #####################################
################################# Apply RMA ####################################
########################## Copied from provided code ###########################

NORMALIZING = 'YES' # If normalizing for the first time select YES, otherwise select NO
setwd(MAIN_DIR) # change working directory 

for (dataset in DATASETS){
    print(dataset)
    
    # Loading the Affy data from Part 1 
    env <- new.env()
    nm <- load(sprintf('%s/%s/%s_Affy.RData', MAIN_DIR, dataset, dataset), env)[1]
    cancer_data <- env[[nm]]
    
    # Normalize with RMA
    if (NORMALIZING == 'YES'){
        
        # Create working directory for RMA
        dir.create(sprintf('%s/%s/rma', MAIN_DIR, dataset))
        
        # Perform RMA and compute expression measures
        cancer_data.rma <- rma(cancer_data)
        
        # Save RMA expression data
        save(cancer_data.rma , file= sprintf('%s/%s/rma/%s.rma.RData', MAIN_DIR, dataset, dataset))
        cancer_data.rma
        
    } else{
        
        # Loading RMA expression dataset from subfolder
        env <- new.env()
        nm <- load(sprintf('%s/%s/rma/%s.rma.RData', MAIN_DIR, dataset, dataset), env)[1]
        cancer_data.rma <- env[[nm]]
        
    }
    
    # Save RMA expression data in matrix form
    cancer_data.matrix <- affy::exprs(cancer_data.rma)
    head(cancer_data.matrix)
    
    # Save expression matrix as RData (.matrix.RData)
    save(cancer_data.matrix, file=sprintf('%s/%s/rma/%s.matrix.RData', MAIN_DIR, dataset, dataset))
    
    # Save expression matrix as csv (.matrix.csv)
    write.csv(cancer_data.matrix, file= sprintf('%s/%s/rma/%s.matrix.csv', MAIN_DIR, dataset, dataset), row.names=TRUE)
}



#################################### Part 3 ####################################
############################ Merge expression sets #############################
####################### Adapted from provided code (HW) ########################

source("util_merge.R")
source("util_misc.R")
source("util_plot.R")

# Merge cel Affy data to combined expression set
esets = list()
for(dataset in DATASETS) {
    env <- new.env()
    filename = sprintf('%s/%s/rma/%s.rma.RData', MAIN_DIR, dataset, dataset)
    print(sprintf('Loading file %s', filename))
    nm <- load(filename, env)[1]
    cancer_data <- env[[nm]]
    pData(cancer_data)['sample'] = dataset
    esets <- append(esets, cancer_data)
}

# eset_COMBAT = inSilicoMerging::merge(esets, method="COMBAT")
eset_COMBAT = merge(esets, method="COMBAT")

# Save merged expression set data
dir.create(sprintf('%s/analysed_datasets', MAIN_DIR))
dir.create(sprintf('%s/analysed_datasets/merged', MAIN_DIR))
save(eset_COMBAT, file=sprintf('%s/analysed_datasets/merged/eset_merged_COMBAT_rma_new.RData', MAIN_DIR))

cancer_data.matrix <- affy::exprs(eset_COMBAT)
head(cancer_data.matrix)

# Generate list of duplicated samples
len = ncol(cancer_data.matrix)
duplicated_names = list()
duplicated_indices = list()
sample_names = colnames(cancer_data.matrix)
j = 0
for(i in 1:len){
    # Retrieve the gene name (GSM...)
    item = unlist(strsplit(unlist(strsplit(colnames(cancer_data.matrix)[i], split='.', fixed=TRUE))[1], split='_', fixed=TRUE))[1]
    if(!item %in% sample_names){
        sample_names[i] = item
    }
    else {
        duplicated_names[j] = sample_names[i]
        duplicated_indices[j] = i
        j = j + 1
    }
}

unlist(duplicated_names)
duplicated_indices = unlist(duplicated_indices)
duplicated_indices

# Remove duplicate samples from gene expression matrix
cancer_data_unique.matrix = cancer_data.matrix[, -duplicated_indices]
len2 = ncol(cancer_data_unique.matrix)
sample_names = colnames(cancer_data_unique.matrix)

for(i in 1:len2){
    item = unlist(strsplit(unlist(strsplit(colnames(cancer_data_unique.matrix)[i], split='.', fixed=TRUE))[1], split='_', fixed=TRUE))[1]
    sample_names[i] = item  
}

colnames(cancer_data_unique.matrix)= sample_names

# Save expression matrix with unique genes as csv (.matrix.csv)
write.csv(cancer_data_unique.matrix, file= sprintf('%s/analysed_datasets/merged/merged_COMBAT_rma_new.matrix.csv',MAIN_DIR), row.names=TRUE)



######################### Map probes to gene symbols ###########################
##################################### (HW) #####################################

# Retrieve list of gene symbols associated with probe ids
BiocManager::install("annotate")
BiocManager::install("hgu133a.db")
BiocManager::install("hgu133b.db")
BiocManager::install("hgu133plus2.db")

library("annotate")
library("hgu133a.db")
library("hgu133b.db")
library("hgu133plus2.db")

set_a <- hgu133aSYMBOL
# Get probe ID mapped to gene symbol
mapped_probes_a <- mappedkeys(set_a)
# Convert to list
probe_list_a <- as.list(set_a[mapped_probes_a])

set_b <- hgu133bSYMBOL
# Get probe ID mapped to gene symbol
mapped_probes_b <- mappedkeys(set_b)
# Convert to list
probe_list_b <- as.list(set_b[mapped_probes_b])

set_plus <- hgu133plus2SYMBOL
# Get probe ID mapped to gene symbol
mapped_probes_plus <- mappedkeys(set_plus)
# Convert to list
probe_list_plus <- as.list(set_plus[mapped_probes_plus])

probe_list_all <- c(probe_list_a, probe_list_b, probe_list_plus)

probes <- as.vector(names(probe_list_all))
genes <- as.vector(unlist(probe_list_all))

probe2gene <- data.frame(probes)
probe2gene$genes = genes

write.csv(probe2gene, 'probe2gene.csv', row.names = FALSE)



############################## Annotation Cleanup ##############################
##################################### (HW) #####################################

# Retrieve annotations for each sample
BiocManager::install("GEOquery")

library("GEOquery")
library(dplyr)

dataset = DATASETS[[17]]

gset <- getGEO(dataset, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


gset_chars <- Biobase::pData(gset) %>% as_tibble()

gset_chars = gset_chars %>% select("geo_accession", contains(":ch1"))
gset_chars = gset_chars %>% rename(sample_id = geo_accession)

colnames(gset_chars) = gsub(":ch1", "", colnames(gset_chars))
View(gset_chars)

write.csv(gset_chars,
          file=sprintf('%s_annotation.csv', dataset),
          row.names = FALSE)



########################### DATASET-SPECIFIC CHANGES ###########################
##################################### (HW) #####################################

# GSE45255
gset_chars <- gset_chars %>%
    select(sample_id,
           death_event=`dfs event (defined as any type of recurrence or death from breast cancer)`,
           death_event_time_yrs=`dfs time`,
           dmfs_event=`dmfs event (defined as distant metastasis or death from breast cancer)`,
           dmfs_event_time_yrs=`dmfs time`,
           dss_event=`dss event (defined as death from breast cancer)`,
           dss_event_time=`DSS time`,
           er=`er status`,
           her2=`her2 status`,
           grade=`histological grade`,
           tissue=histology,
           ln=`ln status`,
           age=`patient age`,
           pgr=`pgr status`,
           size_in_mm=`size (mm)`)

gset_chars$grade <- gsub("G", "", gset_chars$grade)

gset_chars$er[grepl('ER+', gset_chars$er, fixed=TRUE)] <- '1'
gset_chars$er[grepl('ER-', gset_chars$er, fixed=TRUE)] <- '0'

gset_chars$her2[grepl('He+', gset_chars$her2, fixed=TRUE)] <- '1'
gset_chars$her2[grepl('He-', gset_chars$her2, fixed=TRUE)] <- '0'

gset_chars$ln[grepl('LN+', gset_chars$ln, fixed=TRUE)] <- '1'
gset_chars$ln[grepl('LN-', gset_chars$ln, fixed=TRUE)] <- '0'

gset_chars$pgr[grepl('PgR+', gset_chars$pgr, fixed=TRUE)] <- '1'
gset_chars$pgr[grepl('PgR-', gset_chars$pgr, fixed=TRUE)] <- '0'



gset_chars





# GSE78958
gset_chars <- gset_chars %>% select(-c(bmi, `patient ethnicity`))
gset_chars <- gset_chars %>% rename(grade=`tumor grade`, stage=`tumor stage`,
                                    subtype=`tumor subtype (via breastprs)`)

gset_chars$grade[grepl('grade 1', gset_chars$grade, fixed=TRUE)] <- '1'
gset_chars$grade[grepl('grade 2', gset_chars$grade, fixed=TRUE)] <- '2'
gset_chars$grade[grepl('grade 3', gset_chars$grade, fixed=TRUE)] <- '3'

gset_chars$stage <- gsub("Stage ", "", gset_chars$stage)



gset_chars


# GSE7390
gset_chars <- gset_chars %>% select(sample_id,
                      age,
                      e.dmfs, e.os, e.rfs, e.tdm, er,
                      grade,
                      size_in_mm=size,
                      t.dmfs, t.os, t.rfs, t.tdm)

gset_chars$size_in_mm <- as.double(gset_chars$size_in_mm) * 10
gset_chars$t.dmfs <- as.double(gset_chars$t.dmfs) / 365
gset_chars$t.os <- as.double(gset_chars$t.os) / 365
gset_chars$t.rfs <- as.double(gset_chars$t.rfs) / 365
gset_chars$t.tdm <- as.double(gset_chars$t.tdm) / 365

gset_chars



# GSE88770
gset_chars <- gset_chars %>% select(sample_id,
                      death_event=death,
                      e.rfs=drfs_event,
                      t.rfs=drfs_or_last_contact_years,
                      er,
                      gender,
                      grade,
                      her2,
                      ki67,
                      death_event_time_yrs=os_or_last_contact_years,
                      pgr,
                      tissue)

gset_chars[gset_chars == 'Yes'] = '1'
gset_chars[gset_chars == 'Positive'] = '1'
gset_chars[gset_chars == 'No'] = '0'
gset_chars[gset_chars == 'Negative'] = '0'

gset_chars

# GSE87007
gset_chars <- gset_chars %>% select(-c(cellularity))
gset_chars <- gset_chars %>% rename(subtype=`molecular subtype`)

gset_chars

# GSE11001
gset_chars <- gset_chars %>% select(-`tumor cells`)
gset_chars <- gset_chars %>% rename(er=ER, pr=PR, her2=HER2, node=`node stage`,
                                    stage=`tumor stage`)
gset_chars$er[grepl('positive', gset_chars$er, fixed=TRUE)] <- '1'
gset_chars$er[grepl('negative', gset_chars$er, fixed=TRUE)] <- '0'

gset_chars$her2[grepl('positive', gset_chars$her2, fixed=TRUE)] <- '1'
gset_chars$her2[grepl('negative', gset_chars$her2, fixed=TRUE)] <- '0'

gset_chars$pr[grepl('positive', gset_chars$pr, fixed=TRUE)] <- '1'
gset_chars$pr[grepl('negative', gset_chars$pr, fixed=TRUE)] <- '0'

gset_chars



# GSE5460
gset_chars <- gset_chars %>% select(sample_id, grade=`B-R grade`, er=ER,
                      node=`node status`, size_in_mm=`tumor size`)

gset_chars$grade <- as.numeric(as.roman(gset_chars$grade))

gset_chars$er[grepl('pos', gset_chars$er, fixed=TRUE)] <- '1'
gset_chars$er[grepl('neg', gset_chars$er, fixed=TRUE)] <- '0'

gset_chars$node[grepl('pos', gset_chars$node, fixed=TRUE)] <- '1'
gset_chars$node[grepl('neg', gset_chars$node, fixed=TRUE)] <- '0'

gset_chars$size_in_mm <- as.double(gset_chars$size_in_mm) * 10

gset_chars



# GSE50948
gset_chars <- gset_chars %>% select(sample_id, age, bgus.ct, er.ct, er, her2.ct,
                      her2, inflamed_brca=inflammatory.brca, pr.ct, pr,
                      size1=`invasive_tumor_area_size1 [mm]`,
                      size2=`invasive_tumor_area_size2 [mm]`)

gset_chars$er[grepl('ER+', gset_chars$er, fixed=TRUE)] <- '1'
gset_chars$er[grepl('ER-', gset_chars$er, fixed=TRUE)] <- '0'

gset_chars$her2[grepl('HER2+', gset_chars$her2, fixed=TRUE)] <- '1'
gset_chars$her2[grepl('HER2-', gset_chars$her2, fixed=TRUE)] <- '0'

gset_chars$pr[grepl('PR+', gset_chars$pr, fixed=TRUE)] <- '1'
gset_chars$pr[grepl('PR-', gset_chars$pr, fixed=TRUE)] <- '0'

gset_chars$inflamed_brca[grepl('yes', gset_chars$inflamed_brca, fixed=TRUE)] <- '1'
gset_chars$inflamed_brca[grepl('no', gset_chars$inflamed_brca, fixed=TRUE)] <- '0'

gset_chars <- gset_chars %>% mutate(size_in_mm=as.double(size1)+as.double(size2))

gset_chars <- gset_chars %>% select(-c(size1, size2))

gset_chars

# GSE42568
gset_chars <- gset_chars %>% rename(er=er_status, 'death_event'=`overall survival event`,
                      'death_event_time_yrs'=`overall survival time_days`,
                      'relapse_event' = `relapse free survival event`,
                      'relapse_time_yrs'=`relapse free survival time_days`,
                      size_in_mm=size)

gset_chars$size_in_mm <- as.double(gset_chars$size_in_mm) * 10
gset_chars$death_event_time_yrs <- as.double(gset_chars$death_event_time_yrs) / 365
gset_chars$relapse_time_yrs <- as.double(gset_chars$relapse_time_yrs) / 365

gset_chars




# GSE36771
gset_chars <- gset_chars %>% select(-2)
gset_chars <- gset_chars %>% rename(er=`estrogen receptor status`,
                                    pr=`progesterone receptor status`,
                                    grade=`tumour histological grade`)

gset_chars$er[gset_chars$er == 'ER+'] <- '1'
gset_chars$er[gset_chars$er == 'ER-'] <- '0'

gset_chars$pr[gset_chars$pr == 'PgR+'] <- '1'
gset_chars$pr[gset_chars$pr == 'PgR-'] <- '0'

gset_chars$`lymph node metastasis`[gset_chars$`lymph node metastasis` == 'LN+'] <- '1'
gset_chars$`lymph node metastasis`[gset_chars$`lymph node metastasis` == 'LN-'] <- '0'

gset_chars$grade= substring(gset_chars$grade, 2)

gset_chars



# GSE32646
gset_chars <- gset_chars %>% select(-c(3,8,9))
gset_chars <- gset_chars %>% rename(stage=`clinical t stage`,
                      er=`er status ihc`, her2=`her2 status fish`,
                      grade=`histological grade`, pr=`pr status ihc`)
gset_chars$er[grepl('positive', gset_chars$er)] <- '1'
gset_chars$er[grepl('negative', gset_chars$er)] <- '0'

gset_chars$pr[grepl('positive', gset_chars$pr)] <- '1'
gset_chars$pr[grepl('negative', gset_chars$pr)] <- '0'

gset_chars$her2[grepl('positive', gset_chars$her2)] <- '1'
gset_chars$her2[grepl('negative', gset_chars$her2)] <- '0'

gset_chars

# GSE27120
gset_chars <- gset_chars %>% select(-4)
gset_chars <- gset_chars %>% rename(tissue=`cell/tissue type`, er=`er_status (1= positive, 0= negative)`,
                      her2=`her2_status (1= positive, 0= negative)`,
                      pgr=`pgr_status (1= positive, 0= negative)`, size_in_mm=`size_tumor (cm)`,
                      age=`age (year)`)
gset_chars <- gset_chars %>% select(-`ki67 (%)`)
gset_chars$size_in_mm <- as.double(gset_chars$size_in_mm) * 10

gset_chars

# GSE23593
gset_chars <- gset_chars %>% rename(er='er status', pr='pr status', grade='tumor grade')
gset_chars <- gset_chars %>% select(!c('disease state', patient))

gset_chars$er[grepl('Positive', gset_chars$er)] <- '1'
gset_chars$er[grepl('Negative', gset_chars$er)] <- '0'

gset_chars$pr[grepl('Positive', gset_chars$pr)] <- '1'
gset_chars$pr[grepl('Negative', gset_chars$pr)] <- '0'

# GSE20711
gset_chars <- gset_chars %>% select(!c('age (bin)', 'methylation barcode',
                                       'size (bin)', 'node', 'quality control'))
gset_chars <- gset_chars %>% rename(er='er status', her2='her2 status')
gset_chars$t.os <- gsub(" y", "", gset_chars$t.os)

# GSE18864
gset_chars$grade <- as.numeric(as.roman(gset_chars$grade))

gset_chars$er <- sapply(strsplit(gset_chars$`er/pr/her2 status`, '/'), '[', 1)
gset_chars$pr <- sapply(strsplit(gset_chars$`er/pr/her2 status`, '/'), '[', 2)
gset_chars$her2 <- sapply(strsplit(gset_chars$`er/pr/her2 status`, '/'), '[', 3)

gset_chars <- gset_chars %>% select(!'er/pr/her2 status')

# GSE11121
gset_chars$t.dmfs<- as.double(gset_chars$t.dmfs) / 12
gset_chars$size_in_cm <- as.double(gset_chars$size_in_cm) * 10
gset_chars <- gset_chars %>% rename(size_in_mm = size_in_cm)

