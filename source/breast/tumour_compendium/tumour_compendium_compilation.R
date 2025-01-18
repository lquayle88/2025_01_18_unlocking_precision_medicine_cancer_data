## Combine studies to compose tumour compendium
## Lewis Quayle
## 24/02/2021


## load required packages

library("jetset")
library("sva")


## set project directory hierarchy

# home directory

home.dir <- "/Users/Lewis/"

# parent directories

rna.seq.dir <- paste0(home.dir, "Documents/Bioinformatics/rna_seq_analysis/")
parent.dir <- paste0(rna.seq.dir, "public_clinical_data_analysis/")

# data input directory

input.dir <- paste0(parent.dir, "data/tumour_compendium/")

# output directory

output.dir <- paste0(parent.dir, "data/tumour_compendium/tumour_compendium/")


## read-in BrCa core quiescence signature

signature <- read.csv(file = paste0(rna.seq.dir, "BrCa_core_quiescence/results/tabular/BrCa_core_quiescence_signature.csv"),
                      header = TRUE,
                      stringsAsFactors = FALSE)


## read-in sample information for all studies and combine into data frame

# list of files

files <- list.files(path = input.dir,
                    recursive = TRUE,
                    pattern = "sample_info",
                    full.names = TRUE)

# remove the tumour compendium file itself from this list

files <- files[-c(grep(pattern = "tumour_compendium_sample_info", x = files))]

# read-in and combine files into a data frame

combined_sample_info <- data.frame()
batch <- character()

for (i in seq_along(files)) {
  
  # read file
  
  sample.info <- read.csv(file = files[i],
                          header = TRUE,
                          stringsAsFactors = FALSE)
  
  # add to combined data frame
  
  combined_sample_info <- rbind(combined_sample_info, sample.info)
  
  # create batch labels vector for metadata used in adjusting batch effect
  
  study.lab <- sub(pattern = "_.*", replacement = "", x = basename(files[i]))
  study.batch.labels <- rep(study.lab, times = nrow(sample.info))
  batch <- c(batch, study.batch.labels)
  
}

# add study identifier to data frame

combined_sample_info <- data.frame("STUDY" = batch, combined_sample_info)


## read-in gene-level expression data for all studies

# create a gene vector for sig. genes in GPL96 and GPL570

GPL96 <- jmap(chip = "hgu133a", symbol = signature$Symbol)
GPL96 <- GPL96[!is.na(GPL96)]
GPL96 <- names(GPL96)

GPL570 <- jmap(chip = "hgu133plus2", symbol = signature$Symbol)
GPL570 <- GPL570[!is.na(GPL570)]
GPL570 <- names(GPL570)

# list of files

files <- list.files(path = input.dir,
                    recursive = TRUE,
                    pattern = "expression_data",
                    full.names = TRUE)

# remove the tumour compendium file itself from this list

files <- files[-c(grep(pattern = "tumour_compendium_expression_data", x = files))]

# read-in and combine files into a data frame ensuring genes are harmonised

combined_expression_data <- vector(mode = "list", length = length(files))

for (i in seq_along(files)) {
  
  # read file
  
  expression.data <- read.csv(file = files[i],
                              header = TRUE,
                              row.names = 1,
                              stringsAsFactors = FALSE)
  
  # studies using GPL570 reduced to genes represented in GPL96
  
  if (nrow(expression.data) == length(GPL570)) {
    keep <- which(rownames(expression.data) %in% GPL96)
    expression.data <- expression.data[keep, ]
  }
  
  # add to combined data list
  
  combined_expression_data[[i]] <- expression.data
  
}

# convert list to data frame

combined_expression_data <- data.frame(combined_expression_data)

# clean-up environment

rm(list = c("expression.data", "files", "i", "keep", "GPL96", "GPL570",
            "sample.info", "signature", "study.batch.labels", "study.lab"))


## filtration of data

# ensure there are no NAs in the sample data

combined_sample_info <- combined_sample_info[!is.na(combined_sample_info$DMFS_TIME), ]

# remove duplicated patients (those that appear in multiple GEO cohorts)

combined_sample_info <- combined_sample_info[!duplicated(combined_sample_info$SAMPLE_ID), ]

# remove any NAs and duplicated patients in the expression data

keep.patients <- match(combined_sample_info$SAMPLE_ID, colnames(combined_expression_data))
combined_expression_data <- combined_expression_data[ , keep.patients]

# filter sample data to patients with metastasis only

mets.sample <- combined_sample_info[combined_sample_info$DMFS_EVENT == 1, ]

# filter expression data to these patients

keep.patients <- match(mets.sample$SAMPLE_ID, colnames(combined_expression_data))
mets.expr <- combined_expression_data[ , keep.patients]


## adjust for batch effect

# construct the metadata data frame

metadata <- data.frame("SAMPLES" = mets.sample$SAMPLE_ID,
                       "BATCH" = mets.sample$STUDY)

# create model matrix

corr.model <- model.matrix(~1, data = metadata)

# apply ComBat algorithm

harmonised.matrix <- ComBat(dat = mets.expr,
                            batch = metadata$BATCH,
                            mod = corr.model,
                            par.prior = TRUE,
                            prior.plots = FALSE)


## boxplots to observe batch effect before and after removal

#boxplot(mets.expr, outline = FALSE)
#boxplot(harmonised.matrix, outline = FALSE)


## export data

# sample information

write.csv(x = mets.sample,
          file = paste0(output.dir, "tumour_compendium_sample_info.csv"),
          quote = FALSE,
          row.names = FALSE)

# gene level expression data

write.csv(x = harmonised.matrix,
          file = paste0(output.dir, "tumour_compendium_expression_data.csv"),
          quote = FALSE,
          row.names = TRUE)

