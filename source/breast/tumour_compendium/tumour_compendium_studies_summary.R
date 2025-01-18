## Summary table for tumour compendium and constituent data sets
## Lewis Quayle
## 27/02/2021


## load required packages

# NONE


## set project directory hierarchy

# home directory

home.dir <- "/Users/Lewis/"

# parent directories

rna.seq.dir <- paste0(home.dir, "Documents/Bioinformatics/rna_seq_analysis/")
parent.dir <- paste0(rna.seq.dir, "public_clinical_data_analysis/")

# input directory

input.dir <- paste0(parent.dir, "data/")

# output directory

out.dir <- paste0(parent.dir, "results/tabular/")


## read-in tumour compendium sample information

# list of files

files <- list.files(path = input.dir,
                    recursive = TRUE,
                    pattern = ".csv",
                    full.names = TRUE)

# specify file to read

file <- files[grep("tumour_compendium_sample_info", x = files)]

# add to global environment

tc.sample.info <- read.csv(file = file,
                           header = TRUE,
                           stringsAsFactors = FALSE)

# extract vector of study accessions

studies <- unique(tc.sample.info$STUDY)

# append the entire compendium to the list

studies <- c(studies, "tc.sample.info")

# clean-up environment

rm(list = c("file", "files"))


## calculate summary statistics

# storage object

study.stats <- list()

# extract / calculate descriptive stats for each study

for (study in studies) {
  
  # if not the entire compendium then subset to a specific study
  
  if (study != "tc.sample.info") {
    
    study.data <- tc.sample.info[tc.sample.info$STUDY == study, ]
    
  } else  {
    
    study.data <- tc.sample.info
    
  }
  
  # number of cases (all have dist. metastases)
  
  cases <- nrow(study.data)
  
  # number of early (<5 years) distant met. recurrence events
  
  dm.early <- sum(study.data$DMFS_EVENT_EARLY == 1)
  
  # percentage proportion early distant met. recurrence events
  
  dm.early.per.cases <- (dm.early / cases) * 100
  
  # number of late (>5 years) distant met. recurrence events
  
  dm.late <- sum(study.data$DMFS_EVENT_LATE == 1)
  
  # percentage proportion late distant met. recurrence events
  
  dm.late.per.cases <- (dm.late / cases) * 100
  
  # store all extracted metrics
  
  study.stats[[study]] <- c(cases,
                            dm.early,
                            dm.late,
                            dm.early.per.cases,
                            dm.late.per.cases)
  
}


# coerce results to a data frame

study.stats <- as.data.frame(study.stats,
                             row.names = c("No. Cases",
                                           "No. Cases Early Dist. Met. (<5 yrs)",
                                           "No. Cases Late Dist. Met. (>5 yrs)",
                                           "Early Dist. Met. (% of all cases)",
                                           "Late Dist. Met. (% of all cases)"))

# rename the tumour compendium column

colnames(study.stats)[ncol(study.stats)] <- "Compendium"

# transpose the data frame, moving the study ID from rownames to a column

study.stats <- data.frame("Study ID" = colnames(study.stats),
                          t(study.stats),
                          row.names = NULL,
                          check.names = FALSE)

# round the numeric columns to 2 decimal places

study.stats[ , -c(1:4)] <- round(x = study.stats[ , -c(1:4)], digits = 2)


## export

write.csv(x = study.stats,
          file = paste0(out.dir, "tumour_compendium_studies_summary.csv"),
          quote = FALSE,
          row.names = FALSE)

