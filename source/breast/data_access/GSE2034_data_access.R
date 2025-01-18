## GSE2034 data access - tumour data set compendium
## Lewis Quayle (l.quayle@sheffield.ac.uk)
## 2022-07-19


## load libraries

library("dplyr")
library("GEOquery")
library("jetset")
library("org.Hs.eg.db")


## define directory hierarchy

# data set
data_set <- "GSE2034"

# project name
project_name <- "brca_microarray_tumour_compendium"

# project directory
project_dir <- file.path("~", "Documents", "Bioinformatics", project_name)

# output directory
out_dir <- file.path(project_dir, "data", data_set)


## create output directory if it doesn't exist

if (dir.exists(out_dir)) {
  
  print("directory exists")
  
} else {
  
  print("creating directory")
  dir.create(out_dir)
  
}


## access data set

# get GEO object
gse <- getGEO(data_set)

# check data sets present and assign required
length(gse)
gse <- gse[[1]]

# view data and experiment description

#gse
#experimentData(gse)
#abstract(gse)


## extract eData set components

# sample information
sample_info <- read.csv(file = file.path(project_dir, "data", "suppl", "GSE2034_suppl.csv"),
                        header = TRUE,
                        stringsAsFactors = FALSE)

# gene annotation
gene_annot <- fData(gse) 

# gene level expression data
gle_data <- exprs(gse)


## examine the expression data to check normalisation and scale

# scale and normalisation
#summary(gle.data)
#boxplot(gle.data, outline = FALSE)

# log2 transform (data is pre-normalised - see pData(gse))
gle_data <- log2(gle_data)

# check transformation
boxplot(gle_data, outline = FALSE)


## Filtration of the phenotypic data (sample information)

## create the phenotypic data (sample information)

pheno_data <- sample_info %>%
  # select the sample information required
  select(
    "GEO.asscession.number",
    "time.to.relapse.or.last.follow.up..months.",
    "relapse..1.True."
  ) %>%
  # order patient IDs in increasing order
  arrange("GEO.asscession.number") %>%
  # change column names to format used across all DMFS data sets
  set_names(c("SAMPLE_ID", "DMFS_TIME", "DMFS_EVENT")) %>%
  # convert DMFS_TIME from months to years
  mutate("DMFS_TIME" = as.numeric(DMFS_TIME) / 12) %>%
  # create indicators for early (<5 years) and late (>5 years) distant metastasis
  mutate(
    "DMFS_EVENT_EARLY" = if_else(DMFS_EVENT == 1 & DMFS_TIME < 5, true = 1, false = 0),
    "DMFS_EVENT_LATE" = if_else(DMFS_EVENT == 1 & DMFS_TIME >= 5, true = 1, false = 0)
  )


## subset the expression data to patients that exist in the phenotypic data set

expr_data <- gle_data %>%
  # coerce to data frame
  as.data.frame() %>%
  # select only the columns that exist in the photypic data
  select(pheno_data$SAMPLE_ID)


## define jetset probe set for all Entrez IDs in the org.Hs.eg.db database

jetset <- org.Hs.eg.db %>%
  # get all IDs in the database
  keys(keytype = "ENTREZID") %>%
  # map these to the probes in GPL96
  jmap(chip = "hgu133a", eg = .) %>%
  # remove any unmatched probes 
  na.omit()


## filtration of expression matrix to the jetset probe set

expr_data <- expr_data %>%
  # subset to probes in the jetset probe set 
  filter(rownames(expr_data) %in% jetset) %>%
  # change rownames to gene IDs
  set_rownames(names(jetset))


## export data

# phenotypic data
write.csv(
  x = pheno_data,
  file = file.path(out_dir, paste0(data_set, "_pheno.csv")),
  quote = FALSE,
  row.names = FALSE
)

# gene level expression data
write.csv(
  x = expr_data,
  file = file.path(out_dir, paste0(data_set, "_expr.csv")),
  quote = FALSE,
  row.names = TRUE
)
