# 01_prepare_species_data.R
source("R/00_setup.R")
bracken <- read.delim("data/bracken_species.tsv", check.names=FALSE)
s1 <- read.delim("data/S1_bacterial_reads_16.tsv", check.names=FALSE)
dat <- bracken %>%
  left_join(s1 %>% select(fileName,total_reads,bacterial_reads,total_mapped_reads), by="fileName") %>%
  left_join(meta, by="fileName") %>%
  mutate(normalized_count=(new_est_reads/bacterial_reads)*1e6)
stopifnot(n_distinct(dat$fileName)==16)
stopifnot(sum(is.na(dat$bacterial_reads))==0)
write.table(dat,"data/bracken_species_normalized.tsv",sep="\t",quote=FALSE,row.names=FALSE)