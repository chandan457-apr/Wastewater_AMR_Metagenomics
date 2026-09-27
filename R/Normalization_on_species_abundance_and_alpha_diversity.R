# ============================================================
# normalization_on_species_abundance_and_alpha_diversity.R
#
# Normalization of Bracken species abundance
# and calculation of alpha diversity
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(tidyverse)


# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

bracken_file <- "data/bracken_species.tsv"

s1_file <- "data/S1_bacterial_reads_16.tsv"


# ------------------------------------------------------------
# 3. Read Bracken species abundance
# ------------------------------------------------------------

bracken <- read.delim(
  bracken_file,
  check.names = FALSE
)

cat("Bracken dimensions:\n")
print(dim(bracken))

cat("\nBracken columns:\n")
print(colnames(bracken))


# ------------------------------------------------------------
# 4. Read bacterial read information
# ------------------------------------------------------------

s1 <- read.delim(
  s1_file,
  check.names = FALSE
)

cat("\nS1 dimensions:\n")
print(dim(s1))

cat("\nS1 columns:\n")
print(colnames(s1))


# ------------------------------------------------------------
# 5. Define sample metadata
# ------------------------------------------------------------

meta <- tribble(
  ~fileName, ~SRR, ~city_state, ~month_label, ~site,
  
  "OSD104",  "SRR37020310", "Chennai",
  "2022-04", "Perungudi",
  
  "OSD351",  "SRR37020427", "Chennai",
  "2022-10", "Nesapakkam",
  
  "OSD613",  "SRR37020476", "Chennai",
  "2023-04", "Koyambedu",
  
  "OSD1247", "SRR37020297", "Chennai",
  "2023-10", "Alandur",
  
  "OSD17",   "SRR37020554", "Delhi",
  "2022-04", "Najafgarh Drain",
  
  "OSD417",  "SRR37020539", "Delhi",
  "2022-10", "Coronation Pillar",
  
  "OSD572",  "SRR37020546", "Delhi",
  "2023-04", "Coronation Pillar",
  
  "OSD1251", "SRR37020391", "Delhi",
  "2023-10", "Shahadra Drain",
  
  "OSD1",    "SRR37020636", "Kolkata",
  "2022-04", "Kestopur Khal",
  
  "OSD454",  "SRR37020504", "Kolkata",
  "2022-10", "Kestopur Khal",
  
  "OSD563",  "SRR37020650", "Kolkata",
  "2023-04", "Anandapur Khal",
  
  "OSD1196", "SRR37020491", "Kolkata",
  "2023-10", "Kestopur Khal",
  
  "OSD257",  "SRR37020367", "Mumbai",
  "2022-04", "LBS Marg, Bhandup",
  
  "OSD434",  "SRR37020467", "Mumbai",
  "2022-10", "Nehru Nagar Kurla",
  
  "OSD636",  "SRR37020633", "Mumbai",
  "2023-04", "KJ Somaiya Ghatkopar",
  
  "OSD1278", "SRR37020678", "Mumbai",
  "2023-10", "Pratiksha Nagar Wadala"
)


# ------------------------------------------------------------
# 6. Join Bracken with bacterial read information
# ------------------------------------------------------------

dat <- bracken %>%
  
  left_join(
    s1 %>%
      select(
        fileName,
        total_reads,
        bacterial_reads,
        total_mapped_reads
      ),
    by = "fileName"
  ) %>%
  
  left_join(
    meta,
    by = "fileName"
  )


# ------------------------------------------------------------
# 7. Check that every sample received bacterial read data
# ------------------------------------------------------------

cat("\nNumber of samples:\n")
print(n_distinct(dat$fileName))

cat("\nMissing bacterial-read values:\n")
print(sum(is.na(dat$bacterial_reads)))


# ------------------------------------------------------------
# 8. Normalize Bracken abundance
# ------------------------------------------------------------
#
# normalized_count =
#      new_est_reads
#      -------------
#      bacterial_reads
#          × 1,000,000
#
# This gives estimated reads per million bacterial reads.

dat <- dat %>%
  mutate(
    normalized_count =
      (new_est_reads / bacterial_reads) * 1e6
  )


# ------------------------------------------------------------
# 9. Inspect one sample
# ------------------------------------------------------------

cat("\nTop taxa in OSD563:\n")

dat %>%
  filter(fileName == "OSD563") %>%
  select(
    name,
    new_est_reads,
    bacterial_reads,
    normalized_count
  ) %>%
  arrange(desc(normalized_count)) %>%
  head(10) %>%
  print()


# ------------------------------------------------------------
# 10. Save normalized Bracken table
# ------------------------------------------------------------

write.table(
  dat,
  "project_data/processed_data/bracken_species_normalized.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# ALPHA DIVERSITY
# ============================================================


# ------------------------------------------------------------
# 11. Prepare bacterial taxonomic data
# ------------------------------------------------------------

dat_bac <- dat


# ------------------------------------------------------------
# 12. Calculate alpha diversity
# ------------------------------------------------------------

alpha_df <- dat_bac %>%
  
  group_by(
    month_label,
    city_state
  ) %>%
  
  summarise(
    
    # ------------------------------------------
    # Richness
    # Number of taxa with non zero abundance
    # ------------------------------------------
    
    richness =
      sum(normalized_count > 0),
    
    # ------------------------------------------
    # Relative abundance within each sample
    # ------------------------------------------
    
    # Shannon
    shannon =
      -sum(
        (
          normalized_count /
            sum(normalized_count)
        ) *
          log(
            normalized_count /
              sum(normalized_count)
          )
        [
          normalized_count > 0
        ]
      ),
    
    # ------------------------------------------
    # Simpson diversity
    #
    # This is 1 - Simpson dominance
    # Therefore:
    # higher = greater diversity
    # ------------------------------------------
    
    simpson =
      1 -
      sum(
        (
          normalized_count /
            sum(normalized_count)
        )^2
      ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 13. Display alpha diversity results
# ------------------------------------------------------------

cat("\nAlpha-diversity results:\n")

print(alpha_df)


# ------------------------------------------------------------
# 14. Save alpha diversity table
# ------------------------------------------------------------

write.table(
  alpha_df,
  "results/tables/Figure1_alpha_diversity.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Optional summary by city
# ------------------------------------------------------------

city_summary <- alpha_df %>%
  
  group_by(city_state) %>%
  
  summarise(
    mean_richness = mean(richness),
    mean_shannon = mean(shannon),
    mean_simpson = mean(simpson),
    
    .groups = "drop"
  )


cat("\nMean alpha diversity by city:\n")

print(city_summary)