# 00_setup.R
# R packages required to run these codes. 
library(igraph)
library(tibble)
library(tidyverse)
library(dplyr)
library(janitor)
library(reshape2)
library(RColorBrewer)
library(ALDEx2)
library(abdiv)
library(ggalt)
library(paletteer)
library(ggpubr)
library(tidyplots)
library(ggsankey)
library(cowplot)
library(vegan)
library(ecodist)
library(UpSetR)
library(ggplotify)


meta <- tribble(
  ~fileName, ~SRR, ~city_state, ~month_label, ~site,
  "OSD104","SRR37020310","Chennai","2022-04","Perungudi",
  "OSD351","SRR37020427","Chennai","2022-10","Nesapakkam",
  "OSD613","SRR37020476","Chennai","2023-04","Koyambedu",
  "OSD1247","SRR37020297","Chennai","2023-10","Alandur",
  "OSD17","SRR37020554","Delhi","2022-04","Najafgarh Drain",
  "OSD417","SRR37020539","Delhi","2022-10","Coronation Pillar",
  "OSD572","SRR37020546","Delhi","2023-04","Coronation Pillar",
  "OSD1251","SRR37020391","Delhi","2023-10","Shahadra Drain",
  "OSD1","SRR37020636","Kolkata","2022-04","Kestopur Khal",
  "OSD454","SRR37020504","Kolkata","2022-10","Kestopur Khal",
  "OSD563","SRR37020650","Kolkata","2023-04","Anandapur Khal",
  "OSD1196","SRR37020491","Kolkata","2023-10","Kestopur Khal",
  "OSD257","SRR37020367","Mumbai","2022-04","LBS Marg, Bhandup",
  "OSD434","SRR37020467","Mumbai","2022-10","Nehru Nagar Kurla",
  "OSD636","SRR37020633","Mumbai","2023-04","KJ Somaiya Ghatkopar",
  "OSD1278","SRR37020678","Mumbai","2023-10","Pratiksha Nagar Wadala"
)
