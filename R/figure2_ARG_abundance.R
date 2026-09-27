# 04_figure2_ARG_abundance.R
source("R/00_setup.R")
rgi <- read.delim("data/RGI_BWT_16_samples.tsv",check.names=FALSE)
s1 <- read.delim("data/S1_bacterial_reads_16.tsv",check.names=FALSE)
rgi_dat <- rgi %>% filter(`Reference Allele(s) Identity to CARD Reference Protein (%)`>=90) %>% left_join(s1 %>% select(fileName,total_reads),by="fileName") %>% left_join(meta,by="fileName") %>% mutate(normalized_count=`All Mapped Reads`/total_reads*1e6)
rgi_trimmed <- rgi_dat %>% filter(normalized_count>10)
write.table(rgi_trimmed,"results/tables/ARG_abundance.tsv",sep="\t",quote=FALSE,row.names=FALSE)
# Aitchison PCoA
arg_matrix <- rgi_trimmed %>% select(`ARO Term`,fileName,normalized_count) %>% group_by(`ARO Term`,fileName) %>% summarise(normalized_count=sum(normalized_count),.groups="drop") %>% pivot_wider(names_from=fileName,values_from=normalized_count,values_fill=0) %>% column_to_rownames("ARO Term")
arg_pseudocount <- min(arg_matrix[arg_matrix>0])*0.0065
arg_dist <- vegdist(t(arg_matrix),method="aitchison",pseudocount=arg_pseudocount)
arg_pcoa <- ecodist::pco(arg_dist)
arg_percent <- arg_pcoa$values/sum(arg_pcoa$values)*100
arg_pcoa_df <- tibble(fileName=rownames(arg_pcoa$vectors),PCo1=arg_pcoa$vectors[,1],PCo2=arg_pcoa$vectors[,2]) %>% left_join(meta,by="fileName")
write.table(arg_pcoa_df,"results/tables/Figure2A_ARG_Aitchison_coordinates.tsv",sep="\t",quote=FALSE,row.names=FALSE)
# Drug classes: exclude records with missing Drug Class only for class level analysis
betalactam <- c("penam","cephalosporin","carbapenem","cephamycin","penem","monobactam")
MLS <- c("macrolide","lincosamide","streptogramin","streptogramin A","streptogramin B")
drug_class_dat <- rgi_trimmed %>% filter(!is.na(`Drug Class`),str_trim(`Drug Class`)!="") %>% mutate(drug_class=str_replace_all(`Drug Class`," antibiotic",""),drug_class_multidrug=drug_class)
for(i in seq_len(nrow(drug_class_dat))){x<-strsplit(drug_class_dat$drug_class[i],"; ")[[1]];x[x%in%betalactam]<-"Beta-lactam";x[x%in%MLS]<-"MLS";x<-unique(x);drug_class_dat$drug_class_multidrug[i]<-if(length(x)>=2)"multidrug" else x[1]}
rank <- drug_class_dat %>% group_by(city_state,drug_class_multidrug) %>% summarise(total_abundance=sum(normalized_count),.groups="drop")
top10 <- rank %>% group_by(city_state) %>% slice_max(total_abundance,n=10,with_ties=FALSE) %>% ungroup()
top20 <- rank %>% group_by(city_state) %>% slice_max(total_abundance,n=20,with_ties=FALSE) %>% ungroup()
write.table(top10,"results/tables/Figure2B_top10_drug_classes.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(top20,"results/tables/Figure2B_top20_drug_classes.tsv",sep="\t",quote=FALSE,row.names=FALSE)