# figure1C_aitchison.R
source("R/00_setup.R")
dat <- read.delim("data/bracken_species_normalized.tsv", check.names=FALSE)
species_matrix <- dat %>% select(name,fileName,normalized_count) %>% group_by(name,fileName) %>% summarise(normalized_count=sum(normalized_count),.groups="drop") %>% pivot_wider(names_from=fileName,values_from=normalized_count,values_fill=0) %>% column_to_rownames("name")
pseudocount <- min(species_matrix[species_matrix>0])*0.0065
aitchison_dist <- vegdist(t(species_matrix),method="aitchison",pseudocount=pseudocount)
aitchison_pcoa <- ecodist::pco(aitchison_dist)
variance_explained <- aitchison_pcoa$values/sum(aitchison_pcoa$values)*100
pcoa_df <- tibble(fileName=rownames(aitchison_pcoa$vectors),PCo1=aitchison_pcoa$vectors[,1],PCo2=aitchison_pcoa$vectors[,2]) %>% left_join(meta,by="fileName")
write.table(pcoa_df,"results/tables/Figure1C_Aitchison_coordinates.tsv",sep="\t",quote=FALSE,row.names=FALSE)
p1b <- ggplot(pcoa_df,aes(PCo1,PCo2)) + stat_ellipse(aes(colour=city_state),level=0.95,linewidth=0.8,show.legend=FALSE) + geom_point(aes(fill=city_state),shape=21,colour="black",size=4.2,stroke=0.45) + scale_fill_manual(values=city_colors) + scale_colour_manual(values=city_colors) + labs(x=paste0("PCo1 (",round(variance_explained[1],2),"%)"),y=paste0("PCo2 (",round(variance_explained[2],2),"%)"),fill="City") + theme_classic(base_size=14)
ggsave("results/figures/Figure1C_Aitchison_PCoA.png",p1b,width=7.5,height=6.2,units="in",dpi=600,bg="white")