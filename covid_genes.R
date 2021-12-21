gtf <- rtracklayer::import('gtf/Sars_cov_2.ASM985889v3.101.gtf')
gtf_df=as.data.frame(gtf)
covid_genes <- gtf_df$gene_name %>% unique()
saveRDS(covid_genes, "rObjects/covid_genes.rds")
