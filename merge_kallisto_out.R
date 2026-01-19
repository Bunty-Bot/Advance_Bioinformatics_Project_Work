#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
library("tximport")
#library("BiocManager", include.only = c("tximport"))
#see documentation
#browseVignettes("tximport")
#treat kallisto file
#files <- file.path(dir, "kallisto", samples$run, "abundance.tsv.gz")
#txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
#ignoreAfterBar = True is for those genes that have name_gene1 | name_gene_2
#files in this case would be a zipped directory containing abundance.tsv for each of the samples
samples <- list('SRR24630891', 'SRR24630892', 'SRR24630893', 'SRR24630894',
             'SRR24630896', 'SRR24630897', 'SRR24630898', 'SRR24630899',
             'SRR24630900', 'SRR24630902', 'SRR24630903', 'SRR24630901',
             'SRR24630905', 'SRR24630907', 'SRR24630908', 'SRR24630904',
             'SRR24630909', 'SRR24630910', 'SRR24630911', 'SRR24630912',
             'SRR24630913', 'SRR24630914', 'SRR24630915', 'SRR24630916',
             'SRR24630918', 'SRR24630919', 'SRR24630920', 'SRR24630921',
             'SRR24630922', 'SRR24630923')

#901 and 904 had to be run manually bc they produced no abundance.* files when run via the snakefile

files <- list()
for (s in samples){
  dir_name <- paste0(s, "_kallisto")
  file <- file.path(dir_name, "abundance.h5")
  #print(file)
  files <- append(files, file)
}
#i think that 'samples' in the og code refers to the csv containing the meta data,
#and you need to look for the name of the sample in the run column of the meta data

files <- unlist(files)
names(files) <- paste0(samples)
txi.kallisto.tsv <- tximport(files, type = "kallisto",txOut = TRUE, ignoreAfterBar = TRUE)

final_df <- as.data.frame(txi.kallisto.tsv$abundance)
final_df$gene_id <- rownames(final_df)

final_df <- final_df[, c("gene_id", unlist(samples))]

#head(txi.kallisto.tsv$abundance)
write.table(
    final_df,
    file = "kallisto_combined_tpm.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
