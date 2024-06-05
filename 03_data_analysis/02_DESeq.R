
#' ---
#' title: "EV track pilot exp - no trimming - DEseq2"
#' author: "Marieke Vromman"
#' output: 
#'    html_document:
#'       toc: TRUE
#'       toc_float: TRUE
#'       theme: paper
#'       df_print: paged
#'       highlight: tango
#' ---
#' 

#' # File set-up

#' ## Set working directory to current directory
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

#' ## Load standard libraries and resolves conflicts
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer('intersect', 'dplyr')
conflicts_prefer(dplyr::count)

library("DESeq2")
library(ggrepel)


#' ## Set figure theme
mytheme = theme_bw(base_size = 10) + 
  theme(text = element_text(size=10, colour='black'),
        title = element_text(size=10, colour='black'),
        line = element_line(size=0.5),
        axis.title = element_text(size=10, colour='black'),
        axis.text = element_text(size=10, colour='black'),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size=0.5),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.text=element_text(size=10)) 

mytheme_discrete_x = mytheme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#' ## Read data

counts = 
  read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T43.R1/D1515T43.R1_rc.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
  mutate(sample = '250_uM_S4U_1', sample_ID = "D1515T43", conc = 250, seq_read = 'R1', replicate = 1) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T43.R2/D1515T43.R2.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '250_uM_S4U_1', sample_ID = "D1515T43", conc = 250, seq_read = 'R2', replicate = 1)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T44.R1/D1515T44.R1_rc.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '250_uM_S4U_2', sample_ID = "D1515T44", conc = 250, seq_read = 'R1', replicate = 2)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T44.R2/D1515T44.R2.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '250_uM_S4U_2', sample_ID = "D1515T44", conc = 250, seq_read = 'R2', replicate = 2)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T45.R1/D1515T45.R1_rc.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '125_uM_S4U_1', sample_ID = "D1515T45", conc = 125, seq_read = 'R1', replicate = 1)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T45.R2/D1515T45.R2.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '125_uM_S4U_1', sample_ID = "D1515T45", conc = 125, seq_read = 'R2', replicate = 1)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T46.R1/D1515T46.R1_rc.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '125_uM_S4U_2', sample_ID = "D1515T46", conc = 125, seq_read = 'R1', replicate = 2)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T46.R2/D1515T46.R2.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '125_uM_S4U_2', sample_ID = "D1515T46", conc = 125, seq_read = 'R2', replicate = 2)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T47.R1/D1515T47.R1_rc.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '0_uM_S4U_2', sample_ID = "D1515T47", conc = 0, seq_read = 'R1', replicate = 2)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T47.R2/D1515T47.R2.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '0_uM_S4U_2', sample_ID = "D1515T47", conc = 0, seq_read = 'R2', replicate = 2)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T48.R1/D1515T48.R1_rc.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '0_uM_S4U_3', sample_ID = "D1515T48", conc = 0, seq_read = 'R1', replicate = 3)) %>%
  bind_rows(read_tsv('../01_data/20231019_SlamDunk_rc_run/D1515T48.R2/D1515T48.R2.fastq_slamdunk_mapped_filtered_tcount.tsv', skip = 2) %>%
              mutate(sample = '0_uM_S4U_3', sample_ID = "D1515T48", conc = 0, seq_read = 'R2', replicate = 3))

counts

#' remove all transcripts that were not found
counts = counts %>% filter(ReadCount > 0)

counts

#' add annotation

anno = read_tsv('../../../scripts/gencode.v44.chr_patch_hapl_scaff.annotation_transcript_type.bed')

counts = counts %>% 
  rename('chr' = Chromosome,
         'start' = Start,
         'end' = End,
         'transcript_id' = 'Name',
         'strand' = Strand) %>%
  left_join(anno) 

#' simplify gene types

counts %>% count(transcript_type) %>% arrange(desc(n))

counts = counts %>%
  mutate(transcript_type_cat = ifelse(transcript_type == 'protein_coding', 'protein_coding', 'other'),
         transcript_type_cat = ifelse(transcript_type == 'lncRNA', 'lncRNA', transcript_type_cat))
         
counts$transcript_type_cat = factor(counts$transcript_type_cat, levels = c('protein_coding', 'lncRNA', 'other'))

counts

counts %>% count(transcript_type_cat) %>% arrange(desc(n))



#' # Data-analysis


counts_table = counts %>%
  filter(seq_read == 'R2') %>%
  filter(conc == 0 | conc == 250) %>%
  mutate(sample = ifelse(conc == 0, 'untreated', 'treated'),
         sample = paste(sample, replicate, sep = '_')) %>%
  select(transcript_id, sample, ReadCount) %>%
  pivot_wider(names_from = sample, values_from = ReadCount)


counts_table

counts_matrix = as.matrix(counts_table %>% select(-transcript_id))
rownames(counts_matrix) = counts_table$transcript_id

counts_matrix


# remove NAs

counts_matrix[is.na(counts_matrix)] <- 0

# column info

coldata = counts %>%
  filter(conc == 0 | conc == 250) %>%
  select(sample, conc, replicate) %>% unique() %>%
  mutate(condition = ifelse(conc == 0, 'untreated', 'treated'),
         sample = paste(condition, replicate, sep = '_')) %>%
  select(-conc, -replicate) %>%
  mutate(sample = as.factor(sample),
         condition = as.factor(condition))

coldata_matrix = as.matrix(coldata %>% select(-sample))

rownames(coldata_matrix) = coldata$sample

coldata = coldata %>% select(-sample)

#' check if sample names are in same order

all(rownames(coldata_matrix) == colnames(counts_matrix))

dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = coldata_matrix,
                              design = ~ condition)
dds = DESeq(dds)

dds

dds = results(dds)
dds

summary(dds)

diff_exp = as.data.frame(dds) %>%
  rownames_to_column("transcript_id")

diff_exp

# add gene names

diff_exp = diff_exp %>%
  left_join(read_tsv('../../../scripts/gencode.v44.chr_patch_hapl_scaff.annotation_info.txt') %>%
              select(transcript_id, gene_name))


diff_exp = diff_exp %>%
  mutate(diffexpressed = "NO") %>%
  mutate(diffexpressed = ifelse(log2FoldChange > 0.6 & padj < 0.05, "UP", diffexpressed)) %>%
  mutate(diffexpressed = ifelse(log2FoldChange < -0.6 & padj < 0.05, "DOWN", diffexpressed)) %>%
  mutate(label_on = NA) %>%
  mutate(label_on = ifelse(-log10(padj) > 22, gene_name, NA))

diff_exp %>% filter(!is.na(label_on))

#write_tsv(res_df_AGO %>% select(-label_on), "diff_exp/IP_vs_T_all.tsv")

ggplot(data=diff_exp, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label = label_on)) + 
  geom_point() +
  theme_minimal() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("coral", "black", "darkolivegreen3")) +
  geom_text_repel()

counts %>% 
  group_by(sample, seq_read) %>%
  summarise(med_count = median(ReadCount))

diff_exp %>% filter(gene_name == 'HMGCS1')
counts %>% filter(transcript_id == "ENST00000325110.11") %>% view()


