
#' ---
#' title: "EV track pilot exp - no trimming"
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
#' ## How many of the transcripts are labeled?
#' how many transcript have at least one T > C conversion

counts %>% 
  group_by(sample) %>%
  summarise(perc_conversed = 100 * sum(ConversionRate > 0) / n())

100*nrow(counts %>% filter(ConversionRate > 0)) / nrow(counts)


#' bar plot
counts %>% 
  filter(ReadCount >= 5) %>%
  group_by(sample, conc) %>%
  summarise(perc_conversed = sum(ConversionRate >= 0.01) / n()) %>%
  mutate(conc = paste(conc, 'µM', sep = ' ')) %>% 
  ggplot(aes(sample, perc_conversed, fill = conc)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  scale_fill_manual(values = c('grey', 'lightblue2', '#a1d99b')) +
  labs(fill = 'S4U concentration') +
  xlab('') +
  ylab('% of labeled transcripts')


#' bar plot per read 

perc_conversed = counts %>%
  group_by(sample, seq_read, conc) %>%
  summarise(nr_transcripts = n(),
            nr_conversed = sum(ConversionsOnTs > 0),
            perc_conversed = nr_conversed/nr_transcripts)

perc_conversed

perc_conversed %>%
  ggplot(aes(sample, perc_conversed, fill = seq_read)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  facet_wrap(~conc, scales = 'free_x')

#' with filter on minimal conversion rate


perc_conversed = counts %>%
  #filter(ReadCount >= 5) %>%
  group_by(sample, conc) %>%
  summarise(nr_transcripts = n(),
            #nr_conversed = sum(ConversionsOnTs > 0),
            nr_conversed = sum(ConversionRate > 0.01),
            perc_conversed = nr_conversed/nr_transcripts) 

perc_conversed

perc_conversed %>%
  mutate(conc = as.character(conc)) %>%
  ggplot(aes(sample, perc_conversed, fill = conc)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  #facet_wrap(~conc, scales = 'free_x') +
  scale_fill_manual(values = c('grey', 'lightblue2', '#a1d99b')) +
  geom_text(aes(label = sprintf("%0.1f", round(100 * perc_conversed, digits = 4))), vjust = -0.2) +
  labs(fill = 'S4U concentration') +
  xlab('') +
  ylab('% of labeled transcripts')




#' how many labeled per gene category

perc_conversed = counts %>%
  filter(ReadCount >= 5) %>%
  group_by(sample, conc, transcript_type_cat) %>%
  summarise(nr_transcripts = n(),
            nr_conversed = sum(ConversionRate > 0.01),
            perc_conversed = nr_conversed/nr_transcripts) 


perc_conversed

perc_conversed %>%
  mutate(conc = as.character(conc)) %>%
  ggplot(aes(sample, perc_conversed, fill = conc)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c('grey', 'lightblue2', '#a1d99b')) +
  mytheme_discrete_x +
  facet_wrap(~transcript_type_cat, scales = 'free_x') +
  labs(fill = 'S4U concentration') +
  xlab('') +
  ylab('% of labeled transcripts') +
  geom_text(aes(label = sprintf("%0.1f", round(100 * perc_conversed, digits = 4))), vjust = -0.2)

perc_conversed %>% filter(sample == '125_uM_S4U_1') %>%


#' with filter on minimal conversion rate

perc_conversed = counts %>%
  #filter(ReadCount > 99) %>%
  group_by(sample, seq_read, conc, transcript_type_cat) %>%
  summarise(nr_transcripts = n(),
            #nr_conversed = sum(ConversionsOnTs > 0),
            nr_conversed = sum(ConversionRate > 0.01),
            perc_conversed = nr_conversed/nr_transcripts) 


perc_conversed

perc_conversed %>%
  ggplot(aes(sample, perc_conversed, fill = transcript_type_cat)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  facet_wrap(~conc, scales = 'free_x')


# check the other annotations

counts %>%
  filter(ReadCount > 4) %>%
  group_by(sample, conc, transcript_type) %>%
  #filter(n() > 9) %>%
  summarise(nr_transcripts = n(),
            nr_conversed = sum(ConversionRate > 0.01),
            perc_conversed = nr_conversed/nr_transcripts) %>%
  arrange(desc(transcript_type)) %>%
  ggplot(aes(reorder(transcript_type,-nr_transcripts), perc_conversed, fill = transcript_type)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(sample)) + 
  #mytheme_discrete_x +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position = "",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_text(aes(label=nr_transcripts), angle = 90, hjust = 1, vjust = 1,
            y = 1)
  
  

#' ## Conversion rate

#' median conversion rate
counts %>% 
  group_by(sample, seq_read) %>%
  summarise(median_conversionrate = median(ConversionRate))


#' bar plot
counts %>% 
  mutate(conc = as.character(conc)) %>%
  group_by(sample, conc) %>%
  summarise(median_conversionrate = median(ConversionRate)) %>%
  ggplot(aes(sample, median_conversionrate, fill = conc)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  scale_fill_manual(values = c('lightcoral', 'forestgreen', 'deepskyblue')) +
  facet_wrap(~conc, scales = 'free_x') +
  geom_text(aes(label = sprintf("%0.2f", round(100 * median_conversionrate, digits = 4))), vjust = -0.2)

  
#' bar plot with separate seq reads
counts %>% 
  mutate(conc = as.character(conc),
         conc_2 = paste(conc, '_', seq_read, sep = '')) %>%
  group_by(sample, seq_read, conc_2, conc) %>%
  summarise(median_conversionrate = median(ConversionRate)) %>%
  ggplot(aes(sample, median_conversionrate, fill = conc_2)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  scale_fill_manual(values = c('mistyrose1', 'lightcoral', 'lightblue', 'deepskyblue','lightgreen', 'forestgreen')) +
  #geom_text(aes(label = sprintf("%0.2f", round(100 * median_conversionrate, digits = 4))), vjust = -0.2) + 
  facet_wrap(~conc, scales = 'free_x')


#' boxplot
counts %>%
  #filter(ReadCount >= 5) %>%
  mutate(conc = as.character(conc)) %>%
  mutate(conc = paste(conc, 'µM', sep = ' ')) %>% 
  ggplot(aes(sample, ConversionRate, fill = conc)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_y_log10(labels = scales::percent_format()) +
  scale_fill_manual(values = c('grey', 'lightblue2', '#a1d99b')) +
  mytheme_discrete_x +
  labs(fill = 'S4U concentration') +
  xlab('') +
  ylab('conversion rate (log10 scale)') +
  geom_hline(yintercept = 0.01, color = 'red')
  
counts %>% group_by(sample) %>%
  filter(ConversionRate > 0) %>%
  summarise(conv_med = median(ConversionRate))

#' boxplot  with separate seq reads
counts %>%
  mutate(conc = as.character(conc),
         conc_2 = paste(conc, '_', seq_read, sep = '')) %>%
  ggplot(aes(sample, ConversionRate, fill = conc_2)) +
  geom_boxplot() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c('mistyrose1', 'lightcoral', 'lightblue', 'deepskyblue','lightgreen', 'forestgreen')) +
  mytheme_discrete_x


#' density plot
counts %>%
  mutate(conc = as.character(conc)) %>%
  ggplot((aes(ConversionRate, group = sample, color = conc))) +
  geom_density() +
  scale_x_continuous(labels = scales::percent_format())

#' ## link to gene annotation
#' 
#' 
#' barplot
counts %>% 
  mutate(conc = as.character(conc)) %>%
  group_by(sample, conc, transcript_type_cat) %>%
  summarise(median_conversionrate = median(ConversionRate)) %>%
  ggplot(aes(sample, median_conversionrate, fill = transcript_type_cat)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  scale_fill_manual(values = c('lightcoral', 'forestgreen', 'deepskyblue')) +
  facet_wrap(~conc, scales = 'free_x')

#' boxplot

counts %>%
  ggplot(aes(transcript_type_cat, ConversionRate, fill = transcript_type_cat)) +
  geom_boxplot() +
  facet_wrap(~sample, nrow = 3) +
  mytheme_discrete_x +
  scale_y_log10(labels = scales::percent_format())


#' ## correlation read count (CPM) and conversion rate

counts %>%
  mutate(conc = as.character(conc)) %>%
  ggplot(aes(ReadsCPM, ConversionRate, color = conc)) +
  geom_point(alpha = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_log10() +
  facet_wrap(~sample, nrow = 1)

#' ## correlation read count (CPM) and % conversed

counts_cum = counts %>%
  filter(seq_read == "R1") %>%
  select(sample, conc, ConversionRate, ReadCount, transcript_id) %>%
  group_by(sample, conc, ReadCount) %>%
  arrange((ReadCount)) %>%
  summarize(total_n = n(), val_n = sum(ConversionRate > 0.01)) %>%
  ungroup() %>%
  arrange((ReadCount)) %>%
  mutate(total_n_cum = cumsum(total_n), 
         val_n_cum = cumsum(val_n)) %>%
  mutate(perc_val = val_n_cum/total_n_cum)

counts_cum

counts_cum %>% 
  mutate(conc = as.character(conc)) %>%
  ggplot(aes(ReadCount, perc_val, color = conc)) +
  geom_point(size = 1) +
  #geom_smooth(color = '#5AB4E5') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  mytheme +
  facet_wrap(~sample)


count_bin = counts %>%
  mutate(count_bin = ifelse(ReadCount <= 100000, 'bin_100000', 'bin_above'),
         count_bin = ifelse(ReadCount <= 10000, 'bin_10000', count_bin),
         count_bin = ifelse(ReadCount <= 1000, 'bin_1000', count_bin),
         count_bin = ifelse(ReadCount <= 100, 'bin_100', count_bin),
         count_bin = ifelse(ReadCount <= 50, 'bin_50', count_bin),
         count_bin = ifelse(ReadCount <= 10, 'bin_10', count_bin),
         count_bin = ifelse(ReadCount == 1, 'bin_1', count_bin))

count_bin %>% pull(ReadsCPM) %>% quantile()

count_bin$count_bin = factor(count_bin$count_bin,
                             levels = c('bin_1', 'bin_10', 'bin_50', 'bin_100',
                                        'bin_1000', 'bin_10000', 'bin_100000', 
                                        'bin_above'))

count_bin %>%
  group_by(sample, count_bin, conc) %>%
  summarise(perc_conversed = sum(ConversionRate > 0.01)/n()) %>%
  mutate(conc = as.character(conc)) %>%
  ggplot(aes(sample, perc_conversed, fill = count_bin)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  ylab('% of labeled transcripts')


count_bin = counts %>%
  mutate(CPM_bin = ifelse(ReadsCPM <= 10, 'bin_10', 'bin_above'),
         CPM_bin = ifelse(ReadsCPM <= 5, 'bin_5', CPM_bin),
         CPM_bin = ifelse(ReadsCPM <= 1, 'bin_1', CPM_bin),
         CPM_bin = ifelse(ReadsCPM <= 0.5, 'bin_0.5', CPM_bin),
         CPM_bin = ifelse(ReadsCPM <= 0.1, 'bin_0.1', CPM_bin),
         CPM_bin = ifelse(ReadsCPM <= 0.05, 'bin_0.05', CPM_bin),
         CPM_bin = ifelse(ReadsCPM == 0.01, 'bin_0.01', CPM_bin))

count_bin$CPM_bin = factor(count_bin$CPM_bin,
                             levels = c('bin_0.01', 'bin_0.05', 'bin_0.1', 'bin_0.5',
                                        'bin_1', 'bin_5', 'bin_10',
                                        'bin_above'))

count_bin %>%
  group_by(sample, CPM_bin, conc) %>%
  summarise(perc_conversed = sum(ConversionRate > 0.01)/n()) %>%
  mutate(conc = as.character(conc)) %>%
  ggplot(aes(sample, perc_conversed, fill = CPM_bin)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  ylab('% of labeled transcripts')

#' ## choose specific gene of interest

counts %>%
  filter(transcript_type == 'lncRNA') %>%
  group_by(transcript_id) %>%
  filter(n() == 12,
         median(CoverageOnTs) > 1000) %>%
  ungroup() 

counts %>%
  filter(transcript_id == 'ENST00000624927.3') %>%
  ggplot(aes(sample, ConversionRate, color = seq_read, group = seq_read)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x


#' ## distribution counts over gene types

count_bin = counts %>%
  mutate(count_bin = ifelse(ReadCount <= 100000, 'bin_100000', 'bin_above'),
         count_bin = ifelse(ReadCount <= 10000, 'bin_10000', count_bin),
         count_bin = ifelse(ReadCount <= 1000, 'bin_1000', count_bin),
         count_bin = ifelse(ReadCount <= 100, 'bin_100', count_bin),
         count_bin = ifelse(ReadCount <= 50, 'bin_50', count_bin),
         count_bin = ifelse(ReadCount <= 10, 'bin_10', count_bin),
         count_bin = ifelse(ReadCount == 1, 'bin_1', count_bin))


count_bin$count_bin = factor(count_bin$count_bin,
                             levels = c('bin_1', 'bin_10', 'bin_50', 'bin_100',
                                        'bin_1000', 'bin_10000', 'bin_100000', 
                                        'bin_above'))

count_bin %>%
  mutate(conc = as.character(conc)) %>%
  ggplot(aes(transcript_type_cat, fill = count_bin)) +
  geom_bar(position = 'dodge') +
  scale_y_continuous(labels = scales::comma_format()) +
  mytheme_discrete_x +
  ylab('nr of transcripts') 

counts %>%
  ggplot(aes(ReadCount, color = transcript_type_cat)) +
  stat_ecdf() +
  scale_x_log10(labels = scales::comma_format(), 
                breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

#' # % labeled counts plot

counts_cum = counts %>%
  filter(seq_read == 'R1') %>%
  group_by(sample, transcript_type_cat) %>%
  mutate(total_nr_sample = n()) %>%
  group_by(sample, transcript_type_cat, ReadCount, total_nr_sample) %>%
  summarise(total_nr = n(),
            lab_nr = sum(ConversionRate > 0.01)) %>%
  ungroup() %>%
  group_by(sample, transcript_type_cat) %>%
  arrange(desc(ReadCount)) %>%
  mutate(lab_perc_cum = cumsum(lab_nr)/cumsum(total_nr),
         total_perc_sample = cumsum(total_nr) / total_nr_sample) %>%
  ungroup() %>%
  mutate(count_label = 
           ifelse(ReadCount %in% c(1, 3, 5, 10, 100, 1000),
                  sprintf("%d", round(100 * total_perc_sample, digits = 0)),
                  NA),
         count_label = ifelse(ReadCount == 1000 & 
                                transcript_type_cat == 'lncRNA',
                              NA, count_label),
         y_pos = ifelse(sample == '0_uM_S4U_2', 0.2, NA),
         y_pos = ifelse(sample == '0_uM_S4U_3', 0.25, y_pos),
         y_pos = ifelse(sample == '125_uM_S4U_1', 0.55, y_pos),
         y_pos = ifelse(sample == '125_uM_S4U_2', 0.6, y_pos),
         y_pos = ifelse(sample == '250_uM_S4U_1', 0.7, y_pos),
         y_pos = ifelse(sample == '250_uM_S4U_2', 0.75, y_pos)) 

counts_cum

counts_cum_tmp = counts_cum %>%
  filter(transcript_type_cat == 'lncRNA') 

counts_cum_tmp %>%
  ggplot(aes(ReadCount, lab_perc_cum, color = sample)) +
  geom_point() +
  geom_line() +
  #facet_wrap(~transcript_type_cat) +
  scale_x_log10(label = scales::comma_format(),
                breaks = c(1, 3, 5, 10, 100, 1000, 10000, 100000)) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  geom_text(aes(label = count_label),
            vjust = -0.2,
            y = counts_cum_tmp$y_pos,
            size = 3)

counts %>%
  ggplot(aes(ReadCount)) +
  geom_histogram() + 
  facet_wrap(~transcript_type_cat) +
  scale_x_log10(label = scales::comma_format(),
                breaks = c(1, 10, 100, 1000, 10000, 100000))
  
#' # both at the same time?


counts_cum_2 = counts %>%
  filter(seq_read == 'R1') %>%
  group_by(sample, transcript_type_cat) %>%
  mutate(total_nr_sample = n()) %>%
  group_by(sample, transcript_type_cat, ReadCount, total_nr_sample) %>%
  summarise(total_nr = n(),
            lab_nr = sum(ConversionRate > 0.01)) %>%
  ungroup() %>%
  group_by(sample, transcript_type_cat) %>%
  arrange(desc(ReadCount)) %>%
  mutate(lab_perc_cum = cumsum(lab_nr)/cumsum(total_nr),
         perc_sample_cum = cumsum(total_nr) / total_nr_sample) %>%
  ungroup()

counts_cum_2


counts_cum_2 %>%
  ggplot(aes(ReadCount, lab_perc_cum, color = sample)) +
  #geom_point() +
  geom_line() +
  geom_line(aes(ReadCount, perc_sample_cum, group = sample),
            color = 'grey') +
  scale_x_log10(label = scales::comma_format(),
                breaks = c(1, 3, 5, 10, 100, 1000, 10000, 100000)) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  facet_wrap(~transcript_type_cat) +
  ylab('cum % labeled transcripts (color) OR cum % all transcripts (grey)')


# example for in slide

counts %>%
  filter(transcript_id == "ENST00000454562.1",
         seq_read == 'R1',
         sample %in% c('250_uM_S4U_2', '125_uM_S4U_2', 
                       '0_uM_S4U_2')) %>%
  select(transcript_id, sample, CoverageOnTs, 
         TcReadCount, ConversionRate)

counts %>% 
  filter(seq_read == 'R1',
         sample %in% c('250_uM_S4U_2', '125_uM_S4U_2', 
                       '0_uM_S4U_2')) %>%
  group_by(transcript_id) %>%
  filter(mean(ConversionRate) < 0.01) %>% view()


