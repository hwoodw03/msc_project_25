library("GEOquery")
library("tximport")
library("readr")
library("tximeta")
library("ggplot2")
library("dplyr")
library("SummarizedExperiment")
library("DESeq2")
library("fishpond")
library("tidyr")
library("ggrepel")
library("pheatmap")
library("GenomicFeatures")
library("DTUrtle")
library("rtracklayer")
library("org.Hs.eg.db") 

# code for the DESEQ2/Fishpond differential transcript expression and differentila transcript usage testing
# can get metadata using getGEO (GSE243682_series_matrix.txt.gz) replace with file name and cross match to ACP/control samples with metadata from SRA 

bulk_rna <-getGEO(filename = "FILE_NAME.txt.gz") 
metadata <- pData(bulk_rna)

# replace with directory path 
df <- read_csv("./PATH/TO/RUN/TABLE/SraRunTable.csv", col_select = c(Run, AvgSpotLen, Bases, Experiment,'Sample Name', tissue, tumor_full_name, sex, LibraryLayout, Instrument, LibrarySelection, Platform, Organism))

# get desired information columns from metadata and rename as appropriate
colnames(df)
names(df)[names(df) == "tumor_full_name"] <- "condition"
names(df)[names(df) == "Sample Name"] <- "sample_name"

metadata <- metadata[,-c(3:15, 19:43, 47, 49)]
colnames(metadata)
names(metadata)[names(metadata) == "age of_patient_at_time_of_sample_collection_(yrs):ch1"] <- "age_(yrs)"
names(metadata)[names(metadata) == "Sex:ch1"] <- "sex"
names(metadata)[names(metadata) == "tissue:ch1"] <- "tissue"
names(metadata)[names(metadata) == "tumor full_name:ch1"] <- "condition"
names(metadata)[names(metadata) == "title"] <- "sample_description"
names(metadata)[names(metadata) == "molecule_ch1"] <- "Sample_molecule"
names(metadata)[names(metadata) == "extract_protocol_ch1"] <- "rna_preperation"
names(metadata)[names(metadata) == "extract_protocol_ch1.1"] <- "library_preperation"

metadata$sample_name <- rownames(metadata)

dup_cols <- intersect(names(df), names(metadata))
dup_cols <- setdiff(dup_cols, "sample_name")
metadata <- metadata[, !(names(metadata) %in% dup_cols)]

# merge files to get metadata for selected samples with merged experimental runs
metadata <- merge(df, metadata, by="sample_name")

sra <- metadata %>%
  add_count(Experiment, name = "num_runs")

samples <- sra %>%
  group_by(Experiment) %>%
  dplyr::slice(1) %>%  # keep only one row per sample
  ungroup() %>%
  transmute(
    sample_id = ifelse(num_runs > 1, Experiment, Run),
    condition = condition,
    age = `age_(yrs)`, 
    sample_description = sample_description,
    geo_accession = sample_name,
    tissue = tissue,
    sex = sex
  )

# used as coldata for summarised experiment for fishpond and DeSeq2 
coldata <- samples 

#set wd as salmon output directory used  (2024_salmon_quant_final) and check for existance of the salmon quantification files exist in the sample data table 
coldata$names <- coldata$sample_description
coldata$files <- file.path(dir, "DIRECTORY_NAME", paste0(coldata$sample_id, "_quant"), "quant.sf")
all(file.exists(coldata$files))

# check that condition correct for downstream differential analysis + set levels + check file paths 
coldata$condition
class(coldata$condition)
coldata$condition <- factor(coldata$condition, levels = c("normal brain","adamantinomatous craniopharyngioma"), labels = c("control", "ACP"))
class(coldata$condition)
coldata$condition 
coldata$files

# create summarised experiment object with genes and abundance transcript estimates
se <-tximeta(coldata)
se.exons <- addExons(se)
assayNames(se)
colData(se)

# summarise to gene level too 
gse <- summarizeToGene(se)
head(rownames(gse))

# look at CPLX1
gene_of_interest <- "ENSG00000168993.15"
matches <- which(sapply(rowdata$gene_id, function(x) gene_of_interest %in% x))
rowdata[matches, ]
transcript_ids <- rowData(se)$tx_name[matches]
transcript_ids
head(rownames(se))

# all the transcripts of CPLX1 abundance (TPM) in case for visualisation 
expr_abundance <- assay(se, "abundance")[matches, ]
expr_df <- as.data.frame(expr_abundance)
rownames(expr_df) <- transcript_ids

# get CPLX1-201 and CPLX1-203 abundances 
tpm_mat <- assay(se, "abundance")[c("ENST00000304062.11","ENST00000505203.1"), ]
cplx1_names <- c("ENST00000304062.11" = "CPLX1-201","ENST00000505203.1"  = "CPLX1-203")

# some useful plotting with heatmap/expression plot 

tpm_df <- as.data.frame(tpm_mat)
tpm_df$transcript <- rownames(tpm_df)
tpm_long <-tpm_df %>%
  pivot_longer(-transcript, names_to = "sample", values_to = "TPM")
tpm_long <- tpm_long %>%
  left_join(coldata, by = c("sample" = "sample_description"))
tpm_long$transcript_name <- cplx1_names[tpm_long$transcript]
tpm_df$transcript_name <- cplx1_names[tpm_df$transcript]

sample_annotation <- data.frame(condition = coldata$condition)
rownames(sample_annotation) <- coldata$sample_description  # sample names as rownames to match matrix colnames
rownames(tpm_mat) <- cplx1_names[rownames(tpm_mat)]

heatmap <- pheatmap(tpm_mat,
         scale = "row",  # standardize transcripts to mean=0, sd=1 for visualization
         annotation_col = sample_annotation,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Heatmap of CPLX1 Isoform Expression")

isoform_expression_plot <- ggplot(tpm_summary, aes(x = transcript_name, y = mean_logTPM, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_logTPM - se_logTPM, ymax = mean_logTPM + se_logTPM),
                width = 0.2,
                position = position_dodge(width = 0.8)) +
  geom_point(data = tpm_long %>% mutate(logTPM = log2(TPM + 1)),
             aes(x = transcript_name, y = logTPM, color = condition),
             position = position_dodge(width = 0.8),
             size = 2, alpha = 0.7) +
  geom_text_repel(data = tpm_long %>% mutate(logTPM = log2(TPM + 1)),
                  aes(x = transcript_name, y = logTPM, label = sample),
                  color = "black",
                  position = position_dodge(width = 0.8),
                  size = 3,
                  segment.size = 0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Isoform Expression") +
  xlab("Isoform") +
  ggtitle("Mean Expression of CPLX1 Isoforms")

isoform_expression_plot

# DESeQ2 for DTE as se has transcript counts/abundances not summarised to gene level (held in gse)

dds <- DESeqDataSet(se, design = ~ condition)
keep <- rowSums(counts(dds) >= 5) >= 3
dds <- dds[keep, ]
dds <- estimateSizeFactors(dds)

# transform counts to visualise data 
rld <- rlog(dds, blind = FALSE)

levels(dds$condition)

pca_plot <- plotPCA(rld, intgroup = c("condition"))
dds <- DESeq(dds)

# inspect result and can choose what to focus on
res <- results(dds)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)

# focusing on CPLX1
res_df <- as.data.frame(res)
res_df$logPadj <- -log10(res_df$padj)
head(res_df)

# downloaded gencode annotation file and create mapping files
gtf_file <- "PATH/TO/FILE/gencode.v47.annotation.gtf"
gtf <- import(gtf_file)
tx_meta <- gtf[gtf$type == "transcript"]
tx2gene <- data.frame(
  transcript_ids = mcols(tx_meta)$transcript_id,
  gene_ids = mcols(tx_meta)$gene_id,
  gene_name = mcols(tx_meta)$gene_name
)

tx2genemap
gene_map <- unique(tx2genemap[, c("gene_ids", "gene_name")])
gene_map

txmap <- data.frame(
  transcript_ids = mcols(tx_meta)$transcript_id,
  transcript_name = mcols(tx_meta)$transcript_name
)

# merge to get transcript names
res_df$transcript_ids <- rownames(res_df)
res_df <- merge(res_df, txmap, by = "transcript_ids", all.x = TRUE)

# add to dds 
rowData(dds)
matched_names <- res_df$transcript_name[match(rownames(dds), res_df$transcript_ids)]
rowData(dds)$transcript_name <- matched_names
head(rowData(dds))

# define transcript of interest and plot differential transcript expression of what you want

transcripts_of_interest <- c("CPLX1-201", "CPLX1-203")
cplx1_tx_ids <- rownames(dds)[rowData(dds)$transcript_name %in% tx_names]
transcripts_of_interest <- c("CPLX1-201", "CPLX1-203")

# Filter the results data frame
cplx1_df <- res_df %>%
  filter(transcript_name %in% transcripts_of_interest)

Deseq2_dte_plot <- ggplot(cplx1_df, aes(x = transcript_name, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(
    aes(label = paste0("p=", signif(padj, 3))),
    vjust = ifelse(cplx1_df$log2FoldChange > 0, -0.8, 1.5),  # moves label above/below bar
    size = 3,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("lightblue", "salmon"),
    labels = c("Decreased", "Increased")
  ) +
  labs(
    title = "Differential Expression of CPLX1 Isoforms in ACP vs Control",
    x = "Isoform",
    y = "Log2 Fold Change",
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none") 
Deseq2_dte_plot

#################################################################

# fishpond analysis for DTE/DTU 

# plot inferential replicates to see uncertainty in abundance estimation 

plotInfReps(y_filtered_swish, 
            idx = "ENST00000304062.11", 
            x = "condition", 
            main = "CPLX1-201")

plotInfReps(y_filtered_swish, 
            idx = "ENST00000505203.1", 
            x = "condition", 
            main = "CPLX1-203")

txi_counts <- assays(y_filtered_swish)$counts 
summary(rowMeans(txi_counts))
hist(log10(rowMeans(txi_counts) + 1), main = "Transcript Abundance Distribution")

y <- se 
assayNames(y)
rowRanges(y)
gy <- gse 
rowdata <- rowData(se)

y_scale <- scaleInfReps(y)
y_filter_label <- labelKeep(y_scale)
y_filtered <- y_filter_label[mcols(y_filter_label)$keep,]

# set seed 
set.seed(3)
y_filtered_swish <- swish(y_filtered, x="condition")

rowData(y_filtered_swish)
rownames(y_filtered_swish) <- rowData(y_filtered_swish)$gene_id
head(rownames(y_filtered_swish))

# add gene names 
y_filtered_swish <- addIds(y_filtered_swish, "SYMBOL", gene = TRUE)  

# get cplx1 results with direction of change
rownames(y_filtered_swish) <- y_filtered_swish$tx_name
dte_results <- as.data.frame(mcols(y_filtered_swish))
dte_results <- dte_results %>%
  mutate(direction = case_when(
    log2FC > 0 ~ "Increase",
    log2FC < 0 ~ "Decrease",
    TRUE ~ "No change"
  ))
dte_cplx1_result <- dte_results %>% filter(SYMBOL == "CPLX1")
dte_cplx1_result$transcript_name <- cplx1_names[dte_cplx1_result$tx_name]
dte_cplx1_result

# can repeat for gene-level analysis and plotting 

gy_scale <- scaleInfReps(gy)
gy_filter_label <- labelKeep(gy_scale)
gy_filtered <- gy_filter_label[mcols(gy_filter_label)$keep,]
gy_filtered <- addIds(gy_filtered, "SYMBOL", gene = TRUE)

# set genes of interest that can filter by  
genes_of_interest <- c("CPLX1", "CPLX2", "CPLX3", "CPLX4")
gy_filtered_CPLX <- gy_filtered[mcols(gy_filtered)$SYMBOL %in% genes_of_interest, ]

tpm_mat_gene <- assay(gy_filtered_CPLX, "abundance")
tpm_mat_gene_df <- as.data.frame(tpm_mat_gene)
tpm_mat_gene_df$gene_name <- rownames(tpm_mat_gene_df)


tpm_gene_long <-tpm_mat_gene_df %>%
  pivot_longer(-gene_name, names_to = "sample", values_to = "TPM")

tpm_gene_long <- tpm_gene_long %>%
  left_join(coldata, by = c("sample" = "sample_description"))

CPLX_gene_names <- c(
  "ENSG00000168993.15" = "CPLX1",
  "ENSG00000145920.15"  = "CPLX2",
  "ENSG00000213578.6" = "CPLX3"
)

tpm_gene_long$gene_name_abr <- CPLX_gene_names[tpm_gene_long$gene_name]
tpm_mat_gene_df$gene_name_abr <- CPLX_gene_names[tpm_mat_gene_df$gene_name]
tpm_gene_summary <- tpm_gene_long %>%
  mutate(logTPM = log2(TPM + 1)) %>%
  group_by(gene_name_abr, condition) %>%
  summarise(mean_logTPM = mean(logTPM, na.rm = TRUE),
            sd_logTPM = sd(logTPM, na.rm = TRUE),
            n = n(),
            se_logTPM = sd_logTPM / sqrt(n)) %>%
  ungroup()

# plot 
ggplot(tpm_gene_summary, aes(x = gene_name_abr, y = mean_logTPM, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_logTPM - se_logTPM, ymax = mean_logTPM + se_logTPM),
                width = 0.2,
                position = position_dodge(width = 0.8)) + 
  scale_fill_manual(values = c("control" = "salmon", "ACP" = "lightblue")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Isoform Expression") +
  xlab("Gene") +
  ggtitle("Mean Expression of Complexin Genes in ACP vs Control")

# DGE 
set.seed(12)
gy_filtered_swish <- swish(gy_filtered, x="condition")
gy_filtered_swish <- addIds(gy_filtered_swish, "SYMBOL", gene = TRUE)

rownames(gy_filtered_swish) <- gy_filtered_swish$gene_name
dge_results <- as.data.frame(mcols(gy_filtered_swish))

dge_results <- dge_results %>%
  mutate(direction = case_when(
    log2FC > 0 ~ "Increase",
    log2FC < 0 ~ "Decrease",
    TRUE ~ "No change"
  ))

# can choose interactors of interest and used ones closest to on STRING 
interactors_of_interest <-c("CPLX1","STX1A", "SNAP25", "SYT1", "VAMP2", "VAMP3", "SYT2", "STX19", "STX1B","CPLX2", "STX4","SNAP23")
dge_interactors_result <- dge_results %>% filter(SYMBOL %in% interactors_of_interest)
# results of DGE 
dge_interactors_result 

# can plot 
ggplot(dge_interactors_result, aes(x = SYMBOL)) +
  geom_segment(aes(xend = SYMBOL, y = 0, yend = log2FC, color = direction),
               arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  geom_point(aes(y = log2FC, color = direction), size = 4) +
  geom_text(aes(y = log2FC, label = paste0("q=", signif(qvalue, 3))),
            nudge_y = ifelse(dge_interactors_result$log2FC > 0, 0.4, -0.4),
            vjust = 0.5, size = 3, color = "black") +
  scale_color_manual(values = c("Decrease" = "lightblue", "Increase" = "salmon")) +
  labs(title = "Differential Expression of CPLX1 and SNARE Interactors",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

#################################################
iso<- isoformProportions(y_filtered, geneCol = "gene_id", quiet = FALSE)
set.seed(19)
iso_swish <- swish(iso_whole, x="condition") 

# can look at sig transcripts 
sig_transcripts <- iso_swish[which(iso_swish$qvalue < 0.1), ]

# can redefine transcripts of interest if interested in other ones 
cplx1_iso <- iso_swish[mcols(iso_swish)$transcript_name %in% transcripts_of_interest, ]

iso_prop_mat_cplx1 <- assays(cplx1_iso_swish)[["isoProp"]]

res_dtu <- mcols(iso)
cplx1_iso <- res_dtu[rownames(res_dtu) %in% transcripts_of_interest, ]

# get corresponding samples and transcript names 
metadata_df <- as.data.frame(mcols(cplx1_iso_swish_whole))
colnames(txmap) <- c("transcript_id", "transcript_name")
metadata_df$transcript_id <- rownames(metadata_df)
metadata_df <- metadata_df %>%
  left_join(txmap, by = "transcript_id")

df_iso_prop_whole_cplx1 <- as.data.frame(iso_prop_mat_cplx1)
df_iso_prop_whole_cplx1$transcript_id <- rownames(df_iso_prop_whole_cplx1)

df_long_whole_cplx1 <- df_iso_prop_whole_cplx1 %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "sample",
    values_to = "proportion"
  )

df_long_whole_cplx1 <- df_long_whole_cplx1 %>%
  left_join(metadata_df, by = "transcript_id")

sample_data <- as.data.frame(colData(cplx1_iso_swish_whole))
sample_data$sample <- rownames(sample_data)

df_long_whole_cplx1 <- df_long_whole_cplx1 %>%
  left_join(sample_data[, c("sample", "condition")], by = "sample")

# summarise cplx1 results into dataframe with the q value and conditions 

df_whole_cplx1_summary <- df_long_whole_cplx1 %>%
  group_by(transcript_id, transcript_name, condition, pvalue, qvalue) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop")
