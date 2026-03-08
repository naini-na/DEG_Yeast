# Analisis Ekspresi Gen Fase Lipogenik pada Rhodotorula toruloides untuk
# Mengidentifikasi Regulasi Biosintesis Very Long-Chain Fatty Acids (VLCFAs)

# Dataset: GSE241709 (NA04_120h NA04_48h NA14_120h NA14_48h WT_120h WT_48h)  
# Platform: GPL32811	(Illumina NovaSeq 6000 (Rhodotorula toruloides))
# Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) 

# ===============================
# A. PEMANGGILAN LIBRARY 
# ===============================
# library() digunakan agar fungsi di dalam package bisa digunakan  
library(GEOquery) 
library(limma) 
library(pheatmap) 
library(ggplot2) 
library(dplyr) 
library(illuminaHumanv4.db) 
library(AnnotationDbi) 
library(umap)

# ===============================
# B. PENGAMBILAN DATA DARI GEO 
# ===============================
# GEO (Gene Expression Omnibus) adalah database publik milik NCBI 
# getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO 
# GSEMatrix = TRUE -> data diambil dalam format ExpressionSet 
# AnnotGPL  = TRUE -> anotasi gen (Gene Symbol) ikut diunduh 

gset <- getGEO("GSE241709", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#Karena anotasi gen tidak tersedia dalam matrix maka mengambil dari
#supplementary files
getGEOSuppFiles("GSE241709")
list.files("GSE241709")

#Supplementary Files:
#"GSE241709_gene_FPKM.txt.gz" "GSE241709_raw_count.txt.gz"
#Data yang dilanjutkan untuk analisis DEG adalah GSE241709_raw_count.txt.gz
counts <- read.delim("GSE241709/GSE241709_raw_count.txt.gz",
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

dim(counts)
head(counts)
colnames(counts)

#Memisahkan count matrix dari anotasi
count_matrix <- counts[, 1:18]

dim(count_matrix)
head(count_matrix)
head(rownames(count_matrix))

# ===============================
# C. PRE-PROCESSING DATA EKSPRESI  
# ===============================
# Mengubah count_matrix menjadi matriks numerik
# Baris = gen/probe
# Kolom = sampel 
ex <- as.matrix(count_matrix)

#quantile(): menghitung nilai kuantil (persentil) 
#as.numeric(): mengubah hasil quantile (yang berupa named vector) 
#menjadi vektor numerik biasa agar mudah dibandingkan 
qx <- as.numeric(quantile(ex, 
                          c(0, 0.25, 0.5, 0.75, 0.99, 1), 
                          na.rm = TRUE))

if (qx[5] > 100 || (qx[6] - qx[1]) > 50) {
  print("Data likely needs log transformation")
}

#LogTransform adalah variabel logika (TRUE / FALSE) 
#Operator logika: 
#>  : lebih besar dari 
#|| : OR (atau) 
#&& : AND (dan) 
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#IF statement: 
#Jika LogTransform = TRUE, maka lakukan log2 
if (LogTransform) { 
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA 
  ex[ex <= 0] <- NA 
  ex <- log2(ex) 
} 

# ===============================
# D. DEFINISI KELOMPOK SAMPEL 
# ===============================
# 
strain <- factor(c(rep("NA14",6),
                   rep("WT",6),
                   rep("NA04",6)))

time <- factor(rep(c(rep("48h",3),
                     rep("120h",3)),3))

group <- factor(paste(strain, time, sep="_"))

print(group)

# ===============================
# E. DESIGN MATRIX (KERANGKA STATISTIK)
# ===============================

# Buat metadata
count_matrix <- counts[, 1:18]

# Buat faktor strain
strain <- factor(c(rep("NA14",6),
                   rep("WT",6),
                   rep("NA04",6)))

# Buat faktor time
time <- factor(rep(c(rep("48h",3),
                     rep("120h",3)),3))

# Gabungkan jadi metadata
coldata <- data.frame(
  row.names = colnames(count_matrix),
  strain = strain,
  time = time
)

coldata

#Menentukan baseline
coldata$strain <- relevel(coldata$strain, ref = "WT")
coldata$time <- relevel(coldata$time, ref = "48h")

#model.matrix(): 
#Membuat matriks desain untuk model linear 
design <- model.matrix(~ strain + time, data = coldata)
design 

colnames(design) <- make.names(colnames(design))
colnames(design)

# ===============================
# F. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
# ===============================

# 1. Install package edgeR
#edgeR adalah package Bioconductor di R yang digunakan untuk analisis 
#DEG dari data RNA-seq berbasis count
BiocManager::install("edgeR")
library(edgeR)

# Buat DGEList dari raw count
dge <- DGEList(counts = count_matrix)

# Filter low expression
keep <- filterByExpr(dge, group = coldata$strain)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalisasi TMM
dge <- calcNormFactors(dge)

# Voom transform
v <- voom(dge, design, plot=TRUE)

# Linear model
fit <- lmFit(v, design)

# Buat kontras pairwise
# Membandingkan strain dan waktu
contrast_matrix <- makeContrasts(
  
  NA14_vs_WT  = strainNA14,
  NA04_vs_WT  = strainNA04,
  NA14_vs_NA04 = strainNA14 - strainNA04,
  Time120_vs_48 = time120h,
  
  levels = design
)


# Terapkan kontras
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

colnames(fit2$coefficients)

# Ambil hasil DEG
# NA14 vs WT
res_NA14_WT <- topTable(
  fit2,
  coef = "NA14_vs_WT",
  adjust = "fdr",
  number = Inf,
  p.value = 0.05
)

# NA04 vs WT
res_NA04_WT <- topTable(
  fit2,
  coef = "NA04_vs_WT",
  adjust = "fdr",
  number = Inf,
  p.value = 0.05
)

# NA14 vs NA04
res_NA14_NA04 <- topTable(
  fit2,
  coef = "NA14_vs_NA04",
  adjust = "fdr",
  number = Inf,
  p.value = 0.05
)

# Time 120h vs 48h
res_120h_48h <- topTable(
  fit2,
  coef = "Time120_vs_48",
  adjust = "fdr",
  number = Inf,
  p.value = 0.05
)


# Cek jumlah DEG
nrow(res_NA14_WT)
nrow(res_NA04_WT)
nrow(res_NA14_NA04)
nrow(res_120h_48h)
 
# Simpan hasil DEG
write.csv(res_NA14_WT, "DEG_NA14_vs_WT.csv")
write.csv(res_NA04_WT, "DEG_NA04_vs_WT.csv")
write.csv(res_NA14_NA04, "DEG_NA14_vs_NA04.csv")
write.csv(res_120h_48h, "DEG_120h_vs_48h.csv")

# ===============================
# G BOXPLOT DISTRIBUSI NILAI EKSPRESI  
# ===============================

# Data ekspresi untuk visualiasi
expr_voom <- v$E

# Berdasarkan strain
group_colors <- as.numeric(coldata$strain)

boxplot(
  expr_voom,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Distribusi Nilai Ekspresi (Log2-CPM) per Sampel",
  ylab = "Log2 CPM"
)

legend(
  "topright",
  legend = levels(coldata$strain),
  fill = 1:length(levels(coldata$strain)),
  cex = 0.8
)

#Berdasarkan strain dan time
group_combined <- interaction(coldata$strain, coldata$time)
group_colors <- as.numeric(group_combined)

boxplot(
  expr_voom,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Distribusi Log2-CPM berdasarkan Strain dan Waktu",
  ylab = "Log2 CPM"
)

legend(
  "topright",
  legend = levels(group_combined),
  fill = 1:length(levels(group_combined)),
  cex = 0.7
)

# ===============================
# H. DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)
# ===============================

# Warna berdasarkan strain + time
group_combined <- interaction(coldata$strain, coldata$time)
group_colors <- as.numeric(group_combined)

# Plot density sampel pertama dulu
plot(
  density(expr_voom[,1]),
  col = group_colors[1],
  lwd = 2,
  main = "Density Plot Distribusi Log2-CPM",
  xlab = "Log2 CPM"
)

# Tambahkan sampel lain
for (i in 2:ncol(expr_voom)) {
  lines(
    density(expr_voom[,i]),
    col = group_colors[i],
    lwd = 2
  )
}

legend(
  "topright",
  legend = levels(group_combined),
  col = 1:length(levels(group_combined)),
  lwd = 2,
  cex = 0.7
)

# ===============================
# I. UMAP (VISUALISASI DIMENSI RENDAH) 
# ===============================

# Jalankan UMAP
# Transpose (samples jadi baris)
expr_t <- t(expr_voom)

# Jalankan UMAP
umap_result <- umap(expr_t)

# Ambil koordinat
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Tambahkan metadata
umap_df$strain <- coldata$strain
umap_df$time <- coldata$time

# Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, 
                    color = strain, shape = time)) +
  geom_point(size = 4) +
  theme_minimal() +
  ggtitle("UMAP Plot - Transcriptome Profile") +
  theme(text = element_text(size = 14))

# ===============================
# J. VISUALISASI VOLCANO PLOT 
# ===============================

make_volcano <- function(res_table, title_plot) {
  
  volcano_df <- res_table
  
  volcano_df$Significant <- "Not Significant"
  
  volcano_df$Significant[
    volcano_df$adj.P.Val < 0.05 & volcano_df$logFC > 1
  ] <- "Upregulated"
  
  volcano_df$Significant[
    volcano_df$adj.P.Val < 0.05 & volcano_df$logFC < -1
  ] <- "Downregulated"
  
  ggplot(volcano_df,
         aes(x = logFC,
             y = -log10(adj.P.Val),
             color = Significant)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(values = c(
      "Downregulated" = "blue",
      "Not Significant" = "grey70",
      "Upregulated" = "red"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    ggtitle(title_plot) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value") +
    theme(text = element_text(size = 14))
}

# NA14 vs WT
volcano_NA14_WT <- make_volcano(
res_NA14_WT,
"Volcano Plot: NA14 vs WT"
)

volcano_NA14_WT

# NA04 vs WT
volcano_NA04_WT <- make_volcano(
  res_NA04_WT,
  "Volcano Plot: NA04 vs WT"
)

volcano_NA04_WT

# NA14 vs NA04
volcano_NA14_NA04 <- make_volcano(
  res_NA14_NA04,
  "Volcano Plot: NA14 vs NA04"
)

volcano_NA14_NA04


# 120h vs 48h
volcano_120h_48h <- make_volcano(
  res_120h_48h,
  "Volcano Plot: 120h vs 48h"
)

volcano_120h_48h

# ===============================
# K. VISUALISASI HEATMAP 
# ===============================

make_heatmap <- function(res_table, expr_matrix, title_plot) {
  
  # Ambil top 50 berdasarkan adj.P.Val
  top50 <- res_table[order(res_table$adj.P.Val), ]
  top50 <- top50[1:50, ]
  
  genes_top50 <- rownames(top50)
  
  # Subset ekspresi
  heat_data <- expr_matrix[genes_top50, ]
  
  # Scaling per gen (z-score)
  heat_data_scaled <- t(scale(t(heat_data)))
  
  # Annotation kolom
  annotation_col <- data.frame(
    Strain = coldata$strain,
    Time = coldata$time
  )
  
  rownames(annotation_col) <- colnames(expr_matrix)
  
  # Plot heatmap
  pheatmap(
    heat_data_scaled,
    annotation_col = annotation_col,
    show_rownames = FALSE,
    fontsize_col = 10,
    main = title_plot,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete"
  )
}

# NA14 vs WT
make_heatmap(
  res_NA14_WT,
  expr_voom,
  "Heatmap Top 50 DEG: NA14 vs WT"
)

# NA04 vs WT
make_heatmap(
  res_NA04_WT,
  expr_voom,
  "Heatmap Top 50 DEG: NA04 vs WT"
)

# NA14 vs NA04
make_heatmap(
  res_NA14_NA04,
  expr_voom,
  "Heatmap Top 50 DEG: NA14 vs NA04"
)

# 120h vs 48h
make_heatmap(
  res_120h_48h,
  expr_voom,
  "Heatmap Top 50 DEG: 120h vs 48h"
)

# ===============================
# L. MENYIMPAN HASIL  
# ===============================

# write.csv(): menyimpan hasil analisis ke file CSV 
write.csv(topTableResults, "Hasil_GSE241709_DEG.csv") 

message("Analisis selesai. File hasil telah disimpan.") 

# ===============================
# M. GO ENRICHMENT & KEGG PATHWAY ENRICHMENT
# ===============================

# Cek database untuk Rhodosporidium toruloides
search_kegg_organism("toruloides", by="scientific_name")

# Tidak bisa dilakukan enrichment di R, karena
# - Spesies tidak punya OrgDb
# - Tidak punya annotation GO/KEGG
# - Tidak ada file ortholog/mapping

# GO enrichment dan KEGG pathway enrichment akan dilakukan 
# menggunakan IDEP untuk data RNA-Seq

# Menyiapkan tabel anotasi gen
gene_annot <- data.frame(gene_name = counts$gene_name)
gene_annot$gene_id <- rownames(counts)

head(gene_annot)

# Menambahkan ID gen ke DEG time
# Analisis enrichment difokuskan pada perbedaan waktu
# untuk melihat profil transcriptome pada fase lipogenik
res_120h_48h$gene_id <- rownames(res_120h_48h)

# Mapping anotasi ke DEG waktu (time)
res_time_annot <- merge(
  res_120h_48h,
  gene_annot,
  by = "gene_id",
  all.x = TRUE
)

head(res_time_annot)

# Mengambil gen signifikan
deg_time <- res_time_annot[
  res_time_annot$adj.P.Val < 0.05,
]

# Menyimpan data untuk analisis iDEP
write.table(
  deg_time$gene_name,
  "DEG_time_120h_vs_48h.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# karena iDEP membutuhkan data fold changes dan adjust p value, maka perlu
# dibuat file baru berisi gene name, fold changes, dan adjust p value
idep_input <- res_time_annot[, c("gene_name", "logFC", "adj.P.Val")]

colnames(idep_input) <- c("gene_id", "logFC", "adj.P.Val")

write.table(
  idep_input,
  file = "iDEP_time_DEG.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Gene ID tidak terbaca oleh iDEP, maka perlu dikonversi menjadi UniProt
# menggunakan Uniprot ID Mapping (https://www.uniprot.org/id-mapping)

# Menggabungkan ID Uniprot ke DEG
map <- read.delim("idmapping_new.tsv")

res_time_uniprot2 <- merge(
  res_time_annot,
  map,
  by.x = "gene_name",
  by.y = "From",
  all.x = TRUE
)

# Menyimpan hasil 
write.table(
  res_time_uniprot2,
  "for_iDEP.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

idep_input <- res_time_uniprot2[, c("Entry","logFC","adj.P.Val")]
colnames(idep_input) <- c("gene_id","logFC","adj.P.Val")

write.table(
  idep_input,
  "for_iDEP.txt",
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

# Hasil analisis enrichment menggunakan iDEP "GAGAL" karena ID tidak dikenali 
# atau tidak ada dalam database iDEP