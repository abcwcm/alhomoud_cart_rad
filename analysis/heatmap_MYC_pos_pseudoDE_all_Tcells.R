# heatmap for CART rad paper

psedu_DE_res_myc <- read.csv("/SSD/maa7095/scRNAseq/Guzman_CART/scanpy_process/clean_scripts/CART_myc_pos_pseudo_DE.csv", row.names = "X")
log_mat_myc <- read.csv("/SSD/maa7095/scRNAseq/Guzman_CART/scanpy_process/clean_scripts/CART_myc_pos_cells_pseudobulk_log_expr.csv", row.names = "X")


sig_datafram <- psedu_DE_res_myc %>% filter(padj < 0.05)

select_genes_moh <- c("Arid1a", "Ubash3b", "Sdc4", "Sell", "Cd72", "Ndfip2", "Sox4", "Kdm3a", "Areg", "Dapl1", 
                      "Pdcd6", "Pycard", 
                      "Coro1a", "Ciapin1", "Bst2", "Oaz1",
                      "Il12rb2", "Ifi202b", "Mapkapk2", "Atp5c1", "Lat", "Ifitm2", "Tox2", "Bax",
                      "Arap2", "Arhgdib")

#subset mat for genes of interest.
order_samples <- c("out_Mouse_23_T_cell", "out_Mouse_4_T_cell", "out_Mouse_8_T_cell", "out_Mouse_11_T_cell", "out_Mouse_13_T_cell", "out_Mouse_5_T_cell")
sub_log_map <- log_mat_myc[select_genes_moh,order_samples]
#names(sub_log_map) <- c("Cy/mCART19_rep1","Cy/mCART19_rep2","Cy/mCART19_rep3","Cy/LD-TBI/mCART19_rep1","Cy/LD-TBI/mCART19_rep2","Cy/LD-TBI/mCART19_rep3")
names(sub_log_map) <- c("1","2","3","1","2","3")
library(pheatmap)
library(dplyr)

annotation_col <- data.frame(
  Condition = factor(
    c("Cy/mCART19", "Cy/mCART19", "Cy/mCART19", 
      "Cy/LD-TBI/mCART19", "Cy/LD-TBI/mCART19", "Cy/LD-TBI/mCART19"),
    levels = c("Cy/mCART19", "Cy/LD-TBI/mCART19")  # Optional: Specify levels if needed
  )
)

#rownames(annotation_col) <- colnames(sub_log_map)  # Set column names as rownames for annotation


annotation_row <- data.frame(
  Gene_Annotation = factor(c("Inhibition and\nregulation", "Inhibition and\nregulation", "Inhibition and\nregulation", "Inhibition and\nregulation", "Inhibition and\nregulation", "Inhibition and\nregulation", "Inhibition and\nregulation", "Inhibition and\nregulation", "Inhibition and\nregulation", "Inhibition and\nregulation", 
                       "Apoptosis", "Apoptosis", 
                       "Survival and\nproliferation", "Survival and\nproliferation", "Survival and\nproliferation","Survival and\nproliferation",
                       "Activation", "Activation", "Activation","Activation","Activation","Activation","Activation","Activation",
                       "Adhesion", "Adhesion"), 
                     levels = c("Inhibition and\nregulation","Apoptosis","Survival and\nproliferation","Activation","Adhesion"))
)


# Adjust row names to match the new number of rows in annotation_row
#rownames(annotation_row) <- rownames(sub_log_map)  # Ensure rownames match the matrix



z_scores_matrix <- t(scale(t(sub_log_map))) # Scale each row to have mean 0 and sd 1

pdf("myc_pos_pseudoDE_res_publ.pdf", width = 8, height = 10)  # Specify desired width and height for the PDF


Heatmap(
  z_scores_matrix, 
  name = "z-score", 
  show_row_names = TRUE, 
  show_column_names = TRUE, 
  cluster_rows = FALSE,  # No row clustering
  cluster_columns = FALSE,  # No column clustering
  row_split = annotation_row$Gene_Annotation,  # Row split based on gene type annotations
  column_split = annotation_col$Condition,
  top_annotation = HeatmapAnnotation(df = annotation_col,
                                     show_legend = FALSE),  # Column anno??heatmap_legend_sidetations
  left_annotation = rowAnnotation(df = annotation_row,
                                  show_legend = FALSE),  # Row annotations (fix here)
  col = colorRampPalette(c("blue", "white", "red"))(100),  # Color scale
  #column_names_gp = gpar(rot = 90),  # Rotate column names by 90 degrees
  show_heatmap_legend = T, 
)

dev.off()
