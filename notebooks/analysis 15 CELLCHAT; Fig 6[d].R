# 15 CELLCHAT TBI

renv::status()
renv::snapshot()
renv::restore()

install.packages("renv", force=TRUE)
install.packages("devtools")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
BiocManager::install("Biobase")
devtools::install_github("sqjin/CellChat")
devtools::install_github("thomasp85/patchwork")

library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)
library(reticulate)
library(patchwork)
library(grid)
library(renv)
library("ggplotify")
library(cowplot)
library(dichromat)
library(figpatch)
library(ggpubr)
library(RColorBrewer)
options(stringsAsFactors = FALSE)

#setwd("~/Desktop/My projects/[AK002] Thymus aging & regeneration/AK002_R_project/")

# load cellchat objects and rename subsets per dataset (light version) (only secreted)
cellchat.02moD0 <- readRDS("../output/cellchat_io/cd45neg_0147_SLTBI_vdb_02mo_D0_NEW_cellchatDB.rds")
cellchat.02moD1 <- readRDS("../output/cellchat_io/cd45neg_0147_SLTBI_vdb_02mo_D1_NEW_cellchatDB.rds")
cellchat.02moD4 <- readRDS("../output/cellchat_io/cd45neg_0147_SLTBI_vdb_02mo_D4_NEW_cellchatDB.rds")
cellchat.02moD7 <- readRDS("../output/cellchat_io/cd45neg_0147_SLTBI_vdb_02mo_D7_NEW_cellchatDB.rds")
cellchat.18moD0 <- readRDS("../output/cellchat_io/cd45neg_0147_SLTBI_vdb_18mo_D0_NEW_cellchatDB.rds")
cellchat.18moD1 <- readRDS("../output/cellchat_io/cd45neg_0147_SLTBI_vdb_18mo_D1_NEW_cellchatDB.rds")
cellchat.18moD4 <- readRDS("../output/cellchat_io/cd45neg_0147_SLTBI_vdb_18mo_D4_NEW_cellchatDB.rds")
cellchat.18moD7 <- readRDS("../output/cellchat_io/cd45neg_0147_SLTBI_vdb_18mo_D7_NEW_cellchatDB.rds")
levels(cellchat.18moD7@meta$cell_type_subset) 

levels(cellchat.02moD0@meta$cell_type_subset) <- c( "arEC", "capEC",  "venEC", 'capsFB', 'intFB', 'medFB',  "MEC", "vSMCPC", "nmSC", "Fat", "aaTEC1", "aaTEC2",  'cTEC', 'earlyPr', 'mTEC1', 'mTEC prol', 'mTEC2', 'basal', 'tuft', 'neuroendo', 'goblet', 'microfold')
cellchat.02moD0 <- setIdent(cellchat.02moD0, ident.use = "cell_type_subset") # set "cell_type" as default cell identity
levels(cellchat.02moD1@meta$cell_type_subset) <- c( "arEC", "capEC",  "venEC", 'capsFB', 'intFB', 'medFB',  "MEC", "vSMCPC", "nmSC", "Fat", "aaTEC1", "aaTEC2",  'cTEC',  'earlyPr', 'mTEC1', 'mTEC prol', 'mTEC2', 'basal', 'tuft', 'neuroendo', 'goblet', 'microfold')
cellchat.02moD1 <- setIdent(cellchat.02moD1, ident.use = "cell_type_subset") # set "cell_type" as default cell identity
levels(cellchat.02moD4@meta$cell_type_subset) <- c( "arEC", "capEC",  "venEC", 'capsFB', 'intFB', 'medFB',   "MEC", "vSMCPC", "nmSC", "Fat", "aaTEC1", "aaTEC2",  'cTEC', 'earlyPr', 'mTEC1', 'mTEC prol', 'mTEC2', 'basal', 'tuft', 'neuroendo', 'goblet', 'microfold')
cellchat.02moD4 <- setIdent(cellchat.02moD4, ident.use = "cell_type_subset") # set "cell_type" as default cell identity
levels(cellchat.02moD7@meta$cell_type_subset) <- c( "arEC", "capEC",  "venEC", 'capsFB', 'intFB', 'medFB',   "MEC", "vSMCPC", "nmSC", "Fat", "aaTEC1", "aaTEC2",  'cTEC', 'earlyPr' , 'mTEC1', 'mTEC prol', 'mTEC2', 'basal', 'tuft', 'neuroendo', 'goblet', 'microfold')
cellchat.02moD7 <- setIdent(cellchat.02moD7, ident.use = "cell_type_subset") # set "cell_type" as default cell identity
levels(cellchat.18moD0@meta$cell_type_subset) <- c( "arEC", "capEC",  "venEC", 'capsFB', 'intFB', 'medFB', "MEC", "vSMCPC", "nmSC", "Fat", "aaTEC1", "aaTEC2",  'cTEC',  'earlyPr' , 'mTEC1', 'mTEC prol', 'mTEC2', 'basal', 'tuft', 'neuroendo', 'goblet', 'microfold')
cellchat.18moD0 <- setIdent(cellchat.18moD0, ident.use = "cell_type_subset") # set "cell_type" as default cell identity
levels(cellchat.18moD1@meta$cell_type_subset) <-c( "arEC", "capEC",  "venEC", 'capsFB', 'intFB', 'medFB', "MEC", "vSMCPC", "nmSC", "Fat", "aaTEC1", "aaTEC2",  'cTEC', 'earlyPr' ,  'mTEC1', 'mTEC prol', 'mTEC2', 'basal', 'tuft', 'neuroendo', 'goblet', 'microfold')
cellchat.18moD1 <- setIdent(cellchat.18moD1, ident.use = "cell_type_subset") # set "cell_type" as default cell identity
levels(cellchat.18moD4@meta$cell_type_subset) <-c( "arEC", "capEC",  "venEC", 'capsFB', 'intFB', 'medFB', "MEC", "vSMCPC", "nmSC", "Fat", "aaTEC1", "aaTEC2",  'cTEC','earlyPr', 'mTEC1', 'mTEC prol', 'mTEC2',  'basal', 'tuft', 'neuroendo', 'goblet', 'microfold')
cellchat.18moD4 <- setIdent(cellchat.18moD4, ident.use = "cell_type_subset") # set "cell_type" as default cell identity
levels(cellchat.18moD7@meta$cell_type_subset) <- c( "arEC", "capEC",  "venEC", 'capsFB', 'intFB', 'medFB', "MEC", "vSMCPC", "nmSC", "Fat", "aaTEC1", "aaTEC2",  'cTEC', 'earlyPr', 'mTEC1', 'mTEC prol', 'mTEC2', 'basal','tuft', 'neuroendo', 'goblet', 'microfold')
cellchat.18moD7 <- setIdent(cellchat.18moD7, ident.use = "cell_type_subset") # set "cell_type" as default cell identity


#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...
object.list <- list(mo02d0 = cellchat.02moD0, mo18d0 = cellchat.18moD0,
                    mo02d1 = cellchat.02moD1, mo18d1 = cellchat.18moD1,
                    mo02d4 = cellchat.02moD4, mo18d4 = cellchat.18moD4,
                    mo02d7 = cellchat.02moD7, mo18d7 = cellchat.18moD7)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


user_defined_palette =  c('#F6222E', '#3283FE', '#16FF32',  '#BDCDFF', '#AA0DFE','#1CFFCE', '#d62728', '#19c9b3','#FFA5D2',   'grey', '#2ED9FF', '#c1c119', '#8b0000', '#3B00FB', '#FE00FA', "#F8A19F", '#1CBE4F','#B5EFB5', '#BEFFF7','#FEAF16', '#325A9B', '#C075A6')

# Zoom in pathways of interest

# FGF
par(mfrow=c(1,1))
p0 <- netVisual_chord_cell(object.list[[1]],  slot.name = "netP", signaling='FGF', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "FGF (02mo D0)")
p1 <- netVisual_chord_cell(object.list[[7]],  slot.name = "netP", signaling='FGF', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "FGF (02mo D7)")
p2 <- netVisual_chord_cell(object.list[[2]],  slot.name = "netP", signaling='FGF', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "FGF (18mo D0)")
p3 <- netVisual_chord_cell(object.list[[8]],  slot.name = "netP", signaling='FGF', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "FGF (18mo D7)")
(ggdraw(p0,xlim = c(0,1), ylim = c(0,1), clip = "on") | ggdraw(p2,xlim = c(0,1), ylim = c(0,1), clip = "on")) /
(ggdraw(p1, xlim = c(0,1), ylim = c(0,1), clip = "on") | ggdraw(p3, xlim = c(0,1), ylim = c(0,1), clip = "on"))

par(mfrow=c(1,1))
p0 <- netVisual_chord_cell(object.list[[3]],  slot.name = "netP", signaling='FGF', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "FGF (02mo D1)")
p1 <- netVisual_chord_cell(object.list[[5]],  slot.name = "netP", signaling='FGF', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "FGF (02mo D4)")
p2 <- netVisual_chord_cell(object.list[[4]],  slot.name = "netP", signaling='FGF', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "FGF (18mo D1)")
p3 <- netVisual_chord_cell(object.list[[6]],  slot.name = "netP", signaling='FGF', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "FGF (18mo D4)")
(ggdraw(p0,xlim = c(0,1), ylim = c(0,1), clip = "on") | ggdraw(p2,xlim = c(0,1), ylim = c(0,1), clip = "on")) /
  (ggdraw(p1, xlim = c(0,1), ylim = c(0,1), clip = "on") | ggdraw(p3, xlim = c(0,1), ylim = c(0,1), clip = "on"))

# BMP
par(mfrow=c(1,1))
p0 <- netVisual_chord_cell(object.list[[1]],  slot.name = "netP", signaling='BMP', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "BMP (02mo D0)")
p1 <- netVisual_chord_cell(object.list[[7]],  slot.name = "netP", signaling='BMP', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "BMP (02mo D7)")
p2 <- netVisual_chord_cell(object.list[[2]],  slot.name = "netP", signaling='BMP', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "BMP (18mo D0)")
p3 <- netVisual_chord_cell(object.list[[8]],  slot.name = "netP", signaling='BMP', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "BMP (18mo D7)")
(ggdraw(p0,xlim = c(0,1), ylim = c(0,1), clip = "on") | ggdraw(p2,xlim = c(0,1), ylim = c(0,1), clip = "on")) /
(ggdraw(p1, xlim = c(0,1), ylim = c(0,1), clip = "on") | ggdraw(p3, xlim = c(0,1), ylim = c(0,1), clip = "on"))

par(mfrow=c(1,1))
p0 <- netVisual_chord_cell(object.list[[3]],  slot.name = "netP", signaling='BMP', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "BMP (02mo D1)")
p1 <- netVisual_chord_cell(object.list[[5]],  slot.name = "netP", signaling='BMP', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "BMP (02mo D4)")
p2 <- netVisual_chord_cell(object.list[[4]],  slot.name = "netP", signaling='BMP', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "BMP (18mo D1)")
p3 <- netVisual_chord_cell(object.list[[6]],  slot.name = "netP", signaling='BMP', annotationTrackHeight = c(0.1),  lab.cex = 0.5, small.gap = 4, color.use = user_defined_palette, title.name = "BMP (18mo D4)")
(ggdraw(p0,xlim = c(0,1), ylim = c(0,1), clip = "on") | ggdraw(p2,xlim = c(0,1), ylim = c(0,1), clip = "on")) /
  (ggdraw(p1, xlim = c(0,1), ylim = c(0,1), clip = "on") | ggdraw(p3, xlim = c(0,1), ylim = c(0,1), clip = "on"))