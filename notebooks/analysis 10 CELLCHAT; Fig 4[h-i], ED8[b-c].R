# 10 CELLCHAT MAIN

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

p0 <- compareInteractions(cellchat, show.legend = F, color.use = rep(c('#76D6FF', '#FF8072'),4), group = c(1:8), digits = 1, size.text = 10) 
p1 <- compareInteractions(cellchat, show.legend = F, color.use = rep(c('#76D6FF', '#FF8072'),4), group = c(1:8), digits = 1, measure='weight', size.text = 10)
p0|p1

# Circos plots with grouped data
user_defined_palette_groups <- c('#3283FE',  '#aa40fc', '#d62728', '#19c9b3', "#2ED9FF", '#c1c119', "#8b0000",  '#3B00FB',  "#FE00FA", '#1CBE4F', '#ff7f0e')
group.cellType <- c(rep("EC",3), rep("FB",3),   "MEC", "vSMCPC", 'FB', "FB", "aaTEC1", "aaTEC2",  'cTEC',  'earlyPr', 'mTEC1', rep("mTECdiff", 2),  rep("mimetic", 5))
group.cellType <- factor(group.cellType, levels = c("EC", "FB", "MEC", "vSMCPC",  "aaTEC1", "aaTEC2", "cTEC", 'earlyPr', "mTEC1", "mTECdiff", "mimetic"))

user_defined_palette_groups <- c('#ff7f0e', '#BDCDFF', '#aa40fc',  '#1CFFCE', '#d62728', '#19c9b3', "#2ED9FF", '#c1c119', "#8b0000",  '#3B00FB',  "#FE00FA", '#1CBE4F', 'black')
group.cellType <- c(rep("EC",3), 'capsFB', 'intFB', 'medFB',  "MEC", "vSMCPC", 'medFB', "intFB", "aaTEC1", "aaTEC2",  'cTEC',  'earlyPr', 'mTEC1', rep("mTECdiff", 2),  rep("mimetic", 5))
group.cellType <- factor(group.cellType, levels = c("EC", 'capsFB', 'intFB', 'medFB',  "MEC", "vSMCPC",  "aaTEC1", "aaTEC2", "cTEC", 'earlyPr', "mTEC1", "mTECdiff", "mimetic"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged") )

par(mfrow = c(1,1), xpd=TRUE)

p0 <- netVisual_circle(object.list[[1]]@net$count.merged,   weight.scale = T, label.edge = T, edge.weight.max =400,  color.use = user_defined_palette_groups, title.name = paste0("02-mo D0\n"))
p1 <- netVisual_circle(object.list[[2]]@net$count.merged,   weight.scale = T, label.edge = T, edge.weight.max =400,  color.use = user_defined_palette_groups, title.name = paste0("18-mo D0\n"))


user_defined_palette =  c('#F6222E', '#3283FE', '#16FF32',  '#BDCDFF', '#AA0DFE','#1CFFCE', '#d62728', '#19c9b3','#FFA5D2',   'grey', '#2ED9FF', '#c1c119', '#8b0000', '#3B00FB', '#FE00FA', "#F8A19F", '#1CBE4F','#B5EFB5', '#BEFFF7','#FEAF16', '#325A9B', '#C075A6')

par(mfrow=c(1,2), xpd=TRUE)
weight.max <- getMaxWeight(object.list, attribute=c("idents", "count"))
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale =T, label.edge = F, color.use = user_defined_palette,
                   edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("# of interactions - ", names(object.list)[i]))
}

user_defined_palette =  c('#F6222E', '#3283FE', '#16FF32',  '#BDCDFF', '#AA0DFE','#1CFFCE', '#d62728', '#19c9b3','#FFA5D2',   'grey', '#2ED9FF', '#c1c119', '#8b0000', '#3B00FB', '#FE00FA', "#F8A19F", '#1CBE4F','#B5EFB5', '#BEFFF7','#FEAF16', '#325A9B', '#C075A6')

group.new = levels(cellchat.18moD7@meta$cell_type_subset)
cellchat <- liftCellChat(cellchat, group.new) 
p0 <- netVisual_heatmap(cellchat,  comparison = c(1,2), color.heatmap = c("#0000FF", "#00FF00"), color.use = user_defined_palette)

secreted <- c('ANGPTL', 'BMP', 'TWEAK', 'EGF', 'FGF', 'MK', 'PTN', 'IGF', 'VISFATIN', 'NRG', 'SPP1', 'LIFR', 'TGFb', 'IL4', 'GAS', 'CCL', 'CD70', 'MIF', 'GALECTIN', 'PROS')
cell_cell <- c('APP', 'JAM', 'NOTCH', 'EPHA', 'EPHB', 'MPZ', 'CADM', 'CDH', 'DESMOSOME', 'OCLN', 'CD80', 'NRXN', 'PTPRM', 'NCAM')
ecm_R <- c('COLLAGEN', 'FN1', 'LAMININ', 'HSPG', 'TENASCIN', 'VTN')

netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = secreted, sources.use=c(1:8,11:25), targets.use = c("aaTEC1"), annotationTrackHeight = c(0.055), scale = TRUE, lab.cex = 0.0001, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = secreted, sources.use=c(1:8,11:25), targets.use = c("aaTEC2"), annotationTrackHeight = c(0.045), scale = TRUE, lab.cex = 0.0001, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = secreted, sources.use=c(1:8,11:25), targets.use = c("mTEC1"), annotationTrackHeight = c(0.09), scale = TRUE, lab.cex = 0.0001, small.gap = 0.1,  color.use = user_defined_palette, show.legend = F)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = secreted, sources.use=c(1:8,11:25), targets.use = c("earlyPr"), annotationTrackHeight = c(0.07), scale = TRUE, lab.cex = 0.0001, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)

netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = cell_cell, sources.use=c(1:8,11:25),  targets.use = c("aaTEC1"), annotationTrackHeight = c(0.1), scale = TRUE, lab.cex = 0.4, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = cell_cell, sources.use=c(1:8,11:25),  targets.use = c("aaTEC2"), annotationTrackHeight = c(0.1), scale = TRUE, lab.cex = 0.4, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = cell_cell, sources.use=c(1:8,11:25),  targets.use = c("mTEC1"), annotationTrackHeight = c(0.1), scale = TRUE, lab.cex = 0.4, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = cell_cell, sources.use=c(1:8,11:25),  targets.use = c("earlyPr"), annotationTrackHeight = c(0.1), scale = TRUE, lab.cex = 0.4, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)

netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = ecm_R, sources.use=c(1:8,11:25),  targets.use = c("aaTEC1"), annotationTrackHeight = c(0.1), scale = TRUE, lab.cex = 0.4, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = ecm_R, sources.use=c(1:8,11:25),  targets.use = c("aaTEC2"), annotationTrackHeight = c(0.1), scale = TRUE, lab.cex = 0.4, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)
netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = ecm_R, sources.use=c(1:8,11:25),  targets.use = c("mTEC1"), annotationTrackHeight = c(0.1), scale = TRUE, lab.cex = 0.4, small.gap = 0.1, color.use = user_defined_palette, show.legend = F) netVisual_chord_gene(object.list[[2]], slot.name = 'netP', signaling = ecm_R, sources.use=c(1:8,11:25),  targets.use = c("earlyPr"), annotationTrackHeight = c(0.1), scale = TRUE, lab.cex = 0.4, small.gap = 0.1, color.use = user_defined_palette, show.legend = F)
