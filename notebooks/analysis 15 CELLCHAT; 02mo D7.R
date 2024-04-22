# 10 CELLCHAT; 02mo D7

renv::status()
renv::snapshot()
renv::restore()

# Biobase, BiocGenerics, ComplexHeatmap
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
#install.packages("patchwork", force=TRUE)
install.packages("reticulate", force=TRUE)

library(CellChat)
library(patchwork)
library(reticulate)
library(anndata)
library(anndata)
memory.limit(9999999999)

options(stringsAsFactors = FALSE)
adata <- import("anndata", convert = FALSE)
adata_object <- adata$read_h5ad("../output/cellchat_io/cd45neg_0147_SLTBI_vdb_02mo_d7_NEW.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(adata_object$X))
rownames(data.input) <- rownames(py_to_r(adata_object$var))
colnames(data.input) <- rownames(py_to_r(adata_object$obs))

# access meta data
meta.data <- py_to_r(adata_object$obs)
meta <- meta.data

cellchat <- createCellChat(object = data.input)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type_subset")


CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
#dplyr::glimpse(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 20)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

options(repr.plot.width=14, repr.plot.height=7)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netAnalysis_contribution(cellchat, signaling = pathways.show)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

saveRDS(cellchat, file = "../output/cellchat_io/cd45neg_0147_SLTBI_vdb_02mo_D7_NEW_cellchatDB.rds")

sessionInfo()
