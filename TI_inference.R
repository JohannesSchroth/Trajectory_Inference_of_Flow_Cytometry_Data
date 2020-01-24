
library(dyno)
library(tidyverse)


all_fcs <- read.FCS(filename = '/Users/johannesschroth/Desktop/iEMRA clustering stuff/CD8_1.fcs')
colnames(all_fcs)
colnames(all_fcs) <- c(colnames(data)[1:6], 'CD8', 'ADAM28', 'ZOMBIE', 'CD27', 'CD45RA', 'CD28', 'KLRG1', 'CD4', 'CCR7', 'Time', 'Sample')

vars <- colnames(all_fcs)[-grep('FSC|SSC|Time|Sample', colnames(all_fcs))]
all_fcs_transformed <- transform(all_fcs, estimateLogicle(all_fcs, vars))
all_fcs_transformed

df_all <- as.matrix(flowCore::exprs(all_fcs_transformed))[,c(10:13,15)]
rownames(df_all) <- seq(1,length(df_all[,1]),1)

dataset <- wrap_expression(expression = df_all, counts = df_all)
umap <- as.matrix(read.csv('~/Desktop/iEMRA clustering stuff/norm_iEMRAidea.csv')[,8:9])
rownames(umap) <- rownames(df_all)

dataset <- add_dimred(dataset, umap)
guidelines_shiny(dataset)

dataset$expression

dataset <- add_prior_information(dataset, start_id = dataset$cell_ids[5:50])

model <- infer_trajectory(dataset, ti_paga(), verbose = T)
model$pseudotime

install.packages('paga')


plot_dimred(model, color_cells = pseudotime)

