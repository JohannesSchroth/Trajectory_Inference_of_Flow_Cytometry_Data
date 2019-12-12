###12/12/2019
###Johannes Schroth

#Diffusion Pseutotime using the Destiny R Package


packages_1 <- c('cytofkit', 'flowViz', 'ggcyto', 'gridExtra','grid', 'ggplot2', 'flowCore',
                'umap', 'dplyr', 'reshape', 'RColorBrewer', 'FlowSOM', 'Rtsne',
                'ConsensusClusterPlus', 'ggdendro', 'plotly', 'shiny', 'plotly', 'monocle3')
lapply(packages_1, library, character.only=TRUE)

all_fcs <- read.FCS(filename = '/Users/johannesschroth/Desktop/iEMRA clustering stuff/CD4_1.fcs')
colnames(all_fcs)
colnames(all_fcs) <- c(colnames(data)[1:6], 'CD8', 'ADAM28', 'ZOMBIE', 'CD27', 'CD45RA', 'CD28', 'KLRG1', 'CD4', 'CCR7', 'Time', 'Sample')

vars <- colnames(all_fcs)[-grep('FSC|SSC|Time|Sample', colnames(all_fcs))]
all_fcs_transformed <- transform(all_fcs, estimateLogicle(all_fcs, vars))
all_fcs_transformed

df_all <- as.matrix(flowCore::exprs(all_fcs_transformed))[,c(10:13,15)]