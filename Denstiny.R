###12/12/2019
###Johannes Schroth

#Diffusion Pseutotime using the Destiny R Package

packages_1 <- c('cytofkit', 'flowViz', 'ggcyto', 'gridExtra','grid', 'ggplot2', 'flowCore',
                'umap', 'dplyr', 'reshape', 'RColorBrewer', 'FlowSOM', 'Rtsne',
                'ConsensusClusterPlus', 'ggdendro', 'plotly', 'shiny', 'plotly', 'destiny', 'Biobase')
lapply(packages_1, require, character.only=TRUE)

all_fcs <- read.FCS(filename = '/Users/johannesschroth/Desktop/iEMRA clustering stuff/CD8_1.fcs')
colnames(all_fcs)
colnames(all_fcs) <- c(colnames(all_fcs)[1:6], 'CD8', 'ADAM28', 'ZOMBIE', 'CD27', 'CD45RA', 'CD28', 'KLRG1', 'CD4', 'CCR7', 'Time', 'Sample')

vars <- colnames(all_fcs)[-grep('FSC|SSC|Time|Sample', colnames(all_fcs))]
all_fcs_transformed <- transform(all_fcs, estimateLogicle(all_fcs, vars))
all_fcs_transformed

df_all <- as.matrix(flowCore::exprs(all_fcs_transformed))[,c(10:13,15)]

data_es <- as.ExpressionSet(as.data.frame(df_all))

sigmas <- find_sigmas(data_es)

dm <- DiffusionMap(data_es, sigma = optimal_sigma(sigmas))
dpt <- DPT(dm)

reducedDim(data_es, type = )

Clusters <- read.csv('~/Desktop/iEMRA clustering stuff/norm_iEMRAidea.csv')[,12]


ggplot(dm, aes(DC1, DC2, colour = as.factor(Clusters))) +
  geom_density2d(alpha = 0.4) +
  geom_point(aes(colour = as.factor(Clusters)))


plot

ggplot(dm, aes(DC1,DC2,  colour = Clusters)) +
  geom_point() +
  theme_classic() +
  scale_color_gradient2(low="#4575B4", mid= 'grey90' , high='#D73027', midpoint=mean(all_data$Pseudotime), 
                        limits=c(min(dm$CD45RA),max(dm$CD45RA)),
                        guide = guide_colourbar(barwidth = 0.5, barheight = 10))
  


plot(dpt, root = 1, paths_to = c(2,3), col_by = 'branch')

df <- as.data.frame(cbind(DC1 = dm$DC1, DC2 = dm$DC2, DC3 = dm$DC3, Clusters))
df
plot_ly(x = df$DC1, y = df$DC2, z = df$DC3, 
        color = df$Clusters, type = 'scatter3d', mode = 'markers', text = df$Clusters)



