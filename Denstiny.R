###12/12/2019
###Johannes Schroth

#Diffusion Pseutotime using the Destiny R Package

install.packages('ggbeeswarm')
library(ggbeeswarm)

packages_1 <- c('cytofkit', 'flowViz', 'ggcyto', 'gridExtra','grid', 'ggplot2', 'flowCore',
                'umap', 'dplyr', 'reshape', 'RColorBrewer', 'FlowSOM', 'Rtsne',
                'ConsensusClusterPlus', 'ggdendro', 'plotly', 'shiny', 'plotly', 'destiny', 'Biobase', 'slingshot',
                'ggbeeswarm')
lapply(packages_1, require, character.only=TRUE)

all_fcs <- read.FCS(filename = '/Users/johannesschroth/Desktop/iEMRA clustering stuff/CD8_1.fcs')
colnames(all_fcs)
colnames(all_fcs) <- c(colnames(all_fcs)[1:6], 'CD8', 'ADAM28', 'ZOMBIE', 'CD27', 'CD45RA', 'CD28', 'KLRG1', 'CD4', 'CCR7', 'Time', 'Sample')

vars <- colnames(all_fcs)[-grep('FSC|SSC|Time|Sample', colnames(all_fcs))]
all_fcs_transformed <- transform(all_fcs, estimateLogicle(all_fcs, vars))

df_all <- as.matrix(flowCore::exprs(all_fcs_transformed))[,c(10:13,15)]
sce_all <- SingleCellExperiment(as.matrix(t(df_all)))
data_es <- as.ExpressionSet(as.data.frame(df_all))

dm <- DiffusionMap(data_es, sigma = optimal_sigma(find_sigmas(data_es)))
dpt <- DPT(dm)


Clusters <- read.csv('~/Desktop/iEMRA clustering stuff/norm_iEMRAidea.csv')

pca_all <- as.matrix(prcomp(df_all)$x[,1:2])

colData(sce_all)$Clusters <- as.numeric(Clusters$Clusters)
reducedDim(sce_all) <- as.matrix(Clusters[,8:9])

sce_all <- slingshot(sce_all, clusterLabels = 'Clusters', reducedDim = reducedDim(sce_all))
metadata(sce_all)

colData(sce_all)
colData(sce_all)[,2]
rank(colData(sce_all)[,2:5])

colData(sce_all)

dat <- as.data.frame(cbind(Clusters,colData(sce_all)))

d <- reshape::melt(dat, c('slingPseudotime_1','slingPseudotime_2'), c('ADAM28',  'CD27', 'CD45RA', 'CD28', 'KLRG1', 'CCR7'))
d <- na.omit(d)
d <- d[d$slingPseudotime_1 > 0,]
d <- d[d$slingPseudotime_1 < max(d$slingPseudotime_1),]

head(cbind(dpt$DPT1,dpt$DPT2))

head(d)

dat['dpt_pseudotime'] <- rank(dpt$dpt)
dat['branch'] <- dpt$Branch

length(dpt$dpt)

length(rank(na.omit(slingPseudotime(sce_all))))

dpt$Branch
eigenvalues(dm)

ggplot(dat, aes(dat$dpt_pseudotime, as.factor(Clusters), colour = as.factor(Clusters))) +
    geom_quasirandom(groupOnX = FALSE)

ggplot(dat, aes(UMAP1, UMAP2)) +
  geom_point(alpha = 0.1) +
  geom_smooth(data = dat[dat$branch == 1, ], aes(colour = 'red')) +
  geom_smooth(data = dat[dat$branch == 2, ], aes(colour = 'red')) +
  geom_smooth(data = dat[dat$branch == 3, ], aes(colour = 'red'))

names(dpt$branch)


head(dat)


x <- dat$UMAP1[-na]
y <- dat$UMAP2[-na]
dat$UMAP1[which(is.na(dat[,14]),TRUE),]

na <- which(is.na(dat[,14]),TRUE)

dat$UMAP1[-na]

(dat$slingPseudotime_1)


dat$UMAP1[na.omit(dat[,14])]

ggplot(dm, aes(DC1, DC2, colour = as.factor(Clusters))) +
  geom_density(alpha = 0.4) +
  geom_point(UMAP1 ~aes(colour = as.factor(Clusters)))


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



