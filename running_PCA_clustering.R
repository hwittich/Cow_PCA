# loading all necessary packages
library(ade4)
library(ape)
library(data.table)
library(factoextra)
library(ggpubr)
library(ggtree)
library(hierfstat)
library(tidyverse)

# reading genotype data
file <- fread('~/stat488/cattle_genotypes_mod.txt')

# removes first column (population IDs)
file_mod <- file[,-1] %>% as_tibble()

# There are some non-numeric values hidden in the data frame,
# so this should get rid of them.
file_mod <- apply(file_mod, 2, as.numeric) %>% as.data.frame()

# replace NAs by the average genotype because the PCA function used by the authors
# doesn't handle NAs
file_mod <- replace_na(file_mod, as.list(colMeans(file_mod,na.rm=T)))

# removing columns with variance == 0 (that is, columns in which the genotype is the same for everyone)
file_mod[is.nan(file_mod)] <- 0
file_mod <- file_mod[,which(apply(file_mod, 2, var)!=0)]

# runs PCA (using the function from the ade4 package, which the authors used)
res.pca <- dudi.pca(file_mod, # input data frame
                	scannf = FALSE, # do not output scree plot
                	center = TRUE, # center it (they don't mention in the methods, but why not?)
                	scale = FALSE, # genotypes are already in the same scale, so no need to scale it in here
                	nf = 10) # number of components kept in the results

# scree plot
fviz_eig(res.pca)
ggsave('~/stat488/scree_plot.png',width=7,height=5)

# extract different PC columns and convert to data frame
pcaData <- as.data.frame(res.pca$li[, c(1,2)])
pcaData1.3 <- as.data.frame(res.pca$li[, c(1,3)])
pcaData2.3 <- as.data.frame(res.pca$li[, c(2,3)])

# merge with population IDs
pcaData <- cbind(pcaData, file$V1)
colnames(pcaData) <- c('PC1', 'PC2', 'Population')
pcaData1.3 <- cbind(pcaData1.3, file$V1)
colnames(pcaData1.3) <- c('PC1', 'PC3', 'Population')
pcaData2.3 <- cbind(pcaData2.3, file$V1)
colnames(pcaData2.3) <- c('PC2', 'PC3', 'Population')

# plot it
a <- ggplot(pcaData, aes(PC1, PC2, color = Population)) + geom_point()
b <- ggplot(pcaData1.3, aes(PC1, PC3, color = Population)) + geom_point()
c <- ggplot(pcaData2.3, aes(PC2, PC3, color = Population)) + geom_point()
ggarrange(b, c, a, common.legend = TRUE, legend='right')
ggsave('~/stat488/PC_plots.png',width=12,height=10)

##Clustering
#Append population IDs onto file_mod
file_mod.pop <- cbind(file$V1, as.data.frame(file_mod))
colnames(file_mod.pop)[1] <- 'V1'

#Calculate pairwise allele sharing distance for all populations
# Use dist.gene() on SNP-level data to produce distance matrix
# From the ape documentation:
# This function is meant to be very general and accepts different kinds of data
# (alleles, haplotypes, SNP, DNA sequences, ...). The rows of the data matrix represent
# the individuals, and the columns the loci.
ASD_matrix <- dist.gene(x = file_mod.pop,
                    	method = 'pairwise', # distance (d) between two individuals is the number of loci for which they differ
                    	pairwise.deletion = TRUE, # Delete columns with missing values when calculating pairwise distances #shouldn't be a problem but why not
                    	variance = FALSE) # Don't return variances of distances

#Use nj() on distance matrix to build tree
# From the documentation:
# This function performs the neighbor-joining tree estimation of Saitou and Nei (1987).
tree <- nj(ASD_matrix) # One argument: distance matrix

# updating tree leaf-nodes labels to be population IDs
tree_tbl <- as_tibble(tree)
tree_tbl$label[1:nrow(file_mod.pop)] <- file_mod.pop$V1
tree_t <- as.phylo(tree_tbl)

write.tree(tree_t, '~/stat488/STAT488_cattle_ALL_analysis_phylo_tree.txt')

ggtree(tree_t, layout='circular', aes(color=label)) + geom_tiplab(size=5, aes(angle=angle)) +
  labs(color='Population')
ggsave('~/stat488/all_populations_ASD_tree.png',width=15,height=15)

#Calculate pairwise allele sharing distance, filtered down to ten populations
#Should replicate figure 3 in paper
pops <- c('ABO','BIS','CHE','CHF','GUE','HOL','MON','OUL','TAR','TID')
file_filtered <- file_mod.pop %>% filter(V1 %in% pops)
#unique(file_filtered$V1)

ASD_matrix <- dist.gene(x = file_filtered,
                    	method = 'pairwise', # distance (d) between two individuals is the number of loci for which they differ
                    	pairwise.deletion = TRUE, # Delete columns with missing values when calculating pairwise distances #shouldn't be a problem but why not
                    	variance = FALSE) # Don't return variances of distances

#Use nj() on distance matrix to build tree
# From the documentation:
# This function performs the neighbor-joining tree estimation of Saitou and Nei (1987).
tree <- nj(ASD_matrix) # One argument: distance matrix

# updating tree leaf-nodes labels to be population IDs
tree_tbl <- as_tibble(tree)
tree_tbl$label[1:nrow(file_filtered)] <- file_filtered$V1
tree_t <- as.phylo(tree_tbl)

write.tree(tree_t, '~/stat488/STAT488_cattle_filtered_populations_analysis_phylo_tree.txt')

ggtree(tree_t, layout='circular', aes(color=label)) + geom_tiplab(size=5, aes(angle=angle)) +
  labs(color='Population')
ggsave('~/stat488/filtered_populations_ASD_tree.png',width=15,height=15)

#Calculate pairwise Fst between populations
#Use genet.dist() from hierfstat package to produce population pairwise Fst distance matrix
# Fst can be calulcated with the following equation:
# ((avg pairwise ASD between populations)-(avg pairwise ASD within populations)) / (avg pairwise ASD between populations)
Fst_matrix <- genet.dist(dat = file_filtered, #first column is population, rest of columns are genotypes in dosage format
                     	diploid = TRUE, #cattle are diploid organisms
                     	method = 'Fst') #Calculate Fst distance metric

#Use nj() on distance matrix to build tree
# From the documentation:
# This function performs the neighbor-joining tree estimation of Saitou and Nei (1987).
tree <- nj(Fst_matrix) # One argument: distance matrix
write.tree(tree, '~/stat488/STAT488_cattle_filtered_populations_fst_analysis_phylo_tree.txt')

ggtree(tree, layout='circular', aes(color=label)) + geom_tiplab(size=5, aes(angle=angle)) +
  labs(color='Population')
ggsave('~/stat488/filtered_populations_fst_tree.png',width=4,height=4)
