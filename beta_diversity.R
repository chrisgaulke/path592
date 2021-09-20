#Beta Diversity


# "The vegan package provides tools for descriptive community ecology.
# It has most basic functions of diversity analysis, community ordination
# and dissimilarity analysis. Most of its multivariate tools can be
# used for other data types as well." Taken from the vegan help page

library(vegan)
library(ggplot2)
library(ape)


#use random distribution generators to make a neg binomial distributed matrix
#this will be a reasonable-ish approximation of microbiome data

set.seed(731)
big_df1 <- matrix(rnbinom(10000, size = .05, prob = .001),
                 nrow = 100,
                 ncol = 100)

big_df2 <- matrix(rnbinom(10000, size = .05, prob = .001),
                  nrow = 100,
                  ncol = 100)

colnames(big_df1) <- paste("patient", 1:100,sep = "_" )
rownames(big_df1) <- paste("species", 1:100,sep = "_" )

colnames(big_df2) <- paste("patient", 101:200,sep = "_" )
#rownames(big_df2) <- paste("species", 1:100,sep = "_" )
rownames(big_df2) <- paste("species", 51:150,sep = "_" )


big_df <- merge(big_df1,big_df2,by  = 0, all =T )
big_df[is.na(big_df)] <- 0
rownames(big_df) <- big_df$Row.names
big_df$Row.names <- NULL
dim(big_df)

# Rarefy ------------------------------------------------------------------

# Because there were not signficant differences between relative abundance and
# rarefied data abundances we will use rarefied data with the same thresholds
# used in the previous

#remove all samples < 1100

big_df.rare <- big_df[,which(colSums(big_df) >= min(colSums(big_df)))]

#remove all species == 0
big_df.rare <- big_df.rare[which(rowSums(big_df) > 0),]

# rrarefy, note the second "r", requires sampling the rows not colums
# the t() function transposes data (i.e., rows become columns)

big_df.rare <- rrarefy(t(big_df.rare), sample = min(colSums(big_df)))

#check to make sure all
all(colSums(big_df.rare) > 0)


# Beta diversity ----------------------------------------------------------

#You should read the vegdist help file

#Two popular diversity metrics

big_df_rare.bray <- vegdist(x = big_df.rare,
                            method = "bray",
                            binary = FALSE,
                            upper = FALSE,)

big_df_rare.jacc <- vegdist(x = big_df.rare,
                            method = "jaccard",
                            binary = FALSE,
                            upper = FALSE,)

# for these data the distances are similar, but for other data it can make a
# big difference
cor.test(big_df_rare.bray, big_df_rare.jacc)
plot(big_df_rare.bray, big_df_rare.jacc,xlim = c(0,1), ylim=c(0,1))

#View(as.data.frame(as.matrix(big_df_rare.bray)))



# Principle Components Analysis -------------------------------------------


big_df.prcomp <- prcomp(big_df.rare)

#to determine how much variance is explained lets look at the summary
summary(big_df.prcomp)

big_df_prcomp.df <- as.data.frame(big_df.prcomp$x[,1:5])
big_df_prcomp.df$group <- rep(c(1,2), times = c(100,100))


big_df_prcomp.plot <- ggplot(big_df_prcomp.df,
                             aes(x = PC1,
                                 y = PC2,
                                 color = factor(group))
                             )

big_df_prcomp.plot+ geom_point()

# Principle Coordinates Analysis -------------------------------------------

# PCoA wants a matrix of dist object. We will provide both distance objects
# separately, starting with bray

big_df_bray.pcoa <- pcoa(big_df_rare.bray)

#With PCoA the we don't get percent variation but an analogous value of
# brokenstick can be accessed
View(big_df_bray.pcoa$values)

big_df_bray_pcoa.df <- as.data.frame(big_df_bray.pcoa$vectors[,1:5])
big_df_bray_pcoa.df$group <- rep(c(1,2), times = c(100,100))

big_df_bray_pcoa.plot <- ggplot(big_df_bray_pcoa.df,
                             aes(x = Axis.1,
                                 y = Axis.2,
                                 color = factor(group))
)

big_df_bray_pcoa.plot+ geom_point()

#now jaccard

big_df_jacc.pcoa <- pcoa(big_df_rare.jacc)

big_df_jacc_pcoa.df <- as.data.frame(big_df_jacc.pcoa$vectors[,1:5])
big_df_jacc_pcoa.df$group <- rep(c(1,2), times = c(100,100))

big_df_jacc_pcoa.plot <- ggplot(big_df_jacc_pcoa.df,
                                aes(x = Axis.1,
                                    y = Axis.2,
                                    color = factor(group))
)

big_df_jacc_pcoa.plot+ geom_point()



# Non-metric Multidimensional Scaling ------------------------------------

#Now we will try the same approach but with NMDS. For NMDS we have to decide on
# a number of dimensions. The default is 2, but most of the time that will
# yield a very high stress. To attempt to minimize stress we will try a few
# values of k (dimensions) and find the one the results in lowest stress
# our target is around or below .1

(metaMDS(big_df.rare,distance = "bray", k = 2))$stress # too high to be useful
(metaMDS(big_df.rare,distance = "bray", k = 5))$stress # better, but still high
(metaMDS(big_df.rare,distance = "bray", k = 10))$stress # better, but still high
(metaMDS(big_df.rare,distance = "bray", k = 12))$stress # good


big_df_rare_bray.NMDS <- metaMDS(big_df.rare,distance = "bray", k = 12)

# Unlike PCA and PCoA, no % of variance associated with each axis in NMDS

big_df_rare_bray_nmds.df <- as.data.frame(big_df_rare_bray.NMDS$points)
big_df_rare_bray_nmds.df$group <- rep(c(1,2), times = c(100,100))

big_df_rare_bray_nmds.plot <- ggplot(big_df_rare_bray_nmds.df,
                                aes(x = MDS1,
                                    y = MDS2,
                                    color = factor(group))
)

big_df_rare_bray_nmds.plot + geom_point()

#now jaccard
big_df_rare_jacc.NMDS <- metaMDS(big_df.rare,distance = "jaccard", k = 12)
big_df_rare_jacc_nmds.df <- as.data.frame(big_df_rare_jacc.NMDS$points)

big_df_rare_jacc_nmds.df$group <- rep(c(1,2), times = c(100,100))

big_df_rare_jacc_nmds.plot <- ggplot(big_df_rare_jacc_nmds.df,
                                     aes(x = MDS1,
                                         y = MDS2,
                                         color = factor(group))
)

big_df_rare_jacc_nmds.plot + geom_point()





