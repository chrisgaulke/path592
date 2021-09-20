#Alpha Diversity

# "The vegan package provides tools for descriptive community ecology.
# It has most basic functions of diversity analysis, community ordination
# and dissimilarity analysis. Most of its multivariate tools can be
# used for other data types as well." Taken from the vegan help page

library(vegan)

#Make some faux data

df <- data.frame(patient_1 = c(5,5,5,5,5,5),
                 patient_2 = c(0,0,0,0,0,1000),
                 patient_3 = c(0,50,5,100,10,1),
                 patient_4 = c(0,0,5,5,5,0),
                 patient_5 = c(100,100,100,100,100,100))

rownames(df) <- c("species_a",
                  "species_b",
                  "species_c",
                  "species_d",
                  "species_e",
                  "species_f")


#use random distribution generators to make a neg binomial distributed matrix
#this will be a reasonable-ish approximation of microbiome data

set.seed(731)
big_df <- matrix(rnbinom(10000, size = .2, prob = .01),
                 nrow = 100,
                 ncol = 100)

colnames(big_df) <- paste("patient", 1:100,sep = "_" )
rownames(big_df) <- paste("species", 1:100,sep = "_" )


# Rarefy ------------------------------------------------------------------

#Rarefaction is necessary to conduct alpha diversity analysis

rarecurve(x = df,
          step = 10,
          )

#set sample size
sample = min(colSums(big_df)) /10 # you want enough sample points to make smooth

rarecurve(x = big_df,
          step = floor(sample),
          label = FALSE

)

abline(v = 100)
abline(v = 500)
abline(v = 1100) # looks like a few sample

#remove all samples < 1100
big_df.rare <- big_df[,which(colSums(big_df) > 1099)]

#remove all species == 0
big_df.rare <- big_df.rare[which(rowSums(big_df) > 0),]

# rrarefy, note the second "r", requires sampling the rows not columns
# the t() function transposes data (i.e., rows become columns)

big_df.rare <- rrarefy(t(big_df.rare), sample = 1100)

#check to make sure all
all(colSums(big_df.rare) > 0)


big_df.rel <- big_df[,which(colSums(big_df) > 1099)]

big_df.rel <- sweep(x = big_df.rel,
                    MARGIN = 2,
                    STATS = colSums(big_df.rel),
                    FUN = '/'
                    )
big_df.rel <- t(big_df.rel)

#How different are relative abundance and rarefied data tables?

cor.test(matrix(big_df.rare),matrix(big_df.rel))
cor.test(vegdist(big_df.rare),vegdist(big_df.rel))


#df will not be rarefied

# Richness ----------------------------------------------------------------

# The function in R::vegan to determine richness is specnumber

specnumber(x = df,
          MARGIN = 2 #by columns
            )

#unrarefied
big_df.richness <- specnumber(x = big_df,
                    MARGIN = 2 #by columns
                    )

#rarefied
big_df_rare.richness <- specnumber(x = big_df.rare,
                    MARGIN = 1 #by columns
                    )

hist(big_df.richness)
hist(big_df_rare.richness) #rarefied data frame


# Shannon -----------------------------------------------------------------

diversity(df,
          index = "shannon",
          MARGIN = 2
          )

big_df.shannon <- diversity(big_df,
          index = "shannon",
          MARGIN = 2
)

big_df_rare.shannon <- diversity(big_df.rare,
                            index = "shannon",
                            MARGIN = 1
)


hist(big_df.shannon)
hist(big_df_rare.shannon)



# Shannon -----------------------------------------------------------------

diversity(df,
          index = "simpson",
          MARGIN = 2
)

big_df.simpson <- diversity(big_df,
                            index = "simpson",
                            MARGIN = 2
)

big_df_rare.simpson <- diversity(big_df.rare,
                                 index = "simpson",
                                 MARGIN = 1
)


hist(big_df.simpson)
hist(big_df_rare.simpson)

# Richness Correlations ---------------------------------------------------

#Linear relationships cannot be assumed so we will use nonparametic

#Richness vs simpson
plot(big_df_rare.richness, big_df_rare.simpson)
cor.test(big_df_rare.richness, big_df_rare.simpson, method = "spearman")

#Richness vs shannon

plot(big_df_rare.richness, big_df_rare.shannon)
cor.test(big_df_rare.richness, big_df_rare.shannon, method = "spearman")

#simpson vs shannon

plot(big_df_rare.simpson, big_df_rare.shannon)
cor.test(big_df_rare.simpson, big_df_rare.shannon, method = "spearman")





