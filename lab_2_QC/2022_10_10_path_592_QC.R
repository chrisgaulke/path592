# Title and author information --------------------------------------------
#!/usr/bin/R

#################################################
#                                               #
#          2022_10_10_path_592_QC.R             #
#                                               #
#################################################

#Title: Path592 QC for zfish biogeography
#
#Copyright (C) 2022-2023  Christopher A. Gaulke
#author contact: chris.gaulke@gmail.com
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#For a copy of the GNU General Public License see
#<http://www.gnu.org/licenses/>.

#Purpose: The purpose of this script is to demonstrate QC steps availble in
# Dada2

# PACKAGE INSTALLATION ----------------------------------------------------

# If you haven't already, install the required software. The commands in this
# section are run only once and should be commented after you have finished running
# them

# Installing Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

# Installing dada2
BiocManager::install("dada2")

# Installing other packages needed
install.packages("ggplot2")
install.packages("vegan")

# Unpacking and sorting files ---------------------------------------------

# The following commands can be run in a terminal window (Mac) or in Rstudio console
# window. The tar archive should be unpacked in a suitable location
# (i.e., not the desktop). Preferably you have a robust file naming system for
# your research. For example, my research analysis are all kept in
# ~/Documents/research/<name_of_analysis_project>. Each project has the following
# subdirectories: data/, analysis, scripts/. Scripts contains the R project. Other
# directories can be added as needed (e.g., manuscript, reports, ppt, etc.). This
# compartmentalizes the projects and helps keep things organized


#decompress the and untar the archive (Mac)
#tar -xvf path_592_files_biogeo.tar.gz

# SET ENVIRONMENT ---------------------------------------------------------

#load required libraries

library(dada2)
library(ggplot2)
#library(vegan)
#library(reshape2)
#library(pals)
options(stringsAsFactors = F)


# IMPORT DATA -------------------------------------------------------------

path <- "/Users/cgaulke/unsynced_projects/raw_data/2022_10_10_path592_raw/"
filt.path <- "/Users/cgaulke/Documents/research/path592/lab_2_QC/filtered_data/" #filtered file directory make sure to update

#git ignore to ignore this
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# # ANALYSIS: QC ------------------------------------------------------------

#make sure the lengths are the same
length(fnFs) == length(fnRs)

#get sample names, in our case we have to leave some extra on the end
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

#preview
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])


#aggregate all data together now
fnFs_qual.plot <- plotQualityProfile(fnFs,aggregate = T)
fnRs_qual.plot <- plotQualityProfile(fnRs,aggregate = T)

# To adapt code you kinda have to understand what is going on... I can confirm
# my suspicions by profiling the code
plotQualityProfile

#since we have confirmed that they are just wrapping ggplot2 here we can extend
#their code using ggplot2 functions

fnFs_qual.plot + geom_hline(yintercept= 25)
fnFs_qual.plot + geom_hline(yintercept= 25) + geom_vline(xintercept = 250)

fnRs_qual.plot + geom_hline(yintercept= 25)
fnRs_qual.plot + geom_hline(yintercept= 25) + geom_vline(xintercept = 200)

#set up for filtering
filtFs <- file.path(filt.path, "filtered",
                    paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(filt.path, "filtered",
                    paste0(sample.names, "_R_filt.fastq.gz"))

#make these named
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filter and trim
# by looking at the data we see a rapid drop in quality around 250bp R1 (200 R2).
# Since the average drops below ~30 around 250 we will truncate at 200 for the
# reverse. The forward looks better (this is usual) so we will truncate around
# 250. We will also take off about 10 bases on the left as these bases are
# highly skewed.

#note the original tutorial uses generic variable names

filter.out <- filterAndTrim(fnFs, #paths to the input forward reads
                            filtFs, #paths to the output filtered forward reads
                            fnRs, #paths to the input reverse reads
                            filtRs, #paths to the output filtered reverse reads
                            truncLen=c(250,200), #R1 and R2 truncation lengths
                            maxN=0, #max allowable N's (ambiguous bases) in seq
                            maxEE=c(2,2), # max error allowed after filtering
                            truncQ=2, # truncate after first base with this score
                            rm.phix=TRUE, # remove phix reads
                            trimLeft = 10, # number of nt to trim from left
                            compress=TRUE, # gzip files
                            multithread=TRUE # use multiple threads. Turn off on
                                             # windows
                            )
#take a look
View(filter.out)

colMeans(filter.out) #mean number of reads in and out of filtering
mean(1-(filter.out[,2]/filter.out[,1])) #mean % filtered reads

#lets look at these numbers in a little more detail
fivenum(1-(filter.out[,2]/filter.out[,1])) #five number report for % filtered
hist(1-(filter.out[,2]/filter.out[,1]))  # hist #mean % filtered reads

# since some libraries look like there is a higher level of filtration
# than others lets take a closer look at this

sort(1-(filter.out[,2]/filter.out[,1]))

#let's keep this in mind moving forward

