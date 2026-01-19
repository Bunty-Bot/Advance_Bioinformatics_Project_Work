#!/usr/bin/env rscript

#--------------------------------------------------------------------------
# Author: Chris Ambagtsheer (student number: 1216414)
# Description: perform differential pathway analysis
# Usage: rscript pathway_analysis.R statistics_table.tsv output title
# Arguments:
#   pathway_analysis.R: this R script
#   statistics_table.tsv: tab delimited table with pathway ID, pathway 
#     description, P-value (adjusted), query size, term size, overlap size,
#     genes and gene ratio. 
#   output: output directory for figures
#   title: title of figures
# Implementation:
#   import libraries for plotting and colour schemes
#   take top significant pathways
#   make dot plot for significant pathways and relevance
#--------------------------------------------------------------------------

# Define a user library
user_lib <- "~/R/library"

# Create the library if it doesn't exist
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
  message("Created user library at: ", user_lib)
}

# Prepend the user library to the library search path
.libPaths(c(user_lib, .libPaths()))

# List of required packages
packages <- c("ggplot2", "RColorBrewer")

# Check which packages are missing in the user library
installed <- packages %in% installed.packages(lib.loc = user_lib)[, "Package"]

# Install missing packages into the user library
if (any(!installed)) {
  message("Installing missing packages: ", paste(packages[!installed], collapse = ", "))
  install.packages(packages[!installed], lib = user_lib, repos = "https://cloud.r-project.org")
}





# Load all required packages
lapply(packages, library, character.only = TRUE)
args = commandArgs(trailingOnly = TRUE)
df = read.delim(args[1])
output = args[2]
title = args[3]

df = df[order(df$gene_ratio, decreasing = TRUE),]
df = df[df$padj < 0.05,]

df$description <- factor(df$description, 
                             levels = rev(unique(df$description)))

##Creating a figure

#Build output file 
outfile <- file.path(output, paste0(title, "_dotplot.pdf"))

#Open PDF device
pdf(outfile, width = 12, height = 6)

##Define a palette
mypalette <- brewer.pal(3, "Blues")

## Visualizing data with ggplot2
#make a dotplot
ggplot(df) + #an empty plot is created
    
  #the type of plot you want (points/dots) is specified.
  #specify aestetics for this layer
  geom_point(    
    aes(x = gene_ratio, 
        y = description,  
        color = -log10(padj),
        size = pathway_size)) +
  
  #make a theme. change x-axis and title size rel. to default.
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.15)),
        axis.title = element_text(size=rel(1.15))) +
  
  #make labels for the axis
  xlab("Gene ratio") +
  ylab("Top PlantCyc pathways") +
  
  #center the title and make title
  ggtitle(title) +
  theme(plot.title=element_text(hjust=0.5, face = "bold")) +
  
  #define colour gradient colours and name
  scale_color_gradientn(name = "-log10(padj)",
                        colors = mypalette)+
  
  #make legend text bold
  theme(legend.title = element_text(size=rel(1.15),
                                    hjust=0.5, 
                                    face="bold")
)    


dev.off()

  