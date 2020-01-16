##############################
## Code review script - candidate genes
##
## Matt Brachmann (PhDMattyB)
##
## 2020-01-10
##
##############################

setwd('~/PhD/R users group/Gene_region')

library(patchwork)
library(janitor)
library(devtools)
library(skimr)
library(rsed)
library(data.table)
library(sjPlot)
library(tidyverse)

theme_set(theme_bw())

# Other packages to load


# Genome data -------------------------------------------------------------

gff = read_tsv('ref_ASM291031v2_top_level.gff3',
               ## read in the Arctic charr genome data (gff file)
               col_names = F, 
               skip = 9) %>%
  ## there are 9 garbage lines, need to get rid of them
  rename(chr = X1,  
         ## Need to rename the columns to something informative
         source = X2, 
         feature = X3, 
         start = X4, 
         end = X5, 
         score = X6, 
         strand = X7, 
         frame = X8, 
         attribute = X9) %>% 
  filter(feature == 'gene') %>% 
  ## only want the coding gene regions
  arrange(start, 
          end) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = start + (end-start)/2)
## calculate the mid point of each gene from the start and end points


# Our SNP data ------------------------------------------------------------
data = read_tsv('PCAdapt_lake_outliers.txt') %>% 
  ## read in the dataframe of outlier loci within the population
  arrange(CHROME3) %>% 
  ## Arrange the df by chromosome number
  filter(CHROME3 != 'Contigs')
## filter out the unplace contigs

## If you want to view the data:
## View(data)


# gene regions ------------------------------------------------------------
## lets look for gene regions associated with the outlier
## snps on chromosome 8

AC08 = gff %>% 
  filter(chr == 'NC_036848.1') 
## filter out chromosome 8 in the gff file of the Arctic charr genome

AC08_out = data %>% 
  filter(CHROME3 == 'AC08') 
## filter out chromosome 8 in out data set

pos = AC08_out$PDIST.x
## grab the position of outlier snps on chromosome 8

gene_regions = AC08 %>% 
  mutate(hit_dist = abs(mid - pos)) %>% 
  ## find the distance between the mid point and the snp position
  arrange(hit_dist) %>% 
  ## arrange by hit distance
  filter(hit_dist <5000) %>%
  ## lets find genes within 5kb of our snps
  select(chr, 
         start, 
         end, 
         attribute, 
         hit_dist) %>% 
  ## select the data that we need
  pull(attribute)
## pull gets all of the hidden data from the attribute column
## this is the important column, 
## it'll give us the gene ID and Names for the genes on the 
## Arctic charr genome

## to find the actual gene names go to:
## gene 1: https://www.ncbi.nlm.nih.gov/gene/111967237
## gene 2: https://www.ncbi.nlm.nih.gov/gene/111967356

## to find out what they do and what they're associated with:
## gene 1 is the GPC4 gene: https://www.genecards.org/cgi-bin/carddisp.pl?gene=GPC4
## gene 2 is the GIMAP7 gene: https://www.genecards.org/cgi-bin/carddisp.pl?gene=GIMAP7


