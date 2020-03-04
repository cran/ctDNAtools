## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE
)

## ----setup---------------------------------------------------------------
library(ctDNAtools)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)

## Load example data

## example list of predetected mutations
data('mutations',package = 'ctDNAtools')
mutations

## example target file
data('targets', package = 'ctDNAtools')
targets

## Example bam files for middle-treatment and after-treatment
bamT1 <- system.file('extdata', 'T1.bam', package = 'ctDNAtools')
bamT2 <- system.file('extdata', 'T2.bam', package = 'ctDNAtools')

## Eample bam files from samples taken from healthy subjects
bamN1 <- system.file('extdata', 'N1.bam', package = 'ctDNAtools')
bamN2 <- system.file('extdata', 'N2.bam', package = 'ctDNAtools')
bamN3 <- system.file('extdata', 'N3.bam', package = 'ctDNAtools')

## Reference genome
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

## ----basic1--------------------------------------------------------------
test1 <- test_ctDNA(mutations = mutations,
          bam = bamT1,
          reference = BSgenome.Hsapiens.UCSC.hg19, 
          targets = targets,
          informative_reads_threshold = 100)
test1

## ----basic2--------------------------------------------------------------
test2 <- test_ctDNA(mutations = mutations, 
          bam = bamT2,
          reference = BSgenome.Hsapiens.UCSC.hg19,
          targets = targets,
          informative_reads_threshold = 100,
          min_base_quality = 20,
          min_mapq = 30)
test2

## batch runs
## use future_map2_dfr for multi-threading

tests <- map2_dfr(c(bamT1, bamT2),
            list(mutations, mutations), # in case mutations are different
            ~ test_ctDNA(bam = .x, mutations = .y,
                         targets = targets,
                         reference = BSgenome.Hsapiens.UCSC.hg19,
                         informative_reads_threshold = 100))
tests


## ----blacklist-----------------------------------------------------------
## Black list by loci (substition_specific = FALSE)
bg_panel1 <- create_background_panel(bam_list = c(bamN1, bamN2, bamN3), 
          targets = targets, 
          reference = BSgenome.Hsapiens.UCSC.hg19, 
          substitution_specific = FALSE)

black_list1 <- create_black_list(bg_panel1,
          mean_vaf_quantile = 0.99,
          min_samples_one_read = 2,
          min_samples_two_reads = 1)

head(black_list1)

test3 <- test_ctDNA(mutations = mutations, bam = bamT1,
          reference = BSgenome.Hsapiens.UCSC.hg19, targets = targets,
          informative_reads_threshold = 100, black_list = black_list1,
          substitution_specific = FALSE)
test3

## Black list by variants (substition_specific = TRUE)
bg_panel2 <- create_background_panel(bam_list = c(bamN1, bamN2, bamN3),
           targets = targets, 
           reference = BSgenome.Hsapiens.UCSC.hg19,
           substitution_specific = TRUE)

black_list2 <- create_black_list(bg_panel2,
          mean_vaf_quantile = 0.99,
          min_samples_one_read = 2,
          min_samples_two_reads = 1)

head(black_list2)

test4 <- test_ctDNA(mutations = mutations,
          bam = bamT1,
          reference = BSgenome.Hsapiens.UCSC.hg19, 
          targets = targets,
          informative_reads_threshold = 100,
          black_list = black_list2,
          substitution_specific = TRUE)

test4

## ----phasing-------------------------------------------------------------
## Exploiting phased variants in the test
test5 <- test_ctDNA(mutations = mutations, 
          bam = bamT1,
          reference = BSgenome.Hsapiens.UCSC.hg19,
          targets = targets,
          informative_reads_threshold = 100, 
          black_list = black_list2,
          substitution_specific = TRUE, 
          ID_column = "PHASING")
test5

## ----fragments1----------------------------------------------------------

fs1 <- get_fragment_size(bam = bamT1, 
      mapqFilter = 30, 
      isProperPair = NA, 
      min_size = 1, 
      max_size = 400, 
      ignore_trimmed = FALSE, 
      simple_cigar = FALSE, 
      different_strands = TRUE)
head(fs1)

## optional target restriction
fs2 <- get_fragment_size(bam = bamT1, 
      targets = targets)
head(fs2)

## ----fragments2----------------------------------------------------------

fs3 <- get_fragment_size(bam = bamT1, 
      mutations = mutations)
head(fs3)


## ----histograms----------------------------------------------------------
bfs1 <- bin_fragment_size(bam = bamT1,
          min_size = 1, 
          max_size = 400, 
          normalized = TRUE,
          bin_size = 5)

head(bfs1)

## batch execution
## you can use multithreading by using
## furrr::future_map instead of map
bfs <- c(bamT1, bamT2, bamN1, bamN2, bamN3) %>%
  map(bin_fragment_size, bin_size = 5, normalized = TRUE) %>%
  purrr::reduce(inner_join, by = "Breaks")

head(bfs)

bfs %>%
  tidyr::pivot_longer(cols = -Breaks, names_to = "Sample",values_to = "Counts") %>%
  tidyr::separate(Breaks, into = c("start", "end"), sep = "_", convert = T) %>%
  ggplot(aes(x  = start, y = Counts, color = Sample)) +
  geom_line(size = 1) + theme_minimal()


## custom bins
bfs2 <- bin_fragment_size(bam = bamT1,
          normalized = TRUE,
          custom_bins = c(100,200),
          min_size = 1,
          max_size = 400)

bfs2

## restricted targets
bfs2 <- bin_fragment_size(bam = bamT1,
          targets = targets,
          normalized = TRUE,
          custom_bins = c(100,200),
          min_size = 1,
          max_size = 400)

bfs2

## ----profiles------------------------------------------------------------
## create some regions from the targets
 regions <- data.frame(chr = targets$chr,
         start = seq(from = targets$start - 50, to = targets$end + 50, by = 30),
         stringsAsFactors = FALSE) %>%
     mutate(end = start + 30)
     
sfs1 <- summarize_fragment_size(bam = bamT1,
          regions = regions,
          summary_functions = c(Mean = mean, SD = sd))

head(sfs1)

## batch run
sfs <- c(bamT1, bamT2, bamN1, bamN2, bamN3) %>%
  map(summarize_fragment_size, regions = regions, summary_functions = c(Median = median)) %>%
  purrr::reduce(inner_join, by = "Region")

head(sfs)

## ------------------------------------------------------------------------
af <- analyze_fragmentation(bam = bamT1, 
          targets = targets, 
          window_size = 120, 
          step_size = 5, 
          min_size = 120, 
          max_size = 180)

head(af)

af %>% tidyr::pivot_longer(
    cols = c("WPS_adjusted", "n_fragment_ends_adjusted","n_reads"),
    names_to = "Measurement",
    values_to = "value") %>%
  mutate(Measurement = factor(Measurement, 
        levels = c("n_reads", "n_fragment_ends_adjusted", "WPS_adjusted"))) %>%
  ggplot(aes(x = start, y = value, color = Measurement)) +
  geom_line(size = 1) +
  facet_wrap(~Measurement, scales = "free", nrow = 3) +
  theme_linedraw() +
  labs(x = "Chromosome 14", y = "")


