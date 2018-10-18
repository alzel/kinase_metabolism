Figure 2 results
================
Aleksej Zelezniak
2018-10-18

The Deletion of Each Yeast Kinase Triggers a Unique Reconfiguration of Enzyme Expression in the Cell

    ## Warning: package 'tidyverse' was built under R version 3.4.2

    ## -- Attaching packages --------------------------------------------------- tidyverse 1.2.1 --

    ## <U+221A> ggplot2 2.2.1     <U+221A> purrr   0.2.5
    ## <U+221A> tibble  1.4.2     <U+221A> dplyr   0.7.6
    ## <U+221A> tidyr   0.8.1     <U+221A> stringr 1.3.1
    ## <U+221A> readr   1.1.1     <U+221A> forcats 0.3.0

    ## Warning: package 'tibble' was built under R version 3.4.3

    ## Warning: package 'tidyr' was built under R version 3.4.4

    ## Warning: package 'purrr' was built under R version 3.4.4

    ## Warning: package 'dplyr' was built under R version 3.4.4

    ## Warning: package 'stringr' was built under R version 3.4.4

    ## Warning: package 'forcats' was built under R version 3.4.3

    ## -- Conflicts ------------------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    ## Warning: package 'scales' was built under R version 3.4.1

    ## 
    ## Attaching package: 'scales'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard

    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

    ## Warning: package 'gridExtra' was built under R version 3.4.1

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## Warning: package 'reshape2' was built under R version 3.4.3

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

``` r

load("./data/proteins.matrix.sva.0.5.1.RData")
load("./data/proteins.matrix.sva.0.5.1.FC.RData")
load("./data/iMM904._load_.RData")
load("./data/exp_metadata._clean_.RData")
load("./data/orf2name._clean_.RData")
load("./data/gene.annotations._load_.RData")

load("./data/pathway2orf._load_.RData")
UniProt2Reactome <- read.delim("./data/UniProt2Reactome.txt", header=FALSE, stringsAsFactors = F)

protein.matrix <- proteins.matrix.sva.0.5.1
proteins.FC <- proteins.matrix.sva.0.5.1.FC

iMM904[] <- lapply(iMM904, as.character)


kinase_orfs <- unique(as.character(exp_metadata$ORF[exp_metadata$type == "Kinase"]))
uniprot2orf <- gene.annotations %>%  
  filter(V3 ==  "UniProt/Swiss-Prot ID") %>%
  dplyr::select(V1, V4, V6) %>% distinct() %>%
  rename(uniprot_id = V1, ORF = V4,  gene_name = V6)
uniprot2orf[] <- lapply(uniprot2orf[], as.character)

uniprot2orf.kinases <- uniprot2orf %>% filter(ORF %in% kinase_orfs)
```

### Figure 2A

``` r

reference = unique(as.character(proteins.FC$reference))
pval_thr = 0.01
FC_thr = getFC_thr(proteins.matrix=protein.matrix, pval_thr=pval_thr)

proteins.FC.f <- proteins.FC %>% 
  filter(KO %in% unique(exp_metadata$ORF[exp_metadata$type == "Kinase"])) %>%
  mutate(isiMM904 = ORF %in% unique(as.character(iMM904$gene)))

all_proteins <- as.character(unique(proteins.FC.f$ORF))
all_measured_enzymes <- as.vector((proteins.FC.f %>% filter(isiMM904 ==T ) %>% dplyr::select(ORF) %>% distinct())$ORF)

proteins.FC.f.metabolic <- tbl_df(proteins.FC.f) %>% 
  filter(abs(logFC) >= FC_thr, p.value_BH < pval_thr, isiMM904 == T)

proteins.FC.f.metabolic$KO.gene <- orf2name$gene_name[match(proteins.FC.f.metabolic$KO, orf2name$ORF)]
stopifnot(!any(is.na(proteins.FC.f.metabolic$KO.gene)))

FC.f.metabolic.stats <- proteins.FC.f.metabolic %>% 
  group_by(KO.gene) %>% 
  summarise(n = n(),
            n_pos = sum(logFC > 0)/length(all_measured_enzymes),
            n_neg = sum(logFC < 0)/length(all_measured_enzymes)) %>% 
  ungroup() %>% arrange(n)


# all enzymes

x = proteins.FC.f.metabolic
x.wide <- dcast(x, "KO.gene ~ ORF", value.var = "logFC")
x.wide[is.na(x.wide)] <- 0
x.wide.matrix <- x.wide[,-1]
rownames(x.wide.matrix) <- x.wide[,1]
x.wide.matrix <- ifelse(x.wide.matrix != 0, 1, 0)

d.matrix.all <- 1 - as.matrix(dist(x.wide.matrix, method = "binary"))


#upregulated
x = proteins.FC.f.metabolic %>% filter(logFC > 0)
x.wide <- dcast(x, "KO.gene ~ ORF", value.var = "logFC")
x.wide[is.na(x.wide)] <- 0
x.wide.matrix <- x.wide[,-1]
rownames(x.wide.matrix) <- x.wide[,1]
x.wide.matrix <- ifelse(x.wide.matrix != 0, 1, 0)

d.matrix.up <- 1 - as.matrix(dist(x.wide.matrix, method = "binary"))

#downregulated
x = proteins.FC.f.metabolic %>% filter(logFC < 0)
x.wide <- dcast(x, "KO.gene ~ ORF", value.var = "logFC")
x.wide[is.na(x.wide)] <- 0
x.wide.matrix <- x.wide[,-1]
rownames(x.wide.matrix) <- x.wide[,1]
x.wide.matrix <- ifelse(x.wide.matrix != 0, 1, 0)

d.matrix.down <- 1 - as.matrix(dist(x.wide.matrix, method = "binary"))


cl = hclust(dist(d.matrix.all))
cl <- dendextend::rotate(cl, order = as.character(FC.f.metabolic.stats$KO.gene))
d.matrix.all <- d.matrix.all[cl$order,cl$order]


d.matrix.up <- d.matrix.up[rownames(d.matrix.up)[match(rownames(d.matrix.all),rownames(d.matrix.up))], 
                           colnames(d.matrix.up)[match(colnames(d.matrix.all),colnames(d.matrix.up))]]
d.matrix.down <- d.matrix.down[rownames(d.matrix.down)[match(rownames(d.matrix.all),rownames(d.matrix.down))], 
                           colnames(d.matrix.down)[match(colnames(d.matrix.all),colnames(d.matrix.down))]]


zeros.matrix <- matrix(data=0, ncol = ncol(d.matrix.down), nrow = nrow(d.matrix.up))
zeros.matrix[upper.tri(zeros.matrix)] <- d.matrix.up[upper.tri(d.matrix.up)]
zeros.matrix[lower.tri(zeros.matrix)] <- d.matrix.down[lower.tri(d.matrix.down)]*-1





toPlot <- melt(zeros.matrix)
toPlot$x.name <- factor(rownames(d.matrix.all)[toPlot$Var1], levels = rownames(d.matrix.all))
toPlot$y.name <- factor(colnames(d.matrix.all)[toPlot$Var2], levels = colnames(d.matrix.all))

my_breaks <- seq(1, -1, -0.25)
my_colours <- rev(brewer.pal(name = "RdBu", n = length(my_breaks) - 1))
my_colours[c(4,5)] <- "white"
#pheatmap(d.matrix.all, cluster_rows = F, cluster_cols = F)


ggplot(toPlot) +  
  geom_tile(aes(x = x.name, y = y.name, fill = cut(value, breaks = my_breaks)), colour="grey") +
#   scale_fill_gradient2(low="#1F78B4",high="#E31A1C",mid ="white",
#                        breaks = seq(-0.75, 0.75, 0.25),
#                        midpoint=0)  +
  scale_fill_manual(values = my_colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "italic"), 
        aspect.ratio = 1, legend.position = c(0.2, 0.8) ) +
  labs(x="", y = "")
```

<embed src="Figure2_files/figure-markdown_github/heatmap-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
### Figure 2B

``` r

combinations <- combn(unique(proteins.FC.f.metabolic$KO), 2)
combinations.df <- as.data.frame(t(combinations))

KO.genes <- proteins.FC.f.metabolic %>% dplyr::select(KO, ORF) %>% distinct()


overlaps <- plyr::ddply(combinations.df, plyr::.(V1, V2),
                  
                  .fun = function(x) {
                    df1 <- KO.genes %>% filter(KO == x$V1)
                    df2 <- KO.genes %>% filter(KO == x$V2)
                    
                    result <- bind_rows(df1, df2) %>% 
                      group_by(ORF) %>%
                      summarise(gene_n = n()) %>% 
                      summarize(intersection = sum(gene_n == 2),
                                union = n(),
                                overlap = intersection/union)
                    
                    
                    return(result)
                  }  )


overlaps.stats <- overlaps %>% group_by(V1) %>%
  summarize(mean.intersection = mean(intersection, na.rm=T)) %>%
  mutate(gene_name = orf2name$gene_name[match(V1, orf2name$ORF)]) %>%
  left_join(FC.f.metabolic.stats, by = c("gene_name" = "KO.gene"))

toPlot <- FC.f.metabolic.stats %>% 
  dplyr::select(KO.gene, n_pos, n_neg) %>% as.data.frame()%>%
  reshape2::melt(id.vars = "KO.gene") 
toPlot$KO.gene <- factor(toPlot$KO.gene, levels = rownames(d.matrix.all)) 

toPlot.line <- overlaps.stats %>% 
  ungroup() %>% 
  mutate(n_fraction = n/length(all_measured_enzymes),
         mean.intersection_fraction = mean.intersection/length(all_measured_enzymes))

ggplot() + 
  geom_bar(data = toPlot, stat="identity", width=.5, aes(x=KO.gene, y=value, fill=variable)) + 
  geom_line(data = toPlot.line, aes(x = gene_name, y=mean.intersection_fraction, group = 1)) +
  labs(x = "", y = "Fraction of perturbed metabolic network") +
  scale_fill_manual(values = my_colours[c(length(my_colours),1)]) +
  theme_few() +
  theme(legend.position = c(0.1, 0.7),
        axis.text.x = element_text(angle = 90, hjust = 1))
```

<embed src="Figure2_files/figure-markdown_github/barplot-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
### Related Figure 2C and Figure 2D

All the rest plot related to C and D panels as well to Figure S8 amd Figure S9

``` r
#data prep for similarities 
similarities.long <- list()
tmp <- melt(d.matrix.up)
tmp$type = "up"
tmp$dataset = "metabolic"
similarities.long <- lappend(similarities.long, tmp)


tmp <- melt(d.matrix.down)
tmp$type = "down"
tmp$dataset = "metabolic"
similarities.long <- lappend(similarities.long, tmp)

similarities.long <- lapply(similarities.long, 
       FUN = function(x) {
         x %>% mutate(Var1 = orf2name$ORF[match(Var1, orf2name$gene_name)],
                      Var2 = orf2name$ORF[match(Var2, orf2name$gene_name)])
  
       })



proteins.FC.f = proteins.FC %>% filter(KO %in% kinase_orfs)
proteins.FC.f.wide <- dcast(proteins.FC.f, formula = "ORF~KO", value.var = "logFC")
proteins.FC.f.matrix <- as.matrix(proteins.FC.f.wide[,-1])
rownames(proteins.FC.f.matrix) <- proteins.FC.f.wide$ORF

tmp.cor <- cor(proteins.FC.f.matrix)
diag(tmp.cor) = NA
proteins.FC.f.matrix.cor.long <- melt(tmp.cor, varnames = c("Var1", "Var2"))
similarities.long <- lappend(similarities.long, 
                             proteins.FC.f.matrix.cor.long %>% 
                               mutate(type = "pearson_fc",
                                      dataset = "metabolic"))


my_means <- function(proteins.matrix) {
  
  proteins.long = melt(proteins.matrix, id.vars="rownames")
  names(proteins.long) = c("EG.StrippedSequence", "R.Label", "signal")
  proteins.long$ORF = exp_metadata$ORF[match(proteins.long$R.Label, exp_metadata$sample_name)]
  proteins.long.mean = tbl_df(proteins.long) %>% group_by(EG.StrippedSequence, ORF) %>% summarize(mean = mean(signal))
  proteins.mean.df = dcast(proteins.long.mean, formula=EG.StrippedSequence~ORF, value.var="mean")
  
  proteins.mean.matrix = as.matrix(proteins.mean.df[,-1])
  rownames(proteins.mean.matrix) = as.matrix(proteins.mean.df$EG.StrippedSequence)
  return(proteins.mean.matrix)  
}

protein.matrix.mean = my_means(exp(protein.matrix))

tmp.cor <- cor(protein.matrix.mean[,kinase_orfs])
diag(tmp.cor) = NA
protein.matrix.mean.cor.long <- melt(tmp.cor, varnames = c("Var1", "Var2"))

similarities.long <- lappend(similarities.long, 
                             protein.matrix.mean.cor.long %>% mutate(type = "pearson",
                                                                    dataset = "metabolic"))

similarities_dataset <- bind_rows(similarities.long) %>% 
  filter(Var1 != Var2) %>% select(-dataset)

similarities_dataset <- similarities_dataset %>% 
  mutate(X1_uniprot = uniprot2orf$uniprot_id[match(Var1, uniprot2orf$ORF)],
         X2_uniprot = uniprot2orf$uniprot_id[match(Var2, uniprot2orf$ORF)])


similarities_dataset %>%
  ggplot(aes(x = value)) +
    geom_histogram() +
    facet_wrap(~type, scales = "free")
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

<embed src="Figure2_files/figure-markdown_github/unnamed-chunk-1-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
Pathway co-expression
=====================

clearning data for similarity histograms

``` r


names(UniProt2Reactome) <- c("uniprot_id", "reactome_id", "url", "description", "confidence", "species")
UniProt2Reactome <- tbl_df(UniProt2Reactome) %>% 
  filter(species == "Saccharomyces cerevisiae") %>% 
  group_by(reactome_id, uniprot_id) %>% distinct()


uniprot2reactome <- UniProt2Reactome %>% 
  dplyr::filter(uniprot_id %in% uniprot2orf.kinases$uniprot_id) %>% 
  group_by(reactome_id) %>%
  mutate(n = n()) %>% 
  filter(n >2) %>% 
  arrange(reactome_id) %>%
  rename(pathway = reactome_id) %>%
  select(pathway, uniprot_id, n)



uniprot2kegg <- pathway2orf %>% 
  mutate(uniprot_id = uniprot2orf.kinases$uniprot_id[match(ORF, uniprot2orf.kinases$ORF)]) %>%
  filter(!is.na(uniprot_id)) %>%
  group_by(pathway) %>% 
  mutate(n = n()) %>% 
  filter(n >2) %>% arrange(pathway) %>% 
  dplyr::select(pathway, uniprot_id, n)

# in KEGG 4 pathways with more than 2 kinases are mapped 37 kinases
# in Reactome 38 pathways with more than 2 kinases are mapped to 29 kinases 

uniprot2kegg <- uniprot2kegg %>% mutate(pathway_base = "kegg")
uniprot2reactome <- uniprot2reactome %>% mutate(pathway_base = "reactome")

pathway_kinases <- bind_rows(uniprot2kegg, uniprot2reactome)
```

-   In KEGG there are 4 pathways with more than 2 kinase members, in total mapped to 37 kinases
-   In Reactome there are 38 pathways with more than 2 kinase members, in total mapped to 29 kinases
-   In total 49 kinases are mapped to pathways

``` r

getPathwaySimilarity <- function(pathways) {
  ret = plyr::ddply(pathways, plyr::.(pathway, pathway_base), 
        .fun = function(x) {
          return.list = list()
          tmp <<- data.frame(t(combn(unique(as.character(x$uniprot_id)), 2)))
          tmp$n = x$n[1]
          tmp.long <- left_join(tmp, similarities_dataset, by = c("X1" = "X1_uniprot", "X2" = "X2_uniprot"))
          return(tmp.long)
        })
  return(ret)
  
}

set.seed(123)
pathway_kinases_random <- plyr::ddply(pathway_kinases, 
                                    plyr::.(pathway_base),
                                    .fun = function(x) {
                                  # x = pathway_kinases %>% filter(pathway_base == "kegg")
                                  uniprots = unique((pathway_kinases %>% 
                                                    filter(pathway_base == unique(x$pathway_base)) %>%
                                                      dplyr::select(uniprot_id))$uniprot_id)
                                  z <<- plyr::ddply(x, 
                                        plyr::.(pathway), 
                                        .fun = function(y) {
                          #                 y = pathway_kinases %>% filter(pathway == "path:sce04011",
                          #                                                pathway_base == "kegg")
                                          
                                            n = unique(y$n)
                                            tmp = as.matrix(replicate(10, sample(uniprots, n)))
                                            colnames(tmp) = paste0("rep_", 1:ncol(tmp), sep="")
                                            tmp.long <- melt(tmp, varnames = c("X1", "X2"))
                                            tmp.long$n <- n
                                            return(tmp.long)
                                        }
                                        )
                                  
                                          return(data.frame(pathway = paste0(z$pathway,"_",  z$X2),
                                                           uniprot_id = z$value,
                                                           n = z$n))
                                    })
#> Adding missing grouping variables: `pathway`
#> Adding missing grouping variables: `pathway`



set.seed(123)
pathway_similarities <- getPathwaySimilarity(pathways = pathway_kinases)
pathway_similarities$sample_type = "signal" 
pathway_similarities_random <- getPathwaySimilarity(pathways = pathway_kinases_random)
pathway_similarities_random$sample_type = "random"
pathways_dataset <- bind_rows(pathway_similarities, pathway_similarities_random)

toPlot <- pathways_dataset
toPlot.stats <- toPlot %>% group_by(type, pathway_base) %>% 
    summarise(pval = (wilcox.test(value[sample_type == "signal"], value[sample_type == "random"])$'p.value'))
              
  
toPlot.stats$padj <- p.adjust(toPlot.stats$pval, method = "BH")

toPlot.stats.medians <- toPlot %>% 
    group_by(type, pathway_base, sample_type) %>% 
    summarise(median_value = median(value, na.rm = T))
  
toPlot %>% 
  ggplot(aes(x = value)) +
    geom_density(aes(fill = sample_type), alpha = 0.5) +
    facet_wrap(~pathway_base+type, scales = "free", ncol = 3) +
    geom_text(data = toPlot.stats, aes(x=0.5, y = 1, label= paste("p-value=", format(padj, digits=2, scientific=T)))) +
    geom_vline(data = toPlot.stats.medians, aes(xintercept = median_value, colour = sample_type), linetype = 2) +
    # theme_bw() + 
    theme(legend.position = c(0.1, 0.5), aspect.ratio = 5/8)  

p = toPlot %>% filter(type == "pearson_fc") %>%
  ggplot(aes(x = value)) +
    geom_density(aes(fill = sample_type), alpha = 0.5) +
    facet_wrap(~pathway_base+type, scales = "free", ncol = 3) +
    geom_text(data = toPlot.stats %>% filter(type == "pearson_fc"), 
              aes(x=0.5, y = 1, label= paste("p-value=", format(pval, digits=2, scientific=T)))) +
    geom_vline(data = toPlot.stats.medians %>% filter(type == "pearson_fc"), 
               aes(xintercept = median_value, colour = sample_type), linetype = 2) +
    theme_bw() + 
    theme(legend.position = c(0.1, 0.5), aspect.ratio = 5/8)  


file_name = paste("fig2_D", fun_name, sep = ".")
file_path = paste(output_dir, file_name, sep="/")

ggsave(filename = paste(file_path, "pdf", sep = "."), device = "pdf", plot = p) 
#> Saving 6 x 3.71 in image
```

<embed src="Figure2_files/figure-markdown_github/unnamed-chunk-2-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
``` r
pathway_similarities$sorted_pair <- apply(pathway_similarities[ ,c("X1", "X2")] , 1 , 
                                          FUN = function(x) {
                                            paste(sort(c(x[1], x[2])), collapse = "|")
                                          })

pathway_similarities %>% filter(type == "pearson_fc") %>%
  group_by(pathway_base, sorted_pair) %>% 
  summarise(n = n()) %>% arrange(-n) %>%
  ggplot(aes(x = n)) +
    geom_histogram() +
    facet_wrap(~pathway_base, scales = "free")
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.


toPlot <- pathways_dataset %>% filter(!(X1 %in% c("P12688", "P18961")), !(X2 %in% c("P12688", "P18961")))
toPlot.stats <- toPlot %>% group_by(type, pathway_base) %>% 
    summarise(pval = (wilcox.test(value[sample_type == "signal"], value[sample_type == "random"])$'p.value'))
              
toPlot.stats$padj <- p.adjust(toPlot.stats$pval, method = "BH")

toPlot.stats.medians <- toPlot %>% 
    group_by(type, pathway_base, sample_type) %>% 
    summarise(median_value = median(value, na.rm = T))
  
toPlot %>% 
  ggplot(aes(x = value)) +
    geom_density(aes(fill = sample_type), alpha = 0.5) +
    facet_wrap(~pathway_base+type, scales = "free", ncol = 3) +
    geom_text(data = toPlot.stats, aes(x=0.5, y = 1, label= paste("p-value=", format(padj, digits=2, scientific=T)))) +
    geom_vline(data = toPlot.stats.medians, aes(xintercept = median_value, colour = sample_type), linetype = 2) +
    # theme_bw() + 
    theme(legend.position = c(0.1, 0.5), aspect.ratio = 5/8)  
```

<embed src="Figure2_files/figure-markdown_github/no_frequent_pairs-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
<embed src="Figure2_files/figure-markdown_github/no_frequent_pairs-2.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
``` r
sessionInfo()
#> R version 3.4.0 (2017-04-21)
#> Platform: x86_64-apple-darwin15.6.0 (64-bit)
#> Running under: macOS  10.13.5
#> 
#> Matrix products: default
#> BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] C
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] bindrcpp_0.2.2     cluster_2.0.6      ggthemes_3.4.0    
#>  [4] RColorBrewer_1.1-2 reshape2_1.4.3     gridExtra_2.3     
#>  [7] scales_0.5.0       forcats_0.3.0      stringr_1.3.1     
#> [10] dplyr_0.7.6        purrr_0.2.5        readr_1.1.1       
#> [13] tidyr_0.8.1        tibble_1.4.2       ggplot2_2.2.1     
#> [16] tidyverse_1.2.1   
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_0.2.4 haven_1.1.1      lattice_0.20-35  colorspace_1.3-2
#>  [5] htmltools_0.3.6  yaml_2.2.0       rlang_0.2.2      pillar_1.2.1    
#>  [9] foreign_0.8-69   glue_1.3.0       modelr_0.1.1     readxl_1.0.0    
#> [13] bindr_0.1.1      plyr_1.8.4       munsell_0.4.3    gtable_0.2.0    
#> [17] cellranger_1.1.0 rvest_0.3.2      psych_1.8.4      evaluate_0.10.1 
#> [21] knitr_1.20       parallel_3.4.0   broom_0.4.4      Rcpp_0.12.18    
#> [25] backports_1.1.2  jsonlite_1.5     mnormt_1.5-5     hms_0.4.1       
#> [29] digest_0.6.15    stringi_1.2.2    grid_3.4.0       rprojroot_1.3-2 
#> [33] cli_1.0.0        tools_3.4.0      magrittr_1.5     lazyeval_0.2.1  
#> [37] crayon_1.3.4     pkgconfig_2.0.1  xml2_1.2.0       lubridate_1.7.4 
#> [41] assertthat_0.2.0 rmarkdown_1.9    httr_1.3.1       rstudioapi_0.7  
#> [45] R6_2.2.2         nlme_3.1-131.1   compiler_3.4.0
```
