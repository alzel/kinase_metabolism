Figure 3 results
================
Aleksej Zelezniak
2018-10-18

Enzyme Expression Affects Steady-State Metabolism through Redistributing Flux Control

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

``` r

load("./data/exp_metadata._clean_.RData")

#dataset of control coefficients
flux_mca <- read_delim("./data/mca_fluxes_in_kinase_KOs.csv", delim = "\t")
#> Parsed with column specification:
#> cols(
#>   CC = col_double(),
#>   KO = col_character(),
#>   n_changes = col_double(),
#>   parameter = col_character(),
#>   var_type = col_character(),
#>   variable = col_character()
#> )
conc_mca <- read_delim("./data/mca_conc_in_kinase_KOs.csv", delim = "\t")
#> Parsed with column specification:
#> cols(
#>   CC = col_double(),
#>   KO = col_character(),
#>   n_changes = col_double(),
#>   parameter = col_character(),
#>   var_type = col_character(),
#>   variable = col_character()
#> )

# steady-state concentrations
dataset.ss <- read_delim("./data/ss_fluxes_concentrations_in_kinase_KOs.csv", delim = "\t")
#> Parsed with column specification:
#> cols(
#>   variable = col_character(),
#>   type = col_character(),
#>   value = col_double(),
#>   KO = col_character(),
#>   ss_status = col_double(),
#>   n_changes = col_integer()
#> )

# metabolite concentration daata
load("./data/dataPPP_AA.create_datasets.RData")

dataset.mca <- bind_rows(flux_mca, conc_mca) %>% dplyr::rename(VAR = variable) %>% 
  group_by(var_type, KO, VAR, parameter) %>% distinct(.keep_all = T) %>% ungroup()
dataset.mca$KO_name <- as.character(exp_metadata$gene[match(dataset.mca$KO, exp_metadata$ORF)])
dataset.mca <- dataset.mca %>% mutate(CC = ifelse(CC < 10e-06, 0 , CC))


# 
# dataset.ss <- dataset.ss %>% group_by(variable, type) %>% 
#   filter(ss_status < 10-6) %>% 
#   mutate(z_value = (value - mean(value, na.rm=T))/sd(value, na.rm=T))



#don't format the following vectors!!!
parameter_order <- c('HXK2
HXK1
GLK1
PGM1
PGI1
PFK1
PFK2
FBA1
GPD1
GPD2
HOR2
PGK1
GPM1
ENO1
ENO2
CDC19
PDC1
PDC5
PDC6
ADH1
ADH5') %>% str_split("\n") %>% unlist()

reaction_order <- c('HXT
HXK_HXK2
HXK_HXK1
HXK_GLK1
PGM
TPS
TPP
PGI
PFK
FBA
TPI
GPD
GPP
TDH_TDH1
TDH_TDH3
PGK
GPM
ENO_ENO1
ENO_ENO2
PYK_CDC19
ATPase
PDC_PDC1
PDC_PDC5
PDC_PDC6
ADH_ADH1
acetate_branch
UGP
udp_to_utp
AK') %>% str_split("\n") %>% unlist()

metabolite_order <- c('GLC
G6P
F6P
F16bP
GAP
DHAP
BPG
P3G
P2G
PEP
PYR
AcAld
G1P
G3P
UDP
UTP
T6P
ADP
ATP
NAD') %>% str_split("\n") %>% unlist()
```

Figure 3A
=========

``` r


total_kinases = nrow(dataset.mca %>% ungroup() %>% distinct(KO)) - 1
mca_overal <- dataset.mca %>% group_by(var_type, parameter, KO, KO_name) %>% filter( VAR != "AK") %>% dplyr::summarise(overal_CC = sqrt(sum(CC^2, na.rm = T)))

WT.mca_overal <- mca_overal %>% filter(KO == "WT")

mca_overal <- left_join(mca_overal, WT.mca_overal, by = c("var_type" = "var_type", "parameter" = "parameter"))  %>% ungroup() %>%
  mutate(value_change = (overal_CC.x - overal_CC.y)/overal_CC.y) %>% dplyr::select(-matches(".y$")) 
names(mca_overal) <- sub(x = names(mca_overal), pattern = ".x$", replacement = "")

thr = 0.5
value <- mca_overal %>% filter(KO != "WT") %>% 
  filter(abs(value_change) > thr) %>%
  ungroup() %>%
  distinct(KO, var_type) %>% 
  group_by(var_type) %>%
  dplyr::summarise(n = n()/total_kinases)


mca_overal %>%
  ggplot(aes(x = value_change)) +
    geom_histogram(bins = 60) +
    xlim(-1,2) +
    xlab("Change of enzymes overall control of the fluxes/concentrations\nin kinase mutant in contrast to WT") +
    facet_wrap(~var_type, scales = "free") +
    theme_bw(base_family = "Helvetica")
```

<embed src="Figure3_files/figure-markdown_github/overal_FCC-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
Figure 3A(insets)
=================

``` r
## -- Steady-concentration fluxes ----

WT.ss <- dataset.ss %>% filter(KO == "WT")

dataset.ss <- left_join(dataset.ss, WT.ss, by = c("type" = "type", "variable" = "variable"))  %>% 
  mutate(value_change = (value.x - value.y)/value.y) %>% dplyr::select(-matches(".y$")) 
names(dataset.ss) <- sub(x = names(dataset.ss), pattern = ".x$", replacement = "")


dataset.ss <- dataset.ss %>% group_by(variable, type) %>% 
  filter(ss_status < 10-6) %>% mutate(z_value = (value - mean(value, na.rm=T))/sd(value, na.rm=T))

dataset.ss %>% 
  mutate(abs_value_change = ifelse(abs(value_change)>3,3, abs(value_change))) %>%
  ggplot(aes(x = abs_value_change)) +
    geom_density(fill = "black") +
    scale_x_continuous(limits = c(0,3), labels=percent) +
    facet_wrap(~type, scales = "free") +
    labs(x = "Absolute value change compared to WT, %") +
    theme_bw(base_family = "Helvetica") +
    theme(aspect.ratio = 5/8) +
    theme(panel.grid = element_blank())
```

<embed src="Figure3_files/figure-markdown_github/steady-state-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
Figure 3B
=========

``` r
selected_var = "ADH_ADH1"
toPlot <- dataset.mca %>% ungroup() %>%
  filter(VAR == selected_var, var_type == "flux") %>%
  mutate(parameter = factor(parameter, levels = parameter_order))

toPlot <- toPlot %>% 
  group_by(parameter, VAR, var_type) %>%
    dplyr::mutate(CC_scaled = (CC - mean(CC, na.rm = T))/sd(CC, na.rm = T))
toPlot %>%
  ggplot(aes(x = parameter, y = CC)) +
    ylab(paste("Flux control coeficients of", selected_var)) +
    geom_jitter(alpha = 0.25) +
    #ylim(0, 1) +
    geom_point(data = toPlot %>% filter(KO == "WT"), colour = "red", size = 2) +
    theme_bw() +
    theme(legend.position="none")
```

<embed src="Figure3_files/figure-markdown_github/ADH_ADH1-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
Figure 3C
=========

``` r

#small helper to tidy data for pca
tidy_pca = function(x) {
  x.matrix <- x[,-1] %>% as.matrix()
  rownames(x.matrix) <- as.data.frame(x)[,1]
  return(x.matrix)
}


CCflux_matrix <- dataset.mca %>% filter(var_type == "flux", VAR != "AK", !is.na(CC)) %>% 
  dplyr::select(KO, parameter, VAR, CC) %>%  unite(parameter_VAR, parameter, VAR) %>%
  spread(parameter_VAR, CC) %>% tidy_pca()


CCfluxes_pca <- prcomp(CCflux_matrix)


cl <- pam(CCflux_matrix, 4)

toPlot <- CCfluxes_pca

xlabel <- paste("PC1", round(toPlot$sdev[1]/sum(toPlot$sdev),2))
ylabel <- paste("PC2", round(toPlot$sdev[2]/sum(toPlot$sdev),2))


p.pca <- toPlot$x %>% as_tibble() %>% mutate(gene_name = exp_metadata$gene[match(rownames(toPlot$x), exp_metadata$ORF)],
                                                         cluster = cl$clustering[match(rownames(toPlot$x), names(cl$clustering))]) %>% 
  ggplot(aes(x = PC1 , y = PC2)) +
  geom_point(aes(colour = factor(cluster))) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggrepel::geom_text_repel(aes(label = gene_name, x = PC1, y = PC2), size = 3,  segment.alpha = 0.25, segment.size = 0.25) +
  theme(aspect.ratio = 1) +
  xlab(label = xlabel) +
  ylab(label = ylabel) +
  theme_bw(base_family = "Helvetica") +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        aspect.ratio = 1)


sel <- toPlot %>% .$rotation %>% reshape2::melt() %>%
  setNames(c("variable", "component", "value")) %>% 
  filter(component %in% c("PC1", "PC2")) %>%
  ungroup %>% arrange(component, desc(abs(value))) %>% filter(row_number() == 1) %>% dplyr::select(variable)


p.strongest <- dataset.mca %>% filter(var_type == "flux", VAR != "AK", !is.na(CC)) %>% 
  dplyr::select(KO, parameter, VAR, CC) %>% 
  unite(parameter_VAR, parameter, VAR) %>% 
  filter(parameter_VAR %in% sel$variable) %>%
  mutate(cluster = cl$clustering[match(KO, names(cl$clustering))]) %>%
  ggplot(aes(fill = factor(cluster), x = CC)) +
    geom_density(alpha = 0.5) +
    xlab("Control coefficent of HXK2 over HXK_GLK1") +
    theme_bw(base_family = "Helvetica") +
    theme(legend.position = "none", 
          panel.grid = element_blank(), 
          aspect.ratio = 1)

grid.arrange(p.pca, p.strongest, ncol=2)
```

<embed src="Figure3_files/figure-markdown_github/PCA_FCC-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
``` r



dataset.ss_conc <- dataset.ss %>% filter(type == "conc", ss_status < 1e-6)
dataset_metabolites <- dataPPP_AA$metabolites %>% reshape2::melt(id.vars = rownames())

names(dataset_metabolites) <- c("KO", "metabolite", "conc")
dataset_metabolites$gene_name <- as.character(exp_metadata$gene[match(dataset_metabolites$KO, exp_metadata$ORF)])


model2metabolite <- data.frame(model_name =     c("GLC", "G6P", "F6P", "F16bP",      "GAP", "DHAP", "BPG", "P3G",      "P2G", "PEP",     "PYR", "AcAld", "G1P", "G3P", "UDP", "UTP", "T6P", "ADP", "ATP", "NAD"),
                              metabolite_name = c(NA, "G6P...F6P", "F6P", "X1.6.FP.", NA, "DHAP", NA, "X2...3.PG", "X2...3.PG", "PEP", "Pyr", NA,      NA,    "G3P",    NA,    NA ,     NA,  "ADP", "ATP", NA), stringsAsFactors = F)


dataset.ss_conc$metabolite <- model2metabolite$metabolite_name[match(dataset.ss_conc$variable, model2metabolite$model_name)]
dataset <- left_join(dataset.ss_conc, dataset_metabolites, by = c("metabolite", "KO"))

dataset <- dataset %>% group_by(variable) %>% 
                          mutate(z_value = (value - mean(value, na.rm = T))/sd(value, na.rm = T),
                                 z_conc = (conc - mean(conc, na.rm = T))/sd(conc, na.rm = T))


toPlot <- dataset %>% 
   dplyr::group_by(KO) %>%
   filter(!is.na(z_conc), !is.na(z_value)) %>%
   mutate(cor = cor(z_value, z_conc, use = "pairwise.complete.obs"))


p.ss_cor_density <- toPlot %>% distinct(cor) %>%
  ggplot(aes(x = cor))+
    geom_density(fill = "black") +
    theme_bw(base_family = "Helvetica") +
    theme(aspect.ratio = 5/8) +
    theme(panel.grid = element_blank()) +
    xlab("Pearson correlation between steady-state metabolite concentration 
         predictions and measurements per kinase mutant")


toPlot <- toPlot %>% ungroup() %>% filter(dense_rank(-cor) <= 10)
  
stats = with(toPlot, cor.test(z_value, z_conc, use = "pairwise.complete.obs"))  

toPlot.stats <- data.frame(x = -2, 
                           y = 2, 
                           label = paste( "r = ", round(stats$estimate,2),"\n",   "p-value = ", format(stats$p.value, digits = 2, scientific = T), sep = ""))
p.ss_cor <- toPlot %>%
  filter(!is.na(z_conc), !is.na(z_value)) %>% 
  ggplot(aes(x = z_value, y = z_conc)) +
    geom_point() +
    geom_smooth(method = "lm") +
    geom_text(data = toPlot.stats, aes(x = x, y = y , label = label)) +
    theme_bw(base_family = "Helvetica") +
    theme(aspect.ratio = 5/8, legend.position = "none") +
    theme(panel.grid = element_blank()) +
    xlab("Predicted metabolite concentration based on ODE model, standartized value") +
    ylab("Experimetnaly measured metabolite 
         concentration, standartized value")

grid.arrange(p.ss_cor_density, p.ss_cor, ncol=2)
```

<embed src="Figure3_files/figure-markdown_github/steady_state_kinetics-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
``` r

ref = dataset.mca %>% filter(KO_name == "WT")

dataset.mca_FC <- left_join(dataset.mca, ref, by = c( "parameter", "var_type", "VAR")) %>%
  mutate(changeCC = abs(CC.x - CC.y)/CC.y) %>% 
  filter(is.nan(changeCC) != T , is.na(changeCC) != T, is.infinite(changeCC) != T) %>% 
  filter(KO_name.x != "WT")


toPlot <- dataset.mca_FC %>% 
  group_by(KO.x, KO_name.x, var_type) %>%
  dplyr::summarize(median_change_CC = median(changeCC, na.rm = T),
            n_changes = n_changes.x[1]/length(parameter_order)) %>%
  dplyr::mutate(jitter = rnorm(n= 1,  mean = 0,sd = 0.03)) %>%
  dplyr::mutate(cuts = cut(median_change_CC, breaks = c(0, 0.25, 0.5, 0.75, 1)),
                var_type = fct_relevel(var_type, c("flux", "conc")))
  
  
library(ggrepel)
toPlot %>% 
  ggplot(aes(x = cuts, y = n_changes + jitter)) +
    stat_boxplot(geom = "errorbar", width = 0.33) +
    geom_boxplot(width = 0.5, outlier.shape =  NA, fill="lightgrey") +
    geom_point() +
    geom_text_repel(aes(label = KO_name.x), size = 3, segment.alpha = 0.25) +
    facet_wrap(~var_type) +
    ylab("Fraction of changed enzymes in ODE model") +
    xlab("Median change of control coefficient across all\nparameters and reactions per mutant") +
    theme_bw() + theme(aspect.ratio = 5/8) + coord_flip()
```

<embed src="Figure3_files/figure-markdown_github/genes_in_parralel-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
