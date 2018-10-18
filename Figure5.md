Figure 5 results
================
Aleksej Zelezniak
2018-10-18

Machine Learning Regression Predicts the Concentration of Metabolite Pools from Enzyme Abundance

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
#iMM904 yeast metabolic model
load("./data/iMM904._load_.RData")
#experiment metadata
load("./data/exp_metadata._clean_.RData")

#caret models
load("./data/all_final_models.models_summary2.RData")
load("./data/file.list.models_summary2.RData")

#ID maps
metabolite.order <- read.delim("./data/metabolites.txt")
load("./data/metabolite2iMM904._load_.RData")
load("./data/gene.annotations._load_.RData")
load("./data/orf2name._clean_.RData")
```

### Figure 5C

``` r

metabolites.models.long <- all_final_models %>% 
  filter(isImputed == 0, metabolite != "Glu") %>%
  dplyr::select(model, RMSE, Rsquared, normalization, dataset, metabolite, degree, preprocessing) %>% 
  distinct() %>%
  group_by(metabolite, normalization, degree, preprocessing) %>%
  filter(RMSE == min(RMSE,na.rm = T)) %>%
  group_by(metabolite, degree) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))

length(unique(all_final_models$model))
#> [1] 12


metabolite.order = metabolite.order[with(metabolite.order,order(desc(method),pathway,Order, met_name)),]


toPlot <- metabolites.models.long %>% filter(metabolite %in% metabolite.order$metabolite)


toPlot <- toPlot %>% mutate(met_name = metabolite.order$met_name[match(metabolite, metabolite.order$metabolite)],
                            met_name = fct_relevel(met_name, levels=as.character(metabolite.order$met_name)),
                            pathway = metabolite.order$pathway[match(metabolite, metabolite.order$metabolite)],
                            pathway = fct_relevel(pathway,levels = unique(as.character(metabolite.order$pathway))))

tests.pvals <- toPlot %>% 
  group_by(degree) %>% 
  summarise(p.value1 = wilcox.test(toPlot$Rsquared[toPlot$degree == 1], Rsquared)$p.value,
            p.value3 = wilcox.test(toPlot$Rsquared[toPlot$degree == 3], Rsquared)$p.value) %>% 
  as.data.frame() %>% 
  gather(variable, value, -degree) %>% 
  filter(value != 1) %>%
  group_by(degree, variable) %>%
  mutate(Rsquared = jitter(0.7,0.2))


toPlot %>%
  ggplot(aes(x=degree, y = Rsquared, fill = degree)) +
    stat_boxplot(geom ='errorbar', width = 0.5) +
    geom_boxplot() +
    ylab(expression(paste("Cross-validated ", R^2, sep=""))) +
    xlab("Network distance") +
    geom_text(data=tests.pvals, aes(label=format(value, scientific=T))) +
    theme_bw() +
    theme(legend.position="none",
          panel.grid = element_blank())
```

<embed src="Figure5_files/figure-markdown_github/R2_barplot-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
### Figure 5B

``` r
toPlot.points <- toPlot %>% 
  group_by(metabolite) %>%
  filter(degree == 1) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))

toPlot <- toPlot %>% 
  group_by(metabolite) %>%
  filter(Rsquared == max(Rsquared,na.rm = T))
toPlot %>%
  ggplot(aes(x = met_name, y = Rsquared, fill=degree)) +
    geom_bar(position="dodge", stat = "identity") +
    geom_point(data = toPlot.points, aes(x = met_name, y = Rsquared)) +
    ylab(expression(paste("Cross-validated ", R^2, sep=""))) +
    xlab("") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)), limits = c(0,0.8)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.2, 0.7))
  
```

<embed src="Figure5_files/figure-markdown_github/distance-1.pdf" width="70%" style="display: block; margin: auto;" type="application/pdf" />
