## Machine Learning Predicts the Yeast Metabolome from the Quantitative Proteome of Kinase Knockouts
### Zelezniak et al. Cell Systems. 2018.

Code repository for reproducing main analyses and figures of the study

## Viewing code and results
Main results presented in R notebooks with corresponding analyses are named FigureN.Rmd where N is figure number as per manuscript.
To view rendered figures directly on github, simply open on FigureN.md.

## To download code with data
Clone the repository from the command line:
```
$ git clone https://github.com/zelezniak-lab/kinase_metabolism.git
```
Figures R notebooks with corresponding analyses are named FigureN.Rmd where N is figure number as per manuscript.
Open in Rstudio and run the corresponding analyses.

## Dependencies

Most of the figure files plots should work solely with `tidyverse` package. 
There are however places where additional functions called from the following packages: `gridExtra, reshape2, RColorBrewer, ggthemes, cluster, ggrepel`.
Please see `sessionInfo()` output in the end each FigureN.md for more information.

If you have additional questions please contact me via email. 
You can find my email on our group's research page: http://www.zelezniaklab.org



