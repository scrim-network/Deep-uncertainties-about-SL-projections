# Code for *Sources and implications of deep uncertainties surrounding sea-level projections*

## Preamble
This code is intended to accompany the results of

> Bakker, AMR, Louchard, D, and Keller, K. Sources and implications of deep uncertainties surrounding sea-level projections, (conditionally accepted at Climatic Change).

This study explores reasons for different interpretations of estimating sea-level projections and their deep uncertainties. For example, the authors investigate the interpretation of the IPCC's likely range and the role of ice sheets in expert elicitation studies. 

Understanding the reasons for different interpretations *"is important for (i) the design of robust strategies and (ii) the exploration of pathways that may eventually lead to some kind of consensus distributions that are relatively straightforward to interpret."*

## Requirements
### Software
These scripts are written to run in R (tested under R v3.2.1; https://www.r-project.org/). They also require mulitple packages including:  
`reshape2`  
`lhs`  
`ggplot2`  
`grid`  
`gridExtra`

You can install and open the packages in R as shown in the example below, which installs the lhs package.

```R
install.packages("lhs")
library(lhs)
```
* Download the Rfiles
* Download the data. How to obtain the data is stated in the Appendix of the paper.
* Open the `.R` files in R or Rstudio. Change the paths in the `source` and `read.csv` to point to your directory with the appropriate files.
* Source the `main_script.R` file. `main_script.R` executes the scripts in the right order and produces the main figures in the paper.

*Figure may look slightly different than in the paper as a results of an R-update.*

## Contact
Alexander Bakker  
E-mail: alexander.bakker@rws.nl
