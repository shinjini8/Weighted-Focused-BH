# Weighted Focused BH
Implementation of multiple testing procedures introduced in "Leveraging the group structure of hypotheses for more powerful multiple testing with FDR control for the filtered rejection set"
## Dependencies

R is required for the simulation studies and data analyses in Sections 4 and 5 and Wed Appendices G and H.
Please install the following packages in R:
 - `igraph`
 - `multtest`
 - `StructSSI`
 
For data analysis, please install the following package in R:
 - `phyloseq`
 
For replicating the graphs in this paper, the following R packages are required:
 - `ggplot2`
 - `reshape2`
 - `ggpubr`
 - `stringr`

## External Code Dependencies

The simulation studies and data analyses compare the proposed methods with several existing multiple-testing procedures. Implementations of the STAR and Yekutieli procedures were obtained from the original repositories listed below. These source files are not redistributed in this repository.

### STAR

The implementation of the STAR procedure is from the code accompanying Lei et al. (2019), *A General Interactive Framework for False Discovery Rate Control under Structural Constraints*, available at:

https://github.com/lihualei71/STAR

The following files from the original STAR repository are required:

- `STAR_DAG.R`
- `summarize_methods.R`
- `expr_func.R`

These files should be obtained from the original STAR repository and placed in the appropriate working directory before running the corresponding simulation or data-analysis scripts.

### Yekutieli Procedure

The implementation of Yekutieli's hierarchical FDR procedure is from the code accompanying Katsevich, Sabatti, and Bogomolov (2023), *Filtering the Rejection Set while Preserving False Discovery Rate Control*, available in the Focused-BH repository:

https://github.com/ekatsevi/Focused-BH

In particular, the implementation used in our analyses is based on:

`src/methodology/Yekutieli.R`

The original file is available at:

https://github.com/ekatsevi/Focused-BH/blob/master/src/methodology/Yekutieli.R

This file should be obtained from the original Focused-BH repository before running the corresponding simulation or data-analysis scripts. The analysis code additionally uses project-specific utility functions for converting graph representations and summarizing the rejection set.
