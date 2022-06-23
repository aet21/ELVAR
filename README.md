---
title: "Introduction to ELVAR"
author:
- name: "Alok Maity and Andrew E. Teschendorff"
  affiliation: 
  - CAS Key Lab of Computational Biology, PICB, SINH
date: "2022-06-23"
package: ELVAR
output:
  BiocStyle::html_document:
    toc_float: true
---

# Summary

ELVAR is an R-package implementing an Extended Louvain clustering algorithm that takes cell attribute information into acccount when inferring cellular communities from scRNA-seq data. By taking cell attribute information into account, subtle cell clusters representing specific cell-states can be identified. As such, ELVAR can be used to infer differential abundance (DA) of cell-states in relation to disease risk factors or disease itself.

# Installation

To install:

```r
library(devtools)
devtools::install_github("aet21/ELVAR")
```

# References

A generalized Louvain algorithm for improved inference of differential abundance from single-cell RNA-Seq data. Alok K Maity, Andrew E Teschendorff. Submitted 2022.
