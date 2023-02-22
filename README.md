---
title: "Introduction to ELVAR"
author:
- name: "Alok Maity and Andrew E. Teschendorff"
  affiliation: 
  - CAS Key Lab of Computational Biology, PICB, SINH
date: "2023-02-22"
package: ELVAR
output:
  BiocStyle::html_document:
    toc_float: true
---

# Summary

ELVAR is an R-package for differential abundance (DA) testing of cell-types in single-cell RNA-Seq data. It implements an Extended Louvain clustering Algorithm (EVA) that takes cell attribute information into acccount when inferring cellular communities from the cell-cell nearest neighbour graph. By taking cell attribute information into account, improved community detection is possible, which improves the power to detect DA-shifts of cell-types in relation to disease risk factors or disease itself.

# Installation

To install:

```r
library(devtools)
devtools::install_github("aet21/ELVAR")
```

# References

A generalized Louvain algorithm for improved inference of differential abundance from single-cell RNA-Seq data. Alok K Maity, Andrew E Teschendorff. Submitted 2022.
