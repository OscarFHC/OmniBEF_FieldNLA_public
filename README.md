# OmniBEF_FieldNLA

This repository contain the following data and code script for reproducing exactly the same figures and tables (maybe not including the figure and table style) in the maniscript.

* Data in this repository:

(1) raw and scaled NLA dataset of year 2007 and year 2012 (bio07_ or bio12_.csv ).  
(2) coordiante data from 2007 NLA dataset for selecting study sites (LakeforMap.csv).   
(3) raw and scaled data from field experiment (FieldData_.csv).  
(4) zooplankotn count data from the field experiment (Field_zpCount.csv)

* Scripts in this repository:

(1) [`NLA_HM_ModelSelect_brms.R`](https://github.com/OscarFHC/OmniBEF_FieldNLA_public/blob/master/NLA_HM_ModelSelect_brms.R)  
  model selection processes for NLA dataset analyses using [*brms*](https://github.com/paul-buerkner/brms) package).  
(2) [`Figs_FieldNLA.R`](https://github.com/OscarFHC/OmniBEF_FieldNLA_public/blob/master/Figs_FieldNLA.R)  
  Final R script for making all plots and tables in the manuscript.
(3)*NOT used but included here in case anyone is interested in*  
[`NLA_HM_stan.R`](https://github.com/OscarFHC/OmniBEF_FieldNLA_public/blob/master/NLA_HM_stan.R) and [`HM_OmniNoRegion_Reg.stan`](https://github.com/OscarFHC/OmniBEF_FieldNLA_public/blob/master/HM_OmniNoRegion_Reg.stan)  
  R and stan scripts for hierarchical model, where all effects are designed to be region-dependent (like random effect model) and effects of zooplankotn species richness (_zpSR_) is designed to quadratically depend on proportion of omnivores (_Omnip_).  
  
* All figures should be reproducible from the scripts and data in the repository.   
