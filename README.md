# OmniBEF_FieldNLA_public
* Data in this repository:

(1) raw and scaled NLA dataset of year 2007 and year 2012 (bio07_or bio12_ ).  
(2) raw and scaled data from field experiment (FieldData_).  
(3) coordiante data from 2007 NLA dataset for selecting study sites (LakeforMap.csv).   

* Scripts in this repository:

(1) `NLA_HM_ModelSelect_brms.R`  
  model selection processes for NLA dataset analyses using [*brms*](https://github.com/paul-buerkner/brms) package).  
(2) `NLA_HM_stan.R` and `HM_OmniNoRegion_Reg.stan`  
  R and stan scripts for hierarchical model, where all effects are designed to be region-dependent (like random effect model) and effects of zooplankotn species richness (_zpSR_) is designed to quadratically depend on proportion of omnivores (_Omnip_).  
(3) `Figs_FieldNLA.R`  
  Final R script for making all plots and tables in the manuscript.
  
* All figures should be reproducible from the scripts and data in the repository.   
