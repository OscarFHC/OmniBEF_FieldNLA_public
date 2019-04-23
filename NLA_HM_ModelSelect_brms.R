######## Loading necessary libraries ########################################################################################
if (!require(tidyverse)) {
  install.packages("tidyverse", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(tidyverse)
}else{library(tidyverse)}

if (!require(vegan)) {
  install.packages("vegan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(vegan)
}else{library(vegan)}

if (!require(lavaan)) {
  install.packages("lavaan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(lavaan)
}else{library(lavaan)}

if (!require(nlme)) {
  install.packages("nlme", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(nlme)
}else{library(nlme)}

if (!require(piecewiseSEM)) {
  install.packages("piecewiseSEM", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(piecewiseSEM)
}else{library(piecewiseSEM)}

if (!require(brms)) {
  install.packages("brms", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(brms)
}else{library(brms)}

if (!require(blavaan)) {
  install.packages("blavaan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(blavaan)
}else{library(blavaan)}

if (!require(rstanarm)) {
  install.packages("rstanarm", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(rstanarm)
}else{library(rstanarm)}

if (!require(shinystan)) {
  install.packages("shinystan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(shinystan)
}else{library(shinystan)}

# if (!require(rstan)) {
#   install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#   library(rstan)
# }else{library(rstan)
#   rstan_options(auto_write = TRUE)
#   options(mc.cores = parallel::detectCores())}

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

if (!require(loo)) {
  install.packages("loo", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(loo)
}else{library(loo)}

if (!require(ggplot2)) {
  install.packages("ggplot2", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ggplot2)
}else{library(ggplot2)}

if (!require(viridis)) {
  install.packages("viridis", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(viridis)
}else{library(viridis)}

if (!require(GGally)) {
  install.packages("GGally", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(GGally)
}else{library(GGally)}

######## Loading necessary libraries ########################################################################################

##### Reading data from Github ##############################################################################################
bio07_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/NLA_Data/master/Final_cleaned/bio07_scaled.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE) 

bio12_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/NLA_Data/master/Final_cleaned/bio12_scaled.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
##### Reading data from Github ##############################################################################################

#############################################################################################################################
##### Finding the best model for year 2007 ##################################################################################
#############################################################################################################################
##### 2007 #########################################################################################################
### Step 1
fn0.1_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + 
                 (1 + ln_zpSRr*WOmni | ECO3))
fn0.2_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_phySRr + 
                 (1 + ln_zpSRr*WOmni + ln_phySRr | ECO3))
fn0.3_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_zpDen | ECO3))
fn0.4_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_NTL + 
                 (1 + ln_zpSRr*WOmni + ln_NTL | ECO3))
fn0.5_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL | ECO3))
fn0.6_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + 
                 (1 + ln_zpSRr*WOmni + Lat | ECO3))
fn0.7_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_DOC | ECO3))
fn0.8_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_Silica | ECO3))
fn0.9_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + pH + 
                 (1 + ln_zpSRr*WOmni + pH | ECO3))
fn0.10_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Temp + 
                  (1 + ln_zpSRr*WOmni + ln_Temp | ECO3))
fn0.11_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Sol + 
                  (1 + ln_zpSRr*WOmni + ln_Sol | ECO3))
fn0.12_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Cond + 
                  (1 + ln_zpSRr*WOmni + ln_Cond | ECO3))
Mod0.1_07 <- brm(formula = fn0.1_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.2_07 <- brm(formula = fn0.2_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.3_07 <- brm(formula = fn0.3_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.4_07 <- brm(formula = fn0.4_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.5_07 <- brm(formula = fn0.5_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.6_07 <- brm(formula = fn0.6_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.7_07 <- brm(formula = fn0.7_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.8_07 <- brm(formula = fn0.8_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.9_07 <- brm(formula = fn0.9_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.10_07 <- brm(formula = fn0.10_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.11_07 <- brm(formula = fn0.11_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.12_07 <- brm(formula = fn0.12_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo0.1_07 <- loo(Mod0.1_07, reloo = TRUE)
loo0.2_07 <- loo(Mod0.2_07, reloo = TRUE)
loo0.3_07 <- loo(Mod0.3_07, reloo = TRUE)
loo0.4_07 <- loo(Mod0.4_07, reloo = TRUE)
loo0.5_07 <- loo(Mod0.5_07, reloo = TRUE)
loo0.6_07 <- loo(Mod0.6_07, reloo = TRUE)
loo0.7_07 <- loo(Mod0.7_07, reloo = TRUE)
loo0.8_07 <- loo(Mod0.8_07, reloo = TRUE)
loo0.9_07 <- loo(Mod0.9_07, reloo = TRUE)
loo0.10_07 <- loo(Mod0.10_07, reloo = TRUE)
loo0.11_07 <- loo(Mod0.11_07, reloo = TRUE)
loo0.12_07 <- loo(Mod0.12_07, reloo = TRUE)

looic <- rbind(loo0.1_07$estimates["elpd_loo","Estimate"], loo0.2_07$estimates["elpd_loo","Estimate"],
               loo0.3_07$estimates["elpd_loo","Estimate"], loo0.4_07$estimates["elpd_loo","Estimate"],
               loo0.5_07$estimates["elpd_loo","Estimate"], loo0.6_07$estimates["elpd_loo","Estimate"],
               loo0.7_07$estimates["elpd_loo","Estimate"], loo0.8_07$estimates["elpd_loo","Estimate"],
               loo0.9_07$estimates["elpd_loo","Estimate"], loo0.10_07$estimates["elpd_loo","Estimate"],
               loo0.11_07$estimates["elpd_loo","Estimate"], loo0.12_07$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo0.1_07$estimates["elpd_loo","SE"], loo0.2_07$estimates["elpd_loo","SE"],
                  loo0.3_07$estimates["elpd_loo","SE"], loo0.4_07$estimates["elpd_loo","SE"],
                  loo0.5_07$estimates["elpd_loo","SE"], loo0.6_07$estimates["elpd_loo","SE"],
                  loo0.7_07$estimates["elpd_loo","SE"], loo0.8_07$estimates["elpd_loo","SE"],
                  loo0.9_07$estimates["elpd_loo","SE"], loo0.10_07$estimates["elpd_loo","SE"],
                  loo0.11_07$estimates["elpd_loo","SE"], loo0.12_07$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo0.1_07$pointwise[, "elpd_loo"], loo0.2_07$pointwise[, "elpd_loo"],
                                       loo0.3_07$pointwise[, "elpd_loo"], loo0.4_07$pointwise[, "elpd_loo"],
                                       loo0.5_07$pointwise[, "elpd_loo"], loo0.6_07$pointwise[, "elpd_loo"],
                                       loo0.7_07$pointwise[, "elpd_loo"], loo0.8_07$pointwise[, "elpd_loo"],
                                       loo0.9_07$pointwise[, "elpd_loo"], loo0.10_07$pointwise[, "elpd_loo"],
                                       loo0.11_07$pointwise[, "elpd_loo"], loo0.12_07$pointwise[, "elpd_loo"]))
Step1_07 <- round(cbind(looic, looic_SE, stacking_wts), 3)
write.table(Step1_07, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_07_step1.csv", sep = ",",
            row.names = TRUE, col.names = TRUE)
#summary(Mod0.5_07)
### Step 2
fn0.5_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL | ECO3))
fn1.1_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_phySR + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_phySR | ECO3))
fn1.2_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_zpDen | ECO3))
fn1.3_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL | ECO3))
fn1.4_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + Lat +
                 (1 + ln_zpSRr*WOmni + ln_PTL + Lat | ECO3))
fn1.5_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_DOC | ECO3))
fn1.6_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Silica | ECO3))
fn1.7_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + pH + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + pH | ECO3))
fn1.8_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Temp | ECO3))
fn1.9_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Sol + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Sol | ECO3))
fn1.10_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Cond + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Cond | ECO3))
Mod0.5_07 <- brm(formula = fn0.5_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.1_07 <- brm(formula = fn1.1_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.2_07 <- brm(formula = fn1.2_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.3_07 <- brm(formula = fn1.3_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.4_07 <- brm(formula = fn1.4_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.5_07 <- brm(formula = fn1.5_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.6_07 <- brm(formula = fn1.6_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.7_07 <- brm(formula = fn1.7_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.8_07 <- brm(formula = fn1.8_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.9_07 <- brm(formula = fn1.9_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.10_07 <- brm(formula = fn1.10_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo0.5_07 <- loo(Mod0.5_07, reloo = TRUE)
loo1.1_07 <- loo(Mod1.1_07, reloo = TRUE)
loo1.2_07 <- loo(Mod1.2_07, reloo = TRUE)
loo1.3_07 <- loo(Mod1.3_07, reloo = TRUE)
loo1.4_07 <- loo(Mod1.4_07, reloo = TRUE)
loo1.5_07 <- loo(Mod1.5_07, reloo = TRUE)
loo1.6_07 <- loo(Mod1.6_07, reloo = TRUE)
loo1.7_07 <- loo(Mod1.7_07, reloo = TRUE)
loo1.8_07 <- loo(Mod1.8_07, reloo = TRUE)
loo1.9_07 <- loo(Mod1.9_07, reloo = TRUE)
loo1.10_07 <- loo(Mod1.10_07, reloo = TRUE)

looic <- rbind(loo0.5_07$estimates["elpd_loo","Estimate"],
               loo1.1_07$estimates["elpd_loo","Estimate"], loo1.2_07$estimates["elpd_loo","Estimate"],
               loo1.3_07$estimates["elpd_loo","Estimate"], loo1.4_07$estimates["elpd_loo","Estimate"],
               loo1.5_07$estimates["elpd_loo","Estimate"], loo1.6_07$estimates["elpd_loo","Estimate"],
               loo1.7_07$estimates["elpd_loo","Estimate"], loo1.8_07$estimates["elpd_loo","Estimate"],
               loo1.9_07$estimates["elpd_loo","Estimate"], loo1.10_07$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo0.5_07$estimates["elpd_loo","SE"],
                  loo1.1_07$estimates["elpd_loo","SE"], loo1.2_07$estimates["elpd_loo","SE"],
                  loo1.3_07$estimates["elpd_loo","SE"], loo1.4_07$estimates["elpd_loo","SE"],
                  loo1.5_07$estimates["elpd_loo","SE"], loo1.6_07$estimates["elpd_loo","SE"],
                  loo1.7_07$estimates["elpd_loo","SE"], loo1.8_07$estimates["elpd_loo","SE"],
                  loo1.9_07$estimates["elpd_loo","SE"], loo1.10_07$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo0.5_07$pointwise[, "elpd_loo"],
                                       loo1.1_07$pointwise[, "elpd_loo"], loo1.2_07$pointwise[, "elpd_loo"],
                                       loo1.3_07$pointwise[, "elpd_loo"], loo1.4_07$pointwise[, "elpd_loo"],
                                       loo1.5_07$pointwise[, "elpd_loo"], loo1.6_07$pointwise[, "elpd_loo"],
                                       loo1.7_07$pointwise[, "elpd_loo"], loo1.8_07$pointwise[, "elpd_loo"],
                                       loo1.9_07$pointwise[, "elpd_loo"], loo1.10_07$pointwise[, "elpd_loo"]))
Step2_07 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step2_07, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_07_step2.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)

### Step 3
fn1.4_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL | ECO3))
fn2.1_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_phySRr +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_phySRr | ECO3))
fn2.2_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_zpDen +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_zpDen | ECO3))
fn2.3_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL | ECO3))
fn2.4_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_DOC +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_DOC | ECO3))
fn2.5_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_Silica +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_Silica | ECO3))
fn2.6_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + pH +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + pH | ECO3))
fn2.7_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_Temp +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_Temp | ECO3))
fn2.8_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_Sol +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_Sol | ECO3))
fn2.9_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_Cond +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_Cond | ECO3))
Mod1.4_07 <- brm(formula = fn1.4_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.1_07 <- brm(formula = fn2.1_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.2_07 <- brm(formula = fn2.2_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.3_07 <- brm(formula = fn2.3_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.4_07 <- brm(formula = fn2.4_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.5_07 <- brm(formula = fn2.5_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.6_07 <- brm(formula = fn2.6_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.7_07 <- brm(formula = fn2.7_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.8_07 <- brm(formula = fn2.8_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.9_07 <- brm(formula = fn2.9_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo1.4_07 <- loo(Mod1.4_07, reloo = TRUE)
loo2.1_07 <- loo(Mod2.1_07, reloo = TRUE)
loo2.2_07 <- loo(Mod2.2_07, reloo = TRUE)
loo2.3_07 <- loo(Mod2.3_07, reloo = TRUE)
loo2.4_07 <- loo(Mod2.4_07, reloo = TRUE)
loo2.5_07 <- loo(Mod2.5_07, reloo = TRUE)
loo2.6_07 <- loo(Mod2.6_07, reloo = TRUE)
loo2.7_07 <- loo(Mod2.7_07, reloo = TRUE)
loo2.8_07 <- loo(Mod2.8_07, reloo = TRUE)
loo2.9_07 <- loo(Mod2.9_07, reloo = TRUE)

looic <- rbind(loo1.4_07$estimates["elpd_loo","Estimate"],
               loo2.1_07$estimates["elpd_loo","Estimate"], loo2.2_07$estimates["elpd_loo","Estimate"],
               loo2.3_07$estimates["elpd_loo","Estimate"], loo2.4_07$estimates["elpd_loo","Estimate"],
               loo2.5_07$estimates["elpd_loo","Estimate"], loo2.6_07$estimates["elpd_loo","Estimate"],
               loo2.7_07$estimates["elpd_loo","Estimate"], loo2.8_07$estimates["elpd_loo","Estimate"],
               loo2.9_07$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo1.4_07$estimates["elpd_loo","SE"],
                  loo2.1_07$estimates["elpd_loo","SE"], loo2.2_07$estimates["elpd_loo","SE"],
                  loo2.3_07$estimates["elpd_loo","SE"], loo2.4_07$estimates["elpd_loo","SE"],
                  loo2.5_07$estimates["elpd_loo","SE"], loo2.6_07$estimates["elpd_loo","SE"],
                  loo2.7_07$estimates["elpd_loo","SE"], loo2.8_07$estimates["elpd_loo","SE"],
                  loo2.9_07$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo1.4_07$pointwise[, "elpd_loo"],
                                       loo2.1_07$pointwise[, "elpd_loo"], loo2.2_07$pointwise[, "elpd_loo"],
                                       loo2.3_07$pointwise[, "elpd_loo"], loo2.4_07$pointwise[, "elpd_loo"],
                                       loo2.5_07$pointwise[, "elpd_loo"], loo2.6_07$pointwise[, "elpd_loo"],
                                       loo2.7_07$pointwise[, "elpd_loo"], loo2.8_07$pointwise[, "elpd_loo"],
                                       loo2.9_07$pointwise[, "elpd_loo"]))
Step3_07 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step3_07, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_07_step3.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
#summary(Mod2.3_07)

### Step 4
fn2.3_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL +
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL | ECO3))
fn3.1_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr | ECO3))
fn3.2_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_zpDen | ECO3))
fn3.3_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_DOC | ECO3))
fn3.4_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_Silica | ECO3))
fn3.5_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + pH + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + pH | ECO3))
fn3.6_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_Temp | ECO3))
fn3.7_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_Sol + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_Sol | ECO3))
fn3.8_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_Cond + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_Cond | ECO3))
Mod2.3_07 <- brm(formula = fn2.3_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.1_07 <- brm(formula = fn3.1_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.2_07 <- brm(formula = fn3.2_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.3_07 <- brm(formula = fn3.3_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.4_07 <- brm(formula = fn3.4_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.5_07 <- brm(formula = fn3.5_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.6_07 <- brm(formula = fn3.6_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.7_07 <- brm(formula = fn3.7_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.8_07 <- brm(formula = fn3.8_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo2.3_07 <- loo(Mod2.3_07, reloo = TRUE)
loo3.1_07 <- loo(Mod3.1_07, reloo = TRUE)
loo3.2_07 <- loo(Mod3.2_07, reloo = TRUE)
loo3.3_07 <- loo(Mod3.3_07, reloo = TRUE)
loo3.4_07 <- loo(Mod3.4_07, reloo = TRUE)
loo3.5_07 <- loo(Mod3.5_07, reloo = TRUE)
loo3.6_07 <- loo(Mod3.6_07, reloo = TRUE)
loo3.7_07 <- loo(Mod3.7_07, reloo = TRUE)
loo3.8_07 <- loo(Mod3.8_07, reloo = TRUE)

looic <- rbind(loo2.3_07$estimates["elpd_loo","Estimate"],
               loo3.1_07$estimates["elpd_loo","Estimate"], loo3.2_07$estimates["elpd_loo","Estimate"],
               loo3.3_07$estimates["elpd_loo","Estimate"], loo3.4_07$estimates["elpd_loo","Estimate"],
               loo3.5_07$estimates["elpd_loo","Estimate"], loo3.6_07$estimates["elpd_loo","Estimate"],
               loo3.7_07$estimates["elpd_loo","Estimate"], loo3.8_07$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo2.3_07$estimates["elpd_loo","SE"],
                  loo3.1_07$estimates["elpd_loo","SE"], loo3.2_07$estimates["elpd_loo","SE"],
                  loo3.3_07$estimates["elpd_loo","SE"], loo3.4_07$estimates["elpd_loo","SE"],
                  loo3.5_07$estimates["elpd_loo","SE"], loo3.6_07$estimates["elpd_loo","SE"],
                  loo3.7_07$estimates["elpd_loo","SE"], loo3.8_07$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo2.3_07$pointwise[, "elpd_loo"],
                                       loo3.1_07$pointwise[, "elpd_loo"], loo3.2_07$pointwise[, "elpd_loo"],
                                       loo3.3_07$pointwise[, "elpd_loo"], loo3.4_07$pointwise[, "elpd_loo"],
                                       loo3.5_07$pointwise[, "elpd_loo"], loo3.6_07$pointwise[, "elpd_loo"],
                                       loo3.7_07$pointwise[, "elpd_loo"], loo3.8_07$pointwise[, "elpd_loo"]))
Step4_07 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step4_07, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_07_step4.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
#summary(Mod3.1_07)
### Step 5
fn3.1_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + 
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr | ECO3))
fn4.1_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_zpDen +  
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_zpDen | ECO3))
fn4.2_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_DOC +  
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_DOC | ECO3))
fn4.3_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_Silica +  
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_Silica | ECO3))
fn4.4_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + pH +  
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + pH | ECO3))
fn4.5_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_Temp +  
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_Temp | ECO3))
fn4.6_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_Sol +  
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_Sol | ECO3))
fn4.7_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_Cond +  
                 (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + ln_Cond | ECO3))
Mod3.1_07 <- brm(formula = fn3.1_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.1_07 <- brm(formula = fn4.1_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.2_07 <- brm(formula = fn4.2_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.3_07 <- brm(formula = fn4.3_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.4_07 <- brm(formula = fn4.4_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.5_07 <- brm(formula = fn4.5_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.6_07 <- brm(formula = fn4.6_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.7_07 <- brm(formula = fn4.7_07, data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo3.1_07 <- loo(Mod3.1_07, reloo = TRUE)
loo4.1_07 <- loo(Mod4.1_07, reloo = TRUE)
loo4.2_07 <- loo(Mod4.2_07, reloo = TRUE)
loo4.3_07 <- loo(Mod4.3_07, reloo = TRUE)
loo4.4_07 <- loo(Mod4.4_07, reloo = TRUE)
loo4.5_07 <- loo(Mod4.5_07, reloo = TRUE)
loo4.6_07 <- loo(Mod4.6_07, reloo = TRUE)
loo4.7_07 <- loo(Mod4.7_07, reloo = TRUE)

looic <- rbind(loo3.1_07$estimates["elpd_loo","Estimate"],
               loo4.1_07$estimates["elpd_loo","Estimate"], loo4.2_07$estimates["elpd_loo","Estimate"],
               loo4.3_07$estimates["elpd_loo","Estimate"], loo4.4_07$estimates["elpd_loo","Estimate"],
               loo4.5_07$estimates["elpd_loo","Estimate"], loo4.6_07$estimates["elpd_loo","Estimate"],
               loo4.7_07$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo3.1_07$estimates["elpd_loo","SE"],
                  loo4.1_07$estimates["elpd_loo","SE"], loo4.2_07$estimates["elpd_loo","SE"],
                  loo4.3_07$estimates["elpd_loo","SE"], loo4.4_07$estimates["elpd_loo","SE"],
                  loo4.5_07$estimates["elpd_loo","SE"], loo4.6_07$estimates["elpd_loo","SE"],
                  loo4.7_07$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo3.1_07$pointwise[, "elpd_loo"],
                                       loo4.1_07$pointwise[, "elpd_loo"], loo4.2_07$pointwise[, "elpd_loo"],
                                       loo4.3_07$pointwise[, "elpd_loo"], loo4.4_07$pointwise[, "elpd_loo"],
                                       loo4.5_07$pointwise[, "elpd_loo"], loo4.6_07$pointwise[, "elpd_loo"],
                                       loo4.7_07$pointwise[, "elpd_loo"]))
Step5_07 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step5_07, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_07_step5.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
summary(Mod4.4_07)



##### 2007 #########################################################################################################
#############################################################################################################################
##### Finding the best model for year 2007 ##################################################################################
#############################################################################################################################

#############################################################################################################################
##### Finding the best model for year 2012 ##################################################################################
#############################################################################################################################
##### 2012 #########################################################################################################
### Step 1
fn0.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + 
                 (1 + ln_zpSRr*WOmni | ECO3))
fn0.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_phySRr + 
                 (1 + ln_zpSRr*WOmni + ln_phySRr | ECO3))
fn0.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_zpDen | ECO3))
fn0.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_NTL + 
                 (1 + ln_zpSRr*WOmni + ln_NTL | ECO3))
fn0.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL | ECO3))
fn0.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + 
                 (1 + ln_zpSRr*WOmni + Lat | ECO3))
fn0.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_DOC | ECO3))
fn0.8_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_Silica | ECO3))
fn0.9_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + pH + 
                 (1 + ln_zpSRr*WOmni + pH | ECO3))
fn0.10_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Temp + 
                  (1 + ln_zpSRr*WOmni + ln_Temp | ECO3))
fn0.11_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Sol + 
                  (1 + ln_zpSRr*WOmni + ln_Sol | ECO3))
fn0.12_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_Cond + 
                  (1 + ln_zpSRr*WOmni + ln_Cond | ECO3))
Mod0.1_12 <- brm(formula = fn0.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.2_12 <- brm(formula = fn0.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.3_12 <- brm(formula = fn0.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.4_12 <- brm(formula = fn0.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.5_12 <- brm(formula = fn0.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.6_12 <- brm(formula = fn0.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.7_12 <- brm(formula = fn0.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.8_12 <- brm(formula = fn0.8_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.9_12 <- brm(formula = fn0.9_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.10_12 <- brm(formula = fn0.10_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.11_12 <- brm(formula = fn0.11_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0.12_12 <- brm(formula = fn0.12_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo0.1_12 <- loo(Mod0.1_12, reloo = TRUE)
loo0.2_12 <- loo(Mod0.2_12, reloo = TRUE)
loo0.3_12 <- loo(Mod0.3_12, reloo = TRUE)
loo0.4_12 <- loo(Mod0.4_12, reloo = TRUE)
loo0.5_12 <- loo(Mod0.5_12, reloo = TRUE)
loo0.6_12 <- loo(Mod0.6_12, reloo = TRUE)
loo0.7_12 <- loo(Mod0.7_12, reloo = TRUE)
loo0.8_12 <- loo(Mod0.8_12, reloo = TRUE)
loo0.9_12 <- loo(Mod0.9_12, reloo = TRUE)
loo0.10_12 <- loo(Mod0.10_12, reloo = TRUE)
loo0.11_12 <- loo(Mod0.11_12, reloo = TRUE)
loo0.12_12 <- loo(Mod0.12_12, reloo = TRUE)

looic <- rbind(loo0.1_12$estimates["elpd_loo","Estimate"], loo0.2_12$estimates["elpd_loo","Estimate"],
               loo0.3_12$estimates["elpd_loo","Estimate"], loo0.4_12$estimates["elpd_loo","Estimate"],
               loo0.5_12$estimates["elpd_loo","Estimate"], loo0.6_12$estimates["elpd_loo","Estimate"],
               loo0.7_12$estimates["elpd_loo","Estimate"], loo0.8_12$estimates["elpd_loo","Estimate"],
               loo0.9_12$estimates["elpd_loo","Estimate"], loo0.10_12$estimates["elpd_loo","Estimate"],
               loo0.11_12$estimates["elpd_loo","Estimate"], loo0.12_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo0.1_12$estimates["elpd_loo","SE"], loo0.2_12$estimates["elpd_loo","SE"],
                  loo0.3_12$estimates["elpd_loo","SE"], loo0.4_12$estimates["elpd_loo","SE"],
                  loo0.5_12$estimates["elpd_loo","SE"], loo0.6_12$estimates["elpd_loo","SE"],
                  loo0.7_12$estimates["elpd_loo","SE"], loo0.8_12$estimates["elpd_loo","SE"],
                  loo0.9_12$estimates["elpd_loo","SE"], loo0.10_12$estimates["elpd_loo","SE"],
                  loo0.11_12$estimates["elpd_loo","SE"], loo0.12_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo0.1_12$pointwise[, "elpd_loo"], loo0.2_12$pointwise[, "elpd_loo"],
                                       loo0.3_12$pointwise[, "elpd_loo"], loo0.4_12$pointwise[, "elpd_loo"],
                                       loo0.5_12$pointwise[, "elpd_loo"], loo0.6_12$pointwise[, "elpd_loo"],
                                       loo0.7_12$pointwise[, "elpd_loo"], loo0.8_12$pointwise[, "elpd_loo"],
                                       loo0.9_12$pointwise[, "elpd_loo"], loo0.10_12$pointwise[, "elpd_loo"],
                                       loo0.11_12$pointwise[, "elpd_loo"], loo0.12_12$pointwise[, "elpd_loo"]))
Step1_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step1_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step1.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)


### Step 2
fn0.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL | ECO3))
fn1.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_phySRr +
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_phySRr| ECO3))
fn1.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_zpDen | ECO3))
fn1.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL | ECO3))
fn1.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + Lat + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + Lat | ECO3))
fn1.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_DOC | ECO3))
fn1.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Silica | ECO3))
fn1.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + pH + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + pH | ECO3))
fn1.8_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Temp | ECO3))
fn1.9_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Sol +
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_Sol | ECO3))
fn1.10_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_Cond + 
                  (1 + ln_zpSRr*WOmni + ln_PTL + ln_Cond | ECO3))
Mod0.5_12 <- brm(formula = fn0.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.1_12 <- brm(formula = fn1.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.2_12 <- brm(formula = fn1.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.3_12 <- brm(formula = fn1.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.4_12 <- brm(formula = fn1.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.5_12 <- brm(formula = fn1.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.6_12 <- brm(formula = fn1.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.7_12 <- brm(formula = fn1.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.8_12 <- brm(formula = fn1.8_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.9_12 <- brm(formula = fn1.9_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1.10_12 <- brm(formula = fn1.10_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo0.5_12 <- loo(Mod0.5_12, reloo = TRUE)
loo1.1_12 <- loo(Mod1.1_12, reloo = TRUE)
loo1.2_12 <- loo(Mod1.2_12, reloo = TRUE)
loo1.3_12 <- loo(Mod1.3_12, reloo = TRUE)
loo1.4_12 <- loo(Mod1.4_12, reloo = TRUE)
loo1.5_12 <- loo(Mod1.5_12, reloo = TRUE)
loo1.6_12 <- loo(Mod1.6_12, reloo = TRUE)
loo1.7_12 <- loo(Mod1.7_12, reloo = TRUE)
loo1.8_12 <- loo(Mod1.8_12, reloo = TRUE)
loo1.9_12 <- loo(Mod1.9_12, reloo = TRUE)
loo1.10_12 <- loo(Mod1.10_12, reloo = TRUE)

looic <- rbind(loo0.5_12$estimates["elpd_loo","Estimate"],
               loo1.1_12$estimates["elpd_loo","Estimate"], loo1.2_12$estimates["elpd_loo","Estimate"],
               loo1.3_12$estimates["elpd_loo","Estimate"], loo1.4_12$estimates["elpd_loo","Estimate"],
               loo1.5_12$estimates["elpd_loo","Estimate"], loo1.6_12$estimates["elpd_loo","Estimate"],
               loo1.7_12$estimates["elpd_loo","Estimate"], loo1.8_12$estimates["elpd_loo","Estimate"],
               loo1.9_12$estimates["elpd_loo","Estimate"], loo1.10_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo0.5_12$estimates["elpd_loo","SE"],
                  loo1.1_12$estimates["elpd_loo","SE"], loo1.2_12$estimates["elpd_loo","SE"],
                  loo1.3_12$estimates["elpd_loo","SE"], loo1.4_12$estimates["elpd_loo","SE"],
                  loo1.5_12$estimates["elpd_loo","SE"], loo1.6_12$estimates["elpd_loo","SE"],
                  loo1.7_12$estimates["elpd_loo","SE"], loo1.8_12$estimates["elpd_loo","SE"],
                  loo1.9_12$estimates["elpd_loo","SE"], loo1.10_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo0.5_12$pointwise[, "elpd_loo"],
                                       loo1.1_12$pointwise[, "elpd_loo"], loo1.2_12$pointwise[, "elpd_loo"],
                                       loo1.3_12$pointwise[, "elpd_loo"], loo1.4_12$pointwise[, "elpd_loo"],
                                       loo1.5_12$pointwise[, "elpd_loo"], loo1.6_12$pointwise[, "elpd_loo"],
                                       loo1.7_12$pointwise[, "elpd_loo"], loo1.8_12$pointwise[, "elpd_loo"],
                                       loo1.9_12$pointwise[, "elpd_loo"], loo1.10_12$pointwise[, "elpd_loo"]))
Step2_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step2_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step2.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
#summary(Mod1.3_12)
### Step 3
fn1.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL | ECO3))
fn2.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr | ECO3))
fn2.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_zpDen | ECO3))
fn2.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + Lat + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + Lat | ECO3))
fn2.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_DOC | ECO3))
fn2.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Silica | ECO3))
fn2.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + pH + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + pH | ECO3))
fn2.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Temp | ECO3))
fn2.8_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Sol + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Sol | ECO3))
fn2.9_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Cond + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_Cond | ECO3))
Mod1.3_12 <- brm(formula = fn1.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.1_12 <- brm(formula = fn2.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.2_12 <- brm(formula = fn2.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.3_12 <- brm(formula = fn2.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.4_12 <- brm(formula = fn2.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.5_12 <- brm(formula = fn2.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.6_12 <- brm(formula = fn2.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.7_12 <- brm(formula = fn2.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.8_12 <- brm(formula = fn2.8_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod2.9_12 <- brm(formula = fn2.9_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo1.3_12 <- loo(Mod1.3_12, reloo = TRUE)
loo2.1_12 <- loo(Mod2.1_12, reloo = TRUE)
loo2.2_12 <- loo(Mod2.2_12, reloo = TRUE)
loo2.3_12 <- loo(Mod2.3_12, reloo = TRUE)
loo2.4_12 <- loo(Mod2.4_12, reloo = TRUE)
loo2.5_12 <- loo(Mod2.5_12, reloo = TRUE)
loo2.6_12 <- loo(Mod2.6_12, reloo = TRUE)
loo2.7_12 <- loo(Mod2.7_12, reloo = TRUE)
loo2.8_12 <- loo(Mod2.8_12, reloo = TRUE)
loo2.9_12 <- loo(Mod2.9_12, reloo = TRUE)

looic <- rbind(loo1.3_12$estimates["elpd_loo","Estimate"],
               loo2.1_12$estimates["elpd_loo","Estimate"], loo2.2_12$estimates["elpd_loo","Estimate"],
               loo2.3_12$estimates["elpd_loo","Estimate"], loo2.4_12$estimates["elpd_loo","Estimate"],
               loo2.5_12$estimates["elpd_loo","Estimate"], loo2.6_12$estimates["elpd_loo","Estimate"],
               loo2.7_12$estimates["elpd_loo","Estimate"], loo2.8_12$estimates["elpd_loo","Estimate"],
               loo2.9_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo1.3_12$estimates["elpd_loo","SE"],
                  loo2.1_12$estimates["elpd_loo","SE"], loo2.2_12$estimates["elpd_loo","SE"],
                  loo2.3_12$estimates["elpd_loo","SE"], loo2.4_12$estimates["elpd_loo","SE"],
                  loo2.5_12$estimates["elpd_loo","SE"], loo2.6_12$estimates["elpd_loo","SE"],
                  loo2.7_12$estimates["elpd_loo","SE"], loo2.8_12$estimates["elpd_loo","SE"],
                  loo2.9_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo1.3_12$pointwise[, "elpd_loo"],
                                       loo2.1_12$pointwise[, "elpd_loo"], loo2.2_12$pointwise[, "elpd_loo"],
                                       loo2.3_12$pointwise[, "elpd_loo"], loo2.4_12$pointwise[, "elpd_loo"],
                                       loo2.5_12$pointwise[, "elpd_loo"], loo2.6_12$pointwise[, "elpd_loo"],
                                       loo2.7_12$pointwise[, "elpd_loo"], loo2.8_12$pointwise[, "elpd_loo"],
                                       loo2.9_12$pointwise[, "elpd_loo"]))
Step3_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step3_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step3.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
#summary(Mod2.1_12)

### Step 4
fn2.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr | ECO3))
fn3.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_zpDen | ECO3))
fn3.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + Lat + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + Lat | ECO3))
fn3.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC | ECO3))
fn3.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Silica | ECO3))
fn3.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + pH + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + pH | ECO3))
fn3.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Temp | ECO3))
fn3.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Sol + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Sol | ECO3))
fn3.8_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Cond + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_Cond | ECO3))
Mod2.1_12 <- brm(formula = fn2.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.1_12 <- brm(formula = fn3.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.2_12 <- brm(formula = fn3.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.3_12 <- brm(formula = fn3.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.4_12 <- brm(formula = fn3.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.5_12 <- brm(formula = fn3.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.6_12 <- brm(formula = fn3.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.7_12 <- brm(formula = fn3.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod3.8_12 <- brm(formula = fn3.8_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo2.1_12 <- loo(Mod2.1_12, reloo = TRUE)
loo3.1_12 <- loo(Mod3.1_12, reloo = TRUE)
loo3.2_12 <- loo(Mod3.2_12, reloo = TRUE)
loo3.3_12 <- loo(Mod3.3_12, reloo = TRUE)
loo3.4_12 <- loo(Mod3.4_12, reloo = TRUE)
loo3.5_12 <- loo(Mod3.5_12, reloo = TRUE)
loo3.6_12 <- loo(Mod3.6_12, reloo = TRUE)
loo3.7_12 <- loo(Mod3.7_12, reloo = TRUE)
loo3.8_12 <- loo(Mod3.8_12, reloo = TRUE)

looic <- rbind(loo2.1_12$estimates["elpd_loo","Estimate"],
               loo3.1_12$estimates["elpd_loo","Estimate"], loo3.2_12$estimates["elpd_loo","Estimate"],
               loo3.3_12$estimates["elpd_loo","Estimate"], loo3.4_12$estimates["elpd_loo","Estimate"],
               loo3.5_12$estimates["elpd_loo","Estimate"], loo3.6_12$estimates["elpd_loo","Estimate"],
               loo3.7_12$estimates["elpd_loo","Estimate"], loo3.8_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo2.1_12$estimates["elpd_loo","SE"],
                  loo3.1_12$estimates["elpd_loo","SE"], loo3.2_12$estimates["elpd_loo","SE"],
                  loo3.3_12$estimates["elpd_loo","SE"], loo3.4_12$estimates["elpd_loo","SE"],
                  loo3.5_12$estimates["elpd_loo","SE"], loo3.6_12$estimates["elpd_loo","SE"],
                  loo3.7_12$estimates["elpd_loo","SE"], loo3.8_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo2.1_12$pointwise[, "elpd_loo"],
                                       loo3.1_12$pointwise[, "elpd_loo"], loo3.2_12$pointwise[, "elpd_loo"],
                                       loo3.3_12$pointwise[, "elpd_loo"], loo3.4_12$pointwise[, "elpd_loo"],
                                       loo3.5_12$pointwise[, "elpd_loo"], loo3.6_12$pointwise[, "elpd_loo"],
                                       loo3.7_12$pointwise[, "elpd_loo"], loo3.8_12$pointwise[, "elpd_loo"]))
Step4_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step4_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step4.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
summary(Mod3.3_12)
### Step 5
fn3.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC | ECO3))
fn4.1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_zpDen + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_zpDen | ECO3))
fn4.2_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + Lat + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + Lat | ECO3))
fn4.3_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Silica + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Silica | ECO3))
fn4.4_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + pH + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + pH | ECO3))
fn4.5_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Temp + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Temp | ECO3))
fn4.6_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Sol + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Sol | ECO3))
fn4.7_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Cond + 
                 (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + ln_Cond | ECO3))
Mod3.3_12 <- brm(formula = fn3.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.1_12 <- brm(formula = fn4.1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.2_12 <- brm(formula = fn4.2_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.3_12 <- brm(formula = fn4.3_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.4_12 <- brm(formula = fn4.4_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.5_12 <- brm(formula = fn4.5_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.6_12 <- brm(formula = fn4.6_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod4.7_12 <- brm(formula = fn4.7_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
loo3.3_12 <- loo(Mod3.3_12, reloo = TRUE)
loo4.1_12 <- loo(Mod4.1_12, reloo = TRUE)
loo4.2_12 <- loo(Mod4.2_12, reloo = TRUE)
loo4.3_12 <- loo(Mod4.3_12, reloo = TRUE)
loo4.4_12 <- loo(Mod4.4_12, reloo = TRUE)
loo4.5_12 <- loo(Mod4.5_12, reloo = TRUE)
loo4.6_12 <- loo(Mod4.6_12, reloo = TRUE)
loo4.7_12 <- loo(Mod4.7_12, reloo = TRUE)

looic <- rbind(loo3.3_12$estimates["elpd_loo","Estimate"],
               loo4.1_12$estimates["elpd_loo","Estimate"], loo4.2_12$estimates["elpd_loo","Estimate"],
               loo4.3_12$estimates["elpd_loo","Estimate"], loo4.4_12$estimates["elpd_loo","Estimate"],
               loo4.5_12$estimates["elpd_loo","Estimate"], loo4.6_12$estimates["elpd_loo","Estimate"],
               loo4.7_12$estimates["elpd_loo","Estimate"])
looic_SE <- rbind(loo3.3_12$estimates["elpd_loo","SE"],
                  loo4.1_12$estimates["elpd_loo","SE"], loo4.2_12$estimates["elpd_loo","SE"],
                  loo4.3_12$estimates["elpd_loo","SE"], loo4.4_12$estimates["elpd_loo","SE"],
                  loo4.5_12$estimates["elpd_loo","SE"], loo4.6_12$estimates["elpd_loo","SE"],
                  loo4.7_12$estimates["elpd_loo","SE"])
stacking_wts <- stacking_weights(cbind(loo3.3_12$pointwise[, "elpd_loo"],
                                       loo4.1_12$pointwise[, "elpd_loo"], loo4.2_12$pointwise[, "elpd_loo"],
                                       loo4.3_12$pointwise[, "elpd_loo"], loo4.4_12$pointwise[, "elpd_loo"],
                                       loo4.5_12$pointwise[, "elpd_loo"], loo4.6_12$pointwise[, "elpd_loo"],
                                       loo4.7_12$pointwise[, "elpd_loo"]))
Step5_12 <- round(cbind(looic, looic_SE, stacking_wts), 3)
# write.table(Step5_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/ModSelec_12_step5.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
#summary(Mod4.4_12)

##### 2012 #########################################################################################################
#############################################################################################################################
##### Finding the best model for year 2012 ##################################################################################
#############################################################################################################################
