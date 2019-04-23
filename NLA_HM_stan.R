######## Loading necessary libraries ########################################################################################
if (!require(tidyverse)) {
  install.packages("tidyverse", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(tidyverse)
}else{library(tidyverse)}

if (!require(nlme)) {
  install.packages("nlme", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(nlme)
}else{library(nlme)}

if (!require(lavaan)) {
  install.packages("lavaan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(lavaan)
}else{library(lavaan)}

if (!require(piecewiseSEM)) {
  install.packages("piecewiseSEM", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(piecewiseSEM)
}else{library(piecewiseSEM)}

if (!require(brms)) {
  install.packages("brms", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(brms)
}else{library(brms)}

if (!require(ggplot2)) {
  install.packages("ggplot2", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ggplot2)
}else{library(ggplot2)}

if (!require(vegan)) {
  install.packages("vegan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(vegan)
}else{library(vegan)}

if (!require(loo)) {
  install.packages("loo", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(loo)
}else{library(loo)}

if (!require(GGally)) {
  install.packages("GGally", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(GGally)
}else{library(GGally)}

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

if (!require(bayesplot)) {
  install.packages("bayesplot", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(bayesplot)
}else{library(bayesplot)}

# if (!require(rstan)) {
#   install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#   library(rstan)
# }else{library(rstan)
#   rstan_options(auto_write = TRUE)
#   options(mc.cores = parallel::detectCores())}

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

######## Loading necessary libraries ########################################################################################

##### Reading data from Github ##############################################################################################
bio07_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/NLA_Data/master/Final_cleaned/bio07_scaled.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE) 

bio12_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/NLA_Data/master/Final_cleaned/bio12_scaled.csv", 
                 header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
##### Reading data from Github ##############################################################################################

#############################################################################################################################
##### Hierarchical model ####################################################################################################
#############################################################################################################################
##### 2007 ##################################################################
bio_scale <- bio07_scale
bio_07to12 <- bio12_scale[which(bio12_scale$ECO3 %in% bio07_scale$ECO3 == TRUE),]

mm <- model.matrix(~ ln_phySRr + ln_NTL + ln_PTL + Lat, data = bio_scale)
zpSRmm <- model.matrix(~ Omni + Omni2, data = bio_scale)
newmm <- model.matrix(~ ln_phySRr + ln_NTL + ln_PTL + Lat, data = bio_07to12)
newzpSRmm <- model.matrix(~ Omni + Omni2, data = bio_07to12)

stanHM_Omni2No_Reg.dat <- # data list for Not region-denepdent but effects of zpSR depends on Omni analyses
  list(N = nrow(mm), # number of measurements
       K1 = ncol(mm), # number of column of the model matrix for determining phyDen
       X = mm, # number of column of the model matrix for determining phyDen
       J = length(unique(bio_scale[,"ECO3"])), # number of groups (i.e. number of eco-region, it is 48 here)
       region = as.numeric(as.factor(bio_scale[,"ECO3"])), # group (eco-region) id for each entry
       zpSR = bio_scale[,"ln_zpSRr"],
       K2 = ncol(zpSRmm), # number of column of the model matrix for determining the effects of zpSR
       X2 = zpSRmm, # number of column of the model matrix for determining the effects of zpSR
       phyDen = bio_scale[,"ln_phyDen"],
       
       newN = nrow(newmm), # number of measurements
       newK1 = ncol(newmm), # number of column of the model matrix for determining phyDen
       newX = newmm, # number of column of the model matrix for determining phyDen
       newJ = length(unique(bio_07to12[,"ECO3"])), # number of groups (i.e. number of eco-region, it is 48 here)
       newregion = as.numeric(as.factor(bio_07to12[,"ECO3"])), # group (eco-region) id for each entry
       newzpSR = bio_07to12[,"ln_zpSRr"],
       newK2 = ncol(newzpSRmm), # number of column of the model matrix for determining the effects of zpSR
       newX2 = newzpSRmm # number of column of the model matrix for determining the effects of zpSR
  )

stanHM_Omni2No_Reg07 = stan(file = "D:/Research/OmniBEF_NLA/NLA_HM/stan_scripts/HM/HM_OmniNoRegion_Reg_07to12.stan", 
                              data = stanHM_Omni2No_Reg.dat,
                              iter = 12000, warmup = 2000, chains = 4, 
                              control = list(adapt_delta = 0.9, max_treedepth = 15))
check_divergences(stanHM_Omni2No_Reg07)
save(stanHM_Omni2No_Reg07, 
     file = paste("E:/Monthlybackup/Research/OmniBEF_NLA/RData/stanHM_Omni2No_Reg_07to12.RData", sep = ""))
##### 2007 ##################################################################
##### manually summariz paramaters ##########################################
#load(file = "E:/Research_MonthBackup/OmniBEF_NLA/RData/Ver2/stanHM_Omni2No_Reg_07.RData")
param_Omni2No_Reg = extract(stanHM_Omni2No_Reg07)
beta.sum = as.data.frame(param_Omni2No_Reg[["gamma"]]) %>%
  gather(key = "gamma_n", value = "value") %>%
  group_by(gamma_n) %>%
  summarize(
    mean=mean(value),
    lo95=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma"]]))*0.025],
    lo90=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma"]]))*0.05],
    hi90=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma"]]))*0.95],
    hi95=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma"]]))*0.975]
  ) %>%
  mutate(gamma_n = c("Intercept", "ln_phySRr", "ln_NTL", "ln_PTL", "Lat"),
         beta = gamma_n) %>%
  select(-gamma_n)
write.table(beta.sum, file = "D:/Research/OmniBEF_NLA/RData/stanNLA_Omni2No_Reg07_beta.csv", sep = ",", col.names = TRUE, row.names = FALSE)

b.sum = as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]) %>%
  gather(key = "gamma_beta_n", value = "value" ) %>%
  group_by(gamma_beta_n) %>%
  summarize(
    mean=mean(value),
    lo95=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]))*0.025],
    lo90=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]))*0.05],
    hi90=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]))*0.95],
    hi95=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]))*0.975]
  ) %>%
  mutate(gamma_n = c("Intercept", "Omni", "Omni2"),
         beta = gamma_n) %>%
  select(-gamma_n)
write.table(b.sum, file = "D:/Research/OmniBEF_NLA/RData/stanNLA_Omni2No_Reg07_b.csv", sep = ",", col.names = TRUE, row.names = FALSE)

phy_pred = colMeans(as.data.frame(param_Omni2No_Reg[["y_pred"]])) # extracting the estimated means of each measurement
write.table(phy_pred, file = "D:/Research/OmniBEF_NLA/RData/stanNLA_Omni2No_Reg07_pred.csv", sep = ",", col.names = FALSE, row.names = FALSE)

phy_pred_07to12 = colMeans(as.data.frame(param_Omni2No_Reg[["y_12pred"]])) # extracting the "predicted" measurements
write.table(phy_pred_07to12, file = "D:/Research/OmniBEF_NLA/RData/stanNLA_Omni2No_Reg07to12_pred.csv", sep = ",", col.names = FALSE, row.names = FALSE)

rss = sum((phy_pred - bio07_scale[,"ln_phyDen"])^2)
totalrss = sum((mean(bio07_scale[,"ln_phyDen"]) - bio07_scale[,"ln_phyDen"])^2)
R2 = 1 - (rss/totalrss)
R2

bio_07to12 <- bio12_scale[which(bio12_scale$ECO3 %in% bio07_scale$ECO3 == TRUE),]
rss_07to12 = sum((phy_pred_07to12 - bio_07to12[,"ln_phyDen"])^2)
totalrss_07to12 = sum((mean(bio_07to12[,"ln_phyDen"]) - bio_07to12[,"ln_phyDen"])^2)
R2_07to12 = 1 - (rss_07to12/totalrss_07to12)
R2_07to12
##### manually summariz paramaters ##########################################

post <- extract(stanHM_Omni2No_Reg07)

beta0 <- as.data.frame(post[["beta"]][, , 1])
beta0.p <- beta0 %>%
  gather(key = "ECO3", value = "beta_intercept" ) %>%
  ggplot() + 
    geom_boxplot(aes(x = ECO3, y = beta_intercept)) + 
    scale_x_discrete(name = "",
                     breaks = names(beta0),
                     labels = unique(bio07_scale$ECO3)) + 
    scale_y_continuous(name = "")
p2 <- as.data.frame(post[["beta"]][, 21:40, 1]) %>%
  gather(key = "ECO3", value = "beta_intercept" ) %>%
  ggplot() + 
  geom_boxplot(aes(x = ECO3, y = beta_intercept)) + 
  # scale_x_discrete(name = "",
  #                  breaks = paste0("V", seq(21:41)),
  #                  labels = unique(bio07_scale$ECO3)[21:40]) + 
  scale_y_continuous(name = "")

#ppc_dens_overlay(bio07_scale$ln_phyDen, post$y_pred[runif(50, 0, nrow(post$y_pred)),])

##### 2012 ##################################################################
bio_scale <- bio12_scale
mm <- model.matrix(~ ln_phySRr + ln_NTL + ln_PTL + ln_DOC, data = bio_scale)
zpSRmm <- model.matrix(~ Omni + Omni2, data = bio_scale)

stanHM_Omni2No_Reg.dat <- # data list for Not region-denepdent but effects of zpSR depends on Omni analyses
  list(N = nrow(mm), # number of measurements
       K1 = ncol(mm), # number of column of the model matrix for determining phyDen
       X = mm, # number of column of the model matrix for determining phyDen
       J = length(unique(bio_scale[,"ECO3"])), # number of groups (i.e. number of eco-region, it is 48 here)
       region = as.numeric(as.factor(bio_scale[,"ECO3"])), # group (eco-region) id for each entry
       zpSR = bio_scale[,"ln_zpSRr"],
       K2 = ncol(zpSRmm), # number of column of the model matrix for determining the effects of zpSR
       X2 = zpSRmm, # number of column of the model matrix for determining the effects of zpSR
       phyDen = bio_scale[,"ln_phyDen"]
  )

stanHM_Omni2No_Reg12 = stan(file = "D:/Research/OmniBEF_NLA/NLA_HM/stan_scripts/HM/HM_OmniNoRegion_Reg.stan", 
                            data = stanHM_Omni2No_Reg.dat,
                            iter = 12000, warmup = 2000, chains = 4, 
                            control = list(adapt_delta = 0.9, max_treedepth = 15))
check_divergences(stanHM_Omni2No_Reg12)

save(stanHM_Omni2No_Reg12, 
     file = paste("E:/Monthlybackup/Research/OmniBEF_NLA/RData/stanHM_Omni2No_Reg_12_DOC.RData", sep = ""))
##### 2012 ##################################################################
##### manually summariz paramaters ##########################################
#load(file = "E:/Monthlybackup/Research/OmniBEF_NLA/RData/stanHM_Omni2No_Reg_07to12.RData")
param_Omni2No_Reg = extract(stanHM_Omni2No_Reg12)
beta.sum = as.data.frame(param_Omni2No_Reg[["gamma"]]) %>%
  gather(key = "gamma_n", value = "value") %>%
  group_by(gamma_n) %>%
  summarize(
    mean=mean(value),
    lo95=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma"]]))*0.025],
    lo90=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma"]]))*0.05],
    hi90=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma"]]))*0.95],
    hi95=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma"]]))*0.975]
  ) %>%
  mutate(gamma_n = c("Intercept", "ln_phySRr", "ln_NTL", "ln_PTL", "DOC"),
         beta = gamma_n) %>%
  select(-gamma_n)
write.table(beta.sum, file = "D:/Research/OmniBEF_NLA/RData/stanNLA_Omni2No_Reg12_DOC_beta.csv", sep = ",", col.names = TRUE, row.names = FALSE)

b.sum = as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]) %>%
  gather(key = "gamma_beta_n", value = "value" ) %>%
  group_by(gamma_beta_n) %>%
  summarize(
    mean=mean(value),
    lo95=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]))*0.025],
    lo90=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]))*0.05],
    hi90=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]))*0.95],
    hi95=sort(value)[nrow(as.data.frame(param_Omni2No_Reg[["gamma_beta_zpSR"]]))*0.975]
  ) %>%
  mutate(gamma_n = c("Intercept", "Omni", "Omni2"),
         beta = gamma_n) %>%
  select(-gamma_n)
write.table(b.sum, file = "D:/Research/OmniBEF_NLA/RData/stanNLA_Omni2No_Reg12_DOC_b.csv", sep = ",", col.names = TRUE, row.names = FALSE)

# colMeans(extract(stanHM_Omni2No_Reg12, par = "gamma_beta_zpSR")[["gamma_beta_zpSR"]])
phy_pred = colMeans(as.data.frame(param_Omni2No_Reg[["y_pred"]])) # extracting the estimated means of each measurement
write.table(phy_pred, file = "D:/Research/OmniBEF_NLA/RData/stanNLA_Omni2No_Reg12_DOC_pred.csv", sep = ",", col.names = FALSE, row.names = FALSE)

rss = sum((phy_pred - bio12_scale[,"ln_phyDen"])^2)
totalrss = sum((mean(bio12_scale[,"ln_phyDen"]) - bio12_scale[,"ln_phyDen"])^2)
R2 = 1 - (rss/totalrss)
R2
##### manually summariz paramaters ##########################################
#############################################################################################################################
##### Hierarchical model ####################################################################################################
#############################################################################################################################
