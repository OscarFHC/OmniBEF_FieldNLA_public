######## Loading necessary libraries ########################################################################################
if (!require(tidyverse)) {
  install.packages("tidyverse", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(tidyverse)
}else{library(tidyverse)}

if (!require(vegan)) {
  install.packages("vegan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(vegan)
}else{library(vegan)}

if (!require(indicspecies)) {
  install.packages("indicspecies")
  library(indicspecies)
}else{library(indicspecies)}

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

if (!require(loo)) {
  install.packages("loo", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(loo)
}else{library(loo)}

if (!require(blavaan)) {
  install.packages("blavaan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(blavaan)
}else{library(blavaan)}

if (!require(shinystan)) {
  install.packages("shinystan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(shinystan)
}else{library(shinystan)}

if (!require(rstan)) {
  install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
  library(rstan)
}else{library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())}

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

if (!require(ggplot2)) {
  install.packages("ggplot2", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ggplot2)
}else{library(ggplot2)}

if (!require(viridis)) {
  install.packages("viridis", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(viridis)
}else{library(viridis)}

if (!require(ggridges)) {
  install.packages("ggridges", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ggridges)
}else{library(ggridges)}

if (!require(GGally)) {
  install.packages("GGally", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(GGally)
}else{library(GGally)}

if (!require(cowplot)) {
  install.packages("cowplot", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(cowplot)
}else{library(cowplot)}

if (!require(magick)) {
  install.packages("magick")
  library(magick)
}else{library(magick)}

if (!require(png)) {
  install.packages("png")
  library(png)
}else{library(png)}

if (!require(grid)) {
  install.packages("grid")
  library(grid)
}else{library(grid)}

if (!require(ggfortify)) {
  install.packages("ggfortify")
  library(ggfortify)
}else{library(ggfortify)}
######## Loading necessary libraries ########################################################################################

######### NLA analyses #########
##### Reading data from Github ######
bio07_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/bio07_scaled.csv", 
                       header = TRUE, stringsAsFactors = FALSE, fill = TRUE) 

bio12_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/bio12_scaled.csv", 
                       header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
##### Reading data from Github ######

##### 2007 parameterize model #####
fn0_07 <- bf(ln_phyDen ~ ln_phySRr + ln_NTL + ln_PTL + ln_zpSRr)
fn1_07 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + pH +  
               (1 + ln_zpSRr*WOmni + Lat + ln_PTL + ln_NTL + ln_phySRr + pH | ECO3))

Mod0_07 <- brm(formula = fn0_07,# + fncov1_07 + fncov2_07 + fncov3_07 + set_rescor(FALSE),
               data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0_07_loo <- loo(Mod0_07, reloo = TRUE)
Mod1_07 <- brm(formula = fn1_07, #+ fncov1_07 + fncov2_07 + fncov3_07 + set_rescor(FALSE), 
               data = bio07_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1_07_loo <- loo(Mod1_07, reloo = TRUE)


looic <- c(Mod0_07_loo$estimates["looic","Estimate"], Mod1_07_loo$estimates["looic","Estimate"])
elpd <- c(Mod0_07_loo$estimates["elpd_loo","Estimate"], Mod1_07_loo$estimates["elpd_loo","Estimate"])
elpd_SE <- c(Mod0_07_loo$estimates["elpd_loo","SE"], Mod1_07_loo$estimates["elpd_loo","SE"])
lpd_point <- cbind(Mod0_07_loo$pointwise[, "elpd_loo"], Mod1_07_loo$pointwise[, "elpd_loo"])
loo_list <- list(Mod0_07_loo, Mod1_07_loo)
waics <- c(waic(Mod0_07)$estimates[3, 1], waic(Mod1_07)$estimates[3, 1])
dwaic <- waics - min(waics)
waic_wts <- round(exp(-dwaic/2) / sum(exp(-dwaic/2)), 2)
pbma_wts <- pseudobma_weights(lpd_point, BB = FALSE)
pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)

##### 2007 parameterize model #####
##### 2007 R2 of the model #####
bio_07to12 <- bio12_scale[which(bio12_scale$ECO3 %in% bio07_scale$ECO3 == TRUE),]

mod0_pred07_07 <- predict(Mod0_07, newdata = bio07_scale, allow_new_levels = TRUE)
mod1_pred07_07 <- predict(Mod1_07, newdata = bio07_scale, allow_new_levels = TRUE)
mod0_pred07_12 <- predict(Mod0_07, newdata = bio_07to12, allow_new_levels = TRUE)
mod1_pred07_12 <- predict(Mod1_07, newdata = bio_07to12, allow_new_levels = TRUE)

totalrss07 <- sum((mean(bio07_scale[,"ln_phyDen"]) - bio07_scale[,"ln_phyDen"])^2)
totalrss07to12 <- sum((mean(bio_07to12[,"ln_phyDen"]) - bio_07to12[,"ln_phyDen"])^2)

mod0_07R2 <- 1 - ( sum((mod0_pred07_07[,"Estimate"] - bio07_scale[,"ln_phyDen"])^2) / totalrss07)
mod1_07R2 <- 1 - ( sum((mod1_pred07_07[,"Estimate"] - bio07_scale[,"ln_phyDen"])^2) / totalrss07)
mod0_07to12R2 <- 1 - ( sum((mod0_pred07_12[,"Estimate"] - bio_07to12[,"ln_phyDen"])^2) / totalrss07to12)
mod1_07to12R2 <- 1 - ( sum((mod1_pred07_12[,"Estimate"] - bio_07to12[,"ln_phyDen"])^2) / totalrss07to12)

R2 <- cbind(rbind(mod0_07R2, mod1_07R2), rbind(mod0_07to12R2, mod1_07to12R2))
colnames(R2) <- c("exp_R2", "pred_R2")

ic.weight_07 <- round(cbind(elpd, elpd_SE, stacking_wts, R2), 2)
summary(Mod0_07)
# write.table(ic.weight_07, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/IC_weight_07.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
##### 2007 R2 of the model #####

##### 2012 parameterize model #####
fn0_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr)
# fn1_12 <- bf(ln_phyDen ~ ln_NTL + ln_PTL + ln_phySRr + Lat + pH + ln_Silica + ln_zpSRr*WOmni +
#              (1 + ln_NTL + ln_PTL + ln_phySRr + Lat + pH + ln_Silica + ln_zpSRr*WOmni | ECO3))
fn1_12 <- bf(ln_phyDen ~ ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + pH + 
               (1 + ln_zpSRr*WOmni + ln_PTL + ln_NTL + ln_phySRr + ln_DOC + pH | ECO3))

Mod0_12 <- brm(formula = fn0_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod0_12_loo <- loo(Mod0_12, reloo = TRUE)
Mod1_12 <- brm(formula = fn1_12, data = bio12_scale, family = gaussian(), control = list(adapt_delta = 0.9))
Mod1_12_loo <- loo(Mod1_12, reloo = TRUE)


looic <- c(Mod0_12_loo$estimates["looic","Estimate"], Mod1_12_loo$estimates["looic","Estimate"])
elpd <- c(Mod0_12_loo$estimates["elpd_loo","Estimate"], Mod1_12_loo$estimates["elpd_loo","Estimate"])
elpd_SE <- c(Mod0_12_loo$estimates["elpd_loo","SE"], Mod1_12_loo$estimates["elpd_loo","SE"])
lpd_point <- cbind(Mod0_12_loo$pointwise[, "elpd_loo"], Mod1_12_loo$pointwise[, "elpd_loo"])
loo_list <- list(Mod0_12_loo, Mod1_12_loo)
waics <- c(waic(Mod0_12)$estimates[3, 1], waic(Mod1_12)$estimates[3, 1])
dwaic <- waics - min(waics)
waic_wts <- round(exp(-dwaic/2) / sum(exp(-dwaic/2)), 2)
pbma_wts <- pseudobma_weights(lpd_point, BB = FALSE)
pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)
##### 2012 parameterize model #####
##### 2012 R2 of the model #####
bio_12to07 <- bio07_scale[which(bio07_scale$ECO3 %in% bio12_scale$ECO3 == TRUE),]

mod0_pred12_12 <- predict(Mod0_12, newdata = bio12_scale, allow_new_levels = TRUE)
mod1_pred12_12 <- predict(Mod1_12, newdata = bio12_scale, allow_new_levels = TRUE)
mod0_pred12_07 <- predict(Mod0_12, newdata = bio_12to07, allow_new_levels = TRUE)
mod1_pred12_07 <- predict(Mod1_12, newdata = bio_12to07, allow_new_levels = TRUE)

totalrss12 <- sum((mean(bio12_scale[,"ln_phyDen"]) - bio12_scale[,"ln_phyDen"])^2)
totalrss12to07 <- sum((mean(bio_12to07[,"ln_phyDen"]) - bio_12to07[,"ln_phyDen"])^2)

mod0_12R2 <- 1 - ( sum((mod0_pred12_12[,"Estimate"] - bio12_scale[,"ln_phyDen"])^2) / totalrss12)
mod1_12R2 <- 1 - ( sum((mod1_pred12_12[,"Estimate"] - bio12_scale[,"ln_phyDen"])^2) / totalrss12)
mod0_12to07R2 <- 1 - ( sum((mod0_pred12_07[,"Estimate"] - bio_12to07[,"ln_phyDen"])^2) / totalrss12to07)
mod1_12to07R2 <- 1 - ( sum((mod1_pred12_07[,"Estimate"] - bio_12to07[,"ln_phyDen"])^2) / totalrss12to07)

R2 <- cbind(rbind(mod0_12R2, mod1_12R2), rbind(mod0_12to07R2, mod1_12to07R2))
colnames(R2) <- c("exp_R2", "pred_R2")

ic.weight_12 <- round(cbind(elpd, elpd_SE, stacking_wts, R2), 2)
summary(Mod1_12)
# write.table(ic.weight_12, file = "D:/Research/OmniBEF_NLA/NLA_HM/ModResults/brms_Omni/IC_weight_12.csv", sep = ",",
#             row.names = TRUE, col.names = TRUE)
##### 2012 R2 of the model #####
######### NLA analyses #########

######### Field analyses #########
##### Reading data from Github ######
Field.scale <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_scale.csv", 
                          sep = ",", header = TRUE) 
Field.raw <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_raw.csv", 
                        sep = ",", header = TRUE)
##### Reading data from Github ######

##### PCA on envi variables #####
envi_PC <- prcomp(Field.scale[,c("TN", "TP", "Cond", "pH", "Temp", "PAR")], center = TRUE, scale. = TRUE)
#summary(envi_PC)
enviPC1 <- envi_PC$x[,"PC1"]
enviPC2 <- envi_PC$x[,"PC2"]
enviPC3 <- envi_PC$x[,"PC3"]
Field.scale <- Field.scale %>%
  mutate(Envi1 = enviPC1, Envi2 = enviPC2, Envi3 = enviPC3)
##### PCA on envi variables #####

##### running models with GLM #####
glm_f  <- glm(Gzp ~ zpSR*Omnip + phyB + zpD + Envi1 + Envi2 + Envi3, data = Field.scale)
summary(glm_f)
# write.table(
#   as.data.frame(rbind(summary(glm_f)$coefficient[,"Estimate"],
#                       summary(glm_f)$coefficient[,"Estimate"] + summary(glm_f)$coefficient[,"Std. Error"]*qnorm(0.975),
#                       summary(glm_f)$coefficient[,"Estimate"] - summary(glm_f)$coefficient[,"Std. Error"]*qnorm(0.975))),
#   file = "D:/Manuscript/IGP_DivEffects_MS/MS_dissertation/MS_Field_NLA/Field_GLM_ModSum.csv", sep = ",", col.names = TRUE)
##### running models with GLM #####

##### running model with SEM #####
mod <- "
  Gzp ~ zpSR + Omnip + zpSR:Omnip + zpD + phyB
  zpSR ~ Envi1 + Envi2 + Envi3
  zpD ~ Envi1 + Envi2 + Envi3
  phyB ~ Envi1 + Envi2 + Envi3
  
  zpD ~~ phyB
"
mod_lavaan <- sem(mod, data = Field.scale, meanstructure = TRUE)
summary(mod_lavaan)
fitMeasures(mod_lavaan)
##### running model with SEM #####
######### Field analyses #########

######### Fig 3_zpSR and G VS Omni #########
##### preping zp composition data #####
Bio.raw <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_raw.csv", sep = ",", header = TRUE) %>%
  subset(select = c(Samp_ID, zpSR, phySR, zpD, phyB, All_Omni_prop, All_gzp, mzp_gzp, pH, Temp, Cond, PAR, TN_mean, TP_mean)) %>%
  mutate(Samp_ID = as.character(Samp_ID))
colnames(Bio.raw) <- c("Samp_ID", "zpSR", "phySR", "zpD", "phyB", "Omnip", "Gzp", "mzpGzp",
                       "pH", "Temp", "Cond", "PAR", "TN", "TP")

ZP_ID <- read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/Field_zpCount.csv", 
                  sep = ",", header = TRUE, fill = TRUE) %>%
  mutate(Samp_ID = as.character(Samp_ID)) %>%
  subset(Samp_ID %in% Bio.raw$Samp_ID)
# write.table(ZP_ID[which(duplicated(ZP_ID$TSN) == FALSE),] %>% subset(select = -c(SITE_ID, Samp_ID, count, Vol, Den)),
#             file = "D:/Research/OmniBEF_FieldExp/FieldExp_Data/SPcount/Field_zpSPList.csv", 
#             sep = ",", col.names = TRUE, row.names = FALSE)

zpSAD <- ZP_ID %>%
  group_by(FAMILY, FFG) %>%
  summarize(avgDen = mean(Den),
            sumDen = sum(Den))
zpSAD[which(zpSAD$FFG == "Herbivores"),"FG"] <- "H"
zpSAD[which(zpSAD$FFG != "Herbivores"),"FG"] <- "N"
zpSAD_ord <- zpSAD[,which(names(zpSAD) != "FFG")] %>%
  subset(FAMILY != "") %>%
  group_by(FAMILY, FG) %>%
  summarize(avgDen = sum(avgDen)) %>%
  arrange(desc(avgDen)) %>%
  ungroup() %>%
  mutate(FAMILY_ord = paste0(FAMILY, "_", FG)) %>%
  mutate(FAMILY_ord = factor(as.factor(FAMILY_ord), levels = FAMILY_ord)) 

ZPcom <- as.data.frame(matrix(0, length(unique(ZP_ID$Samp_ID)), length(unique(ZP_ID$TSN)) + 1))
colnames(ZPcom) <- c("Samp_ID", as.character(unique(ZP_ID$TSN)))

for (i in 1 : length(unique(ZP_ID$Samp_ID))){
  ZPcom[i, "Samp_ID"] <- unique(ZP_ID$Samp_ID)[i]
  temp <- ZP_ID[which(ZP_ID$Samp_ID == unique(ZP_ID$Samp_ID)[i]),]
  
  for (j in colnames(ZPcom)[which(colnames(ZPcom) %in% unique(temp$TSN))]){
    ZPcom[i, which(colnames(ZPcom) == j)] <- sum(temp[which(temp$TSN == j), "Den"]) * 1000
  }
}

ZPcom <- ZPcom %>%
  right_join(Bio.raw, by = "Samp_ID") %>%
  right_join(read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/Field_zpCount.csv", 
                      sep = ",", header = TRUE, fill = TRUE) %>%
               mutate(Samp_ID = as.character(Samp_ID)) %>%
               subset(Samp_ID %in% Bio.raw$Samp_ID) %>%
               group_by(Samp_ID) %>%
               summarize(
                 Omni_pt = sum(count[which(FFG != "Herbivores")]) / sum(count),
                 Omni_abspt = length(which(FFG != "Herbivores"))
               ), 
             by = "Samp_ID")
# write.table(ZPcom, file = "D:/Research/OmniBEF_FieldExp/FieldExp_Data/SPcount/Field_zpComm.csv",
#             sep = ",", col.names = TRUE, row.names = FALSE)
##### preping zp composition data #####
##### panel a Omnip ~ zpSR #####
A_zpSR_Omni <- ZPcom %>%
  ggplot() + 
  geom_point(aes(x = Omnip*100, y = zpSR), size = 4) + 
  geom_smooth(aes(x = Omnip*100, y = zpSR), se = FALSE, method = "lm", color = "black", linetype = "dashed", size = 3) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-2.5, 75)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(5, 26)) + 
  labs(x = "", #expression(atop("Omnivorous consumption", "(% microzooplankton density consumed)")),
       y = expression("Zooplankton taxonomic richness (" * italic(zpSR) * ")")) +
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 42, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        axis.title.y = element_text(size = 42, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        plot.margin = margin(t = 24, r = 36, b = 0, l = 12, "pt"))
##### panel a Omnip ~ zpSR #####
##### panel b Omnip ~ G #####
B_Gzp_Omni <- ZPcom %>%
  ggplot() + 
  geom_point(aes(x = Omnip*100, y = Gzp), size = 4) + 
  geom_smooth(aes(x = Omnip*100, y = Gzp), se = FALSE, method = "lm", color = "black", linetype = "dashed", size = 3) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-2.5, 75)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
  labs(x = "",
       y = expression(atop("% Phytoplankton biomass consumed", 
                           "by entire zooplankton community ( " * italic("G") * ")"))) +
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 42, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        axis.title.y = element_text(size = 42, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        plot.margin = margin(t = 24, r = 36, b = 0, l = 12, "pt"))
##### panel b Omnip ~ G #####
Omni.p <- plot_grid(A_zpSR_Omni, B_Gzp_Omni, ncol = 2, labels = c("a.", "b."), label_size = 30) + 
  draw_label(label = expression("Omnivorous consumption (" * italic(Omni) * "; % microzooplankton density consumed)"),
             x = 0.5, y = 0.02, hjust = 0.5, vjust = 0.5, colour = "black", size = 42)
# ggsave(plot = Omni.p,
#        file = "D:/Manuscript/IGP_DivEffects_MS/MS_dissertation/MS_Field_NLA/Figs/Fig3_Omni.tiff",
#        width = 60, height = 32, units = c("cm"),
#        dpi = 600)
######### Fig 3_zpSR and G VS Omni #########

######### Fig 4_ ZP taxa (Cyclops, Daphnia) VS Omni #########
##### preping zp composition data #####
Bio.raw <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_raw.csv", sep = ",", header = TRUE) %>%
  subset(select = c(Samp_ID, zpSR, phySR, zpD, phyB, All_Omni_prop, All_gzp, mzp_gzp, pH, Temp, Cond, PAR, TN_mean, TP_mean)) %>%
  mutate(Samp_ID = as.character(Samp_ID))
colnames(Bio.raw) <- c("Samp_ID", "zpSR", "phySR", "zpD", "phyB", "Omnip", "Gzp", "mzpGzp",
                       "pH", "Temp", "Cond", "PAR", "TN", "TP")

ZP_ID <- read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/Field_zpCount.csv", 
                  sep = ",", header = TRUE, fill = TRUE) %>%
  mutate(Samp_ID = as.character(Samp_ID)) %>%
  subset(Samp_ID %in% Bio.raw$Samp_ID)
# write.table(ZP_ID[which(duplicated(ZP_ID$TSN) == FALSE),] %>% subset(select = -c(SITE_ID, Samp_ID, count, Vol, Den)),
#             file = "D:/Research/OmniBEF_FieldExp/FieldExp_Data/SPcount/Field_zpSPList.csv", 
#             sep = ",", col.names = TRUE, row.names = FALSE)

zpSAD <- ZP_ID %>%
  group_by(FAMILY, FFG) %>%
  summarize(avgDen = mean(Den),
            sumDen = sum(Den))
zpSAD[which(zpSAD$FFG == "Herbivores"),"FG"] <- "H"
zpSAD[which(zpSAD$FFG != "Herbivores"),"FG"] <- "N"
zpSAD_ord <- zpSAD[,which(names(zpSAD) != "FFG")] %>%
  subset(FAMILY != "") %>%
  group_by(FAMILY, FG) %>%
  summarize(avgDen = sum(avgDen)) %>%
  arrange(desc(avgDen)) %>%
  ungroup() %>%
  mutate(FAMILY_ord = paste0(FAMILY, "_", FG)) %>%
  mutate(FAMILY_ord = factor(as.factor(FAMILY_ord), levels = FAMILY_ord)) 

ZPcom <- as.data.frame(matrix(0, length(unique(ZP_ID$Samp_ID)), length(unique(ZP_ID$TSN)) + 1))
colnames(ZPcom) <- c("Samp_ID", as.character(unique(ZP_ID$TSN)))

for (i in 1 : length(unique(ZP_ID$Samp_ID))){
  ZPcom[i, "Samp_ID"] <- unique(ZP_ID$Samp_ID)[i]
  temp <- ZP_ID[which(ZP_ID$Samp_ID == unique(ZP_ID$Samp_ID)[i]),]
  
  for (j in colnames(ZPcom)[which(colnames(ZPcom) %in% unique(temp$TSN))]){
    ZPcom[i, which(colnames(ZPcom) == j)] <- sum(temp[which(temp$TSN == j), "Den"]) * 1000
  }
}

ZPcom <- ZPcom %>%
  right_join(Bio.raw, by = "Samp_ID") %>%
  right_join(read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/Field_zpCount.csv", 
                      sep = ",", header = TRUE, fill = TRUE) %>%
               mutate(Samp_ID = as.character(Samp_ID)) %>%
               subset(Samp_ID %in% Bio.raw$Samp_ID) %>%
               group_by(Samp_ID) %>%
               summarize(
                 Omni_pt = sum(count[which(FFG != "Herbivores")]) / sum(count),
                 Omni_abspt = length(which(FFG != "Herbivores"))
               ), 
             by = "Samp_ID")
# write.table(ZPcom, file = "D:/Research/OmniBEF_FieldExp/FieldExp_Data/SPcount/Field_zpComm.csv",
#             sep = ",", col.names = TRUE, row.names = FALSE)
##### preping zp composition data #####
##### zp abund rank #####
zp_rank.p <- zpSAD_ord %>%
  ggplot() + 
  geom_bar(aes(x = FAMILY_ord, y = avgDen), stat = "identity", width = 0.8,  position = "dodge2") + 
  scale_fill_manual(name = "") + 
  labs(x = "Zooplankton family",
       y = expression("Average density (ind./mL)")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, size = 24),
        axis.text.y = element_text(size = 24),
        axis.title = element_text(size = 32, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        #axis.title.y.right = element_text(size = 32, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        plot.margin = margin(t = 24, r = 36, b = 24, l = 36, "pt"),
        legend.position = c(0, 0.95), 
        legend.title = element_text(size = 24), 
        legend.text = element_text(size = 24))
# ggsave(zp_rank.p,
#        file = "D:/Manuscript/IGP_DivEffects_MS/MS_dissertation/MS_Field_NLA/Figs/SupFigs/FigS3_zp_AbundRank.tiff",
#        width = 58, height = 32, units = c("cm"),
#        dpi = 600)
zpSAD_ord <- zpSAD_ord %>% mutate(prop = avgDen/sum(avgDen))
sum(zpSAD_ord$prop[1:10])
##### zp abund rank #####
##### panel a, b Cyclops, Daphnia ~ Omni #####
for(i in 1:10){
  domTSN <- as.character(unique(ZP_ID[which(ZP_ID$FAMILY %in% as.character(zpSAD_ord$FAMILY[i])),"TSN"]))
  if (length(domTSN)>1){
    ZPcom[, paste0("domzp", i)] <- rowSums(ZPcom[, domTSN])/(ZPcom[, "zpD"]*1000)
  }else{ZPcom[, paste0("domzp", i)] <- ZPcom[, domTSN]/(ZPcom[, "zpD"]*1000)}
}
names(ZPcom)

summary(lm(Omnip ~ domzp1 + domzp2 + domzp3 + domzp4 + domzp5 + domzp6 + domzp7 + domzp8 + domzp9 + domzp10, data = ZPcom))
max(ZPcom$domzp3)
zpSAD_ord$FAMILY[3]
zpSAD_ord$FAMILY[8]

A_Cyclops_Omni <- ZPcom %>%
  ggplot() + 
  geom_point(aes(x = domzp3*100, y = Omnip*100), size = 4) +
  geom_smooth(aes(x = domzp3*100, y = Omnip*100), se = FALSE, method = "lm", color = "black", size = 3, linetype = "solid") + 
  scale_x_continuous(expand = c(0, 0), limits = c(-2.5, 100*max(ZPcom$domzp3)+2.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) +
  labs(x = expression(atop("% of " * italic(Cyclopidae), "") ),
       y = expression(atop("Omnivorous consumption (" * italic(Omni) * ";", "% microzooplankton density consumed)"))) +
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 42, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        axis.title.y = element_text(size = 42, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        plot.margin = margin(t = 42, r = 12, b = 24, l = 12, "pt"))
B_Daphnia_Omni <- ZPcom %>%
  ggplot() + 
  geom_point(aes(x = domzp8*100, y = Omnip*100), size = 4) +
  geom_smooth(aes(x = domzp8*100, y = Omnip*100), se = FALSE, method = "lm", color = "black", size = 3, linetype = "solid") + 
  scale_x_continuous(expand = c(0, 0), limits = c(-2.5, 100*max(ZPcom$domzp8)+2.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) +
  labs(x = expression(atop("% of " * italic(Daphniidae), "") ),
       y = "") +  #expression(atop("Omnivorous consumption (" * italic(Omni) * ";", "% microzooplankton density consumed)"))) +
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 42, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        axis.title.y = element_text(size = 42, margin = margin(t = 24, r = 24, b = 24, l = 24, "pt")),
        plot.margin = margin(t = 42, r = 12, b = 24, l = 12, "pt"))
Omni_zp <- plot_grid(A_Cyclops_Omni, B_Daphnia_Omni, ncol = 2, labels = c("a.", "b."), label_size = 30)
# ggsave(plot = Omni_zp,
#        file = "D:/Manuscript/IGP_DivEffects_MS/MS_dissertation/MS_Field_NLA/Figs/Fig4_Omni_zp.jpeg",
#        width = 60, height = 32, units = c("cm"),
#        dpi = 600)
##### panel a, b Cyclops, Daphnia ~ Omni #####
######### Fig 4_ ZP taxa (Cyclops, Daphnia) VS Omni #########


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

if (!require(cowplot)) {
  install.packages("cowplot", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(cowplot)
}else{library(cowplot)}

if (!require(ggfortify)) {
  install.packages("ggfortify")
  library(ggfortify)
}else{library(ggfortify)}

######## Loading necessary libraries ########################################################################################

######### Fig S1_NLA correlation plots ##########
bio07_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/bio07_scaled.csv", 
                       header = TRUE, stringsAsFactors = FALSE, fill = TRUE) 
bio12_scale = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/bio12_scaled.csv", 
                       header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

bio07_raw = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/bio07_unscaled_raw.csv", 
                     header = TRUE, stringsAsFactors = FALSE, fill = TRUE) %>%
  mutate(phyB = ln_phyDen, phySR = ln_phySRr, zpDen = ln_zpDen, zpSR = ln_zpSRr, Omni = WOmni,
         Lat = LAT_DD, TN = NTL, TP = PTL, Cond = ln_Cond, DOC = ln_DOC, Silica = ln_Silica, Temp = ln_Temp, Light = ln_Sol)
bio12_raw = read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/bio12_unscaled_raw.csv", 
                     header = TRUE, stringsAsFactors = FALSE, fill = TRUE) %>%
  mutate(phyB = ln_phyDen, phySR = ln_phySRr, zpDen = ln_zpDen, zpSR = ln_zpSRr, Omni = WOmni,
         Lat = LAT_DD, TN = NTL, TP = PTL, Cond = ln_Cond, DOC = ln_DOC, Silica = ln_Silica, Temp = ln_Temp, Light = ln_Sol)

NLA07_Corr <- plot_grid(
  ggmatrix_gtable(ggpairs(bio07_scale,
                          column = c("ln_phyDen", "ln_phySR", "ln_zpDen", "ln_zpSR", "WOmni"),
                          columnLabels = c("Log[phyB]", "Log[zpSR]", "log[phySR]", "Log[zpDen]", "Omnip"),
                          upper = list(continuous = wrap("cor", size = 4)),
                          lower = list(continuous = "smooth"))),
  ggmatrix_gtable(ggpairs(bio07_scale,
                          column = c("ln_NTL", "ln_PTL", "Lat", "ln_DOC", "ln_Silica", "pH", "ln_Temp", "ln_Sol", "ln_Cond"),
                          columnLabels = c("Log(TN)", "Log(TP)", "Lat", "Log(DOC)", "Log(Silica)", "pH", "Log(Temp)", 
                                           "Log(Sol)", "Log(Cond)"),
                          upper = list(continuous = wrap("cor", size = 4)),
                          lower = list(continuous = "smooth"))),
  ncol = 1, labels = c("a.", "b."), label_size = 30, scale = 0.95)
# ggsave(NLA07_Corr,
#        file = "D:/Manuscript/IGP_DivEffects_MS_Figs/Field_NLA_Figs/SupFigs/FieldNLA_fs1_1_NLA07_Corr.pdf", 
#        width = 40, height = 40, units = c("cm"), dpi = 600)

NLA12_Corr <- plot_grid(
  ggmatrix_gtable(ggpairs(bio12_scale,
                          column = c("ln_phyDen", "ln_phySR", "ln_zpDen", "ln_zpSR", "WOmni"),
                          columnLabels = c("Log[phyB]", "Log[zpSR]", "log[phySR]", "Log[zpDen]", "Omnip"),
                          upper = list(continuous = wrap("cor", size = 4)),
                          lower = list(continuous = "smooth"))),
  ggmatrix_gtable(ggpairs(bio12_scale,
                          column = c("ln_NTL", "ln_PTL", "Lat", "ln_DOC", "ln_Silica", "pH", "ln_Temp", "ln_Sol", "ln_Cond"),
                          columnLabels = c("Log(TN)", "Log(TP)", "Lat", "Log(DOC)", "Log(Silica)", "pH", "Log(Temp)", 
                                           "Log(Sol)", "Log(Cond)"),
                          upper = list(continuous = wrap("cor", size = 4)),
                          lower = list(continuous = "smooth"))),
  ncol = 1, labels = c("a.", "b."), label_size = 30, scale = 0.95)
# ggsave(NLA12_Corr,
#        file = "D:/Manuscript/IGP_DivEffects_MS_Figs/Field_NLA_Figs/SupFigs/FieldNLA_fs1_2_NLA12_Corr.pdf", 
#        width = 40, height = 40, units = c("cm"), dpi = 600)
######### Fig S1_NLA correlation plots ##########



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

if (!require(cowplot)) {
  install.packages("cowplot", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(cowplot)
}else{library(cowplot)}

if (!require(ggfortify)) {
  install.packages("ggfortify")
  library(ggfortify)
}else{library(ggfortify)}

if (!require(png)) {
  install.packages("png")
  library(png)
}else{library(png)}

if (!require(grid)) {
  install.packages("grid")
  library(grid)
}else{library(grid)}

######## Loading necessary libraries ########################################################################################

######### Site maps ##########
Lake <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/LakeforMap.csv", 
                   sep = ",", header = TRUE)
img <- readPNG("D:/Research/OmniBEF_FieldExp/flowCAM pics/home.png")
home <- rasterGrob(img, interpolate = TRUE) 

map.MI <- map_data("state", c("michigan:south"), boundary = TRUE, interior = FALSE)
map.OH <- map_data("state", c("ohio"), boundary = TRUE, interior = FALSE)
map.IN <- map_data("state", c("indiana"), boundary = TRUE, interior = FALSE)

SiteMap <- ggplot() + 
  geom_polygon(data = map.MI, aes(x = long, y = lat), fill = NA, color = "black") + 
  geom_polygon(data = map.OH, aes(x = long, y = lat), fill = NA, color = "black") + 
  geom_polygon(data = map.IN, aes(x = long, y = lat), fill = NA, color = "black") + 
  coord_cartesian(xlim = c( (min(Lake$LON_DD) - 0.7), (max(Lake$LON_DD) + 0.7) ),
                  ylim = c( (min(Lake$LAT_DD) - 0.5), (max(Lake$LAT_DD) + 0.5) )) + 
  geom_point(data = Lake, aes(x = LON_DD, y = LAT_DD), size = 5) + 
  annotation_custom(home, xmin = -83.49, xmax = -83.99, ymin = 42.19, ymax = 42.39) + 
  geom_point(aes(x = -83.7430, y = 42.2808), color = "red", size = 7, shape = 64) +
  labs(x = "Longitude", y = "Latitude") + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, "pt"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 36, margin = margin(t = 6, r = 0, b = 6, l = 0, "pt"), hjust = 0.5),
        axis.title.y = element_text(size = 36, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")),
        axis.text = element_text(size = 42))
# ggsave(plot = SiteMap, width = 5.37*5, height = 5*5, units = "cm",
#        file = "D:/Manuscript/IGP_DivEffects_MS_Figs/Field_NLA_Figs/SiteMap.png", dpi = 1200)

map.US <- map_data("usa", boundary = TRUE, interior = FALSE)
USMap <- ggplot() + 
  geom_polygon(data = map.US, aes(x = long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3) +
  labs(x = "Longitude", y = "Latitude") + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, "pt"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 24, margin = margin(t = 6, r = 0, b = 6, l = 0, "pt"), hjust = 0.5),
        axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")),
        axis.text = element_text(size = 28))
# ggsave(plot = USMap, width = 5.37*5, height = 5*5, units = "cm",
#        file = "D:/Manuscript/IGP_DivEffects_MS_Figs/Field_NLA_Figs/USMap.png", dpi = 1200)

# ggdraw() +
#   draw_plot(USMap + theme(legend.justification = "bottom"), 0, 0, 1, 1) +
#   draw_plot(SiteMap +
#               theme(legend.justification = "top"), 0.5, 0.52, 0.5, 0.4)
######### Site maps ##########

######### Fig S2_Field Bio correlation plots ##########
Dat.scale <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_scale.csv", 
                        sep = ",", header = TRUE) 
Bio.raw <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_raw.csv", 
                      sep = ",", header = TRUE) %>%
  subset(select = c(Samp_ID, zpSR, phySR, zpD, phyB, All_Omni_prop, All_gzp, mzp_gzp, pH, Temp, Cond, PAR, TN_mean, TP_mean)) %>%
  mutate(Samp_ID = as.character(Samp_ID))
colnames(Bio.raw) <- c("Samp_ID", "zpSR", "phySR", "zpD", "phyB", "Omnip", "Gzp", "mzpGzp",
                       "pH", "Temp", "Cond", "PAR", "TN", "TP")
panel.plot <- function(x, y) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  ct <- cor.test(x,y)
  sig <- symnum(ct$p.value, corr = FALSE, na = FALSE,
                cutpoints = c(0, 0.1, 1),
                symbols = c("*", " "))
  r <- round(ct$estimate, 2)
  #rt <- format(r, digits=2)[1]
  if(cor.test(x, y)$p.value < 0.05){
    text(0.5, 0.5, r, cex = 3, font = 2)
    text(0.75, 0.75, sig, cex = 3)
  }else {text(.5, .5, r, cex = 2)}
}
panel.smooth <- function (x, y) {
  points(x, y, cex = 1.5, pch = 16)
  if(cor.test(x, y)$p.value < 0.05){
    abline(lm(y~x), col = "black", lwd = 2)
  }else {abline(lm(y~x), col = "black", lty = "dashed")}
}

par(cex.axis = 2.5)
pairs(cbind(log(Bio.raw[,c("phyB", "phySR", "zpD", "zpSR")]), Bio.raw[, "Omnip"]), 
      labels = c(expression(atop(Phytoplankton, "biovolume*")), expression(atop(Phytoplankton, "species richness*")), 
                 expression(atop(Zooplankton, "density*")), expression(atop(Zooplankton, "species richness*")), 
                 expression(atop(Omnivorous, consumption))),
      cex.labels = 3,
      lower.panel = panel.smooth, upper.panel = panel.plot)
######### Fig S2_Field Bio correlation plots ##########

######### Fig S2_Field Envi correlation plots ##########
Dat_scale <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_scale.csv", 
                        sep = ",", header = TRUE)
Dat_raw <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_raw.csv", 
                      sep = ",", header = TRUE) %>%
  mutate(CM = log(All_gzp),
         phyB = log(phyB),
         phySR = log(phySR),
         zpD = log(zpD),
         zpSR = log(zpSR),
         Omni = log(-All_Omni), 
         Temp = log(Temp),
         Cond = log(Cond),
         ph = log(pH),
         PAR = log(PAR),
         TN = log(TN_mean),
         TP = log(TP_mean))
EnviCorr <- ggpairs(Dat_raw, 
                    column = c("Temp", "Cond", "pH", "PAR", "TN", "TP"),
                    columnLabels = c("Log(Temp)", "Log(Cond)", "Log(pH)", "Log(PAR)", "Log(TN)", "Log(TP)"),
                    upper = list(continuous = wrap("cor", size = 8)),
                    lower = list(continuous = "smooth")) + 
            theme(axis.line = element_line(colour = "black"),
                  # axis.title.x = element_text(size = 36, margin = margin(t = 12, r = 12, b = 0, l = 0, "pt"), hjust = 0.5),
                  # axis.title.y = element_text(size = 36, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")),
                  axis.text = element_text(size = 20))
# ggsave(EnviCorr, file = "D:/Manuscript/IGP_DivEffects_MS/MS_dissertation/MS_FieldExp_Ch3/Figs/FigS2_1_EnviCorr.pdf", 
#        width = 38, height = 24, units = c("cm"), dpi = 600)

envi_PC <- prcomp(Dat_scale[,c("TN", "TP", "Cond", "pH", "Temp", "PAR")], center = TRUE, scale. = TRUE)
EnviPCA12 <- autoplot(envi_PC, data = Dat.scale, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 10) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 12, r = 12, b = 12, l = 12, "pt"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 24, margin = margin(t = 12, r = 12, b = 0, l = 0, "pt"), hjust = 0.5),
        axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")),
        axis.text = element_text(size = 20))
EnviPCA13 <- autoplot(envi_PC, x = 1, y = 3, data = Dat.scale, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 10) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 12, r = 12, b = 12, l = 12, "pt"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 24, margin = margin(t = 12, r = 12, b = 0, l = 0, "pt"), hjust = 0.5),
        axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")),
        axis.text = element_text(size = 20))
EnviPCA23 <- autoplot(envi_PC, x = 2, y = 3, data = Dat.scale, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 10) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 12, r = 12, b = 12, l = 12, "pt"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 24, margin = margin(t = 12, r = 12, b = 0, l = 0, "pt"), hjust = 0.5),
        axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")),
        axis.text = element_text(size = 20))
EnviPCA <- plot_grid(EnviPCA12, EnviPCA13, EnviPCA23, scale = 0.95, ncol = 3, labels = c("b.", "c.", "d."), label_size = 30)
FinalCorr <- plot_grid(ggmatrix_gtable(EnviCorr), EnviPCA, ncol = 1, rel_heights = c(3, 1), labels = c("a.",""), label_size = 30)
# ggsave(FinalCorr, 
#        file = "D:/Manuscript/IGP_DivEffects_MS_Figs/Field_NLA_Figs/SupFigs/FieldNLA_fs2_2_Field_EnviCorr.pdf", 
#        width = 40, height = 40, units = c("cm"), dpi = 600)
######### Fig S2_Field Envi correlation plots ##########




