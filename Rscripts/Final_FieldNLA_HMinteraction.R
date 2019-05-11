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

######## Preping the data ###################################################################################################
Field.scale <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_scale.csv", 
                          sep = ",", header = TRUE) 
envi_PC <- prcomp(Field.scale[,c("TN", "TP", "Cond", "pH", "Temp", "PAR")], center = TRUE, scale. = TRUE)
#summary(envi_PC)
enviPC1 <- envi_PC$x[,"PC1"]
enviPC2 <- envi_PC$x[,"PC2"]
enviPC3 <- envi_PC$x[,"PC3"]

Bio.raw <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_raw.csv", sep = ",", header = TRUE) %>%
  subset(select = c(Samp_ID, zpSR, phySR, zpD, phyB, All_Omni_prop, All_gzp, mzp_gzp, pH, Temp, Cond, PAR, TN_mean, TP_mean)) %>%
  mutate(Samp_ID = as.character(Samp_ID))
colnames(Bio.raw) <- c("Samp_ID", "zpSR", "phySR", "zpD", "phyB", "Omnip", "Gzp", "mzpGzp",
                       "pH", "Temp", "Cond", "PAR", "TN", "TP")

Bio.scale <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/FieldDat_scale.csv", 
                        sep = ",", header = TRUE) %>%
  mutate(Envi1 = enviPC1, Envi2 = enviPC2, Envi3 = enviPC3) %>%
  subset(select = c(Samp_ID, zpSR, phySR, zpD, phyB, All_Omni_prop, All_gzp, mzp_gzp, 
                    pH, Temp, Cond, PAR, TN_mean, TP_mean, Envi1, Envi2, Envi3)) %>%
  mutate(Samp_ID = as.character(Samp_ID))
colnames(Bio.scale) <- c("Samp_ID", "zpSR", "phySR", "zpD", "phyB", "Omnip", "Gzp", "mzpGzp",
                         "pH", "Temp", "Cond", "PAR", "TN", "TP", "Envi1", "Envi2", "Envi3")

ZP_ID <- read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/Field_zpCount.csv", 
                  sep = ",", header = TRUE, fill = TRUE) %>%
  mutate(Samp_ID = as.character(Samp_ID)) %>%
  subset(Samp_ID %in% Bio.scale$Samp_ID)
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
  right_join(Bio.scale, by = "Samp_ID") %>%
  right_join(read.csv(file = "https://raw.githubusercontent.com/OscarFHC/OmniBEF_FieldNLA_public/master/Field_zpCount.csv", 
                      sep = ",", header = TRUE, fill = TRUE) %>%
               mutate(Samp_ID = as.character(Samp_ID)) %>%
               subset(Samp_ID %in% Bio.scale$Samp_ID) %>%
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
    ZPcom[, paste0("domzp", i)] <- rowSums(ZPcom[, domTSN])/(Bio.raw[, "zpD"]*1000)
  }else{ZPcom[, paste0("domzp", i)] <- ZPcom[, domTSN]/(Bio.raw[, "zpD"]*1000)}
}

ZPcom <- ZPcom %>%
  mutate(interact = zpSR*Omnip)
######## Preping the data ###################################################################################################

######## HM model ###########################################################################################################
##### running HM model #####
mm_Gzp <- model.matrix(~ zpSR + phyB + zpD + Envi1 + Envi2 + Envi3, data = ZPcom)
mm_b_intn = model.matrix(~ domzp1 + domzp2 + domzp3 + domzp4 + domzp5 + domzp6 + domzp7 + domzp8 + domzp9 + domzp10, data = ZPcom)

HM_Omni.dat <- # data list for Not region-denepdent but effects of zpSR depends on Omni analyses
  list(N = nrow(mm_Gzp), # number of obs
       K_Gzp = ncol(mm_Gzp), # number of column of the model matrix for determining Gz
       X_Gzp = mm_Gzp, # model matrix determining  Gz
       Gzp = ZPcom[,"Gzp"],
       
       intn = ZPcom[,"zpSR"],
       K_bintn = ncol(mm_b_intn), # number of column of the model matrix for determining the effects of zpSR
       X_bintn = mm_b_intn # model matrix determining the effects of zpSR on Gz
       )

HM_Omni = stan(file = "D:/Research/OmniBEF_FieldNLA_public/stan/Field_Omni_201905.stan", 
                       data = HM_Omni.dat,
                       iter = 6000, warmup = 2000, chains = 4,
                       control = list(adapt_delta = 0.95, max_treedepth = 15))
##### running HM model #####
##### manually summariz paramaters #####
param_Omni = extract(HM_Omni)
beta.sum = as.data.frame(param_Omni[["gamma_Gzp"]]) %>%
  gather(key = "gamma_Gzp_n", value = "value") %>%
  group_by(gamma_Gzp_n) %>%
  summarize(
    mean=mean(value),
    lo95=sort(value)[nrow(as.data.frame(param_Omni[["gamma_Gzp"]]))*0.025],
    lo90=sort(value)[nrow(as.data.frame(param_Omni[["gamma_Gzp"]]))*0.05],
    hi90=sort(value)[nrow(as.data.frame(param_Omni[["gamma_Gzp"]]))*0.95],
    hi95=sort(value)[nrow(as.data.frame(param_Omni[["gamma_Gzp"]]))*0.975]
  ) %>%
  mutate(gamma_Gzp_n = c("Intercept", "ln_zpSR", "Omni", "ln_phyB", "ln_zpD", "Envi1", "Envi2", "Envi3", "interaction"),
         beta = gamma_Gzp_n) %>%
  select(-gamma_Gzp_n)
write.table(beta.sum, file = "D:/Research/OmniBEF_FieldNLA_public/stanFieldNLA_Omni_beta.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)

b.sum = as.data.frame(param_Omni[["gamma_beta_intn"]]) %>%
  gather(key = "gamma_beta_intn", value = "value" ) %>%
  group_by(gamma_beta_intn) %>%
  summarize(
    mean=mean(value),
    lo95=sort(value)[nrow(as.data.frame(param_Omni[["gamma_beta_intn"]]))*0.025],
    lo90=sort(value)[nrow(as.data.frame(param_Omni[["gamma_beta_intn"]]))*0.05],
    hi90=sort(value)[nrow(as.data.frame(param_Omni[["gamma_beta_intn"]]))*0.95],
    hi95=sort(value)[nrow(as.data.frame(param_Omni[["gamma_beta_intn"]]))*0.975]
  ) %>%
  mutate(gamma_beta_intn = c("Intercept", "domzp1", "domzp2", "domzp3", "domzp4", "domzp5", 
                             "domzp6", "domzp7", "domzp8", "domzp9", "domzp10"),
         beta = gamma_beta_intn) %>%
  select(-gamma_beta_intn)
write.table(b.sum, file = "D:/Research/OmniBEF_FieldNLA_public/stanFieldNLA_Omni_b.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)

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
######## HM model ###########################################################################################################

glm_f  <- glm(Gzp ~ zpSR*(domzp1 + domzp2 + domzp3 + domzp4 + domzp5 + domzp6 + domzp7 + domzp8 + domzp9 + domzp10) + 
                    phyB + zpD + Envi1 + Envi2 + Envi3, data = ZPcom)
summary(glm_f)

ZPcom <-  ZPcom %>% mutate(domzp10q = domzp10^2)
glm_single  <- glm(Gzp ~ zpSR*(domzp10) + phyB + zpD + Envi1 + Envi2 + Envi3, data = ZPcom)
summary(glm_single)

plot(Gzp ~ domzp10, data = ZPcom)

mod <- "
  Gzp ~ zpSR + domzp10 + zpSR:domzp10 + zpD + phyB
  zpSR ~ Envi1 + Envi2 + Envi3
  zpD ~ Envi1 + Envi2 + Envi3
  phyB ~ Envi1 + Envi2 + Envi3
  
  zpD ~~ phyB
"
mod_lavaan <- sem(mod, data = ZPcom, meanstructure = TRUE)
summary(mod_lavaan)
fitMeasures(mod_lavaan)
ZPcom$domzp10
