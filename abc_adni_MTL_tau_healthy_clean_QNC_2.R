library(tidyverse)
library(lubridate)
library("ggpubr")
library(pwr)
library(car)
library(corrplot)
library(lmtest)

#setting working directory
setwd("C:/Users/maron/OneDrive/Documents/Penn Documents/RotationIMcMillanLab/ABC_CSVs")

#PET-tau data
abc_mri <- read_csv('stats_lr_cleanup_sorted_2.13.20.csv')
abc_mri <- abc_mri%>%
  rename(INDDID = 'RID')%>%
  mutate(INDDID = as.character(INDDID))%>%
  filter(!is.na(left_HIPP_tau))%>%
  filter(left_HIPP_tau > 0)

abc_mri$INDDID <- gsub('_', '.', abc_mri$INDDID)

#diagnosis data (filtering normal)
abc_DX <- read_csv('DX_ABC_2_23_2020.csv')
abc_DX <- abc_DX%>%
  filter(grepl('Normal', GlobalDx))%>%
  mutate(INDDID = as.character(INDDID))

#merging data
full_data <- left_join(abc_mri, abc_DX, by = 'INDDID')
full_data <- full_data%>%
  filter(!is.na(GlobalDx))

#amyloid classification data (filtering negative)
abc_abeta <- read_csv('AmyloidPET_VisualReads_daw_IN.csv')
abc_abeta <- abc_abeta%>%
  mutate(INDDID = as.character(INDDID))%>%
  rename(Amyloid_Status = 'Amyloid Status')%>%
  filter(Amyloid_Status == 0)
abc_abeta <- abc_abeta[!duplicated(abc_abeta$INDDID), ]

#merging data
full_data <- left_join(full_data, abc_abeta, by = 'INDDID')
full_data <- full_data%>%
  filter(!is.na(Amyloid_Status))

#epigenetic age data
abc_mdna <- read_csv('ABCdnaage_18Feb2020.csv')
abc_mdna <- abc_mdna%>%
  mutate(DNAdate = mdy(DNAdate))%>%
  mutate(DNAmAge = as.numeric(DNAmAge))%>%
  mutate(INDDID = as.character(INDDID))%>%
  filter(!is.na(DNAdate))
abc_mdna <- abc_mdna[!duplicated(abc_mdna$INDDID), ]

#merging data
full_data <- left_join(full_data, abc_mdna, by = 'INDDID')
full_data <- full_data%>%
  mutate(INDDID  = as.character(INDDID ))%>%
  filter(!is.na(left_HIPP_tau))%>%
  filter(!is.na(DNAmAge))%>%
  mutate(Sex = as.factor(Sex))

#mini mental state exam data
abc_MMSE <- read_csv('MMSE_ABC_2_23_2020.csv')
abc_MMSE <- abc_MMSE%>%
  mutate(INDDID = as.character(INDDID))%>%
  mutate(TestDate = mdy(TestDate))

#merging and filtering for scores equal to or more than 27
full_data <- left_join(full_data, abc_MMSE, by = 'INDDID')
full_data <- full_data%>%
  filter(!is.na(MMSETotal))%>%
  mutate(datediff = abs(TestDate - DNAdate))%>%
  group_by(INDDID) %>%
  slice(which.min(datediff))%>%
  ungroup(INDDID)%>%
  filter(MMSETotal >= 27)

full_data$BothEXTHIPP_no36_tau <- rowMeans(full_data[c('left_EXTHIPPOno36_tau', 'right_EXTHIPPOno36_tau')], na.rm=TRUE)

full_data$Amy_class <- ifelse(full_data$Amyloid_Status == 0, 'Negative', 'Positive')
full_data <- full_data%>%
  mutate(Amy_class = as.factor(Amy_class))

#age linear model
abc_lm_age <- lm(DNAmAge ~ ActualAge + Sex, data = full_data)

#age discordance
smage <- summary(abc_lm_age)
full_data <- full_data%>%
  mutate(age_discordance = smage$residuals)

#removing outliers
full_data <- full_data%>%
  filter(DNAmAge > mean(DNAmAge)-3*sd(DNAmAge) & DNAmAge < mean(DNAmAge)+ 3*sd(DNAmAge))%>%
  filter(ActualAge > mean(ActualAge)-3*sd(ActualAge) & ActualAge < mean(ActualAge)+ 3*sd(ActualAge))%>%
  filter(left_HIPP_tau > mean(left_HIPP_tau)- 3 * sd(left_HIPP_tau) & left_HIPP_tau < mean(left_HIPP_tau)+ 3 * sd(left_HIPP_tau))

full_data_cormat <- full_data%>%
  select(DNAmAge, EXTHIPPno36both_tau)

#data visualization
my_graph <- ggplot(data = full_data_edited) +
  geom_point(aes(x = DNAmAge, y =full_data_edited$BothEXTHIPP_no36_tau), size = 3)
my_graph +
  scale_fill_brewer(palette="Set1")+
  xlab('Epigenetic Age') + 
  ylab('SUVR of the Hippocampus')+
  stat_smooth(aes(x = DNAmAge, y =full_data_edited$BothEXTHIPP_no36_tau), method = 'lm', se = TRUE)+
  theme_classic()+
  theme (axis.text=element_text(size=25))+
  theme(axis.title.y = element_text(size = 30), axis.title.x = element_text(size = 30))

########################################################################################
#statistical analysis
full_data_lm <- full_data%>%
  select(contains("tau"))%>%
  select(c(1:32, 157), )
full_data_lm <- data.frame(full_data_lm)
predict_value <- full_data%>%
  select(DNAmAge)#exchange with DNAmAge/age_discordance
predict_value <- as.numeric(unlist(predict_value))
my_lms <- lapply(1:33, function(x) lm(full_data_lm[,x]~predict_value + full_data$Sex))

#p value function
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

r_squared <- vector("double", ncol(full_data_lm))
p_value <- vector("double", ncol(full_data_lm))
p_value_DNAm_effect <- vector("double", ncol(full_data_lm))
name <- vector("character", ncol(full_data_lm))
r_coeff <- vector("double", ncol(full_data_lm))
f2_num <- vector("double", ncol(full_data_lm))

for (i in seq_along(full_data_lm)) {
  name[[i]] <- names(full_data_lm)[[i]]
  p_value[[i]] <- lmp(my_lms[[i]])
  r_squared[[i]] <- summary(my_lms[[i]])$r.squared
  p_value_DNAm_effect[[i]] <- summary(my_lms[[i]])$coefficients[2,4]
  f2_num[[i]] <- summary(my_lms[[i]])$r.squared / (1 - summary(my_lms[[i]])$r.squared)
  
}
#power analysis
pwr <- list()
n_num <- vector("double", ncol(full_data_lm))
  
for (i in seq_along(f2_num)) {
  pwr[[i]] <- pwr.f2.test(u = 2, f2 = f2_num[[i]] , sig.level = 0.05, power = 0.8)
  n_num[[i]] <- pwr[[i]]$v + 3
  
}

r_coeff <- data.frame(r_cor = apply(full_data_lm, 2, function(y) cor(y, predict_value)))
r_coeff <- as.numeric(unlist(r_coeff))

stats_df_MTL <- data.frame(name, r_squared, p_value, r_coeff, p_value_DNAm_effect, n_num)
stats_df_MTL <- stats_df_MTL%>%
  mutate(p_value_adj = p.adjust(stats_df$p_value, method = "fdr", n = length(p_value)))%>%
  mutate(p_value_DNAm_effect_adj = p.adjust(stats_df$p_value_DNAm_effect, method = "fdr", n = length(p_value_DNAm_effect)))

#power analysis
pwr.f2.test(u = 2, f2 = 0.0917 , sig.level = 0.05, power = 0.8)

setwd("C:/Users/maron/OneDrive/Documents/Penn Documents/RotationIMcMillanLab/ADNI_CSVs")

#cleaning data
adni_tau <- read_csv('MRI_cog_data_1-14-19-1_noNA_nocomma.csv')
adni_tau <- adni_tau%>%
  filter(!is.na(NUM))%>%
  filter(DXCURREN == 1)%>% 
  filter(MMSCORE >= 27)%>%
  mutate(TAUSCANDATE = mdy(TAUSCANDATE))%>%
  mutate(RID = as.character(RID))%>%
  mutate(DXCURREN = as.factor(DXCURREN))%>%
  mutate(PTDOB = mdy(PTDOB))%>%
  filter(!is.na(left_CA1_tau))%>%
  filter(!is.na(Right_FO_frontal_operculum_thickness))

adni_mdna <- read_csv('ADNI_mDNAage.csv')
adni_mdna <- adni_mdna%>%
  filter(!is.na(RID))%>%
  mutate(mDNA_Date = mdy(mDNA_Date))%>%
  select(mDNA_Date, RID, DNAmAge)%>%
  mutate(DNAmAge = as.numeric(DNAmAge))%>%
  mutate(RID = as.character(RID))%>%
  filter(!is.na(mDNA_Date))

full_data <- left_join(adni_mdna, adni_tau, by = 'RID')
full_data <- full_data%>%
  mutate(datediff = abs(TAUSCANDATE - mDNA_Date))%>%
  group_by(RID) %>%
  slice(which.min(datediff))%>%
  ungroup(RID)%>%
  mutate(RID = as.character(RID))%>%
  filter(!is.na(left_CA1_tau))

adni_AV45 <- read_csv('UCBERKELEYAV45_08_27_19.csv')
adni_AV45 <- adni_AV45%>%
  mutate(EXAMDATE = mdy(EXAMDATE))%>%
  select(EXAMDATE, RID, SUMMARYSUVR_WHOLECEREBNORM, SUMMARYSUVR_COMPOSITE_REFNORM)%>%
  mutate(RID = as.character(RID))

full_data <- left_join(full_data, adni_AV45, by = 'RID')
full_data <- full_data%>%
  mutate(datediff_mDNA = abs(EXAMDATE - mDNA_Date))%>%
  group_by(RID) %>%
  slice(which.min(datediff_mDNA))%>%
  ungroup(RID)%>%
  mutate(RID = as.character(RID))

adni_demo <- read_csv("ADNIMERGE.csv")
adni_demo <- adni_demo%>%
  select(RID, PTGENDER, PTEDUCAT, EXAMDATE, EXAMDATE_bl, AV45_bl)%>%
  mutate(EXAMDATE = ymd(EXAMDATE))%>%
  group_by(RID) %>%
  slice(which.max(EXAMDATE))%>%
  ungroup(RID)%>%
  mutate(RID = as.character(RID))%>%
  mutate(PTGENDER = as.factor(PTGENDER))%>%
  mutate(Sex = as.numeric (PTGENDER))

full_data <- left_join(full_data, adni_demo, by = 'RID')
full_data <- full_data%>%
  mutate(age_calc = (mDNA_Date-PTDOB)/365.25)%>%
  mutate(age_calc = as.numeric(age_calc))%>%
  filter(!is.na(age_calc))%>%
  mutate(datediff_bl_tau = abs(EXAMDATE_bl - TAUSCANDATE))%>%
  mutate(datediff_bl_mDNA = abs(EXAMDATE_bl - mDNA_Date))%>%
  filter(SUMMARYSUVR_WHOLECEREBNORM <= 1.11)%>%
  filter(!is.na(MMSCORE))%>%
  mutate(Sex = as.factor(Sex))

#age linear model
adni_lm_age <- lm(DNAmAge ~ age_calc + Sex + datediff, data = full_data)

#age discordance
smage <- summary(adni_lm_age)
full_data <- full_data%>%
  mutate(age_discordance = smage$residuals)

#removing outliers
full_data <- full_data%>%
  filter(cereb_tau > mean(cereb_tau)-3*sd(cereb_tau) & cereb_tau < mean(cereb_tau)+ 3*sd(cereb_tau))%>%
  filter(DNAmAge > mean(DNAmAge)-3*sd(DNAmAge) & DNAmAge < mean(DNAmAge)+ 3*sd(DNAmAge))%>%
  filter(age_calc > mean(age_calc)-3*sd(age_calc) & age_calc < mean(age_calc)+ 3*sd(age_calc))
full_data_lm <- full_data%>%
  select(contains("thickness"), contains("_tau"), -Slice_Thickness, DNAmAge, age_calc, Sex, datediff, age_discordance, PTEDUCAT, SUMMARYSUVR_WHOLECEREBNORM)
full_data_lm <- full_data_lm%>%
  mutate(cort_thk_mean = rowMeans(full_data_lm[,1:104]))%>%
  filter(cort_thk_mean > mean(cort_thk_mean)-3*sd(cort_thk_mean) & cort_thk_mean < mean(cort_thk_mean)+ 3*sd(cort_thk_mean))
full_data <- full_data_lm%>%
  mutate(datediff_year = datediff/365.25)

full_data_cormat <- full_data%>%
  select(DNAmAge, datediff)%>%
  mutate(datediff = as.numeric(datediff))

#statistical analysis
full_data_lm <- full_data%>%
  select(contains("_tau"), -Left_Hippocampus_tau, -Right_Hippocampus_tau)
full_data_lm <- data.frame(full_data_lm)

#slecting hipppocampal subfields
full_data_lm <- full_data_lm%>%
  select(1:35)
predict_value <- full_data%>%
  select(DNAmAge)#exchange with DNAmAge/age_discordance
predict_value <- as.numeric(unlist(predict_value))
my_lms <- lapply(1:35, function(x) lm(full_data_lm[,x]~predict_value + full_data$Sex +full_data$datediff))

#p value function
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

r_squared <- vector("double", ncol(full_data_lm))
p_value <- vector("double", ncol(full_data_lm))
p_value_DNAm_effect <- vector("double", ncol(full_data_lm))
name <- vector("character", ncol(full_data_lm))
r_coeff <- vector("double", ncol(full_data_lm))

for (i in seq_along(full_data_lm)) {
  name[[i]] <- names(full_data_lm)[[i]]
  p_value[[i]] <- lmp(my_lms[[i]])
  r_squared[[i]] <- summary(my_lms[[i]])$r.squared
  p_value_DNAm_effect[[i]] <- summary(my_lms[[i]])$coefficients[2,4]
}

r_coeff <- data.frame(r_cor = apply(full_data_lm, 2, function(y) cor(y, predict_value)))
r_coeff <- as.numeric(unlist(r_coeff))

stats_df <- data.frame(name, r_squared, p_value, r_coeff, p_value_DNAm_effect)
stats_df <- stats_df%>%
  mutate(p_value_adj = p.adjust(stats_df$p_value, method = "fdr", n = length(p_value)))%>%
  mutate(p_value_DNAm_effect_adj = p.adjust(stats_df$p_value_DNAm_effect, method = "fdr", n = length(p_value_DNAm_effect)))

hipp_lm <- lm(EXTHIPPno36both_tau ~ DNAmAge + Sex +datediff, data = full_data)

#qq plot
qqPlot(hipp_lm$residuals)

#correlation matrix
M <- cor(full_data_cormat)
corrplot(M, method="color")

res.aov <- aov(full_data$datediff ~ full_data$Sex, data = full_data)
# Summary of the analysis
summary(res.aov)

res.aov <- aov(full_data$DNAmAge ~ full_data$Sex, data = full_data)
# Summary of the analysis
summary(res.aov)

#durbin-watson test
dwtest(hipp_lm)
