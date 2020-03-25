library(tidyverse)
library(lubridate)
library("ggpubr")
library(pwr)

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
  mutate(Sex = as.factor(Sex))%>%
  filter(!is.na(DNAmAge))%>%
  mutate(Sex = as.numeric (Sex))

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

#data visualization
my_graph <- ggplot(data = full_data) +
  geom_point(aes(x = DNAmAge, y =full_data$BothEXTHIPP_no36_tau), size = 3)
my_graph +
  scale_fill_brewer(palette="Set1")+
  xlab('Epigenetic Age') + 
  ylab('SUVR of the Hippocampus')+
  stat_smooth(aes(x = DNAmAge, y =full_data$BothEXTHIPP_no36_tau), method = 'lm', se = FALSE, color = 'black')+
  theme_classic()+
  theme (axis.text=element_text(size=25))+
  theme(axis.title.y = element_text(size = 30), axis.title.x = element_text(size = 30))

ggplot(full_data, aes(x = full_data$DNAmAge))+
  geom_histogram(bins = 10)

########################################################################################
#statistical analysis
full_data_lm <- full_data%>%
  select(contains("tau"))%>%
  select(1:32, )
full_data_lm <- data.frame(full_data_lm)
predict_value <- full_data%>%
  select(DNAmAge)#exchange with DNAmAge/age_discordance
predict_value <- as.numeric(unlist(predict_value))
my_lms <- lapply(1:32, function(x) lm(full_data_lm[,x]~predict_value + full_data$Sex + full_data$MMSETotal + full_data$Amyloid_Status))

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
pwr <- list()
n_num <- vector("double", ncol(full_data_lm))
  
for (i in seq_along(f2_num)) {
  pwr[[i]] <- pwr.f2.test(u = 2, f2 = f2_num[[i]] , sig.level = 0.05, power = 0.8)
  n_num[[i]] <- pwr[[i]]$v + 3
  
}

r_coeff <- data.frame(r_cor = apply(full_data_lm, 2, function(y) cor(y, predict_value)))
r_coeff <- as.numeric(unlist(r_coeff))

stats_df <- data.frame(name, r_squared, p_value, r_coeff, p_value_DNAm_effect, n_num)
stats_df <- stats_df%>%
  mutate(p_value_adj = p.adjust(stats_df$p_value, method = "fdr", n = length(p_value)))%>%
  mutate(p_value_DNAm_effect_adj = p.adjust(stats_df$p_value_DNAm_effect, method = "fdr", n = length(p_value_DNAm_effect)))

pwr <- pwr.f2.test(u = 2, f2 = 0.157, sig.level = 0.05, power = 0.8)