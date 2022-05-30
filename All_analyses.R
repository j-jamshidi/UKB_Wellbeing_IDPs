###Load data of participants with ICV available (ukb N# 42809)====
setwd("Z:/UKB/p_analysis/R/Brain_analysis/")
ukb <- read.csv("brain_available_samples_DK_plusVol.csv", header = TRUE)
ukb <- as.data.frame(ukb)
library(tidyverse)
library(psychometric)
library(fmsb) #calculate the 95% intervals around R2
library(lm.beta) #to standardize the estimate (in case!)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggseg)
library(ggseg3d)
library(plotly)
library(gtsummary)
library(viridis)
#deal with missing data&change the codes
ukb[ukb == -1] <- NA
ukb[ukb == -3] <- NA
ukb[ukb == -818] <- NA
ukb[ukb == -121] <- NA
ukb[ukb == -9] <- NA
ukb[ukb == " "] <- NA
#----Remove participant with confounding diseases, both at baseline or imaging visit; 1814 removed (ukb N#40993)====
#non-cancer
ukb[46:108][is.na(ukb[46:108])] <- 0
for (i in c(46:108)) {
  ukb <- ukb[ukb[i]  != 1081 & ukb[i]  != 1082 & ukb[i]  != 1083  & ukb[i]  != 1086  & ukb[i]  != 1240  & ukb[i]  != 1244 &
               ukb[i]  != 1245  & ukb[i]  != 1246  & ukb[i]  != 1247 & ukb[i]  != 1258 & ukb[i]  != 1259 & ukb[i]  != 1261 &
               ukb[i]  != 1262 & ukb[i]  != 1236  & ukb[i]  != 1264 & ukb[i]  != 1266 & ukb[i]  != 1425 & ukb[i]  != 1433 &
               ukb[i]  != 1434 & ukb[i]  != 1491 & ukb[i]  != 1524 & ukb[i]  != 1583 & ukb[i]  != 1659 & ukb[i]  != 1683 , ]
}
#cancer
ukb[109:119][is.na(ukb[109:119])] <- 0
for (i in c(109:119)) {
  ukb <- ukb[ukb[i]  != 1031 & ukb[i]  != 1032 & ukb[i]  != 1029, ]
}

#----Remove outliers -+3SD for ICV (ukb N#40885)====
ukb$ICV <- scale(ukb$ICV)
ukb$ICV[ukb$ICV >= 3] <- NA
ukb$ICV[ukb$ICV  <= (-3)] <- NA
ukb <- ukb[complete.cases(ukb$ICV),]
#----Remove people with missing covariates (ukbb N#40087)====
ukbb <- ukb
#remove missing education N#40345--
ukbb <- ukbb[complete.cases(ukbb$Education_2),]
#remove missing ethnicity N#40255---
ukbb <- ukbb[complete.cases(ukbb$Ethnicity_W_nW),]
#remove missing deprivation index N#40217--
ukbb <- ukbb[complete.cases(ukbb$Townsend_deprivation_index_at_recruitment),]
#Smoking status--N#40096---
ukbb <- ukbb[complete.cases(ukbb$Smoking_status_img),]
#Alcohol intake freq--#4087---
ukbb <- ukbb[complete.cases(ukbb$Alcohol_intake_frequency_img),]
#set categorical variables
ukbb$Education_2 <- as.factor(ukbb$Education_2)
ukbb$Ethnicity_W_nW <- as.factor(ukbb$Ethnicity_W_nW)
ukbb$Sex <- as.factor(ukbb$Sex)
ukbb$cat_BMI_imp <- as.factor(ukbb$cat_BMI_imp)
ukbb$Smoking_status_img <- as.factor(ukbb$Smoking_status_img)
ukbb$Alcohol_intake_frequency_img <- as.factor(ukbb$Alcohol_intake_frequency_img)
ukbb$Assessment_centre_Imaging <- as.factor(ukbb$Assessment_centre_Imaging)
#remove missing for BMI--
#######---------------------------------------------------------------------------------------------##
#----standardize wb_index2 and, remove 3SD outliers, select participants with wb_index 2 available (ukbbw N#38982)======
ukbb$wb_index_2 <- as.numeric(scale(ukbb$wb_index_2))
ukbb$wb_index_2[ukbb$wb_index_2 > 3] <- NA
ukbb$wb_index_2[ukbb$wb_index_2 < (-3)] <- NA
ukbbw <- ukbb[complete.cases(ukbb$wb_index_2),]  
ukbbw$wb_index_2 <- as.numeric(scale(ukbbw$wb_index_2))

#----standardize prs and, remove 3SD outliers, select participants with prs available (ukbbp N#19987)======
ukbb$PRS_Brain <- as.numeric(scale(ukbb$PRS_Brain))
ukbb$PRS_Brain[ukbb$PRS_Brain > 3] <- NA
ukbb$PRS_Brain[ukbb$PRS_Brain < (-3)] <- NA
ukbbp <- ukbb[complete.cases(ukbb$PRS_Brain),] 
ukbb$PRS_Brain <- scale(ukbb$PRS_Brain)
################################################## Wellbeing index analysis ###########################################.=====
#----standardize brain values and mark 3SD outliers as NA & plot outliers-----
ukbbw[161:205] <- as.numeric(scale(ukbbw[161:205])) #volumes
ukbbw[327:528] <- as.numeric(scale(ukbbw[327:528])) #SA and CT and volume from DK
ukbbw[161:205][ukbbw[161:205] > 3] <- NA
ukbbw[161:205][ukbbw[161:205] < (-3)] <- NA
ukbbw[327:528][ukbbw[327:528] > 3] <- NA
ukbbw[327:528][ukbbw[327:528] < (-3)] <- NA
#----plot the number of outliers (volume)
ss<- as.data.frame(summary(is.na(ukbbw[161:205]))) 
ss$Var1 <-NULL
ss<- ss[!grepl("Mode", ss$Freq),]
ss<- ss[!grepl("FALSE", ss$Freq),]
ss$Freq <- str_remove_all(ss$Freq, "TRUE :")
colnames(ss) <- c("variable", "N_missing")
ss$N_missing <- as.numeric(ss$N_missing)
volume_3SDw_Otliers <- ss %>% 
  mutate_at("variable", str_replace, "vgm_", "") %>%
  mutate_at("variable", str_replace, "Volume_of_", "") %>%
  mutate(across("variable", str_replace_all, "_", " ")) %>%
  mutate_at("variable", str_replace, "left","lh") %>%
  mutate_at("variable", str_replace, "right","rh") %>%
  mutate_at("variable", str_replace, "grey","Total grey") %>%
  mutate_at("variable", str_replace, "white","Total white") 
#plot missing and put them in a dataframe
volume_3SDw_Otliers$variable <- as.factor(volume_3SDw_Otliers$variable)
volume_3SDw_Otliers$variable <- factor(volume_3SDw_Otliers$variable, levels = volume_3SDw_Otliers$variable[order(volume_3SDw_Otliers$N_missing)])
V_3SDw_Otliers<-ggplot(volume_3SDw_Otliers , aes(y=variable , x = N_missing )) +
  geom_bar(stat = "identity", fill=rgb(0.2,0.4,0.6), width = .6) + theme_classic() +
  scale_x_continuous(limits = c(0,410), expand = c(0, 0))+ 
  theme(panel.border = element_blank(), axis.line.y = element_blank()) + labs(x="#N +/-3SD", title = "Volume")
#plot the number of outliers (SA)
ss<- as.data.frame(summary(is.na(ukbbw[327:393]))) 
ss$Var1 <-NULL
ss<- ss[!grepl("Mode", ss$Freq),]
ss<- ss[!grepl("FALSE", ss$Freq),]
ss$Freq <- str_remove_all(ss$Freq, "TRUE :")
colnames(ss) <- c("variable", "N_missing")
ss$N_missing <- as.numeric(ss$N_missing)
SA_3SDw_Otliers <- ss %>%
  mutate_at("variable", str_replace, "SA_", "") %>%
  mutate(across("variable", str_replace_all, "_", " "))
#plot missing and put them in a dataframe
SA_3SDw_Otliers$variable <- as.factor(SA_3SDw_Otliers$variable)
SA_3SDw_Otliers$variable <- factor(SA_3SDw_Otliers$variable, levels = SA_3SDw_Otliers$variable[order(SA_3SDw_Otliers$N_missing)])
S_3SDw_Otliers<-ggplot(SA_3SDw_Otliers , aes(y=variable , x = N_missing )) +
  geom_bar(stat = "identity", fill=rgb(0.2,0.4,0.6), width = .6) + theme_classic() +
  scale_x_continuous(limits = c(0,810), expand = c(0, 0))+ 
  theme(panel.border = element_blank(), axis.title.y = element_blank(),axis.line.y = element_blank())+
  labs(x="#N +/-3SD", title = "Surface area")
#plot the number of outliers (CT)
ss<- as.data.frame(summary(is.na(ukbbw[396:462]))) 
ss$Var1 <-NULL
ss<- ss[!grepl("Mode", ss$Freq),]
ss<- ss[!grepl("FALSE", ss$Freq),]
ss$Freq <- str_remove_all(ss$Freq, "TRUE :")
colnames(ss) <- c("variable", "N_missing")
ss$N_missing <- as.numeric(ss$N_missing)
CT_3SDw_Otliers <- ss %>%
  mutate_at("variable", str_replace, "CT_", "") %>%
  mutate(across("variable", str_replace_all, "_", " "))
#plot missing and put them in a dataframe
CT_3SDw_Otliers$variable <- as.factor(CT_3SDw_Otliers$variable)
CT_3SDw_Otliers$variable <- factor(CT_3SDw_Otliers$variable, levels = CT_3SDw_Otliers$variable[order(CT_3SDw_Otliers$N_missing)])
c_3SDw_Otliers<-ggplot(CT_3SDw_Otliers , aes(y=variable , x = N_missing )) +
  geom_bar(stat = "identity", fill=rgb(0.2,0.4,0.6), width = .6) + theme_classic() +
  scale_x_continuous(limits = c(0,810), expand = c(0, 0))+ 
  theme(panel.border = element_blank(), axis.title.y = element_blank(),axis.line.y = element_blank())+
  labs(x="#N +/-3SD", title = "Cortical thickness")
#plot the number of outliers (VOL_DK)
ss<- as.data.frame(summary(is.na(ukbbw[463:528]))) 
ss$Var1 <-NULL
ss<- ss[!grepl("Mode", ss$Freq),]
ss<- ss[!grepl("FALSE", ss$Freq),]
ss$Freq <- str_remove_all(ss$Freq, "TRUE :")
colnames(ss) <- c("variable", "N_missing")
ss$N_missing <- as.numeric(ss$N_missing)
Vo_3SDw_Otliers <- ss %>%
  mutate_at("variable", str_replace, "V_", "") %>%
  mutate(across("variable", str_replace_all, "_", " "))
#plot missing and put them in a dataframe
Vo_3SDw_Otliers$variable <- as.factor(Vo_3SDw_Otliers$variable)
Vo_3SDw_Otliers$variable <- factor(Vo_3SDw_Otliers$variable, levels = Vo_3SDw_Otliers$variable[order(Vo_3SDw_Otliers$N_missing)])
Vo_3SDw_Otliers<-ggplot(Vo_3SDw_Otliers , aes(y=variable , x = N_missing )) +
  geom_bar(stat = "identity", fill=rgb(0.2,0.4,0.6), width = .6) + theme_classic() +
  scale_x_continuous(limits = c(0,810), expand = c(0, 0))+ 
  theme(panel.border = element_blank(), axis.title.y = element_blank(),axis.line.y = element_blank())+
  labs(x="#N +/-3SD", title = "Cortical volume")
#put all plots together
ggarrange(V_3SDw_Otliers,Vo_3SDw_Otliers,c_3SDw_Otliers,S_3SDw_Otliers,
          ncol = 4, nrow = 1, widths = c(1.3,2,2,2))
#----restandardise values------
ukbbw$wb_index_2 <- as.numeric(scale(ukbbw$wb_index_2))
ukbbw[161:205] <- as.numeric(scale(ukbbw[161:205]))
ukbbw[327:528] <- as.numeric(scale(ukbbw[327:528]))
#Quick check of normality===== 
par(mfrow= c(5,6))
for (i in c(161:205)){
  hist(ukbbw[,i], xlab = colnames(ukbbw[i]), breaks= 40, main = "", ylab = "")
} 
par(mfrow= c(5,6))
for (i in c(327:528)){
  hist(ukbbw[,i], xlab = colnames(ukbbw[i]), breaks= 40, main = "", ylab = "")
}
###--------------------------Linear models for basic covariates------------------
##linear models for  Brain VOLUME from  DK,FIRST, FAST-------
bv_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(bv_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(161:205,463:528)) {
  fitm <-lm.beta(lm( ukbbw[,i]~ wb_index_2+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+ 
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ ICV
                     , ukbbw)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbw[,i])) # Sample size of the test
  bv_ROIs[nrow(bv_ROIs)+1,] <- c(colnames(ukbbw)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"vol" ) #save in a dataframe
}
bv_ROIs$pFDR <- p.adjust(bv_ROIs$pval, method = "BH", n = nrow(bv_ROIs))
bv_ROIs$pbonf <- p.adjust(bv_ROIs$pval, method = "bonferroni", 245)
##linear models for CT-------
ct_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(ct_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval", "group")
#loop for models
for (i in c(396:462)) {
  fitm <-lm.beta(lm( ukbbw[,i]~ wb_index_2+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+ 
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ICV
                     , ukbbw)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbw[,i])) # Sample size of the test
  ct_ROIs[nrow(ct_ROIs)+1,] <- c(colnames(ukbbw)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"CT" ) #save in a dataframe
}
ct_ROIs$pFDR <- p.adjust(ct_ROIs$pval, method = "BH", n = nrow(ct_ROIs))
ct_ROIs$pbonf <- p.adjust(ct_ROIs$pval, method = "bonferroni", 245)

##linear models for  SA-------
sa_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(sa_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(327:393)) {
  fitm <-lm.beta(lm( ukbbw[,i]~ wb_index_2+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+ 
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ ICV
                     , ukbbw)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbw[,i])) # Sample size of the test
  sa_ROIs[nrow(sa_ROIs)+1,] <- c(colnames(ukbbw)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"SA" ) #save in a dataframe
}
sa_ROIs$pFDR <- p.adjust(sa_ROIs$pval, method = "BH", n = nrow(sa_ROIs))
sa_ROIs$pbonf <- p.adjust(sa_ROIs$pval, method = "bonferroni", 245)

###--------------------------Linear models with full covariates------------------
##linear models for  Brain VOLUMES from DK,FIRST,FAST-------
bvf_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(bvf_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(161:205,463:528)) {
  fitm <-lm.beta(lm( ukbbw[,i]~ wb_index_2+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+cat_BMI_imp+
                       Education_2+Ethnicity_W_nW+Townsend_deprivation_index_at_recruitment+Smoking_status_img+Alcohol_intake_frequency_img+
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ ICV
                     , ukbbw)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbw[,i])) # Sample size of the test
  bvf_ROIs[nrow(bvf_ROIs)+1,] <- c(colnames(ukbbw)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"vol" ) #save in a dataframe
}
bvf_ROIs$pFDR <- p.adjust(bvf_ROIs$pval, method = "BH", n = nrow(bvf_ROIs))
bvf_ROIs$pbonf <- p.adjust(bvf_ROIs$pval, method = "bonferroni", 245)

##linear models for CT-------
ctf_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(ctf_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(396:462)) {
  fitm <-lm.beta(lm( ukbbw[,i]~ wb_index_2+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+cat_BMI_imp+
                       Education_2+Ethnicity_W_nW+Townsend_deprivation_index_at_recruitment+Smoking_status_img+Alcohol_intake_frequency_img+
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ ICV
                     , ukbbw)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbw[,i])) # Sample size of the test
  ctf_ROIs[nrow(ctf_ROIs)+1,] <- c(colnames(ukbbw)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"CT" ) #save in a dataframe
}
ctf_ROIs$pFDR <- p.adjust(ctf_ROIs$pval, method = "BH", n = nrow(ctf_ROIs))
ctf_ROIs$pbonf <- p.adjust(ctf_ROIs$pval, method = "bonferroni", 245)

##linear models for  SA-------
saf_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(saf_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(327:393)) {
  fitm <-lm.beta(lm( ukbbw[,i]~ wb_index_2+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+cat_BMI_imp+
                       Education_2+Ethnicity_W_nW+Townsend_deprivation_index_at_recruitment+Smoking_status_img+Alcohol_intake_frequency_img+
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ ICV
                     , ukbbw)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbw[,i])) # Sample size of the test
  saf_ROIs[nrow(saf_ROIs)+1,] <- c(colnames(ukbbw)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"SA" ) #save in a dataframe
}
saf_ROIs$pFDR <- p.adjust(saf_ROIs$pval, method = "BH", n = nrow(saf_ROIs))
saf_ROIs$pbonf <- p.adjust(saf_ROIs$pval, method = "bonferroni", 245)

#Put all analysis in a single dataframe-------------
#basic models
all_ROIs <- bv_ROIs
all_ROIs[112:178,] <- ct_ROIs
all_ROIs[179:245,] <- sa_ROIs
#with all covariates
allf_ROIs <- bvf_ROIs
allf_ROIs[112:178,] <- ctf_ROIs
allf_ROIs[179:245,] <- saf_ROIs
###Plotting time!=======================================================
#BASIC models; all regions that survive FDR-P-value of 0.0166 (0.05/3)=======
s <- all_ROIs
s[2:8,10:11] <- as.numeric(unlist(s[2:8,10:11]))
s <- s[s$pFDR < 0.0166,]
s$estimate <- as.numeric(s$estimate)
s$lci_estimate <- as.numeric(s$lci_estimate)
s$uci_estimate <- as.numeric(s$uci_estimate)
plot_all_ROIs<- ggplot() +
  geom_pointrange(data=s, aes(x=variable, y=estimate, ymin=lci_estimate, ymax=uci_estimate,color=group),shape =15, size=.8)  +
  ggtitle("Brain ~ WB +sex+age+age2+assescentre+ICV+scanerXYZ")+
  geom_hline(yintercept=0, lty=2, color ="gray20") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("standardised estimate (95% CI)") + xlab("")+
  theme_classic()+
  theme(text =element_text(size=14),axis.text = element_text(size=14),plot.title = element_text(color="gray20", size=12, face="italic"))
plot_all_ROIs
#FULL models; all regions that survive FDR-P-value of 0.0166 (0.05/3)=======
sf <- allf_ROIs
sf[2:8,10:11] <- as.numeric(unlist(sf[2:8,10:11]))
sf <- sf[sf$pFDR < 0.0166,]
sf$estimate <- as.numeric(sf$estimate)
sf$lci_estimate <- as.numeric(sf$lci_estimate)
sf$uci_estimate <- as.numeric(sf$uci_estimate)
plot_allf_ROIs<- ggplot() +
  geom_pointrange(data=sf, aes(x=variable, y=estimate, ymin=lci_estimate, ymax=uci_estimate,color=group),shape =15, size=.8)  +
  ggtitle("Brain ~ WB +sex+age+age2+assescentre+ICV+scanerXYZ+ethn+edu+depindex+alchol+smoking+BMI")+
  geom_hline(yintercept=0, lty=2, color ="gray20") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("standardised estimate (95% CI)") + xlab("")+
  theme_classic()+
  theme(text =element_text(size=14),axis.text = element_text(size=14),plot.title = element_text(color="gray20", size=12, face="italic"))
plot_allf_ROIs

#A plot with shows which variables still significant after adding all covariates===== 
#make a data frame with basic and full model FDR pvalue for plotting 
all_ROIs_both <- all_ROIs
all_ROIs_both$pFDR_full <-allf_ROIs$pFDR
all_ROIs_both <- all_ROIs_both %>% 
  mutate(sig.with.covars = if_else(pFDR_full <0.0166, "YES", "NO"))
#plot it!
sb <- all_ROIs_both
sb[2:8,10:12] <- as.numeric(unlist(sb[2:8,10:12]))
sb <- sb[sb$pFDR < 0.0166,]
sb$estimate <- as.numeric(sb$estimate)
sb$lci_estimate <- as.numeric(sb$lci_estimate)
sb$uci_estimate <- as.numeric(sb$uci_estimate)
plot_allb_ROIs<- ggplot() +
  geom_pointrange(data=sb, aes(x=variable, y=estimate, ymin=lci_estimate, ymax=uci_estimate,color=sig.with.covars), size=.8, shape=15)  +
  ggtitle("Brain ~ WB +sex+age+age2+assescentre+ICV+scanerXYZ+(ethn+edu+depindex+alchol+smoking+BMI)")+
  geom_hline(yintercept=0, lty=2, color ="gray20") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("standardised estimate (95% CI)") + xlab("")+
  theme_classic()+
  theme(text=element_text(size=14, family = "sans"),axis.text.y = element_text(size=14),axis.line.y = element_blank(),
        axis.title.x= element_text(size=14),legend.title = element_text(size=14),axis.text.x = element_text(size=14))+
  scale_color_manual(values=c("grey70", "black"))+
  scale_x_discrete(limits = c(
    "SA_rh_supramarginal","SA_lh_supramarginal","SA_rh_superiorparietal","SA_lh_superiorparietal","SA_rh_precuneus","SA_lh_precuneus",
    "SA_lh_precentral","SA_lh_postcentral","SA_rh_pericalcarine","Total_SA","CT_rh_superiorparietal","CT_lh_rostralanteriorcingulate",
    "CT_lh_posteriorcingulate","CT_rh_postcentral","CT_rh_pericalcarine","CT_rh_lingual","CT_rh_lateraloccipital","CT_lh_lateraloccipital",
    "CT_lh_insula","CT_rh_cuneus","CT_lh_cuneus","V_lh_precuneus","V_rh_lateraloccipital","V_rh_supramarginal","vgm_X_Cerebellum_vermis","vgm_IX_Cerebellum_right",
    "vgm_IX_Cerebellum_left","vgm_VIIIb_Cerebellum_vermis","vgm_VIIIb_Cerebellum_right","vgm_VIIIb_Cerebellum_left",
    "vgm_VIIIa_Cerebellum_right","vgm_VIIIa_Cerebellum_left","vgm_VIIb_Cerebellum_right","vgm_VIIb_Cerebellum_left",
    "vgm_Crus_II_Cerebellum_right","vgm_Crus_II_Cerebellum_left",
    "vgm_Crus_I_Cerebellum_right","vgm_Crus_I_Cerebellum_left","Volume_of_thalamus_right","Volume_of_thalamus_left","Volume_of_caudate_left",
    "Volume_of_accumbens_right","Volume_of_accumbens_left","vgm_BrainStem","Volume_of_grey_matter"),
    
    labels= c("Supramarginal right","Supramarginal left","Superior parietal right","Superior parietal left","Precuneus right","Precuneus left",
              "Precentral left","Postcentral left","Pericalcarine right","Total SA ","Superior parietal right","Rostral anterior cingulate left",
              "Posterior cingulate left","Postcentral right","Pericalcarine right","Lingual right","Lateral occipital right","Lateral occipital left",
              "Insula left","Cuneus right","Cuneus left","Precuneus left",	"Lateral occipital right",	"Supramarginal right","X Cerebellum vermis","IX Cerebellum right","IX Cerebellum left",
              "VIIIb Cerebellum vermis","VIIIb Cerebellum right","VIIIb Cerebellum left","VIIIa Cerebellum right","VIIIa Cerebellum left",
              "VIIb Cerebellum right","VIIb Cerebellum left","Crus II Cerebellum right",
              "Crus II Cerebellum left","Crus I Cerebellum right","Crus I Cerebellum left","Thalamus right","Thalamus left","Caudate left",
              "Accumbens right","Accumbens left","Brain stem","Total grey matter")
  )+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

plot_allb_ROIs










#------Plot Vol, CT and SA on brain using ggseg  ============
#make the readable data frame for atlas 
bv_ROIs_dk <- bv_ROIs[46:111,]
ctsa <- ct_ROIs
ctsa[68:134,] <- sa_ROIs
ctsa[135:200,] <-bv_ROIs_dk
ctsa <- ctsa[ctsa$pFDR <0.0166,]
ctsa <- ctsa %>%
  mutate_at("variable", str_replace, "CT_", "") %>%
  mutate_at("variable", str_replace, "SA_", "") %>%
  mutate_at("variable", str_replace, "V_", "") %>%
  rename(label=variable)
ctsa$t <- as.numeric(ctsa$t)
#change labels
grp.labs <- c("Cortical thickness", "Surface area", "Cortical volume")
names(grp.labs) <- c("CT", "SA", "vol")
#---plot using ggseg-----------------.
CT_SA_vOl_brainplot <- ctsa %>% group_by(group) %>%
  ggseg(position = "stacked",
        mapping = aes(fill = t, color=t),
        show.legend = T, color="gray20") +
  facet_wrap(~group, ncol= 3, labeller = labeller(group=grp.labs))+
  theme(axis.text = element_text(size = 15, family = "sans"),axis.title =element_text(size = 15, family = "sans"),
        strip.text = element_text(size = 17, family = "sans", face = "bold"), 
        legend.text = element_text(size = 13, family = "sans"),
        legend.title = element_text(size=15, family = "sans")
  )+
  labs(title = "")+
  scale_fill_gradient2(low ="blue",
                       high = "darkred")+
  guides(fill = guide_colourbar(barheight = 8)) 


CT_SA_vOl_brainplot
#------Plot CT and SA on brain using ggseg3d ============
#separate CT and SA
ct <- ctsa[ctsa$group=="CT",]
sa <- ctsa[ctsa$group=="SA",]
##---------------------CT---------
#ct right lateral
ct_rl<- ggseg3d(.data = ct,
                atlas = dk_3d,
                hemisphere = "right",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("right lateral") %>%
  remove_axes()
ct_rl
#ct right medial
ct_rm<- ggseg3d(.data = ct,
                atlas = dk_3d,
                hemisphere = "right",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("right medial") %>%
  remove_axes()
ct_rm

#ct left  lateral
ct_ll<- ggseg3d(.data = ct,
                atlas = dk_3d,
                hemisphere = "left",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("left lateral") %>%
  remove_axes()
ct_ll

#ct left  medial
ct_lm<- ggseg3d(.data = ct,
                atlas = dk_3d,
                hemisphere = "left",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("left medial") %>%
  remove_axes()
ct_lm
##---------------------SA-------
#sa right lateral
sa_rl<- ggseg3d(.data = sa,
                atlas = dk_3d,
                hemisphere = "right",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("right lateral") %>%
  remove_axes()
sa_rl
#sa right medial
sa_rm<- ggseg3d(.data = sa,
                atlas = dk_3d,
                hemisphere = "right",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("right medial") %>%
  remove_axes()
sa_rm

#sa left  lateral
sa_ll<- ggseg3d(.data = sa,
                atlas = dk_3d,
                hemisphere = "left",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("left lateral") %>%
  remove_axes()
sa_ll

#sa left  medial
sa_lm<- ggseg3d(.data = sa,
                atlas = dk_3d,
                hemisphere = "left",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("left medial") %>%
  remove_axes()
sa_lm


subplot(sa_rl,sa_ll,ct_rl,ct_ll, nrows = 2)

dev.off()
ct_rl

####++++++----=====A table with the sample demographics======-----++++++----
ukbbwA <-ukbbw
levels(ukbbwA$Education_2) <- c("No college degree", "College degree")
levels(ukbbwA$Sex) <-c("Female", "Male")
levels(ukbbwA$Ethnicity_W_nW) <- c("Not White","White")
levels(ukbbwA$Assessment_centre_Imaging) <- c("Bristol", "Cheadle", "Newcastle","Reading")
levels(ukbbwA$Smoking_status_img) <- c("Never", "Previous", "Current")
levels(ukbbwA$Alcohol_intake_frequency_img) <- c("Daily or almost daily","Three or four times a week","Once or twice a week",
                                                 "One to three times a month","Special occasions only","Never")
levels(ukbbwA$cat_BMI_imp) <- c("Underweight","Normal", 'Overweight', "Obese")
sumr <- ukbbwA %>% select(Age_mri,Sex,Education_2,Ethnicity_W_nW, Assessment_centre_Imaging,Smoking_status_img, 
                          Alcohol_intake_frequency_img,cat_BMI_imp,Townsend_deprivation_index_at_recruitment,wb_index_2) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 2,
    label = list(Age_mri ~ "Age", Education_2 ~ "Education" , Ethnicity_W_nW ~ "Ethnicity", 
                 Assessment_centre_Imaging ~ "Assessment centre", Smoking_status_img ~ "Smoking status",
                 Alcohol_intake_frequency_img ~"Alcohol intake frequency", cat_BMI_imp ~ "BMI (categorical)",
                 Townsend_deprivation_index_at_recruitment ~ "Townsend deprivation index", wb_index_2 ~ "Wellbeing index")) %>%
  modify_caption("Participants Characteristics (phenotype)") %>%
  bold_labels()
sumr


####################################################### PRS analysis #####################################################.=====
#----standardize brain values and mark 3SD outliers as NA & plot outliers-----
ukbbp[161:205] <- as.numeric(scale(ukbbp[161:205])) #volumes
ukbbp[327:528] <- as.numeric(scale(ukbbp[327:528])) #SA and CT and cortical volume
ukbbp[161:205][ukbbp[161:205] > 3] <- NA
ukbbp[161:205][ukbbp[161:205] < (-3)] <- NA
ukbbp[327:528][ukbbp[327:528] > 3] <- NA
ukbbp[327:528][ukbbp[327:528] < (-3)] <- NA
#----plot the number of outliers (volume)
ss<- as.data.frame(summary(is.na(ukbbp[161:205]))) 
ss$Var1 <-NULL
ss<- ss[!grepl("Mode", ss$Freq),]
ss<- ss[!grepl("FALSE", ss$Freq),]
ss$Freq <- str_remove_all(ss$Freq, "TRUE :")
colnames(ss) <- c("variable", "N_missing")
ss$N_missing <- as.numeric(ss$N_missing)
volumep_3SDw_Otliers <- ss %>%
  mutate_at("variable", str_replace, "vgm_", "") %>%
  mutate_at("variable", str_replace, "Volume_of_", "") %>%
  mutate(across("variable", str_replace_all, "_", " ")) %>%
  mutate_at("variable", str_replace, "left","lh") %>%
  mutate_at("variable", str_replace, "right","rh") %>%
  mutate_at("variable", str_replace, "grey","Total grey") %>%
  mutate_at("variable", str_replace, "white","Total white")
#plot missing and put them in a dataframe
volumep_3SDw_Otliers$variable <- as.factor(volumep_3SDw_Otliers$variable)
volumep_3SDw_Otliers$variable <- factor(volumep_3SDw_Otliers$variable, levels = volumep_3SDw_Otliers$variable[order(volumep_3SDw_Otliers$N_missing)])
Vp_3SDw_Otliers <- ggplot(volumep_3SDw_Otliers , aes(y=variable , x = N_missing )) +
  geom_bar(stat = "identity", fill=rgb(0.2,0.4,0.6), width = .6) + theme_classic() +
  scale_x_continuous(limits = c(0,255), expand = c(0, 0))+ 
  theme(panel.border = element_blank(), axis.line.y = element_blank()) + labs(x="#N +/-3SD", title = "Volume")
#plot the number of outliers (SA)
ss<- as.data.frame(summary(is.na(ukbbp[327:393]))) 
ss$Var1 <-NULL
ss<- ss[!grepl("Mode", ss$Freq),]
ss<- ss[!grepl("FALSE", ss$Freq),]
ss$Freq <- str_remove_all(ss$Freq, "TRUE :")
colnames(ss) <- c("variable", "N_missing")
ss$N_missing <- as.numeric(ss$N_missing)
SAp_3SDw_Otliers <- ss %>%
  mutate_at("variable", str_replace, "SA_", "") %>%
  mutate(across("variable", str_replace_all, "_", " "))
#plot missing and put them in a dataframe
SAp_3SDw_Otliers$variable <- as.factor(SAp_3SDw_Otliers$variable)
SAp_3SDw_Otliers$variable <- factor(SAp_3SDw_Otliers$variable, levels = SAp_3SDw_Otliers$variable[order(SAp_3SDw_Otliers$N_missing)])
Sp_3SDw_Otliers<- ggplot(SAp_3SDw_Otliers , aes(y=variable , x = N_missing )) +
  geom_bar(stat = "identity", fill=rgb(0.2,0.4,0.6), width = .6) + theme_classic() +
  scale_x_continuous(limits = c(0,450), expand = c(0, 0))+ 
  theme(panel.border = element_blank(), axis.title.y = element_blank(),axis.line.y = element_blank())+
  labs(x="#N +/-3SD", title = "Surface area")
#plot the number of outliers (CT)
ss<- as.data.frame(summary(is.na(ukbbp[396:462]))) 
ss$Var1 <-NULL
ss<- ss[!grepl("Mode", ss$Freq),]
ss<- ss[!grepl("FALSE", ss$Freq),]
ss$Freq <- str_remove_all(ss$Freq, "TRUE :")
colnames(ss) <- c("variable", "N_missing")
ss$N_missing <- as.numeric(ss$N_missing)
CTp_3SDw_Otliers <- ss %>%
  mutate_at("variable", str_replace, "CT_", "") %>%
  mutate(across("variable", str_replace_all, "_", " "))
#plot missing and put them in a dataframe
CTp_3SDw_Otliers$variable <- as.factor(CTp_3SDw_Otliers$variable)
CTp_3SDw_Otliers$variable <- factor(CTp_3SDw_Otliers$variable, levels = CTp_3SDw_Otliers$variable[order(CTp_3SDw_Otliers$N_missing)])
Cp_3SDw_Otliers<- ggplot(CTp_3SDw_Otliers , aes(y=variable , x = N_missing )) +
  geom_bar(stat = "identity", fill=rgb(0.2,0.4,0.6), width = .6) + theme_classic() +
  scale_x_continuous(limits = c(0,450), expand = c(0, 0))+ 
  theme(panel.border = element_blank(), axis.title.y = element_blank(),axis.line.y = element_blank())+
  labs(x="#N +/-3SD", title = "Cortical thickness")
#plot the number of outliers (VOL_DK)
ss<- as.data.frame(summary(is.na(ukbbp[463:528]))) 
ss$Var1 <-NULL
ss<- ss[!grepl("Mode", ss$Freq),]
ss<- ss[!grepl("FALSE", ss$Freq),]
ss$Freq <- str_remove_all(ss$Freq, "TRUE :")
colnames(ss) <- c("variable", "N_missing")
ss$N_missing <- as.numeric(ss$N_missing)
Vop_3SDw_Otliers <- ss %>%
  mutate_at("variable", str_replace, "V_", "") %>%
  mutate(across("variable", str_replace_all, "_", " "))
#plot missing and put them in a dataframe
Vop_3SDw_Otliers$variable <- as.factor(Vop_3SDw_Otliers$variable)
Vop_3SDw_Otliers$variable <- factor(Vop_3SDw_Otliers$variable, levels = Vop_3SDw_Otliers$variable[order(Vop_3SDw_Otliers$N_missing)])
Vop_3SDw_Otliers<-ggplot(Vop_3SDw_Otliers , aes(y=variable , x = N_missing )) +
  geom_bar(stat = "identity", fill=rgb(0.2,0.4,0.6), width = .6) + theme_classic() +
  scale_x_continuous(limits = c(0,450), expand = c(0, 0))+ 
  theme(panel.border = element_blank(), axis.title.y = element_blank(),axis.line.y = element_blank())+
  labs(x="#N +/-3SD", title = "Cortical volume")
#put all plots together
ggarrange(Vp_3SDw_Otliers,Vop_3SDw_Otliers,Cp_3SDw_Otliers,Sp_3SDw_Otliers,
          ncol = 4, nrow = 1, widths = c(1.4,2,2,2))
#----restandardise values------
ukbbp$PRS_Brain <- as.numeric(scale(ukbbp$PRS_Brain))
ukbbp[162:205] <- as.numeric(scale(ukbbp[162:205]))
ukbbp[327:528] <- as.numeric(scale(ukbbp[327:528]))
#Quick check of normality===== 
par(mfrow= c(5,6))
for (i in c(161:205)){
  hist(ukbbp[,i], xlab = colnames(ukbbp[i]), breaks= 40, main = "", ylab = "")
} 
par(mfrow= c(5,6))
for (i in c(327:528)){
  hist(ukbbp[,i], xlab = colnames(ukbbp[i]), breaks= 40, main = "", ylab = "")
}
###--------------------------Linear models for basic covariates------------------
##linear models for  Brain VOLUMES DK,FIRST and FAST-------
bvp_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(bvp_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(161:205,463:528)) {
  fitm <-lm.beta(lm( ukbbp[,i]~ PRS_Brain+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+ 
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ICV+
                       genotype_array+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10
                     , ukbbp)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbp[,i])) # Sample size of the test
  bvp_ROIs[nrow(bvp_ROIs)+1,] <- c(colnames(ukbbp)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"vol" ) #save in a dataframe
}
bvp_ROIs$pFDR <- p.adjust(bvp_ROIs$pval, method = "BH", n = nrow(bvp_ROIs))
bvp_ROIs$pbonf <- p.adjust(bvp_ROIs$pval, method = "bonferroni", 245)

##linear models for CT-------
ctp_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(ctp_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval", "group")
#loop for models
for (i in c(396:462)) {
  fitm <-lm.beta(lm( ukbbp[,i]~ PRS_Brain+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+ 
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ICV+
                       genotype_array+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10
                     , ukbbp)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbp[,i])) # Sample size of the test
  ctp_ROIs[nrow(ctp_ROIs)+1,] <- c(colnames(ukbbp)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"CT" ) #save in a dataframe
}
ctp_ROIs$pFDR <- p.adjust(ctp_ROIs$pval, method = "BH", n = nrow(ctp_ROIs))
ctp_ROIs$pbonf <- p.adjust(ctp_ROIs$pval, method = "bonferroni", 245)

##linear models for  SA-------
sap_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(sap_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(327:393)) {
  fitm <-lm.beta(lm( ukbbp[,i]~ PRS_Brain+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+ 
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ ICV+
                       genotype_array+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10
                     , ukbbp)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbp[,i])) # Sample size of the test
  sap_ROIs[nrow(sap_ROIs)+1,] <- c(colnames(ukbbp)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"SA" ) #save in a dataframe
}
sap_ROIs$pFDR <- p.adjust(sap_ROIs$pval, method = "BH", n = nrow(sap_ROIs))
sap_ROIs$pbonf <- p.adjust(sap_ROIs$pval, method = "bonferroni", 245)

###--------------------------Linear models with full covariates------------------
##linear models for  Brain VOLUMES-------
bvpf_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(bvpf_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(161:205,463:528)) {
  fitm <-lm.beta(lm( ukbbp[,i]~ PRS_Brain+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+cat_BMI_imp+
                       Education_2+Townsend_deprivation_index_at_recruitment+Smoking_status_img+Alcohol_intake_frequency_img+
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ ICV+
                       genotype_array+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10
                     , ukbbp)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbp[,i])) # Sample size of the test
  bvpf_ROIs[nrow(bvpf_ROIs)+1,] <- c(colnames(ukbbp)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"vol" ) #save in a dataframe
}
bvpf_ROIs$pFDR <- p.adjust(bvpf_ROIs$pval, method = "BH", n = nrow(bvpf_ROIs))
bvpf_ROIs$pbonf <- p.adjust(bvpf_ROIs$pval, method = "bonferroni", 245)

##linear models for CT-------
ctpf_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(ctpf_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(396:462)) {
  fitm <-lm.beta(lm( ukbbp[,i]~ PRS_Brain+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+cat_BMI_imp+
                       Education_2+Townsend_deprivation_index_at_recruitment+Smoking_status_img+Alcohol_intake_frequency_img+
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ ICV+
                       genotype_array+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10
                     , ukbbp)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbp[,i])) # Sample size of the test
  ctpf_ROIs[nrow(ctpf_ROIs)+1,] <- c(colnames(ukbbp)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"CT" ) #save in a dataframe
}
ctpf_ROIs$pFDR <- p.adjust(ctpf_ROIs$pval, method = "BH", n = nrow(ctpf_ROIs))
ctpf_ROIs$pbonf <- p.adjust(ctpf_ROIs$pval, method = "bonferroni", 245)

##linear models for  SA-------
safp_ROIs <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(safp_ROIs) <- c("variable", "estimate", "lci_estimate", "uci_estimate", "stderr_estimate","t","N","pval","group")
#loop for models
for (i in c(327:393)) {
  fitm <-lm.beta(lm( ukbbp[,i]~ PRS_Brain+Sex+Age_mri+Age_mri_squared+Assessment_centre_Imaging+cat_BMI_imp+
                       Education_2+Townsend_deprivation_index_at_recruitment+Smoking_status_img+Alcohol_intake_frequency_img+
                       Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position+ ICV+
                       genotype_array+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10
                     , ukbbp)) #full model
  fit <- summary(fitm)
  p_val <- fit$coefficients[2,5] #p-value
  t <- fit$coefficients[2,4] #t-value
  estimate <- fit$coefficients[2,2] # standardized coefficient of the regression or estimate
  lci_estimate <- confint(fitm)[2,1]# 95% confidence interval of the estimate lower band
  uci_estimate <- confint(fitm)[2,2]# 95% confidence interval of the estimate upper band
  std_error <- fit$coefficients[2,3]# standard error of estimate
  n <- length(na.omit(ukbbp[,i])) # Sample size of the test
  safp_ROIs[nrow(safp_ROIs)+1,] <- c(colnames(ukbbp)[i],estimate,lci_estimate,uci_estimate,std_error,t,n,p_val,"SA" ) #save in a dataframe
}
safp_ROIs$pFDR <- p.adjust(safp_ROIs$pval, method = "BH", n = nrow(safp_ROIs))
safp_ROIs$pbonf <- p.adjust(safp_ROIs$pval, method = "bonferroni", 245)

#Put all analysis in a single dataframe-------------
#basic models
allp_ROIs <- bvp_ROIs
allp_ROIs[112:178,] <- ctp_ROIs
allp_ROIs[179:245,] <- sap_ROIs
#with all covariates
allpf_ROIs <- bvpf_ROIs
allpf_ROIs[112:178,] <- ctpf_ROIs
allpf_ROIs[179:245,] <- safp_ROIs
###Plotting time!=======================================================
#BASIC models; all regions that survive FDR-P-value of 0.0166 (0.05/3)=======
s <- allp_ROIs
s[2:8,10:11] <- as.numeric(unlist(s[2:8,10:11]))
s <- s[s$pFDR < 0.0166,]
s$estimate <- as.numeric(s$estimate)
s$lci_estimate <- as.numeric(s$lci_estimate)
s$uci_estimate <- as.numeric(s$uci_estimate)
plot_allp_ROIs<- ggplot() +
  geom_pointrange(data=s, aes(x=variable, y=estimate, ymin=lci_estimate, ymax=uci_estimate,color=group),shape =15, size=.8)  +
  ggtitle("Brain ~ WB_PRS +sex+age+age2+assescentre+ICV+scanerXYZ")+
  geom_hline(yintercept=0, lty=2, color ="gray20") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("standardised estimate (95% CI)") + xlab("")+
  theme_classic()+
  theme(text =element_text(size=14),axis.text = element_text(size=14),plot.title = element_text(color="gray20", size=12, face="italic"))
plot_allp_ROIs
#FULL models; all regions that survive FDR-P-value of 0.0166 (0.05/3)=======
sf <- allpf_ROIs
sf[2:8,10:11] <- as.numeric(unlist(sf[2:8,10:11]))
sf <- sf[sf$pFDR < 0.0166,]
sf$estimate <- as.numeric(sf$estimate)
sf$lci_estimate <- as.numeric(sf$lci_estimate)
sf$uci_estimate <- as.numeric(sf$uci_estimate)
plot_allpf_ROIs<- ggplot() +
  geom_pointrange(data=sf, aes(x=variable, y=estimate, ymin=lci_estimate, ymax=uci_estimate,color=group),shape =15, size=.8)  +
  ggtitle("Brain ~ WB_PRS +sex+age+age2+assescentre+ICV+scanerXYZ+ethn+edu+depindex+alchol+smoking+BMI")+
  geom_hline(yintercept=0, lty=2, color ="gray20") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("standardised estimate (95% CI)") + xlab("")+
  theme_classic()+
  theme(text =element_text(size=14),axis.text = element_text(size=14),plot.title = element_text(color="gray20", size=12, face="italic"))
plot_allpf_ROIs

#A plot with shows which variables still significant after adding all covariates===== 
#make a data frame with basic and full model FDR pvalue for plotting 
allp_ROIs_both <- allp_ROIs
allp_ROIs_both$pFDR_full <-allpf_ROIs$pFDR
allp_ROIs_both <- allp_ROIs_both %>% 
  mutate(sig.with.covars = if_else(pFDR_full <0.0166, "YES", "NO"))
#plot it!
sb <- allp_ROIs_both
sb[2:8,10:12] <- as.numeric(unlist(sb[2:8,10:12]))
sb <- sb[sb$pFDR < 0.0166,]
sb$estimate <- as.numeric(sb$estimate)
sb$lci_estimate <- as.numeric(sb$lci_estimate)
sb$uci_estimate <- as.numeric(sb$uci_estimate)
plot_allb_ROIs<- ggplot() +
  geom_pointrange(data=sb, aes(x=variable, y=estimate, ymin=lci_estimate, ymax=uci_estimate,color=sig.with.covars,shape=group), size=.8)  +
  ggtitle("Brain ~ WB_PRS +sex+age+age2+assescentre+ICV+scanerXYZ+(ethn+edu+depindex+alchol+smoking+BMI)")+
  geom_hline(yintercept=0, lty=2, color ="gray20") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("standardised estimate (95% CI)") + xlab("")+
  theme_classic()+
  theme(text =element_text(size=14),axis.text = element_text(size=14),plot.title = element_text(color="gray20", size=12, face="italic"))+
  scale_color_manual(values=c("grey70", "black"))+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

plot_allb_ROIs


#------Plot CT and SA on brain using ggseg  ============
#make the readable data frame for atlas 
ctsa <- ctp_ROIs
ctsa[68:134,] <- sap_ROIs
ctsa <- ctsa[ctsa$pFDR <0.0166,]
ctsa <- ctsa %>%
  mutate_at("variable", str_replace, "CT_", "") %>%
  mutate_at("variable", str_replace, "SA_", "") %>%
  rename(label=variable)
ctsa$t <- as.numeric(ctsa$t)
#---plot using ggseg-----------------.
CT_SA_brainplot <- ctsa %>% group_by(group) %>%
  ggseg(position = "stacked",
        mapping = aes(fill = t, color=t),
        show.legend = T, color="black") +
  facet_wrap(~group, ncol= 2)+
  labs(title = "")+
  scale_fill_gradient2(low ="darkred",
                       high = "#132B45")+
  guides(fill = guide_colourbar(barheight = 10))  

CT_SA_brainplot
#------Plot CT and SA on brain using ggseg3d ============
#separate CT and SA
ct <- ctsa[ctsa$group=="CT",]
sa <- ctsa[ctsa$group=="SA",]
##---------------------CT---------
#ct right lateral
ct_rl<- ggseg3d(.data = ct,
                atlas = dk_3d,
                hemisphere = "right",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("right lateral") %>%
  remove_axes()
ct_rl
#ct right medial
ct_rm<- ggseg3d(.data = ct,
                atlas = dk_3d,
                hemisphere = "right",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("right medial") %>%
  remove_axes()
ct_rm

#ct left  lateral
ct_ll<- ggseg3d(.data = ct,
                atlas = dk_3d,
                hemisphere = "left",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("left lateral") %>%
  remove_axes()
ct_ll

#ct left  medial
ct_lm<- ggseg3d(.data = ct,
                atlas = dk_3d,
                hemisphere = "left",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("left medial") %>%
  remove_axes()
ct_lm
##---------------------SA-------
#sa right lateral
sa_rl<- ggseg3d(.data = sa,
                atlas = dk_3d,
                hemisphere = "right",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("right lateral") %>%
  remove_axes()
sa_rl
#sa right medial
sa_rm<- ggseg3d(.data = sa,
                atlas = dk_3d,
                hemisphere = "right",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("right medial") %>%
  remove_axes()
sa_rm

#sa left  lateral
sa_ll<- ggseg3d(.data = sa,
                atlas = dk_3d,
                hemisphere = "left",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("left lateral") %>%
  remove_axes()
sa_ll

#sa left  medial
sa_lm<- ggseg3d(.data = sa,
                atlas = dk_3d,
                hemisphere = "left",
                label = "label",
                colour = "t", text = "t",
                palette = c("darkred" = -6, "white" =
                              0, "#132B45" = 5.2)
) %>%
  pan_camera("left medial") %>%
  remove_axes()
sa_lm


subplot(sa_rl,sa_ll,ct_rl,ct_ll, nrows = 2)

dev.off()
ct_rl

####++++++----=====A table with the sample demographics======-----++++++----
ukbbpA <- ukbbp
levels(ukbbpA$Education_2) <- c("No college degree", "College degree")
levels(ukbbpA$Sex) <-c("Female", "Male")
levels(ukbbpA$Ethnicity_W_nW) <- c("Not White","White")
levels(ukbbpA$Assessment_centre_Imaging) <- c("Bristol", "Cheadle", "Newcastle","Reading")
levels(ukbbpA$Smoking_status_img) <- c("Never", "Previous", "Current")
levels(ukbbpA$Alcohol_intake_frequency_img) <- c("Daily or almost daily","Three or four times a week","Once or twice a week",
                                                 "One to three times a month","Special occasions only","Never")
levels(ukbbpA$cat_BMI_imp) <- c("Underweight","Normal", 'Overweight', "Obese")
sumrp <- ukbbpA %>% select(Age_mri,Sex,Education_2,Ethnicity_W_nW, Assessment_centre_Imaging,Smoking_status_img, 
                           Alcohol_intake_frequency_img,cat_BMI_imp,Townsend_deprivation_index_at_recruitment,wb_index_2) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 2,
    label = list(Age_mri ~ "Age", Education_2 ~ "Education" , Ethnicity_W_nW ~ "Ethnicity", 
                 Assessment_centre_Imaging ~ "Assessment centre", Smoking_status_img ~ "Smoking status",
                 Alcohol_intake_frequency_img ~"Alcohol intake frequency", cat_BMI_imp ~ "BMI (categorical)",
                 Townsend_deprivation_index_at_recruitment ~ "Townsend deprivation index", wb_index_2 ~ "Wellbeing index")) %>%
  modify_caption("Participants Characteristics (PRS)") %>%
  bold_labels()
sumrp







####################################################### Joint Wellbeing index and PRSanalysis #############################.=====
#------ ploting PRS and wellbeing significant regions on the same plot-----
#merge pheno and prs
all_ROIs_both$type <- "Wellbeing index"
allp_ROIs_both$type <- "Wellbeing index PGS"
all_phen_prs <- all_ROIs_both
all_phen_prs[246:490,] <- allp_ROIs_both

sb <- all_phen_prs
sb[2:8,10:12] <- as.numeric(unlist(sb[2:8,10:12]))
sb <- sb[sb$pFDR < 0.0166,]
sb$estimate <- as.numeric(sb$estimate)
sb$lci_estimate <- as.numeric(sb$lci_estimate)
sb$uci_estimate <- as.numeric(sb$uci_estimate)
#plot using ggplot
plot_allb_ROIs<- ggplot() +
  geom_pointrange(data=sb, aes(x=variable, y=estimate, ymin=lci_estimate, ymax=uci_estimate,color=sig.with.covars,shape=group), size=1)  +
  ggtitle("")+
  geom_hline(yintercept=0, lty=2, color ="gray20") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  facet_grid(cols = vars(type))+
  ylab("Standardised estimate (95% CI)") + xlab("")+
  theme_bw()+
  theme(text=element_text(size=14, family = "sans"),axis.text.y = element_text(size=13),axis.line.y = element_blank(),
        axis.title.x= element_text(size=12),legend.title = element_text(size=12))+
  scale_color_manual(values=c("grey70", "black"))+
  scale_x_discrete(limits = c(
    "SA_rh_supramarginal","SA_lh_supramarginal","SA_rh_superiorparietal","SA_lh_superiorparietal","SA_rh_precuneus","SA_lh_precuneus",
    "SA_lh_precentral","SA_lh_postcentral","SA_rh_pericalcarine","Total_SA","CT_rh_superiorparietal","CT_lh_rostralanteriorcingulate",
    "CT_lh_posteriorcingulate","CT_rh_postcentral","CT_rh_pericalcarine","CT_rh_lingual","CT_rh_lateraloccipital","CT_lh_lateraloccipital",
    "CT_lh_insula","CT_rh_cuneus","CT_lh_cuneus","V_lh_precuneus","V_rh_lateraloccipital","V_rh_supramarginal","vgm_X_Cerebellum_vermis","vgm_IX_Cerebellum_right",
    "vgm_IX_Cerebellum_left","vgm_VIIIb_Cerebellum_vermis","vgm_VIIIb_Cerebellum_right","vgm_VIIIb_Cerebellum_left",
    "vgm_VIIIa_Cerebellum_right","vgm_VIIIa_Cerebellum_left","vgm_VIIb_Cerebellum_right","vgm_VIIb_Cerebellum_left",
    "vgm_Crus_II_Cerebellum_right","vgm_Crus_II_Cerebellum_left",
    "vgm_Crus_I_Cerebellum_right","vgm_Crus_I_Cerebellum_left","Volume_of_thalamus_right","Volume_of_thalamus_left","Volume_of_caudate_left",
    "Volume_of_accumbens_right","Volume_of_accumbens_left","vgm_BrainStem","Volume_of_grey_matter"),
    
    labels= c("Supramarginal right","Supramarginal left","Superior parietal right","Superior parietal left","Precuneus right","Precuneus left",
              "Precentral left","Postcentral left","Pericalcarine right","Total SA ","Superior parietal right","Rostral anterior cingulate left",
              "Posterior cingulate left","Postcentral right","Pericalcarine right","Lingual right","Lateral occipital right","Lateral occipital left",
              "Insula left","Cuneus right","Cuneus left","Precuneus left",	"Lateral occipital right",	"Supramarginal right","X Cerebellum vermis","IX Cerebellum right","IX Cerebellum left",
              "VIIIb Cerebellum vermis","VIIIb Cerebellum right","VIIIb Cerebellum left","VIIIa Cerebellum right","VIIIa Cerebellum left",
              "VIIb Cerebellum right","VIIb Cerebellum left","Crus II Cerebellum right",
              "Crus II Cerebellum left","Crus I Cerebellum right","Crus I Cerebellum left","Thalamus right","Thalamus left","Caudate left",
              "Accumbens right","Accumbens left","Brain stem","Total grey matter")
  )
plot_allb_ROIs
#------ wellbeing and PRS samples demographics together!========
#joining the two dataframes 
ukbbwA$type <- "phenotype"
ukbbpA$type <- "PRS"
ukbbwpA <- ukbbwA
ukbbwpA[38984:58970, ] <- ukbbpA
ukbbwpA <- ukbbwpA[complete.cases(ukbbwpA$type),]
summary(as.factor(ukbbwpA$type))
##the demographic table
sumrpw <- ukbbwpA %>% select(Age_mri,Sex,Education_2,Ethnicity_W_nW, Assessment_centre_Imaging,Smoking_status_img, 
                             Alcohol_intake_frequency_img,cat_BMI_imp,Townsend_deprivation_index_at_recruitment,wb_index_2,type) %>%
  tbl_summary(by = type,
              statistic = list(all_continuous() ~ "{mean} ({sd})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 2,
              label = list(Age_mri ~ "Age", Education_2 ~ "Education" , Ethnicity_W_nW ~ "Ethnicity", 
                           Assessment_centre_Imaging ~ "Assessment centre", Smoking_status_img ~ "Smoking status",
                           Alcohol_intake_frequency_img ~"Alcohol intake frequency", cat_BMI_imp ~ "BMI (categorical)",
                           Townsend_deprivation_index_at_recruitment ~ "Townsend deprivation index", wb_index_2 ~ "Wellbeing index")) %>%
  modify_caption("Participants Characteristics") %>%
  bold_labels()
sumrpw

#------ Plots for comparing wellbeing index in covariates ---------------
#create appropriate dataset for the plot
covw <- ukbbwA %>% select(Sex,Education_2,Ethnicity_W_nW, Assessment_centre_Imaging,Smoking_status_img, 
                          Alcohol_intake_frequency_img,cat_BMI_imp,wb_index_2) %>%
  rename("Education"=Education_2 ,"Ethnicity"=Ethnicity_W_nW, 
         "Assessment centre"= Assessment_centre_Imaging,"Smoking status"=Smoking_status_img ,
         "Alcohol intake frequency"=Alcohol_intake_frequency_img , "BMI (categorical)"=cat_BMI_imp 
         ,"Wellbeing index"= wb_index_2)%>%
  gather(variable, measurement, Sex:"BMI (categorical)", factor_key = T) 
#order levels
covw$measurement <- factor(covw$measurement , levels=c("Bristol","Cheadle","Newcastle","Reading","College degree",
                                                       "No college degree","Never","Previous","Current","Female","Male","Special occasions only",
                                                       "One to three times a month","Once or twice a week","Three or four times a week","Daily or almost daily",
                                                       "Underweight","Normal","Overweight","Obese","White","Not White")) 
#make vector for statistical comparisions
co_temp <- ukbbwA %>% select(Sex,Education_2,Ethnicity_W_nW, Assessment_centre_Imaging,Smoking_status_img, 
                             Alcohol_intake_frequency_img,cat_BMI_imp,wb_index_2) %>%
  gather(variable, measurement, Sex:cat_BMI_imp, factor_key = T)
#run tests
covw_copm <- compare_means( wb_index_2 ~ measurement ,co_temp, group.by = "variable" ,  method = "kruskal.test") 
#rename variables accordingly
covw_copm$variable <- factor(c("Sex","Education","Ethnicity","Assessment centre","Smoking status","Alcohol intake frequency","BMI (categorical)"),
                             levels = c("Sex","Education","Ethnicity","Assessment centre","Smoking status","Alcohol intake frequency","BMI (categorical)"))
covw_copm<- as.data.frame(covw_copm)

#plot it 
ggplot(covw, aes(x=`Wellbeing index`,y=measurement, fill=variable, color=variable)) + 
  geom_violin(trim=FALSE)+
  scale_fill_viridis_d(alpha = .2)+
  scale_colour_viridis_d(alpha = .5)+
  geom_boxplot(width=.4, fill="white", color="black", outlier.shape = NA, lwd=.8)+
  theme_bw()+
  geom_vline(xintercept = median(covw$`Wellbeing index`), linetype = 2, color="grey40" )+
  theme(text = element_text(size = 15))+
  facet_col(vars(variable), scales = "free_y", space = "free", strip.position = "top")+
  theme(axis.text.y = element_text( hjust = 1, vjust = .5, size = 13), axis.title.y = element_blank(),axis.title.x = element_text(size=13, hjust = .53), 
        legend.position = "none", strip.text = element_text(hjust = 0, face = "italic"))+
  geom_text(x=3.6,y=c(2.5,2.5,2.5,4,3,6,4), aes(label =paste("P =",formatC(covw_copm$p,digits = 2))), 
            data = covw_copm, color="black", angle = 270, nudge_y = 2, hjust=0)+
  scale_x_continuous(breaks = c(-2,0,2),labels = c(-2,0,2))




write.csv(bv_ROIs, "bv.csv")
write.csv(bvf_ROIs , "bvf.csv")
write.csv(ctp_ROIs , "ctp.csv")
write.csv(ctpf_ROIs , "ctpf.csv")
write.csv(sap_ROIs , "sap.csv")
write.csv(safp_ROIs , "sapf.csv")
