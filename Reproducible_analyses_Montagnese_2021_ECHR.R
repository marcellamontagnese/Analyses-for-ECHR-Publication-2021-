#### Analyses for paper by Montagnese et al., ####
## "Cognition, hallucination severity and hallucination-specific insight in neurodegenerative disorders and eye disease"

####################### 
#### Prerequisites ####
#install.packages('lme4','haven','tidyverse','RColorBrewer','lmerTest','arm','ggplot2')
library(lme4) # for the analysis
library(haven) # to load the SPSS .sav file
library(tidyverse) # needed for data manipulation.
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(lmerTest)# to get p-value estimations that are not part of the standard lme4 packages
library(arm)
library(ggplot2)
library(emmeans)
library(dplyr)
library(rlang)
library(sjPlot)

#### Dataset preparation ####
dataset <- read.csv('/Users/marcellamontagnese/Google Drive/PhD KCL/Publications in progress/Prabitha ECHR/dataset_4.csv')
dataset_4 <- dataset
#### Scale explanatory variables of interest 
dataset_4$Fright_x_freq <- scale(dataset_4$Fright_x_freq, center = TRUE, scale = TRUE) # scaling explanatory variables
dataset_4$Believe_real_combined_MM_updated <- scale(dataset_4$Believe_real_combined_MM_updated, center = TRUE, scale = TRUE) # scaling explanatory variables
dataset_4$PNEVHI_How_often <- scale(dataset_4$PNEVHI_How_often, center = TRUE, scale = TRUE) # scaling explanatory variables
dataset_4$PNEVHI_How_long <- scale(dataset_4$PNEVHI_How_long, center = TRUE, scale = TRUE) # scaling explanatory variables
dataset_4$PNEVHI_Frightening <- scale(dataset_4$PNEVHI_Frightening, center = TRUE, scale = TRUE) # scaling explanatory variables
dataset_4$Age <- scale(dataset_4$Age, center = TRUE, scale = TRUE) # scaling explanatory variable
dataset_4$MMSE_tot <- scale(dataset_4$MMSE_tot, center = TRUE, scale = TRUE) # scaling explanatory variable
#### Set constrasts
dataset_4$Diagnose_Prabs <- relevel(dataset_4$Diagnose_Prabs, ref = "EyeDisease")
contrasts(dataset_4$Diagnose_Prabs) # to get contrasts. When you have categorical variables as predictors, R uses one of the levels as the reference.


#### Check Null Models ####
Null<-lmer(Believe_real_combined_MM_updated ~ 1 # Degree of insight as predictor
           +(1|Diagnose_Prabs), # each group gets its own intercept 
           data=dataset, REML = FALSE)
summary(Null)

Null<-lmer(MMSE_tot ~ 1 # Degree of insight as predictor
           +(1|Diagnose_Prabs), # each group gets its own intercept 
           data=dataset, REML = FALSE)
summary(Null)

# Check the intra-class correlation (ICC) to determine if multi-level modeling is the correct choice for our analysis. 
# If my ICC is greater than 0, it is a multi-level study.
ICC.Model<-function(Model.Name) {
  tau.Null<-as.numeric(lapply(summary(Model.Name)$varcor, diag))
  sigma.Null <- as.numeric(attr(summary(Model.Name)$varcor, "sc")^2)
  ICC.Null <- tau.Null/(tau.Null+sigma.Null)
  return(ICC.Null)
}
ICC.Model(Null) # --> it is a multilevel study as the ICC is greater than zero! (0.62)


#### Linear mixed models ####
## Testing Hypothesis 1 
## Model 1 and 2 look at degree of insight as a predictor ## 
## Model 1 - severity as one measure
mixed.lmer_random_1_ML<- lmer(Believe_real_combined_MM_updated ~ Fright_x_freq * MMSE_tot +  PNEVHI_How_long + Age
                              + (1|Diagnose_Prabs), data = dataset_4, REML = FALSE)
summary(mixed.lmer_random_1_ML)
mixed.lmer_random_table_1_ML <- tab_model(mixed.lmer_random_1_ML) # to get HTML summary table
mixed.lmer_random_table_1_ML
## Model 2 - severity split into its component measures of frequency and distress
mixed.lmer_random_2_ML<- lmer(Believe_real_combined_MM_updated ~ PNEVHI_How_often + PNEVHI_Frightening + MMSE_tot + PNEVHI_How_long + Age
                              + (1|Diagnose_Prabs), data = dataset_4, REML = FALSE)
summary(mixed.lmer_random_2_ML)
mixed.lmer_random_table_2_ML <- tab_model(mixed.lmer_random_2_ML) # to get HTML summary table
mixed.lmer_random_table_2_ML

## Testing Hypothesis 2
## Model 3 and 4 look at total MMSE score as a predictor
# firstly with severity + duration
mixed.lmer_random_1_ML_MMSE<- lmer(MMSE_tot ~ Fright_x_freq * Believe_real_combined_MM_updated +  PNEVHI_How_long + Age
                                   + (1|Diagnose_Prabs), data = dataset_4, REML = FALSE)
summary(mixed.lmer_random_1_ML_MMSE)
mixed.lmer_random_table_1_ML_MMSE <- tab_model(mixed.lmer_random_1_ML_MMSE)# to get HTML summary table
mixed.lmer_random_table_1_ML_MMSE
# then again but this time with severity split into its components of distress and frequency + the usual duration aside
mixed.lmer_random_2_ML_MMSE<- lmer(MMSE_tot ~  Believe_real_combined_MM_updated +  PNEVHI_How_often + PNEVHI_Frightening + PNEVHI_How_long + Age
                                   + (1|Diagnose_Prabs), data = dataset_4, REML = FALSE)
summary(mixed.lmer_random_2_ML_MMSE)
mixed.lmer_random_table_2_ML_MMSE <- tab_model(mixed.lmer_random_2_ML_MMSE)# to get HTML summary table
mixed.lmer_random_table_2_ML_MMSE


#### Spearman Correlation analyses with Bonferroni correction ####
library(ggplot2)
library(reshape2)
library(Hmisc)
library(stats)

## Subset dataset by group
DLB <- dataset_4%>% filter(Diagnose_Prabs == "DLB")
PD <- dataset_4%>% filter(Diagnose_Prabs == "PD")
PDD <- dataset_4%>% filter(Diagnose_Prabs == "PDD")
EY <- dataset_4%>% filter(Diagnose_Prabs == "EyeDisease")

## Open this function before starting - as it's called later on
abbreviateSTR <- function(value, prefix){  # format string more concisely
  lst = c()
  for (item in value) {
    if (is.nan(item) || is.na(item)) { # if item is NaN return empty string
      lst <- c(lst, '')
      next
    }
    item <- round(item, 2) # round to two digits
    if (item == 0) { # if rounding results in 0 clarify
      item = '<.001'
    }
    item <- as.character(item)
    item <- sub("(^[0])+", "", item)    # remove leading 0: 0.05 -> .05
    item <- sub("(^-[0])+", "-", item)  # remove leading -0: -0.05 -> -.05
    lst <- c(lst, paste(prefix, item, sep = ""))
  }
  return(lst)
}

#### For the **DLB** group
dataset_corr <- DLB[, c("Fright_x_freq", "Believe_real_combined_MM_updated", "MMSE_tot", "PNEVHI_How_long", "PNEVHI_How_often","PNEVHI_Frightening")] # "categorical_fluency_Animal"
d <- dataset_corr
cormatrix = rcorr(as.matrix(d[,c(1,2,3,4,5,6)]), type='spearman') # Type of correlation
cormatrix$adj = p.adjust(cormatrix$P, method = "bonferroni") # Which multiple comparison correction to apply
cordata = melt(cormatrix$r) 
cordata$labelr = abbreviateSTR(melt(cormatrix$r)$value, 'r')
cordata$labelP = abbreviateSTR(melt(cormatrix$P)$value, 'P non-adjusted')
cordata$labelPadj = abbreviateSTR(melt(cormatrix$adj)$value, 'P')
cordata$label = paste(cordata$labelr, "\n", 
                      cordata$labelPadj, sep = "")
cordata$strike = ""
cordata$strike[cormatrix$adj > 0.05] = "X" # changed this to only show the significant ones
txtsize <- par('din')[2] / 2
# Plot
DLB_corr <- ggplot(cordata, aes(x=Var1, y=Var2, fill=value), type = "upper") + geom_tile() + 
  #theme(axis.text.x = element_text(angle=90, hjust=TRUE, size=16)) +
  xlab("") + ylab("") + labs(title="* DLB group") +
  geom_text(label=cordata$label, size=txtsize*1.3, color="black") + 
  geom_text(label=cordata$strike, size=txtsize*3 , color="black", alpha=0.2) +
  scale_fill_gradient2(low = "#00688B", high = "#EE5C42", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 13, hjust = 1))+
  theme(axis.text.y = element_text(angle = 0, vjust = 1, 
                                   size = 13, hjust = 1)) +
  scale_x_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot', 'Hallucinations  \nduration', 'Hallucinations  \nfrequency','Hallucinations  \ndistress')) + # 'Categorical fluency score'
  scale_y_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot', 'Hallucinations  \nduration', 'Hallucinations  \nfrequency','Hallucinations  \ndistress'))
DLB_corr


### PD
dataset_corr <- PD[, c("Fright_x_freq", "Believe_real_combined_MM_updated", "MMSE_tot","edu_year_Prabs","PNEVHI_How_long", "PNEVHI_How_often","PNEVHI_Frightening")] # "categorical_fluency_Animal"
d <- dataset_corr
cormatrix = rcorr(as.matrix(d[,c(1,2,3,4,5,6,7)]), type='spearman')
cormatrix$adj = p.adjust(cormatrix$P, method = "bonferroni") 
cordata = melt(cormatrix$r) 

cordata$labelr = abbreviateSTR(melt(cormatrix$r)$value, 'r')
cordata$labelP = abbreviateSTR(melt(cormatrix$P)$value, 'P non-adjusted')
cordata$labelPadj = abbreviateSTR(melt(cormatrix$adj)$value, 'P')
cordata$label = paste(cordata$labelr, "\n", 
                      cordata$labelPadj, sep = "")
cordata$strike = ""
cordata$strike[cormatrix$adj > 0.05] = "X" # changed this to only show the significant ones
txtsize <- par('din')[2] / 2
# Plot
PD_corr <- ggplot(cordata, aes(x=Var1, y=Var2, fill=value), type = "upper") + geom_tile() + 
  #theme(axis.text.x = element_text(angle=90, hjust=TRUE, size=16)) +
  xlab("") + ylab("") + labs(title="PD group") +
  geom_text(label=cordata$label, size=txtsize*1.3, color="black") + 
  geom_text(label=cordata$strike, size=txtsize*3 , color="black", alpha=0.2) +
  scale_fill_gradient2(low = "#00688B", high = "#EE5C42", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 13, hjust = 1))+
  theme(axis.text.y = element_text(angle = 0, vjust = 1, 
                                   size = 13, hjust = 1))  +
  scale_x_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot','Education \n(years)','Hallucinations  \nduration','Hallucinations  \nfrequency','Hallucinations  \ndistress')) + # 'Categorical fluency score'
  scale_y_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot','Education \n(years)','Hallucinations  \nduration','Hallucinations  \nfrequency','Hallucinations  \ndistress'))
PD_corr


### PDD
dataset_corr <- PDD[, c("Fright_x_freq", "Believe_real_combined_MM_updated", "MMSE_tot","edu_year_Prabs","PNEVHI_How_long", "PNEVHI_How_often","PNEVHI_Frightening")]
d <- dataset_corr
cormatrix = rcorr(as.matrix(d[,c(1,2,3,4,5,6,7)]), type='spearman')
cormatrix$adj = p.adjust(cormatrix$P, method = "bonferroni") 
cordata = melt(cormatrix$r) 

cordata$labelr = abbreviateSTR(melt(cormatrix$r)$value, 'r')
cordata$labelP = abbreviateSTR(melt(cormatrix$P)$value, 'P non-adjusted')
cordata$labelPadj = abbreviateSTR(melt(cormatrix$adj)$value, 'P')
cordata$label = paste(cordata$labelr, "\n", 
                      cordata$labelPadj, sep = "")
cordata$strike = ""
cordata$strike[cormatrix$adj > 0.05] = "X" # changed this to only show the significant ones
txtsize <- par('din')[2] / 2
# Plot
PDD_corr <- ggplot(cordata, aes(x=Var1, y=Var2, fill=value), type = "upper") + geom_tile() + 
  #theme(axis.text.x = element_text(angle=90, hjust=TRUE, size=16)) +
  xlab("") + ylab("") + labs(title="PDD group") +
  geom_text(label=cordata$label, size=txtsize*1.3, color="black") + 
  geom_text(label=cordata$strike, size=txtsize*3 , color="black", alpha=0.3) +
  scale_fill_gradient2(low = "#00688B", high = "#EE5C42", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 13, hjust = 1))+
  theme(axis.text.y = element_text(angle = 0, vjust = 1, 
                                   size = 13, hjust = 1))  +
  scale_x_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot','Education \n(years)','Hallucinations  \nduration','Hallucinations  \nfrequency','Hallucinations  \ndistress')) +
  scale_y_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot','Education \n(years)','Hallucinations  \nduration','Hallucinations  \nfrequency','Hallucinations  \ndistress'))
PDD_corr

### EY
dataset_corr <- EY[, c("Fright_x_freq", "Believe_real_combined_MM_updated", "MMSE_tot","edu_year_Prabs","PNEVHI_How_long", "PNEVHI_How_often","PNEVHI_Frightening")] # "categorical_fluency_Animal"
d <- dataset_corr
cormatrix = rcorr(as.matrix(d[,c(1,2,3,4,5,6,7)]), type='spearman')
cormatrix$adj = p.adjust(cormatrix$P, method = "bonferroni") 
cordata = melt(cormatrix$r) 
cordata$labelr = abbreviateSTR(melt(cormatrix$r)$value, 'r')
cordata$labelP = abbreviateSTR(melt(cormatrix$P)$value, 'P non-adjusted')
cordata$labelPadj = abbreviateSTR(melt(cormatrix$adj)$value, 'P')
cordata$label = paste(cordata$labelr, "\n", 
                      cordata$labelPadj, sep = "")
cordata$strike = ""
cordata$strike[cormatrix$adj > 0.05] = "X" # changed this to only show the significant ones
txtsize <- par('din')[2] / 2
# Plot
EY_corr <- ggplot(cordata, aes(x=Var1, y=Var2, fill=value), type = "upper") + geom_tile() + 
  #theme(axis.text.x = element_text(angle=90, hjust=TRUE, size=16)) +
  xlab("") + ylab("") + labs(title="ED group") +
  geom_text(label=cordata$label, size=txtsize*1.3, color="black") + 
  geom_text(label=cordata$strike, size=txtsize*3 , color="black", alpha=0.3) +
  scale_fill_gradient2(low = "#00688B", high = "#EE5C42", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 13, hjust = 1))+
  theme(axis.text.y = element_text(angle = 0, vjust = 1, 
                                   size = 13, hjust = 1))  +
  scale_x_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot','Education \n(years)','Hallucinations  \nduration','Hallucinations  \nfrequency','Hallucinations  \ndistress')) + # 'Categorical fluency score'
  scale_y_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot','Education \n(years)','Hallucinations  \nduration','Hallucinations  \nfrequency','Hallucinations  \ndistress'))
EY_corr


### All four groups as one group
dataset_corr <- dataset_4[, c("Fright_x_freq", "Believe_real_combined_MM_updated", "MMSE_tot","edu_year_Prabs","PNEVHI_How_long","PNEVHI_How_often","PNEVHI_Frightening")] # "categorical_fluency_Animal"
d <- dataset_corr
cormatrix = rcorr(as.matrix(d[,c(1,2,3,4,5,6,7)]), type='spearman')
cormatrix$adj = p.adjust(cormatrix$P, method = "bonferroni") 
cordata = melt(cormatrix$r) 
cordata$labelr = abbreviateSTR(melt(cormatrix$r)$value, 'r')
cordata$labelP = abbreviateSTR(melt(cormatrix$P)$value, 'P non-adjusted')
cordata$labelPadj = abbreviateSTR(melt(cormatrix$adj)$value, 'P')
cordata$label = paste(cordata$labelr, "\n", 
                      cordata$labelPadj, sep = "")
cordata$strike = ""
cordata$strike[cormatrix$adj > 0.05] = "X" # changed this to only show the significant ones
txtsize <- par('din')[2] / 2
# Plot
everyone_corr <- ggplot(cordata, aes(x=Var1, y=Var2, fill=value), type = "upper") + geom_tile() + 
  #theme(axis.text.x = element_text(angle=90, hjust=TRUE, size=16)) +
  xlab("") + ylab("") + labs(title="All four patients groups") +
  geom_text(label=cordata$label, size=txtsize*1.3, color="black") + 
  geom_text(label=cordata$strike, size=txtsize*3 , color="black", alpha=0.3) +
  scale_fill_gradient2(low = "#00688B", high = "#EE5C42", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 13, hjust = 1))+
  theme(axis.text.y = element_text(angle = 0, vjust = 1, 
                                   size = 13, hjust = 1))  +
  scale_x_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot','*Education \n(years)','Hallucinations  \nduration','Hallucinations  \nfrequency','Hallucinations  \ndistress')) + # 'Categorical fluency score'
  scale_y_discrete(labels = c('Hallucinations \nseverity','Degree of \ninsight', 'MMSE Tot','*Education \n(years)','Hallucinations  \nduration','Hallucinations  \nfrequency','Hallucinations  \ndistress'))
everyone_corr

## Bring all together in one figure
library(ggpubr)
# For all 4 groups separately but next to each other with shared legend
setwd("~/Google Drive/PhD KCL/Publications in progress/Prabitha ECHR/Figures")
tiff("all_corr_groups_scaled.tiff", units="in", width=13, height=10.5, res=300)
all_corr_groups <- ggarrange(DLB_corr, PDD_corr, PD_corr,EY_corr, ncol = 2, nrow = 2, common.legend = TRUE)
all_corr_groups
dev.off()
# For the overall group 
setwd("~/Google Drive/PhD KCL/Publications in progress/Prabitha ECHR/Figures")
tiff("overall_group.tiff", units="in", width=9, height=7.5, res=300)
everyone_corr
dev.off()


#### Supplementary analyses ####
#### Mutivariate regression analyses

## Subset dataset into different diagnostic groups
DLB <- dataset_4%>% filter(Diagnose_Prabs == "DLB")
PD <- dataset_4%>% filter(Diagnose_Prabs == "PD")
PDD <- dataset_4%>% filter(Diagnose_Prabs == "PDD")
EY <- dataset_4%>% filter(Diagnose_Prabs == "EyeDisease")

## Run analyses per each group
# DLB
model_DLB <- lm(Believe_real_combined_MM_updated ~ MMSE_tot * Fright_x_freq 
                + PNEVHI_How_long + PNEVHI_How_often + PNEVHI_Frightening
                + Age,
                data = DLB)
summary(model_DLB)

# PDD 
model_PDD <- lm(Believe_real_combined_MM_updated ~ MMSE_tot * Fright_x_freq 
                + PNEVHI_How_long + PNEVHI_How_often + PNEVHI_Frightening
                + Age, 
                data = PDD)
summary(model_PDD)

# PD 
model_PD <- lm(Believe_real_combined_MM_updated ~ MMSE_tot * Fright_x_freq 
               + PNEVHI_How_long + PNEVHI_How_often + PNEVHI_Frightening
               + Age,
               data = PD)
summary(model_PD)

# EY 
model_EY_visual_acc <- lm(Believe_real_combined_MM_updated ~ MMSE_tot * Fright_x_freq 
                          + PNEVHI_How_long + PNEVHI_How_often + PNEVHI_Frightening 
                          + Age, 
                          data = EY)
summary(model_EY_visual_acc)

## Create regression tables
library(gtsummary)
tbl_regression(model_DLB, intercept=TRUE)
tbl_regression(model_PDD, intercept=TRUE)
tbl_regression(model_PD, intercept=TRUE)
tbl_regression(model_EY_visual_acc, intercept=TRUE)



####################### 