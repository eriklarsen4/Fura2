
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(purrr)
library(reshape2)




## Ca Img Response Statistics Script Developed by Erik Larsen ##

  ## Enter the data for each agonist
##### LY #####
LY344864 = data.frame(genotype = c(rep("WT", 5), rep("Mut", 5)),
                      id = rep(1:5,2),
                      nneurons = c(255, 332, 80, 233, 340, 217, 632, 139, 101, 154),
                      nresponses = c(15, 18, 1, 2, 11, 6, 9, 7, 3, 3),
                      mag_avg = c(0.3199307, 0.6069476, 1.793108, 0.4647599, 0.1927499, 0.4432441, 0.4922991, 0.7120788, 0.3830038, 0.4321384))
LY344864$phat = LY344864$nresponses/LY344864$nneurons
LY344864$phat.as.percent = LY344864$nresponses/LY344864$nneurons*100

##### CYM #####
CYM5442 = data.frame(genotype = c(rep("WT", 4), rep("Mut", 4)),
                     id = rep(1:4,2),
                     nneurons = c(352,420,116,79,124,597,100,232),
                     nresponses = c(17,29,6,5,7,12,2,13),
                     mag_avg = c(0.3555906, 0.4063199, 0.3325119, 0.2322338, 0.3074834, 0.3576473, 0.7396424, 0.5005961))
CYM5442$phat = CYM5442$nresponses/CYM5442$nneurons
CYM5442$phat.as.percent = CYM5442$nresponses/CYM5442$nneurons*100

##### B #####
Beta_alanine = data.frame(genotype = c(rep("WT", 3), rep("Mut", 4)),
                          id = c(1:3,1:4),
                          nneurons = c(79,226,744,301,176,404,431),
                          nresponses = c(20,57,92,49,41,43,79),
                          mag_avg = c(0.5043806, 0.505493, 0.3508911, 0.3281273, 0.5383118, 0.373826, 0.4885986))
Beta_alanine$phat = Beta_alanine$nresponses/Beta_alanine$nneurons
Beta_alanine$phat.as.percent = Beta_alanine$nresponses/Beta_alanine$nneurons*100

##### AITC #####
AITC = data.frame(genotype = c(rep("WT", 5), rep("Mut", 4)),
                  id = c(1:5,1:4),
                  nneurons = c(144,226,115,238,417,86,292,517,346),
                  nresponses = c(18,45,29,28,184,4,114,41,170),
                  mag_avg = c(0.4248014, 0.5948393, 0.3601463, 0.4433002, 0.536738, 0.6128877, 0.6831767, 0.4383469, 0.8302011))
AITC$phat = AITC$nresponses/AITC$nneurons
AITC$phat.as.percent = AITC$nresponses/AITC$nneurons*100

##### Capsaicin #####
Capsaicin = data.frame(genotype = c(rep("WT", 4), rep("Mut", 4)),
                       id = rep(1:4,2),
                       nneurons = c(270, 139, 220, 476, 124, 515, 91, 491),
                       nresponses = c(163, 112, 129, 258, 82, 216, 42, 254),
                       mag_avg = c(0.5626785, 1.10915, 1.312774, 0.5619838, 1.059688, 0.6714705, 1.358373, 1.112391))
Capsaicin$phat = Capsaicin$nresponses/Capsaicin$nneurons
Capsaicin$phat.as.percent = Capsaicin$nresponses/Capsaicin$nneurons*100

##### CQ #####
Chloroquine = data.frame(genotype = c(rep("WT", 3), rep("Mut", 3)),
                         id = rep(1:3,2),
                         nneurons = c(295,232,230,219,336,196),
                         nresponses = c(27,4,22,23,65,21))
Chloroquine$phat = Chloroquine$nresponses/Chloroquine$nneurons
Chloroquine$phat.as.percent = Chloroquine$nresponses/Chloroquine$nneurons*100

##### IL31 #####
IL31 = data.frame(genotype = c(rep("WT", 3), rep("Mut", 4)),
                  id = c(c(1:3), c(1:4)),
                  nneurons = c(244, 189, 230, 404, 298, 436, 294),
                  nresponses = c(9, 5, 6, 3, 29, 62, 16))
IL31$phat = IL31$nresponses/IL31$nneurons
IL31$phat.as.percent = IL31$nresponses/IL31$nneurons*100



##### Analyze data to determine significance based on pooling data and using contingency table, Welch's correction on Fisher's exact t test #####


  ## Define a function that will create a contingency table of pooled responses to an agonist and return whether there is a significant difference between genotypes using Fisher's Exact Test
ContingencyTableTest = function(Agonist, WTN, MutN){

  ## Create the contingency table of pooled responses across mice to a particular agonist between genotypes
  contingency.pooled = matrix(c(sum(Agonist[1:WTN,3]) - sum(Agonist[1:WTN,4]), sum(Agonist[(WTN+1):(WTN+MutN),3]) - sum(Agonist[(WTN+1):(WTN+MutN),4]), sum(Agonist[1:WTN,4]), sum(Agonist[(WTN+1):(WTN+MutN),4])), ncol = 2)
  
  print(contingency.pooled)
  
  ## Perform the Fisher's Exact Test on the data
  fisher.test(contingency.pooled)
  return(fisher.test(contingency.pooled))
}
  ## Perform the analysis for a given dataset
#ContingencyTableTest(Agonist = LY344864, WTN = 5, MutN = 5)



    ## Assumes the true proportion of responding neurons is the same among all mice of each genotype
      ## Ignores large variations due to sampling bias and each mouse's intrinsic individuality/variation
        ## --> estimate the proportion of responding neurons from each mouse; calculate the variance of each estimate and weight the estimate by its variance

##### Analyze counts data to determine whether pooling is a valid strategy; if not for one, not for any #####

PoolingTest = function(Agonist, Genotype){
  set.seed(2000)
  if (Genotype == "WT"){
    contingency = matrix(c(subset(Agonist, genotype == "WT")$nneurons - subset(Agonist, genotype == "WT")$nresponses, subset(Agonist, genotype == "WT")$nresponses), ncol = 2)}
  else if (Genotype == "Mut"){
    contingency = matrix(c(subset(Agonist, genotype == "Mut")$nneurons - subset(Agonist, genotype == "Mut")$nresponses, subset(Agonist, genotype == "Mut")$nresponses), ncol = 2)}
    print(contingency)
    return(fisher.test(contingency, simulate.p.value = T))
}

#attach(LY344864)
#PoolingTest(Agonist = LY344864, Genotype = "WT")
#PoolingTest(Agonist = LY344864, Genotype = "Mut")
#detach(Capsaicin)


##### Analyze counts data to determine significance based on alternative options to pooling #####
  ## Options are:
  ## Unpaired t/linear regression (same)
  ## Weighted linear regression (weighting proportion estimates by the number of neurons per mouse)
  ## Transformed linear regression (to account for low-count biases)
    ## Use weighted linear regression

AlternativesToPoolingCounts = function(Agonist, Analysis_Type, WTN, MutN){
  ## Determine the weight of data from each mouse: number of neurons examined in each mouse scaled by the total number of neurons examined across all mice in which the given agonist was used to challenge; multiply by 10 to give an "effective mouse weight" where integers are given in fractions of whole mice
  mouse_weights = nneurons/sum(nneurons)*10
  mouse_weights <<- mouse_weights
  
  if (Analysis_Type == "t test"){
    print(t.test(phat.as.percent ~ genotype, Agonist, var.equal = F))
  }
  if (Analysis_Type == "linear regression"){
    print(summary(lm(phat.as.percent ~ genotype, Agonist)))
  }
  if (Analysis_Type == "weighted linear regression"){
      ## Create the model as an object to extract the fitted mean and rmse
    weighted_linear_model = lm(phat.as.percent ~ genotype, Agonist, weights = mouse_weights)
      ## Create the summary of the model as an object to extract the standard error of the fitted estimate
    sum_weighted_linear_model = summary(weighted_linear_model)
      ## Extract the fitted mean and store it in the data.frame for graphing
    Agonist$fitted_mean = as.numeric(c(rep(weighted_linear_model$fitted.values[1], WTN), rep(weighted_linear_model$fitted.values[(WTN+1)], MutN)))
      ## Extract the se and store it in the data.frame for graphing
    Agonist$se = as.numeric(c(rep(sum_weighted_linear_model$coefficients[,2][2], WTN), rep(sum_weighted_linear_model$coefficients[,2][1], MutN)))
    
    
    print(summary(lm(phat.as.percent ~ genotype, Agonist, weights = mouse_weights)))
  }
  if (Analysis_Type == "transformed weighted linear regression"){
      ## Create the model as an object to extract the fitted mean and rmse
    weighted_linear_model = lm(phat.as.percent ~ genotype, Agonist, weights = mouse_weights)
      ## Create the summary of the model as an object to extract the standard error of the fitted estimate
    sum_weighted_linear_model = summary(weighted_linear_model)
      ## Extract the fitted mean and store it in the data.frame for graphing
    Agonist$fitted_mean = as.numeric(c(rep(weighted_linear_model$fitted.values[1], WTN), rep(weighted_linear_model$fitted.values[(WTN+1)], MutN)))
      ## Extract the se and store it in the data.frame for graphing
    Agonist$se = as.numeric(c(rep(sum_weighted_linear_model$coefficients[,2][2], WTN), rep(sum_weighted_linear_model$coefficients[,2][1], MutN)))
    
      ## Transform the percent estimate to account for low count bias (negative binomial distribution)
    Agonist$ptransf = asin(sqrt(Agonist$phat))
      ## Move the decimal over
    Agonist$ptransf.as.percent = Agonist$ptransf*100
    print(summary(lm(ptransf.as.percent ~ genotype, Agonist, weights = mouse_weights)))
    
      ## Create the model as an object to extract the fitted mean and rmse
    weighted_linear_model = lm(ptransf.as.percent ~ genotype, Agonist, weights = mouse_weights)
      ## Create the summary of the model as an object to extract the standard error of the fitted estimate
    sum_weighted_linear_model = summary(weighted_linear_model)
    
      ## Extract the fitted mean and store it in the data.frame for graphing
    Agonist$ptransf.fitted_mean = as.numeric(c(rep(weighted_linear_model$fitted.values[1], WTN), rep(weighted_linear_model$fitted.values[(WTN+1)], MutN)))
      ## Extract the se and store it in the data.frame for graphing
    Agonist$ptransf.se = as.numeric(c(rep(sum_weighted_linear_model$coefficients[,2][2], WTN), rep(sum_weighted_linear_model$coefficients[,2][1], MutN)))
    
    
  }
  
  Agonist <<- Agonist
  sum_weighted_linear_model <<- sum_weighted_linear_model
  weighted.p <<- Agonist$phat*mouse_weights
  ptransf <<- Agonist$ptransf
  
}

#attach(Chloroquine)
#AlternativesToPoolingCounts(Agonist = Chloroquine, Analysis_Type = "t test", WTN = 4, MutN = 4)
#attach(Chloroquine)
#AlternativesToPoolingCounts(Agonist = Chloroquine, Analysis_Type = "linear regression", WTN = 4, MutN = 4)
attach(Capsaicin)
AlternativesToPoolingCounts(Agonist = Capsaicin, Analysis_Type = "weighted linear regression", WTN = 4, MutN = 4)
#attach(Capsaicin)
#AlternativesToPoolingCounts(Agonist = IL31, Analysis_Type = "transformed weighted linear regression", WTN = 3, MutN = 4)
detach(Capsaicin)




##### Analyze magnitudes data to determine significance based on alternative options to pooling #####
  ## Options are:
  ## Unpaired t/linear regression (same)
  ## Weighted linear regression (weighting proportion estimates by the number of neurons per mouse)
  ## Transformed linear regression (to account for low-count biases)
  ## Use weighted linear regression
AlternativesToPoolingMags = function(Agonist, Analysis_Type, WTN, MutN){
  ## Determine the weight of data from each mouse: number of neurons examined in each mouse scaled by the total number of neurons examined across all mice in which the given agonist was used to challenge; multiply by 10 to give an "effective mouse weight" where integers are given in fractions of whole mice
  mouse_weights = nresponses/sum(nresponses)*10
  mouse_weights <<- mouse_weights
  
  if (Analysis_Type == "t test"){
    print(t.test(mag_avg ~ genotype, Agonist, var.equal = F))
  }
  if (Analysis_Type == "linear regression"){
    print(summary(lm(mag_avg ~ genotype, Agonist)))
  }
  if (Analysis_Type == "weighted linear regression"){
    ## Create the model as an object to extract the fitted mean and rmse
    weighted_linear_model = lm(mag_avg ~ genotype, Agonist, weights = mouse_weights)
    ## Create the summary of the model as an object to extract the standard error of the fitted estimate
    sum_weighted_linear_model = summary(weighted_linear_model)
    ## Extract the fitted mean and store it in the data.frame for graphing
    Agonist$fitted_mean = as.numeric(c(rep(weighted_linear_model$fitted.values[1], WTN), rep(weighted_linear_model$fitted.values[(WTN+1)], MutN)))
    ## Extract the se and store it in the data.frame for graphing
    Agonist$se = as.numeric(c(rep(sum_weighted_linear_model$coefficients[,2][2], WTN), rep(sum_weighted_linear_model$coefficients[,2][1], MutN)))
    
    
    print(summary(lm(mag_avg ~ genotype, Agonist, weights = mouse_weights)))
  }
  
  Agonist <<- Agonist
  sum_weighted_linear_model <<- sum_weighted_linear_model
  weighted.p <<- Agonist$phat*mouse_weights
  ptransf <<- Agonist$ptransf
  
}
attach(LY344864)
AlternativesToPoolingMags(Agonist = LY344864, Analysis_Type = "weighted linear regression", WTN = 5, MutN = 5)
detach(LY344864)




##### Plot the results #####

GGBar = ggplot(data = Agonist, aes(x = genotype, y = phat.as.percent, fill = genotype)) +
  #geom_boxplot(width = 0.5) +
  #geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  scale_shape_manual(values = c(19,1)) +
  geom_jitter(position = position_jitter(0.1)) +
  geom_segment(aes(x = 1, y = Agonist$fitted_mean[1], xend = 2, yend = Agonist$fitted_mean[nrow(Agonist)]), color = "darkgoldenrod") +
  scale_x_discrete(limits = c("WT", "Mut")) +
  geom_point(aes(y = fitted_mean), color = "darkgoldenrod", size = 3) +
  geom_errorbar(aes(y = fitted_mean, ymin = (fitted_mean - se), ymax = (fitted_mean + se), width = 0.1), color = "darkgoldenrod") +
  labs(title = expression(paste(italic("In vitro"), " Capsaicin Responses")), x = "Genotype", y = "% Responsive Neurons\nPer Mouse") +
  scale_fill_manual(values = c("gray", "black")) +
  #geom_line(data = weighted_linear_model, aes(y = predwlm)) + 
  coord_cartesian(ylim = c(0.0,100)) +
  #geom_smooth(method = "lm", color = "yellow") + 
  geom_signif(stat = "identity", aes(x = 1.1, xend = 1.9, y = 85, yend = 85, annotation = "p = 0.1279")) +
  #geom_line(yintercept = 83, aes(x = genotype)) +
  geom_text(x = 0.8, y = 59.9, label = "59.9", color = "darkgoldenrod", size = 4) +
  geom_text(x = 2.21, y = 48.6, label = "48.6", color = "darkgoldenrod", size = 4) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.9), panel.grid = element_blank(), legend.position = "none")

GGBar


ggsave(filename = "Capsaicin responses w sem", plot = last_plot(), device = "tiff", dpi = "print", units = "in", width = 3, height = 3, path = "M:/Erik/Data/Calcium Imaging Analysis")




