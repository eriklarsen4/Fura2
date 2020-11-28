
library(tidyverse)
library(ggplot2)
library(purrr)
library(reshape2)
library(dplyr)
library(ecp)
library(strucchange)

##### Data acquisition and prep for perfusion experiment analysis. Developed by Erik Larsen #####

##### Fluorescence Value Acquisition and Prep #####

  ## Change the working directory to the directory containing all the files and scripts
setwd("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/")
  ## Import the fluorescence files if necessary
files = list.files(pattern = ".csv")
  ## Import the fluorescence scripts tailored specifically to each run
scripts = list.files(pattern = ".R")
  ## Character vector containing WT experimental dates
WT = c("2019_10_30", "2019_11_8", "2019_11_15", "2019_12_4", "2019_12_11")
  ## Character vector containing Mut experimental dates
Mut = c("2019_11_6", "2019_11_11", "2019_11_20", "2019_12_6", "2019_12_13")

  ## Subset the scripts (and files, if necessary)
WTscripts = c(scripts[which( !is.na(str_extract(scripts, regex(pattern = WT[1])) ))], scripts[which( !is.na(str_extract(scripts, regex(pattern = WT[2])) ))], scripts[which( !is.na(str_extract(scripts, regex(pattern = WT[3])) ))], scripts[which( !is.na(str_extract(scripts, regex(pattern = WT[4])) ))], scripts[which( !is.na(str_extract(scripts, regex(pattern = WT[5])) ))])
  ## Subset the Mut scripts (and files, if necessary)
Mutscripts = c(scripts[which( !is.na(str_extract(scripts, regex(pattern = Mut[1])) ))], scripts[which( !is.na(str_extract(scripts, regex(pattern = Mut[2])) ))], scripts[which( !is.na(str_extract(scripts, regex(pattern = Mut[3])) ))], scripts[which( !is.na(str_extract(scripts, regex(pattern = Mut[4])) ))], scripts[which( !is.na(str_extract(scripts, regex(pattern = Mut[5])) ))])


WTfiles = c(files[which( !is.na(str_extract(files, regex(pattern = WT[1])) ))], files[which( !is.na(str_extract(files, regex(pattern = WT[2])) ))], files[which( !is.na(str_extract(files, regex(pattern = WT[3])) ))], files[which( !is.na(str_extract(files, regex(pattern = WT[4])) ))], files[which( !is.na(str_extract(files, regex(pattern = WT[5])) ))])

Mutfiles = c(files[which( !is.na(str_extract(files, regex(pattern = Mut[1])) ))], files[which( !is.na(str_extract(files, regex(pattern = Mut[2])) ))], files[which( !is.na(str_extract(files, regex(pattern = Mut[3])) ))], files[which( !is.na(str_extract(files, regex(pattern = Mut[4])) ))], files[which( !is.na(str_extract(files, regex(pattern = Mut[5])) ))])

  ## Open the desired script
  ## WT Capsaicin: 2,3,4,6,  8,11,12,  21,22,26,27,33,34,36,39
  ## Mutant Capsaicin: 1,2,7,  13,14,16,  29,31,  36,37,38,40,42,44,45
  ## WT AITC: 1,6,7,  8,11,13,  17,18,20,  25,28,29,30,  32,35,37,38
  ## Mutant AITC: 1,2,  18,19,21,23,24,25,  26,27,30,33,34,35,  36,39,41,46
  ## WT Beta: 17,20,  21,24,25,28,30,  31,33,35,37,38,39
  ## Mutant Beta: 11,16,  18,19,22,23,  31,32,34,35,  37,38,40,41,42,46
  ## WT CYM5442: 1,2,3,4,7,  9,10,13,14,  15,16,18,19,  23,29
  ## Mutant CYM5442: 4,5,6,  13,14,15,17  21,25,  26,27,28,29,33
  ## WT LY344864: 1,3,4,7,  9,10,12,14,  15,16,19,  22,23,24,26,27,  31,32,34,36
  ## Mutant LY344864: 3,4,5,6,7,  8,13,15,17,  20,24,  28,30,  39,45

#file.edit(WTscripts[37])
#file.edit(Mutscripts[16])


Find_Matching_Files = function(Script_Number, Genotype){
  ## For finding files that are analyzed by a particular script
  Matching_Files = vector()
  ## Find the matching WT files analyzed by a given script (indexed by script; change to the desired index)
  if (Genotype == "WT"){
    Matching_Files = WTfiles[
      which( !is.na(
        str_extract(WTfiles, regex(
          pattern = paste(
            gsub(" script.R", "", WTscripts[Script_Number], fixed = TRUE), ".+$", sep = ""))
        )
      )
      )
    ]
  }
  ## Find the matching Mut files analyzed by a given script (indexed by script; change to the desired index)
  if (Genotype == "Mut"){
    Matching_Files = Mutfiles[
      which( !is.na(
        str_extract(Mutfiles, regex(
          pattern = paste(
            gsub(" script.R", "", Mutscripts[Script_Number], fixed = TRUE), ".+$", sep = ""))
        )
      )
      )
    ]
  }
  Matching_Files <<- Matching_Files
  return(Matching_Files)
}
Find_Matching_Files(Script_Number = 1, Genotype = "Mut")

##### Find Peak Response magnitudes and times #####
  ## Create a function that determines from the ROI number (+1; time vector is the first column in the data frame):
    ## 1) The maximum smooth response for each agonist
    ## 2) The frame when that/those responses occur

Find_Peak_Times_And_Mags = function (Cell_Number_Plus_One){
  ## Concatenate the max smoothed responses and their respective frame numbers (when they occur) for non-HiK agonists
    ## Initialize variables
  Peak_Mags = as.numeric(vector())
  Peak_Frames = as.numeric(vector())
  HiK_Peak_Mag = as.numeric(vector())
  HiK_Peak_Frame = as.numeric(vector())
  for (j in 1:(nrow(response_frame)-1) ){
      ## Concatenate the max smoothed responses to non-HiK agonists into a variable
      Peak_Mags[j] = cbind(max(
        SmoothFLUO[ response_frame[j,1]:response_frame[j,2], Cell_Number_Plus_One ] ))
      ## Concatenate the respective frame numbers when these responses occur
      Peak_Frames[j] = cbind(which(SmoothFLUO[, Cell_Number_Plus_One] == max( SmoothFLUO[ response_frame[j,1]:response_frame[j,2], Cell_Number_Plus_One ] ) ))
  }
  Peak_Mags = Peak_Mags[!is.na(Peak_Mags) == TRUE]
  Peak_Frames = Peak_Frames[!is.na(Peak_Frames) == TRUE]
  
      ## Create a variable for the max smoothed response to HiK
      HiK_Peak_Mag = max(
        SmoothFLUO[ response_frame[nrow(response_frame),1]:response_frame[nrow(response_frame),2], Cell_Number_Plus_One ]  )
      ## Create a variable for when the max HiK response occurs
      HiK_Peak_Frame = which(SmoothFLUO[, Cell_Number_Plus_One] == max( SmoothFLUO[ response_frame[nrow(response_frame),1]:response_frame[nrow(response_frame),2], Cell_Number_Plus_One ] ) )
      
    ## Export all these values for a given ROI into the global environment for plotting
  Peak_Mags <<- Peak_Mags
  Peak_Frames <<- Peak_Frames
  HiK_Peak_Mag <<- HiK_Peak_Mag
  HiK_Peak_Frame <<- HiK_Peak_Frame
  return(print(c("Peak Agonist Frames", Peak_Frames, "Peak HiK Frame", HiK_Peak_Frame)))
}

  ##### Smooth the FLUO data (span = 10 to retain trace shape but reduce noise) acquired from each run's script #####
  ## Create a function that will perform local weighted regression (smoothing) on the FLUO data
    ## The majority of the ROIs do not need smoothing; however, it is useful for noisy ROIs, particularly quantifying magnitudes at precise points
    ## Smoothing may lessen transient response magnitudes

loess.filter = function (x, span){
  loess(formula = paste(x, "t", sep = "~"), data = FLUO, degree = 1, span = span)$fitted
}

    ## Append the time vector to the FLUO data
FLUO = cbind(t, FLUO)
    ## Extract the strings of the columns/ROIs
ROIs = c(colnames(FLUO[c(2:ncol(FLUO))]))
  ## Apply the smoothing to all ROIs in the FLUO data
SmoothFLUO = as.data.frame(lapply(ROIs, loess.filter, span = 10/nrow(FLUO)), col.names = colnames(FLUO[c(2:ncol(FLUO))]))
  ## Append the time vector to the smoothed FLUO data
SmoothFLUO = cbind(t, SmoothFLUO)
  ## Find, store the smoothed ROI FLUO slope change (change in smooth 340/380 / change in time) in a matrix
Slope_Change_DF = as.data.frame(diff(
  as.matrix(SmoothFLUO[ , c(2:ncol(SmoothFLUO)) ]
  )
)
)
  ## Retain the same column names as from the FLUO and SmoothFLUO data frames
colnames(Slope_Change_DF) = colnames(SmoothFLUO[c(2:ncol(SmoothFLUO))])
options(scipen = 999)

##### Perform Analysis #####

  ## Define a function that will compute (smoothed or raw) responses based on response thresholds used in the Itch field (10% above baseline, also following a 2% increase in slope). Should align with video responses, where ROIs turn green. Count data for the imaging run are exported and returned. Exported raw traces for all ROIs are commented out. Also generates variables used for ggplotting the smoothed data

Analysis = function(Slope_Change_Threshold){
  
    ##### Create a matrix to house the response calls to all agonists #####
  RESPONSE_CALLS = matrix(nrow = length(agonist[seq(2,length(agonist),2)]), ncol = ncol(SmoothFLUO)-1)
  RESPONSE_CALLS = as.data.frame(RESPONSE_CALLS)
  colnames(RESPONSE_CALLS) = colnames(SmoothFLUO)[2:ncol(SmoothFLUO)]
  rownames(RESPONSE_CALLS) = rownames(response_frame)
    ## Create a matrix that houses where each ROI experiences a ratiometric increase in slope above a threshold
      ## Loop through an populate the slope matrix
  Slope_idx_DF = matrix(nrow = nrow(SmoothFLUO), ncol = (ncol(SmoothFLUO)-1))
  for (i in 1:ncol(Slope_idx_DF)) {
    for (j in 1:nrow(Slope_idx_DF)) {
      Slope_idx_DF[j,i] = paste(which(Slope_Change_DF[j,i] > as.numeric(Slope_Change_Threshold)), collapse = "")
    }
  }
    ## Create a variable that stores the ROI indeces that experience an increase in slope at any point in the run
      ## Necessary for conditional statements in the triple-for-loop assigning calls
  response_idx = as.numeric(vector())
  for (k in 1:ncol(Slope_idx_DF)){
    response_idx[k] = if_else( 
      any(which(Slope_idx_DF[,k] > 0)) == TRUE, 1, 0)
  }
  #which(response_idx == 0)
  
    ##### Assign response calls based on the "10% above baseline" criteria, but also occurring after a/n (arbitrary) 2% increase in slope #####
  for (i in 1:(ncol(SmoothFLUO)-1) ){
    for (j in 1:(nrow(response_frame)-1) ){
      for ( k in 1:nrow(SmoothFLUO) ){
        RESPONSE_CALLS[j,c(which(response_idx == 0))] = "Non-responder"
        if (
          k == c(which(response_idx == 1)) & any(response_frame[j,1]:response_frame[j,2] %in% which(Slope_idx_DF[ , i] > 0)  == TRUE) ) {
          RESPONSE_CALLS[j,i] = if_else(max( SmoothFLUO[ response_frame[j,1]:response_frame[j,2], i+1 ] ) >
                                          (max( SmoothFLUO[ exposure_frame[j,1]:exposure_frame[j,2], i+1 ] )*Agonists_df[j,2] + max( SmoothFLUO[ exposure_frame[j,1]:exposure_frame[j,2], i+1 ] ) ), "Responder", "Non-responder")
        } else if (
          k == c(which(response_idx == 1)) & any(response_frame[j,1]:response_frame[j,2] %in% which(Slope_idx_DF[ , i] > 0)  == FALSE) ) {
          RESPONSE_CALLS[j,i] = "Non-responder"
        }
      }
      RESPONSE_CALLS[nrow(response_frame),i] = if_else(max( SmoothFLUO[ response_frame[nrow(response_frame),1]:response_frame[nrow(response_frame),2], i+1 ] ) >
                                                         (max( SmoothFLUO[ exposure_frame[ length(agonist)-1,1]:exposure_frame[ length(agonist)-1,2], i+1 ] )*Agonists_df[nrow(Agonists_df),2] + max( SmoothFLUO[ exposure_frame[ length(agonist)-1,1]:exposure_frame[ length(agonist)-1,2], i+1 ] ) ), "Neuron", "Non-neuron")
    }
  }

  
    ##### Create a dataframe to store response stats to each agonist #####
  Response_Stats_df = matrix(nrow = length(agonist[seq(2,length(agonist),2)]), ncol = 5)
    ## Rename columns and rows
  rownames(Response_Stats_df) = c(agonist[seq(2,length(agonist),2)])
  colnames(Response_Stats_df) = c("# of Responders", "% of Neurons responding", "% Neurons of ROIs", "Max Response Mag.", "Avg. Peak Response")
    ## Remove NAs
  Response_Stats_df[is.na(Response_Stats_df)] = 0
  
    ## Fill it with the number of responders, neurons, percent of neurons responding to an agonist, the max response of each agonist, and the percent of neurons out of all the ROIs. Also printed below
  
    ## See how many ROIs are classified as "responders", which ones they are, and the percentages of responses
      ## How many responders?
  for (i in 1:(nrow(RESPONSE_CALLS)-1)){
    if (length(which(RESPONSE_CALLS[i,] == "Responder")) == 0){
      Response_Stats_df[i,1] = 0
    } else {
      Response_Stats_df[i,1] = length(which(RESPONSE_CALLS[i,] == "Responder")) 
    }
  }
    ## How many neurons?
  for (i in nrow(RESPONSE_CALLS)){
    Response_Stats_df[nrow(Response_Stats_df),1] = length(which(RESPONSE_CALLS[i,] == "Neuron"))
  }
    ## Which are responders?
  for (i in 1:(nrow(RESPONSE_CALLS)-1)){
    print(which(RESPONSE_CALLS[i,] == "Responder"))
  }
    ## Which are neurons?
  for (i in nrow(RESPONSE_CALLS)){
    print(which(RESPONSE_CALLS[i,] == "Neuron")) 
  }
    ## Of all responders, what were the maximum responses to each respective agonist?
  for (j in 1:(nrow(RESPONSE_CALLS)-1)) {
    Response_Stats_df[j,4] = max(SmoothFLUO[ response_frame[j,1]:response_frame[j,2] , c(2:ncol(SmoothFLUO))])
  }
  
    ## Of all ROIs, what were the maximum responses to HiK?
    Response_Stats_df[nrow(RESPONSE_CALLS),4] = max(
      SmoothFLUO[ c(response_frame[nrow(response_frame),1]:response_frame[nrow(response_frame),2]) , c(2:ncol(SmoothFLUO)) ] )
    ## Determine percentages
      ## Of neurons responding to each agonist
  for (i in 1:(nrow(RESPONSE_CALLS)-1)){
    if (length(which(RESPONSE_CALLS[i,] == "Responder")) == 0){
      Response_Stats_df[i,2] = 0
    } else {
      Response_Stats_df[i,2] = (length(which(RESPONSE_CALLS[i,] == "Responder"))/length(which(RESPONSE_CALLS[nrow(RESPONSE_CALLS),] == "Neuron")))*100 
    }
  }
    Response_Stats_df[nrow(Response_Stats_df),2] = (length(which(RESPONSE_CALLS[nrow(RESPONSE_CALLS),] == "Neuron")) / length(which(RESPONSE_CALLS[nrow(RESPONSE_CALLS),] == "Neuron")))*100
  
  for (i in 1:(nrow(RESPONSE_CALLS)-1)){
    Response_Stats_df[i,3] = NA
  }
    ## Find the percent of neurons in the field
  for (i in 1:ncol(Response_Stats_df)){
    for (j in 1:nrow(Response_Stats_df)){
      Response_Stats_df[nrow(Response_Stats_df),3] =  (Response_Stats_df[nrow(Response_Stats_df),1]/ncol(SmoothFLUO))*100
    }
  }
    ## Show the percent of neurons in the field
  for (i in nrow(RESPONSE_CALLS)){
    print((length(which(RESPONSE_CALLS[i,] == "Neuron"))/length(SmoothFLUO))*100)
  }
    ## Create the variables that will store peak agonist responses and average peak agonist responses
    First_Ag_Max_Responses = vector()
    First_Ag_Max_Avg = vector()
    Second_Ag_Max_Responses = vector()
    Second_Ag_Max_Avg = vector()
    HiK_Max_Responses = vector()
    HiK_Max_Avg = vector()
    
    if ( length(agonist) > 6 ){
      Third_Ag_Max_Responses = vector()
      Third_Ag_Max_Avg = vector()
    }
    
      ## Find and store those values for each agonist
  for (i in 1:length(which(RESPONSE_CALLS[1,] == "Responder"))){
      if (length(which(RESPONSE_CALLS[1,] == "Responder")) == 0) {
        First_Ag_Max_Responses = 0
      } else if (length(which(RESPONSE_CALLS[1,] == "Responder")) == 1){
        First_Ag_Max_Responses[i] = max( SmoothFLUO[ response_frame[1,1]:response_frame[1,2], which(RESPONSE_CALLS[1,] == "Responder")+1][] ) 
      } else if (length(which(RESPONSE_CALLS[1,] == "Responder")) > 1) {
        First_Ag_Max_Responses[i] = max( SmoothFLUO[ response_frame[1,1]:response_frame[1,2], which(RESPONSE_CALLS[1,] == "Responder")+1][i] )
      }
  }
    for (i in 1:length(which(RESPONSE_CALLS[1,] == "Responder"))){
      if (length(which(RESPONSE_CALLS[1,] == "Responder")) == 0){
        First_Ag_Max_Avg = 0
      } else {
        First_Ag_Max_Avg = (sum( First_Ag_Max_Responses) / length(which(RESPONSE_CALLS[1,] == "Responder")))
      }
    }
      
    for (i in 1:length(which(RESPONSE_CALLS[2,] == "Responder"))){
      if (length(which(RESPONSE_CALLS[2,] == "Responder")) == 0){
        Second_Ag_Max_Responses = 0
      } else if (length(which(RESPONSE_CALLS[2,] == "Responder")) == 1) {
        Second_Ag_Max_Responses[i] = max( SmoothFLUO[ response_frame[2,1]:response_frame[2,2], which(RESPONSE_CALLS[2,] == "Responder")+1][] )
      } else if (length(which(RESPONSE_CALLS[2,] == "Responder")) > 1) {
        Second_Ag_Max_Responses[i] = max( SmoothFLUO[ response_frame[2,1]:response_frame[2,2], which(RESPONSE_CALLS[2,] == "Responder")+1][i] )
      }
    }
    
  for  (i in 1:length(which(RESPONSE_CALLS[2,] == "Responder"))){
    if (length(which(RESPONSE_CALLS[2,] == "Responder")) == 0){
      Second_Ag_Max_Avg = 0
    } else {
      Second_Ag_Max_Avg = (sum( Second_Ag_Max_Responses) / length(which(RESPONSE_CALLS[2,] == "Responder")))
    }
  }
    
    if ( length(agonist) > 6) {
      for (i in 1:length(which(RESPONSE_CALLS[3,] == "Responder"))){
        if (length(which(RESPONSE_CALLS[3,] == "Responder")) == 0) {
          Third_Ag_Max_Responses = 0
        } else if (length(which(RESPONSE_CALLS[3,] == "Responder")) == 1){
          Third_Ag_Max_Responses[i] = max( SmoothFLUO[ response_frame[3,1]:response_frame[3,2], which(RESPONSE_CALLS[3,] == "Responder")+1][] ) 
        } else if (length(which(RESPONSE_CALLS[3,] == "Responder")) > 1) {
          Third_Ag_Max_Responses[i] = max( SmoothFLUO[ response_frame[3,1]:response_frame[3,2], which(RESPONSE_CALLS[3,] == "Responder")+1][i] )
        }
      }
      for (i in 1:length(which(RESPONSE_CALLS[3,] == "Responder"))){
        if (length(which(RESPONSE_CALLS[3,] == "Responder")) == 0){
          Third_Ag_Max_Avg = 0
        } else {
          Third_Ag_Max_Avg = (sum( Third_Ag_Max_Responses) / length(which(RESPONSE_CALLS[3,] == "Responder")))
        }
      }
    }
      
  for (i in 1:length(which(RESPONSE_CALLS[nrow(RESPONSE_CALLS),] == "Neuron"))){
      if (length(which(RESPONSE_CALLS[nrow(RESPONSE_CALLS),] == "Neuron")) == 0) {
        HiK_Max_Responses[i] = max( SmoothFLUO[ response_frame[nrow(RESPONSE_CALLS),1]:response_frame[nrow(RESPONSE_CALLS),2], which(RESPONSE_CALLS[nrow(RESPONSE_CALLS),] == "Neuron")+1][] )
      } else if (length(which(RESPONSE_CALLS[nrow(RESPONSE_CALLS),] == "Neuron")) > 1) {
        HiK_Max_Responses[i] = max( SmoothFLUO[ response_frame[nrow(RESPONSE_CALLS),1]:response_frame[nrow(RESPONSE_CALLS),2], which(RESPONSE_CALLS[nrow(RESPONSE_CALLS),] == "Neuron")+1][i] )
      }
    }
    HiK_Max_Avg = (sum( HiK_Max_Responses) / length(which(RESPONSE_CALLS[nrow(RESPONSE_CALLS),] == "Neuron")))
      
    ## Put them in the Response Stats DF
  Response_Stats_df[1,5] = First_Ag_Max_Avg
  Response_Stats_df[2,5] = Second_Ag_Max_Avg
  if ( length(agonist) > 6 ){
    Response_Stats_df[3,5] = Third_Ag_Max_Avg
  }
  Response_Stats_df[nrow(RESPONSE_CALLS),5] = HiK_Max_Avg

  
  # Make individual plots as separate png files; these are response profiles for one cell (ROI), each. These give rough ideas of responses by each ROI
  #  dir.create("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fura2 Images\\2019_11_15\\A1plots")
  #  for (i in seq_along(FLUO[,1:ncol(FLUO)])) {
  #    #tempname = paste("cell_",i,".png", sep="")
  #    outpath = file.path("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fura2 Images\\2019_11_15\\A1plots", paste("cell_",i,".png", sep=""))
  #    png(filename = outpath)
  #    plot(t, FLUO[,i], xlim = c(t[1],t[length(t)]), ylim = c(0,4))
  #    dev.off()
  #  }
  First_Ag_Max_Responses <<- First_Ag_Max_Responses
  First_Ag_Max_Avg <<- First_Ag_Max_Avg
  Second_Ag_Max_Responses <<- Second_Ag_Max_Responses
  Second_Ag_Max_Avg <<- Second_Ag_Max_Avg
  if ( length(agonist) > 6 ) {
    Third_Ag_Max_Responses <<- Third_Ag_Max_Responses
    Third_Ag_Max_Avg <<- Third_Ag_Max_Avg
  }
  HiK_Max_Responses <<- HiK_Max_Responses
  HiK_Max_Avg <<- HiK_Max_Avg
  FLUO <<- FLUO
  SmoothFLUO <<- SmoothFLUO
  RESPONSE_CALLS <<- RESPONSE_CALLS
  Slope_Change_DF <<- Slope_Change_DF
  Slope_idx_DF <<- Slope_idx_DF
  response_idx <<- response_idx
  Response_Stats_df <<- Response_Stats_df
  
  
  return(Response_Stats_df)
  
  #for (i in 1:nrow(Agonists_df)){
  # Determine percentages
  # Of responders to each agonist; % responders of all neurons, and percent of neurons in the field
  #print((length(which(Agonist_responses[i,] == "Responder"))/length(which(Agonist_responses[nrow(Agonists_df),] == "Neuron")))*100)
  #print((length(which(Agonist_responses[i,] == "Neuron"))/length(FLUO))*100)
  #}
}

Analysis(Slope_Change_Threshold = 0.02)



##### GGPlots #####
  ## Check for structural change in the data with e.divisive function?
    #e.divisive(as.matrix(Slope_Change_DF[,59]), min.size = 10, alpha = 1)
    #DIV = e.divisive(diff(as.matrix(SmoothFLUO[,110])), min.size = 10, alpha = 1)
    #t[DIV$estimates[-c(1,150)]]

Find_Peak_Times_And_Mags(Cell_Number_Plus_One = 96)

  ## Make plots of perfusion data
CaImaging = ggplot(data = FLUO, aes(x = FLUO$t, y = FLUO$MEAN_GREEN.1.95)) +
  geom_rect(data = FLUO, xmin = t[exposure_frame[2,1]], xmax = t[exposure_frame[2,2]], ymin = -Inf, ymax = Inf, fill = "rosybrown1", alpha = 0.01) +
  geom_rect(data = FLUO, xmin = t[exposure_frame[4,1]], xmax = t[exposure_frame[4,2]], ymin = -Inf, ymax = Inf, fill = "peachpuff2", alpha = 0.01) +
  geom_rect(data = FLUO, xmin = t[exposure_frame[6,1]], xmax = t[exposure_frame[6,2]], ymin = -Inf, ymax = Inf, fill = "darkseagreen", alpha = 0.01) +
  geom_text(x = t[exposure_frame[2,1]]+(t[exposure_frame[2,2]]-t[exposure_frame[2,1]])/2, y  = Response_Stats_df[nrow(RESPONSE_CALLS),4]-.05, label = "\u03B2-alanine\napplication", color = "rosybrown1", size = 5, fontface = "bold") +
  geom_text(x = t[exposure_frame[4,1]]+(t[exposure_frame[4,2]]-t[exposure_frame[4,1]])/2, y  = Response_Stats_df[nrow(RESPONSE_CALLS),4]-.05, label = "AITC\napplication", color = "peachpuff2", size = 5, fontface = "bold") +
  geom_text(x = t[exposure_frame[6,1]]+(t[exposure_frame[6,2]]-t[exposure_frame[6,1]])/2, y  = Response_Stats_df[nrow(RESPONSE_CALLS),4]-.05, label = "140mM K\napplication", color = "darkseagreen", size = 3.7, fontface = "bold") +
  #geom_vline(xintercept = t[which(Slope_Change_DF[,95] > 0.02)], color = "blue") +
  #geom_vline(xintercept = t[Peak_Frames], color = "red", size = 0.9) +
  #geom_vline(xintercept = t[DIV$estimates][-c(1,150)]], color = "brown") +
  #geom_point(color = "black") +
  geom_line(color = "black", size = 1.1) +
  #geom_smooth(method = "loess", span = 10/nrow(FLUO), method.args = list(degree = 1), aes(x = FLUO$t, y = FLUO$MEAN_GREEN.1.95), color = "navy") +
  coord_cartesian(xlim = c(0,t[length(t)]), ylim = c(1,Response_Stats_df[nrow(RESPONSE_CALLS),4])) +
  labs(title = "Mut DRG Responses to\n5mM \u03B2-alanine, 200\u03BCM AITC, 140mM K", x = "Time (s)", y = "340:380\nNormalized Fluorescence (a.u.)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank()) +
  geom_hline(yintercept = 1.0, color = "black", linetype = "solid")

CaImaging





B_Peak_Mags_WT = vector()
B_Peak_Mags_Mut = vector()
LY_Peak_Mags_WT = vector()
LY_Peak_Mags_Mut = vector()
CYM_Peak_Mags_WT = vector()
CYM_Peak_Mags_Mut = vector()
AITC_Peak_Mags_WT = vector()
AITC_Peak_Mags_Mut = vector()
Cap_Peak_Mags_WT = vector()
Cap_Peak_Mags_Mut = vector()

  ## Create a function that determines baseline ratiometric fluorescences before increases in fluorescence in "responder" ROIs
Find_Baselines = function(SmoothFLUO, agonist, response_frame){
  ## Determine baseline ratiometric fluorescences for smoothed fluorescence data frame
  temp = vector()
  Baseline_DF = matrix(nrow = (length(agonist[seq(2,length(agonist),2)]))-1, ncol = ncol(SmoothFLUO)-1)
  Baseline_DF = as.data.frame(Baseline_DF)
  colnames(Baseline_DF) = colnames(SmoothFLUO)[2:ncol(SmoothFLUO)]
  rownames(Baseline_DF) = rownames(response_frame)[1:((length(agonist[seq(2,length(agonist),2)]))-1)]
  
  for (i in 1:ncol(RESPONSE_CALLS)){
    for (j in 1:(nrow(RESPONSE_CALLS)-1)){
      if (RESPONSE_CALLS[j,i] == "Responder"){
        temp = diff(as.numeric(Slope_idx_DF[ , i]))
        Baseline_DF[j,i] = SmoothFLUO[ which(temp == 0)[which(which(temp == 0 ) %in% response_frame[j,1]:response_frame[j,2]) [1]]-1, i+1]
      }
    }
  }
  Baseline_DF <<- Baseline_DF
}
  
  ## Find magnitudes of ratiometric fluorescence increases in "responder" ROIs
Mag_Quant = function(Responders_List, Agonist_Number_From_response_frame, Agonist, Genotype){
  ## Function will return the baseline ratiometric and background-subtracted fluorescence values of a given list of regions of interest and any that undetectable. NAs represent ROIs that don't reach 2% slope increase or 10% FLUO increase
  ## Exported to the Global Environment are the net fluorescence magnitudes (relative to baseline)
  ## Arg "Genotype" is either "WT" or "Mut"
  ## Arg "Agonist" is "LY344864", "CYM5442", "Beta-alanine", "AITC", or "Capsaicin"
  ## Arg "Agonist_Number_From_response_frame" is 1 or 2 or 3 (in a couple cases)
  ## Arg "Responders_List" is a concatenated list of integers of the ROIs of a given experiment deemed to be "responders"
  if (Agonist == "LY344864"){
    if (Genotype == "WT"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
      ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        LY_Peak_Mags_WT = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        LY_Peak_Mags_WT = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      LY_Peak_Mags_WT <<- LY_Peak_Mags_WT
    }
    else if (Genotype == "Mut"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
      ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        LY_Peak_Mags_Mut = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        LY_Peak_Mags_Mut = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      LY_Peak_Mags_Mut <<- LY_Peak_Mags_Mut
    }
  }
  if (Agonist == "CYM5442"){
    if (Genotype == "WT"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
      ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        CYM_Peak_Mags_WT = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        CYM_Peak_Mags_WT = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      CYM_Peak_Mags_WT <<- CYM_Peak_Mags_WT
    }
    else if (Genotype == "Mut"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
      ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        CYM_Peak_Mags_Mut = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        CYM_Peak_Mags_Mut = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      CYM_Peak_Mags_Mut <<- CYM_Peak_Mags_Mut
    }
  }
  if (Agonist == "Beta-alanine"){
    if (Genotype == "WT"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
      ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        B_Peak_Mags_WT = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        B_Peak_Mags_WT = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      B_Peak_Mags_WT <<- B_Peak_Mags_WT
    }
    else if (Genotype == "Mut"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
      ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        B_Peak_Mags_Mut = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        B_Peak_Mags_Mut = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      B_Peak_Mags_Mut <<- B_Peak_Mags_Mut
    }
  }
  if (Agonist == "AITC"){
    if (Genotype == "WT"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
      ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        AITC_Peak_Mags_WT = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        AITC_Peak_Mags_WT = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      AITC_Peak_Mags_WT <<- AITC_Peak_Mags_WT
    }
    else if (Genotype == "Mut"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
      ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        AITC_Peak_Mags_Mut = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        AITC_Peak_Mags_Mut = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      AITC_Peak_Mags_Mut <<- AITC_Peak_Mags_Mut
    }
  }
  if (Agonist == "Capsaicin"){
    if (Genotype == "WT"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
      ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        Cap_Peak_Mags_WT = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        Cap_Peak_Mags_WT = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      Cap_Peak_Mags_WT <<- Cap_Peak_Mags_WT
    }
    else if (Genotype == "Mut"){
      values = vector()
      for (i in 1:length(c(Responders_List))){
        values[i] = max( SmoothFLUO[ response_frame[Agonist_Number_From_response_frame,1]:response_frame[Agonist_Number_From_response_frame,2], c(Responders_List)+1][i] )
      }
        ## Find the baselines for responding ROIs
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
        ## Find the responders that don't reach a 2% slope or 10% magnitude increase
      colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
        ## Find baselines of all responding ROIs that meet response criteria
      as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
        ## Find the magnitudes of responses (maxima of responders - baseline of those responders), accounting for responses that don't meet mag or slope criteria
      values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
        ## Find the net magnitudes of responses for responding ROIs
      values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      ## Rename the vector with the values, subtracted by the baseline
      if ( any(is.na(Baseline_DF[Agonist_Number_From_response_frame, c(Responders_List)])) == TRUE){
        Cap_Peak_Mags_Mut = values[which(!is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))] - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][-which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])
      } else {
        Cap_Peak_Mags_Mut = values - as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)])
      }
      Cap_Peak_Mags_Mut <<- Cap_Peak_Mags_Mut
    }
  }
  return( c(as.numeric(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]), colnames(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)][which(is.na(Baseline_DF[Agonist_Number_From_response_frame,c(Responders_List)]))])))
  
}

##### WT LY #####
  ## 10/30 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(17,28,36,43,54,61,65,68,87,91), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in A4 from 10/30
A4_10_30_LY_Mags = as.numeric(c(LY_Peak_Mags_WT, 0.2143427, 0.2324647, 0.1986785))
  ## Find the avg response magnitude
A4_10_30_LY_avg = sum(A4_10_30_LY_Mags) / length(A4_10_30_LY_Mags)

  ## 10/30 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(35,57,64,121), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in B2 from 10/30
B2_10_30_LY_Mags = LY_Peak_Mags_WT
  ## Find the avg response magnitude
B2_10_30_LY_avg = sum(B2_10_30_LY_Mags) / length(B2_10_30_LY_Mags)

  ## 10/30 B3
  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in B3 from 10/30
#B3_10_30_LY_Mags = LY_Peak_Mags_WT
  ## Find the avg response magnitude
#B3_10_30_LY_avg = sum(B3_10_30_LY_Mags) / length(B3_10_30_LY_Mags)

  ## 10/30 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(40), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in C1 from 10/30
  ## Manually enter #40
C1_10_30_LY_Mags = as.numeric(c(0.2395488))
  ## Find the avg response magnitude
C1_10_30_LY_avg = sum(C1_10_30_LY_Mags) / length(C1_10_30_LY_Mags)

  ## 11/8 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(13,14,38,42,84,102,114,139,160,165,196,226), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in A3 from 11/8
A3_11_8_LY_Mags = LY_Peak_Mags_WT
  ## Find the avg response magnitude
A3_11_8_LY_avg = sum(A3_11_8_LY_Mags) / length(A3_11_8_LY_Mags)

  ## 11/8 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(17), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in A4 from 11/8
  ## Manually enter #17
A4_11_8_LY_Mags = as.numeric(0.1028318)
  ## Find the avg response magnitude
A4_11_8_LY_avg = sum(A4_11_8_LY_Mags) / length(A4_11_8_LY_Mags)

  ## 11/8 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,2,4), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in B2 from 11/8
B2_11_8_LY_Mags = LY_Peak_Mags_WT
  ## Find the avg response magnitude
B2_11_8_LY_avg = sum(B2_11_8_LY_Mags) / length(B2_11_8_LY_Mags)

  ## 11/8 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(7,9), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in B4 from 11/8
B4_11_8_LY_Mags = LY_Peak_Mags_WT
  ## Find the avg response magnitude
B4_11_8_LY_avg = sum(B2_11_8_LY_Mags) / length(B2_11_8_LY_Mags)

  ## 11/15 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(17), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in A1 from 11/15
  ## Manually enter #17
A1_11_15_LY_Mags = as.numeric(1.1793108)
  ## Find the avg response magnitude
A1_11_15_LY_avg = sum(A1_11_15_LY_Mags) / length(A1_11_15_LY_Mags)

  ## 11/15 A3
  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in A3 from 11/15
#A3_11_15_LY_Mags = LY_Peak_Mags_WT
  ## Find the avg response magnitude
#A3_11_15_LY_avg = sum(A3_11_15_LY_Mags) / length(A3_11_15_LY_Mags)

  ## 11/15 B3
  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in B3 from 11/15
#B3_11_15_LY_Mags = LY_Peak_Mags_WT
  ## Find the avg response magnitude
#B3_11_15_LY_avg = sum(B3_11_15_LY_Mags) / length(B3_11_15_LY_Mags)

  ## 12/4 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(35), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in A3 from 12/4
  ## Manually enter #35
A3_12_4_LY_Mags = as.numeric(0.6717732)
  ## Find the avg response magnitude
A3_12_4_LY_avg = sum(A3_12_4_LY_Mags) / length(A3_12_4_LY_Mags)

  ## 12/4 A4
  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in A4 from 12/4
#A4_12_4_LY_Mags = LY_Peak_Mags_WT
## Find the avg response magnitude
#A4_12_4_LY_avg = sum(A4_12_4_LY_Mags) / length(A4_12_4_LY_Mags)

  ## 12/4 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(21), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in B1 from 12/4
B1_12_4_LY_Mags = as.numeric(0.2577467)
  ## Find the avg response magnitude
B1_12_4_LY_avg = sum(B1_12_4_LY_Mags) / length(B1_12_4_LY_Mags)

#  ## 12/4 B3
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
#  ## Store the response net magnitude increases in fluorescence in response to LY in B3 from 12/4
#B3_12_4_LY_Mags = LY_Peak_Mags_WT
#  ## Find the avg response magnitude
#B3_12_4_LY_avg = sum(B3_12_4_LY_Mags) / length(B3_12_4_LY_Mags)

  ## 12/4 B3_2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in B3_2 from 12/4
B3_2_12_4_LY_Mags = LY_Peak_Mags_WT
  ## Find the avg response magnitude
B3_2_12_4_LY_avg = sum(B3_2_12_4_LY_Mags) / length(B3_2_12_4_LY_Mags)

  ## 12/11 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(25,95,99), Agonist_Number_From_response_frame = 2, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in A1 from 12/11
A1_12_11_LY_Mags = as.numeric(c(LY_Peak_Mags_WT, 0.1182547))
  ## Find the avg response magnitude
A1_12_11_LY_avg = sum(A1_12_11_LY_Mags) / length(A1_12_11_LY_Mags)

  ## 12/11 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(41,90), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in A3 from 12/11
  ## Manually enter #41, #90
A3_12_11_LY_Mags = as.numeric(c(0.1731586, 0.1979932))
  ## Find the avg response magnitude
A3_12_11_LY_avg = sum(A3_12_11_LY_Mags) / length(A3_12_11_LY_Mags)

  ## 12/11 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(34,43,73,112), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in B2 from 12/11
B2_12_11_LY_Mags = as.numeric(c(LY_Peak_Mags_WT, 0.104087))
  ## Find the avg response magnitude
B2_12_11_LY_avg = sum(B2_12_11_LY_Mags) / length(B2_12_11_LY_Mags)

  ## 12/11 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(29,43), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to LY in B4 from 12/11
B4_12_11_LY_Mags = LY_Peak_Mags_WT
  ## Find the avg response magnitude
B4_12_11_LY_avg = sum(B4_12_11_LY_Mags) / length(B4_12_11_LY_Mags)

##### Mut LY #####
  ## 11/6 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(10,24), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to LY in A4 from 11/6
A4_11_6_LY_Mags = LY_Peak_Mags_Mut
  ## Find the avg response magnitude
A4_11_6_LY_avg = sum(A4_11_6_LY_Mags) / length(A4_11_6_LY_Mags)

  ## 11/6 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(21,23,50,66), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to LY in B1 from 11/6
B1_11_6_LY_Mags = LY_Peak_Mags_Mut
  ## Find the avg response magnitude
B1_11_6_LY_avg = sum(B1_11_6_LY_Mags) / length(B1_11_6_LY_Mags)

#  ## 11/6 B2
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to LY in B2 from 11/6
#B2_11_6_LY_Mags = LY_Peak_Mags_Mut
#  ## Find the avg response magnitude
#B2_11_6_LY_avg = sum(B2_11_6_LY_Mags) / length(B2_11_6_LY_Mags)

#  ## 11/6 B3
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to LY in B3 from 11/6
#B3_11_6_LY_Mags = LY_Peak_Mags_Mut
#  ## Find the avg response magnitude
#B3_11_6_LY_avg = sum(B3_11_6_LY_Mags) / length(B3_11_6_LY_Mags)

#  ## 11/6 B4
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to LY in B4 from 11/6
#B4_11_6_LY_Mags = LY_Peak_Mags_Mut
#  ## Find the avg response magnitude
#B4_11_6_LY_avg = sum(B4_11_6_LY_Mags) / length(B4_11_6_LY_Mags)

  ## 11/11 A1
  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to LY in B1 from 11/11
  ## Manually add #33,34
A1_11_11_LY_Mags = as.numeric(c( 0.2795854, 0.3046652))
  ## Find the avg response magnitude
A1_11_11_LY_avg = sum(A1_11_11_LY_Mags) / length(A1_11_11_LY_Mags)

#  ## 11/11 A4
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to LY in A4 from 11/11
#A4_11_11_LY_Mags = LY_Peak_Mags_Mut
#  ## Find the avg response magnitude
#A4_11_11_LY_avg = sum(A4_11_11_LY_Mags) / length(A4_11_11_LY_Mags)

  ## 11/11 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(70,80,85,140,144,145), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to LY in B2 from 11/11
B2_11_11_LY_Mags = LY_Peak_Mags_Mut
  ## Find the avg response magnitude
B2_11_11_LY_avg = sum(B2_11_11_LY_Mags) / length(B2_11_11_LY_Mags)

  ## 11/11 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(32), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to LY in B4 from 11/11
  ## Manually add 32
B4_11_11_LY_Mags = as.numeric(0.6752444)
  ## Find the avg response magnitude
B4_11_11_LY_avg = sum(B4_11_11_LY_Mags) / length(B4_11_11_LY_Mags)

  ## 11/20 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(10,13), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to LY in A4 from 11/20
A4_11_20_LY_Mags = LY_Peak_Mags_Mut
  ## Find the avg response magnitude
A4_11_20_LY_avg = sum(A4_11_20_LY_Mags) / length(A4_11_20_LY_Mags)

#  ## 11/20 B2
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to LY in B2 from 11/20
#B2_11_20_LY_Mags = LY_Peak_Mags_Mut
#  ## Find the avg response magnitude
#B2_11_20_LY_avg = sum(B2_11_20_LY_Mags) / length(B2_11_20_LY_Mags)

  ## 11/20 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(35,62,69,77,85), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to LY in B4 from 11/20
B4_11_20_LY_Mags = LY_Peak_Mags_Mut
  ## Find the avg response magnitude
B4_11_20_LY_avg = sum(B4_11_20_LY_Mags) / length(B4_11_20_LY_Mags)

#  ## 12/6 A4
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to LY in A4 from 12/6
#A4_12_6_LY_Mags = LY_Peak_Mags_Mut
#  ## Find the avg response magnitude
#A4_12_6_LY_avg = sum(A4_12_6_LY_Mags) / length(A4_12_6_LY_Mags)

  ## 12/6 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(45,66,68), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to LY in B2 from 12/6
B2_12_6_LY_Mags = LY_Peak_Mags_Mut
  ## Find the avg response magnitude
B2_12_6_LY_avg = sum(B2_12_6_LY_Mags) / length(B2_12_6_LY_Mags)

#  ## 12/13 B1
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to LY in B1 from 12/13
#B1_12_13_LY_Mags = LY_Peak_Mags_Mut
#  ## Find the avg response magnitude
#B1_12_13_LY_avg = sum(B1_12_13_LY_Mags) / length(B1_12_13_LY_Mags)

  ## 12/13 C2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(8,51,91), Agonist_Number_From_response_frame = 1, Agonist = "LY344864", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to LY in C2 from 12/13
C2_12_13_LY_Mags = LY_Peak_Mags_Mut
  ## Find the avg response magnitude
C2_12_13_LY_avg = sum(C2_12_13_LY_Mags) / length(C2_12_13_LY_Mags)

##### WT CYM #####
  ## 10/30 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(31,40,64:66,68,81,87,90), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in A4 from 10/30
A4_10_30_CYM_Mags = CYM_Peak_Mags_WT
  ## Find the avg response magnitude
A4_10_30_CYM_avg = sum(A4_10_30_CYM_Mags) / length(A4_10_30_CYM_Mags)

  ## 10/30 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(57,151), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B1 from 10/30
B1_10_30_CYM_Mags = CYM_Peak_Mags_WT
  ## Find the avg response magnitude
B1_10_30_CYM_avg = sum(B1_10_30_CYM_Mags) / length(B1_10_30_CYM_Mags)

  ## 10/30 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(57,121), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B2 from 10/30
B2_10_30_CYM_Mags = CYM_Peak_Mags_WT
  ## Find the avg response magnitude
B2_10_30_CYM_avg = sum(B2_10_30_CYM_Mags) / length(B2_10_30_CYM_Mags)

  ## 10/30 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,9,39), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B3 from 10/30
B3_10_30_CYM_Mags = CYM_Peak_Mags_WT
  ## Find the avg response magnitude
B3_10_30_CYM_avg = sum(B3_10_30_CYM_Mags) / length(B3_10_30_CYM_Mags)

  ## 10/30 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(32), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in C1 from 10/30
C1_10_30_CYM_Mags = as.numeric(0.3335494)
  ## Find the avg response magnitude
C1_10_30_CYM_avg = sum(C1_10_30_CYM_Mags) / length(C1_10_30_CYM_Mags)

  ## 11/8 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(15,21,34,43,81,84,94,103,107,114,126,129,139,152,160,165,227), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in A3 from 11/8
A3_11_8_CYM_Mags = CYM_Peak_Mags_WT
  ## Find the avg response magnitude
A3_11_8_CYM_avg = sum(A3_11_8_CYM_Mags) / length(A3_11_8_CYM_Mags)

  ## 11/8 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = as.numeric(17), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in A4 from 11/8
A4_11_8_CYM_Mags = as.numeric(0.3953106)
  ## Find the avg response magnitude
A4_11_8_CYM_avg = sum(A4_11_8_CYM_Mags) / length(A4_11_8_CYM_Mags)

  ## 11/8 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(15,16,26,37,65,66,92), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B3 from 11/8
B3_11_8_CYM_Mags = as.numeric(c(CYM_Peak_Mags_WT, 0.1718266))
  ## Find the avg response magnitude
B3_11_8_CYM_avg = sum(B3_11_8_CYM_Mags) / length(B3_11_8_CYM_Mags)

  ## 11/8 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(5,7,9,19), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B4 from 11/8
B4_11_8_CYM_Mags = CYM_Peak_Mags_WT
  ## Find the avg response magnitude
B4_11_8_CYM_avg = sum(B4_11_8_CYM_Mags) / length(B4_11_8_CYM_Mags)

  ## 11/15 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,19), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in A1 from 11/15
A1_11_15_CYM_Mags = CYM_Peak_Mags_WT
  ## Find the avg response magnitude
A1_11_15_CYM_avg = sum(A1_11_15_CYM_Mags) / length(A1_11_15_CYM_Mags)

  ## 11/15 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(5), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in A3 from 11/15
A3_11_15_CYM_Mags = as.numeric(0.3079227)
  ## Find the avg response magnitude
A3_11_15_CYM_avg = sum(A3_11_15_CYM_Mags) / length(A3_11_15_CYM_Mags)

  ## 11/15 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(40,42), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B2 from 11/15
B2_11_15_CYM_Mags = CYM_Peak_Mags_WT
  ## Find the avg response magnitude
B2_11_15_CYM_avg = sum(B2_11_15_CYM_Mags) / length(B2_11_15_CYM_Mags)

  ## 11/15 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1), Agonist_Number_From_response_frame = 3, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B3 from 11/15
B3_11_15_CYM_Mags = as.numeric(0.1854545)
  ## Find the avg response magnitude
B3_11_15_CYM_avg = sum(B3_11_15_CYM_Mags) / length(B3_11_15_CYM_Mags)

  ## 12/4 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,44), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in A4 from 12/4
A4_12_4_CYM_Mags = as.numeric(c(CYM_Peak_Mags_WT, 0.1027383))
  ## Find the avg response magnitude
A4_12_4_CYM_avg = sum(A4_12_4_CYM_Mags) / length(A4_12_4_CYM_Mags)

  ## 12/4 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(9,11,12), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to CYM in C1 from 12/4
C1_12_4_CYM_Mags = CYM_Peak_Mags_WT
  ## Find the avg response magnitude
C1_12_4_CYM_avg = sum(C1_12_4_CYM_Mags) / length(C1_12_4_CYM_Mags)

##### Mut CYM #####
  ## 11/6 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,10,36), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B1 from 11/6
B1_11_6_CYM_Mags = CYM_Peak_Mags_Mut
  ## Find the avg response magnitude
B1_11_6_CYM_avg = sum(B1_11_6_CYM_Mags) / length(B1_11_6_CYM_Mags)

  ## 11/6 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,9), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B2 from 11/6
B2_11_6_CYM_Mags = CYM_Peak_Mags_Mut
  ## Find the avg response magnitude
B2_11_6_CYM_avg = sum(B2_11_6_CYM_Mags) / length(B2_11_6_CYM_Mags)

  ## 11/6 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(17,29), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B3 from 11/6
B3_11_6_CYM_Mags = as.numeric(c(CYM_Peak_Mags_Mut, 0.1703919))
  ## Find the avg response magnitude
B3_11_6_CYM_avg = sum(B3_11_6_CYM_Mags) / length(B3_11_6_CYM_Mags)

#  ## 11/11 A4
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to CYM in A4 from 11/11
#A4_11_11_CYM_Mags = CYM_Peak_Mags_Mut
#  ## Find the avg response magnitude
#A4_11_11_CYM_avg = sum(A4_11_11_CYM_Mags) / length(A4_11_11_CYM_Mags)

  ## 11/11 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(15,18,42,58,102,131), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B1 from 11/11
B1_11_11_CYM_Mags = as.numeric(c(CYM_Peak_Mags_Mut, 0.2068902, 0.2661255))
  ## Find the avg response magnitude
B1_11_11_CYM_avg = sum(B1_11_11_CYM_Mags) / length(B1_11_11_CYM_Mags)

  ## 11/11 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(20,37,96,145), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B2 from 11/11
B2_11_11_CYM_Mags = CYM_Peak_Mags_Mut
  ## Find the avg response magnitude
B2_11_11_CYM_avg = sum(B2_11_11_CYM_Mags) / length(B2_11_11_CYM_Mags)

  ## 11/11 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(37,94), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B4 from 11/11
B4_11_11_CYM_Mags = CYM_Peak_Mags_Mut
  ## Find the avg response magnitude
B4_11_11_CYM_avg = sum(B4_11_11_CYM_Mags) / length(B4_11_11_CYM_Mags)

#  ## 11/20 A4
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to CYM in A4 from 11/20
#A4_11_20_CYM_Mags = CYM_Peak_Mags_Mut
#  ## Find the avg response magnitude
#A4_11_20_CYM_avg = sum(A4_11_20_CYM_Mags) / length(A4_11_20_CYM_Mags)

#  ## 11/20 B1
#  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
#  ## Find response magnitudes
#Mag_Quant(Responders_List = c(16,32), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "Mut")
#  ## Store the response net magnitude increases in fluorescence in response to CYM in B1 from 11/20
#B1_11_20_CYM_Mags = CYM_Peak_Mags_Mut
#  ## Find the avg response magnitude
#B1_11_20_CYM_avg = sum(B1_11_20_CYM_Mags) / length(B1_11_20_CYM_Mags)

  ## 11/20 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(16,32), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in C1 from 11/20
C1_11_20_CYM_Mags = CYM_Peak_Mags_Mut
  ## Find the avg response magnitude
C1_11_20_CYM_avg = sum(C1_11_20_CYM_Mags) / length(C1_11_20_CYM_Mags)

  ## 12/6 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in A1 from 12/6
  ## Add #4 manually
A1_12_6_CYM_Mags = as.numeric(0.3401297)
  ## Find the avg response magnitude
A1_12_6_CYM_avg = sum(A1_12_6_CYM_Mags) / length(A1_12_6_CYM_Mags)

  ## 12/6 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(6,21), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in A3 from 12/6
  ## Add #6,21 manually
A3_12_6_CYM_Mags = as.numeric(c(0.5920739, 0.5334424))
  ## Find the avg response magnitude
A3_12_6_CYM_avg = sum(A3_12_6_CYM_Mags) / length(A3_12_6_CYM_Mags)

  ## 12/6 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(5), Agonist_Number_From_response_frame = 2, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in A4 from 12/6
  ## Add #5 manually
A4_12_6_CYM_Mags = as.numeric(0.4747043)
  ## Find the avg response magnitude
A4_12_6_CYM_avg = sum(A4_12_6_CYM_Mags) / length(A4_12_6_CYM_Mags)

  ## 12/6 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(20,23,27,33,67), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in B1 from 12/6
B1_12_6_CYM_Mags = CYM_Peak_Mags_Mut
  ## Find the avg response magnitude
B1_12_6_CYM_avg = sum(B1_12_6_CYM_Mags) / length(B1_12_6_CYM_Mags)

  ## 12/6 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(28,33,55,57), Agonist_Number_From_response_frame = 1, Agonist = "CYM5442", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to CYM in C1 from 12/6
C1_12_6_CYM_Mags = CYM_Peak_Mags_Mut
  ## Find the avg response magnitude
C1_12_6_CYM_avg = sum(C1_12_6_CYM_Mags) / length(C1_12_6_CYM_Mags)

##### WT B #####
  ## 11/15 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(5,8:10,15:18,21,23,34,36,40,43), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in A4 from 11/15
A4_11_15_B_Mags = B_Peak_Mags_WT
  ## Find the avg response magnitude
A4_11_15_B_avg = sum(A4_11_15_B_Mags) / length(A4_11_15_B_Mags)

  ## 11/15 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,4,20,21,31,39), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B4 from 11/15
B4_11_15_B_Mags = as.numeric(c(B_Peak_Mags_WT, 0.09809681))
  ## Find the avg response magnitude
B4_11_15_B_avg = sum(B4_11_15_B_Mags) / length(B4_11_15_B_Mags)

  ## 12/4 A1
  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in A1 from 12/4
#A1_12_4_B_Mags = B_Peak_Mags_WT
  ## Find the avg response magnitude
#A1_12_4_B_avg = sum(A1_12_4_B_Mags) / length(A1_12_4_B_Mags)

  ## 12/4 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,3,6,9,18,19), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B1 from 12/4
B1_12_4_B_Mags = B_Peak_Mags_WT
  ## Find the avg response magnitude
B1_12_4_B_avg = sum(B1_12_4_B_Mags) / length(B1_12_4_B_Mags)

  ## 12/4 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1:3,5,10,17,18,20,26,27,29,32,35,40:43), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B2 from 12/4
B2_12_4_B_Mags = B_Peak_Mags_WT
  ## Find the avg response magnitude
B2_12_4_B_avg = sum(B2_12_4_B_Mags) / length(B2_12_4_B_Mags)

  ## 12/4 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,2,5,14,19,21,22,24,25,27,28,31,34,39,41,50,52,53,65,66,70,71,79,90,93), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B4 from 12/4
B4_12_4_B_Mags = as.numeric(c(B_Peak_Mags_WT, 0.1498191, 0.103728))
  ## Find the avg response magnitude
B4_12_4_B_avg = sum(B4_12_4_B_Mags) / length(B4_12_4_B_Mags)

  ## 12/4 C2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,7,15,28,34,39,58,62,64), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in C2 from 12/4
C2_12_4_B_Mags = as.numeric(c(B_Peak_Mags_WT, 0.3026034))
  ## Find the avg response magnitude
C2_12_4_B_avg = sum(C2_12_4_B_Mags) / length(C2_12_4_B_Mags)

  ## 12/11 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,4,8,19,22,27,48,56,57,77,81,92,94,99,103,115,118), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in A1 from 12/11
A1_12_11_B_Mags = as.numeric(c(B_Peak_Mags_WT, 0.1424865))
  ## Find the avg response magnitude
A1_12_11_B_avg = sum(A1_12_11_B_Mags) / length(A1_12_11_B_Mags)

  ## 12/11 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,3,10,11,16,22,25,28,38,44,46,47,55,56,66,69,72,77,84,94,97,98,112,132,148,157,157,169,170,185,197,209,210,215,236,237,241), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B1 from 12/11
B1_12_11_B_Mags = as.numeric(c(B_Peak_Mags_WT, 0.1020502))
  ## Find the avg response magnitude
B1_12_11_B_avg = sum(B1_12_11_B_Mags) / length(B1_12_11_B_Mags)

  ## 12/11 B3
## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
## Find response magnitudes
Mag_Quant(Responders_List = c(16,35,38,43,46,50,53,54,62,65,67,78), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
## Store the response net magnitude increases in fluorescence in response to Beta in B3 from 12/11
B3_12_11_B_Mags = as.numeric(c(B_Peak_Mags_WT, 0.1222118))
## Find the avg response magnitude
B3_12_11_B_avg = sum(B3_12_11_B_Mags) / length(B3_12_11_B_Mags)

  ## 12/11 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(9,18,19,28,30,43,44,69,74,76,81,85,88,91,97,109,118,143,149,173,181), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in C1 from 12/11
C1_12_11_B_Mags = as.numeric(c(B_Peak_Mags_WT, 0.2257973))
  ## Find the avg response magnitude
C1_12_11_B_avg = sum(C1_12_11_B_Mags) / length(C1_12_11_B_Mags)

  ## 12/11 C2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(28,30), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in C1 from 12/11
C2_12_11_B_Mags = B_Peak_Mags_WT
  ## Find the avg response magnitude
C2_12_11_B_avg = sum(C2_12_11_B_Mags) / length(C2_12_11_B_Mags)

  ## 12/11 C3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(53,71,131), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Beta in C3 from 12/11
C3_12_11_B_Mags = B_Peak_Mags_WT
  ## Find the avg response magnitude
C3_12_11_B_avg = sum(C3_12_11_B_Mags) / length(C3_12_11_B_Mags)

##### Mut B #####
  ## 11/11 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,6,22,23,30,37,38,85,90,129,131), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in A3 from 11/11
A3_11_11_B_Mags = as.numeric(c(B_Peak_Mags_Mut, 0.449383, 0.9216661, 0.2235338, 0.5278941, 0.6031466, 0.2686332, 0.2896009, 0.4520882, 0.1837269))
  ## Find the avg response magnitude
A3_11_11_B_avg = sum(C3_12_11_B_Mags) / length(C3_12_11_B_Mags)

  ## 11/11 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,5,27,30,33,55,57,59,66,73,78:81,88,92,94,103,104,113,127,138,140,142,148,154,156,158,167,187,190,194,208,212,215,216,228,229), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B3 from 11/11
B3_11_11_B_Mags = as.numeric(c(B_Peak_Mags_Mut, 0.0963641, 0.1208852, 0.1876179))
  ## Find the avg response magnitude
B3_11_11_B_avg = sum(B3_11_11_B_Mags) / length(B3_11_11_B_Mags)

  ## 11/20 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(11,28), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in A1 from 11/20
A1_11_20_B_Mags = B_Peak_Mags_Mut
  ## Find the avg response magnitude
A1_11_20_B_avg = sum(A1_11_20_B_Mags) / length(A1_11_20_B_Mags)

  ## 11/20 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(6,12,13,15,22,26,38,49), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in A3 from 11/20
A3_11_20_B_Mags = B_Peak_Mags_Mut
  ## Find the avg response magnitude
A3_11_20_B_avg = sum(A3_11_20_B_Mags) / length(A3_11_20_B_Mags)

  ## 11/20 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,3,16,18,22,23,24,27,28,31), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B2 from 11/20
B2_11_20_B_Mags = B_Peak_Mags_Mut
  ## Find the avg response magnitude
B2_11_20_B_avg = sum(B2_11_20_B_Mags) / length(B2_11_20_B_Mags)

  ## 11/20 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(7,15,18,20,21,23,25,29,31,34,36,48,50,52,61,66,68,71,75,79,81), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B3 from 11/20
B3_11_20_B_Mags = as.numeric(c(B_Peak_Mags_Mut, 0.2101837))
  ## Find the avg response magnitude
B3_11_20_B_avg = sum(B3_11_20_B_Mags) / length(B3_11_20_B_Mags)

  ## 12/6 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,8,15,16,30,33,35), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B3 from 12/6
B3_12_6_B_Mags = B_Peak_Mags_Mut
  ## Find the avg response magnitude
B3_12_6_B_avg = sum(B3_12_6_B_Mags) / length(B3_12_6_B_Mags)

  ## 12/6 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,3,10,12,25,29,40,49,52), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B4 from 12/6
B4_12_6_B_Mags = B_Peak_Mags_Mut
  ## Find the avg response magnitude
B4_12_6_B_avg = sum(B4_12_6_B_Mags) / length(B4_12_6_B_Mags)

  ## 12/6 C2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,13,19,26,27,32,36,42,45,65,69,70,73,76:78,82,85,97,126,129,131), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in C2 from 12/6
C2_12_6_B_Mags = B_Peak_Mags_Mut
  ## Find the avg response magnitude
C2_12_6_B_avg = sum(C2_12_6_B_Mags) / length(C2_12_6_B_Mags)

  ## 12/6 C3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(26,49,74,164,227), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in C3 from 12/6
C3_12_6_B_Mags = as.numeric(c(B_Peak_Mags_Mut, 0.1895349))
  ## Find the avg response magnitude
C3_12_6_B_avg = sum(C3_12_6_B_Mags) / length(C3_12_6_B_Mags)

  ## 12/13 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(6,36,37,42,43,49,53:55,60,62,73), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in A3 from 12/13
A3_12_13_B_Mags = as.numeric(c(B_Peak_Mags_Mut, 0.2785107))
  ## Find the avg response magnitude
A3_12_13_B_avg = sum(A3_12_13_B_Mags) / length(A3_12_13_B_Mags)

  ## 12/13 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,8,17,26,32), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in A4 from 12/13
A4_12_13_B_Mags = as.numeric(c(B_Peak_Mags_Mut, 0.2934113))
  ## Find the avg response magnitude
A4_12_13_B_avg = sum(A4_12_13_B_Mags) / length(A4_12_13_B_Mags)

  ## 12/13 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,13,14,16,17,24,34,39), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B2 from 12/13
B2_12_13_B_Mags = B_Peak_Mags_Mut
  ## Find the avg response magnitude
B2_12_13_B_avg = sum(B2_12_13_B_Mags) / length(B2_12_13_B_Mags)

  ## 12/13 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,3,7,8,11,22,33,43,44,51,58,59,62,70,76,86,87), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B3 from 12/13
B3_12_13_B_Mags = as.numeric(c(B_Peak_Mags_Mut, 0.1292285))
  ## Find the avg response magnitude
B3_12_13_B_avg = sum(B3_12_13_B_Mags) / length(B3_12_13_B_Mags)

  ## 12/13 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(10,20,26,29,43,55,58,61,63,64,66:68,72,76,78,79,82,93,100:102,109,110), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in B4 from 12/13
B4_12_13_B_Mags = B_Peak_Mags_Mut
  ## Find the avg response magnitude
B4_12_13_B_avg = sum(B4_12_13_B_Mags) / length(B4_12_13_B_Mags)

  ## 12/13 C3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(20,21,28,31,35,38,59,67,76,82,87,90,95), Agonist_Number_From_response_frame = 1, Agonist = "Beta-alanine", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Beta in C3 from 12/13
C3_12_13_B_Mags = as.numeric(c(B_Peak_Mags_Mut, 0.1800381))
  ## Find the avg response magnitude
C3_12_13_B_avg = sum(C3_12_13_B_Mags) / length(C3_12_13_B_Mags)

##### WT AITC #####
  ## 10/30 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(37,40,45,50,68,74), Agonist_Number_From_response_frame = 3, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A4 from 10/30
A4_10_30_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.1579174))
  ## Find the avg response magnitude
A4_10_30_AITC_avg = sum(A4_10_30_AITC_Mags) / length(A4_10_30_AITC_Mags)

  ## 10/30 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(6,9,15,19,28,30), Agonist_Number_From_response_frame = 1, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B4 from 10/30
B4_10_30_AITC_Mags = AITC_Peak_Mags_WT
  ## Find the avg response magnitude
B4_10_30_AITC_avg = sum(A4_10_30_AITC_Mags) / length(A4_10_30_AITC_Mags)

  ## 10/30 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(10,12,17,33,49,66), Agonist_Number_From_response_frame = 3, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C1 from 10/30
C1_10_30_AITC_Mags = AITC_Peak_Mags_WT
  ## Find the avg response magnitude
C1_10_30_AITC_avg = sum(C1_10_30_AITC_Mags) / length(C1_10_30_AITC_Mags)

  ## 11/8 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,4,22,28,34,36,42,54,58,68,82), Agonist_Number_From_response_frame = 1, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A1 from 11/8
A1_11_8_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.2264347))
  ## Find the avg response magnitude
A1_11_8_AITC_avg = sum(A1_11_8_AITC_Mags) / length(A1_11_8_AITC_Mags)

  ## 11/8 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,4,6,7,17,19,41), Agonist_Number_From_response_frame = 1, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B1 from 11/8
B1_11_8_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.1850825))
  ## Find the avg response magnitude
B1_11_8_AITC_avg = sum(B1_11_8_AITC_Mags) / length(B1_11_8_AITC_Mags)

  ## 11/8 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,3,6,9,15,18,19,23,26,27,30,32,33,57,64:66,73,77,94,103,106,107,109), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B3 from 11/8
B3_11_8_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.2058071, 0.2555837))
  ## Find the avg response magnitude
B3_11_8_AITC_avg = sum(B3_11_8_AITC_Mags) / length(B3_11_8_AITC_Mags)

  ## 11/15 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(16:18,21,23,26,29,34,40,43), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A4 from 11/15
A4_11_15_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.3709418, 0.1953214))
  ## Find the avg response magnitude
A4_11_15_AITC_avg = sum(A4_11_15_AITC_Mags) / length(A4_11_15_AITC_Mags)

  ## 11/15 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,5,8,11,15,21,23,26,27,28,31,35,36,42), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B2 from 11/15
B2_11_15_AITC_Mags = AITC_Peak_Mags_WT
  ## Find the avg response magnitude
B2_11_15_AITC_avg = sum(B2_11_15_AITC_Mags) / length(B2_11_15_AITC_Mags)

  ## 11/15 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,20,34), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B4 from 11/15
B4_11_15_AITC_Mags = AITC_Peak_Mags_WT
  ## Find the avg response magnitude
B4_11_15_AITC_avg = sum(B4_11_15_AITC_Mags) / length(B4_11_15_AITC_Mags)

  ## 12/4 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,7,8,20,32,33,35,37), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B2 from 12/4
B2_12_4_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.2847332))
  ## Find the avg response magnitude
B2_12_4_AITC_avg = sum(B2_12_4_AITC_Mags) / length(B2_12_4_AITC_Mags)

  ## 12/4 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(12,24,27,28,34,36,56,72,74,88,90), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B4 from 12/4
B4_12_4_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.4823356))
  ## Find the avg response magnitude
B4_12_4_AITC_avg = sum(B4_12_4_AITC_Mags) / length(B4_12_4_AITC_Mags)

  ## 12/4 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,11,18,19), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C1 from 12/4
C1_12_4_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.5817458))
  ## Find the avg response magnitude
C1_12_4_AITC_avg = sum(C1_12_4_AITC_Mags) / length(C1_12_4_AITC_Mags)

  ## 12/4 C2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(19,21), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C2 from 12/4
C2_12_4_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.3173787, 0.4101486))
  ## Find the avg response magnitude
C2_12_4_AITC_avg = sum(C2_12_4_AITC_Mags) / length(C2_12_4_AITC_Mags)

  ## 12/11 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2:4,7,8,10:13,15,17,19,21,23,25,32,34,35,37,39,47,53,54,57,60,61,65,66,68,74,77,79,81,84), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A3 from 12/11
A3_12_11_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.3136015))
  ## Find the avg response magnitude
A3_12_11_AITC_avg = sum(A3_12_11_AITC_Mags) / length(A3_12_11_AITC_Mags)

  ## 12/11 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1:3,6,12,16,19,20,23,24,28,32:39,41,43,44,46,50,53:55,57,59,61,62,65,67,69,70,75,78,80,87,91,100), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B3 from 12/11
B3_12_11_AITC_Mags = AITC_Peak_Mags_WT
  ## Find the avg response magnitude
B3_12_11_AITC_avg = sum(B3_12_11_AITC_Mags) / length(B3_12_11_AITC_Mags)

  ## 12/11 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2:4,6:10,12,18,19,21,25,28,35,43,44,46,48:50,54,56,58,60,65,66,71,73:77,80:89,91,92,97,100,101,103,107:110,113:116,118,120:124,132:135,137,140,143,146,147,149,151,154,156,163,173:177,181,183:187,192,193,195:199,204), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C1 from 12/11
C1_12_11_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.1629207))
  ## Find the avg response magnitude
C1_12_11_AITC_avg = sum(C1_12_11_AITC_Mags) / length(C1_12_11_AITC_Mags)

  ## 12/11 C2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,5:7,11,13,14,22,23,26,28:30,34,35,40:42,45,47,49,50,55:57,59,62:69,71,75,76,78,91,98,102), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C2 from 12/11
C2_12_11_AITC_Mags = as.numeric(c(AITC_Peak_Mags_WT, 0.1125343))
  ## Find the avg response magnitude
C2_12_11_AITC_avg = sum(C2_12_11_AITC_Mags) / length(C2_12_11_AITC_Mags)


##### Mut AITC #####
  ## 11/6 A1
  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 1, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A1 from 11/6
A1_11_6_AITC_Mags = as.numeric(0.3673884)
  ## Find the avg response magnitude
A1_11_6_AITC_avg = sum(A1_11_6_AITC_Mags) / length(A1_11_6_AITC_Mags)

  ## 11/6 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,23,32), Agonist_Number_From_response_frame = 1, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A3 from 11/6
A3_11_6_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
A3_11_6_AITC_avg = sum(A3_11_6_AITC_Mags) / length(A3_11_6_AITC_Mags)

  ## 11/20 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(10,15,24,25,27,28,40), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A1 from 11/20
A1_11_20_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
A1_11_20_AITC_avg = sum(A1_11_20_AITC_Mags) / length(A1_11_20_AITC_Mags)

  ## 11/20 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,6,12:15,25,38,39,43,45,48), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A3 from 11/20
A3_11_20_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
A3_11_20_AITC_avg = sum(A3_11_20_AITC_Mags) / length(A3_11_20_AITC_Mags)

  ## 11/20 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,5,7:10,14:17,26,28:31,35,40,42), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B1 from 11/20
B1_11_20_AITC_Mags = as.numeric(c(AITC_Peak_Mags_Mut, 0.2065447))
  ## Find the avg response magnitude
B1_11_20_AITC_avg = sum(B1_11_20_AITC_Mags) / length(B1_11_20_AITC_Mags)

  ## 11/20 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(6,11,15,20,21,25,26,28,31,34,36,43,49,50,52,58,65,68,69,71,75,81), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B3 from 11/20
B3_11_20_AITC_Mags = as.numeric(c(AITC_Peak_Mags_Mut, 0.2276222, 0.2188026))
  ## Find the avg response magnitude
B3_11_20_AITC_avg = sum(B3_11_20_AITC_Mags) / length(B3_11_20_AITC_Mags)

  ## 11/20 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1:3,9,11,12,15,18,21:23,25,29,30,32:36,40:42,44,47:49,51:53,55,58,59,66,67,69,75,80,81,83:88,90,92), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B4 from 11/20
B4_11_20_AITC_Mags = as.numeric(c(AITC_Peak_Mags_Mut, 0.1831757))
  ## Find the avg response magnitude
B4_11_20_AITC_avg = sum(B4_11_20_AITC_Mags) / length(B4_11_20_AITC_Mags)

  ## 11/20 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(5,7,9,12,14,15,18,27,29), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C1 from 11/20
C1_11_20_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
C1_11_20_AITC_avg = sum(C1_11_20_AITC_Mags) / length(C1_11_20_AITC_Mags)

  ## 12/6 A1
  ## Find baselines
#Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
#Mag_Quant(Responders_List = c(), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A1 from 12/6
#A1_12_6_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
#A1_12_6_AITC_avg = sum(A1_12_6_AITC_Mags) / length(A1_12_6_AITC_Mags)

  ## 12/6 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(5,6,8,21,22,31), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A3 from 12/6
A3_12_6_AITC_Mags = as.numeric(c(AITC_Peak_Mags_Mut, 0.1303312, 0.2332332))
  ## Find the avg response magnitude
A3_12_6_AITC_avg = sum(A3_12_6_AITC_Mags) / length(A3_12_6_AITC_Mags)

  ## 12/6 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(30,43,49,51,54,60), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B2 from 12/6
B2_12_6_AITC_Mags = as.numeric(c(AITC_Peak_Mags_Mut, 0.1673817))
  ## Find the avg response magnitude
B2_12_6_AITC_avg = sum(B2_12_6_AITC_Mags) / length(B2_12_6_AITC_Mags)

  ## 12/6 C1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,6,8,12,20,22,24,27,28,43,46,47,49,52,54,57), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C1 from 12/6
C1_12_6_AITC_Mags = as.numeric(c(AITC_Peak_Mags_Mut, 0.09975966))
  ## Find the avg response magnitude
C1_12_6_AITC_avg = sum(C1_12_6_AITC_Mags) / length(C1_12_6_AITC_Mags)

  ## 12/6 C2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,6,13,19,26,27,32,40,66,78,82,112,127), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C2 from 12/6
C2_12_6_AITC_Mags = as.numeric(c(AITC_Peak_Mags_Mut, 0.3231445, 0.3659826, 0.3037502))
  ## Find the avg response magnitude
C2_12_6_AITC_avg = sum(C2_12_6_AITC_Mags) / length(C2_12_6_AITC_Mags)

  ## 12/6 C3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(52,60,93,130,143,160,164,167,180,183,210,227,244,245,297), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C3 from 12/6
C3_12_6_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
C3_12_6_AITC_avg = sum(C3_12_6_AITC_Mags) / length(C3_12_6_AITC_Mags)

  ## 12/13 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1:5,7:10,13,17,19,21:27,29:31,38,39,44,52:55,59:61,63,64,69:80,84:90,92:94,101,102,106), Agonist_Number_From_response_frame = 1, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in A1 from 12/13
A1_12_13_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
A1_12_13_AITC_avg = sum(A1_12_13_AITC_Mags) / length(A1_12_13_AITC_Mags)

  ## 12/13 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,9,10,11,16,19,20,21,23,24,26,27,31:33,35:39,43,44,47,50,52,54:57,59,61,64), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B1 from 12/13
B1_12_13_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
B1_12_13_AITC_avg = sum(B1_12_13_AITC_Mags) / length(B1_12_13_AITC_Mags)

  ## 12/13 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1:8,10:12,14,16,17,21,24,27,29,31:34,37,44:47,49:51,53,54,56,58,62,65,68,70,71,75,76,87,89), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in B3 from 12/13
B3_12_13_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
B3_12_13_AITC_avg = sum(B3_12_13_AITC_Mags) / length(B3_12_13_AITC_Mags)

  ## 12/13 C3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,6,8:11,14,15,17,19,20,28,35,44,46,49,51,56,57,59,67,70,74,76,78,88:92,95,104,105,107,109), Agonist_Number_From_response_frame = 2, Agonist = "AITC", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to AITC in C3 from 12/13
C3_12_13_AITC_Mags = AITC_Peak_Mags_Mut
  ## Find the avg response magnitude
C3_12_13_AITC_avg = sum(C3_12_13_AITC_Mags) / length(C3_12_13_AITC_Mags)

##### WT Cap #####
  ## 10/30 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,5,7,14,17,18,20,24,30,31,34,39,43,47,51,69,70,76,88,96,99,107,110,115,117,121,125,128,136,138,144,152,168,171,179,195,206,228,235,239,254,258,261,264,265,271,273), Agonist_Number_From_response_frame = 3, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B1 from 10/30
B1_10_30_Cap_Mags = as.numeric(c(Cap_Peak_Mags_WT, 0.2476478))
  ## Find the avg response magnitude
B1_10_30_Cap_avg = sum(B1_10_30_Cap_Mags) / length(B1_10_30_Cap_Mags)

  ## 10/30 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,5,6,10:12,15,16,18,24,27,28,30,35,39,48,50:52,54,55,57,63,70:74,79,80,83,84,87:91,93,97:102,106,109,112:115,120:124), Agonist_Number_From_response_frame = 3, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B2 from 10/30
B2_10_30_Cap_Mags = as.numeric(c(Cap_Peak_Mags_WT, 0.4451265, 0.1691926))
  ## Find the avg response magnitude
B2_10_30_Cap_avg = sum(B2_10_30_Cap_Mags) / length(B2_10_30_Cap_Mags)

  ## 10/30 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,2,4,5,7:13,15:21,23:26,29:35,37,40:44,46:49), Agonist_Number_From_response_frame = 4, Agonist = "Capsaicin", Genotype = "WT")
## Store the response net magnitude increases in fluorescence in response to Capsaicin in B3 from 10/30
B3_10_30_Cap_Mags = Cap_Peak_Mags_WT
## Find the avg response magnitude
B3_10_30_Cap_avg = sum(B3_10_30_Cap_Mags) / length(B3_10_30_Cap_Mags)

  ## 10/30 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(6,7,9:13,15,16,18,19,21,23,24,26,28,29,31:36), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B4 from 10/30
B4_10_30_Cap_Mags = Cap_Peak_Mags_WT
  ## Find the avg response magnitude
B4_10_30_Cap_avg = sum(B4_10_30_Cap_Mags) / length(B4_10_30_Cap_Mags)

  ## 11/8 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1:4,6:10,12:14,16,19,20,21,23,24,26,27,29:31,34:43,45,46,48,50,52,53,55,57:59,63:65,67,69:80,82,84), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A1 from 11/8
A1_11_8_Cap_Mags = Cap_Peak_Mags_WT
  ## Find the avg response magnitude
A1_11_8_Cap_avg = sum(A1_11_8_Cap_Mags) / length(A1_11_8_Cap_Mags)

  ## 11/8 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,6:12,14:21,23:29,31,33,35:38,40:42), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B1 from 11/8
B1_11_8_Cap_Mags = Cap_Peak_Mags_WT
  ## Find the avg response magnitude
B1_11_8_Cap_avg = sum(B1_11_8_Cap_Mags) / length(B1_11_8_Cap_Mags)

  ## 11/8 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2:4,6,8:10,13:15,17:21,24,26,28,29), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B2 from 11/8
B2_11_8_Cap_Mags = Cap_Peak_Mags_WT
  ## Find the avg response magnitude
B2_11_8_Cap_avg = sum(B2_11_8_Cap_Mags) / length(B2_11_8_Cap_Mags)

  ## 12/4 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1:5,8,11,12,16:18,20:22,25,28,32:36,38,40:42), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A1 from 12/4
A1_12_4_Cap_Mags = Cap_Peak_Mags_WT
  ## Find the avg response magnitude
A1_12_4_Cap_avg = sum(A1_12_4_Cap_Mags) / length(A1_12_4_Cap_Mags)

  ## 12/4 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3:6,8,10,11,14:16,18,19,21:23,29:31,33,35,37,39,41,44:48,51:53,56,58,60,62,65,68,69,76,79,82,87,89,92,96,97), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A3 from 12/4
A3_12_4_Cap_Mags = Cap_Peak_Mags_WT
  ## Find the avg response magnitude
A3_12_4_Cap_avg = sum(A3_12_4_Cap_Mags) / length(A3_12_4_Cap_Mags)

  ## 12/4 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,4,6,8,9,13:19,22,23,25,27:32,34,36:38,40,42,44), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B3 from 12/4
B3_12_4_Cap_Mags = Cap_Peak_Mags_WT
  ## Find the avg response magnitude
B3_12_4_Cap_avg = sum(B3_12_4_Cap_Mags) / length(B3_12_4_Cap_Mags)

  ## 12/4 B3.2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,4,5,8,9,11:13,18,20,21,24,25,27:31,33,36,37,40:46,48), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B3.2 from 12/4
B3.2_12_4_Cap_Mags = as.numeric(c(Cap_Peak_Mags_WT, 0.2859188 ))
  ## Find the avg response magnitude
B3.2_12_4_Cap_avg = sum(B3.2_12_4_Cap_Mags) / length(B3.2_12_4_Cap_Mags)

  ## 12/11 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2:5,14,16,18:20,22,23,25,34,38,39,44,46,47,49,52,54:57,60,61,65,66,68:70,72,75,77,78,83,84,87,93,95,98,106,108,112,114,121,122,135,144,148,149,153:155,157:159,161,166:169,171,175:178,185,186,195,196,198,209:212,215,216,219,221:226,231,232,235:238,244), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B1 from 12/11
B1_12_11_Cap_Mags = as.numeric(c(Cap_Peak_Mags_WT, 0.1576741, 0.213067))
  ## Find the avg response magnitude
B1_12_11_Cap_avg = sum(B1_12_11_Cap_Mags) / length(B1_12_11_Cap_Mags)

  ## 12/11 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3:6,12,14,16,19,22,23,25:27,34,36:40,44,47:52,56,57,60,62:64,69,73,78:80,93,96,90,91,93:96,98,99,101,105,108,110:112,115:117,119:122,124,128,129,132), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B2 from 12/11
B2_12_11_Cap_Mags = as.numeric(c(Cap_Peak_Mags_WT, 0.1919255))
  ## Find the avg response magnitude
B2_12_11_Cap_avg = sum(B2_12_11_Cap_Mags) / length(B2_12_11_Cap_Mags)

  ## 12/11 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1:5,13,14,17,25,29,31:35,38,39,43:45,48,53:55,57,59,61,64:69), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B4 from 12/11
B4_12_11_Cap_Mags = Cap_Peak_Mags_WT
  ## Find the avg response magnitude
B4_12_11_Cap_avg = sum(B4_12_11_Cap_Mags) / length(B4_12_11_Cap_Mags)

  ## 12/11 C3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,8,10:12,15,18,28:32,36,38,40,43:46,49,50,52:55,57,58,60:63,71,74,76:79,83,87,92,93,100,104:107,109:111,113,115:117,120,121,124,125,131,133:137), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "WT")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B4 from 12/11
C3_12_11_Cap_Mags = as.numeric(c(Cap_Peak_Mags_WT, 0.447322))
  ## Find the avg response magnitude
C3_12_11_Cap_avg = sum(C3_12_11_Cap_Mags) / length(C3_12_11_Cap_Mags)

##### Mut Cap #####
  ## 11/6 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,3,6,7,9,10,15:22,28,31,32,35,36,39,41:46,48,51,52), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A1 from 11_6
A1_11_6_Cap_Mags = Cap_Peak_Mags_Mut
  ## Find the avg response magnitude
A1_11_6_Cap_avg = sum(A1_11_6_Cap_Mags) / length(A1_11_6_Cap_Mags)

  ## 11/6 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1:4,6,8:10,12,13,16,18,19,21,22,23,25:37,41), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A3 from 11_6
A3_11_6_Cap_Mags = Cap_Peak_Mags_Mut
A3_11_6_Cap_avg = sum(A3_11_6_Cap_Mags) / length(A3_11_6_Cap_Mags)

  ## 11/6 B4
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,4,6,8,11,12,14:17,19,23,25:31,33,37:39,42), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A1 from 11_6
B4_11_6_Cap_Mags = as.numeric(c(Cap_Peak_Mags_Mut, 0.2478465 ) )
B4_11_6_Cap_avg = sum(B4_11_6_Cap_Mags) / length(B4_11_6_Cap_Mags)

  ## 11/11 B3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,6,8,12,14,18,19,25,30,32,37,38,44,47,49:51,55,57,59,62,65,68,69,71,73,75,77,78,82,84,87,90,93:98,100,101,104,106,109,111,114,116,117,118,121,127,130,132,137,138,140,142,143,144,148,149,151,154,155,160,165,166,168,173,174,179,190,192,194,197:199,201,203,205,206,208:210,212,213,219,221,226,230,233), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B3 from 11_11
B3_11_11_Cap_Mags = Cap_Peak_Mags_Mut
B3_11_11_Cap_avg = sum(B3_11_11_Cap_Mags) / length(B3_11_11_Cap_Mags)

  ## 11/11 B1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,2,6,13,15:18,22,24,28,29,31:33,35,36,40:43,47,51,54,57:59,62,63,65,69,78,83,85,104,106,109,110,114,115,118,119,121,126:128,130,131,133,136), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B1 from 11_11
B1_11_11_Cap_Mags = Cap_Peak_Mags_Mut
B1_11_11_Cap_avg = sum(B1_11_11_Cap_Mags) / length(B1_11_11_Cap_Mags)

  ## 11/11 A4
  ## Find baselines; adjust the function for this experiment
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(6,16,17,18,24,30,34:36,51,60,61,63:65,69,92,93,96,99,100,103,104,109,110,112,113,115,120:123,130:132,134,137,139,140,144,145,147,149,155,159,160,161,170,173,174,179,182,188,192,195,197,198,200:202,205,206,208,209,213,219,222,230,232:235,237,238,240), Agonist_Number_From_response_frame = 3, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A4 from 11_11
A4_11_11_Cap_Mags = as.numeric(c(Cap_Peak_Mags_Mut, 0.1544591))
A4_11_11_Cap_avg = sum(A4_11_11_Cap_Mags) / length(A4_11_11_Cap_Mags)

  ## 12/6 B1
  ## Find baselines; adjust the function for this experiment
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(4,5,8,12,14,20,22,23,24,30:33,35,38,40,41,44,47,49,51:53,56,59,61,65:67,70,71), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B1 from 12_6
B1_12_6_Cap_Mags = Cap_Peak_Mags_Mut
B1_12_6_Cap_avg = sum(B1_12_6_Cap_Mags) / length(B1_12_6_Cap_Mags)

  ## 12/6 B3
  ## Find baselines; adjust the function for this experiment
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(8,15,18,26,27,30:35), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B3 from 12_6
B3_12_6_Cap_Mags = Cap_Peak_Mags_Mut
B3_12_6_Cap_avg = sum(B3_12_6_Cap_Mags) / length(B3_12_6_Cap_Mags)

  ## 12/13 A1
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(8,9,13,15,17,22,24,27,30,32,35,37:39,41:43,54,57,60,61,66,70,74,77,82,86,96,98,100:105), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A1 from 12_13
A1_12_13_Cap_Mags = Cap_Peak_Mags_Mut
A1_12_13_Cap_avg = sum(A1_12_13_Cap_Mags) / length(A1_12_13_Cap_Mags)

  ## 12/13 A3
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(1,2,6,8:12,14:16,18,19,21,33,35:37,39:43,45,46,49,55,57,62,63,66,67,71,73), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A3 from 12_13
A3_12_13_Cap_Mags = as.numeric(c(Cap_Peak_Mags_Mut, 0.7656966))
A3_12_13_Cap_avg = sum(A3_12_13_Cap_Mags) / length(A3_12_13_Cap_Mags)

  ## 12/13 A4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,4,6,7,8,16,17,19,23,24,26,28:30,34,36,37), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in A4 from 12_13
A4_12_13_Cap_Mags = Cap_Peak_Mags_Mut
A4_12_13_Cap_avg = sum(A4_12_13_Cap_Mags) / length(A4_12_13_Cap_Mags)

  ## 12/13 B2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,4,8,9,14,17,19,21,22,25,26,30:32,36,38:43,45), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B2 from 12_13
B2_12_13_Cap_Mags = Cap_Peak_Mags_Mut
B2_12_13_Cap_avg = sum(B2_12_13_Cap_Mags) / length(B2_12_13_Cap_Mags)

  ## 12/13 B4
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,10,11,12,16:19,21,22,24,25,28:31,33,35,40,42:44,46,50,53,55,58,61:66,68,69,71:74,76,80,82,84:91,93,95,96,99:102,104:106,108:110,112), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B4 from 12_13
B4_12_13_Cap_Mags = as.numeric(c(Cap_Peak_Mags_Mut, 0.7306436))
B4_12_13_Cap_avg = sum(B4_12_13_Cap_Mags) / length(B4_12_13_Cap_Mags)

### Run "Analysis" and "Find_Baselines" functions in C1late script ###
  ## 12/13 C1late
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(3,4,7,12,17,21,22,24:27,29,38,40,41,44:48,52,54,60,62,64:66,68:70,73:76,78,79,81:83), Agonist_Number_From_response_frame = 1, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B4 from 12_13
C1late_12_13_Cap_Mags = Cap_Peak_Mags_Mut
C1late_12_13_Cap_avg = sum(C1late_12_13_Cap_Mags) / length(C1late_12_13_Cap_Mags)

  ## 12/13 C2
  ## Find baselines
Find_Baselines(SmoothFLUO = SmoothFLUO, agonist = agonist, response_frame = response_frame)
  ## Find response magnitudes
Mag_Quant(Responders_List = c(2,6,8,17,20,22,23,25,27,28,30,32,33,39,40,45:47,50,51,53,59,61,63,65,66,70,74,86,88,90:95,97,99,102,107,108,111,113), Agonist_Number_From_response_frame = 2, Agonist = "Capsaicin", Genotype = "Mut")
  ## Store the response net magnitude increases in fluorescence in response to Capsaicin in B4 from 12_13
C2_12_13_Cap_Mags = Cap_Peak_Mags_Mut
C2_12_13_Cap_avg = sum(C2_12_13_Cap_Mags) / length(C2_12_13_Cap_Mags)


sum(c(A1_11_6_Cap_Mags, A3_11_6_Cap_Mags, B4_11_6_Cap_Mags, B3_11_11_Cap_Mags, B1_11_11_Cap_Mags, A4_11_11_Cap_Mags, B1_12_6_Cap_Mags, B3_12_6_Cap_Mags, A1_12_13_Cap_Mags, A3_12_13_Cap_Mags, A4_12_13_Cap_Mags, B4_12_13_Cap_Mags, C1late_12_13_Cap_Mags, C2_12_13_Cap_Mags)) / length(c(A1_11_6_Cap_Mags, A3_11_6_Cap_Mags, B4_11_6_Cap_Mags, B3_11_11_Cap_Mags, B1_11_11_Cap_Mags, A4_11_11_Cap_Mags, B1_12_6_Cap_Mags, B3_12_6_Cap_Mags, A1_12_13_Cap_Mags, A3_12_13_Cap_Mags, A4_12_13_Cap_Mags, B4_12_13_Cap_Mags, C1late_12_13_Cap_Mags, C2_12_13_Cap_Mags))

##### Notes #####
  ## Code to export raw traces: manually change the final subdirectory title in the file destination path when analying a new script/data
    ## Make individual plots as separate png files; these are response profiles for one cell (ROI), each. These give rough ideas of responses by each ROI
#dir.create("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fura2 Images\\2019_10_30\\A4plots")
#for (i in seq_along(FLUO[,1:ncol(FLUO)])) {
#  #tempname = paste("cell_",i,".png", sep="")
#  outpath = file.path("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fura2 Images\\Paper Data Copies\\2019_10_30\\A4plots", paste("cell_",i,".png", sep=""))
#  png(filename = outpath)
#  plot(t, FLUO[,i], xlim = c(t[1],t[length(t)]), ylim = c(0,3))
#  dev.off()
#}