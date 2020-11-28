## 12_13 dissection script



library(tidyverse)
library(ggplot2)
library(purrr)
library(reshape2)
library(dplyr)

##### Data acquisition and prep for perfusion experiment analysis. Developed by Erik Larsen #####

##### Fluorescence Value Acquisition and Prep #####
  ## Import the fluorescence value CSVs.
C1late = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/2019_12_13 C1late.csv")


  ## Remove the extra column that gives the wrong (?) ratio
C1late = C1late[,-5]
  ## Replace NAs with the cell number for each cell
C1late[,1] = rep(1:max(C1late$Field., na.rm = T), each = max(C1late$Object.))

  ## Reshape the data into something that is usable for plotting for green
Green = reshape(data = C1late, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")

  ## Reshape the data into something that is usable for plotting for red
Red = reshape(data = C1late, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")

  ## Rename the "Fields" column to "t" (time)
names(Green)[1] = "t"
names(Red)[1] = "t"

  ## Remove any extra frames. Varies by experiment
  ## Create a vector that has values in the corresponding range
  ## This will cut extra fluorescence values acquired beyond e.g. 210 frames (after experiment conclusion)
trim = seq(1:150)

  ## Apply the trim to the 340, 380 measurements
Green = Green[trim,]
Red = Red[trim,]

  ## Remove the extra column from the dataframes. Depending on run/experiment, hiccups in imaging may need to be cut
Green = Green[-c(1:78),-1]
Red = Red[-c(1:78),-1]


  ## Save the time as a variable to append after computing the green/red ratio
t = ((1:nrow(Green))*3 + 234)


##### Background Fluorescence Values Acquisition and Prep #####

  ## Repeat above steps for the background ROIs
C1latebkgd = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/2019_12_13 C1latebkgd.csv")

C1latebkgd = C1latebkgd[,-5]
  ## Replace NAs with the ROI number for each ROI
C1latebkgd[,1] = rep(1:max(C1latebkgd$Field., na.rm = T), each = max(C1latebkgd$Object.))
  ## Reshape the data into something that is usable for plotting for green
Greenbkgd = reshape(data = C1latebkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")
  ## Reshape the data into something that is usable for plotting for red
Redbkgd = reshape(data = C1latebkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")

  ## Rename the "Fields" column to "t" (time)
names(Greenbkgd)[1] = "t"
names(Redbkgd)[1] = "t"

  ## Create a vector that has values from 1 to whatever is supposed to be the end of the experiment (if late in stopping the recording, extra frames will be added).
trim = seq(1:150)
  ## Apply trim to the 380/340 measurements
Greenbkgd = Greenbkgd[trim,]
Redbkgd = Redbkgd[trim,]

  ## Remove unusable columns from the dataframes; also remove times if experiment contained hiccups or loss of the image field
Greenbkgd = Greenbkgd[-c(1:78),-1]
Redbkgd = Redbkgd[-c(1:78),-1]

  ## Take the mean of the bkgd for each background ROI, then use that mean to normalize the fluorescence in each ROI of the fluorescence data frame
AvgGbkgd = rowMeans(Greenbkgd)
AvgRbkgd = rowMeans(Redbkgd)

  ## Subtract the background fluorescence measurements from those of the ROIs.
C1lateGreen = Green - AvgGbkgd
C1lateRed = Red - AvgRbkgd

  ## Compute the ratio of intracellular Calcium bound to the Fura dye at rest to Calcium bound to Fura dye after agonist induces IC Calcium influx
FLUO = C1lateGreen/C1lateRed

  ## Set the measured ROI fluorescences to a value of 1 in the first image for each ROI. This normalizes every frame to the first (1) for all ROI/columns
for (i in 1:ncol(FLUO)){
  FLUO[,i] = (FLUO[,i]) / (FLUO[1,i]) 
}




##### Plotting prep #####


# Find the minimum, maximum for plotting scales and removing artefacts
min(FLUO[,])
max(FLUO[,])
lower_lim = min(FLUO[,])
upper_lim = max(FLUO[,])
which(FLUO == min(FLUO), arr.ind = TRUE)
which(FLUO == max(FLUO), arr.ind = TRUE)

# careful with "dev.off()". Avoid when able. COMMENTED OUT
#dev.off()
#plot(t, FLUO[,206], xlim = c(0,331), ylim = c(0,21))

##### Run conditions #####


  ## Create a vector list of agonists in the order they were added in the experiment.
  ## There is a rough, 15s (or ~5 frame) delay for liquid to reach the coverslip from its bottle near the microscope. Once one agonist was turned on, it was the only agonist to which all the cells were exposed.

agonist = c("1uM Capsaicin", "Ringer's", "HiK")

  ## Create a matrix containing the frame number at which the agonist was added and turned off
exposure_frame = cbind(c(1, 22, 52), c(22, 52, 72))
  ## Rename the start and end columns
colnames(exposure_frame) = c("Start", "End")
  ## Rename the rows with the agonists
rownames(exposure_frame) = c(agonist)

  ## Create a second vector to fill with exposure times without Ringer's physiological solution
response_frame = cbind(c(1, 53), c(52, 72))
  ## Rename the start and end columns
colnames(response_frame) = c("Start", "End")
  ## Rename the rows with the agonists
rownames(response_frame) = c(agonist[seq(1,length(agonist),2)])
  ## For an example, find the fluorescence recorded over the given interval of one of the agonists in a given cell
FLUO[exposure_frame[1,1]:exposure_frame[1,2],5]


  ## Create a data frame to house all the agonists of each experiment along with the 0.10 magnitude threshold necessary to classify responders
Agonists_df = data.frame(c(agonist[seq(1,length(agonist),2)]), c(0.10,0.10))
  ## Rename the columns
colnames(Agonists_df) = c("Agonist", "Percent Response Threshold")



##### Perform Analysis #####
  ## Define the function that will compute smoothed responses based on response thresholds used in the Itch field (10% above baseline, also following a 2% increase in slope). Aligns with video responses, where ROIs turn green. Count data for the imaging run are returned, though verified manually. Exported raw traces for all ROIs are commented out. Also generates variables used for ggplotting the smoothed data
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


  ## Unique to this imaging run
Analysis = function(Agonists_information, Slope_Change_Threshold, Fluorescence_Information){
  
  ##### Create a matrix to house the response calls to all agonists #####
  RESPONSE_CALLS = matrix(nrow = length(agonist[seq(1,length(agonist),2)]), ncol = ncol(SmoothFLUO)-1)
  RESPONSE_CALLS = as.data.frame(RESPONSE_CALLS)
  colnames(RESPONSE_CALLS) = colnames(SmoothFLUO)[2:ncol(SmoothFLUO)]
  rownames(RESPONSE_CALLS) = rownames(response_frame)
  ## Create a matrix that houses where each ROI experiences a ratiometric increase in slope above a threshold
  ## Loop through an populate the slope matrix
  Slope_idx_DF = matrix(nrow = nrow(SmoothFLUO), ncol = (ncol(SmoothFLUO)-1))
  for (i in 1:ncol(Slope_idx_DF)) {
    for (j in 1:nrow(Slope_idx_DF)) {
      Slope_idx_DF[j,i] = paste(which(Slope_Change_DF[j,i] > 0.02), collapse = "")
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
  Response_Stats_df = matrix(nrow = length(agonist[seq(1,length(agonist),2)]), ncol = 5)
  ## Rename columns and rows
  rownames(Response_Stats_df) = c(agonist[seq(1,length(agonist),2)])
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
  
  if ( length(agonist) > 3){
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
  if ( length(agonist) > 3){
    Response_Stats_df[length(agonist[seq(1,length(agonist),2)]),5] = Second_Ag_Max_Avg
  }
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
  if ( length(agonist) > 3 ) {
    Second_Ag_Max_Responses <<- Second_Ag_Max_Responses
    Second_Ag_Max_Avg <<- Second_Ag_Max_Avg
  }
  
  if ( length(agonist) > 6 ) {
    Third_Ag_Max_Responses <<- Third_Ag_Max_Responses
    Third_Ag_Max_Avg <<- Third_Ag_Max_Avg
  }
  HiK_Max_Responses <<- HiK_Max_Responses
  HiK_Max_Avg <<- HiK_Max_Avg
  RESPONSE_CALLS <<- RESPONSE_CALLS
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

  ## Run the function
Analysis(Agonists_information = Agonists_df, Slope_Change_Threshold = 0.02, Fluorescence_Information = SmoothFLUO)
  ## Create the function that returns baseline values that are used for response calls; unique to this experiment
Find_Baselines = function(SmoothFLUO, agonist, response_frame){
  ## Determine baseline ratiometric fluorescences for smoothed fluorescence data frame
  temp = vector()
  Baseline_DF = matrix(nrow = (length(agonist[seq(1,length(agonist),2)]))-1, ncol = ncol(SmoothFLUO)-1)
  Baseline_DF = as.data.frame(Baseline_DF)
  colnames(Baseline_DF) = colnames(SmoothFLUO)[2:ncol(SmoothFLUO)]
  rownames(Baseline_DF) = rownames(response_frame)[1:((length(agonist[seq(1,length(agonist),2)]))-1)]
  
  for (i in 1:ncol(RESPONSE_CALLS)){
    if (RESPONSE_CALLS[1,i] == "Responder"){
      temp = diff(as.numeric(Slope_idx_DF[ , i]))
      Baseline_DF[i] = SmoothFLUO[ which(temp == 0)[which(which(temp == 0 ) %in% response_frame[1,1]:response_frame[1,2]) [1]]-1, i+1]
    }
  }
  Baseline_DF <<- Baseline_DF
}

