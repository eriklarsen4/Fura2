## 10_30 dissection script



library(tidyverse)
library(ggplot2)
library(purrr)
library(reshape2)
library(dplyr)

##### Data acquisition and prep for perfusion experiment analysis. Developed by Erik Larsen #####

##### Fluorescence Value Acquisition and Prep #####
  ## Import the fluorescence value CSVs.
B4 = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/2019_10_30 B4late.csv")
  ## Remove the extra column that gives the wrong (?) ratio
B4 = B4[,-5]
  ## Replace NAs with the ROI# for each ROI
B4[,1] = rep(1:max(B4$Field., na.rm = T), each = max(B4$Object.))
  ## Reshape the data into something that is usable for plotting for green
Green = reshape(data = B4, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")
  ## Reshape the data into something that is usable for plotting for red
Red = reshape(data = B4, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")
  ## Rename the "Fields" column to "t" (time)
names(Green)[1] = "t"
names(Red)[1] = "t"

  ## Remove any extra frames. Varies by experiment, so check notes
  ## Create a vector that has values in the corresponding range
  ## This will cut extra fluorescence values acquired beyond e.g. 210 frames
trim = seq(1:230)

  ## Apply the trim to the 340, 380 measurements
Green = Green[trim,]
Red = Red[trim,]
  ## Remove the extra column from the dataframes. Depending on run/experiment, hiccups in imaging may need to be cut
Green = Green[-c(1:142),-1]
Red = Red[-c(1:142),-1]

  ## Save the time as a variable to add after computing the green/red ratio
t = ((1:nrow(Green))*3 -3)


##### Background Fluorescence Values Acquisition and Prep #####
  ## Repeat above steps for the background ROIs
B4bkgd = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/2019_10_30 B4bkgd.csv")
B4bkgd = B4bkgd[,-5]
  ## Replace NAs with the ROI# for each ROI
B4bkgd[,1] = rep(1:max(B4bkgd$Field., na.rm = T), each = max(B4bkgd$Object.))
  ## Reshape the data into something that is usable for plotting for green
Greenbkgd = reshape(data = B4bkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")
  ## Reshape the data into something that is usable for plotting for red
Redbkgd = reshape(data = B4bkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")
  ## Rename the "Fields" column to "t" (time)
names(Greenbkgd)[1] = "t"
names(Redbkgd)[1] = "t"

  ## Create a vector that has values from 1 to whatever is supposed to be the end of the experiment (if late in stopping the recording, extra frames will be added).
trim = seq(1:230)

  ## Apply the trim to the 380/340 measurements
Greenbkgd = Greenbkgd[trim,]
Redbkgd = Redbkgd[trim,]
  ## Remove unusable columns from the dataframes; also remove times if experiment contained hiccups or loss of the image field
Greenbkgd = Greenbkgd[-c(1:142),-1]
Redbkgd = Redbkgd[-c(1:142),-1]
  ## Take the mean of the bkgd for each of three ROIs, then use that mean to normalize the fluorescence in each ROI of the fluorescence data frame
AvgGbkgd = rowMeans(Greenbkgd)
AvgRbkgd = rowMeans(Redbkgd)
  ## Subtract the background fluorescence measurements from those of the ROIs.
B4Green = Green - AvgGbkgd
B4Red = Red - AvgRbkgd

  ## Compute the ratio of intracellular Calcium bound to the Fura dye at rest to Calcium bound to Fura dye after agonist induces IC Calcium influx
FLUO = B4Green/B4Red

  ## Set the measured ROI fluorescences to a value of 1 at time 0 for each ROI. This normalizes every frame to the first (1) for all ROI/columns
for (i in 1:ncol(FLUO)){
  FLUO[,i] = (FLUO[,i]) / (FLUO[1,i]) 
}

##### Plotting prep #####

  ## Find the minimum, maximum for plotting scales and removing artefacts
min(FLUO[,])
max(FLUO[,])
lower_lim = min(FLUO[,])
upper_lim = max(FLUO[,])
which(FLUO == min(FLUO), arr.ind = TRUE)
which(FLUO == max(FLUO), arr.ind = TRUE)

# careful with "dev.off()". Avoid when able.
#dev.off()
#plot(t, FLUO[,206], xlim = c(0,331), ylim = c(0,21))

##### Run conditions #####

  ## Create a vector list of agonists in the order they were added in the experiment.
  ## There is a rough, 8s (or ~5 frame) delay for liquid to reach the coverslip from its bottle near the microscope. Once one agonist was turned on, it was the only agonist to which all the cells were exposed.

agonist = c("200uM AITC", "Ringer's", "1uM Cap", "Ringer's", "HiK")
  ## Create a matrix containing the frame number at which the agonist was added and turned off
exposure_frame = cbind(c(1, 15, 35, 55, 75), c(15, 35, 55, 75, 88))
  ## Rename the start and end columns
colnames(exposure_frame) = c("Start", "End")
  ## Rename the rows with the agonists
rownames(exposure_frame) = c(agonist)

  ## Create a second vector to fill with exposure times without Ringer's
response_frame = cbind(c(1,36,76), c(35,75,88))
  ## Rename the start and end columns
colnames(response_frame) = c("Start", "End")
  ## Rename the rows with the agonists
rownames(response_frame) = c(agonist[seq(1,length(agonist),2)])

  ## For an example, find the fluorescence recorded over the given interval of one of the agonists in a given cell
FLUO[exposure_frame[1,1]:exposure_frame[1,2],5]

  ## Create a data frame to house all the agonists of each experiment along with the 0.10 magnitude threshold necessary to classify responders
Agonists_df = data.frame(c(agonist[seq(1,length(agonist),2)]), c(0.10, 0.10, 0.10))
  ## Rename the columns
colnames(Agonists_df) = c("Agonist", "Percent Response Threshold")


##### Smooth the FLUO data (span = 10) #####
  ## Append the time vector to the FLUO data
FLUO = cbind(t, FLUO)
  ## Extract the strings of the columns/ROIs
ROIs = c(colnames(FLUO[c(2:ncol(FLUO))]))
  ## Create a function that will perform local weighted regression (smoothing) on the FLUO data
loess.filter = function (x, span){
  loess(formula = paste(x, "t", sep = "~"), data = FLUO, degree = 1, span = span)$fitted
}
  ## Apply the smoothing to all ROIs in the FLUO data
SmoothFLUO = as.data.frame(lapply(ROIs, loess.filter, span = 10/nrow(FLUO)), col.names = colnames(FLUO[c(2:ncol(FLUO))]))

  ## Append the time vector to the smoothed FLUO data
SmoothFLUO = cbind(t, SmoothFLUO)
##### Find, store the smoothed ROI FLUO slope change (change in smooth 340/380 / change in time) in a matrix #####
Slope_Change_DF = as.data.frame(diff(
  as.matrix(SmoothFLUO[ , c(2:ncol(SmoothFLUO)) ]
  )
)
)
  ## Retain the same column names as from the FLUO and SmoothFLUO data frames
colnames(Slope_Change_DF) = colnames(SmoothFLUO[c(2:ncol(SmoothFLUO))])
  ## Turn off scientific notation
options(scipen = 999)
## Check for structural change in the data
#e.divisive(as.matrix(Slope_Change_DF[,59]), min.size = 10, alpha = 1)

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



  ## Define the function that will compute smoothed responses based on response thresholds used in the Itch field (10% above baseline, also following a 2% increase in slope). Should align with video responses, where ROIs turn green. Count data for the imaging run are exported. Exported raw traces for all ROIs are commented out. Also generates variables used for ggplotting the smoothed data
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
  
  
  ## Which are responders?
  for (i in 1:(nrow(RESPONSE_CALLS)-1)){
    print(which(RESPONSE_CALLS[i,] == "Responder"))
  }
  ## Which are neurons?
  for (i in nrow(RESPONSE_CALLS)){
    print(which(RESPONSE_CALLS[i,] == "Neuron")) 
  }
  ## Show the percent of neurons in the field
  for (i in nrow(RESPONSE_CALLS)){
    print((length(which(RESPONSE_CALLS[i,] == "Neuron"))/length(SmoothFLUO))*100)
  }
 
  RESPONSE_CALLS <<- RESPONSE_CALLS
  Slope_idx_DF <<- Slope_idx_DF
  response_idx <<- response_idx

}


  ## Run the function
Analysis(Agonists_information = Agonists_df, Slope_Change_Threshold = 0.02, Fluorescence_Information = SmoothFLUO)

