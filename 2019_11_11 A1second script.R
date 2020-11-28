## 11_10 dissection script

##### DON'T USE THIS SCRIPT #####

library(tidyverse)
library(ggplot2)
library(purrr)
library(reshape2)
library(dplyr)

##### Data acquisition and prep for perfusion experiment analysis. Developed by Erik Larsen #####

##### Fluorescence Value Acquisition and Prep #####
  ## Import the fluorescence value CSVs.
A1second = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/2019_11_11 A1second.csv")
  ## Remove the extra column
A1second = A1second[,-5]
  ## Replace NAs with the ROI# for each ROI
A1second[,1] = rep(1:max(A1second$Field., na.rm = T), each = max(A1second$Object.))
  ## Reshape the data into something that is usable for plotting for green
Green = reshape(data = A1second, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")
  ## Reshape the data into something that is usable for plotting for red
Red = reshape(data = A1second, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")
  ## Rename the "Fields" column to "t" (time)
names(Green)[1] = "t"
names(Red)[1] = "t"
  ## Remove any extra frames. Varies by experiment, so check notes
  ## Create a vector that has values in the corresponding range
  ## This will cut extra fluorescence values acquired beyond e.g. 210 frames
trim = seq(1:125)
  ## Apply the trim to the 340, 380 measurements
Green = Green[trim,]
Red = Red[trim,]
  ## Remove the extra column from the dataframes. Depending on run/experiment, hiccups in imaging may need to be cut
Green = Green[-c(1:76),-1]
Red = Red[-c(1:76),-1]
  ## Save the time as a variable to add after computing the green/red ratio since computing math cannot be done in specified columns of dataframes in R
t = ((1:nrow(Green))*3 - 3)


##### Background Fluorescence Values Acquisition and Prep #####
  ## Repeat above steps for the background ROIs
A1secondbkgd = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/2019_11_11 A1secondbkgd.csv")
A1secondbkgd = A1secondbkgd[,-5]
  ## Replace NAs with the ROI# for each ROI
A1secondbkgd[,1] = rep(1:max(A1secondbkgd$Field., na.rm = T), each = max(A1secondbkgd$Object.))
  ## Reshape the data into something that is usable for plotting for green
Greenbkgd = reshape(data = A1secondbkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")
  ## Reshape the data into something that is usable for plotting for red
Redbkgd = reshape(data = A1secondbkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")
  ## Rename the "Fields" column to "t" (time)
names(Greenbkgd)[1] = "t"
names(Redbkgd)[1] = "t"
  ## Create a vector that has values from 1 to whatever is supposed to be the end of the experiment (if late in stopping the recording, extra frames will be added).
trim = seq(1:125)
  ## Apply the trim to the 380/340 measurements
Greenbkgd = Greenbkgd[trim,]
Redbkgd = Redbkgd[trim,]
  ## Remove unusable columns from the dataframes; also remove times if experiment contained hiccups or loss of the image field
Greenbkgd = Greenbkgd[-c(1:76),-1]
Redbkgd = Redbkgd[-c(1:76),-1]
  ## Take the mean of the bkgd for each background ROI, then use that mean to normalize the fluorescence in each ROI of the fluorescence data frame
AvgGbkgd = rowMeans(Greenbkgd)
AvgRbkgd = rowMeans(Redbkgd)
  ## Subtract the background fluorescence measurements from those of the ROIs.
A1secondGreen = Green - AvgGbkgd
A1secondRed = Red - AvgRbkgd
  ## Compute the ratio of intracellular Calcium bound to the Fura dye at rest to Calcium bound to Fura dye after agonist induces IC Calcium influx
FLUO = A1secondGreen/A1secondRed
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
  ## There is a rough, 15s (or ~5 frame) delay for liquid to reach the coverslip from its bottle near the microscope. Once one agonist was turned on, it was the only agonist to which all the cells were exposed.

agonist = c("Ringer's", "25uM CYM5442")
  ## Create a matrix containing the frame number at which the agonist was added and turned off
exposure_frame = cbind(c(1, 19), c(19,49))
  ## Rename the start and end columns
colnames(exposure_frame) = c("Start", "End")
  ## Rename the rows with the agonists
rownames(exposure_frame) = c(agonist)
  ## Create a second vector to fill with exposure times without Ringer's
response_frame = cbind(c(14), c(49))
  ## Rename the start and end columns
colnames(response_frame) = c("Start", "End")
  ## Rename the rows with the agonists
rownames(response_frame) = c(agonist[seq(2,length(agonist),2)])
  ## For an example, find the fluorescence recorded over the given interval of one of the agonists in a given cell
FLUO[exposure_frame[1,1]:exposure_frame[1,2],5]

  ## Create a data frame to house all the agonists of each experiment along with the 0.10 magnitude threshold necessary to classify responders
Agonists_df = data.frame(c(agonist[seq(2,length(agonist),2)]), c(0.10))
  ## Rename the columns
colnames(Agonists_df) = c("Agonist", "Percent Response Threshold")


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
    for (j in 1:(nrow(response_frame)) ){
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
  for (i in 1:(nrow(RESPONSE_CALLS))){
    if (length(which(RESPONSE_CALLS[i,] == "Responder")) == 0){
      Response_Stats_df[i,1] = 0
    } else {
      Response_Stats_df[i,1] = length(which(RESPONSE_CALLS[i,] == "Responder")) 
    }
  }
  ## Which are responders?
  for (i in 1:(nrow(RESPONSE_CALLS))){
    print(which(RESPONSE_CALLS[i,] == "Responder"))
  }
  ## Of all responders, what were the maximum responses to each respective agonist?
  for (j in 1:(nrow(RESPONSE_CALLS))) {
    Response_Stats_df[j,4] = max(SmoothFLUO[ response_frame[j,1]:response_frame[j,2] , c(2:ncol(SmoothFLUO))])
  }
  
  ## Create the variables that will store peak agonist responses and average peak agonist responses
  First_Ag_Max_Responses = vector()
  First_Ag_Max_Avg = vector()
  
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
  
  
  ## Put them in the Response Stats DF
  Response_Stats_df[1,5] = First_Ag_Max_Avg
  Response_Stats_df[nrow(RESPONSE_CALLS),5] = HiK_Max_Avg
  
  FLUO <<- FLUO
  SmoothFLUO <<- SmoothFLUO
  RESPONSE_CALLS <<- RESPONSE_CALLS
  Slope_Change_DF <<- Slope_Change_DF
  Slope_idx_DF <<- Slope_idx_DF
  response_idx <<- response_idx
  Response_Stats_df <<- Response_Stats_df
  
  
  return(Response_Stats_df)

}
Analysis(Slope_Change_Threshold = 0.02)

##### GGplotting #####
Find_Peak_Times_And_Mags = function (Cell_Number_Plus_One){
  ## Concatenate the max smoothed responses and their respective frame numbers (when they occur) for non-HiK agonists
  ## Initialize variables
  Peak_Mags = as.numeric(vector())
  Peak_Frames = as.numeric(vector())
  #HiK_Peak_Mag = as.numeric(vector())
  #HiK_Peak_Frame = as.numeric(vector())
  for (j in 1:(nrow(response_frame)) ){
    ## Concatenate the max smoothed responses to non-HiK agonists into a variable
    Peak_Mags[j] = cbind(max(
      SmoothFLUO[ response_frame[j,1]:response_frame[j,2], Cell_Number_Plus_One ] ))
    ## Concatenate the respective frame numbers when these responses occur
    Peak_Frames[j] = cbind(which(SmoothFLUO[, Cell_Number_Plus_One] == max( SmoothFLUO[ response_frame[j,1]:response_frame[j,2], Cell_Number_Plus_One ] ) ))
  }
  Peak_Mags = Peak_Mags[!is.na(Peak_Mags) == TRUE]
  Peak_Frames = Peak_Frames[!is.na(Peak_Frames) == TRUE]
  
  ## Create a variable for the max smoothed response to HiK
  #HiK_Peak_Mag = max(
  #  SmoothFLUO[ response_frame[nrow(response_frame),1]:response_frame[nrow(response_frame),2], Cell_Number_Plus_One ]  )
  ## Create a variable for when the max HiK response occurs
  #HiK_Peak_Frame = which(SmoothFLUO[, Cell_Number_Plus_One] == max( SmoothFLUO[ response_frame[nrow(response_frame),1]:response_frame[nrow(response_frame),2], Cell_Number_Plus_One ] ) )
  
  ## Export all these values for a given ROI into the global environment for plotting
  Peak_Mags <<- Peak_Mags
  Peak_Frames <<- Peak_Frames
  #HiK_Peak_Mag <<- HiK_Peak_Mag
  #HiK_Peak_Frame <<- HiK_Peak_Frame
  return(print(c("Peak Agonist Frames", Peak_Frames))) #"Peak HiK Frame", #HiK_Peak_Frame)))
}



Find_Peak_Times_And_Mags(Cell_Number_Plus_One = 62)

## Make plots of perfusion data
CaImaging = ggplot(data = FLUO, aes(x = FLUO$t, y = FLUO$MEAN_GREEN.1.61)) +
  geom_rect(data = FLUO, xmin = t[exposure_frame[2,1]], xmax = t[exposure_frame[2,2]], ymin = -Inf, ymax = Inf, fill = "rosybrown1", alpha = 0.01) +
  #geom_rect(data = FLUO, xmin = t[exposure_frame[4,1]], xmax = t[exposure_frame[4,2]], ymin = -Inf, ymax = Inf, fill = "peachpuff2", alpha = 0.01) +
  #geom_rect(data = FLUO, xmin = t[exposure_frame[6,1]], xmax = t[exposure_frame[6,2]], ymin = -Inf, ymax = Inf, fill = "darkseagreen", alpha = 0.01) +
  geom_text(x = t[exposure_frame[2,1]]+(t[exposure_frame[2,2]]-t[exposure_frame[2,1]])/2, y  = Response_Stats_df[nrow(RESPONSE_CALLS),4]-.05, label = "CYM5442\napplication", color = "rosybrown1", size = 5, fontface = "bold") +
  #geom_text(x = t[exposure_frame[4,1]]+(t[exposure_frame[4,2]]-t[exposure_frame[4,1]])/2, y  = Response_Stats_df[nrow(RESPONSE_CALLS),4]-.05, label = "CYM5442\napplication", color = "peachpuff2", size = 5, fontface = "bold") +
  #geom_text(x = t[exposure_frame[6,1]]+(t[exposure_frame[6,2]]-t[exposure_frame[6,1]])/2, y  = Response_Stats_df[nrow(RESPONSE_CALLS),4]-.05, label = "140mM K\napplication", color = "darkseagreen", size = 3.75, fontface = "bold") +
  geom_vline(xintercept = t[which(Slope_Change_DF[,61] > 0.02)], color = "blue") +
  geom_vline(xintercept = t[Peak_Frames], color = "red", size = 0.9) +
  #geom_vline(xintercept = t[DIV$estimates][-c(1,150)]], color = "brown") +
  geom_point(color = "black") +
  geom_line(color = "black", size = 1.1) +
  geom_smooth(method = "loess", span = 10/nrow(FLUO), method.args = list(degree = 1), aes(x = FLUO$t, y = FLUO$MEAN_GREEN.1.61), color = "navy") +
  coord_cartesian(xlim = c(0,t[length(t)]), ylim = c(0.5,4)) +
  labs(title = "Mut DRG Response to\n25\u03BCM CYM5442", x = "Time (s)", y = "340:380\nNormalized Fluorescence (a.u.)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank()) +
  geom_hline(yintercept = 1.0, color = "black", linetype = "solid")

CaImaging