## 11_10 dissection script



library(tidyverse)
library(ggplot2)
library(purrr)
library(reshape2)
library(dplyr)

##### Data acquisition and prep for perfusion experiment analysis. Developed by Erik Larsen #####

##### Fluorescence Value Acquisition and Prep #####
# Import the fluorescence value CSVs.
# Make sure the csv's are trimmed appropriately. Get rid of summary rows and columns (the first set of values), and delete extra rows so that the first row containing text is the first row of the file.
A3late = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/2019_11_11/A3late.csv")


# Remove the extra column that gives the wrong (?) ratio
A3late = A3late[,-5]
# Replace NAs with the cell number for each cell
A3late[,1] = rep(1:max(A3late$Field., na.rm = T), each = max(A3late$Object.))

# Reshape the data into something that is usable for plotting for green
Green = reshape(data = A3late, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")

# Reshape the data into something that is usable for plotting for red
Red = reshape(data = A3late, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")

# Rename the "Fields" column to "t" (time)
names(Green)[1] = "t"
names(Red)[1] = "t"

# Remove any extra frames. Varies by experiment, so check notes
# Create a vector that has values in the corresponding range
# This will cut extra fluorescence values acquired beyond e.g. 210 frames
trim = seq(1:160)

# Apply it to the 340, 380 measurements
Green = Green[trim,]
Red = Red[trim,]

# Remove the extra column from the dataframes, along with the first two frames (rows 1 and 2) because the microscope/filter wheel shifts the entire image, depending on exposure settings, and finally settles on a baseline image in the third frame
Green = Green[-c(1:82),-1]
Red = Red[-c(1:82),-1]


# Save the time as a variable to add after computing the green/red ratio since computing math cannot be done in specified columns of dataframes in R
t = ((1:nrow(Green))*3)


##### Background Fluorescence Values Acquisition and Prep #####

# Repeat above steps for the background ROIs
A3latebkgd = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/2019_11_11/A3earlybkgd.csv")

A3latebkgd = A3latebkgd[,-5]

# Replace NAs with the cell number for each cell
A3latebkgd[,1] = rep(1:max(A3latebkgd$Field., na.rm = T), each = max(A3latebkgd$Object.))

# Reshape the data into something that is usable for plotting for green
Greenbkgd = reshape(data = A3latebkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")

# Reshape the data into something that is usable for plotting for red
Redbkgd = reshape(data = A3latebkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")

# Rename the "Fields" column to "t" (time)
names(Greenbkgd)[1] = "t"
names(Redbkgd)[1] = "t"

# Remove any extra frames. May vary by experiment.
# Create a vector that has values from 1 to whatever is supposed to be the end of the experiment (If late in stopping the recording, extra frames will be added).
# This will cut extra fluorescence values acquired.
trim = seq(1:160)

# Apply it to the 380/340 measurements
Greenbkgd = Greenbkgd[trim,]
Redbkgd = Redbkgd[trim,]

# Create the first column as a vector to add to the data frame.
# (Finding the ratio of green/red cannot be done across the data frame with the time column still in the data frame from the original data. 

# Remove the column from the dataframes
Greenbkgd = Greenbkgd[-c(1:82),-1]
Redbkgd = Redbkgd[-c(1:82),-1]

# Take the mean of the bkgd for each background ROI, then use that mean to normalize the fluorescence in each ROI of the fluorescence data frame
AvgGbkgd = rowMeans(Greenbkgd)
AvgRbkgd = rowMeans(Redbkgd)

# Subtract the background fluorescence measurements from those of the ROIs.
A3lateGreen = Green - AvgGbkgd
A3lateRed = Red - AvgRbkgd

# Compute the ratio of intracellular Calcium bound to the Fura dye at rest to Calcium bound to Fura dye after agonist induces IC Calcium influx
FLUO = A3lateGreen/A3lateRed

# Set the measured ROI fluorescences to a value of 1 at time 0 for each ROI. This normalizes every frame to the first (1) for all ROI/columns
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

##### Responses Output #####

# For each experiment, these variables must be adjusted to perform analysis:

# Which agonists were used/the order in which they were used
# The time when each agonist was turned on and off (exposure)
# The threshold above normalized baseline which will determine a valid response to agonist exposure


# Create a vector list of agonists in the order they were added in the experiment.

# There is a rough, 15s (or ~5 frame) delay for liquid to reach the coverslip from its bottle near the microscope. Once one agonist was turned on, it was the only agonist to which all the cells were exposed.

agonist = c("Ringer's", "2.5uM LY344862", "Ringer's", "HiK")

# Create a matrix containing the frame number at which the agonist was added and turned off
exposure_frame = cbind(c(1, 13, 38, 58), c(13, 38, 58, 78))

# Rename the start and end columns
colnames(exposure_frame) = c("Start", "End")
# Rename the rows with the agonists
rownames(exposure_frame) = c(agonist)


# Create a second vector to fill with exposure times without Ringer's
response_frame = cbind(c(13,58), c(58,78))
# Rename the start and end columns
colnames(response_frame) = c("Start", "End")
# Rename the rows with the agonists
rownames(response_frame) = c(agonist[seq(2,length(agonist),2)])

# For an example, find the fluorescence recorded over the given interval of one of the agonists in a given cell
FLUO[exposure_frame[1,1]:exposure_frame[1,2],5]


# Create a data frame to house all the agonists of each experiment along with the percent threshold we deem is worthy of a salient response
Agonists_df = data.frame(c(agonist[seq(2,length(agonist),2)]), c(0.50, 1.5))
# Rename the columns
colnames(Agonists_df) = c("Agonist", "Percent Response Threshold")



# Define the function that will plot responses for each individual ROI and compute responses based on arbitrary computational thresholds that may align with video responses, where ROIs clearly turn green
Analysis = function(Agonists_information, Exposure_Information, Fluorescence_Information){
  
  ##### Establish Baselines, Thresholds, and Determine Responses for Classification #####
  
  # These values depend on the onset of agonist in each experiment and subsequent responses
  # Set an arbitrary baseline average FL (pre-agonist) for all cells prior to addition of the first agonist in frame indicated by notes
  # This has been revised to include "wash out", where Ringer's solution is added between agonists, and is therefore included in agonist information
  # Ringer's solution is therefore the baseline for agonist addition immediately following.
  
  # Create the baseline fluorescence df
  baseline_df = matrix(nrow = length(agonist)/2, ncol = length(FLUO))
  # Fill with the first agonist
  baseline_df[1,] = colSums(FLUO[exposure_frame[1,1]:exposure_frame[1,2], 1:length(FLUO)], dims = 1)/(exposure_frame[1,2]-exposure_frame[1,1]+1)
  # Fill with the next agonist
  baseline_df[nrow(baseline_df)-1,] = colSums(FLUO[exposure_frame[length(agonist)-3,1]:exposure_frame[length(agonist)-3,2], 1:length(FLUO)], dims = 1)/(exposure_frame[length(agonist)-3,2]-exposure_frame[length(agonist)-3,1])
  # Fill with the Hi K
  baseline_df[nrow(baseline_df),] = colSums(FLUO[exposure_frame[length(agonist)-1,1]:exposure_frame[length(agonist)-1,2], 1:length(FLUO)], dims = 1)/(exposure_frame[length(agonist)-1,2]-exposure_frame[length(agonist)-1,1])
  
  
  # Calculate the response threshold for the first agonist
  # Include the thresholds for each agonist in the Agonists_df
  Agonist_threshold_df = matrix(nrow = length(agonist)/2, ncol = length(FLUO))
  
  # Keep the same cell #
  rownames(Agonist_threshold_df) = c(agonist[seq(2,length(agonist),2)])
  
  # Append the threshold values of each cell
  for (i in 1:nrow(baseline_df))
    Agonist_threshold_df[i,] = baseline_df[i,]*Agonists_df[i,2] + baseline_df[i,]
  
  
  # Determine responders to the first agonist
  # Create an empty vector to fill with whether each cell (ROI) responded to an agonist
  Agonist_responses = data.frame()
  
  # Fill the vector with ROIs that responded to the agonist or not, based on the pre-determined, arbitrary threshold of what constitutes a "response"
  for (i in 1:ncol(Agonist_threshold_df)){
    for (j in 1:nrow(Agonist_threshold_df)){
      Agonist_responses[j,i] = if_else(any(FLUO[response_frame[j,1]:response_frame[j,2],i] > Agonist_threshold_df[1,i]), "Responder", "Non-responder")
    }
  }
  
  # Make response calls to High Potassium addition; replace "responders" as "Neurons", and call non-responders to High Potassium "Non-neurons"
  for (i in 1:ncol(Agonist_threshold_df)){
    Agonist_responses[nrow(Agonists_df),i] = if_else(Agonist_responses[nrow(Agonists_df),i] == "Responder", "Neuron", 
                                                     if_else(FLUO[71,i] - baseline_df[nrow(Agonists_df),i] > Agonist_threshold_df[nrow(Agonists_df),i], "Neuron", "Non-neuron"))
  }
  
  # Replace false-positive responses to an agonist; be careful.. some ROIs may already have had their gradients altered by agonist addition, so they may not respond to High Potassium after depleting their internal calcium stores..
  for (i in 1:nrow(Agonist_threshold_df)) {
    Agonist_responses[i,which(Agonist_responses[i,] == "Responder" & Agonist_responses[length(agonist),] == "Non-neuron")] = "Non-responder"
  }
  
  
  
  # Create the response dataframe to store response stats to each agonist
  Response_Stats_df = matrix(nrow = length(agonist)/2, ncol = 4)
  # Rename columns and rows
  rownames(Response_Stats_df) = c(agonist[seq(2,length(agonist),2)])
  colnames(Response_Stats_df) = c("# of Responders", "% of Neurons responding", "% Neurons of ROIs", "Maximum Response Magnitude")
  # Remove NAs
  Response_Stats_df[is.na(Response_Stats_df)] = 0
  
  # Fill it with the number of responders, neurons, percent of neurons responding to an agonist, the max response of each agonist, and the percent of neurons out of all the ROIs. Also printed below
  
  # See how many ROIs are classified as "responders", which ones they are, and the percentages of responses
  for (i in 1:(nrow(Agonist_threshold_df)-1)){
    # How many responders?
    Response_Stats_df[i,1] = length(which(Agonist_responses[i,] == "Responder"))
  }
  
  for (i in nrow(Agonist_threshold_df)){
    # How many neurons?
    Response_Stats_df[nrow(Response_Stats_df),1] = length(which(Agonist_responses[i,] == "Neuron"))
  }
  
  for (i in 1:(nrow(Agonist_threshold_df)-1)){
    # Which are responders?
    print(which(Agonist_responses[i,] == "Responder"))
  }
  
  for (i in nrow(Agonists_df)){
    # Which are neurons?
    print(which(Agonist_responses[i,] == "Neuron")) 
  }
  
  for (i in 1:(nrow(Agonist_threshold_df))){
    # What were the maximum responses to each agonist?
    Response_Stats_df[i,4] = max(FLUO[response_frame[i,1]:response_frame[i,2],])
  }
  
  
  for (i in 1:(nrow(Agonist_threshold_df)-1)){
    # Determine percentages
    # Of neurons responding to each agonist
    Response_Stats_df[i,2] = (length(which(Agonist_responses[i,] == "Responder"))/length(which(Agonist_responses[nrow(Agonists_df),] == "Neuron")))*100
  }
  
  
  for (i in 1:ncol(Response_Stats_df)){
    for (j in 1:nrow(Response_Stats_df)){
      # Find the percent of neurons in the field
      Response_Stats_df[nrow(Response_Stats_df),3] =  (Response_Stats_df[nrow(Response_Stats_df),1]/ncol(FLUO))*100
    }
  }
  
  
  for (i in nrow(Agonists_df)){
    # Show the percent of neurons in the field
    print((length(which(Agonist_responses[i,] == "Neuron"))/length(FLUO))*100)
  }
  
  
  # Make individual plots as separate png files; these are response profiles for one cell (ROI), each. These give rough ideas of responses by each ROI
  dir.create("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fura2 Images\\2019_11_11\\A3lateplots")
  for (i in seq_along(FLUO[,1:ncol(FLUO)])) {
    #tempname = paste("cell_",i,".png", sep="")
    outpath = file.path("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fura2 Images\\2019_11_11\\A3lateplots", paste("cell_",i,".png", sep=""))
    png(filename = outpath)
    plot(t, FLUO[,i], xlim = c(t[1],t[length(t)]), ylim = c(0,3))
    dev.off()
  }
  
  return(Response_Stats_df)
  
  #for (i in 1:nrow(Agonists_df)){
  # Determine percentages
  # Of responders to each agonist; % responders of all neurons, and percent of neurons in the field
  #print((length(which(Agonist_responses[i,] == "Responder"))/length(which(Agonist_responses[nrow(Agonists_df),] == "Neuron")))*100)
  #print((length(which(Agonist_responses[i,] == "Neuron"))/length(FLUO))*100)
  #}
}

# Run the function
Analysis(Agonists_information = Agonists_df, Exposure_Information = exposure_frame, Fluorescence_Information = FLUO)

##### GGPlots #####


# Add the time vector again
FLUO = cbind(t, FLUO)

# Make plots of perfusion data
CaImaging = ggplot(data = FLUO, aes(x = FLUO$t, y = FLUO$MEAN_GREEN.1.102)) + geom_vline(xintercept = 30, color = "forest green") + geom_vline(xintercept = 120, color = "forest green") + geom_vline(xintercept = 210, color = "red") + geom_vline(xintercept = 300, color = "red") + geom_vline(xintercept = 390, color = "purple") + geom_line(color = "blue") + coord_cartesian(xlim = c(0,450), ylim = c(0,3)) + labs(x = "Time (s)", y = "340:380 Normalized Fluorescence (a.u.)", title = "In vitro WT DRG Responses to 2.5uM LY344862, 25uM CYM544, 140mM K") + theme(plot.title = element_text(hjust = 0.5)) + theme_light() + geom_hline(yintercept = 1.0, color = "black", linetype = "solid")

CaImaging