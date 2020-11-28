## 12_5 dissection script



library(tidyverse)
library(ggplot2)
library(purrr)
library(reshape2)
library(dplyr)

##### Data acquisition and prep for perfusion experiment analysis. Developed by Erik Larsen #####

##### Fluorescence Value Acquisition and Prep #####
  ## Import the fluorescence value CSVs.
B3 = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/2019_12_6 B3.csv")
  ## Remove the extra column
B3 = B3[,-5]
  ## Replace NAs with the ROI# for each ROI
B3[,1] = rep(1:max(B3$Field., na.rm = T), each = max(B3$Object.))

  ## Reshape the data into something that is usable for plotting for green
Green = reshape(data = B3, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")
  ## Reshape the data into something that is usable for plotting for red
Red = reshape(data = B3, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")
  ## Rename the "Fields" column to "t" (time)
names(Green)[1] = "t"
names(Red)[1] = "t"

  ## Remove any extra frames. Varies by experiment, so check notes
  ## Create a vector that has values in the corresponding range
  ## This will cut extra fluorescence values acquired beyond e.g. 210 frames
trim = seq(1:150)

  ## Apply the trim to the 340, 380 measurements
Green = Green[trim,]
Red = Red[trim,]

  ## Remove the extra column from the dataframes. Depending on run/experiment, hiccups in imaging may need to be cut
Green = Green[,-1]
Red = Red[,-1]

  ## Save the time as a variable to add after computing the green/red ratio
t = ((1:nrow(Green))*3 - 3)


##### Background Fluorescence Values Acquisition and Prep #####

  ## Repeat above steps for the background ROIs
B3bkgd = read.csv("M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/2019_12_6 B3bkgd.csv")
B3bkgd = B3bkgd[,-5]
  ## Replace NAs with the cell number for each cell
B3bkgd[,1] = rep(1:max(B3bkgd$Field., na.rm = T), each = max(B3bkgd$Object.))
  ## Reshape the data into something that is usable for plotting for green
Greenbkgd = reshape(data = B3bkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_RED.1")
  ## Reshape the data into something that is usable for plotting for red
Redbkgd = reshape(data = B3bkgd, idvar = "Field.", timevar = "Object.", direction = "wide", drop = "MEAN_GREEN.1")

  ## Rename the "Fields" column to "t" (time)
names(Greenbkgd)[1] = "t"
names(Redbkgd)[1] = "t"

  ## Create a vector that has values from 1 to whatever is supposed to be the end of the experiment (if late in stopping the recording, extra frames will be added).
trim = seq(1:150)

  ## Apply it to the 380/340 measurements
Greenbkgd = Greenbkgd[trim,]
Redbkgd = Redbkgd[trim,]

  ## Remove unuseable columns from the dataframes; also remove times if experiment contained hiccups or loss of the image field
Greenbkgd = Greenbkgd[,-1]
Redbkgd = Redbkgd[,-1]
  ## Take the mean of the bkgd for each background ROI, then use that mean to normalize the fluorescence in each ROI of the fluorescence data frame
AvgGbkgd = rowMeans(Greenbkgd)
AvgRbkgd = rowMeans(Redbkgd)

  ## Subtract the background fluorescence measurements from those of the ROIs.
B3Green = Green - AvgGbkgd
B3Red = Red - AvgRbkgd

  ## Compute the ratio of intracellular Calcium bound to the Fura dye at rest to Calcium bound to Fura dye after agonist induces IC Calcium influx
FLUO = B3Green/B3Red

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

# careful with "dev.off()". Avoid when able. COMMENTED OUT
#dev.off()
#plot(t, FLUO[,206], xlim = c(0,331), ylim = c(0,21))

##### Run conditions #####

  ## Create a vector list of agonists in the order they were added in the experiment.
  ## There is a rough, 15s (or ~5 frame) delay for liquid to reach the coverslip from its bottle near the microscope. Once one agonist was turned on, it was the only agonist to which all the cells were exposed.

agonist = c("Ringer's", "5mM Beta-Alanine", "Ringer's", "1uM Capsaicin", "Ringer's", "HiK")

  ## Create a matrix containing the frame number at which the agonist was added and turned off
exposure_frame = cbind(c(1, 15, 45, 75, 105, 135), c(15, 45, 75, 105, 135, 150))
  ## Rename the start and end columns
colnames(exposure_frame) = c("Start", "End")
  ## Rename the rows with the agonists
rownames(exposure_frame) = c(agonist)


  ## Create a second vector to fill with exposure times without Ringer's
response_frame = cbind(c(16,76,136), c(75,135,150))
  ## Rename the start and end columns
colnames(response_frame) = c("Start", "End")
  ## Rename the rows with the agonists
rownames(response_frame) = c(agonist[seq(2,length(agonist),2)])

  ## For an example, find the fluorescence recorded over the given interval of one of the agonists in a given cell
FLUO[exposure_frame[1,1]:exposure_frame[1,2],5]

  ## Create a data frame to house all the agonists of each experiment along with the 0.10 magnitude threshold necessary to classify responders
Agonists_df = data.frame(c(agonist[seq(2,length(agonist),2)]), c(0.10, 0.10, 0.10))
  ## Rename the columns
colnames(Agonists_df) = c("Agonist", "Percent Response Threshold")

