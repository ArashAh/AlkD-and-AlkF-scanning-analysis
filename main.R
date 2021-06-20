
##### Packages #####

library(magrittr)
library(tidyverse)
library(viridis)
library(broom)
library(plotly)
library(e1071)
library(randomForest)
library(h2o)
library(party)
library(data.table)
library(microbenchmark)
library(dplyr)
library(DT)
library(shiny)
library(stringr)
library(scales)
library(Rcpp)
library(RcppRoll)
library(zoo)
library(XML)
library(ggplot2)
library(jsonlite)
library(caret)
library(R.matlab)
source("functions.R")



############### 1. Trajectory detection and noise exclusion #################

# This step is performed separately for each protein

# Direct to the spatially filtered data:

filtered.data.path <- "processed/1_AlkD_filtered_data.rds"

collect.filtered.data <- readRDS(filtered.data.path)

collect.noise.excluded <- NULL
collect.trajectory.address <- NULL

for(i in collect.filtered.data %>%
    filter(spatial.filter == "Yes") %>%
    "$"(data.set.name) %>% 
    unique()){
  
 #  Separate the data sets into raw and filtered
  
  raw.data.set <- collect.filtered.data %>% 
    filter(data.set.name == i)
  
  filtered.data.set <- raw.data.set %>% 
    filter(spatial.filter == "Yes")

filtered.matix <- MakeItMatrix(RawData = raw.data.set, 
             InputData = filtered.data.set)  
  
detected.trajectories <- FindTrajectory(FilteredDataSet = filtered.matix,
                                               dxmax = 600, dxmin = 0.1, 
                                               dymax = 200) %>% 
  MakeLongForm() %>% 
  CorrectBlinking() %>% 
  MakeZeroNA()


# Noise exlusion

noise.excluded <- detected.trajectories %>% 
  ExcludeNoise()

collect.noise.excluded <- bind_rows(collect.noise.excluded, noise.excluded)

}

#  Direct to save the noise excluded data: 
noise.excluded.path <- "processed/2_AlkD_detected_noise_excluded.rds"

saveRDS(collect.noise.excluded, noise.excluded.path)


# Number of data sets with existing data, frames and observations after displacement filtering. 


post.displ.filt.summary <- collect.noise.excluded %>%
  filter(displacement.filter == "Yes", duration.filter == "Yes", !is.na(data.set.name)) %>% 
  group_by(data.set.name) %>%
  summarise(frames = length(unique(frame.number)), observations = n())

post.displ.filt.summary %>% 
  summarise(data.sets = n(), frames = sum(frames), observations = sum(observations)) 


post.larg.displ.filt.summary <- collect.noise.excluded %>%
  filter(large.displacement.filter == "Yes", duration.filter == "Yes",!is.na(data.set.name)) %>% 
  group_by(data.set.name) %>%
  summarise(frames = length(unique(frame.number)), observations = n())

post.larg.displ.filt.summary %>%
  summarise(data.sets = n(), frames = sum(frames), observations = sum(observations)) 


#############################2. Data transformation ###############################

# Direct to noise excluded data: 
noise.excluded.path <- "processed/2_AlkD_detected_noise_excluded.rds"

collect.noise.excluded <- readRDS(noise.excluded.path)


# Based on the file names with which data sets are saved TransformData() function might differ 
# There are two versions of this function since we had two different style of data set names  

transformed.data <- collect.noise.excluded %>% 
  TransformData() %>%
  mutate(enzyme = "AlkD") %>% 
  arrange(trajectory.unique.id, frame.number) 

# A summary of the transformed data

transformed.data %>%
  filter(duration.filter == "Yes", !is.na(data.set.name)) %>%
  group_by(displacement.filter, large.displacement.filter) %>%
  summarise(length(unique(trajectory.unique.id)))

# Exporting the address of trajectories for cutting the frames in FIJI

trajectory.address <- transformed.data %>% 
  filter(enzyme == "AlkD") %>% 
MakeTrajectoryAddress()


# Path of file names for saving the trajectory address files 

trajectory.address.path1 <- ""

write.table(trajectory.address,
            trajectory.address.path1,
            sep = "\t", col.names = TRUE)

to.cut.data.sets <- trajectory.address$file.name


# Path of file names for saving the trajectory address files number 2  

trajectory.address.path2 <- ""

write.table(to.cut.data.sets,trajectory.address.path2, 
            quote = FALSE, row.names = FALSE)



# Visual inspection of the detected trajectories: 
# Inspect the trajectories in the output video files of the FIJI script 
# Import the name of trajectories that are not visually approved 

remove.trajectories <-c()
q <- "Yes"
while(q=="Yes")
{
  remove.trajectories <- c(remove.trajectories, 
                         as.integer(readline("insert the trajectory id: ")))
}


remove.trajectories <- remove.trajectories[!is.na(remove.trajectories)]

a <- transformed.data %>%
  filter(duration.filter == "Yes", 
         large.displacement.filter == "Yes", 
         !is.na(data.set.name)) %>% 
  .$trajectory.unique.id %>% unique()

remove.trajectories1 <- setdiff(a, remove.trajectories)

  transformed.data %<>% 
  group_by(trajectory.unique.id) %>% 
  mutate(duration.filter = 
           ifelse(trajectory.unique.id %in% remove.trajectories1, "No", duration.filter))

# Removing the noise from the data
  
transformed.data %<>%  
  filter(duration.filter == "Yes", !is.na(data.set.name), displacement.filter == "Yes")



# Fitting the DNA regression line and setting on.DNA filter

fit.on.trajectory <- transformed.data %>% 
  filter(large.displacement.filter == "Yes", 
         !is.na(data.set.name), 
         duration.filter == "Yes") %>% 
  group_by(data.set.name) %>% 
do(tidy(lm(Y ~ X, data = .))) %>% 
  mutate(fit.slope = estimate[term == "X"], 
         fit.intercept = estimate[term == "(Intercept)"]) %>%
  filter(term == "X") %>%  
  select(data.set.name, fit.slope, fit.intercept)


transformed.data %<>% 
  filter(duration.filter == "Yes", !is.na(data.set.name)) %>% 
  left_join(.,fit.on.trajectory) %>% 
  rowwise %>% 
  mutate(corrected_X = RotatePoints(X, Y, -fit.slope)[1], 
         corrected_Y = RotatePoints(X, Y, -fit.slope)[2]) %>% 
  ungroup() 


on.dna.filter.set <- transformed.data %>% 
  filter(large.displacement.filter == "Yes") %>% 
  group_by(data.set.name) %>% 
  summarise(mid.y = mean(corrected_Y)) 

 
transformed.data <-  left_join(transformed.data, on.dna.filter.set, by = "data.set.name")


transformed.data %<>%  mutate(on.dna.filter = ifelse(corrected_Y > mid.y-200 & 
                                corrected_Y < mid.y+200, "Yes", "No")) %>% 
  filter(!is.na(fit.slope))

transformed.data %>% 
  filter(large.displacement.filter == "Yes") %>% 
  group_by(data.set.name) %>% 
  summarise(length(unique(trajectory.unique.id))) %>% view()

# Corrrecting for datasets with low number of large.displacement filter trajectories


correct.fit <- transformed.data %>% 
  filter(large.displacement.filter == "Yes") %>% 
  group_by(data.set.name) %>% 
  summarise(large.dis.in.traj = length(unique(trajectory.unique.id))) %>%
  mutate(correction.filter = ifelse(large.dis.in.traj <10 , "No" , "Yes")) %>% 
  select(-large.dis.in.traj)

transformed.data <- left_join(transformed.data, correct.fit)

transformed.data %<>%
  mutate(displacement.filter = ifelse(correction.filter == "No" & 
                                        large.displacement.filter == "No", "No", displacement.filter))

# Overview of number of trajectories 

transformed.data %>%
  filter(duration.filter == "Yes", !is.na(data.set.name)) %>%
  group_by(displacement.filter, large.displacement.filter) %>%
  summarise(length(unique(trajectory.unique.id)))


# This part of the analysis is performed separately for each proteins 

# The data for wild-type AlkD and AlkF are saved together 
# and the data for two rounds of AlkF mutant are saved separately 

# Direct to transformed data path and files names 

transformed.data.path <- "processed/3_AlkD_AlkF_transfromed_data.rds"

saveRDS(transformed.data, transformed.data.path)
transformed.data <- readRDS(transformed.data.path)



###################3.Data collecting #############################################

# Direct to path and file names of the different proteins that are saved from previous step 

mut.Alkf1.path <- "processed/3_mut-AlkF_transfromed_data.rds"
mut.Alkf2.path <- "processed/3_mut-AlkF2_transfromed_data.rds"
Alks.path <- "processed/3_AlkD_AlkF_transfromed_data.rds"
Endogg.path <- "processed/3_AlignedDataBlinkingCorrected.2017-04-07.rds"

# Import all data from previous and current experiments 

mut.Alkf1.data <- readRDS(mut.Alkf1.path)
mut.Alkf2.data <- readRDS(mut.Alkf2.path)
Alks.data <- readRDS(Alks.path)
Endogg.data <- readRDS(Endogg.path)


# See the naming of the columns from different data sets 

names(Alks.data)
names(Endogg.data)
names(mut.Alkf1.data)
names(mut.Alkf2.data)

# Make unique trajectory ID and arrange the columns for a row bind 


mut.Alkf1.data %<>% 
  mutate(trajectory.unique.id = 
           paste(enzyme, "1", trajectory.unique.id, sep = "_")) %>% 
  select(-duration.filter, -trajectory.id, -fit.slope, -fit.intercept, 
         -mid.y) %>% 
  mutate(displacement.filter = ifelse(displacement.filter == "Yes", 
                                      TRUE, FALSE),
         large.displacement.filter = ifelse(large.displacement.filter == "Yes", 
                                            TRUE, FALSE),
         on.dna.filter = ifelse(on.dna.filter == "Yes", 
                                TRUE, FALSE)) 

## fix the frame interval for AlkF data where exposure time was written instead of frame interval 

mut.Alkf2.data$frame.interval[mut.Alkf2.data$frame.interval == 10] = 13.5


mut.Alkf2.data %<>% 
  mutate(trajectory.unique.id = 
           paste(enzyme, "2", trajectory.unique.id, sep = "_")) %>% 
  select(-duration.filter, -trajectory.id, -fit.slope, -fit.intercept, 
         -mid.y) %>% 
  mutate(displacement.filter = ifelse(displacement.filter == "Yes", 
                                      TRUE, FALSE),
         large.displacement.filter = ifelse(large.displacement.filter == "Yes", 
                                            TRUE, FALSE),
         on.dna.filter = ifelse(on.dna.filter == "Yes", 
                                TRUE, FALSE))%>% 
  filter(large.displacement.filter) %>% 
  mutate(corrected_X = X, corrected_Y = Y)


Alks.data %<>% 
  mutate(trajectory.unique.id = 
                        paste(enzyme, trajectory.unique.id, sep = "_")) %>% 
  select(-duration.filter, -trajectory.id, -fit.slope, -fit.intercept, 
         -mid.y) %>% 
  mutate(displacement.filter = ifelse(displacement.filter == "Yes", 
                                      TRUE, FALSE),
         large.displacement.filter = ifelse(large.displacement.filter == "Yes", 
                                            TRUE, FALSE),
         on.dna.filter = ifelse(on.dna.filter == "Yes", 
                                            TRUE, FALSE))

 Endogg.data %<>% 
   ungroup() %>% 
   mutate(Unique_trajectory_ID =
            paste(Enzyme, Unique_trajectory_ID, sep = "_")) %>% 
   select(Unique_trajectory_ID, File_name, Frame_number, X, Y, 
          `Delta_X > 300`, Visual_confirmation, Intensity, Enzyme, NaCl, 
          Frame_interval, corrected_X, corrected_Y, On_DNA) %>% 
   set_colnames(c("trajectory.unique.id", "data.set.name", "frame.number",
                  "X", "Y", "displacement.filter", "large.displacement.filter",
                  "intensity", "enzyme", "salt.concentration", "frame.interval",
                  "corrected_X", "corrected_Y", "on.dna.filter")) %>% 
   mutate(displacement.filter = ifelse(displacement.filter == "Yes", 
                                       TRUE, FALSE),
          large.displacement.filter = ifelse(large.displacement.filter == "Yes", 
                                       TRUE, FALSE))

# Combining all data sets 
 
tidy.data <- bind_rows(mut.Alkf1.data, mut.Alkf2.data, Alks.data, Endogg.data) %>% 
  arrange(trajectory.unique.id, frame.number)

# Renaming the proteins 

 tidy.data$enzyme <-  factor(tidy.data$enzyme, 
                             levels = c("hOgg1", "EndoV", "mEndoV", "AlkD", "AlkF", "mut-AlkF"))
 levels(tidy.data$enzyme) <-c("hOGG1","EndoV", "mut-EndoV", "AlkD", "AlkF", "mut-AlkF")
 

 # Smoothening the effect of the on DNA filter to avoid chopping trajectories 
 
 tidy.data %<>%  
   group_by(trajectory.unique.id) %>% 
   mutate(correc.on.dna.filter = sum(on.dna.filter)/length(trajectory.unique.id))

 tidy.data %<>% 
   mutate(on.dna.filter = ifelse(correc.on.dna.filter > 0.5 | 
                                   large.displacement.filter, TRUE, FALSE )) %>% 
   select(-correction.filter, -correc.on.dna.filter)
 
 
 # Summary: 
 
 tidy.data %>% 
   filter(large.displacement.filter, on.dna.filter) %>%  
   group_by(enzyme) %>% 
   summarise(num.larg.dis.traj = length(unique(trajectory.unique.id)),
             num.dis.large.frames = n()) %>% 
   left_join(., tidy.data %>% 
               filter(displacement.filter, on.dna.filter) %>%  
               group_by(enzyme) %>% 
               summarise(num.dis.traj = length(unique(trajectory.unique.id)) ,
                         num.dis.frames = n())) %>% 
   view()
 
 
 # Direct to path and file name to save the collected data that contains all proteins 
 collected.data.path <- "processed/4_collected_data_2020_03_04.rds"
 
 saveRDS(tidy.data, collected.data.path)
 
 

############## 4.Instantaneous diffusion rate #####################################

tidy.data <- readRDS(collected.data.path)


sourceCpp('MSDsComplete.cpp')

# localMSDs are given in nm2 we need D= MSD/2t in terms of  micrometer2/s
# therefore we divide the values of MSD by (2000 * frame.interval (in ms)) 

local.MSD <- tidy.data %>% 
  group_by(trajectory.unique.id) %>% 
  mutate(time=frame.number-min(frame.number)) %>% 
  arrange(time) %>% 
  mutate(local.msd.05 = localMSDcomplete(corrected_X, 5)/ (2000*frame.interval), 
         local.msd.07 = localMSDcomplete(corrected_X, 7)/ (2000*frame.interval), 
         local.msd.10 = localMSDcomplete(corrected_X, 10)/ (2000*frame.interval), 
         local.msd.15 = localMSDcomplete(corrected_X, 15)/ (2000*frame.interval)) 

# Plot instantaneous diffusion rate distribution for different widths of the moving window

local.msd.path <- "processed/6_local_MSD_2020_03_04.rds"
saveRDS(local.MSD, local.msd.path)
local.MSD <- readRDS(local.msd.path)

# width of moving window = 5

local.MSD%>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.05!=0) %>% 
  group_by(enzyme) %>%
  ggplot(aes(x=local.msd.05)) +
  facet_wrap(~enzyme, ncol = 3) +
  geom_histogram(bins=65, position = "stack",
                 aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                 fill = "black") + 
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Fractions")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.025,0.05,0.075))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm")) +
    ggtitle("Moving window d = 5")


# width of moving window = 5

local.MSD%>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.07!=0) %>% 
  group_by(enzyme) %>%
  ggplot(aes(x=local.msd.07)) +
  facet_wrap(~enzyme, ncol = 3) +
  geom_histogram(bins=65, position = "stack",
                 aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                 fill = "black") + 
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Fractions")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.025,0.05,0.075))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm")) +
  ggtitle("Moving window d = 7")


# width of moving window = 10 

local.MSD%>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.10!=0) %>% 
  group_by(enzyme) %>%
  ggplot(aes(x=local.msd.10)) +
  facet_wrap(~enzyme, ncol = 3) +
  geom_histogram(bins=65, position = "stack",
                 aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                 fill = "black") + 
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Fractions")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.025,0.05,0.075))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm")) +
  ggtitle("Moving window d = 10")

# width of moving window = 15

local.MSD%>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.15!=0) %>% 
  group_by(enzyme) %>%
  ggplot(aes(x=local.msd.15)) +
  facet_wrap(~enzyme, ncol = 3) +
  geom_histogram(bins=65, position = "stack",
                 aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                 fill = "black") + 
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Fractions")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.025,0.05,0.075))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm")) +
  ggtitle("Moving window d = 15")


# Overlay of density plots for all proteins 

local.MSD %>% 
  filter(local.msd.05!=0, displacement.filter, on.dna.filter) %>% 
  ggplot(aes(x=local.msd.05, color = enzyme)) + 
  geom_density( size = 3) + 
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Density")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.5, 1))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm"))+ 
  scale_fill_discrete(name = "Energy barrier")


## Calculation of spatial precision 

spatial.precision <- readRDS("processed/spatial.precision.data.rds")

spatial.precision %>% filter(data == "spatial.precision") %>% 
  ungroup() %>% summarise(length(unique(trajectory.unique.id))) 

spatial.precision %>% filter(data == "spatial.precision") %>% 
  group_by(trajectory.unique.id) %>% summarise(sdx= sd(X), sdy= sd(Y)) %>% ungroup() %>% 
  summarise(meansdx= mean(sdx), sdsdx = sd(sdx))

spatial.precision %>% filter(local.msd.05!=0, data == "spatial.precision" )%>% 
  ungroup() %>%  
  summarise(ave.diffusion.rate = mean(local.msd.05, na.rm = T), 
            sd.diffusion.rate = sd(local.msd.05, na.rm = T)) 



################# 5.Simulation ###################################################

# calculation of average diffusion rate of Each proteins

local.MSD %>% filter(local.msd.15!=0, displacement.filter, on.dna.filter)%>% 
  group_by(enzyme) %>% 
  summarise(ave.diffusion.rate = mean(local.msd.15, na.rm = T), 
            sd.diffusion.rate = sd(local.msd.15, na.rm = T)) 




## simulation 

##hOgg1

sim.hOgg1 <- NULL
for (i in 1:4000) {
  D = 0.13
  t = 0.0235
  SD <- sqrt(2 * D * t)
  x0 = 0
  n = sample(x = 5:100, size = 1)
  Step <- c(x0, rnorm(n = (n - 1), mean = 0, sd = SD))
  corrected_X <- cumsum(Step)
  Traj <- data.frame(frame.number = seq(1, n, 1), corrected_X, 
                     Step, 
                     trajectory.unique.id = paste("sim-hOgg1", i, sep = "-"))
  sim.hOgg1 <- rbind(sim.hOgg1, Traj)
  print(i)
}

sim.hOgg1 %<>% mutate(enzyme = "hOgg1", frame.interval = 23.5, 
                      corrected_X = corrected_X*1000)

## Endov 
sim.EndoV <- NULL
for (i in 1:4000) {
  D = 0.65
  t = 0.0075
  SD <- sqrt(2 * D * t)
  x0 = 0
  n = sample(x = 5:100, size = 1)
  Step <- c(x0, rnorm(n = (n - 1), mean = 0, sd = SD))
  corrected_X <- cumsum(Step)
  Traj <- data.frame(frame.number = seq(1, n, 1), corrected_X, 
                     Step, 
                     trajectory.unique.id = paste("sim-EndoV", i, sep = "-"))
  sim.EndoV <- rbind(sim.EndoV, Traj)
  print(i)
}

sim.EndoV %<>% mutate(enzyme = "EndoV", frame.interval = 7.5, 
                      corrected_X = corrected_X*1000)


## mEndoV

sim.mEndoV <- NULL
for (i in 1:4000) {
  D = 0.69
  t = 0.0075
  SD <- sqrt(2 * D * t)
  x0 = 0
  n = sample(x = 5:100, size = 1)
  Step <- c(x0, rnorm(n = (n - 1), mean = 0, sd = SD))
  corrected_X <- cumsum(Step)
  Traj <- data.frame(frame.number = seq(1, n, 1), corrected_X, 
                     Step, 
                     trajectory.unique.id = paste("sim-mEndoV", i, sep = "-"))
  sim.mEndoV <- rbind(sim.mEndoV, Traj)
  print(i)
}

sim.mEndoV %<>% mutate(enzyme = "mEndoV", frame.interval = 7.5, 
                       corrected_X = corrected_X*1000)

## AlkD

sim.AlkD <- NULL
for (i in 1:4000) {
  D = 0.08
  t = 0.0335
  SD <- sqrt(2 * D * t)
  x0 = 0
  n = sample(x = 5:100, size = 1)
  Step <- c(x0, rnorm(n = (n - 1), mean = 0, sd = SD))
  corrected_X <- cumsum(Step)
  Traj <- data.frame(frame.number = seq(1, n, 1), corrected_X, 
                     Step, 
                     trajectory.unique.id = paste("sim-AlkD", i, sep = "-"))
  sim.AlkD <- rbind(sim.AlkD, Traj)
  print(i)
}

sim.AlkD %<>% mutate(enzyme = "AlkD", frame.interval = 33.5, 
                       corrected_X = corrected_X*1000)


## AlkF

sim.AlkF <- NULL
for (i in 1:4000) {
  D = 0.13
  t = 0.0335
  SD <- sqrt(2 * D * t)
  x0 = 0
  n = sample(x = 5:100, size = 1)
  Step <- c(x0, rnorm(n = (n - 1), mean = 0, sd = SD))
  corrected_X <- cumsum(Step)
  Traj <- data.frame(frame.number = seq(1, n, 1), corrected_X, 
                     Step, 
                     trajectory.unique.id = paste("sim-AlkF", i, sep = "-"))
  sim.AlkF <- rbind(sim.AlkF, Traj)
  print(i)
}

sim.AlkF %<>% mutate(enzyme = "AlkF", frame.interval = 33.5, 
                     corrected_X = corrected_X*1000)


## AlkF mnutant 

sim.mut.AlkF <- NULL

for (i in 1:4000) {
  D = 0.55
  t = 0.0335
  SD <- sqrt(2 * D * t)
  x0 = 0
  n = sample(x = 5:100, size = 1)
  Step <- c(x0, rnorm(n = (n - 1), mean = 0, sd = SD))
  corrected_X <- cumsum(Step)
  Traj <- data.frame(frame.number = seq(1, n, 1), corrected_X, 
                     Step, 
                     trajectory.unique.id = paste("sim-mut-AlkF", i, sep = "-"))
  sim.mut.AlkF <- rbind(sim.mut.AlkF, Traj)
  print(i)
}

sim.mut.AlkF %<>% mutate(enzyme = "mut-AlkF", frame.interval = 33.5, 
                     corrected_X = corrected_X*1000)



sim.data <- bind_rows( sim.hOgg1, sim.EndoV, sim.mEndoV ,sim.AlkD, sim.AlkF, sim.mut.AlkF)

sim.data$enzyme <-  factor(sim.data$enzyme, levels = c("hOgg1", "EndoV", 
                                                         "mEndoV", "AlkD", "AlkF", "mut-AlkF"))
levels(sim.data$enzyme) <-c("hOGG1","EndoV", "mut-EndoV", "AlkD","AlkF",  "mut-AlkF")


## localMSDs of simulations 

local.MSD.sim <- sim.data %>% 
  group_by(trajectory.unique.id) %>% 
  mutate(time=frame.number-min(frame.number)) %>% 
  arrange(time) %>% 
  mutate(local.msd.05 = localMSDcomplete(corrected_X, 5)/ (2000*frame.interval), 
         local.msd.07 = localMSDcomplete(corrected_X, 7)/ (2000*frame.interval), 
         local.msd.10 = localMSDcomplete(corrected_X, 10)/ (2000*frame.interval), 
         local.msd.15 = localMSDcomplete(corrected_X, 15)/ (2000*frame.interval)) 


# Control by recalculating the average diffusion rate of the simulated data  

local.MSD.sim %>% filter(local.msd.15!=0)%>% 
  group_by(enzyme) %>% 
  summarise(ave.diffusion.rate = mean(local.msd.15, na.rm = T), 
            sd.diffusion.rate = sd(local.msd.15, na.rm = T)) 

# Overlay of the simulation and the real data 

local.MSD %>%
  filter(local.msd.05!=0, displacement.filter, on.dna.filter) %>% 
  ggplot(aes(x=local.msd.05)) + 
  facet_wrap(~enzyme, ncol = 3) + 
  geom_histogram(bins=100, position = "stack",
                 aes(y=20*(..count..)/
                       tapply(..count..,..PANEL..,sum)[..PANEL..])) +
  geom_density(data =local.MSD.sim %>% filter(local.msd.05!=0) ,
               aes(x=local.msd.05), 
               color = "black")+
  geom_density(color = "red")+
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Density")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.5, 1))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm"))+ 
  scale_fill_discrete(name = "Diffusion mode")

simulated.data.path <- "Processed/7_local_MSD_sim_2020_03_04.rds"

saveRDS(local.MSD.sim, simulated.data.path)
local.MSD.sim <-readRDS(simulated.data.path)

# Extract the upper limit of diffusion for simulated random walk corresponding to each real protein
# The upper limit is defined as max of diffusion rate where density distribution of instantaneous 
# diffusion rate drops to 0.01


# AlkD

to.plot.info <- local.MSD.sim %>%
  filter(local.msd.05!=0, enzyme == "AlkD") %>% ungroup() %>% 
  select(local.msd.05) 

p.info <- ggplot(to.plot.info) + 
  geom_density(aes(x=local.msd.05))+ 
  scale_x_log10()



p.data <- ggplot_build(p.info)

density.data <- p.data[[1]][1] %>% as.data.frame()
density.data %<>% mutate(limit = ifelse(scaled < 0.01, "No", "Yes")) %>% 
  mutate(limit2 = ifelse(limit==lead(limit, n=1), 0,10^x))

density.data %>% ggplot()+ geom_point(aes(x =x, y =scaled, color = limit))


# upper limit of AlkD 

max(density.data$limit2, na.rm = T)

# AlkF

to.plot.info <- local.MSD.sim %>%
  filter(local.msd.05!=0, enzyme == "AlkF") %>% ungroup() %>% 
  select(local.msd.05) 

p.info <- ggplot(to.plot.info) + 
  geom_density(aes(x=local.msd.05))+ 
  scale_x_log10()



p.data <- ggplot_build(p.info)

density.data <- p.data[[1]][1] %>% as.data.frame()
density.data %<>% mutate(limit = ifelse(scaled < 0.01, "No", "Yes")) %>% 
  mutate(limit2 = ifelse(limit==lead(limit, n=1), 0,10^x))

density.data %>% ggplot()+ geom_point(aes(x =x, y =scaled, color = limit))

max(density.data$limit2, na.rm = T)


# Diffusion plot with upper limit in place 

local.MSD %>% filter(enzyme %in% c("AlkD", "AlkF", "mut-AlkF")) %>% 
  filter(local.msd.05!=0, displacement.filter, on.dna.filter) %>% 
  ggplot(aes(x=local.msd.05)) + 
  facet_wrap(~enzyme, ncol = 3) + 
  geom_histogram(bins=100, position = "stack",
                 aes(y=20*(..count..)/
                       tapply(..count..,..PANEL..,sum)[..PANEL..])) +
  geom_density(data =local.MSD.sim %>% filter(local.msd.05!=0,
                                              enzyme %in% c("AlkD", "AlkF", "mut-AlkF")) ,
               aes(x=local.msd.05), 
               color = "black")+
  geom_density(color = "red")+
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Density")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.5, 1))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm"))+ 
  scale_fill_discrete(name = "Diffusion mode")+ 
  geom_segment(data=data.frame(x=c(0.49, 0.8, 0.8), 
                               y= c(0.7, 0.7, 0.7), 
                               enzyme=c("AlkD", "AlkF", "mut-AlkF")), 
               aes(x= x, y= 0, xend= x ,yend=y), inherit.aes=FALSE, size= 0.6)



#################### 6. Markov classification ####################################################

local.MSD <- readRDS(local.msd.path)
local.MSD.sim <- readRDS(simulated.data.path)

# Producing text files of trajectories to be imported to MATLAB for HMM analysis

# AlkD-23.5

traj.id <- local.MSD %>%
  filter(enzyme=="AlkD", frame.interval== 23.5, displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkD", frame.interval== 23.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X, corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/AlkD/23.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}


# AlkD-33.5

traj.id <- local.MSD %>%
  filter(enzyme=="AlkD", frame.interval== 33.5, displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkD", frame.interval== 33.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X, corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/AlkD/33.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}

# AlkF-13.5

traj.id <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 13.5,displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 13.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X,corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/AlkF/13.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}


# AlkF-23.5

traj.id <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 23.5,displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 23.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X,corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/AlkF/23.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}


# AlkF-33.5

traj.id <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 33.5,displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 33.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X,corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/AlkF/33.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}


# mut-AlkF-13.5

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 13.5,displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 13.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X,corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/mutAlkF/13.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}




# mut-AlkF-15.5

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 15.5,displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 15.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X,corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/mutAlkF/15.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}



# mut-AlkF-18.5

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 18.5,displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 18.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X,corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/mutAlkF/18.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}


# mut-AlkF-23.5

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 23.5,displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 23.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X,corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/mutAlkF/23.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}


# mut-AlkF-33.5

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 33.5,displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 33.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    select(corrected_X,corrected_Y) 
  
  write.table(to.save,
              paste("processed/matlab2/to_markov/mutAlkF/33.5/" ,j,".txt", 
                    sep = ""), sep="\t",col.names =F,row.names = F ,quote= F)
  print(j) 
  
}

#

# Importing the output of HMM analysis from MATLAB

# AlkD-23.5

data <- readMat('Processed/matlab2/from_markov/AlkD_23.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]


traj.id <- local.MSD %>%
  filter(enzyme=="AlkD", frame.interval== 23.5, displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkD", frame.interval== 23.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.AlkD.23.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.AlkD.23.5 <- bind_rows(to.out.AlkD.23.5, to.save)
  print(j) 
}



# AlkD-33.5

data <- readMat('Processed/matlab2/from_markov/AlkD_33.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]


traj.id <- local.MSD %>%
  filter(enzyme=="AlkD", frame.interval== 33.5, displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkD", frame.interval== 33.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.AlkD.33.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.AlkD.33.5 <- bind_rows(to.out.AlkD.33.5, to.save)
  print(j) 
}

# AlkF-13.5

data <- readMat('Processed/matlab2/from_markov/AlkF_13.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]

traj.id <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 13.5, displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 13.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.AlkF.13.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.AlkF.13.5 <- bind_rows(to.out.AlkF.13.5, to.save)
  print(j) 
}



# AlkF-23.5

data <- readMat('Processed/matlab2/from_markov/AlkF_23.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]

traj.id <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 23.5, displacement.filter, on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 23.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.AlkF.23.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.AlkF.23.5 <- bind_rows(to.out.AlkF.23.5, to.save)
  print(j) 
}



# AlkF-33.5

data <- readMat('Processed/matlab2/from_markov/AlkF_33.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]

traj.id <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 33.5, displacement.filter,on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="AlkF", frame.interval== 33.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.AlkF.33.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.AlkF.33.5 <- bind_rows(to.out.AlkF.33.5, to.save)
  print(j) 
}


# mut-AlkF-13.5

data <- readMat('Processed/matlab2/from_markov/mut_AlkF_13.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 13.5, displacement.filter,on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 13.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.mut.AlkF.13.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.mut.AlkF.13.5 <- bind_rows(to.out.mut.AlkF.13.5, to.save)
  print(j) 
}


# mut-AlkF-15.5

data <- readMat('Processed/matlab2/from_markov/mut_AlkF_15.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 15.5, displacement.filter,on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 15.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.mut.AlkF.15.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.mut.AlkF.15.5 <- bind_rows(to.out.mut.AlkF.15.5, to.save)
  print(j) 
}



# mut-AlkF-18.5

data <- readMat('Processed/matlab2/from_markov/mut_AlkF_18.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 18.5, displacement.filter,on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 18.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.mut.AlkF.18.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.mut.AlkF.18.5 <- bind_rows(to.out.mut.AlkF.18.5, to.save)
  print(j) 
}




# mut-AlkF-23.5

data <- readMat('Processed/matlab2/from_markov/mut_AlkF_23.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 23.5, displacement.filter,on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 23.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.mut.AlkF.23.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.mut.AlkF.23.5 <- bind_rows(to.out.mut.AlkF.23.5, to.save)
  print(j) 
}



# mut-AlkF-33.5

data <- readMat('Processed/matlab2/from_markov/mut_AlkF_33.5_modes_2_06_Mar_2020_HMMresult.mat')
a1 <- data[[2]][[10]][[5]]

traj.id <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 33.5, displacement.filter,on.dna.filter) %>%
  extract2("trajectory.unique.id") %>% 
  unique()

to.save2 <- local.MSD %>%
  filter(enzyme=="mut-AlkF", frame.interval== 33.5, displacement.filter, on.dna.filter) %>%
  ungroup()

j=0
to.out.mut.AlkF.33.5 <- NULL

for(i in traj.id){
  j=j+1
  to.save <- to.save2 %>% 
    filter(trajectory.unique.id == i) %>% 
    add.col(a1[[j]] %>% unlist())
  
  to.out.mut.AlkF.33.5 <- bind_rows(to.out.mut.AlkF.33.5, to.save)
  print(j) 
}





to.out <- bind_rows(to.out.AlkD.23.5, to.out.AlkD.33.5,
                    to.out.AlkF.13.5,to.out.AlkF.23.5, to.out.AlkF.33.5, 
                    to.out.mut.AlkF.13.5, to.out.mut.AlkF.15.5, to.out.mut.AlkF.18.5, 
                    to.out.mut.AlkF.23.5, to.out.mut.AlkF.33.5)

segmented.local.MSD2 <- left_join(local.MSD %>%
                                   filter(displacement.filter,on.dna.filter,
                                          enzyme %in% c("AlkD", "AlkF", "mut-AlkF"))%>%
                                   mutate(upper.limit = ifelse(enzyme == "AlkD",  0.49, 0.8),
                                          energy.barrier =
                                            ifelse(log(upper.limit/(local.msd.05)) <= 2 ,
                                                   "Ea < 2KbT", "Ea > 2KbT")),
                                 to.out) %>%
  arrange(trajectory.unique.id, frame.number) %>%
  ungroup()

segmented.data.path <- "processed/9_segmented_2_states_local_MSD_2020_04_08.rds"

saveRDS(segmented.local.MSD2, segmented.data.path)


segmented.local.MSD2<- readRDS(segmented.data.path)



# Ploting 

# Density 

segmented.local.MSD2 %>% 
  filter(local.msd.05!=0, displacement.filter, on.dna.filter, 
         !is.na(Markov.state)) %>% 
  ggplot() + 
  facet_wrap(~enzyme, ncol = 1) + 
  geom_density(aes(x=local.msd.05, color= as.factor(Markov.state)), size =1.5)+
  geom_density(data =local.MSD.sim %>% filter(enzyme %in% c("AlkD", "AlkF", "mut-AlkF"),
               local.msd.05!=0),
               aes(x=local.msd.05), 
               color = "black")+
    scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Density")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.5, 1))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm"))+ 
  scale_fill_discrete(name = "Markov state")



# histograms stack

segmented.local.MSD2%>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.05!=0) %>% 
  group_by(enzyme) %>%
  ggplot(aes(x=local.msd.05)) +
  facet_wrap(~enzyme, ncol = 1) +
  geom_histogram(bins=100, position = "stack",
                 aes(y=21*(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..], 
                     fill= as.factor(Markov.state))) + 
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Fractions")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.5,1))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm")) + 
  ggtitle("Moving window d = 5") + geom_density()+
  geom_density(data =local.MSD.sim %>% filter(enzyme %in% c("AlkD", "AlkF", "mut-AlkF"),
                                              local.msd.05!=0),
               aes(x=local.msd.05), 
               color = "red", size = 1.4)


# histograms dodge 

segmented.local.MSD2%>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.05!=0) %>% 
  group_by(enzyme) %>%
  ggplot(aes(x=local.msd.05)) +
  facet_wrap(~enzyme, ncol = 1) +
  geom_histogram(bins=100, position = "dodge",
                 aes(y=21*(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..], 
                     fill= as.factor(Markov.state))) + 
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Fractions")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.5,1))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm")) + 
  ggtitle("Moving window d = 5") + geom_density()+
  geom_density(data =local.MSD.sim %>% filter(enzyme %in% c("AlkD", "AlkF", "mut-AlkF"),
                                              local.msd.05!=0),
               aes(x=local.msd.05), 
               color = "red", size = 1.4)



# Energy barrier plot

segmented.local.MSD2%>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.05!=0) %>% 
  group_by(enzyme) %>%
  ggplot(aes(x=local.msd.05)) + 
  facet_wrap(~enzyme, ncol = 1) +
  geom_histogram(bins=100, position = "dodge",
                 aes(y=21*(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..], 
                     fill= as.factor(energy.barrier))) + 
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  ylab("Fractions")+ 
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.5,1))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm")) + 
  ggtitle("Moving window d = 5") + geom_density()+
  geom_density(data =local.MSD.sim %>% filter(enzyme %in% c("AlkD", "AlkF", "mut-AlkF"),
                                              local.msd.05!=0),
               aes(x=local.msd.05), 
               color = "red", size = 1.4)



# Markov-energy comparison 

markov.result <- segmented.local.MSD2 %>% 
  filter(local.msd.05!=0, !is.na(Markov.state)) %>% 
  group_by(enzyme, Markov.state) %>% 
  summarise(diffusion.rate = mean(local.msd.05), 
            diffusion.rate.error = sd(local.msd.05)/sqrt(n()),
            occupancy = n()) %>% 
  group_by(enzyme) %>% 
  mutate(l= sum(occupancy)) %>% 
  mutate(occupancy = occupancy/l) %>% 
  select(-l) %>%
  mutate(Segmentation = "Markov states") %>% 
 mutate(upper.limit =  1) %>% 
   mutate(energy.barrier =
            ifelse(log(upper.limit/(diffusion.rate)) <= 2 , 
                   "Ea < 2KbT","Ea > 2KbT")) %>% 
   select(-upper.limit)


energy.result <- segmented.local.MSD2 %>% 
  filter(local.msd.05!=0, !is.na(Markov.state)) %>% 
  group_by(enzyme, energy.barrier) %>% 
  summarise(diffusion.rate = mean(local.msd.05), 
            diffusion.rate.error = sd(local.msd.05)/sqrt(n()), 
            occupancy = n()) %>% 
  group_by(enzyme) %>% 
  mutate(l= sum(occupancy)) %>% 
  mutate(occupancy = occupancy/l) %>% select(-l) %>%
  mutate(Segmentation = "Energy barrier")

diffusion.occupancy <- bind_rows( energy.result, markov.result)

diffusion.occupancy$energy.barrier <- 
  factor(diffusion.occupancy$energy.barrier,
         levels = c('Ea > 2KbT', 'Ea < 2KbT'))

diffusion.occupancy %>% ggplot() + 
  geom_point(aes(x = diffusion.rate,
                 y = occupancy, 
                 color = energy.barrier , 
                 shape = Segmentation), size =3)+ 
  facet_wrap(~enzyme)


# Calculating a confusion matrix for the two models

segmented.local.MSD22 <- segmented.local.MSD2
library(plyr)
segmented.local.MSD2$energy.barrier <- 
  revalue(segmented.local.MSD2$energy.barrier, 
          c("Ea < 2KbT" = "2", "Ea > 2KbT"= "1"))

detach("package:plyr", unload=TRUE)


confusionMatrix(as.factor(segmented.local.MSD2$energy.barrier), 
                as.factor(segmented.local.MSD2$Markov.state)) 


# Calculating average diffusion rates of the states 

segmented.local.MSD2 %>% filter(local.msd.05!=0, displacement.filter, on.dna.filter)%>% 
  group_by(enzyme, Markov.state) %>% 
  summarise(ave.diffusion.rate = mean(local.msd.05, na.rm = T))


############################ 7. Salt Dependence and hopping ###################################

segmented.local.MSD2 <- readRDS(segmented.data.path)


## 2-state
segmented.local.MSD2$energy.barrier <- 
  factor(segmented.local.MSD2$energy.barrier,
         levels = c('Ea > 2KbT', 'Ea < 2KbT'))

# salt dependence based on Markov states 

segmented.local.MSD2 %>% 
  ungroup() %>% 
  filter(local.msd.05!=0, displacement.filter, salt.concentration %in% c(10,50,100), on.dna.filter) %>% 
  group_by(enzyme, salt.concentration, Markov.state ) %>% 
  summarise(meanMSD = mean(local.msd.05), 
            errorMSD= sd(local.msd.05)/sqrt(n())) %>% 
  ggplot(aes(x = salt.concentration, y = meanMSD, color = as.factor(Markov.state)))+ 
  geom_point()+
  geom_errorbar(aes(x= salt.concentration, 
                    ymax = meanMSD + errorMSD, 
                    ymin = meanMSD - errorMSD,
                    color = as.factor(Markov.state))) +
  facet_wrap(~enzyme)+ 
  scale_x_log10()
  
# salt dependence based on Energy barrier states 


segmented.local.MSD2 %>% 
  ungroup() %>% 
  filter(local.msd.05!=0, large.displacement.filter, salt.concentration %in% c(10,50,100), on.dna.filter) %>% 
  group_by(enzyme, salt.concentration, energy.barrier ) %>% 
  summarise(meanMSD = mean(local.msd.05), 
            errorMSD= sd(local.msd.05)/sqrt(n())) %>% 
  ggplot(aes(x = salt.concentration, y = meanMSD, color = as.factor(energy.barrier)))+ 
  geom_point()+
  geom_errorbar(aes(x= salt.concentration, 
                    ymax = meanMSD + errorMSD, 
                    ymin = meanMSD - errorMSD,
                    color = energy.barrier)) +
  facet_wrap(~enzyme)+ 
  scale_x_log10()



  
################## 8.Stepping analysis #############################################

real.data <- readRDS(collected.data.path)

levels(real.data$enzyme)

# number of trajectories 

real.data %>% 
  filter(on.dna.filter, displacement.filter) %>%
  group_by(enzyme) %>% 
  summarise(length(trajectory.unique.id %>% unique()))

#  binding lifetime 

real.data %>% 
  filter(salt.concentration %in% c(10,50,100), on.dna.filter, large.displacement.filter) %>% 
  group_by(enzyme, salt.concentration, trajectory.unique.id) %>%
  summarize(number.of.frames = n(), 
            duration = number.of.frames * unique(frame.interval)/1000) %>% 
  group_by(enzyme, salt.concentration) %>% 
  summarise(binding.life.time = mean(duration), 
            standard.error = sd(duration)/sqrt(length(duration))) %>% 
  ggplot()+ 
  geom_point(aes(x= salt.concentration, y= binding.life.time, color = enzyme), 
             size = 1.2)+ 
  geom_errorbar(aes(x = salt.concentration, 
                                  ymin = binding.life.time - standard.error, 
                                  ymax = binding.life.time + standard.error, 
                    color = enzyme), width = 0.1) +
  scale_x_log10(breaks = c(1,10,50,100,150))



# average binding lifetime and displacement per trajectory 

real.data %>% 
  filter(salt.concentration>0.1, on.dna.filter, displacement.filter) %>% 
  group_by(enzyme, trajectory.unique.id) %>%
  summarize(number.of.frames = n(), 
            duration = number.of.frames * unique(frame.interval)/1000, 
            displacement = (max(corrected_X)- min(corrected_X))/0.34) %>% 
  group_by(enzyme) %>% 
  summarise(mean(duration), mean(displacement))
  

#  Scanning speed distribution  

real.data%>%
  filter(on.dna.filter, displacement.filter) %>% 
  group_by(trajectory.unique.id) %>%
  mutate(Step=abs((lead(corrected_X, n=1)- corrected_X)/frame.interval)) %>% 
  group_by(enzyme) %>%
  ggplot(aes(x=Step, y=..density..))+
  geom_histogram( bins = 300)+ 
  facet_wrap(~enzyme,ncol = 3)+ 
  coord_cartesian(xlim =c(0,50)) + 
  xlab("Scanning speed (nm/ms) ") + 
  ylab("Density") 

real.data%>%
  filter(on.dna.filter, displacement.filter) %>% 
  group_by(trajectory.unique.id) %>%
  mutate(Step=abs((lead(corrected_X, n=1)- corrected_X)/frame.interval)) %>% 
  filter(Step<50) %>% 
  group_by(enzyme) %>% 
  ggplot()+ 
  geom_density(aes(x=Step , color = enzyme), size= 1)+ 
  coord_cartesian(xlim =c(0,50))+ 
  xlab("Scanning speed (nm/ms) ") + 
  ylab("Density") 

# Redundancy-efficiency analysis 

real.data %>%  filter(on.dna.filter, displacement.filter) %>% 
  group_by(trajectory.unique.id) %>% 
  mutate(duration = (max(frame.number)-min(frame.number))*frame.interval,  
         speed = (lead(corrected_X, n=1)- corrected_X)/frame.interval, 
         walked.length = mean(abs(speed), na.rm = T) *duration,
         covered.length = max(corrected_X)-min(corrected_X), 
         walking.speed = mean(abs(speed), na.rm = T),
         covering.speed = (max(corrected_X)-min(corrected_X))/duration) %>% 
  group_by(enzyme) %>% 
  summarise(average.walked.length = mean(walked.length, na.rm = T), 
            average.covered.length = mean(covered.length, na.rm = T),
            average.walking.speed = mean(walking.speed, na.rm = T), 
            average.covering.speed = mean(covering.speed, na.rm = T)) %>% 
  ggplot()+
  geom_point(aes(x = average.covering.speed/0.34, 
                 y = average.walked.length/average.covered.length, 
                 color = enzyme), size = 2)+
  scale_color_manual(values = c("magenta",
                                "dark blue", "dodger blue", "light sky blue",
                                "seagreen","green"))+
  ylab("Redundancy") + 
  xlab("Efficiency (bp/ms)") 


##### Data file paths ##### 

filtered.data.path <- "processed/1_AlkD_filtered_data.rds"
noise.excluded.path <- "processed/2_AlkD_detected_noise_excluded.rds"
transformed.data.path <- "processed/3_AlkD_AlkF_transfromed_data.rds"
collected.data.path <- "processed/4_collected_data_2020_03_04.rds"
local.msd.path <- "processed/6_local_MSD_2020_03_04.rds"
simulated.data.path <- "Processed/7_local_MSD_sim_2020_03_04.rds"
segmented.data.path <- "processed/9_segmented_2_states_local_MSD_2020_04_08.rds"




##### Fig. 1 ##### 

collected.data <- readRDS(collected.data.path)


theme_white2 <- theme(panel.background = element_blank(),
                      legend.key = element_blank(),
                      legend.background = element_blank(),
                      strip.background = element_blank(),
                      plot.background = element_blank(),
                      panel.grid = element_blank(),
                      strip.text = element_blank(),
                      axis.text =element_text(size = 6,color = "black"), 
                      axis.line = element_blank(),
                      axis.ticks = element_line(color = "black"),
                      axis.title.y = element_text(size = 6,color = "black"),
                      axis.title.x = element_text(size = 6,color = "black"), 
                      panel.border = element_rect(colour = "black", fill = NA, 
                                                  size= 1))

collected.data %>% filter(trajectory.unique.id == "AlkD_2691") %>% 
  mutate(time = (frame.number - min(frame.number))*frame.interval, 
         displacement = corrected_X- corrected_X[1]) %>% 
  ggplot(aes(x = displacement/1000, y = time/1000))+
  geom_path(color = "red", alpha = 0.5)+ geom_point(color= "red", size = 0.8)+ 
  theme_white2+ ylab("Time (s)") +
  xlab(expression(Displacement~along~DNA~(mu*m)))+ 
  scale_x_continuous(position = "top", limits = c(-1.0,1.1))+ 
  scale_y_reverse()

ggsave(filename = "Fig.1.svg" ,path= "../../Figures/" , 
       width=8, height= 6, units= "cm", dpi=300)






##### Fig. 2d #####

segmented.local.MSD <- readRDS(segmented.data.path)
local.MSD.sim <- readRDS(simulated.data.path)



theme_white4 <- theme(panel.background = element_blank(),
                      legend.key = element_blank(),
                      legend.background = element_blank(),
                      strip.background = element_blank(),
                      plot.background = element_blank(),
                      axis.line = element_line(color = "black"),
                      axis.ticks = element_line(color = "black"),
                      strip.text = element_text(size = 8, color = "black"),
                      axis.title.y = element_text(size = 7,color = "black"),
                      axis.title.x = element_text(size = 7,color = "black"),
                      axis.text = element_text(color = "black", size= 7),
                      legend.position="top", 
                      legend.direction = "horizontal",
                      legend.text = element_text(size=7), 
                      legend.title = element_text(size=7), 
                      panel.spacing = unit(1, "mm"))


segmented.local.MSD%>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.05!=0) %>% 
  group_by(enzyme) %>%
  ggplot(aes(x=local.msd.05)) +
  facet_wrap(~enzyme, ncol = 3) +
  geom_histogram(bins=99, position = "dodge",
                 aes(y=21*(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..], 
                     fill= as.factor(Markov.state)))  +
  scale_x_log10(breaks= c(0.01,0.1,1,10))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm")) + 
  geom_density(color = "red")+
  ylab("Density")+ 
  geom_density(data =local.MSD.sim %>% filter(local.msd.05!=0,
                                              enzyme %in% c("AlkD", "AlkF", "mut-AlkF")) ,
               aes(x=local.msd.05), 
               color = "black")+
  geom_density(data =local.MSD.sim %>% filter(local.msd.05!=0,
                                              enzyme == "AlkF") %>% mutate(enzyme = "mut-AlkF"),
               aes(x=local.msd.05), linetype="dashed")+
  geom_segment(data=data.frame(x=c(0.49/exp(2), 0.8/exp(2), 0.8/exp(2)), 
                               y= c(0.76, 0.44, 0.44), 
                               enzyme=c("AlkD", "AlkF", "mut-AlkF")), 
               aes(x= x, y= 0, xend= x ,yend=y), inherit.aes=FALSE, linetype="11")+
  xlab( expression(Instantaneous~diffusion~rate~(mu*m^2/s))) +
  scale_y_continuous(breaks = c(0.5, 1), limits = c(0,1.25))+
  annotation_logticks(side= "b",
                      short = unit(0.3,"mm"), 
                      mid = unit(0.6,"mm"), 
                      long = unit(1,"mm"))+ 
  scale_fill_manual(values = c("#D55E00", "#0072B2"),
  name = "Markov classification \n of protein trajectories", 
                        labels= c("Slow", "Fast"))+ 
  theme_white4

ggsave(filename = "Fig.2d.svg" ,path= "../../Figures/" , 
       width=17.5, height= 6, units= "cm", dpi=300)

##### Fig. 2d inset ##### 

energy <- readRDS("processed/energy.differnce.data.rds")

theme_white0 <- theme(panel.background = element_blank(),
                      legend.key = element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.background = element_blank(),
                      strip.background = element_blank(),
                      plot.background = element_blank(),
                      panel.grid = element_blank(),
                      strip.text = element_blank(),
                      axis.text.y = element_blank(),
                      axis.text.x =element_text(size = 6,color = "black"), 
                      axis.line = element_blank(),
                      axis.ticks = element_line(color = "black"),
                      axis.title.x = element_text(size = 6,color = "black"),
                      axis.title.y = element_blank(), 
                      panel.border = element_rect(colour = "black", fill = NA, 
                                                  size= 1))


energy %>% ggplot()+ 
  geom_segment(aes(x = 0, y = xi, xend = energy.difference, 
                   yend = xi), lineend = "round", linejoin = "round", size =0.5, color = "gray")+  
  geom_text (aes(x = 0.5, y =xi + 0.1  , label = enzyme), size =2) + 
  geom_segment(aes(x= 1.5, xend = 1.5, y= 0, yend = 0.9), linetype='dotted', col = 'red', size = 0.5)+
  geom_segment(aes(x = energy.difference, y= xi -0.07, xend= energy.difference, yend = xi + 0.07), 
               lineend = "round", linejoin = "round", size =0.5)+
  xlim(c(0, 2.8)) + ylim(c(0,0.9)) + theme_white0 + xlab(expression(Energy~barrier~difference(~k[B]*T)))


ggsave(filename = "Fig.2d.inset.svg" ,path= "../../Figures/" , 
       width=4, height= 2.5, units= "cm", dpi=300)


##### Fig. 2e ##### 

segmented.local.MSD <- readRDS(segmented.data.path)

markov.result <- segmented.local.MSD %>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.05!=0, !is.na(Markov.state)) %>% 
  group_by(enzyme, Markov.state) %>% 
  summarise(diffusion.rate = mean(local.msd.05), 
            diffusion.rate.error = sd(local.msd.05)/sqrt(n()),
            occupancy = n()) %>% 
  group_by(enzyme) %>% 
  mutate(l= sum(occupancy)) %>% 
  mutate(occupancy = occupancy/l) %>% 
  select(-l) %>%
  mutate(Segmentation = "Markov states")
  

energy.result <- segmented.local.MSD %>% 
  filter(displacement.filter, on.dna.filter) %>% 
  filter(local.msd.05!=0, !is.na(Markov.state)) %>% 
  group_by(enzyme, energy.barrier) %>% 
  summarise(diffusion.rate = mean(local.msd.05), 
            diffusion.rate.error = sd(local.msd.05)/sqrt(n()), 
            occupancy = n()) %>% 
  group_by(enzyme) %>% 
  mutate(l= sum(occupancy)) %>% 
  mutate(occupancy = occupancy/l) %>% select(-l) %>%
  mutate(Segmentation = "Energy barrier")


theme_white4 <- theme(panel.background = element_blank(),
                      legend.key = element_blank(), 
                      legend.background = element_blank(),
                      strip.background = element_blank(), 
                      plot.background = element_blank(), 
                      panel.grid.major = element_line(color = "gray"),
                      axis.line = element_line(color = "black"),
                      axis.ticks.y = element_line(color = "black"), 
                      axis.ticks.x = element_blank(), 
                      strip.text = element_blank(), 
                      axis.title.y = element_text(size = 7, color = "black"),
                      axis.title.x = element_text(size = 7, color = "black"), 
                      axis.text.x = element_text (size=7, color= "black"), 
                      axis.text.y = element_text (size=7, color= "black"), 
                      legend.position="top", legend.direction = "horizontal", 
                      legend.text = element_text(size=7), 
                      legend.title = element_text(size=7), panel.spacing = unit(10, "mm"))
##plot
markov.result %>% 
  ggplot() + geom_point(aes(x = diffusion.rate, y = occupancy, color = as.factor(Markov.state) , 
                            shape = Segmentation), size =3)+ 
  geom_point(data = energy.result, aes(x = diffusion.rate, y = occupancy, color = energy.barrier , 
                                       shape = Segmentation), size =3)+
  facet_wrap(~enzyme)+ 
  facet_wrap(~enzyme,ncol=3, scales = "free_y")+ 
  theme_white4+ scale_x_log10(breaks= c(0.05, 0.1, 0.2, 0.4, 0.8))+ 
  ylab("Mode occupancy")+ 
  xlab(expression(Average~diffusion~rate~(mu*m^2/s))) + 
  annotation_logticks(side= "b", 
                      short = unit(0.3,"mm"), mid = unit(0.6,"mm"), 
                      long = unit(1,"mm"))+ 
  scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9", "#E69F00"))+
  scale_shape_discrete(name = "Classification method")+ 
  scale_y_continuous(limits=c(0,0.8), breaks = c(0.25,0.5,0.75))

ggsave(filename = "Fig.2e.svg" ,path= "../../Figures/" , 
       width=17.5, height= 6, units= "cm", dpi=300)


##### Fig. 3a ##### 

segmented.local.MSD <- readRDS(segmented.data.path)

segmented.local.MSD$energy.barrier <- 
  factor(segmented.local.MSD$energy.barrier, levels = c('Ea > 2KbT', 'Ea < 2KbT'))


theme_white6 <- theme(panel.background = element_blank(),
                      legend.key = element_blank(),
                      legend.background = element_blank(),
                      strip.background = element_blank(),
                      plot.background = element_blank(),
                      panel.grid = element_blank(),
                      axis.line = element_line(color = "black"),
                      axis.ticks = element_line(color = "black"),
                      strip.text = element_text(size = 8,
                                                color = "black"),
                      axis.title.y = element_text(size = 7,
                                                  color = "black"),
                      axis.title.x = element_text(size = 7,
                                                  color = "black"),
                      axis.text = element_text(color = "black",
                                               size= 7),
                      legend.text = element_text(size=7), 
                      legend.position = "top", 
                      legend.title = element_text(size=7), 
                      panel.spacing = unit(3, "mm"))

segmented.local.MSD %>% 
  ungroup() %>% 
  filter(local.msd.05!=0, 
         displacement.filter, 
         on.dna.filter,
         salt.concentration !=150,
         salt.concentration != 1) %>% 
  group_by(enzyme, salt.concentration, energy.barrier ) %>% 
  summarise(meanMSD = mean(local.msd.05), 
            errorMSD= sd(local.msd.05)/sqrt(n())) %>% 
  ggplot(aes(x = salt.concentration, y = meanMSD, 
             color = as.factor(energy.barrier)))+ 
  geom_point()+
  geom_errorbar(aes(x= salt.concentration, 
                    ymax = meanMSD + errorMSD, 
                    ymin = meanMSD - errorMSD,
                    color = as.factor(energy.barrier))) +
  facet_wrap(~enzyme, scales = "free_y")+ 
  xlab("Salt concentration (mM)")+
  ylab(expression(Average~diffusion~rate~(mu*m^2/s)))+ 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), 
                     name = "Energy barrier", 
                      labels= c(expression(E[a]>~2*~k[B]*T), 
                                expression(E[a]<2*~k[B]*T)))+ 
  annotation_logticks(side= "b",short = unit(0,"mm"),
                      mid = unit(0.6,"mm"), long = unit(1,"mm"))+ 
  scale_x_log10()+
  theme_white6

ggsave(filename = "Fig.3a.svg" ,path= "../../Figures/" , 
       width=8, height= 6, units= "cm", dpi=300)


##### Fig. 3b ##### 

real.data <- readRDS(collected.data.path)

real.data$enzyme <- 
  factor(real.data$enzyme, levels = c("hOGG1", "AlkD", "AlkF", "mut-AlkF", "EndoV", "mut-EndoV" ))


theme_white7 <- theme(panel.background = element_blank(),
                      legend.key = element_blank(),
                      legend.background = element_blank(),
                      strip.background = element_blank(),
                      plot.background = element_blank(),
                      axis.line = element_line(color = "black"),
                      panel.grid.major = element_line(color = "gray"),
                      axis.ticks = element_line(color = "black"),
                      strip.text = element_text(size = 8, color = "black"),
                      axis.title.y = element_text(size = 7,color = "black"),
                      axis.title.x = element_text(size = 7,color = "black"),
                      axis.text = element_text(color = "black", size= 7),
                      legend.position="none", 
                      legend.direction = "vertical",
                      legend.text = element_text(size=6), 
                      legend.title = element_text(size=7), 
                      panel.spacing = unit(1, "mm"))




real.data %>%  filter(on.dna.filter, displacement.filter) %>% 
  group_by(trajectory.unique.id) %>% 
  mutate(duration = (max(frame.number)-min(frame.number))*frame.interval,  
         speed = (lead(corrected_X, n=1)- corrected_X)/frame.interval, 
         walked.length = mean(abs(speed), na.rm = T) *duration,
         covered.length = max(corrected_X)-min(corrected_X), 
         walking.speed = mean(abs(speed), na.rm = T),
         covering.speed = (max(corrected_X)-min(corrected_X))/duration) %>% 
  group_by(enzyme) %>% 
  summarise(average.walked.length = mean(walked.length, na.rm = T), 
            average.covered.length = mean(covered.length, na.rm = T),
            average.walking.speed = mean(walking.speed, na.rm = T), 
            average.covering.speed = mean(covering.speed, na.rm = T)) %>% 
  ggplot()+
  geom_point(aes(x = average.covering.speed/0.34, 
                 y = average.walked.length/average.covered.length, 
                 color = enzyme), size = 2)+
  scale_color_manual(values = c("magenta",
                                "dark blue", "dodger blue", "light sky blue",
                                "seagreen","green"))+
  ylab("Redundancy") + 
  xlab("Efficiency (bp/ms)") + 
  theme_white7 


ggsave(filename = "Fig.3b.svg" ,path= "../../Figures/" , 
       width=8, height= 5, units= "cm", dpi=300)



