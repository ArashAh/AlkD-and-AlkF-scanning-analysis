#####
# construction of the matrix from outputs of thunderSTORM

MakeItMatrix <- function(RawData,InputData)
{
  # creation of a matrix from raw data, this matrix can be used
  # for visualization of the raw data parallel to the filtered
  # and detected data
  row = max(RawData$frame) # max number of frame 
 
  col <- InputData %>% 
    group_by(frame.number) %>% 
    summarise(nSignals = n()) %>% 
    select(nSignals) %>% 
    max()
  
  DataMatx <- matrix(0,row,col)
  DataMaty <- matrix(0,row,col)
  DataMatf <- matrix(0,row,1)
  k = 0
  p = 0
  
  for(i in 1:row) 
  {
    DataMatf[i,1] <- i
    
    l <- length(InputData$frame.number[InputData$frame.number==i])
    if (l==0)
    {
      DataMatx[i,] <- 0
      DataMaty[i,] <- 0
      k <- k+1
    }
    else 
    {
      for(j in 1:l)
      {
        DataMatx[i,j] <- InputData$X[i-k+p+j-1]
        DataMaty[i,j] <- InputData$Y[i-k+p+j-1]
        
      }
      p <- p+l-1
    }
  }
  DataMat <- data.frame(DataMatx, DataMaty)
  for (i in 1:col)
  {
    names(DataMat)[i]<-paste("x",i, sep = "")
    names(DataMat)[col+i]<-paste("y",i, sep = "")
  }
  return(DataMat)
}

#####



#####
# Finding the trajectories 

FindTrajectory <- function(FilteredDataSet, dxmax, dxmin, dymax)
  {
  info <- 0
  row <- nrow(FilteredDataSet)
  col <- ncol(FilteredDataSet)/2
  # detection loop 
  for(i in 1:(row-1)) 
  {
    for(j in 1:(col))
    {
      if (FilteredDataSet[i,j]!=0)
      {
        if(j!=col)
        {
          if (sqrt((FilteredDataSet[i,j])^2 + 
                   (FilteredDataSet[i,(j+col)])^2) -
              sqrt((FilteredDataSet[i,(j+1)])^2 + 
                   (FilteredDataSet[i,(j+col+1)])^2) < 200 &
              abs(FilteredDataSet[i,j]- FilteredDataSet[i,(j+1)]) <200) 
            # if two signal has been detected within one PSF in one frame, 
            #consider the mean of the two 
          {
            FilteredDataSet[i,j] <- (FilteredDataSet[i,j]+
                                       FilteredDataSet[i,(j+1)])/2
            FilteredDataSet[i,j+col] <- (FilteredDataSet[i,j+col]+
                                           FilteredDataSet[i,(j+col+1)])/2
            FilteredDataSet[i,j+1] <- 0
            FilteredDataSet[i,j+col+1] <- 0
          }
        }
        deltax<- abs(FilteredDataSet[i,j]-FilteredDataSet[i+1,1:col])   
        deltay<- abs(FilteredDataSet[i,(j+col)]-
                       FilteredDataSet[i+1,(col+1):(2*col)])   
        min <- which.min(deltax) 
        # guarantees that we follow the closest signali in the next frame 
        
        if(deltax[min]<dxmax & deltax[min]>dxmin & deltay[min]<dymax)
        {
          if (i!=1)
          {
            if (FilteredDataSet[i-1,j]!=0) 
            {# this is for separation of two consacutive trajectories 
              if (abs(FilteredDataSet[i,j]-FilteredDataSet[i-1,j])> dxmax |
                  abs(FilteredDataSet[i,j]-FilteredDataSet[i-1,j])< dxmin |
                  abs(FilteredDataSet[i,(j+col)]-
                      FilteredDataSet[i-1,(j+col)])> dymax)
              { # basically, if this is the first point 
                #then make the previous point zero, make sure it is 
                #not following another trajectory 
                
                # print("first point detected")
                # print(i)
                # print(j)
                
                FilteredDataSet[i-1,j] <- 0
                FilteredDataSet[i-1,(j+col)] <- 0
              }
            }
          }
          if (abs(FilteredDataSet[i+1,min]-FilteredDataSet[i,j]) ==
              min(abs(FilteredDataSet[i,1:col]-FilteredDataSet[i+1,min])))
          { # in case there are two signal within the dxmax and dymax, 
            #this help to choose the clsoest signal to the signal in 
            #the last frame 
            ax <- FilteredDataSet[i+1,j]
            FilteredDataSet[i+1,j] <- FilteredDataSet[i+1,min]  
            FilteredDataSet[i+1,min] <- ax
            ay <- FilteredDataSet[i+1,(j+col)]
            FilteredDataSet[i+1,(j+col)] <- FilteredDataSet[i+1,(min+col)]  
            FilteredDataSet[i+1,(min+col)] <- ay
          } 
          else
          {
            # print("it didn't move")
            # print(i)
            # print(j)
          }
        } 
        else if (i!=1)
        {
          if (FilteredDataSet[i-1,j]!=0)
          {
            if (abs(FilteredDataSet[i,j]-FilteredDataSet[i-1,j])<dxmax &
                abs(FilteredDataSet[i,j]-FilteredDataSet[i-1,j])> dxmin &
                abs(FilteredDataSet[i,(j+col)]-
                    FilteredDataSet[i-1,(j+col)])<dymax)
            { # this is how the last point of the trajectory is checked and 
              #preserved 
              # print("last point or single point detected")
              # print(i)
              # print(j)
            }
            else 
            {
              FilteredDataSet[i,j] <- 0
              FilteredDataSet[i,(j+col)] <- 0
            }
          }
          else 
          {
            FilteredDataSet[i,j] <- 0
            FilteredDataSet[i,(j+col)] <- 0
          }
        }
        else 
        {
          FilteredDataSet[i,j] <- 0
          FilteredDataSet[i,(j+col)] <- 0
        }
      }
      else 
      {  
        FilteredDataSet[i,j] <- 0
        FilteredDataSet[i,(j+col)] <- 0
      }
    }
  }
  return(FilteredDataSet)
}

#####




##### 

MakeZeroNA <- function(InputDataSet)
{
  InputDataSet[InputDataSet==0] <- NA
  return(InputDataSet)
}

##### 




##### 

# Make long form 

MakeLongForm <- function(InputMatrix) {
  LongForm <- data_frame(frame.number = rep(1:(nrow(InputMatrix)), (ncol(InputMatrix)/2)), 
                         X = InputMatrix[, 1:(ncol(InputMatrix)/2)] %>%
                           unlist() %>% 
                           as.vector(), 
                         Y = InputMatrix[, (ncol(InputMatrix)/2 + 1):
                                           ncol(InputMatrix)] %>% 
                           unlist() %>% 
                           as.vector(), 
                         track.register = rep(LETTERS[1:(ncol(InputMatrix)/2)], 
                                 each = nrow(InputMatrix)))
  return(LongForm)
  
}
#####




#####

CorrectBlinking <- function(InputData){
  output <- InputData %>% 
  mutate(pick = ifelse(X == 0 & lead(X, n=1)!= 0 & lag(X, n=1)!=0, 1, 0)) %>% 
  mutate(setx = ifelse(pick ==1 & abs(lead(X, n=1)-lag(X, n=1)) < 600 &
                         abs(lead(Y, n=1)-lag(Y, n=1)) < 600, 
                       (lead(X, n=1)+ lag(X, n=1))/2,0),
         sety= ifelse(pick ==1 & abs(lead(X, n=1)-lag(X, n=1)) < 600 &
                        abs(lead(Y, n=1)-lag(Y, n=1)) < 600, 
                      (lead(Y, n=1) + lag(Y, n=1))/2,0)) %>% 
  mutate(X = ifelse(setx == 0 , X, setx),
         Y = ifelse(sety == 0 , Y, sety), 
         blinked =ifelse(setx != 0 & sety != 0, "Yes", "No")) %>% 
  select(-c(pick, setx, sety))
}
#####

add.col<-function(df, Markov.state) {n.row<-dim(df)[1]
length(Markov.state)<-n.row
cbind(df, Markov.state)
}






#####

ExcludeNoise <- function(DetectedTrajectories){

  output <- DetectedTrajectories %>%
  mutate(see = ifelse(is.na(X), 1,0), see2 = cumsum(see)) %>% 
  mutate(trajectory.id = ifelse(see == 1, 0,see2)) %>% 
  group_by(trajectory.id) %>% 
  mutate(delta.t = max(frame.number)-min(frame.number), 
         delta.x = max(X)-min(X)) %>% 
  mutate(duration.filter = ifelse(delta.t >= 5 & !is.na(X),"Yes", "No"),
         displacement.filter = ifelse(delta.x > 300 & !is.na(X), "Yes", "No"),
         large.displacement.filter = ifelse(delta.x > 600 &
                                              !is.na(X), "Yes", "No")) %>% 
  ungroup() %>% filter(!is.na(X)) %>%
  select(-c(see, see2, delta.t, delta.x))

new.id <- output %>%  
  group_indices(trajectory.id) 

output$trajectory.id <- new.id

output <- left_join(output, raw.data.set) %>% 
  select(data.set.name, trajectory.id, frame.number:Y, 
         duration.filter:large.displacement.filter, intensity)

output$data.set.name[is.na(output$data.set.name)] <- output$data.set.name[1]
output$intensity[is.na(output$intensity)] <- 0

return(output)
}

#####







##### 

MakeTrajectoryAddress <- function(Input){
output <- Input %>% 
  filter(large.displacement.filter == "Yes", duration.filter == "Yes") %>% 
  group_by(trajectory.unique.id) %>% 
  summarise(file.name = data.set.name %>% unique(),
            start.frame = min(frame.number), 
            end.frame = max(frame.number))
}

#####







#####

TransformData <- function(Input){
  Output <- Input %<>% 
  mutate(n = data.set.name) %>% 
  separate(n, c("enzyme","salt.concentration","date", 
                "frame.interval","dilution", "id"), sep = "_") %>% 
  select(-date, -dilution, -id) %>% 
  mutate(frame.interval = str_replace(frame.interval, "mspf", "") %>% 
           as.numeric()) %>% 
  mutate(salt.concentration = str_replace(salt.concentration, 
                                          "mMNaCl", ""),
         salt.concentration = str_replace(salt.concentration, 
                                          "mMKCl", ""), 
         salt.concentration = as.numeric(salt.concentration)) %>% 
  unite(trajectory.unique.id, data.set.name, trajectory.id, remove = FALSE) %>% 
  mutate(trajectory.unique.id = as.factor(trajectory.unique.id) %>% 
           as.numeric()) 
}

TransformData2 <- function(Input){
  Output <- Input %<>% 
    mutate(n = data.set.name) %>% 
    separate(n, c("enzyme","del1","del2","salt.concentration","del3","frame.interval"), sep = "_") %>% 
    select(-del1, -del2, -del3)  %>% 
    mutate(salt.concentration = str_replace(salt.concentration, "mM", "") %>% 
           as.numeric()) %>% 
    mutate(frame.interval = str_replace(frame.interval, "ms", "") %>% 
             as.numeric()) %>%  
    unite(trajectory.unique.id, data.set.name, trajectory.id, remove = FALSE) %>% 
    mutate(trajectory.unique.id = as.factor(trajectory.unique.id) %>% 
             as.numeric()) 
}
#####






#####

OverviewData <- function(Input){
  
  Output1 <- Input %>% 
    group_by(enzyme, salt.concentration) %>% 
    summarise(number.of.data.sets = length(unique(data.set.name))) 
  
  write.table(Output1,"../owerview_tracking_output/Overview1.txt")

  Output2 <- Input %>% 
    group_by(enzyme, duration.filter, displacement.filter, 
             large.displacement.filter) %>% 
    summarise(nubmer.of.trajectories = length(unique(trajectory.unique.id))) 

  write.table(Output2,"../owerview_tracking_output/Overview2.txt")
    
  Output3 <- Input %>% 
    filter(duration.filter == "Yes") %>% 
    group_by(enzyme, salt.concentration, large.displacement.filter) %>% 
    summarise(number.of.trajectories = length(unique(trajectory.unique.id))) 
  
  write.table(Output3,"../owerview_tracking_output/Overview3.txt")
  
  return(c(print(Output1),print(Output2),print(Output3)))
}

#####






#####

RotatePoints <- function(xPos, yPos, estimate) {
  coordinates <- matrix(c(xPos, yPos), ncol = 1)
  mat <- matrix(c(cos(atan(estimate)), sin(atan(estimate)), 
                  -sin(atan(estimate)), cos(atan(estimate))), ncol = 2)
  (mat %*% coordinates) %>% as.vector()
}

#####















