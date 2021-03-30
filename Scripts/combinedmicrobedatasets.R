library(tidyverse)

#load data to merge

FCM_Fixed_clean <- read_csv("Data/Nyssa_FCM_Fixed_clean.csv")
Tidepools_fDOM_Final <- read_csv("Data/Silbiger_Tidepools_fDOM_Final.csv")

combinedFCMandfDOMdata<-left_join(FCM_Fixed_clean,Tidepools_fDOM_Final)

write.csv(combinedFCMandfDOMdata, file = 'Data/combined_FCMandfDOMdata.csv')
