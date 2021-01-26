#Biogeochem for microbial samples data processing
## By: Jenn Fields
## Last updated: 12.23.2020
###########################
## clear workspace
rm(list=ls())
####load libraries####
library(seacarb) # for carbonate chem
library(oce) # for carbonate chem
library(vegan) #for multi variate analysis
library(devtools)
library(AUC)
library(factoextra)
library(plotrix) #for SE
library(MASS) #for stat analysis assumptions
library(car)#for stat analysis assumptions
library(lubridate) #for data/time
library(hms) #date/time
library(ggrepel) #for pcas
library(heplots)
library(zoo) #rolling averages
library(broom) #for tidy function
library(tidyverse) #for all things %>%

# load data
#source scripts
source("scripts/tidepoolphysicalparameters.R")#load physical parameters

#Load Data
Nutrients<-read_csv("Data/Biogeochem/RawNutrientData.csv")
CarbChem<-read_csv("Data/Biogeochem/ChemData.csv")
Salinity<-read_csv("Data/Biogeochem/SampleSalinityData_Adjusted.csv")
TA<- read_csv("Data/Biogeochem/TASamples_Adjusted.csv")
## bring in pH calibration files
pHcalib<-read.csv('Data/Biogeochem/SummerTrispHcalibration.csv')

# First calculate salinity from conductivity
CarbChem$SalCal<-swSCTp(as.numeric(CarbChem$Conductivity)/1000, CarbChem$Temp.pool, pressure=rep(0, each=nrow(CarbChem)),
                        conductivityUnit="mS/cm")

#plot salinity straight from instrument versus calculated from conuctivity and temp
#plot(CarbChem$SalCal, CarbChem$Salinity, xlab='Salinity from conductivity', ylab="Salinity")
#looks very similar...good

#########Cleaning Light and Temp logger data#################
#Clean up time code:
#change date/time data for lubridate package
#first column date and time in same column

CarbChem$Sampling_Day<-mdy(CarbChem$Sampling_Day, quiet=FALSE, tz="America/Los_Angeles", truncated=0)
CarbChem$Date_Time <- paste(CarbChem$Sampling_Day, CarbChem$Sampling_time)
CarbChem$Date_Time<-ymd_hms(CarbChem$Date_Time, quiet=FALSE, tz="America/Los_Angeles", truncated=0) #format date and time

CarbChemAll<-CarbChem


####Microbe data sets (all time points)#####
#remove time 5 from after period since no microbes taken
CarbChemAll<-CarbChemAll %>%
  filter(!(Before_After =="After" & Time_Point == '5')) #filter out time 5 from after period because no microbes

#create data set with just the tp id and start and stop times
Start<- CarbChemAll %>%
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 1 ) %>%
  dplyr::mutate(StartTime = Date_Time[Time_Point== 1]) %>%
  dplyr::select(PoolID, StartTime)

Stopbefore<- CarbChemAll %>%
  filter(Foundation_spp != 'Ocean') %>%
  filter(Before_After =='Before') %>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 5) %>%
  dplyr::mutate(StopTime = Date_Time[Time_Point== 5]) %>%
  dplyr::select(PoolID, StopTime)

stopnightbefore<-CarbChemAll%>%
  filter(Foundation_spp != 'Ocean') %>%
  filter(Before_After =='Before' & Day_Night =="Night" & Sampling_Day == "2019-07-12") %>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 4) %>%
  dplyr::mutate(StopTime = Date_Time[Time_Point== 4]) %>%
  dplyr::select(PoolID, StopTime)

stopnightbefore<-rbind(Stopbefore,stopnightbefore)

Stopafter<-CarbChemAll %>%
  filter(Foundation_spp != 'Ocean') %>%
  filter(Before_After =="After")%>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 4) %>%
  dplyr::mutate(StopTime = Date_Time[Time_Point== 4]) %>%
  dplyr::select(PoolID, StopTime)

Stopall<-rbind(stopnightbefore,Stopafter)

StartStopall<-left_join(Start,Stopall) #combine both start and stop

temp.data <-
  list.files(path = 'Data/LightandTemp',pattern = ".csv", full.names = TRUE) %>% #all temp and light files in the is folder
  # list files in directory following a particular pattern
  set_names(.) %>% # get the column names
  map_dfr(read.csv, .id = "file.ID") %>% # join all files together in one data frame by file ID
  group_by(file.ID) 

temp.data$Date_Time<- mdy_hm(temp.data$Date.Time, quiet=FALSE, tz="America/Los_Angeles", truncated=0) #format date and time 
#View(temp.data)

#convert lux to par
#parameters from Long 2012 Field exp values
#A1 = -4924.7
#t1 = 20992.9
#y0 =  4929.0
#PARLICOR = A1 e^(–HOBO/t1) + y0
A1 <- -4924.7
t1 <- 20992.9
y0 <- 4929.0

#covert lux to numeric and get rids of pesky , in lux values
temp.data$Intensity.lux <- as.numeric(gsub(",","",temp.data$Intensity.lux))
x <- temp.data$Intensity.lux #create vector for just lux column
x[x == 0] <- NA #convert 0 to NA
#x[x >= 	19671] <- NA #get rid of values over 	19671 lux/30000PAR (not natural?)
#View(x)

# enter into equation: PARLICOR = A1e(–HOBO/t1) + y
#assign to new column for Par
temp.data$Par<- (A1* (exp( -x / t1))) + y0

temp.data$Par[is.na(temp.data$Par)]<- 0
temp.data$PoolID<-as.character(temp.data$PoolID)

StartandStopall<-left_join(StartStopall,temp.data) #join with temp data

TempandLightSumall<-StartandStopall %>% #create new data frame finding temp and light at start and stop time of water sampling
  dplyr::filter(!(Date_Time >= StopTime | Date_Time <= StartTime)) %>%
  dplyr::group_by(PoolID,Before_After,Day_Night,by60=cut(Date_Time, "60 min")) %>% #hourly means of pools per ~time point
  dplyr::summarise(Temp.max = max(Temp.C, na.rm=T), 
                   Temp.mean = mean(Temp.C, na.rm = T),
                   Par.mean = mean(Par, na.rm = T),
                   Par.max = max(Par,na.rm = T),
                   Par.sum = sum(Par,na.rm = T))

###########Nutrient Data Cleaning##########
#Missing 8 samples from nutritent analysis out of 612 samples
#N_17_5_BC was lost before taking nutrients \
# D_24_1_BC, N_24_1_BC, N2_Ocean_5_BC, D_1_2_AI, D_16_2_AI, D_4_3_AI, N_20_3_AI were missing
#from lab analyses
#We removed time points that nutrients were missing
#missing data was somewhat randomly distributed in dataframe

CarbChemAll<-left_join(CarbChemAll,Nutrients) #join with carbchem data

#Remove rows missing nutrient values from further analysis
CarbChemAll <- CarbChemAll %>%
  filter(Id_code != 'N_17_5_BC' & Id_code != 'D_24_1_BC' & Id_code != 'N_24_1_BC' &
           Id_code != 'N2_Ocean_5_BC' & Id_code != 'D_1_2_AI' & Id_code != 'D_16_2_AI'&
           Id_code != 'D_4_3_AI' & Id_code != 'N_20_3_AI') 

#Conversions from ug/l to µmol/l
#1 µg P/l = 1/MW P = 0.032285 µmol/l
#1 µg/l NO3 = 1/ MW NO3 µg/l = 0.016128 µmol/l
#1 µg/l NO2 = 1/ MW NO2 µg/l = 0.021736 µmol/l
# NO3 and NO2 = 0.016128 + 0.021736
#1 µg/l NH4 = 1/ MW NH4 µg/l = 0.055437 µmol/l

CarbChemAll$PO_ug_L<-as.numeric(CarbChemAll$PO_ug_L)
CarbChemAll$NN_ug_L<-as.numeric(CarbChemAll$NN_ug_L)
CarbChemAll$NH4_ug_L<-as.numeric(CarbChemAll$NH4_ug_L)
#Convert ug/L to umol/L for phosphate, nitrate, ammonium 
CarbChemAll$PO_umol_L <- CarbChemAll$PO_ug_L * 0.032285
CarbChemAll$NN_umol_L <- CarbChemAll$NN_ug_L * (0.016128 + 0.021736)
CarbChemAll$NH4_umol_L <- CarbChemAll$NH4_ug_L * 0.055437

#Converting negative values that were below detection to 0 values
CarbChemAll$PO_umol_L[which(CarbChemAll$PO_umol_L < 0)] <- 0
CarbChemAll$NN_umol_L [which(CarbChemAll$NN_umol_L  < 0)] <- 0
CarbChemAll$NH4_umol_L[which(CarbChemAll$NH4_umol_L < 0)] <- 0

####Adding TA data and lab salinity data####

CarbChemAll<-left_join(CarbChemAll,TA) #join with TA samples
CarbChemAll<-left_join(CarbChemAll,Salinity) #join with salinity in lab

###########CO2 Calculations##############
##correct the TA for the calibration error
CarbChemAll$TA_CRM<-CarbChemAll$TA_Raw-(2207.03*(CarbChemAll$CRM_off/100)) #compared to TA of batch #169 of CRM used for all samples

## take the mV calibration files by each date and use them to calculate pH and pH insitu
CarbChemAll<-pHcalib %>%
  dplyr::nest_by(date_triscal)%>%
  dplyr::mutate(fitpH = list(lm(mVTris~TTris, data = data))) %>% # linear regression of mV and temp of the tris
  dplyr::summarise(broom::tidy(fitpH)) %>% # make the output tidy
  dplyr::select(date_triscal, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  left_join(CarbChemAll,.) %>% # join with the pH sample data
  mutate(mVTris = Temp.in*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  mutate(pH = pH(Ex=mVolts,Etris=mVTris,S=Salinity,T=Temp.in)) %>% # calculate pH of the samples using the pH seacarb function
  mutate(pH_insitu = pHinsi(pH = pH, ALK = TA_CRM/1000000, Tinsi = Temp.pool, Tlab = Temp.in, S=Salinity, 
                            Pt=PO_umol_L/1000000, Sit=0, k1k2 = "x",  kf = "x", ks="d", pHscale = "T", b="u74")) 

###########DO regression#################
# The incorrect multiparameter sensor was used during Time point 1 of Day TP 1-16 in After period
# We were missing Dissolved Oxygen (DO mg/L and %) values for Ocean, TP 1-15
# We ran a regression of Dissolved Oxygen with pH and Temperature from the Before period to replace DO values
# Due to the strong correlation between DO, Temp, and pH
# We used the regression values to back calculate DO for the missing tide pools
# For the Ocean sample, we took the global average of the DO values of that day from time points 2-5 to replace
# the ocean value from time point 1

#Time point 1 for days (TP 1-32) from the Before period without ocean sample
DayTime1Chem<- CarbChemAll %>%
  filter(PoolID != "Ocean") %>%
  filter(Day_Night == "Day") %>%
  filter(Time_Point == "1") %>%
  filter(Before_After == "Before")
#View(DayTime1Chem)

DOTime1regression<-lm(DO_mg_L~pH_insitu + Temp.pool, data = DayTime1Chem) #DO mg/L,pH, and temp for time point 1 of days
summary(DOTime1regression)
#r2 of 0.9083
#B0 = -36.4751
#BpH =  5.4970  
#B temp  = 0.1504
DOpercentTime1regression<-lm(DO_p~pH_insitu + Temp.pool, data = DayTime1Chem) #%DO,pH, and temp for time point 1 of days
summary(DOpercentTime1regression)
#r2 of 0.9039
#B0 = -470.450
# BpH = 66.208 
# BTemp = 3.935

resid(DOTime1regression) #better to name these as something
DOregressionres<-resid(DOTime1regression)
#Now we can test the normality of the residuals
qqp(DOregressionres, "norm")
#normality
#We can also call the fitted y values as:
fitted(DOTime1regression)

#To test for homogeneity of variance, we want to plot the fitted (predicted) values
#against the residuals
plot(DOregressionres~fitted(DOTime1regression))
#good so can use model for DO values

#plot(DO~pH + Temperature, data=DayTime1Chem)
#abline(DOTime1regression, col="blue")

#select from dataset tide pools 1-15 for time point one on 08052019 sampling day during after period since missing DO values
tpday1AI<-which(CarbChemAll$Time_Point == "1" & CarbChemAll$Before_After == "After" & CarbChemAll$Sampling_Day == "2019-08-05" & 
                  CarbChemAll$PoolID != "Ocean" & CarbChemAll$PoolID != "16")
#View(tpday1AI)

#DO mg/L regression equation
#r2 of 0.9083
#B0 = -36.4751
#BpH =  5.4970  
#B temp  = 0.1504
CarbChemAll$DO_mg_L[tpday1AI]<- -36.4751 + (5.4970 * CarbChemAll$pH_insitu[tpday1AI]) + (0.1504 * CarbChemAll$Temp.pool[tpday1AI])
#View(CarbChem)

#DO % regression equation
#r2 of 0.9039
#B0 = -470.450
# BpH = 66.208 
# BTemp = 3.935
CarbChemAll$DO_p[tpday1AI]<- -470.450 + (66.208  *CarbChemAll$pH_insitu[tpday1AI]) + (3.935 * CarbChemAll$Temp.pool[tpday1AI])

#Now for ocean sample--select from data set time one on 0805209 sampling day
Oceanday1AI<-which(CarbChemAll$Time_Point == "1" & CarbChemAll$Before_After == "After" & CarbChemAll$Sampling_Day == "2019-08-05" & 
                     CarbChemAll$PoolID == "Ocean")
#View(Oceanday1AI)
CarbChemAll$DO_mg_L[Oceanday1AI]<- -36.4751 + (5.4970 * CarbChemAll$pH_insitu[Oceanday1AI]) + (0.1504 * CarbChemAll$Temp.pool[Oceanday1AI])
CarbChemAll$DO_p[Oceanday1AI]<- -470.450 + (66.208  * CarbChemAll$pH_insitu[Oceanday1AI]) + (3.935 * CarbChemAll$Temp.pool[Oceanday1AI])

########NEC and NEP rates##########
#calculate all the CO2Sys params
PoolCO2all<-carb(flag=8, CarbChemAll$pH_insitu, CarbChemAll$TA_CRM/1000000, S=CarbChemAll$Salinity, T=CarbChemAll$Temp.pool, Patm=1, P=0, 
                 Pt=CarbChemAll$PO_umol_L/1000000, Sit=0,
                 k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propagation
er<-errors(flag=8, CarbChemAll$pH_insitu, CarbChemAll$TA_CRM/1000000, 
           S=CarbChemAll$Salinity, T=CarbChemAll$Temp.pool, 
           Patm=1, P=0,Pt=CarbChemAll$PO_umol_L/1000000,
           Sit=0,evar1 = 0.01, evar2 = 5e-6) 

#average error for DIC based on pH and TA
mean(er$DIC*1000000)
#7.091473
sd(er$DIC*1000000)/sqrt(nrow(er))
#0.05899179

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
PoolCO2all[,c("CO2","HCO3","CO3","DIC","ALK")]<-PoolCO2all[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#combine all the data
CarbChemAll<-cbind(CarbChemAll,PoolCO2all[-c(1:6,10:13)])

#Normalize the TA and DIC to salinity

#TA normalized to constant salinity, but not nutrients here
CarbChemAll$TA_NormSal<-CarbChemAll$TA_CRM*(CarbChemAll$Salinity_Lab/34)

#Normalize DIC constant salinity
CarbChemAll$DIC_Norm<-CarbChemAll$DIC*(CarbChemAll$Salinity_Lab/34)

MicrobeSamplingchem<-CarbChemAll%>%
  filter(Time_Point == 1 | Time_Point == 4 | Time_Point == 5)
#write.csv(MicrobeSamplingchem,file="Output/MicrobeCarbChem.csv")
####NEC and NEP microbes####
#for TA, DIC, Sampling Time to find change over time
AllTPsamples<- CarbChemAll %>%
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  dplyr::filter(Id_code != 'D_24_2_BC'|Id_code != 'D_24_3_BC'| Id_code != 'D_24_4_BC'|Id_code == 'N_24_5_BC'|
                  Id_code != 'N_24_2_BC'|Id_code != 'N_24_3_BC'| Id_code != 'N_24_4_BC'|Id_code == 'D_24_5_BC'|
                  Id_code != 'D_1_2_AI'|Id_code != 'D_16_2_AI'|Id_code != 'D_4_3_AI'| 
                  Id_code != 'N_20_3_AI') %>% #remove samples that don't have nutrients
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::mutate(TADeltaTime = TA_NormSal - lag(TA_NormSal, default = first(TA_NormSal)), # subtracts time point 1&2,2&3,3&4
                DICDeltaTime = DIC_Norm-lag(DIC_Norm, default = first(DIC_Norm)),
                DeltaTime = (difftime(Date_Time,lag(Date_Time,default = first(Date_Time)))/3600), #sampling diff btwn time pts into hrs
                DeltaTime = ifelse(DeltaTime < 0 , DeltaTime + 24, DeltaTime), #if Delta Time is less than 0 then +24 hrs
                #converted to account from time went from 23/24 to 0:00 and added 24 hrs %>%
                DeltaNN = NN_umol_L - lag(NN_umol_L, default = first(NN_umol_L)),
                DeltaNH4 = NH4_umol_L - lag(NH4_umol_L, default = first (NH4_umol_L)),
                DeltaPO = PO_umol_L - lag(PO_umol_L, default =first(PO_umol_L)),
                DeltapH = pH_insitu -lag(pH_insitu, default = first (pH_insitu)),
                DeltapCO2 = pCO2 - lag(pCO2, default = first(pCO2)),
                DeltaDO = DO_mg_L - lag(DO_mg_L, default = first(DO_mg_L))) %>%
  dplyr::filter(TADeltaTime != 0) %>% #filters out row where deltas = 0
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,TADeltaTime, DICDeltaTime,
                DeltaTime,DeltaNN,DeltaNH4,DeltaPO,DeltapH,DeltapCO2,DeltaDO)

#because tp 24 in Before period did not have time pt 1 so took differences btw time point 2,3 and 4
TP24all <- CarbChemAll %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  dplyr::filter(Id_code == 'D_24_2_BC' | Id_code == 'D_24_3_BC'| Id_code == 'D_24_4_BC'|Id_code == 'D_24_5_BC'|
                  Id_code == 'N_24_2_BC'|Id_code == 'N_24_3_BC'|Id_code == 'N_24_4_BC'|Id_code == 'N_24_5_BC') %>%
  dplyr::filter(Time_Point == 2 | Time_Point == 3 | Time_Point == 4 | Time_Point == 5) %>%
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::mutate(TADeltaTime = TA_NormSal - lag(TA_NormSal, default = first(TA_NormSal)), # subtracts time point 2&3,3&4
                DICDeltaTime = DIC_Norm-lag(DIC_Norm, default = first(DIC_Norm)),
                DeltaTime = (difftime(Date_Time,lag(Date_Time,default = first(Date_Time)))/3600), #sampling diff btwn time pts into hrs
                DeltaTime = ifelse(DeltaTime < 0 , DeltaTime + 24, DeltaTime), #if Delta Time is less than 0 then +24 hrs
                #converted to account from time went from 23/24 to 0:00 and added 24 hrs %>%
                DeltaNN = NN_umol_L - lag(NN_umol_L, default = first(NN_umol_L)),
                DeltaNH4 = NH4_umol_L - lag(NH4_umol_L, default = first (NH4_umol_L)),
                DeltaPO = PO_umol_L - lag(PO_umol_L, default =first(PO_umol_L)),
                DeltapH = pH_insitu -lag(pH_insitu, default = first (pH_insitu)),
                DeltapCO2 = pCO2 - lag(pCO2, default = first(pCO2)),
                DeltaDO = DO_mg_L - lag(DO_mg_L, default = first(DO_mg_L))) %>%
  filter(TADeltaTime != 0) %>% #filters out row where deltas = 0
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,TADeltaTime, DICDeltaTime,
         DeltaTime,DeltaNN,DeltaNH4,DeltaPO,DeltapH,DeltapCO2,DeltaDO)

#combine dataframes
AllTPsamples<-rbind(AllTPsamples,TP24all)

#bring in physical parameters for nec and nep calculations
PhysicalParameters<-TidePooldes[, c("PoolID", "Before_After", "SurfaceArea", "Vol", "SAtoV", "TideHeight")] #pull out the necessary columns and treatment 
PhysicalParameters$PoolID<-as.character(PhysicalParameters$PoolID)
PhysicalParameters$Before_After<-as.character(PhysicalParameters$Before_After)

#join with delta samples
AllTPsamples<-left_join(AllTPsamples,PhysicalParameters)


#create an ifelse statement for sampling volumes
#Each time point more water was taken out of each pool and need to account for change in volume
#for NEC calculation
AllTPsamples <- AllTPsamples %>%
  mutate(SamplingVolume = ifelse(Time_Point == '2',Vol - 0.00055, #took out 550mL-->0.00055m3 
                                 #at time point 1 so volume at time pt 2 is -0.00055 m3
                                 ifelse(Time_Point == '3', Vol - 0.00095,
                                        #took out 400mL at time point 2 so now volume @time3
                                        #is - 950mL-->0.00095m3
                                        #took out 400 mL at time point 3 so volume @time4
                                        #is 550+400+400-->0.00135m3
                                        ifelse(Time_Point == '4', Vol - 0.00135,
                                               #took out 400 ml at time point 4
                                               #+0.0004=0.00175
                                               ifelse(Time_Point =='5', Vol - 0.00175,Vol))))) %>%
  
  #Some sampling times volume wasn't exact these are take from the note section in the data set 
  #note: these are for the time point after the note was taken
  mutate(AdjSamplingVolume = ifelse(Id_code == 'D_4_2_BC', SamplingVolume - 0.00055,
                                    ifelse(Id_code == 'D_18_2_BC', SamplingVolume - 0.00004,
                                           ifelse(Id_code == 'D_20_2_BC', SamplingVolume - 0.00008,
                                                  ifelse(Id_code == 'D_28_2_BC', SamplingVolume -0.0001,
                                                         ifelse(Id_code == 'D_5_3_BC', SamplingVolume -0.00005,
                                                                ifelse(Id_code == 'D_6_3_BC', SamplingVolume-0.00001,
                                                                       ifelse(Id_code == 'D_17_3_BC', SamplingVolume-0.00015,
                                                                              ifelse(Id_code == 'D_7_3_AI', SamplingVolume-0.00004,
                                                                                     ifelse(Id_code == 'D_9_3_AI', SamplingVolume-0.00004,
                                                                                            ifelse(Id_code == 'N_2_3_AI', SamplingVolume-0.0005,
                                                                                                   ifelse(Id_code == 'N_9_2_AI', SamplingVolume + 0.00001,
                                                                                                          ifelse(Id_code == 'N_1_3_AI', SamplingVolume - 0.00002,
                                                                                                                 ifelse(Id_code == 'N_2_3_AI', SamplingVolume -0.000018,
                                                                                                                        ifelse(Id_code == 'N_4_3_AI', SamplingVolume -0.000015,
                                                                                                                               ifelse(Id_code == 'N_5_3_AI',SamplingVolume -0.00001,
                                                                                                                                      ifelse(Id_code == 'N_6_3_AI',SamplingVolume -0.00004,
                                                                                                                                             ifelse(Id_code == 'N_12_3_AI', SamplingVolume -0.000035,
                                                                                                                                                    ifelse(Id_code == 'D_4_4_BC', SamplingVolume - 0.00055,
                                                                                                                                                           ifelse(Id_code == 'D_18_4_BC', SamplingVolume - 0.00004,
                                                                                                                                                                  ifelse(Id_code == 'D_20_4_BC', SamplingVolume - 0.00008,
                                                                                                                                                                         ifelse(Id_code == 'D_28_4_BC', SamplingVolume -0.0001,
                                                                                                                                                                                ifelse(Id_code == 'D_5_4_BC', SamplingVolume -0.00005,
                                                                                                                                                                                       ifelse(Id_code == 'D_6_4_BC', SamplingVolume-0.00001,
                                                                                                                                                                                              ifelse(Id_code == 'D_17_4_BC', SamplingVolume-0.00015,
                                                                                                                                                                                                     ifelse(Id_code == 'D_7_4_AI', SamplingVolume-0.00004,
                                                                                                                                                                                                            ifelse(Id_code == 'D_9_4_AI', SamplingVolume-0.00004,
                                                                                                                                                                                                                   ifelse(Id_code == 'N_2_4_AI', SamplingVolume-0.0005,
                                                                                                                                                                                                                          ifelse(Id_code == 'N_9_4_AI', SamplingVolume + 0.00001,
                                                                                                                                                                                                                                 ifelse(Id_code == 'N_1_4_AI', SamplingVolume - 0.00002,
                                                                                                                                                                                                                                        ifelse(Id_code == 'N_2_4_AI', SamplingVolume -0.000018,
                                                                                                                                                                                                                                               ifelse(Id_code == 'N_4_4_AI', SamplingVolume -0.000015,
                                                                                                                                                                                                                                                      ifelse(Id_code == 'N_5_4_AI',SamplingVolume -0.00001,
                                                                                                                                                                                                                                                             ifelse(Id_code == 'N_6_4_AI',SamplingVolume -0.00004,
                                                                                                                                                                                                                                                                    ifelse(Id_code == 'N_12_4_AI', SamplingVolume -0.000035,
                                                                                                                                                                                                                                                                           ifelse(Id_code == 'N_5_4_AI', SamplingVolume -0.00002,
                                                                                                                                                                                                                                                                                  ifelse(Id_code == 'N_9_4_AI', SamplingVolume -0.00001, SamplingVolume)))))))))))))))))))))))))))))))))))))



#Normalizing change in TA to change in nutrients based on Wolf-Gladrow et al. 2007
AllTPsamples$DeltaTA_N_Norm<-AllTPsamples$TADeltaTime - (AllTPsamples$DeltaNN) - (2*AllTPsamples$DeltaPO) + (AllTPsamples$DeltaNH4)
#for every mol Nitrate and phosophate--increases TA by one and two moles respectively
#so subtract to normalise TA where as TA decreases with every mole of Nh4 so + Nh4 back

AllTPsamples$DeltaTA_N_Norm<- -1 * AllTPsamples$DeltaTA_N_Norm #to change to positive calcification and negative dissolution
AllTPsamples$DICDeltaTime <- -1 * (AllTPsamples$DICDeltaTime) #to change to positive photosynthesis and neg respiration 


AllTPsamples$DeltaTime<-as.numeric(AllTPsamples$DeltaTime) #as numeric since attributes from difftime function were interferring 
#with NEC/NEP calcs

#Equation for NEC rates
#NEC delta TA/2 * seawater density * (volume/ surface area) / time  Divided by 1000 so mmols m2 hr 
AllTPsamples$NEC.mmol.m2.hr<- ((AllTPsamples$DeltaTA_N_Norm)/2) * (1023) * (AllTPsamples$AdjSamplingVolume/AllTPsamples$SurfaceArea) * (1/AllTPsamples$DeltaTime) * (1/1000)
#without nutrients:
#DeltaSamples$NECnotNNorm<- ((DeltaSamples$TADeltaTime)/2) * (1023) * (DeltaSamples$AdjSamplingVolume/DeltaSamples$SurfaceArea) * (1/DeltaSamples$DeltaTime) * (1/1000)
#Function for calculation of air-sea CO2 flux
#adapted from MatLab Code by Cecilia Chapa Balcorta copyright @ 2015

#input:
#pCO2_agua= seawater pCO2 (uatm)
#pCO2_atm=  atmospheric pCO2 (uatm)
#T=  Temperature (Celsius)
#S=  Salinity 
#u = Wind speed (m/s)

#get average salinity, temp and pco2 between time points
AirSeaFluxTPall<- CarbChemAll %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::filter(Id_code != 'D_24_2_BC'|Id_code != 'D_24_3_BC'| Id_code != 'D_24_4_BC'|
                  Id_code != 'N_24_2_BC'|Id_code != 'N_24_3_BC'| Id_code != 'N_24_4_BC'|Id_code == 'D_24_5_BC'|
                  Id_code != 'D_1_2_AI' |Id_code != 'D_16_2_AI'|Id_code != 'D_4_3_AI'|Id_code == 'N_24_5_BC'|
                  Id_code != 'N_20_3_AI') %>% #take out these ids since not following 4-1 rule
  #dplyr::filter(Time_Point == 1 | Time_Point == 4) %>%
  dplyr::mutate(Tempmean = rollmean(Temp.pool,2, na.pad=TRUE, align="right"), #takes rolling means every two pts
                pCO2mean = rollmean(pCO2,2, na.pad=TRUE, align="right"), #avg btw time 1 & 2, time 2&3 etc. 
                Salinitymean = rollmean(Salinity,2, na.pad=TRUE, align="right"),
                Windmean = rollmean(Wind_Speed.m.s,2, na.pad=TRUE, align="right"),
                pHmean = rollmean(pH_insitu, 2, na.pad=TRUE, align="right"),
                DOmean = rollmean(DO_mg_L, 2, na.pad=TRUE, align="right"),
                NNmean = rollmean(NN_umol_L, 2, na.pad=TRUE, align="right"),
                NH4mean = rollmean(NH4_umol_L, 2, na.pad=TRUE, align="right"),
                POmean = rollmean(PO_umol_L, 2, na.pad=TRUE, align="right")) %>% 
  filter(Tempmean != 'NA') %>% #removes row that values are NA (time 1:time1) 
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,
         Tempmean,pCO2mean,Salinitymean,Windmean,pHmean,DOmean, NNmean, NH4mean,POmean) #selects the columns needed for air-sea flux equation

#TP24
TP24AirseafluxTPall<- CarbChemAll %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::filter(Id_code == 'D_24_2_BC' | Id_code == 'D_24_3_BC'| Id_code == 'D_24_4_BC'|Id_code == 'D_24_5_BC'|
                  Id_code == 'N_24_2_BC'|Id_code == 'N_24_3_BC'|Id_code == 'N_24_4_BC'|Id_code == 'N_24_5_BC') %>%
  dplyr::filter(Time_Point == 2 |Time_Point == 3| Time_Point == 4 | Time_Point ==5) %>% #no time point 1
  dplyr::mutate(Tempmean = rollmean(Temp.pool,2, na.pad=TRUE, align="right"), #takes rolling means every two pts
                pCO2mean = rollmean(pCO2,2, na.pad=TRUE, align="right"), #avg btw time 2& 3, time 3&4 etc. 
                Salinitymean = rollmean(Salinity,2, na.pad=TRUE, align="right"),
                Windmean = rollmean(Wind_Speed.m.s,2, na.pad=TRUE, align="right"),
                pHmean = rollmean(pH_insitu, 2, na.pad=TRUE, align="right"),
                DOmean = rollmean(DO_mg_L, 2, na.pad=TRUE, align="right"),
                NNmean = rollmean(NN_umol_L, 2, na.pad=TRUE, align="right"),
                NH4mean = rollmean(NH4_umol_L, 2, na.pad=TRUE, align="right"),
                POmean = rollmean(PO_umol_L, 2, na.pad=TRUE, align="right")) %>% 
  filter(Tempmean != 'NA') %>% #removes row that values are NA (time 1:time1) 
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,
         Tempmean,pCO2mean,Salinitymean,Windmean,pHmean,DOmean,NNmean, NH4mean,POmean) #selects the columns needed for air-sea flux equation


#combine dataframes
AirSeaFluxall<-rbind(AirSeaFluxTPall,TP24AirseafluxTPall)

#combine with delta samples dataframe
AllTPsamples<-left_join(AllTPsamples,AirSeaFluxall)

#Setting parameters to make it easier in functions
T<-AllTPsamples$Tempmean
S<-AllTPsamples$Salinitymean
u<-AllTPsamples$Windmean

pCO2_water<-AllTPsamples$pCO2mean
pCO2_atm<- 410.86 #average uatm CO2 in atms btween July (411.77) and August 2019 (409.95) from Mauna Loa Station


#Air-sea CO2 is calculated as follows:

#FCO2 = K * a * (dpCO2) 

#Where K=is the transfer velocity according to Wanninkhof (1992).
#a = CO2 solibility constant according to Weiss (1974)
#dpCO2 is the difference of air and seawater pCO2 


#Schmidt Number function
#For water of salinity=35 and temperature range 0-30deg C
# A,B,C,D are constants
Sc<-function(T) {
  A <- 2073.1
  B <- 125.62
  C <- 3.6276
  D <- 0.043219
  Sc <- A - (B * T) + (C *(T^2)) - (D * (T^3))
  return(Sc)
}

#Solibuility constant (Weiss, 1974) 
Ko_weiss<- function(T,S) {
  A <- c(-60.2409, 93.4517, 23.3585) #mol/kg.atm
  B <- c(0.023517, -0.023656, 0.0047036)  #mol/kg.atm
  T <- T + 273.15 #Conversion from Celsius degrees to Kelvins
  Ln_Ko = A[1] + (A[2] * (100/T)) + (A[3] * (log(T/100))) + S * (B[1] + (B[2] * (T/100)) + (B[3] * (T/100)^2))
  Ko <- exp(Ln_Ko)
  return(Ko)
}

#CO2 Transfer velocity calculation 
slowwind<-0.31*(u^2) * ((Sc(T)/660)^-0.5) #for wind <=6 m/s
highwind<-0.39*(u^2) * ((Sc(T)/660)^-0.5) #for wind > 6 m/s

AllTPsamples$K<- ifelse(u <= 6, slowwind, highwind) #if else statement for K [transfer velocity according to Wanninkhof (1992)]

AllTPsamples$dpCO2 <- pCO2_water - pCO2_atm #difference of air and seawater pCO2

AllTPsamples$a <- Ko_weiss(T,S) #CO2 solibility constant according to Weiss (1974) solubility in mmol L^-1 atm^-1 or mmol m^-3 uatm^-1

AllTPsamples$F_CO2 <- ((0.24 * AllTPsamples$K * AllTPsamples$a * AllTPsamples$dpCO2)/ 24) #CO2 flux (mmol m^-2 hr^-1) (rate is by day so divide by 24 to convert per hour)


#NEP with FCO2
# NEP = (delta DIC) * SW density * (volume/SA) / time  -  NEC - Fugosity of CO2) *1000 to convert to mmol 
AllTPsamples$NEP.mmol.m2.hr<-((AllTPsamples$DICDeltaTime) * (1023) * (AllTPsamples$AdjSamplingVolume/AllTPsamples$SurfaceArea) * (1/AllTPsamples$DeltaTime) * (1/1000)) - (AllTPsamples$NEC.mmol.m2.hr) - (AllTPsamples$F_CO2)

AllTPsamples$NEC.mmol.m2.hr<-as.numeric(AllTPsamples$NEC.mmol.m2.hr) 
AllTPsamples$NEP.mmol.m2.hr<-as.numeric(AllTPsamples$NEP.mmol.m2.hr)

Microbesnecnep<-AllTPsamples%>%
  filter(Time_Point ==4|Time_Point ==5)
#write.csv(Microbesnecnep, file="Output/MicrobesTime4and5NECNEP.csv")

####Integrated measurements over low tide 1-4 or 1-5#####

#Dataframe of microbe samples were taken at time 4
Time4CarbChem<-CarbChemAll%>%
  filter(Before_After =="After" | Sampling_Day == "2019-07-12")
#Dataframe when microbe samples were taken at time 5
Time5CarbChem<-CarbChemAll%>%
  filter(Before_After =="Before" & Sampling_Day != "2019-07-12")

#for TA, DIC, Sampling Time to find change over time
Time5Delta<- Time5CarbChem%>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  dplyr::filter(Id_code != 'D_24_2_BC',Id_code != 'D_24_3_BC', Id_code != 'D_24_4_BC',Id_code != 'D_24_5_BC',
                Id_code != 'N_24_2_BC',Id_code != 'N_24_3_BC', Id_code != 'N_24_4_BC',Id_code != 'N_24_5_BC') %>%
  dplyr::filter(Time_Point == 1 |Time_Point ==5) %>% #taking the integrated value over low tide
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::mutate(TADeltaTime = TA_NormSal - lag(TA_NormSal, default = first(TA_NormSal)), # subtracts time point 4 from 1
                DICDeltaTime = DIC_Norm-lag(DIC_Norm, default = first(DIC_Norm)),
                DeltaTime = (difftime(Date_Time,lag(Date_Time,default = first(Date_Time)))/3600), #sampling diff btwn time pts into hrs
                DeltaTime = ifelse(DeltaTime < 0 , DeltaTime + 24, DeltaTime), #if Delta Time is less than 0 then +24 hrs
                #converted to account from time went from 23/24 to 0:00 and added 24 hrs %>%
                DeltaNN = NN_umol_L - lag(NN_umol_L, default = first(NN_umol_L)),
                DeltaNH4 = NH4_umol_L - lag(NH4_umol_L, default = first (NH4_umol_L)),
                DeltaPO = PO_umol_L - lag(PO_umol_L, default =first(PO_umol_L)),
                DeltapH = pH_insitu -lag(pH_insitu, default = first (pH_insitu)),
                DeltapCO2 = pCO2 - lag(pCO2, default = first(pCO2)),
                DeltaDO = DO_mg_L - lag(DO_mg_L, default = first(DO_mg_L))) %>%
  dplyr::filter(TADeltaTime != 0) %>% #filters out row where deltas = 0
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,TADeltaTime, DICDeltaTime,
                DeltaTime,DeltaNN,DeltaNH4,DeltaPO,DeltapH,DeltapCO2,DeltaDO)

#because tp 24 in Before period did not have time pt 1 so took differences btw time point 2 and 4
TP24Time5 <- Time5CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  dplyr::filter(Id_code == 'D_24_2_BC' | Id_code == 'D_24_3_BC'| Id_code == 'D_24_4_BC',Id_code != 'D_24_5_BC',
                  Id_code == 'N_24_2_BC'|Id_code == 'N_24_3_BC'|Id_code == 'N_24_4_BC',Id_code != 'N_24_5_BC') %>%
  dplyr::filter(Time_Point == 2 |Time_Point == 5) %>%
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::mutate(TADeltaTime = TA_NormSal - lag(TA_NormSal, default = first(TA_NormSal)), # subtracts time point 1&2,2&3,3&4
                DICDeltaTime = DIC_Norm-lag(DIC_Norm, default = first(DIC_Norm)),
                DeltaTime = (difftime(Date_Time,lag(Date_Time,default = first(Date_Time)))/3600), #sampling diff btwn time pts into hrs
                DeltaTime = ifelse(DeltaTime < 0 , DeltaTime + 24, DeltaTime), #if Delta Time is less than 0 then +24 hrs
                #converted to account from time went from 23/24 to 0:00 and added 24 hrs %>%
                DeltaNN = NN_umol_L - lag(NN_umol_L, default = first(NN_umol_L)),
                DeltaNH4 = NH4_umol_L - lag(NH4_umol_L, default = first (NH4_umol_L)),
                DeltaPO = PO_umol_L - lag(PO_umol_L, default =first(PO_umol_L)),
                DeltapH = pH_insitu -lag(pH_insitu, default = first (pH_insitu)),
                DeltapCO2 = pCO2 - lag(pCO2, default = first(pCO2)),
                DeltaDO = DO_mg_L - lag(DO_mg_L, default = first(DO_mg_L))) %>%
  filter(TADeltaTime != 0) %>% #filters out row where deltas = 0
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,TADeltaTime, DICDeltaTime,
         DeltaTime,DeltaNN,DeltaNH4,DeltaPO,DeltapH,DeltapCO2,DeltaDO)

#combine dataframes
Time5Delta<-rbind(Time5Delta,TP24Time5)

#for TA, DIC, Sampling Time to find change over time 
Time4Delta<- Time4CarbChem%>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  dplyr::filter(Time_Point == 1 |Time_Point ==4) %>% #taking the integrated value over low tide
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::mutate(TADeltaTime = TA_NormSal - lag(TA_NormSal, default = first(TA_NormSal)), # subtracts time point 4 from 1
                DICDeltaTime = DIC_Norm-lag(DIC_Norm, default = first(DIC_Norm)),
                DeltaTime = (difftime(Date_Time,lag(Date_Time,default = first(Date_Time)))/3600), #sampling diff btwn time pts into hrs
                DeltaTime = ifelse(DeltaTime < 0 , DeltaTime + 24, DeltaTime), #if Delta Time is less than 0 then +24 hrs
                #converted to account from time went from 23/24 to 0:00 and added 24 hrs %>%
                DeltaNN = NN_umol_L - lag(NN_umol_L, default = first(NN_umol_L)),
                DeltaNH4 = NH4_umol_L - lag(NH4_umol_L, default = first (NH4_umol_L)),
                DeltaPO = PO_umol_L - lag(PO_umol_L, default =first(PO_umol_L)),
                DeltapH = pH_insitu -lag(pH_insitu, default = first (pH_insitu)),
                DeltapCO2 = pCO2 - lag(pCO2, default = first(pCO2)),
                DeltaDO = DO_mg_L - lag(DO_mg_L, default = first(DO_mg_L))) %>%
  dplyr::filter(TADeltaTime != 0) %>% #filters out row where deltas = 0
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,TADeltaTime, DICDeltaTime,
                DeltaTime,DeltaNN,DeltaNH4,DeltaPO,DeltapH,DeltapCO2,DeltaDO)

AllTimePts<-rbind(Time5Delta,Time4Delta)
#bring in physical parameters for nec and nep calculations
PhysicalParameters<-TidePooldes[, c("PoolID", "Before_After", "SurfaceArea", "Vol", "SAtoV", "TideHeight")] #pull out the necessary columns and treatment 
PhysicalParameters$PoolID<-as.character(PhysicalParameters$PoolID)
PhysicalParameters$Before_After<-as.character(PhysicalParameters$Before_After)

#join with delta samples
AllTimePts<-left_join(AllTimePts,PhysicalParameters)

#create an ifelse statement for sampling volumes
#Each time point more water was taken out of each pool and need to account for change in volume
#for NEC calculation
AllTimePts <- AllTimePts %>%
  mutate(SamplingVolume = #ifelse(Time_Point == '2',Vol - 0.00055, #took out 550mL-->0.00055m3 
           #at time point 1 so volume at time pt 2 is -0.00055 m3
           #ifelse(Time_Point == '3', Vol - 0.00095,
           #took out 400mL at time point 2 so now volume @time3
           #is - 950mL-->0.00095m3
           ifelse(Time_Point == '4', Vol - 0.00135, 
                  ifelse(Time_Point == '5', Vol -0.00175,Vol)))%>%
  #took out 400 mL at time point 3 so volume @time4
  #is 550+400+400-->0.00135m3; took 400ml out at time 4 so time 5 is 550+400+400+400-->0.00175m3
  #Some sampling times volume wasn't exact these are take from the note section in the data set 
  #note: Time 5 for all before and time 4 for afters
  mutate(AdjSamplingVolume = ifelse(Id_code == 'D_4_5_BC', SamplingVolume - 0.00055,
ifelse(Id_code == 'D_18_5_BC', SamplingVolume - 0.00004,
ifelse(Id_code == 'D_20_5_BC', SamplingVolume - 0.00008,
ifelse(Id_code == 'D_28_5_BC', SamplingVolume -0.0001,
ifelse(Id_code == 'D_5_5_BC', SamplingVolume -0.00005,
ifelse(Id_code == 'D_6_5_BC', SamplingVolume-0.00001,
ifelse(Id_code == 'D_17_5_BC', SamplingVolume-0.00015,
ifelse(Id_code == 'D_7_4_AI', SamplingVolume-0.00004,
ifelse(Id_code == 'D_9_4_AI', SamplingVolume-0.00004,
ifelse(Id_code == 'N_2_4_AI', SamplingVolume-0.0005,
ifelse(Id_code == 'N_9_4_AI', SamplingVolume + 0.00001,
ifelse(Id_code == 'N_1_4_AI', SamplingVolume - 0.00002,
ifelse(Id_code == 'N_2_4_AI', SamplingVolume -0.000018,
ifelse(Id_code == 'N_4_4_AI', SamplingVolume -0.000015,
ifelse(Id_code == 'N_5_4_AI',SamplingVolume -0.00001,
ifelse(Id_code == 'N_6_4_AI',SamplingVolume -0.00004,
ifelse(Id_code == 'N_12_4_AI', SamplingVolume -0.000035,
ifelse(Id_code == 'N_5_4_AI', SamplingVolume -0.00002,
ifelse(Id_code == 'N_9_4_AI', SamplingVolume -0.00001, SamplingVolume)))))))))))))))))))) 



#Normalizing change in TA to change in nutrients based on Wolf-Gladrow et al. 2007
AllTimePts$DeltaTA_N_Norm<-AllTimePts$TADeltaTime - (AllTimePts$DeltaNN) - (2*AllTimePts$DeltaPO) + (AllTimePts$DeltaNH4)
#for every mol Nitrate and phosophate--increases TA by one and two moles respectively
#so subtract to normalise TA where as TA decreases with every mole of Nh4 so + Nh4 back

AllTimePts$DeltaTA_N_Norm<- -1 * AllTimePts$DeltaTA_N_Norm #to change to positive calcification and negative dissolution
AllTimePts$DICDeltaTime <- -1 * (AllTimePts$DICDeltaTime) #to change to positive photosynthesis and neg respiration 


AllTimePts$DeltaTime<-as.numeric(AllTimePts$DeltaTime) #as numeric since attributes from difftime function were interferring 
#with NEC/NEP calcs

#Equation for NEC rates
#NEC delta TA/2 * seawater density * (volume/ surface area) / time  Divided by 1000 so mmols m2 hr 
AllTimePts$NEC.mmol.m2.hr<- ((AllTimePts$DeltaTA_N_Norm)/2) * (1023) * (AllTimePts$AdjSamplingVolume/AllTimePts$SurfaceArea) * (1/AllTimePts$DeltaTime) * (1/1000)
#without nutrients:
#DeltaSamples$NECnotNNorm<- ((DeltaSamples$TADeltaTime)/2) * (1023) * (DeltaSamples$AdjSamplingVolume/DeltaSamples$SurfaceArea) * (1/DeltaSamples$DeltaTime) * (1/1000)
#Function for calculation of air-sea CO2 flux
#adapted from MatLab Code by Cecilia Chapa Balcorta copyright @ 2015

#input:
#pCO2_agua= seawater pCO2 (uatm)
#pCO2_atm=  atmospheric pCO2 (uatm)
#T=  Temperature (Celsius)
#S=  Salinity 
#u = Wind speed (m/s)

#get average salinity, temp and pco2 between time points
AirSeaFluxTime5<- Time5CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::filter(Id_code != 'D_24_2_BC',Id_code != 'D_24_3_BC', Id_code != 'D_24_4_BC',Id_code != 'D_24_5_BC',
                Id_code != 'N_24_2_BC',Id_code != 'N_24_3_BC', Id_code != 'N_24_4_BC',Id_code != 'N_24_5_BC') %>% #take these out because not follow 1-5 rule
  dplyr::filter(Time_Point == 1 | Time_Point == 5) %>%
  dplyr::mutate(Tempmean = rollmean(Temp.pool,2, na.pad=TRUE, align="right"), #takes rolling means every two pts
                pCO2mean = rollmean(pCO2,2, na.pad=TRUE, align="right"), #avg btw time 1 & 2, time 2&3 etc. 
                Salinitymean = rollmean(Salinity,2, na.pad=TRUE, align="right"),
                Windmean = rollmean(Wind_Speed.m.s,2, na.pad=TRUE, align="right"),
                pHmean = rollmean(pH_insitu, 2, na.pad=TRUE, align="right"),
                DOmean = rollmean(DO_mg_L, 2, na.pad=TRUE, align="right"),
                NNmean = rollmean(NN_umol_L, 2, na.pad=TRUE, align="right"),
                NH4mean = rollmean(NH4_umol_L, 2, na.pad=TRUE, align="right"),
                POmean = rollmean(PO_umol_L, 2, na.pad=TRUE, align="right")) %>% 
  filter(Tempmean != 'NA') %>% #removes row that values are NA (time 1:time1) 
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,
         Tempmean,pCO2mean,Salinitymean,Windmean,pHmean,DOmean, NNmean, NH4mean,POmean) #selects the columns needed for air-sea flux equation

#TP24 time 5
TP24time5Airflux<- Time5CarbChem%>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::filter(Id_code == 'D_24_2_BC' | Id_code == 'D_24_3_BC'| Id_code == 'D_24_4_BC'|Id_code == 'D_24_5_BC'|
                  Id_code == 'N_24_2_BC'|Id_code == 'N_24_3_BC'|Id_code == 'N_24_4_BC'|Id_code == 'N_24_5_BC') %>%
  dplyr::filter(Time_Point == 2 | Time_Point == 5) %>%
  dplyr::mutate(Tempmean = rollmean(Temp.pool,2, na.pad=TRUE, align="right"), #takes rolling means every two pts
                pCO2mean = rollmean(pCO2,2, na.pad=TRUE, align="right"), #avg btw time 1 & 2, time 2&3 etc. 
                Salinitymean = rollmean(Salinity,2, na.pad=TRUE, align="right"),
                Windmean = rollmean(Wind_Speed.m.s,2, na.pad=TRUE, align="right"),
                pHmean = rollmean(pH_insitu, 2, na.pad=TRUE, align="right"),
                DOmean = rollmean(DO_mg_L, 2, na.pad=TRUE, align="right"),
                NNmean = rollmean(NN_umol_L, 2, na.pad=TRUE, align="right"),
                NH4mean = rollmean(NH4_umol_L, 2, na.pad=TRUE, align="right"),
                POmean = rollmean(PO_umol_L, 2, na.pad=TRUE, align="right")) %>% 
  filter(Tempmean != 'NA') %>% #removes row that values are NA (time 1:time1) 
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,
         Tempmean,pCO2mean,Salinitymean,Windmean,pHmean,DOmean,NNmean, NH4mean,POmean) #selects the columns needed for air-sea flux equation

#combine dataframes
AirSeaFluxTime5<-rbind(AirSeaFluxTime5,TP24time5Airflux)

AirSeaFluxTime4<- Time4CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::filter(Time_Point == 1 | Time_Point == 4) %>%
  dplyr::mutate(Tempmean = rollmean(Temp.pool,2, na.pad=TRUE, align="right"), #takes rolling means every two pts
                pCO2mean = rollmean(pCO2,2, na.pad=TRUE, align="right"), #avg btw time 1 & 2, time 2&3 etc. 
                Salinitymean = rollmean(Salinity,2, na.pad=TRUE, align="right"),
                Windmean = rollmean(Wind_Speed.m.s,2, na.pad=TRUE, align="right"),
                pHmean = rollmean(pH_insitu, 2, na.pad=TRUE, align="right"),
                DOmean = rollmean(DO_mg_L, 2, na.pad=TRUE, align="right"),
                NNmean = rollmean(NN_umol_L, 2, na.pad=TRUE, align="right"),
                NH4mean = rollmean(NH4_umol_L, 2, na.pad=TRUE, align="right"),
                POmean = rollmean(PO_umol_L, 2, na.pad=TRUE, align="right")) %>% 
  filter(Tempmean != 'NA') %>% #removes row that values are NA (time 1:time1) 
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,
         Tempmean,pCO2mean,Salinitymean,Windmean,pHmean,DOmean, NNmean, NH4mean,POmean) #selects the columns needed for air-sea flux equation

AllAirseaflux<-rbind(AirSeaFluxTime5,AirSeaFluxTime4)
#combine with delta samples dataframe
AllTimesamples<-left_join(AllTimePts,AllAirseaflux)

#Setting parameters to make it easier in functions
T<-AllTimesamples$Tempmean
S<-AllTimesamples$Salinitymean
u<-AllTimesamples$Windmean

pCO2_water<-AllTimesamples$pCO2mean
pCO2_atm<- 410.86 #average uatm CO2 in atms btween July (411.77) and August 2019 (409.95) from Mauna Loa Station


#Air-sea CO2 is calculated as follows:

#FCO2 = K * a * (dpCO2) 

#Where K=is the transfer velocity according to Wanninkhof (1992).
#a = CO2 solibility constant according to Weiss (1974)
#dpCO2 is the difference of air and seawater pCO2 


#Schmidt Number function
#For water of salinity=35 and temperature range 0-30deg C
# A,B,C,D are constants
Sc<-function(T) {
  A <- 2073.1
  B <- 125.62
  C <- 3.6276
  D <- 0.043219
  Sc <- A - (B * T) + (C *(T^2)) - (D * (T^3))
  return(Sc)
}

#Solibuility constant (Weiss, 1974) 
Ko_weiss<- function(T,S) {
  A <- c(-60.2409, 93.4517, 23.3585) #mol/kg.atm
  B <- c(0.023517, -0.023656, 0.0047036)  #mol/kg.atm
  T <- T + 273.15 #Conversion from Celsius degrees to Kelvins
  Ln_Ko = A[1] + (A[2] * (100/T)) + (A[3] * (log(T/100))) + S * (B[1] + (B[2] * (T/100)) + (B[3] * (T/100)^2))
  Ko <- exp(Ln_Ko)
  return(Ko)
}

#CO2 Transfer velocity calculation 
slowwind<-0.31*(u^2) * ((Sc(T)/660)^-0.5) #for wind <=6 m/s
highwind<-0.39*(u^2) * ((Sc(T)/660)^-0.5) #for wind > 6 m/s

AllTimesamples$K<- ifelse(u <= 6, slowwind, highwind) #if else statement for K [transfer velocity according to Wanninkhof (1992)]

AllTimesamples$dpCO2 <- pCO2_water - pCO2_atm #difference of air and seawater pCO2

AllTimesamples$a <- Ko_weiss(T,S) #CO2 solibility constant according to Weiss (1974) solubility in mmol L^-1 atm^-1 or mmol m^-3 uatm^-1

AllTimesamples$F_CO2 <- ((0.24 * AllTimesamples$K * AllTimesamples$a * AllTimesamples$dpCO2)/ 24) #CO2 flux (mmol m^-2 hr^-1) (rate is by day so divide by 24 to convert per hour)


#NEP with FCO2
# NEP = (delta DIC) * SW density * (volume/SA) / time  -  NEC - Fugosity of CO2) *1000 to convert to mmol 
AllTimesamples$NEP.mmol.m2.hr<-((AllTimesamples$DICDeltaTime) * (1023) * (AllTimesamples$AdjSamplingVolume/AllTimesamples$SurfaceArea) * (1/AllTimesamples$DeltaTime) * (1/1000)) - (AllTimesamples$NEC.mmol.m2.hr) - (AllTimesamples$F_CO2)

AllTimesamples$NEC.mmol.m2.hr<-as.numeric(AllTimesamples$NEC.mmol.m2.hr) 
AllTimesamples$NEP.mmol.m2.hr<-as.numeric(AllTimesamples$NEP.mmol.m2.hr)

#write.csv(AllTimesamples, "Output/Integratedtime1thru4or5necnep.csv")
