# TidePoolMicrobes
Microbial communities of mussel, surfgrass, and removal Oregon tide pools at Otter Rock Marine Reserve

Contents

**[Protocols](https://github.com/jenniferfields/TidePoolMicrobes/tree/master/Protocols)**

Contains all protocols for [biogeochemistry and microbial sampling](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Protocols/TidePoolSampling_SOP.md), [community composition](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Protocols/CommunityComposition_SOP.md), and [tide pool volume](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Protocols/Dye_Method_Protocol.md)

**[Scripts](https://github.com/jenniferfields/TidePoolMicrobes/tree/master/Scripts)**
* [Biogeochemistry Script](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Scripts/MicrobeCommBiogeochem.R) - script for processing all carbonate and biogeochemistry data and calculation of net ecosystem production and calcification
* [Community Composition Script](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Scripts/TidePoolCommunityComposition.R) - script for processing community composition data for sessile and mobile organisms
* [Tide Pool Physical Parameters Script](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Scripts/tidepoolphysicalparameters.R) - script for processing physical parameters of the tide pools  
* [Combined microbe dataset](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Scripts/combinedmicrobedatasets.R) - script to combine fDOM and FCM datasets

**[Data](https://github.com/jenniferfields/TidePoolMicrobes/tree/master/Data)**
Contains all clean and raw data from project.
1. *[Biogeochemistry Data](https://github.com/jenniferfields/TidePoolMicrobes/tree/master/Data/Biogeochem)*

	* [Integrated NEP and NEC](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/Biogeochem/Integratedtime1thru4or5necnep.csv) - calculated net ecosystem production and net ecosystem calcification over the low tide period (Time 1-4/5)
	*  [NEP and NEC for last time point](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/Biogeochem/MicrobesTime4and5NECNEP.csv) - NEP and NEC for only last time point of the series (either time point 4 or time point 5)
	* [Microbe Biogeochm](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/Biogeochem/MicrobeCarbChem.csv) - Contains all carbonate and biogeochemistry over the microbe sampling time points
		* Columns to use: Dissolved oxygen (DO) (mg/L & %), Conductivity, Salinity, PO_umol_L, NN_umol_L, NH4_umol_L, pH_insitu (this is the corrected pH to use), CO2, fCO2, pCO2, HCO3, CO3, DIC, ALK, OmegaAragonite, OmegaCalcite, TA_NormSal & DIC_Norm (normalized to salinity in lab)
	*  [NEP and NEC for last time point](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/Biogeochem/MicrobesTime4and5NECNEP.csv) - NEP and NEC for only last time point of the series (either time point 4 or time point 5)
	*  Other data files are needed to run script
  
3. *[Community Composition Data](https://github.com/jenniferfields/TidePoolMicrobes/tree/master/Data/CommunityComposition)*

	*	[Sessile species list](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/CommunityComposition/Sessilesspplist.csv) - all sessile community with data normalized to the sum of the total cover (since it can be greater than 100%)
	*	[Sessile functional groups](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/CommunityComposition/Sessilefunctionalgroups.csv) - sessile data switched to long format with functional group attached to each species on the list
 
	*	[Mobile species list](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/CommunityComposition/Mobilespplist.csv) - all mobile species 

	*	[Mobile functional groups](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/CommunityComposition/Mobilespplist.csv) - mobile data switched to long format with functional group attached to each species on the list
	* Others are the raw datasets
4. *[Light and Temperature Logger Data](https://github.com/jenniferfields/TidePoolMicrobes/tree/master/Data/LightandTemp)*
	* [Temp and light by hour](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/LightandTemp/LightandTempbysamplinghour.csv) - Temperature and light values from HOBO® Pendant loggers averaged for each hour of sampling (n = 3-4 values averaged)
	* Others are for script to run
5. *[Microbe Data](https://github.com/jenniferfields/TidePoolMicrobes/tree/master/Data/Microbe_Clean)*
	* [fDOM and FCM data](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/Microbe_Clean/combined_FCMandfDOMdata.csv) - Combined fDOM and FCM datasets
	* Others are fDOM and FCM as separate datasets 
6. *[Pool Physical Parameters Data](https://github.com/jenniferfields/TidePoolMicrobes/tree/master/Data/PoolPhysicalParameters)*
	*	[Physical parameter measurements](https://github.com/jenniferfields/TidePoolMicrobes/blob/master/Data/PoolPhysicalParameters/TidePoolDes.csv)
	*	Other data files are needed to run script

7. [Microbe datasets we edited for merging with rest of data](https://github.com/jenniferfields/TidePoolMicrobes/tree/master/Data/Microbes_withnotes)

