### Data visualization and analysis for fDOM and microbe datasets ###
### By Nyssa Silbiger #####
### Created on 8/16/2023 #####


### Load libraries   

library(here)
library(tidyverse)
library(janitor)
library(ggtext)
library(patchwork)
library(calecopal)
library(lubridate)
library(lme4)
library(lmerTest)
library(vegan)
library(broom)
library(broom.mixed)

## read in the data ####

microbedata<-read_csv(here("Data","Microbe_Clean","combined_FCMandfDOMdata.csv"))

biogeodata<-read_csv(here("Data","Biogeochem","MicrobeCarbChem.csv"))

BenthicData<-read_csv(here("Data","Microbe_Clean","CommunityData.csv"))

MetaData<-read_csv(here("Data","Microbe_Clean","TidePoolDescriptions.csv"))

# Pull out the NEC and NEP rates
MetabRates<-read_csv(here("Output","MicrobesTime4and5NECNEP.csv")) %>%
  filter(Time_Point == 4) %>% # use th3 4 hour calculations
  select(PoolID, Foundation_spp, Before_After, Removal_Control, Day_Night, NEC.mmol.m2.hr, NEP.mmol.m2.hr) %>%
  clean_names()%>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))
  
######

head(microbedata)
head(biogeodata)

data_all<-left_join(biogeodata, microbedata, by = "Id_code") %>%
  dplyr::select(PoolID, Foundation_spp, Day_Night, Sampling_group, Before_After, Removal_Control, Group, Time_Point,
         Sampling_Day, Sampling_time, DO_mg_L, Salinity, Temp.pool,PO_umol_L, NN_umol_L,
         NH4_umol_L, ph = pH_insitu, DIC, ALK, OmegaAragonite, `Heterotrophic Bacterioplankton/μL`,`Synechoococcus/μL`,
         `Autotrophic PicoEukaryotes/μL`,`M:C`,BIX,HIX,FI, `Ultra Violet Humic-like`,
         `Marine Humic-like`,`Visible Humic-like`,`Tryptophan-like`,`Tyrosine-like`,`Phenylalanine-like`
) %>%
  clean_names()%>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))


### export the joined data file ####
write_csv(data_all, here("Data","Microbe_Clean","joinedData.csv"))

### Make a few visuals ####


data_all<-data_all %>%
#  filter(sampling_day !=ymd("2019-07-11"), time_point != 5)%>%
  group_by(before_after, time_point, day_night, sampling_group)%>%
  mutate(do_change = do_mg_l - do_mg_l[foundation_spp == "Ocean"],
         temp_change = temp_pool - temp_pool[foundation_spp == "Ocean"],
         po_change = po_umol_l - po_umol_l[foundation_spp == "Ocean"],
         nn_change = nn_umol_l - nn_umol_l[foundation_spp == "Ocean"],
         nh4_change = nh4_umol_l - nh4_umol_l[foundation_spp == "Ocean"],
         ph_change = ph - ph[foundation_spp == "Ocean"],
         hetero_change = heterotrophic_bacterioplankton_m_l - heterotrophic_bacterioplankton_m_l[foundation_spp == "Ocean"],
         syn_change = synechoococcus_m_l - synechoococcus_m_l[foundation_spp == "Ocean"],
         mc_change = m_c - m_c[foundation_spp == "Ocean"],
         bix_change = bix - bix[foundation_spp == "Ocean"],
         hix_change = hix - hix[foundation_spp == "Ocean"],
         fi_change = fi - fi[foundation_spp == "Ocean"],
         uv_humic_change = ultra_violet_humic_like - ultra_violet_humic_like[foundation_spp == "Ocean"],
         marine_change = marine_humic_like - marine_humic_like[foundation_spp == "Ocean"],
         visible_change = visible_humic_like - visible_humic_like[foundation_spp == "Ocean"],
         trypto_change = tryptophan_like - tryptophan_like[foundation_spp == "Ocean"],
         tyrosine_change = tyrosine_like - tyrosine_like[foundation_spp == "Ocean"],
         phenyl_change = phenylalanine_like - phenylalanine_like[foundation_spp == "Ocean"]
  ) 
  

# boxplots showing only control pools for timepoints 1 and 4 in the after period.  Values are differences from the ocean   
data_all  %>%
  filter(foundation_spp !="Ocean", before_after !="Before", 
         #time_point !=5, 
         removal_control == "Control") %>%
  dplyr::select(pool_id, foundation_spp, day_night, time_point, do_change:phenyl_change)%>%
  mutate(nh4_change = ifelse(nh4_change>200, NA, nh4_change))%>% # drop the 2 crazy outliers
  pivot_longer(cols = do_change:phenyl_change, names_to = "parameter", values_to = "values") %>%
  mutate(found_time = paste(foundation_spp,time_point))%>%
  ggplot(aes(x = found_time, y = values, fill = day_night))+
  geom_boxplot()+
  labs(y = "difference from ocean sample")+
  facet_wrap(~parameter, scales = "free")


pcadata<-data_all  %>%
  filter(foundation_spp !="Ocean", before_after !="Before", time_point !=5, removal_control == "Control") %>%
  ungroup()%>%
  dplyr::select(do_change:phenyl_change) 
  
# run a pca
pca<-prcomp(pcadata, scale. = TRUE, center = TRUE)

# calculate percent explained by each PC
perc.explained<-round(100*pca$sdev/sum(pca$sdev),1)

# Extract the scores and loadings
PC_scores <-as_tibble(pca$x[,1:2])

PC_loadings<-as_tibble(pca$rotation)%>%
  bind_cols(labels = rownames(pca$rotation))

# compare timepoint 1 to 4
# bind the data together
data_pca<-data_all  %>%
  filter(foundation_spp !="Ocean", before_after !="Before", time_point !=5, removal_control == "Control") %>%
  ungroup()%>%
  bind_cols(PC_scores)


p1<-data_pca %>%
  mutate(time_point = factor(time_point))%>%
  ggplot(aes(x = PC1, y = PC2, color = foundation_spp, shape = time_point))+
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
 # scale_shape_manual(values = c(1, 22,15,16))+
  scale_colour_manual(values = c("#D64550","#EA9E8D"))+
  scale_fill_manual(values = c("#D64550","#EA9E8D"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(fill = foundation_spp, label = paste(time_point, foundation_spp), color =foundation_spp),
    alpha = .35, show.legend = FALSE,  label.buffer = unit(1, "mm"), con.cap=0, tol = 0.05)+
  geom_point(size = 2) +
  labs(
     x = paste0("PC1 ","(",perc.explained[1],"%)"),
     y = paste0("PC2 ","(",perc.explained[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_blank(),
  #      strip.text = element_blank()
        )+
  facet_wrap(~day_night)


p2<-PC_loadings %>%
  ggplot(aes(x=PC1, y=PC2, label=labels))+
  geom_richtext(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA) +
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),size = 1.2,
               arrow=arrow(length=unit(0.1,"cm")))+
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(
       y = "",
       x = "")+
  #y = paste0("PC2 ","(",perc.explained_both[2],"%)"))+
#  scale_color_manual(values = wes_palette("Darjeeling1"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.75, 0.75),
        legend.position = "none",
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))


p1/plot_spacer()/p2+plot_layout(heights = c(2,-0.1, 1))

## only do timepoint 4
pcadata<-data_all  %>%
  filter(foundation_spp !="Ocean", before_after !="Before", time_point ==4, removal_control == "Control") %>%
  ungroup()%>%
  dplyr::select(do_change:phenyl_change) 

# run a pca
pca<-prcomp(pcadata, scale. = TRUE, center = TRUE)

# calculate percent explained by each PC
perc.explained<-round(100*pca$sdev/sum(pca$sdev),1)

# Extract the scores and loadings
PC_scores <-as_tibble(pca$x[,1:2])

PC_loadings<-as_tibble(pca$rotation)%>%
  bind_cols(labels = rownames(pca$rotation))

# bind the data together
data_pca<-data_all  %>%
  filter(foundation_spp !="Ocean", before_after !="Before", time_point ==4, removal_control == "Control") %>%
  ungroup()%>%
  bind_cols(PC_scores)


p1<-data_pca %>%
  ggplot(aes(x = PC1, y = PC2, color = foundation_spp, shape = day_night))+
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  # scale_shape_manual(values = c(1, 22,15,16))+
  scale_colour_manual(values = c("#D64550","#EA9E8D"))+
  scale_fill_manual(values = c("#D64550","#EA9E8D"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(fill = foundation_spp, label = paste(day_night, foundation_spp), color =foundation_spp),
    alpha = .35, show.legend = FALSE,  label.buffer = unit(1, "mm"), con.cap=0, tol = 0.05)+
  geom_point(size = 2) +
  labs(
    x = paste0("PC1 ","(",perc.explained[1],"%)"),
    y = paste0("PC2 ","(",perc.explained[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_blank(),
        #      strip.text = element_blank()
  )


p2<-PC_loadings %>%
  ggplot(aes(x=PC1, y=PC2, label=labels))+
  geom_richtext(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA) +
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),size = 1.2,
               arrow=arrow(length=unit(0.1,"cm")))+
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(
    y = "",
    x = "")+
  #y = paste0("PC2 ","(",perc.explained_both[2],"%)"))+
  #  scale_color_manual(values = wes_palette("Darjeeling1"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.75, 0.75),
        legend.position = "none",
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))


p1/p2

### calculate growth and uptake rates ####
#### There is a mix between samples that were collected at timepoint 4 and 5 as the endpoint and most missing samples are from the "before" period.

### add in the SA and vol meta data
data_all<-data_all %>%
  left_join(MetaData %>%
              clean_names() %>%
              mutate(pool_id = as.character(pool_id)))


## plot the raw ocean data of before and after
data_all %>%
  ungroup()%>%
  filter(foundation_spp == "Ocean") %>%
  filter(day_night == "Day") %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)"))) %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like) %>%
  select(month, do_mg_l,heterotrophic_bac_m_l=heterotrophic_bacterioplankton_m_l,prot,temperature = temp_pool, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  group_by(month, name)%>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val= sd(value, na.rm = TRUE)/sqrt(n())) %>%
  ggplot(aes(x = month, y = mean_val))+
  geom_point(size = 3)+
  geom_errorbar(aes(x = month, ymin = mean_val - se_val, ymax = mean_val+se_val), width = 0.1)+
  facet_wrap(~name, scales = "free_y", strip.position = "left", ncol = 1)+
  labs(x = "",
       y = "",
       title = "Mean Ocean values")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside")
  
ggsave(here("Output","OceanVals.png"), width = 4, height = 9)

# Before Day nutrients were collected at 1,4, and 5 and microbes at 1,5
# Before night, nutrients are 1,4, and microbes at 1,5
# After day and night both sets are at 1 and 4

data_rates_before<-data_all%>%
  ungroup()%>%
  mutate(sampling_datetime = mdy_hms(paste(sampling_day, sampling_time))) %>% 
  filter(before_after == "Before")%>%
  reframe(difftime_hours = as.numeric(sampling_datetime[time_point==4]-sampling_datetime[time_point==1]),
          difftime_hours5 = as.numeric(sampling_datetime[time_point==5]-sampling_datetime[time_point==1]),
          do_mg_l_rate = (do_mg_l[time_point==4]-do_mg_l[time_point == 1])/difftime_hours,
          temp_rate = (temp_pool[time_point==4]-temp_pool[time_point == 1])/difftime_hours,
          po_rate = (po_umol_l[time_point==4]-po_umol_l[time_point == 1])/difftime_hours,
          nn_rate = (nn_umol_l[time_point==4]-nn_umol_l[time_point == 1])/difftime_hours,
          nh4_rate = (nh4_umol_l[time_point==4]-nh4_umol_l[time_point == 1])/difftime_hours,
          ph_rate = (ph[time_point==4]-ph[time_point == 1])/difftime_hours,
          hetero_rate = (heterotrophic_bacterioplankton_m_l[time_point==5]-heterotrophic_bacterioplankton_m_l[time_point == 1])/difftime_hours5,
          syn_rate = (synechoococcus_m_l[time_point==5]-synechoococcus_m_l[time_point == 1])/difftime_hours5,
          auto_rate = (autotrophic_pico_eukaryotes_m_l[time_point==5]-autotrophic_pico_eukaryotes_m_l[time_point == 1])/difftime_hours5,
          uv_rate = (ultra_violet_humic_like[time_point==5] - ultra_violet_humic_like[time_point == 1])/difftime_hours5,
          mh_rate = (marine_humic_like[time_point==5] - marine_humic_like[time_point == 1])/difftime_hours5,
          vh_rate = (visible_humic_like[time_point==5] - visible_humic_like[time_point == 1])/difftime_hours5,
          tryp_rate = (tryptophan_like[time_point==5] - tryptophan_like[time_point == 1])/difftime_hours5,
          tyro_rate = (tyrosine_like[time_point==5] - tyrosine_like[time_point == 1])/difftime_hours5,
          phenyl_rate = (phenylalanine_like[time_point==5] - phenylalanine_like[time_point == 1])/difftime_hours5,
          mc_rate = (m_c[time_point==5] - m_c[time_point == 1])/difftime_hours5,
          bix_rate = (bix[time_point==5] - bix[time_point == 1])/difftime_hours5,
          hix_rate = (hix[time_point==5] - hix[time_point == 1])/difftime_hours5,
          fi_rate = (fi[time_point==5] - fi[time_point == 1])/difftime_hours5,
          .by = c(foundation_spp, removal_control, day_night, pool_id, sampling_group))%>%
  drop_na(foundation_spp) %>%
  mutate(before_after = "Before")




data_rates_after<-data_all%>%
  ungroup()%>%
  mutate(sampling_datetime = mdy_hms(paste(sampling_day, sampling_time)),
         humics = ultra_violet_humic_like+marine_humic_like+ visible_humic_like,
         prot = tryptophan_like + tyrosine_like + phenylalanine_like) %>% 
  filter(before_after == "After")%>%
  reframe(difftime_hours = as.numeric(sampling_datetime[time_point==4]-sampling_datetime[time_point==1]),
          do_mg_l_rate = (do_mg_l[time_point==4]-do_mg_l[time_point == 1])/difftime_hours,
          temp_rate = (temp_pool[time_point==4]-temp_pool[time_point == 1])/difftime_hours,
          po_rate = (po_umol_l[time_point==4]-po_umol_l[time_point == 1])/difftime_hours,
          nn_rate = (nn_umol_l[time_point==4]-nn_umol_l[time_point == 1])/difftime_hours,
          nh4_rate = (nh4_umol_l[time_point==4]-nh4_umol_l[time_point == 1])/difftime_hours,
          ph_rate = (ph[time_point==4]-ph[time_point == 1])/difftime_hours,
          hetero_rate = (heterotrophic_bacterioplankton_m_l[time_point==4]-heterotrophic_bacterioplankton_m_l[time_point == 1])/difftime_hours,
          syn_rate = (synechoococcus_m_l[time_point==4]-synechoococcus_m_l[time_point == 1])/difftime_hours,
          auto_rate = (autotrophic_pico_eukaryotes_m_l[time_point==4]-autotrophic_pico_eukaryotes_m_l[time_point == 1])/difftime_hours,
          humics_rate = (humics[time_point==4] - humics[time_point == 1])/difftime_hours,
          prot_rate = (prot[time_point==4] - prot[time_point == 1])/difftime_hours,
          uv_rate = (ultra_violet_humic_like[time_point==4] - ultra_violet_humic_like[time_point == 1])/difftime_hours,
          mh_rate = (marine_humic_like[time_point==4] - marine_humic_like[time_point == 1])/difftime_hours,
          vh_rate = (visible_humic_like[time_point==4] - visible_humic_like[time_point == 1])/difftime_hours,
          tryp_rate = (tryptophan_like[time_point==4] - tryptophan_like[time_point == 1])/difftime_hours,
          tyro_rate = (tyrosine_like[time_point==4] - tyrosine_like[time_point == 1])/difftime_hours,
          phenyl_rate = (phenylalanine_like[time_point==4] - phenylalanine_like[time_point == 1])/difftime_hours,
          mc_rate = (m_c[time_point==4] - m_c[time_point == 1])/difftime_hours,
          bix_rate = (bix[time_point==4] - bix[time_point == 1])/difftime_hours,
          hix_rate = (hix[time_point==4] - hix[time_point == 1])/difftime_hours,
          fi_rate = (fi[time_point==4] - fi[time_point == 1])/difftime_hours,
          .by = c(foundation_spp, removal_control, day_night, pool_id, sampling_group),)%>%
  drop_na(foundation_spp)%>%
  mutate(before_after = "After")

# bring before and after together
data_rates <-
  bind_rows(data_rates_before, data_rates_after) %>%
  mutate(nh4_rate = ifelse(nh4_rate>10, NA, nh4_rate)) # remove crazy Nh4 outlier


## add in NEC and NEP
data_rates<-data_rates %>%
  left_join(MetabRates %>%
              mutate(pool_id = as.character(pool_id)))

### calculate true rates normalized to tide pool volume and SA while also controling for the change in the ocean
ocean_rates<-data_rates %>%
  filter(pool_id == "Ocean") %>%
  select(day_night, before_after, sampling_group, do_mg_l_rate:auto_rate, humics_rate: fi_rate) %>%
  rename_with(~paste0(., "_ocean"), -c(day_night:sampling_group)) # rename all the rates to say ocean

# bring to ocean rates with the pool rates for easier normalization
# calculate rates normalized to ocean change, tide pool volume and SA 

poolrates<-data_rates%>%
  filter(pool_id != "Ocean") %>%
  select(-c(difftime_hours, difftime_hours5)) %>% # remove things I dont need
  left_join(ocean_rates) %>% # join with the ocean rates
  left_join(MetaData %>%
              clean_names() %>%
              mutate(pool_id = as.character(pool_id)))%>%
  mutate(do_mg_m2_hr_ocean = (do_mg_l_rate - do_mg_l_rate_ocean)*(vol/surface_area), # normalized to ocean
         do_mg_m2_hr = (do_mg_l_rate)*(vol/surface_area), # not ocean
         po_umol_m2_hr_ocean = (po_rate  - po_rate_ocean)*(vol/surface_area),
         po_umol_m2_hr = (po_rate)*(vol/surface_area),
         nn_umol_m2_hr_ocean = (nn_rate  - nn_rate_ocean)*(vol/surface_area),
         nn_umol_m2_hr = (nn_rate)*(vol/surface_area),
         nh4_umol_m2_hr_ocean = (nh4_rate  - nh4_rate_ocean)*(vol/surface_area),
         nh4_umol_m2_hr = (nh4_rate)*(vol/surface_area),
         ph_m2_hr_ocean = (ph_rate  - ph_rate_ocean)/(surface_area),
         ph_m2_hr = (ph_rate/surface_area),
         hetero_counts_m2_hr_ocean = (hetero_rate  - hetero_rate_ocean)*(1000*vol/surface_area),# bacteria are in counts per ml
         hetero_counts_m2_hr = (hetero_rate)*(vol/surface_area),
         syn_counts_m2_hr_ocean = (syn_rate  - syn_rate_ocean)*(1000*vol/surface_area),# bacteria are in counts per ml
         syn_counts_m2_hr = (syn_rate)*(vol/surface_area),
         auto_counts_m2_hr_ocean = (auto_rate  - auto_rate_ocean)*(1000*vol/surface_area),# bacteria are in counts per ml
         auto_counts_m2_hr = (auto_rate)*(vol/surface_area),
         prot_raman_m2_hr_ocean = (prot_rate  - prot_rate_ocean)/surface_area,# 
         prot_raman_m2_hr = (prot_rate/surface_area),
         humics_raman_m2_hr_ocean = (humics_rate  - humics_rate_ocean)/surface_area,# 
         humics_raman_m2_hr = (humics_rate/surface_area),
         uv_raman_m2_hr_ocean = (uv_rate  - uv_rate_ocean)/surface_area,# 
         uv_raman_m2_hr = (uv_rate/surface_area),
         mh_raman_m2_hr_ocean = (mh_rate  - mh_rate_ocean)/surface_area,# 
         mh_raman_m2_hr = (mh_rate/surface_area),
         vh_raman_m2_hr_ocean = (vh_rate  - vh_rate_ocean)/surface_area,# 
         vh_raman_m2_hr = (vh_rate/surface_area),
         tryp_raman_m2_hr_ocean = (tryp_rate  - tryp_rate_ocean)/surface_area,# 
         tryp_raman_m2_hr = (tryp_rate/surface_area),
         tyro_raman_m2_hr_ocean = (tyro_rate  - tyro_rate_ocean)/surface_area,# 
         tyro_raman_m2_hr = (tyro_rate/surface_area),
         phenyl_raman_m2_hr_ocean = (phenyl_rate  - phenyl_rate_ocean)/surface_area,# 
         phenyl_raman_m2_hr = (phenyl_rate/surface_area),
         mc_hr_ocean = (mc_rate - mc_rate_ocean),# this is a ratio so not normalized to anything
         mc_hr = mc_rate,
         bix_hr_ocean = (bix_rate - bix_rate_ocean),# this is a ratio so not normalized to anything
         bix_hr = bix_rate,
         hix_hr_ocean = (hix_rate - hix_rate_ocean),# this is a ratio so not normalized to anything
         hix_hr = hix_rate,
         fi_hr_ocean = (fi_rate - fi_rate_ocean),# this is a ratio so not normalized to anything
         fi_hr = fi_rate
  ) %>%
  select(foundation_spp:pool_id,before_after,do_mg_m2_hr_ocean:fi_hr)
    

## Add in the eco metab rates to the pool rates
poolrates <- poolrates %>%
  left_join(MetabRates %>%
              mutate(pool_id = as.character(pool_id)))


# make a bunch of boxplots
data_rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  pivot_longer(cols= do_mg_l_rate:auto_rate) %>%
  dplyr::select(-difftime_hours, -difftime_hours5, -sampling_group) %>%
  ggplot(aes(x = before_after, y = value, fill = removal_control))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_boxplot()+
  facet_wrap(name*foundation_spp~day_night, scales = "free")

## Make some reaction norms
mean_rates<-data_rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  pivot_longer(cols= c(do_mg_l_rate:auto_rate, humics_rate,prot_rate)) %>%
 # filter(!pool_id %in% 1:16)%>% ## uneven sample sizes... only 16 pools have everything
  group_by(foundation_spp, removal_control, day_night, before_after, name)%>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE)/sqrt(n()))


mean_poolrates<- poolrates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
 # select(-contains("_ocean")) %>%# remove the ocean control
  pivot_longer(cols= c(do_mg_m2_hr:fi_hr)) %>%
  group_by(foundation_spp, removal_control, day_night, before_after, name)%>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE)/sqrt(n()))


mean_poolrates %>%
  filter(day_night != "Night") %>%
  filter(name %in% c("auto_counts_m2_hr_ocean","syn_counts_m2_hr_ocean",
                   "hetero_counts_m2_hr_ocean","do_mg_m2_hr_ocean",
                   "nh4_umol_m2_hr_ocean","nn_umol_m2_hr_ocean")) %>%
  ggplot(aes(x = before_after, y = mean_value, color = removal_control, group = removal_control))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point()+
  geom_errorbar(aes(x = before_after, y = mean_value, ymin = mean_value-se_value, ymax = mean_value+se_value), width = 0.01)+
  geom_line()+
  labs(x = "time period",
       y = "Rate per hour")+
  facet_wrap(name~foundation_spp, scales = "free", ncol = 2)+
  theme_bw()
  


## with only the control pools, but across both timepoints
mean_rates_control<-data_rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  filter(removal_control != "Removal")%>%
  pivot_longer(cols= c(do_mg_l_rate:auto_rate, humics_rate,prot_rate)) %>%
  # filter(!pool_id %in% 1:16)%>% ## uneven sample sizes... only 16 pools have everything
  group_by(foundation_spp, day_night, name)%>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE)/sqrt(n()))


# rename the foundation species for the ocean to have one group for mussels and one for phyllo for easier plotting
myt_ocean<-mean_rates %>%
  filter(foundation_spp == "Ocean") %>%
  mutate(foundation_spp = "Mytilus")

phylo_ocean<-mean_rates %>%
  filter(foundation_spp == "Ocean") %>%
  mutate(foundation_spp = "Phyllospadix")

mean_plot<-mean_rates %>%
  filter(foundation_spp != "Ocean")%>%
  bind_rows(myt_ocean) %>%
  bind_rows(phylo_ocean) %>%
mutate(nicenames = case_when(
  name == "auto_rate" ~"Autotrophic Bacteria <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
  name == "do_mg_l_rate" ~ "DO <br> (mg L<sup>-1</sup> hr<sup>-1</sup>)",
  name == "hetero_rate" ~ "Heterotrophic Bacteria <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
  name == "nh4_rate" ~ "Ammonium <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
  name == "nn_rate" ~ "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>hr<sup>-1</sup>)"))

# make a set of reaction norm plots for the day


day_plot<-mean_plot %>%
  filter(day_night == "Day")%>%
  filter(name %in% c(
         "hetero_rate","do_mg_l_rate",
         "nh4_rate","nn_rate"))%>%
  mutate(nicenames = factor(nicenames, levels = c( "DO <br> (mg L<sup>-1</sup> hr<sup>-1</sup>)","Heterotrophic Bacteria <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                                                   "Ammonium <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)","Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>hr<sup>-1</sup>)")))%>%
  ggplot(aes(x = before_after, y = mean_value, color = removal_control, group = removal_control, shape = removal_control))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(x = before_after, y = mean_value, ymin = mean_value-se_value, ymax = mean_value+se_value), width = 0.01)+
  geom_line()+
  labs(x = "",
       y = "", 
       color = "",
       shape = "")+
  scale_color_manual(values = c("grey30","#79ACBD","grey30"))+
  scale_shape_manual(values = c(16,16,1))+
  #facet_nested_wrap(vars(nicenames,foundation_sp), scales = "free", ncol = 2)
  facet_wrap(nicenames~foundation_spp, scales = "free_y", ncol = 2, strip.position = "left")+
  facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits=c(-1, 5)),
      scale_y_continuous(limits = c(-75, 150)),
      scale_y_continuous(limits = c(-1, 4)),
      scale_y_continuous(limits = c(-1.5, 2))
      ), each = 2))+
  theme_bw()+
  theme(strip.text.y.left  = ggtext::element_markdown(size = 16),
        strip.text.x = element_blank(),
        strip.placement = "outside",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.background = element_blank())
  
day_plot

ggsave(here("Output","BACI_plot.pdf"), width = 8, height = 8)

## run some two-way ANOVAs

## I need to layer then so that ocean is included as removal control for both M and P

data_anova_M<-data_rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  select(foundation_spp, before_after,removal_control, day_night, pool_id, "auto_rate",
         "hetero_rate","do_mg_l_rate",
         "nh4_rate","nn_rate")%>%
  pivot_longer(cols= c("auto_rate",
                       "hetero_rate","do_mg_l_rate",
                       "nh4_rate","nn_rate")) %>%
  filter(foundation_spp == "Mytilus")

data_anova_P<-data_rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  select(foundation_spp, before_after,removal_control, day_night, pool_id, "auto_rate",
         "hetero_rate","do_mg_l_rate",
         "nh4_rate","nn_rate")%>%
  pivot_longer(cols= c("auto_rate",
                       "hetero_rate","do_mg_l_rate",
                       "nh4_rate","nn_rate")) %>%
  filter(foundation_spp == "Phyllospadix")


data_anova_om <-data_rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  select(foundation_spp, before_after,removal_control, day_night, pool_id, "auto_rate",
         "hetero_rate","do_mg_l_rate",
         "nh4_rate","nn_rate")%>%
  pivot_longer(cols= c("auto_rate",
                       "hetero_rate","do_mg_l_rate",
                       "nh4_rate","nn_rate")) %>%
  filter(foundation_spp == "Ocean")%>%
  mutate(foundation_spp = "Mytilus")

data_anova_op <-data_rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  select(foundation_spp, before_after,removal_control, day_night, pool_id, "auto_rate",
         "hetero_rate","do_mg_l_rate",
         "nh4_rate","nn_rate")%>%
  pivot_longer(cols= c("auto_rate",
                       "hetero_rate","do_mg_l_rate",
                       "nh4_rate","nn_rate")) %>%
  filter(foundation_spp == "Ocean")%>%
  mutate(foundation_spp = "Phyllospadix")

data_anova <-data_anova_M %>%
  bind_rows(data_anova_om)%>%
  bind_rows(data_anova_P)%>%
  bind_rows(data_anova_op)%>%
  filter(day_night == "Day")%>%
  select(!day_night)%>%
  nest_by(foundation_spp, name) %>%
  mutate(mod = list(lmer(value~removal_control*before_after+(1|pool_id), data = data)),
         modstat = list(broom::glance(mod)),
         res =  list(broom::tidy(mod)),
         ano = list(broom::tidy(anova(mod))))

data_anova %>%
  select(foundation_spp, name, modstat, res, ano) %>%
  unnest(ano) %>%
  filter(term != "Residuals") %>%
  filter(p.value <= 0.05) # pull out just the significant values


### do the anova without the ocean sample since too low sample size
data_rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  select(foundation_spp, before_after,removal_control, day_night, pool_id, "auto_rate",
         "hetero_rate","do_mg_l_rate",
         "nh4_rate","nn_rate")%>%
  pivot_longer(cols= c("auto_rate",
                       "hetero_rate","do_mg_l_rate",
                       "nh4_rate","nn_rate")) %>%
  filter(foundation_spp != "Ocean") %>%
  filter(day_night == "Day")%>%
  select(!day_night)%>%
  nest_by(foundation_spp, name) %>%
  mutate(mod = list(lmer(value~removal_control*before_after+(1|pool_id), data = data)),
         modstat = list(broom::glance(mod)),
         res =  list(broom::tidy(mod)),
         ano = list(broom::tidy(anova(mod))))%>%
select(foundation_spp, name, modstat, res, ano) %>%
  unnest(ano) %>%
  filter(term != "Residuals") %>%
  filter(p.value <= 0.05)



# make a set of reaction norm plots for the night

night_plot<-mean_rates %>%
  filter(foundation_spp != "Ocean")%>%
  bind_rows(myt_ocean) %>%
  bind_rows(phylo_ocean) %>%
  filter(day_night == "Night")%>%
  ggplot(aes(x = before_after, y = mean_value, color = removal_control, group = removal_control))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point()+
  geom_errorbar(aes(x = before_after, y = mean_value, ymin = mean_value-se_value, ymax = mean_value+se_value), width = 0.01)+
  geom_line()+
  labs(x = "time period",
       y = "Rate per hour")+
  facet_wrap(name~foundation_spp, scales = "free", ncol = 2)+
  theme_bw()

night_plot




### just look at controls in the after period 
#averaged across both the before and after period for the controls
mean_rates_control %>%
  mutate(name = case_when(name  == "auto_rate"~ "Autotrophic pico counts uL-1 hr-1",
                          name == "do_mg_l_rate"~ "DO mg L-1 hr-1",
                          name == "hetero_rate"~"Heterotrophic counts uL-1 hr-1",
                          name == "nh4_rate"~"NH4 umol L-1 hr-1",
                          name == "nn_rate"~"NN umol L-1 hr-1",
                          name == "ph_rate"~ "pH unit hr-1",
                          name == "po_rate"~"PO umol L-1 hr-1",
                          name == "syn_rate"~"Synechoococcus counts uL-1 hr-1",
                          name == "temp_rate"~"Temperature (deg C) hr-1",
                          name  == "humics_rate"~"Humics",
                          name == "prot_rate"~"Prot"))%>%
  mutate(name = factor(name, levels=c("Heterotrophic counts uL-1 hr-1","Synechoococcus counts uL-1 hr-1",
                                      "Autotrophic pico counts uL-1 hr-1","NN umol L-1 hr-1",
                                      "NH4 umol L-1 hr-1","PO umol L-1 hr-1",
                                      "DO mg L-1 hr-1", "pH unit hr-1", "Temperature (deg C) hr-1",
                                      "Humics","Prot")))%>%
#  filter(removal_control %in% c("Control","Ocean"), before_after == "After") %>%
  ggplot(aes(x = day_night, color = foundation_spp, y = mean_value, group = foundation_spp))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point()+
  geom_errorbar(aes(x = day_night, y = mean_value, ymin = mean_value-se_value, ymax = mean_value+se_value), width = 0.01)+
  geom_line()+
  scale_color_manual(values = c("grey30","#79ACBD","#79B38F"))+
  labs(x = "",
       y = "Mean rate (unit per hour)",
       color = "")+
  facet_wrap(~name, scales = "free", ncol = 3)+
  theme_bw()

ggsave(here("Output","mean_rates_2mo.png"), width = 8, height = 8)


## just the after period alone
mean_rates%>%
  mutate(name = case_when(name  == "auto_rate"~ "Autotrophic pico counts uL-1 hr-1",
                          name == "do_mg_l_rate"~ "DO mg L-1 hr-1",
                          name == "hetero_rate"~"Heterotrophic counts uL-1 hr-1",
                          name == "nh4_rate"~"NH4 umol L-1 hr-1",
                          name == "nn_rate"~"NN umol L-1 hr-1",
                          name == "ph_rate"~ "pH unit hr-1",
                          name == "po_rate"~"PO umol L-1 hr-1",
                          name == "syn_rate"~"Synechoococcus counts uL-1 hr-1",
                          name == "temp_rate"~"Temperature (deg C) hr-1",
                          name  == "humics_rate"~"Humics",
                          name == "prot_rate"~"Prot"))%>%
  mutate(name = factor(name, levels=c("Heterotrophic counts uL-1 hr-1","Synechoococcus counts uL-1 hr-1",
                                      "Autotrophic pico counts uL-1 hr-1","NN umol L-1 hr-1",
                                      "NH4 umol L-1 hr-1","PO umol L-1 hr-1",
                                      "DO mg L-1 hr-1", "pH unit hr-1", "Temperature (deg C) hr-1",
                                      "Humics","Prot")))%>%
  filter(removal_control %in% c("Control","Ocean"), before_after == "After") %>%
  ggplot(aes(x = day_night, color = foundation_spp, y = mean_value, group = foundation_spp))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point()+
  geom_errorbar(aes(x = day_night, y = mean_value, ymin = mean_value-se_value, ymax = mean_value+se_value), width = 0.01)+
  geom_line()+
  scale_color_manual(values = c("grey30","#79ACBD","#79B38F"))+
  labs(x = "",
       y = "Mean rate (unit per hour)",
       color = "")+
  facet_wrap(~name, scales = "free", ncol = 3)+
  theme_bw()

ggsave(here("Output","mean_rates_1mo.png"), width = 8, height = 8)


### PCA of fDOM from just the control and after period ####

pca_after<-data_all %>%
  filter(removal_control != "Removal", before_after == "After") %>%
  filter(tryptophan_like < 1.3)%>% # there is a big mussel outlier
  ungroup()%>%
  mutate(time_point = factor(time_point))%>%
  #dplyr::select(foundation_spp, day_night, time_point, ultra_violet_humic_like, marine_humic_like, visible_humic_like, tryptophan_like, tyrosine_like, phenylalanine_like)%>%
  dplyr::select(foundation_spp, day_night,time_point, m_c, bix, hix, fi,ultra_violet_humic_like, marine_humic_like, visible_humic_like, tryptophan_like, tyrosine_like, phenylalanine_like)%>%
  drop_na()

# run a pca
pca<-prcomp(pca_after[,4:13], scale. = TRUE, center = TRUE)

# calculate percent explained by each PC
perc.explained<-round(100*pca$sdev/sum(pca$sdev),1)

# Extract the scores and loadings
PC_scores <-as_tibble(pca$x[,1:2])

PC_loadings<-as_tibble(pca$rotation)%>%
  bind_cols(labels = rownames(pca$rotation))

# bind the data together
data_pca_after<-pca_after  %>%
    bind_cols(PC_scores)

p1<-data_pca_after %>%
  ggplot(aes(x = PC1, y = PC2, color = foundation_spp, shape = time_point))+
  coord_cartesian(xlim = c(-12, 12), ylim = c(-10, 10)) +
  # scale_shape_manual(values = c(1, 22,15,16))+
  scale_color_manual(values = c("grey30","#79ACBD","#79B38F"))+
  scale_fill_manual(values = c("grey30","#79ACBD","#79B38F"))+
    geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(fill = foundation_spp, label = paste(time_point), color =foundation_spp),
    alpha = .35, show.legend = FALSE,  label.buffer = unit(1, "mm"), con.cap=0, tol = 0.05)+
  geom_point(size = 2) +
  labs(
    x = "",
    #x = paste0("PC1 ","(",perc.explained[1],"%)"),
    y = paste0("PC2 ","(",perc.explained[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_blank(),
        #      strip.text = element_blank()
  )+
  facet_wrap(day_night~foundation_spp)

# loadings
p2<-PC_loadings %>%
  ggplot(aes(x=PC1, y=PC2, label=labels))+
  geom_richtext(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA) +
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),size = 1.2,
               arrow=arrow(length=unit(0.1,"cm")))+
  coord_cartesian(xlim = c(-12, 12), ylim = c(-10, 10)) +
  labs(
    y = "",
    x = paste0("PC1 ","(",perc.explained[1],"%)"))+
  #y = paste0("PC2 ","(",perc.explained_both[2],"%)"))+
  #  scale_color_manual(values = wes_palette("Darjeeling1"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.75, 0.75),
        legend.position = "none",
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))



p1/plot_spacer()/p2+plot_layout(heights = c(2,-0.1, 1))
ggsave(here("Output","pca_fdom.png"), width = 8, height = 8)



## Make a plot of the benthic data

BenthicData %>%
  select(Foundation_spp, Removal_Control, PoolID, Before_After, macroalgae, Diatoms, consumers, allCCA, AdjMusselCover, AdjSurfgrassCover)%>%
  mutate(rock_sand = (100-(macroalgae+consumers+allCCA+AdjMusselCover+AdjSurfgrassCover)),
         macroalgae = macroalgae - Diatoms,
         Before_After = factor(Before_After, levels = c("Before","After"))) %>% # the macroalgae includes diatom cover and I want to see the difference
  pivot_longer(cols = macroalgae:rock_sand) %>%
  ggplot(aes(x = Before_After, y = value, fill = name))+
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cal_palette("tidepool", n = 7, type = "continuous"))+
  facet_wrap(Before_After~Removal_Control, scale = "free_x")+
  theme_classic()+
  facet_wrap(Removal_Control~PoolID, ncol = 8)

BenthicData %>%
  filter(Before_After == "After")%>%
  select(PoolID, Before_After, macroalgae, Diatoms, consumers, allCCA, AdjMusselCover, AdjSurfgrassCover)%>%
  mutate(PoolID =  as.factor(PoolID)) %>%
  mutate(rock_sand = (100-(macroalgae+consumers+allCCA+AdjMusselCover+AdjSurfgrassCover)),
         macroalgae = macroalgae - Diatoms,
         PoolID = fct_reorder2(PoolID,rock_sand, AdjMusselCover),) %>% # the macroalgae includes diatom cover and I want to see the difference
  pivot_longer(cols = macroalgae:rock_sand) %>%
  ggplot(aes(x = PoolID, y = value, fill = name))+
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cal_palette("tidepool", n = 7, type = "continuous"))+
  theme_classic()

#### Take a regression approach with and w/o the ocean normalization
Benthic<-BenthicData %>%
  select(-c(Removal_Control, Foundation_spp))%>%
  rename(pool_id = PoolID, before_after = Before_After) %>%
  mutate(pool_id = as.character(pool_id))

PoolRates_combined <-poolrates %>%
  left_join(Benthic, by = c("pool_id", "before_after"))             


PoolRates_combined %>%
  filter(day_night == "Day",
         before_after == "After") %>%
  ggplot(aes(x = allproddom, y = do_mg_m2_hr_ocean))+
  geom_point()
  

PoolRates_combined %>%
  filter(day_night == "Day",
         before_after == "After") %>%
  ggplot(aes(x = allproddom, y = nh4_umol_m2_hr_ocean))+
  geom_point()


PoolRates_combined %>%
  filter(day_night == "Day",
         before_after == "After",
         mc_hr > -0.3) %>%
  ggplot(aes(x = do_mg_m2_hr, y = mc_hr))+
  geom_point()+
  geom_smooth(method = "lm")

anova(lm(data = PoolRates_combined %>%
           filter(day_night == "Day",
                  before_after == "After",
                  mc_hr > -0.3), mc_hr~do_mg_m2_hr))

PoolRates_combined %>%
  filter(day_night == "Day",
         before_after == "After",
         mc_hr > -0.3) %>%
  ggplot(aes(x =  tryp_raman_m2_hr, y = hetero_counts_m2_hr))+
  geom_point()+
  geom_smooth(method = "lm")

anova(lm(data = PoolRates_combined %>%
           filter(day_night == "Day",
                  before_after == "After",
                  mc_hr > -0.3), tryp_raman_m2_hr~hetero_counts_m2_hr))

### Take a regression approach #####

##### STOPPED HERE!####

## combine with community comp data

Day_rates<-BenthicData  %>%
  rename(pool_id = PoolID, removal_control = Removal_Control, before_after = Before_After, foundation_spp = Foundation_spp) %>%
  mutate(pool_id = as.character(pool_id))%>%
  right_join(data_rates) %>%
  filter(day_night == "Day") %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)")))

Night_rates<-BenthicData  %>%
  rename(pool_id = PoolID, removal_control = Removal_Control, before_after = Before_After, foundation_spp = Foundation_spp) %>%
  mutate(pool_id = as.character(pool_id))%>%
  right_join(data_rates) %>%
  filter(day_night == "Night",
         do_mg_l_rate > -4) %>% # one crazy outlier
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)")))

P_Nuts<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  filter(month == "August (Upwelling)")%>%
  filter(po_rate> -2)%>% # remove one outlier
 # filter(foundation_spp == "Phyllospadix")%>%
#  select(month,prodphyllodom,nn_rate,nh4_rate, po_rate)%>%
  pivot_longer(cols = c(nn_rate,nh4_rate, hetero_rate, do_mg_l_rate)) %>%
 # filter(removal_control == "Control")%>%
ggplot(aes(x = prodphyllodom, y = value))+
  geom_vline(aes(xintercept  = 0))+
  geom_hline(aes(yintercept  = 0))+
  geom_point(aes(color =  month))+
  geom_smooth(method = "lm")+
  labs(x = "Producer-dominance",
       y = "umol L-1 hr-1",
       color = "")+
  theme_bw()+
  facet_wrap(~name, scales = "free_y")

## Make plot just for prod and DO, and Het and DO
weird <- scales::trans_new("signed_log",
                           transform=function(x) sign(x)*log(abs(x)),
                           inverse=function(x) sign(x)*exp(abs(x)))

P_DO_Prod<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  filter(month == "August (Upwelling)")%>%
  ggplot(aes(x = allproddom, y = do_mg_l_rate))+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_hline(aes(yintercept  = log(1)), lty = 2)+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", color = "black")+
#  coord_trans(y = function(x)log(x+1))+
  geom_text(aes(x = 50, y = 6, label = "p <0.001"))+
  geom_text(aes(x = 50, y = 5.5, label = "R2 = 0.31"))+
  geom_text(aes(x = 50, y = -0.5, label = "Producer dominated"))+
  geom_text(aes(x = -50, y = -0.5, label = "Consumer dominated"))+
  labs(x = "Producer-dominance (% Producers - % Consumers)",
       y = "Productivity (DO mg L-1 hr-1)",
       color = "")+
  theme_bw()+
  theme(panel.grid = element_blank())

anova(lm(data = Day_rates %>%  
           filter(foundation_spp != "Ocean") %>%
           filter(month == "August (Upwelling)"),do_mg_l_rate~ allproddom))

summary(lm(data = Day_rates %>%  
     filter(foundation_spp != "Ocean") %>%
     filter(month == "August (Upwelling)"),do_mg_l_rate~ allproddom))

P_DO_het<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  filter(month == "August (Upwelling)")%>%
  ggplot(aes(x = hetero_rate, y = do_mg_l_rate))+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", color = "black", lty = 2)+
  geom_text(aes(x = 50, y = 6, label = "p = 0.19"))+
  labs(x = "Heterotrophic Bacteria production (counts L-1 hr-1)",
       y = "Productivity (DO mg L-1 hr-1)",
       color = "")+
  theme_bw()+
  theme(panel.grid = element_blank())

P_DO_Prod+P_DO_het
ggsave(here("Output","Community_DO.png"), width = 8, height = 4)

Prod_DO_mod<-lm(data = Day_rates %>%  
                filter(foundation_spp != "Ocean",
                       before_after == "After"),do_mg_l_rate~hetero_rate+allproddom )
anova(Prod_DO_mod)
summary(Prod_DO_mod)

# mod of producer dominance and nuts
NN_prod_mod<-lm(data = Day_rates %>%  
                 filter(foundation_spp != "Ocean")%>%
                   filter(before_after == "After"),nn_rate~poly(allproddom,2))
anova(NN_prod_mod)
summary(NN_prod_mod)

NH4_prod_mod<-lm(data = Day_rates %>%  
                  filter(foundation_spp != "Ocean")%>%
                    filter(before_after == "After"),nh4_rate~poly(prodphyllodom,2))
anova(NH4_prod_mod)

# producers and H back
Hetero_prod_mod<-lmer(data = Day_rates %>%  
                   filter(foundation_spp != "Ocean"),hetero_rate~prodphyllodom+(1|pool_id)  )
anova(Hetero_prod_mod)

Nuts_prod<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
#  filter(nh4_rate < 4)%>% # remove one outlier
 # select(month,nn_rate,nh4_rate, po_rate, hetero_rate, syn_rate, auto_rate)%>%
  rename(`N+N` = nn_rate, NH4 = nh4_rate )%>%
  pivot_longer(cols = c(`N+N`,NH4)) %>%
   filter(before_after == "After")%>%
  ggplot(aes(y = value, x = allproddom))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", color = "black", formula = 'y~poly(x,2)')+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  labs(x = "Producer-dominance (% Producers - % Consumers)",
       y = "umol L-1 hr-1",
       color = "")+
  theme_bw()+
  facet_wrap(~name, ncol = 1, scales = "free_y")+
  theme(panel.grid = element_blank())


NN_prod<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  #  filter(nh4_rate < 4)%>% # remove one outlier
  # select(month,nn_rate,nh4_rate, po_rate, hetero_rate, syn_rate, auto_rate)%>%
  rename(`N+N` = nn_rate, NH4 = nh4_rate )%>%
#  pivot_longer(cols = c(`N+N`,NH4)) %>%
  filter(before_after == "After")%>%
  ggplot(aes(y = `N+N`, x = allproddom))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", color = "black", formula = 'y~poly(x,2)')+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_text(aes(x = -50, y = 2.5, label = "p = 0.01"))+
  labs(x = "Producer-dominance (% Producers - % Consumers)",
       y = "N+N (umol L-1 hr-1)",
       color = "Producer-dominance (% Producers - % Consumers)")+
  theme_bw()+
  theme(panel.grid = element_blank())

NH4_prod<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  #  filter(nh4_rate < 4)%>% # remove one outlier
  # select(month,nn_rate,nh4_rate, po_rate, hetero_rate, syn_rate, auto_rate)%>%
  rename(`N+N` = nn_rate, NH4 = nh4_rate )%>%
  #  pivot_longer(cols = c(`N+N`,NH4)) %>%
  filter(before_after == "After")%>%
  ggplot(aes(y = NH4, x = allproddom))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", color = "black", formula = 'y~poly(x,2)')+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_text(aes(x = -50, y = 4, label = "p < 0.001"))+
  labs(x = "Producer-dominance (% Producers - % Consumers)",
       y = "NH4 (umol L-1 hr-1)",
       color = "")+
  theme_bw()+
  theme(panel.grid = element_blank())

NN_prod/NH4_prod

ggsave(here("Output","Prod_N.png"), height = 8, width = 4)
# run a model of NN and NH4~ heterotrophic Bac
NN_DO_mod<-lmer(data = Day_rates %>%  
                 filter(foundation_spp != "Ocean"),nn_rate~do_mg_l_rate+(1|pool_id) )
anova(NN_DO_mod)

NN_DO_mod<-lm(data = Day_rates %>%  
                  filter(foundation_spp != "Ocean",
                         before_after == "After"),nn_rate~do_mg_l_rate)
anova(NN_DO_mod)


NH_DO_mod<-lm(data = Day_rates %>%  
                filter(foundation_spp != "Ocean",
                       before_after == "After"),nh4_rate~allproddom)
anova(NH_DO_mod)

NH4_DO_mod<-lmer(data = Day_rates %>%  
                 filter(foundation_spp != "Ocean"),nh4_rate~do_mg_l_rate+(1|pool_id))
anova(NH4_DO_mod)


Prot_Hetero<-Day_rates %>%  
  filter(foundation_spp != "Ocean")%>%
  filter(prot_rate < 0.3)%>% # remove three outlier... figure out which..
 #  pivot_longer(cols = c(humics_rate, prot_rate)) %>%
  # filter(removal_control == "Control")%>%
  ggplot(aes(y = prot_rate, x = hetero_rate))+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5)+
  geom_text(aes(x = 175, y = 0.1, label = "p = 0.036"))+
  ylim(-.3,.3)+
  geom_smooth(method = "lm", color = "black")+
  labs(y = "Proteinaceaous fDOM (raman units hr-1)",
       x = "heterotrophic bacteria (counts mL-1 hr-1)",
       color = "")+
  theme_bw()+
  theme(panel.grid = element_blank())

Prot_NH4<-Day_rates %>%  
  filter(foundation_spp != "Ocean",
         before_after == "After")%>%
  filter(nh4_rate<5,prot_rate < 0.3)%>%
  #filter(prot_rate < 0.3, nh4_rate<4)%>% # remove three outlier... figure out which..
  #  pivot_longer(cols = c(humics_rate, prot_rate)) %>%
  # filter(removal_control == "Control")%>%
  ggplot(aes(x = nh4_rate, y = prot_rate))+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5)+
  ylim(-.3,.3)+
  geom_text(aes(x = 2.5, y = -0.1, label = "p = 0.027"))+
  geom_smooth(method = "lm", color = "black", lty = 2)+
  labs(y = "Proteinaceaous fDOM (raman units hr-1)",
       x = "NH4 rate (umol L-1 hr-1)",
       color = "")+
  theme_bw()+
  theme(panel.grid = element_blank())

Het_prot_mod<-lm(data = Day_rates %>% 
                   filter(nh4_rate<5,
                          prot_rate < 0.3)%>%
                   drop_na(prot_rate,nh4_rate), prot_rate~nh4_rate+hetero_rate)

anova(Het_prot_mod)



Prot_Hetero|Prot_NH4
ggsave(here("Output","Prot_het_nh4.png"), width = 8, height = 4)



## bring all the plots together
(Prot_Hetero|Prot_NH4)/(NN_prod|NH4_prod)/(P_DO_Prod|P_DO_het)

ggsave(here("Output","regressionmods.png"), height = 10, width = 8)

## production driving prod fdom
Day_rates %>%  
  filter(foundation_spp != "Ocean",
         before_after == "After")%>%
  filter(mc_rate > -0.3)%>%
  #filter(prot_rate < 0.3, nh4_rate<4)%>% # remove three outlier... figure out which..
  #  pivot_longer(cols = c(humics_rate, prot_rate)) %>%
  # filter(removal_control == "Control")%>%
      ggplot(aes(x = do_mg_l_rate, y = mc_rate))+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5)+
#  ylim(-.3,.3)+
 # geom_text(aes(x = 2.5, y = -0.1, label = "p = 0.027"))+
  geom_smooth(method = "lm", color = "black", lty = 2)+
  labs(y = "",
       x = "DO (mg L-1 hr-1)",
       color = "")+
  theme_bw()+
  theme(panel.grid = element_blank())


mc_DO_mod<-lm(data = Day_rates %>% 
                   filter(mc_rate > -0.3)%>%
                   drop_na(mc_rate), mc_rate~do_mg_l_rate)

anova(mc_DO_mod)


## Do a multivariate community analysis on the x
BD<- Day_rates %>%  
  filter(foundation_spp != "Ocean")%>%
  filter(before_after == "After")  %>%
  drop_na(nn_rate, nh4_rate, do_mg_l_rate)%>%
  select(AdjMusselCover, AdjSurfgrassCover, Diatoms, allCCA, consumers, macroalgae)

CCA_mod <- cca(BD ~ nn_rate + nh4_rate +  hetero_rate+prot_rate+humics_rate +do_mg_l_rate, data = Day_rates %>%  
                 filter(foundation_spp != "Ocean")%>%
                 filter(before_after == "After")%>%
                 drop_na(nn_rate, nh4_rate, do_mg_l_rate))
plot(CCA_mod)

summary(CCA_mod)

adonis2(BD ~ nn_rate + nh4_rate +  hetero_rate+prot_rate+humics_rate+humics_rate +do_mg_l_rate , data = Day_rates %>%  
         filter(foundation_spp != "Ocean")%>%
         filter(before_after == "After")%>%
         drop_na(nn_rate, nh4_rate, do_mg_l_rate))


# vectors
ccavectors <- as.matrix(scores(CCA_mod, display = "bp", scaling = "species")*2) %>% 
  as.data.frame() %>%
  mutate(nicenames = c("Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>hr<sup>-1</sup>)",
                       "Ammonium <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
                       "Heterotrophic Bacteria <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                       "Humic-like fDOM <br> (raman units hr<sup>-1</sup>)",
                       "Proteinaceous-like fDOM <br> (raman units hr<sup>-1</sup>)",
                       "DO <br> (mg L<sup>-1</sup> hr<sup>-1</sup>)"
                       )
  )
           

# site coordinates
site_data <- scores(CCA_mod, display = "sites") %>% 
  as.data.frame() %>% 
  bind_cols(Day_rates %>%  
              filter(foundation_spp != "Ocean")%>%
              filter(before_after == "After")%>%
              drop_na(nn_rate, nh4_rate, do_mg_l_rate))

# species coordinates
species_data <- scores(CCA_mod, display = "species") %>% 
  as.data.frame() %>%
  mutate(names = c("% Mussels","% Surfgrass", "% Diatoms"," % CCA", "% Non-Mussel Inverts", "% Fleshy Macroalgae"))

# plotting
plot_cca <- ggplot(site_data) +
  geom_point(aes(x = CCA1, y = CCA2), shape = 19, size = 2, alpha = 0.2) +
  coord_fixed() +
  geom_segment(data = ccavectors, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#  geom_text(data = species_data, aes(x = CCA1, y = CCA2, label = names),  size = 3, color = "slateblue") +
  scale_x_continuous(limits = c(-2,2)) +
#  scale_y_continuous(limits = c(-3, 12)) +
  geom_richtext(data = ccavectors, aes(x = CCA1, y = CCA2, label = nicenames), nudge_x = 0.2, nudge_y = 0.2, fill = NA, label.color = NA) +
  labs(title = "Canonical Correspondence Analysis",
       x = "CCA 1 (44.5%)",
       y = "CCA 2 (19.7%)")+
  theme_bw()+
  theme(panel.grid = element_blank())

plot_cca

plot_cca2 <- ggplot(site_data) +
  geom_point(aes(x = CCA1, y = CCA2), shape = 19, size = 2, alpha = 0.2) +
  coord_fixed() +
#  geom_segment(data = ccavectors, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_text(data = species_data, aes(x = CCA1, y = CCA2, label = names),  size = 3.5, color = "slateblue") +
    scale_x_continuous(limits = c(-2, 2)) +
  #  scale_y_continuous(limits = c(-3, 12)) +
#  geom_richtext(data = ccavectors, aes(x = CCA1, y = CCA2, label = nicenames), nudge_x = 0.2, nudge_y = 0.2, fill = NA, label.color = NA) +
  labs(title = "Canonical Correspondence Analysis",
       x = "CCA 1 (44.5%)",
       y = "CCA 2 (19.7%)")+
  theme_bw()+
  theme(panel.grid = element_blank())

plot_cca|plot_cca2

ggsave(here("Output","CCA_plot.png"), width = 12, height = 8)

# comm<-Communitymetrics %>% 
#   clean_names() %>%
#   mutate(pool_id = as.character(pool_id))%>%
#   left_join(Day_rates) %>%
#   filter(day_night == "Day",
#          before_after == "After") %>%
#   drop_na(nn_rate, nh4_rate, do_mg_l_rate)%>%
#     select(rock:stylantheca_spp) 
#   
#   
#   
# comm
# 
#   birdCCA <- cca(comm ~ nn_rate + nh4_rate +  hetero_rate+ mc_rate , data = Day_rates %>%  
#                    filter(foundation_spp != "Ocean")%>%
#                    filter(before_after == "After")%>%
#                    drop_na(nn_rate, nh4_rate, do_mg_l_rate))
# 
#   plot(birdCCA)
#   