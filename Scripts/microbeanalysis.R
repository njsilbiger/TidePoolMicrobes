### Data visualization and analysis for fDOM and microbe datasets ###
### By Nyssa Silbiger #####
### Created on 8/16/2023 #####


### Load libraries

library(here)
library(tidyverse)
library(janitor)
library(ggtext)
library(patchwork)


## read in the data ####

#biogeodata<-read_csv(here("Output","MicrobesTime4and5NECNEP.csv"))
microbedata<-read_csv(here("Data","Microbe_Clean","combined_FCMandfDOMdata.csv"))

biogeodata<-read_csv(here("Data","Biogeochem","MicrobeCarbChem.csv"))


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
  filter(foundation_spp !="Ocean", before_after !="Before", time_point !=5, removal_control == "Control") %>%
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
          .by = c(foundation_spp, removal_control, day_night, pool_id))%>%
  drop_na(foundation_spp) %>%
  mutate(before_after = "Before")

data_rates_after<-data_all%>%
  ungroup()%>%
  mutate(sampling_datetime = mdy_hms(paste(sampling_day, sampling_time))) %>% 
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
          .by = c(foundation_spp, removal_control, day_night, pool_id, sampling_group))%>%
  drop_na(foundation_spp)%>%
  mutate(before_after = "After")

# bring before and after together
data_rates <-
  bind_rows(data_rates_before, data_rates_after) %>%
  mutate(nh4_rate = ifelse(nh4_rate>10, NA, nh4_rate)) # remove crazy Nh4 outlier

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
  pivot_longer(cols= do_mg_l_rate:auto_rate) %>%
 # filter(!pool_id %in% 1:16)%>% ## uneven sample sizes... only 16 pools have everything
  group_by(foundation_spp, removal_control, day_night, before_after, name)%>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE)/sqrt(n()))

# rename the foundation species for the ocean to have one group for mussels and one for phyllo for easier plotting
myt_ocean<-mean_rates %>%
  filter(foundation_spp == "Ocean") %>%
  mutate(foundation_spp = "Mytilus")

phylo_ocean<-mean_rates %>%
  filter(foundation_spp == "Ocean") %>%
  mutate(foundation_spp = "Phyllospadix")

# make a set of reaction norm plots for the day
day_plot<-mean_rates %>%
  filter(foundation_spp != "Ocean")%>%
  bind_rows(myt_ocean) %>%
  bind_rows(phylo_ocean) %>%
  filter(day_night == "Day")%>%
  ggplot(aes(x = before_after, y = mean_value, color = removal_control, group = removal_control))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point()+
  geom_errorbar(aes(x = before_after, y = mean_value, ymin = mean_value-se_value, ymax = mean_value+se_value), width = 0.01)+
  geom_line()+
  labs(x = "time period",
       y = "Rate per hour")+
  facet_wrap(name~foundation_spp, scales = "free", ncol = 2)+
  theme_bw()

day_plot

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
mean_rates %>%
  filter(removal_control == "Control", before_after == "After")
