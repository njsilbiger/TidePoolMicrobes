### Using the cleaned dataset with only start and end data ###


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
library(ggh4x)


### read in data #########
BenthicData<-read_csv(here("Data","Microbe_Clean","CommunityData.csv"))

MetaData<-read_csv(here("Data","Microbe_Clean","TidePoolDescriptions.csv"))

data_all<-read_csv(here("Data","Microbe_Clean","joinedData_edited.csv"))


### show difference from the ocean
data_all<-data_all %>%
  mutate(nh4_umol_l = ifelse(nh4_umol_l>200, NA, nh4_umol_l),
         synechoococcus_m_l = ifelse(synechoococcus_m_l>50, NA, synechoococcus_m_l)) %>% # drop the 2 crazy outliers
  mutate(time_point = factor(time_point, levels = c("start","end")))%>%
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


## add in the metadata
data_all<-data_all %>%
  left_join(MetaData %>%
              clean_names() %>%
              mutate(pool_id = as.character(pool_id)))%>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)"))) 


## make some plots
data_all  %>%
  filter(foundation_spp !="Ocean", before_after !="Before") %>%
  dplyr::select(pool_id, foundation_spp, removal_control,day_night, time_point, do_change:phenyl_change)%>%
  mutate(nh4_change = ifelse(nh4_change>200, NA, nh4_change))%>% # drop the 2 crazy outliers
  pivot_longer(cols = do_change:phenyl_change, names_to = "parameter", values_to = "values") %>%
  mutate(found_time = paste(foundation_spp,time_point))%>%
  ggplot(aes(x = found_time, y = values, fill = day_night))+
  geom_boxplot()+
  labs(y = "difference from ocean sample")+
  facet_wrap(removal_control~parameter, scales = "free")

## make a pca of all the data
pcadata<-data_all  %>%
  filter(foundation_spp !="Ocean", 
         #before_after !="Before", 
         #time_point !="end", 
         #removal_control == "Control"
         ) %>%
  ungroup()%>%
  dplyr::select(do_change:phenyl_change) %>%
  drop_na()

# run a pca
pca<-prcomp(pcadata, scale. = TRUE, center = TRUE)

# calculate percent explained by each PC
perc.explained<-round(100*pca$sdev/sum(pca$sdev),1)

# Extract the scores and loadings
PC_scores <-as_tibble(pca$x[,1:2])

PC_loadings<-as_tibble(pca$rotation)%>%
  bind_cols(labels = rownames(pca$rotation))


data_pca<-data_all  %>%
  filter(foundation_spp !="Ocean", 
         #before_after !="Before", 
         #time_point !="end", 
         #removal_control == "Control"
  ) %>%
  ungroup()%>%
  drop_na()%>%
  bind_cols(PC_scores)


p1<-data_pca %>%
  #mutate(time_point = factor(time_point))%>%
  ggplot(aes(x = PC1, y = PC2, color = removal_control, shape = time_point))+
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  # scale_shape_manual(values = c(1, 22,15,16))+
  scale_colour_manual(values = c("#D64550","#EA9E8D"))+
  scale_fill_manual(values = c("#D64550","#EA9E8D"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(fill = removal_control, label = paste(time_point, removal_control), color =removal_control),
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
  facet_wrap(foundation_spp~day_night)


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

## plot the raw ocean data of before and after
data_all %>%
  ungroup()%>%
  filter(foundation_spp == "Ocean") %>%
 # filter(day_night == "Day")  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like) %>%
  select(month, day_night, do_mg_l,heterotrophic_bac_m_l=heterotrophic_bacterioplankton_m_l,prot,temperature = temp_pool, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  group_by(month, name, day_night)%>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val= sd(value, na.rm = TRUE)/sqrt(n())) %>%
  ggplot(aes(x = month, y = mean_val, color = day_night))+
  geom_point(size = 3)+
  geom_errorbar(aes(x = month, ymin = mean_val - se_val, ymax = mean_val+se_val), width = 0.1)+
  facet_wrap(~name, scales = "free_y", strip.position = "left", ncol = 1)+
  labs(x = "",
       y = "",
       title = "Mean Ocean values")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside")


data_long<-data_all %>%
  ungroup()%>%
  filter(removal_control != "Removal") %>%
  filter(day_night == "Day")  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like) %>%
  select(month, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot, humic, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l)

data_all %>%
  ungroup()%>%
  filter(removal_control != "Removal") %>%
   filter(day_night == "Day")  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like) %>%
  select(month, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot, humic, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  group_by(month, name, foundation_spp,time_point, day_night, removal_control)%>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val= sd(value, na.rm = TRUE)/sqrt(n())) %>%
  ggplot(aes(x = time_point, y = mean_val, color = foundation_spp, group = foundation_spp)
             #group = interaction(foundation_spp, removal_control))
         )+
  geom_point(size = 3)+
 # geom_point(data = data_long, aes(x = time_point, y = value))+
  geom_line()+
  geom_errorbar(aes(x = time_point, ymin = mean_val - se_val, ymax = mean_val+se_val), width = 0.1)+
  facet_wrap(name~month, scales = "free_y", strip.position = "left", ncol = 2)+
  labs(x = "",
       y = "",
       title = "Mean values")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside")

### create a function that calculated difference between start and end values for each tide pool

change_val<-function(x, time){
  val<-x[time == "end"]-x[time == "start"]
  val<-as.numeric(val)
  return(val)
}


# first get the times 
## calculate the change in time for each tidepool to join
times<-data_all%>%
  ungroup()%>%
  mutate(sampling_datetime = mdy_hms(paste(sampling_day, sampling_time))) %>% 
  select(pool_id,before_after,day_night, time_point, sampling_datetime,sampling_group) %>%
  group_by(pool_id, before_after, day_night,sampling_group) %>%
  reframe(diff_time = change_val(sampling_datetime, time_point)) 

# now calculate all the rates per hour
Rates<-data_all%>%
  ungroup()%>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like) %>%
  #mutate(sampling_datetime = as.numeric(mdy_hms(paste(sampling_day, sampling_time)))) %>% 
  select(pool_id:removal_control, time_point, do_mg_l, po_umol_l:nh4_umol_l, heterotrophic_bacterioplankton_m_l:fi, prot, humic) %>%
  # mutate(heterotrophic_bacterioplankton_m_l = heterotrophic_bacterioplankton_m_l *1000, # convert to per L
  #        autotrophic_pico_eukaryotes_m_l =  autotrophic_pico_eukaryotes_m_l*1000,
  #        synechoococcus_m_l = synechoococcus_m_l*1000
  #        )%>%
  pivot_longer(cols = do_mg_l:humic)%>%
  group_by(pool_id, before_after, removal_control,day_night, foundation_spp, sampling_group, name) %>%
  reframe(change = change_val(value, time_point)) %>% # calculate the difference between start and end
  left_join(times)%>%# join with the times
  left_join(MetaData %>%
              clean_names %>%
              mutate(pool_id = as.character(pool_id)) %>%
              select(pool_id, before_after, surface_area,vol))%>% # add in the tide pool info
  mutate(rate_hr = change/diff_time,# difference in value per hour
         rate_m2_hr = rate_hr*vol/surface_area
         ) %>%
  mutate(rate_hr = ifelse(name == "po_umol_l" & rate_hr < -2 & day_night == "Day", NA, rate_hr))%>%
  mutate(rate_hr = ifelse(name == "synechoococcus_m_l" & rate_hr > 0.6 & day_night == "Day", NA, rate_hr))%>%
  mutate(rate_hr = ifelse(name == "autotrophic_pico_eukaryotes_m_l" & rate_hr < -0.6 & day_night == "Day", NA, rate_hr))
  


# make a bunch of boxplots
Rates %>%
  filter(day_night == "Day")%>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  ggplot(aes(x = before_after, y = rate_hr, fill = removal_control))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_boxplot()+
  facet_wrap(name~foundation_spp, scales = "free")
  
## Make some reaction norms
mean_rates<-Rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  group_by(foundation_spp, removal_control, day_night, before_after, name)%>%
  summarise(mean_value = mean(rate_hr, na.rm = TRUE),
            se_value = sd(rate_hr, na.rm = TRUE)/sqrt(n()),
            mean_value_m2 = mean(rate_m2_hr, na.rm = TRUE),
            se_value_m2 = sd(rate_m2_hr, na.rm = TRUE)/sqrt(n()))


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
    name == "autotrophic_pico_eukaryotes_m_l" ~"Autotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
    name == "do_mg_l" ~ "DO <br> (mg L<sup>-1</sup> hr<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
    name == "synechoococcus_m_l" ~ "Synechoococcus <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "Ammonium <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
    name == "po_umol_l" ~ "Phosphate <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
    name == "nn_umol_l" ~ "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>hr<sup>-1</sup>)",
    name == "humic" ~ "Humic_like <br> (Raman units hr<sup>-1</sup>)",
name == "prot" ~ "Proteinaceous <br> (Raman units hr<sup>-1</sup>)")
)

# make a set of reaction norm plots for the day
day_plot<-mean_plot %>%
  filter(day_night == "Day",
         !name %in% c("fi","hix","m_c","bix") )%>%
  # filter(name %in% c(
  #   "hetero_rate","do_mg_l_rate",
  #   "nh4_rate","nn_rate"))%>%
  mutate(nicenames = factor(nicenames, levels = c("DO <br> (mg L<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Ammonium <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>hr<sup>-1</sup>)",
                                                  "Phosphate <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Autotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Synechoococcus <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Heterotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Humic_like <br> (Raman units hr<sup>-1</sup>)",
                                                  "Proteinaceous <br> (Raman units hr<sup>-1</sup>)"
                                                  
    
  )))%>%
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
      scale_y_continuous(limits=c(-0.5, 4)),
      scale_y_continuous(limits=c(-1, 5)),
      scale_y_continuous(limits=c(-1.5, 2.5)),
      scale_y_continuous(limits = c(-0.2, 0.4)),
      scale_y_continuous(limits=c(-0.25, 0.25)),
      scale_y_continuous(limits = c(-0.3, 0.2)),
      scale_y_continuous(limits = c(-75, 150)),
      scale_y_continuous(limits=c(-0.01, 0.08)),
      scale_y_continuous(limits=c(-0.1, 0.2))
      
    ), each = 2))+
  theme_bw()+
  theme(strip.text.y.left  = ggtext::element_markdown(size = 16),
        strip.text.x = element_blank(),
        strip.placement = "outside",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "left")

day_plot
ggsave(plot = day_plot, filename = here("Output","AllRates.pdf"), width = 8, height = 18)

# make a set of reaction norm plots for the night
night_plot<-mean_plot %>%
  filter(day_night == "Night",
         !name %in% c("fi","hix","m_c","bix") )%>%
  # filter(name %in% c(
  #   "hetero_rate","do_mg_l_rate",
  #   "nh4_rate","nn_rate"))%>%
  mutate(nicenames = factor(nicenames, levels = c("DO <br> (mg L<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Ammonium <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>hr<sup>-1</sup>)",
                                                  "Phosphate <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Autotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Synechoococcus <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Heterotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Humic_like <br> (Raman units hr<sup>-1</sup>)",
                                                  "Proteinaceous <br> (Raman units hr<sup>-1</sup>)"
                                                  
                                                  
  )))%>%
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
      scale_y_continuous(limits=c(-2, 0.5)),
      scale_y_continuous(limits=c(-2, 3)),
      scale_y_continuous(limits=c(-1, 0.5)),
      scale_y_continuous(limits = c(-0.4, 0.8)),
      scale_y_continuous(limits=c(-0.5, 0.25)),
      scale_y_continuous(limits = c(-2, 1)),
      scale_y_continuous(limits = c(-100, 150)),
      scale_y_continuous(limits=c(-0.1, 0.1)),
      scale_y_continuous(limits=c(-0.2, 0.15))
      
    ), each = 2))+
  theme_bw()+
  theme(strip.text.y.left  = ggtext::element_markdown(size = 16),
        strip.text.x = element_blank(),
        strip.placement = "outside",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "left")

night_plot
ggsave(plot = night_plot, filename = here("Output","AllRates_night.pdf"), width = 8, height = 18)

## Make a plot of the benthic data

BenthicData %>%
  select(Foundation_spp, Removal_Control, PoolID, Before_After, macroalgae, Diatoms, consumers, allCCA, AdjMusselCover, AdjSurfgrassCover)%>%
  mutate(rock_sand = (100-(macroalgae+consumers+allCCA+AdjMusselCover+AdjSurfgrassCover)),
         macroalgae = macroalgae - Diatoms,
         Before_After = factor(Before_After, levels = c("Before","After"))) %>% # the macroalgae includes diatom cover and I want to see the difference
  pivot_longer(cols = macroalgae:rock_sand) %>%
  ggplot(aes(x = PoolID, y = value, fill = name))+
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cal_palette("tidepool", n = 7, type = "continuous"))+
  facet_wrap(Before_After~Removal_Control, scale = "free_x")+
  theme_classic()

BenthicData %>%
 # filter(Before_After == "After")%>%
  mutate(Before_After = factor(Before_After, levels = c("Before","After")))%>%
  select(PoolID, Before_After,  Macroalgae = macroalgae , Diatoms, `Non-mussel Consumers` = consumers, CCA = allCCA, Mussels = AdjMusselCover, Surfgrass = AdjSurfgrassCover)%>%
  mutate(PoolID =  as.factor(PoolID)) %>%
  mutate(`Rock/Sand` = (100-(Macroalgae+`Non-mussel Consumers`+CCA+Mussels+Surfgrass)),
         Macroalgae = Macroalgae - Diatoms,
         PoolID = fct_reorder2(PoolID,`Rock/Sand`, Mussels),) %>% # the macroalgae includes diatom cover and I want to see the difference
  pivot_longer(cols = Macroalgae:`Rock/Sand`) %>%
  mutate(name = factor(name, levels = c("Surfgrass","Macroalgae","Mussels","Non-mussel Consumers", "Diatoms","CCA","Rock/Sand")))%>%
  ggplot(aes(x = PoolID, y = value, fill = name))+
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cal_palette("tidepool", n = 7, type = "continuous"))+
  labs(x = "",
       y = "% Cover",
       fill = "")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        axis.text.x = element_blank())+
  facet_wrap(~Before_After, ncol = 1)

ggsave(here("Output","BenthicCover.png"), width = 8, height = 8)

#### Take a regression approach #############

## combine with community comp data

Day_rates<-BenthicData  %>%
  rename(pool_id = PoolID, removal_control = Removal_Control, before_after = Before_After, foundation_spp = Foundation_spp) %>%
  mutate(pool_id = as.character(pool_id))%>%
  select(pool_id, before_after, AdjMusselCover:allproddom) %>%
  right_join(Rates, by = c("pool_id", "before_after")) %>%
  filter(day_night == "Day") %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)")))

Night_rates<-BenthicData  %>%
  rename(pool_id = PoolID, removal_control = Removal_Control, before_after = Before_After, foundation_spp = Foundation_spp) %>%
  mutate(pool_id = as.character(pool_id))%>%
  select(pool_id, before_after, AdjMusselCover:allproddom) %>%
  right_join(Rates, by = c("pool_id", "before_after")) %>%
  filter(day_night == "Night") %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)")))

## make some plots
P_Nuts<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  filter(!name %in% c("fi","bix","hix","m_c") )%>%
  mutate(nicenames = case_when(
    name == "autotrophic_pico_eukaryotes_m_l" ~"Autotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
    name == "do_mg_l" ~ "DO <br> (mg L<sup>-1</sup> hr<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
    name == "synechoococcus_m_l" ~ "Synechoococcus <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "Ammonium <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
    name == "po_umol_l" ~ "Phosphate <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
    name == "nn_umol_l" ~ "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>hr<sup>-1</sup>)",
    name == "humic" ~ "Humic-like <br> (Raman units hr<sup>-1</sup>)",
    name == "prot" ~ "Proteinaceous <br> (Raman units hr<sup>-1</sup>)") )%>%
  mutate(nicenames = factor(nicenames, levels = c("DO <br> (mg L<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Ammonium <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>hr<sup>-1</sup>)",
                                                  "Phosphate <br> (&mu;mol L<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Humic-like <br> (Raman units hr<sup>-1</sup>)",
                                                  "Proteinaceous <br> (Raman units hr<sup>-1</sup>)",
                                                  "Autotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Synechoococcus <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
                                                  "Heterotrophic <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)"
    
  )))%>%
  ggplot(aes(x = allproddom, y = rate_hr))+
  geom_vline(aes(xintercept  = 0))+
  geom_hline(aes(yintercept  = 0))+
  geom_point(aes(color =  month), alpha = 0.5)+
  geom_smooth(method = "lm",
              formula = 'y~poly(x,2)', color = "black")+
  labs(x = "Producer-dominance (%)",
       y = "",
       color = "")+
  scale_color_manual(values = cal_palette(name = "chaparral1"))+
  theme_bw()+
  facet_wrap(~nicenames, scales = "free_y", strip.position = "left")+
  theme(strip.text.y.left = element_markdown(size = 16),
        strip.placement = "outside",
         strip.background = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.position = "top",
        legend.direction = "horizontal")

ggsave(here("Output","Regressionplots.png"), width = 12, height = 8)

# community and DO
P_DO_Prod<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  filter(month == "August (Upwelling)",
         name == "do_mg_l")%>%
  ggplot(aes(x = allproddom, y = rate_hr))+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_hline(aes(yintercept  = log(1)), lty = 2)+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", color = "black", formula = "y~poly(x,2)")+
  #  coord_trans(y = function(x)log(x+1))+
#  geom_text(aes(x = 50, y = 6, label = "p <0.001"))+
#  geom_text(aes(x = 50, y = 5.5, label = "R2 = 0.31"))+
  geom_text(aes(x = 50, y = -0.5, label = "Producer dominated"))+
  geom_text(aes(x = -50, y = -0.5, label = "Consumer dominated"))+
  labs(x = "Producer-dominance (% Producers - % Consumers)",
       y = "Productivity (DO mg L-1 hr-1)",
       color = "")+
  theme_bw()+
  theme(panel.grid = element_blank())

summary(lm(data = Day_rates %>%  
           filter(foundation_spp != "Ocean",
                  name == "do_mg_l") %>%
           filter(month == "August (Upwelling)"),rate_hr~ allproddom))

# heterotrophic bacteria and DO

P_DO_het<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  select(-c(change:vol, rate_m2_hr))%>%
  pivot_wider(values_from = rate_hr, names_from = name)%>%
  filter(month == "August (Upwelling)")%>%
  ggplot(aes(x = heterotrophic_bacterioplankton_m_l, y = do_mg_l))+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", color = "black", lty = 2)+
  geom_text(aes(x = 50, y = 6, label = "p = 0.19"))+
  labs(x = "Heterotrophic Bacteria production (counts mL-1 hr-1)",
       y = "",
       color = "")+
  theme_bw()+
  theme(panel.grid = element_blank())

P_DO_Prod+P_DO_het
ggsave(here("Output","Community_DO.png"), width = 8, height = 4)


## create a wide dataset
Day_rates_wide<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  select(-c(change:vol, rate_m2_hr))%>%
  pivot_wider(values_from = rate_hr, names_from = name)

Prod_DO_mod<-lm(data = Day_rates_wide %>%  
                  filter(before_after == "After"),do_mg_l~heterotrophic_bacterioplankton_m_l+allproddom )
anova(Prod_DO_mod)
summary(Prod_DO_mod)

# mod of producer dominance and nuts
NN_prod_mod<-lm(data = Day_rates_wide %>%  
                  filter(before_after == "After"),nn_umol_l~poly(allproddom,2))
anova(NN_prod_mod)
summary(NN_prod_mod)

po_prod_mod<-lm(data = Day_rates_wide %>%  
                  filter(before_after == "After"),po_umol_l~poly(allproddom,2))
anova(po_prod_mod)
summary(po_prod_mod)

nh4_prod_mod<-lm(data = Day_rates_wide %>%  
                  filter(before_after == "After"),nh4_umol_l~poly(allproddom,2))
anova(nh4_prod_mod)
summary(nh4_prod_mod)

### prod dom ~ bacterial counts
Het_prod_mod<-lm(data = Day_rates_wide %>%  
                  filter(before_after == "After"),heterotrophic_bacterioplankton_m_l~poly(allproddom,2))
anova(Het_prod_mod)
summary(Het_prod_mod)

Syn_prod_mod<-lm(data = Day_rates_wide %>%  
                   filter(before_after == "After"),synechoococcus_m_l~allproddom)
anova(Syn_prod_mod)
summary(Syn_prod_mod)

Auto_prod_mod<-lm(data = Day_rates_wide %>%  
                   filter(before_after == "After"),autotrophic_pico_eukaryotes_m_l~allproddom)
anova(Auto_prod_mod)
summary(Auto_prod_mod)

### plot do ~ nh4, nn, po, Hback, Syn, Auto
Day_rates_wide %>%
  select(do_mg_l, nh4_umol_l, nn_umol_l, po_umol_l, heterotrophic_bacterioplankton_m_l, synechoococcus_m_l, autotrophic_pico_eukaryotes_m_l, month) %>%
  pivot_longer(cols = heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l)%>%
  ggplot(aes(x = nh4_umol_l, y = value, color = month))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~name, scales = "free")
  
### fDom and bacteria
Prot_Hetero<-Day_rates_wide %>%  
  filter(foundation_spp != "Ocean")%>%
  ggplot(aes(y = prot, x = heterotrophic_bacterioplankton_m_l, color = month))+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5)+
  geom_text(aes(x = 175, y = 0.1, label = "p = 0.036"))+
  ylim(-.3,.3)+
  geom_smooth(method = "lm")+
  labs(y = "Proteinaceaous fDOM (raman units hr-1)",
       x = "heterotrophic bacteria (counts mL-1 hr-1)")+
  theme_bw()+
  theme(panel.grid = element_blank())

Prot_NH4<-Day_rates_wide %>%  
  filter(foundation_spp != "Ocean")%>%
  ggplot(aes(y = prot, x = nh4_umol_l, color = month))+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5)+
 # geom_text(aes(x = 175, y = 0.1, label = "p = 0.036"))+
  ylim(-.3,.3)+
  geom_smooth(method = "lm")+
  labs(y = "Proteinaceaous fDOM (raman units hr-1)",
       x = "NH4 rate (umol L-1 hr-1)")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  facet_wrap(~month, scale = "free")

Het_prot_mod<-lm(data = Day_rates_wide %>%
                   filter(month != "July"), prot~nh4_umol_l+heterotrophic_bacterioplankton_m_l)

anova(Het_prot_mod)

### DO and fDOM

DO_Humic<-Day_rates_wide %>%  
  filter(foundation_spp != "Ocean")%>%
  ggplot(aes(y = humic, x = do_mg_l, color = month))+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5)+
  # geom_text(aes(x = 175, y = 0.1, label = "p = 0.036"))+
  ylim(-.3,.3)+
  geom_smooth(method = "lm", color = "black")+
  labs(y = "Humic fDOM (raman units hr-1)",
       x = "DO rate (mg L-1 hr-1)")+
  theme_bw()+
  theme(panel.grid = element_blank())
#  facet_wrap(~month, scale = "free")

DO_humic_mod<-lm(data = Day_rates_wide %>%
                   filter(month == "July"), humic~do_mg_l)

anova(DO_humic_mod)
