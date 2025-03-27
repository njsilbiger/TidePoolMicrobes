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
library(brms)


### read in data #########
BenthicData<-read_csv(here("Data","Microbe_Clean","CommunityData.csv"))

MetaData<-read_csv(here("Data","Microbe_Clean","TidePoolDescriptions.csv"))

data_all<-read_csv(here("Data","Microbe_Clean","joinedData_edited.csv"))


### show difference from the ocean
data_all<-data_all %>%
  mutate(nh4_umol_l = ifelse(nh4_umol_l>60, NA, nh4_umol_l),
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

# total fDOM in the ocean
data_all %>%
  filter(removal_control == "Ocean",
         day_night == "Day") %>%
  select(sampling_group, before_after, time_point,ultra_violet_humic_like:phenylalanine_like) %>%
  mutate(totalfDOM = rowSums(across(ultra_violet_humic_like:phenylalanine_like))) %>%
  group_by(before_after)%>%
  summarise(mean_fdom = mean(totalfDOM),
            se_fDOM = sd(totalfDOM, na.rm = TRUE))

### ADD TOTAL FDOM TO ANALYSIS AND SEE IF THERE IS AN EFFECT

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
  filter(foundation_spp == "Ocean",
         day_night == "Day") %>%
 # filter(day_night == "Day")  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = prot+humic) %>%
  select(month, day_night, do_mg_l,heterotrophic_bac_m_l=heterotrophic_bacterioplankton_m_l,total_fDOM,temperature = temp_pool, nn_umol_l, nh4_umol_l) %>%
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
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = prot+humic) %>%
  select(month, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot, humic,total_fDOM, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l)

### MAKE THIS PLOT PRETTY #######

# All pools raw data, then control in Aug and removal in Aug
data_long_day<-data_all %>%
  ungroup()%>%
  filter(day_night == "Day")  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = prot+humic) %>%
  select(month, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l, nn_umol_l, nh4_umol_l, bix, hix, m_c, total_fDOM) %>%
  pivot_longer(cols = do_mg_l:total_fDOM) %>%
  group_by(month, name, foundation_spp,time_point, day_night, removal_control)%>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val= sd(value, na.rm = TRUE)/sqrt(n())) %>%
  mutate(time_point_clean = ifelse(time_point == "start","Early", "Late"))%>%
  mutate(nicenames = case_when(
   # name == "autotrophic_pico_eukaryotes_m_l" ~"Autotrophic <br> (# mL<sup>-1</sup>)",
    name == "do_mg_l" ~ "DO <br> (mg L<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic <br> (# mL<sup>-1</sup>)",
  #  name == "synechoococcus_m_l" ~ "Synechoococcus <br> (# mL<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
    name == "po_umol_l" ~ "Phosphate <br> (&mu;mol L<sup>-1</sup>)",
    name == "nn_umol_l" ~ "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
  #  name == "humic" ~ "Humic-like <br> (Raman units)",
  name == "total_fDOM" ~"Total fDOM <br> (Raman units)",
  #  name == "prot" ~ "Proteinaceous <br> (Raman units)",
    name == "hix"~"HIX",
    name == "bix"~"BIX",
    name == "m_c"~"M:C")
  ) %>%
  mutate(nicenames = factor(nicenames, levels = c("DO <br> (mg L<sup>-1</sup>)",
                                                  "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Phosphate <br> (&mu;mol L<sup>-1</sup>)",
                                                #  "Autotrophic <br> (# mL<sup>-1</sup>)",
                                                #  "Synechoococcus <br> (# mL<sup>-1</sup>)",
                                                  "Heterotrophic <br> (# mL<sup>-1</sup>)",
                                                  # "Humic-like <br> (Raman units)",
                                                  # "Proteinaceous <br> (Raman units)",
                                                  "HIX","BIX","M:C"
                                                  
                                                  
  )))

#### FIX STARTING HERE

P_july<-data_all %>%
  ungroup()%>%
  filter(day_night == "Day")  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = prot+humic) %>%
  select(month, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l, nn_umol_l, nh4_umol_l, bix, hix, m_c, total_fDOM) %>%
  pivot_longer(cols = do_mg_l:total_fDOM) %>%
  group_by(month, name, foundation_spp,time_point, day_night)%>% ### all the pools together
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val= sd(value, na.rm = TRUE)/sqrt(n())) %>%
  mutate(time_point_clean = ifelse(time_point == "start","Early", "Late"))%>%
  mutate(nicenames = case_when(
  #  name == "autotrophic_pico_eukaryotes_m_l" ~"Autotrophic <br> (# mL<sup>-1</sup>)",
    name == "do_mg_l" ~ "DO <br> (mg L<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic <br> (# mL<sup>-1</sup>)",
 #   name == "synechoococcus_m_l" ~ "Synechoococcus <br> (# mL<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
    name == "po_umol_l" ~ "Phosphate <br> (&mu;mol L<sup>-1</sup>)",
    name == "nn_umol_l" ~ "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
  #  name == "humic" ~ "Humic-like <br> (Raman units)",
  #  name == "prot" ~ "Proteinaceous <br> (Raman units)",
    name == "total_fDOM" ~"Total fDOM <br> (Raman units)",
    name == "hix"~"HIX",
    name == "bix"~"BIX",
    name == "m_c"~"M:C")
  ) %>%
  mutate(nicenames = factor(nicenames, levels = c("DO <br> (mg L<sup>-1</sup>)",
                                                  "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Phosphate <br> (&mu;mol L<sup>-1</sup>)",
                                            #      "Autotrophic <br> (# mL<sup>-1</sup>)",
                                            #      "Synechoococcus <br> (# mL<sup>-1</sup>)",
                                                  "Heterotrophic <br> (# mL<sup>-1</sup>)",
                                                 # "Humic-like <br> (Raman units)",
                                                #  "Proteinaceous <br> (Raman units)",
                                                  "HIX","BIX","M:C", "Total fDOM <br> (Raman units)"
                                                  
                                                  
  )))%>%
  filter(month == "July",
#         foundation_spp != "Ocean"
  ) %>%
  ggplot(aes(x = time_point_clean, y = mean_val, color = foundation_spp, group = foundation_spp)
             #group = interaction(foundation_spp, removal_control))
         )+
  geom_point(size = 3)+
 # geom_point(data = data_long, aes(x = time_point, y = value))+
  geom_line()+
  geom_errorbar(aes(x = time_point_clean, ymin = mean_val - se_val, ymax = mean_val+se_val), width = 0.1)+
  facet_wrap(~nicenames, scales = "free_y", strip.position = "left", ncol = 1)+
  facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits=c(6, 25)),
      scale_y_continuous(limits=c(0, 30)),
      scale_y_continuous(limits=c(0, 16)),
         #   scale_y_continuous(limits = c(0, 3)),
   #   scale_y_continuous(limits=c(0, 6)),
      scale_y_continuous(limits = c(0, 1000)),
      scale_y_continuous(limits = c(0.5, 1.5)),
      scale_y_continuous(limits=c(0.8, 1.2)),
      scale_y_continuous(limits=c(0.9, 1.35)),
      scale_y_continuous(limits=c(0, 2))

    ), each = 1))+
  labs(x = "",
       y = "",
       color = "",
       title = "Before \n (July)")+
  scale_color_manual(values = c("grey30",
                                "#79ACBD",
                                "#567d46"))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_markdown(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 14)
          )

# August control
P_Aug_control<-data_long_day %>%
  filter(month != "July",
         removal_control != "Removal",
#         foundation_spp != "Ocean"
  ) %>%
  ggplot(aes(x = time_point_clean, y = mean_val, color = foundation_spp, group = foundation_spp)
         #group = interaction(foundation_spp, removal_control))
  )+
  geom_point(size = 3)+
  # geom_point(data = data_long, aes(x = time_point, y = value))+
  geom_line()+
  geom_errorbar(aes(x = time_point_clean, ymin = mean_val - se_val, ymax = mean_val+se_val), width = 0.1)+
  facet_wrap(~nicenames, scales = "free_y", strip.position = "left", ncol = 1)+
  facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits=c(6, 25)),
      scale_y_continuous(limits=c(0, 30)),
      scale_y_continuous(limits=c(0, 16)),
  #    scale_y_continuous(limits = c(0, 3)),
  #    scale_y_continuous(limits=c(0, 6)),
      scale_y_continuous(limits = c(0, 1000)),
  scale_y_continuous(limits = c(0.5, 1.5)),
  scale_y_continuous(limits=c(0.8, 1.2)),
  scale_y_continuous(limits=c(0.9, 1.35)),
  scale_y_continuous(limits=c(0, 2))
    ), each = 1))+
  labs(x = "",
       y = "",
       color = "",
       title = "After Control \n (August upwelling)")+
  scale_color_manual(values = c("grey30",
                                "#79ACBD",
                                "#567d46"))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 14)
  )

# August impact

P_Aug_impact<-data_long_day %>%
  filter(month != "July",
         removal_control != "Control",
 #        foundation_spp != "Ocean"
  ) %>%
  ggplot(aes(x = time_point_clean, y = mean_val, color = foundation_spp, group = foundation_spp)
         #group = interaction(foundation_spp, removal_control))
  )+
  geom_point(size = 3)+
  # geom_point(data = data_long, aes(x = time_point, y = value))+
  geom_line()+
  geom_errorbar(aes(x = time_point_clean, ymin = mean_val - se_val, ymax = mean_val+se_val), width = 0.1)+
  facet_wrap(~nicenames, scales = "free_y", strip.position = "left", ncol = 1)+
  facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits=c(6, 25)),
      scale_y_continuous(limits=c(0, 30)),
      scale_y_continuous(limits=c(0, 16)),
   #   scale_y_continuous(limits = c(0, 3)),
  #    scale_y_continuous(limits=c(0, 6)),
      scale_y_continuous(limits = c(0, 1000)),
  scale_y_continuous(limits = c(0.5, 1.5)),
  scale_y_continuous(limits=c(0.8, 1.2)),
  scale_y_continuous(limits=c(0.9, 1.35)),
  scale_y_continuous(limits=c(0, 2))
    ), each = 1))+
  labs(x = "",
       y = "",
       color = "",
       title = "After Impact \n (August upwelling)")+
  scale_color_manual(values = c("grey30",
                                "#79ACBD",
                                "#567d46"))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 14)
  )


P_july+ theme(legend.position = "none")|P_Aug_control+ theme(legend.position = "none")+plot_layout(guides = "collect")|P_Aug_impact&theme(legend.position = "bottom")
ggsave(here("Output","mean_chem.pdf"), width = 10, height = 16, device = cairo_pdf)

####### Look at it the BACI way #####

data_all %>%
  ungroup()%>%
  filter(day_night == "Day",
         foundation_spp != "Ocean")  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = prot+humic) %>%
  select(month,pool_id, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,total_fDOM, humic, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  group_by(foundation_spp, pool_id,removal_control,name, time_point) %>%
  reframe(value_diff = value[month == "August (Upwelling)"] - value[month == "July"]) %>%
  ungroup()%>%
  group_by(foundation_spp, removal_control, name, time_point)%>%
  summarise(value_diff_mean = mean(value_diff, na.rm = TRUE),
            value_diff_se = sd(value_diff, na.rm = TRUE)/sqrt(n()),
  )%>%
  ggplot(aes(x = removal_control, y = value_diff_mean, color = foundation_spp, group = foundation_spp))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin = value_diff_mean-value_diff_se,ymax = value_diff_mean+value_diff_se), width = 0.1)+
  facet_wrap(time_point~name, scales = "free_y")


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
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = prot+humic) %>%
  #mutate(sampling_datetime = as.numeric(mdy_hms(paste(sampling_day, sampling_time)))) %>% 
  select(pool_id:removal_control, time_point, do_mg_l, po_umol_l:nh4_umol_l, heterotrophic_bacterioplankton_m_l:fi, prot, humic,total_fDOM) %>%
   mutate(heterotrophic_bacterioplankton_m_l = heterotrophic_bacterioplankton_m_l *1000, # convert to per L
          autotrophic_pico_eukaryotes_m_l =  autotrophic_pico_eukaryotes_m_l*1000,
          synechoococcus_m_l = synechoococcus_m_l*1000
          )%>%
  pivot_longer(cols = do_mg_l:total_fDOM)%>%
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
  

## run the BACI design
Rates %>%
  filter(day_night == "Day",
         foundation_spp != "Ocean")%>%
  group_by(foundation_spp, pool_id,removal_control,name) %>%
  reframe(Rate_diff = rate_hr[before_after == "After"] - rate_hr[before_after == "Before"]) %>%
  ungroup()%>%
  group_by(foundation_spp, removal_control, name)%>%
  summarise(rate_diff_mean = mean(Rate_diff, na.rm = TRUE),
            rate_diff_se = sd(Rate_diff, na.rm = TRUE)/sqrt(n()),
  )%>%
  ggplot(aes(x = removal_control, y = rate_diff_mean, color = foundation_spp, group = foundation_spp))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin = rate_diff_mean-rate_diff_se,ymax = rate_diff_mean+rate_diff_se), width = 0.1)+
  facet_wrap(~name, scales = "free_y")


# make a bunch of boxplots
Rates %>%
  filter(day_night == "Day")%>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  ggplot(aes(x = before_after, y = rate_hr, fill = removal_control))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_boxplot()+
  facet_wrap(name~foundation_spp, scales = "free")
  
# To make the arrows in the tidepool plot

Rates %>%
  filter(name == "do_mg_l",
         day_night == "Day") %>%
  group_by(foundation_spp, before_after, removal_control) %>%
  summarise(mean_DO = mean(rate_m2_hr, na.rm = TRUE)) %>%
  pivot_wider(names_from = before_after, values_from = mean_DO) %>%
  mutate(change = (After - Before)*10)


## Make some reaction norms
mean_rates<-Rates %>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  group_by(foundation_spp, removal_control, day_night, before_after, name)%>%
  summarise(mean_value = mean(rate_hr, na.rm = TRUE),
            se_value = sd(rate_hr, na.rm = TRUE)/sqrt(n()),
            mean_value_m2 = mean(rate_m2_hr, na.rm = TRUE),
            se_value_m2 = sd(rate_m2_hr, na.rm = TRUE)/sqrt(n()))

## make a plot without the ocean data ####

mean_rates %>%
  mutate(month = ifelse(before_after == "Before", "July", "August (Upwelling)"),
         month = factor(month, levels = c("July","August (Upwelling)"))) %>%
  filter(day_night == "Day",
         foundation_spp != "Ocean") %>%
  ggplot(aes(x = month, y = mean_value_m2, color = foundation_spp, shape = removal_control, group = removal_control))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin = mean_value_m2 - se_value_m2, ymax = mean_value_m2+se_value_m2), width = 0.1)+
  facet_wrap(foundation_spp~name, scale = "free_y")

### Table with just ocean data

mean_rates %>%
  mutate(month = ifelse(before_after == "Before", "July", "August (Upwelling)"),
         month = factor(month, levels = c("July","August (Upwelling)"))) %>%
  filter(day_night == "Day",
         foundation_spp == "Ocean") 
  
### have all the before data be averaged out
mean_rates_before<-Rates %>%
  filter(before_after == "Before")%>%
  mutate(before_after = factor(before_after, levels = c("Before","After")))%>%
  group_by(foundation_spp, day_night, before_after,name)%>%
  summarise(mean_value = mean(rate_hr, na.rm = TRUE),
            se_value = sd(rate_hr, na.rm = TRUE)/sqrt(n()),
            mean_value_m2 = mean(rate_m2_hr, na.rm = TRUE),
            se_value_m2 = sd(rate_m2_hr, na.rm = TRUE)/sqrt(n()))



# repeat it so that removal and control are on there twice to make it easier to plot
mean_rates_before<-mean_rates_before %>%
  mutate(removal_control = ifelse(foundation_spp == "Ocean", "Ocean","Control"))%>%
  bind_rows(mean_rates_before %>%
              mutate(removal_control = ifelse(foundation_spp == "Ocean", "Ocean","Removal")))

# bring before and after together
mean_rates<-mean_rates %>%
  filter(before_after != "Before") %>%
  bind_rows(mean_rates_before)

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
    name == "autotrophic_pico_eukaryotes_m_l" ~"&Delta;Autotrophic <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "do_mg_l" ~ "&Delta;DO <br> (mg m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "synechoococcus_m_l" ~ "Synechoococcus <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "&Delta; Ammonium <br> (&mu;mol m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "po_umol_l" ~ "&Delta;Phosphate <br> (&mu;mol m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite <br> (&mu;mol m<sup>-2</sup>hr<sup>-1</sup>)",
    name == "humic" ~ "&Delta;Humic-like <br> (Raman units hr<sup>-1</sup>)",
    name == "prot" ~ "&Delta;Proteinaceous <br> (Raman units hr<sup>-1</sup>)",
    name == "bix"~"&Delta;BIX <br> (# m<sup>-2</sup> hr<sup>-1</sup>)" ,
    name == "hix"~"&Delta;HIX <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "m_c"~"&Delta;M:C <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "fi"~"&Delta;FI <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "total_fDOM"~"&Delta;fDOM <br> (Raman units hr<sup>-1</sup>)")
)

## reaction norm for m2 and w/o ocean ####
mean_plotdata<-mean_plot %>%
  filter(day_night == "Day",
         foundation_spp != "Ocean",
         removal_control != "Ocean",
         name %in% c("m_c","bix","hix","total_fDOM","do_mg_l", "nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l") )%>%
  mutate(month = ifelse(before_after == "Before", "July", "August (Upwelling)"),
         month = factor(month, levels = c("July","August (Upwelling)"))) %>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;DO <br> (mg m<sup>-2</sup> hr<sup>-1</sup>)",
                                                  "&Delta; Ammonium <br> (&mu;mol m<sup>-2</sup> hr<sup>-1</sup>)",
                                                  "&Delta;Nitrate+Nitrite <br> (&mu;mol m<sup>-2</sup>hr<sup>-1</sup>)",
                                                  "&Delta;Heterotrophic <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
                                                  "&Delta;BIX <br> (# m<sup>-2</sup> hr<sup>-1</sup>)" ,
                                                  "&Delta;M:C <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
                                                  "&Delta;HIX <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
                                                  "&Delta;fDOM <br> (Raman units hr<sup>-1</sup>)"
                                                  
    
  )))


mussel_rates<-mean_plotdata %>% 
 filter(foundation_spp == "Mytilus") %>%
  ggplot(aes(x = month, y = mean_value_m2, color = foundation_spp, 
             group = interaction(removal_control, foundation_spp), 
             shape = removal_control))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(x = month, y = mean_value_m2, ymin = mean_value_m2-se_value_m2, ymax = mean_value_m2+se_value_m2), width = 0.01)+
  geom_line()+
  labs(x = "",
       y = "", 
       color = "",
       shape = "",
      title = "Mussels"
       )+
  scale_color_manual(values = c("grey30"), guide = "none")+
  #scale_color_manual(values = c("grey30","#567d46"))+
  
  scale_shape_manual(values = c(16,1))+
  facet_wrap(~nicenames, scales = "free_y", ncol = 1, strip.position = "left")+
   facetted_pos_scales(
     y = rep(list(
       scale_y_continuous(limits=c(-0.1, 0.4)),
       scale_y_continuous(limits=c(0, 0.3)),
       scale_y_continuous(limits=c(-.15,0.2)),
        scale_y_continuous(limits = c(-5100, 10000)),
       scale_y_continuous(limits=c(-0.009, 0.002)),
       scale_y_continuous(limits = c(-0.015,0.005)),
       scale_y_continuous(limits = c(-0.0075,0.012))
      ), each = 1))+
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
        legend.justification = "left", 
        plot.title = element_text(hjust = 0.5, size = 14))
  
  
surfgrass_rates<-mean_plotdata %>% 
  filter(foundation_spp == "Phyllospadix") %>%
  ggplot(aes(x = month, y = mean_value_m2, color = foundation_spp, 
             group = interaction(removal_control, foundation_spp), 
             shape = removal_control))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(x = month, y = mean_value_m2, ymin = mean_value_m2-se_value_m2, ymax = mean_value_m2+se_value_m2), width = 0.01)+
  geom_line()+
  labs(x = "",
       y = "", 
       color = "",
       shape = "",
       title = "Surfgrass"
  )+
  scale_color_manual(values = c("#567d46"), guide = "none")+
  #scale_color_manual(values = c("grey30","#567d46"))+
  
  scale_shape_manual(values = c(16,1), guide = "none")+
  facet_wrap(~nicenames, scales = "free_y", ncol = 1, strip.position = "left")+
  facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits=c(-0.1, 0.4)),
      scale_y_continuous(limits=c(0, 0.3)),
      scale_y_continuous(limits=c(-.15,0.2)),
      scale_y_continuous(limits = c(-5100, 10000)),
      scale_y_continuous(limits=c(-0.009, 0.002)),
      scale_y_continuous(limits = c(-0.015,0.005)),
      scale_y_continuous(limits = c(-0.0075,0.012))
    ), each = 1))+
  theme_bw()+
  theme(
        strip.text.y.left  = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "left", 
        plot.title = element_text(hjust = 0.5, size = 14))


mussel_rates|surfgrass_rates+plot_layout(guides = "collect")&theme(legend.position = "bottom")

ggsave(filename = here("Output","AllRates.pdf"), width = 8, height = 18, device = cairo_pdf)


#### run two-way anovas
mods<-Rates %>%
  filter(day_night == "Day",
         foundation_spp != "Ocean",
         removal_control != "Ocean")%>%
#  group_by(foundation_spp, pool_id,removal_control,name) %>%
  
  # reframe(rate_diff = rate_hr[before_after == "After"] - rate_hr[before_after == "Before"]) %>%
  # ungroup() %>%
  ### removed the pools that were impacted for this analysis to see the true effect of community and upwelling
  mutate(together = paste(removal_control, before_after))%>%
  filter(together != "Removal After"
  #  removal_control == "Control"
    )%>%
   filter(name %in% c("m_c","bix", "nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "hix") )%>%
#  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l") )%>%
  group_by(name)%>%
  mutate(#rate_diff_scale = as.numeric(scale(rate_diff, scale = TRUE,center = TRUE)),
         rate_hr_scale = as.numeric(scale(rate_m2_hr, scale = TRUE,center = TRUE))) %>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) {
                       lm(rate_hr_scale  ~ foundation_spp*before_after,
                           data = df)
                     })) %>%
  
  # mutate(model = map(data, 
  #                    function(df) {
  #                      brm(rate_diff_scale  ~ removal_control,
  #                          data = df, chains = 3, iter = 10000, warmup = 5000, thin = 2)
  #                    })) %>%
  mutate(
    tidy = map(model, tidy),
    glance = map(model, glance)
  ) %>%
  mutate(nicenames = case_when(
    name == "autotrophic_pico_eukaryotes_m_l" ~"&Delta;Autotrophic <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "do_mg_l" ~ "&Delta;DO <br> (mg m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic bacteria",
    name == "synechoococcus_m_l" ~ "Synechoococcus <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "&Delta;Ammonium",
    name == "po_umol_l" ~ "&Delta;Phosphate <br> (&mu;mol m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite",
    name == "humic" ~ "&Delta;Humic-like <br> (Raman units hr<sup>-1</sup>)",
    name == "prot" ~ "&Delta;Proteinaceous <br> (Raman units hr<sup>-1</sup>)",
    name == "bix"~"&Delta;BIX" ,
    name == "hix"~"&Delta;HIX",
    name == "m_c"~"&Delta;M:C",
    name == "fi"~"&Delta;FI <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "total_fDOM"~"&Delta;fDOM <br> (Raman units hr<sup>-1</sup>)")
  )%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;Ammonium",
                                                  "&Delta;Nitrate+Nitrite",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic bacteria"
                                                  
  )))


r1<-mods%>%
  unnest(tidy)%>%
  filter(#p.value <=  0.05,
  #  effect == "fixed",
    term != "(Intercept)",
   # term != "foundation_sppPhyllospadix:before_afterBefore" # there is no sig interaction so remove
   #  name !="heterotrophic_bacterioplankton_m_l"
  ) %>%
   mutate(term = case_when(term == "before_afterBefore"~"Upwelling Effect",
                          term == "foundation_sppPhyllospadix"~"Foundation spp. Effect",
                          term =="foundation_sppPhyllospadix:before_afterBefore"~"Interaction"))%>%
  mutate(term = factor(term,levels = c("Foundation spp. Effect",
                                       "Upwelling Effect",
                                       "Interaction")))%>%
  mutate(alpha = ifelse(p.value <= 0.05, 1, 0.5))%>%
  ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha))+
  geom_point(size = 3)+
  geom_errorbar(aes(xmin = estimate-std.error, xmax = estimate+std.error), width = 0.1)+
  geom_vline(xintercept=0)+
  labs(x = "Standardized effect size <br> (rates m<sup>-2</sup> hr<sup>-1</sup>)",
       y = "")+
  facet_wrap(~term, scale = "free_y", ncol = 1)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_markdown(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_markdown(size = 14),
        legend.position = "none")
#ggsave(here("Output","Effects_ControlOnly.png"), width = 6, height = 8)

### Do this again with the real BACI design
mods_BACI<-Rates %>%
  filter(day_night == "Day",
         foundation_spp != "Ocean",
         removal_control != "Ocean")%>%
  group_by(foundation_spp, pool_id,removal_control,name) %>%
   reframe(rate_diff = rate_m2_hr[before_after == "After"] - rate_m2_hr[before_after == "Before"]) %>%
    ungroup() %>%
  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "m_c","bix","hix") )%>%
  group_by(foundation_spp, name)%>%
  mutate(rate_diff_scale = as.numeric(scale(rate_diff, scale = TRUE,center = TRUE))
  #  rate_hr_scale = as.numeric(scale(rate_m2_hr, scale = TRUE,center = TRUE))
  ) %>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) {
                       lm(rate_diff_scale  ~ removal_control,
                           data = df)
                     })) %>%
    mutate(
    tidy = map(model, tidy),
    glance = map(model, glance)
  ) %>%
  mutate(nicenames = case_when(
    name == "autotrophic_pico_eukaryotes_m_l" ~"&Delta;Autotrophic <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "do_mg_l" ~ "&Delta;DO <br> (mg m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria",
    name == "synechoococcus_m_l" ~ "Synechoococcus <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "&Delta; Ammonium",
    name == "po_umol_l" ~ "&Delta;Phosphate <br> (&mu;mol m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite",
    name == "humic" ~ "&Delta;Humic-like <br> (Raman units hr<sup>-1</sup>)",
    name == "prot" ~ "&Delta;Proteinaceous <br> (Raman units hr<sup>-1</sup>)",
    name == "bix"~"&Delta;BIX" ,
    name == "hix"~"&Delta;HIX",
    name == "m_c"~"&Delta;M:C",
    name == "fi"~"&Delta;FI <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "allfDOM"~"&Delta;fDOM <br> (Raman units hr<sup>-1</sup>)")
  )%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta; Ammonium",
                                                  "&Delta;Nitrate+Nitrite",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic Bacteria"
                                                  
  )))


r2<-mods_BACI%>%
  unnest(tidy)%>%
  filter(#p.value <=  0.05,
         #  effect == "fixed",
         term != "(Intercept)",
         #  name !="heterotrophic_bacterioplankton_m_l"
  ) %>%
  mutate(alpha = ifelse(p.value<= 0.05,1, 0.5))%>%
  ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha))+
  geom_point(size = 3)+
  geom_errorbar(aes(xmin = estimate-std.error, xmax = estimate+std.error), width = 0.1)+
  geom_vline(xintercept=0)+
  labs(x = "Standardized effect size <br> (rates m<sup>-2</sup> hr<sup>-1</sup>)",
       y = "")+
  facet_wrap(~foundation_spp, scale = "free_y", ncol = 1)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_markdown(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_markdown(size = 14),
        legend.position = "none")
#ggsave(here("Output","Effects_BACI.png"), width = 6, height = 8)

### run ANOVA with the raw data

data_end<-data_all %>%
  ungroup()%>%
  filter(day_night == "Day",
         foundation_spp != "Ocean"
         )  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like) %>%
  select(month,pool_id, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot, humic,m_c, bix, hix,fi, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  
  group_by(foundation_spp, pool_id,removal_control,name, time_point) %>%
   reframe(value_diff = value[month == "August (Upwelling)"] - value[month == "July"]) %>%
   ungroup() %>%

  
  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "m_c","bix","hix") )%>%
  #filter(month == "August (Upwelling)")%>%
#  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l") )%>%
  group_by(foundation_spp, name)%>%
  #mutate(value_diff_scale = as.numeric(scale(value_diff, scale = TRUE,center = TRUE))) %>%
  mutate(value_scale = as.numeric(scale(value_diff, scale = TRUE,center = TRUE))) %>%
  ungroup()%>%
   mutate(nicenames = case_when(
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria",
    name == "nh4_umol_l" ~ "&Delta;Ammonium",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite",
    name == "bix"~"&Delta;BIX" ,
    name == "hix"~"&Delta;HIX",
    name == "m_c"~"&Delta;M:C",
    )
  )%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;Ammonium",
                                                  "&Delta;Nitrate+Nitrite",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic Bacteria"
                                                  
  )))


  #filter(time_point == "start")

mods2<-data_end %>%
  filter(time_point == "start")%>%
#  mutate(value_scale_log = sign(value_scale)*sqrt(abs(value_scale)))%>%
  group_by(foundation_spp,nicenames, time_point)%>%
  nest() %>%
  # mutate(model = map(data, 
  #                    function(df) {
  #                      brm(value_scale  ~ removal_control,
  #                          data = df, chains = 3, iter = 10000, warmup = 5000, thin = 2)
  #                    })) %>%
  mutate(model = map(data,
                     function(df) {
                       lm(value_scale~ removal_control, data = df) #log transformed log(abs(value_diff_scale))*sign(value_diff_scale)
                     })) %>%
  mutate(
    tidy = map(model, tidy),
    glance = map(model, glance)
  )


conc2<-mods2%>%
  unnest(tidy)%>%
  filter(
    #p.value <=  0.05,
 #   effect == "fixed",
    term != "(Intercept)",
    #  name !="heterotrophic_bacterioplankton_m_l"
  ) %>%
  mutate(alpha = ifelse(p.value <= 0.055, 1, 0.5))%>%
  ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha))+
  geom_point(size = 3)+
  geom_errorbar(aes(xmin = estimate-std.error, xmax = estimate+std.error), width = 0.1)+
  geom_vline(xintercept = 0)+
  labs(x = "Standardized effect size <br> (concentration or density)",
       y = "")+
  facet_wrap(~foundation_spp, ncol = 1)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_markdown(size = 14),
        legend.position = "none")


r2|conc2
ggsave(here("Output","BACIEffects"), width = 8, height = 8, device = cairo_pdf)


#### concentration only for control pools 
data_end2<-data_all %>%
  ungroup()%>%
  filter(day_night == "Day",
         foundation_spp != "Ocean"
  )  %>%
  dplyr::select(before_after,pool_id, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,m_c, bix, hix,total_fDOM, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  group_by(foundation_spp, pool_id,removal_control,name, time_point) %>%
  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "m_c","bix","hix") )%>%
  mutate(together = paste(removal_control, before_after))%>%
  filter(together != "Removal After"
         #  removal_control == "Control"
  )%>%
  #filter(month == "August (Upwelling)")%>%
  #  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l") )%>%
  group_by(name)%>%
  #mutate(value_diff_scale = as.numeric(scale(value_diff, scale = TRUE,center = TRUE))) %>%
  mutate(value_scale = as.numeric(scale(value, scale = TRUE,center = TRUE))) %>%
  ungroup()%>%
  mutate(nicenames = case_when(
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic Bacteria",
    name == "nh4_umol_l" ~ "Ammonium",
    name == "nn_umol_l" ~ "Nitrate+Nitrite",
    name == "bix"~"BIX" ,
    name == "hix"~"HIX",
    name == "m_c"~"M:C",
  )
  )%>%
  mutate(nicenames = factor(nicenames, levels = c("Ammonium",
                                                  "Nitrate+Nitrite",
                                                  "BIX" ,
                                                  "M:C",
                                                  "HIX",
                                                  "Heterotrophic Bacteria"
                                                  
  )))



mods3<-data_end2 %>%
  filter(time_point == "start")%>%
  group_by(nicenames)%>%
  nest() %>%
   mutate(model = map(data,
                     function(df) {
                       lm(value_scale~ foundation_spp*before_after, data = df) #log transformed log(abs(value_diff_scale))*sign(value_diff_scale)
                     })) %>%
  mutate(
    tidy = map(model, tidy),
    glance = map(model, glance)
  )



conc1<-mods3%>%
  unnest(tidy)%>%
  filter(#p.value <=  0.05,
    #  effect == "fixed",
    term != "(Intercept)",
  #  term != "foundation_sppPhyllospadix:before_afterBefore" # there is no sig interaction so remove
    #  name !="heterotrophic_bacterioplankton_m_l"
  ) %>%
  mutate(term = case_when(term == "before_afterBefore"~"Upwelling Effect",
                          term == "foundation_sppPhyllospadix"~"Foundation spp. Effect",
                          term =="foundation_sppPhyllospadix:before_afterBefore"~"Interaction"))%>%
  mutate(term = factor(term,levels = c("Foundation spp. Effect",
                                       "Upwelling Effect",
                                       "Interaction")))%>%
  mutate(alpha = ifelse(p.value <= 0.05, 1, 0.5))%>%
  ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha))+
  geom_point(size = 3)+
  geom_errorbar(aes(xmin = estimate-std.error, xmax = estimate+std.error), width = 0.1)+
  geom_vline(xintercept=0)+
  labs(x = "Standardized effect size <br> (concentration or density)",
       y = "")+
  facet_wrap(~term, scale = "free_y", ncol = 1)+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_markdown(size = 14),
        legend.position = "none")


r1|conc1
ggsave(here("Output","ControlOnlyEffects.pdf"), width = 8, height = 8, device = cairo_pdf)

### THEN JUST THE AFTER SET WITH AND W/O BACI

### TRY JUST LOOKING AT CONTROL POOLS FOR SPECIES AND UPWELLING IMPACTS, THEN JUST THE AFTER SET WITH AND W/O BACI


### normalize to the ocean ####
data_norm<-data_all %>%
  ungroup()%>%
  filter(day_night == "Day",
          foundation_spp == "Ocean"
  )  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like) %>%
  select(month,pool_id, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot,bix, m_c, hix, fi, humic, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  group_by(name, month, time_point) %>%
  summarise(ocean_mean = mean(value, na.rm = TRUE)) %>%
  left_join(data_all %>%
              ungroup()%>%
              filter(day_night == "Day",
                     foundation_spp != "Ocean"
              )  %>%
              mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
                     humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like) %>%
              select(month,pool_id, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot,bix, m_c, hix, fi, humic, nn_umol_l, nh4_umol_l) %>%
              pivot_longer(cols = do_mg_l:nh4_umol_l)) %>%
  mutate(value_norm = value - ocean_mean)


### make the rates but normalized to the ocean
rate_norm<-data_all%>%
  ungroup()%>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         allfDOM = prot+humic) %>%
  #mutate(sampling_datetime = as.numeric(mdy_hms(paste(sampling_day, sampling_time)))) %>% 
  select(pool_id:removal_control, time_point, do_mg_l, po_umol_l:nh4_umol_l, heterotrophic_bacterioplankton_m_l:fi, prot, humic,allfDOM) %>%
  mutate(heterotrophic_bacterioplankton_m_l = heterotrophic_bacterioplankton_m_l *1000, # convert to per L
         autotrophic_pico_eukaryotes_m_l =  autotrophic_pico_eukaryotes_m_l*1000,
         synechoococcus_m_l = synechoococcus_m_l*1000
  )%>%
  pivot_longer(cols = do_mg_l:allfDOM)%>%
  group_by(name, before_after, time_point) %>%
  summarise(ocean_mean = mean(value, na.rm = TRUE)) %>%
  left_join(data_all%>%
              ungroup()%>%
              mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
                     humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
                     allfDOM = prot+humic) %>%
              #mutate(sampling_datetime = as.numeric(mdy_hms(paste(sampling_day, sampling_time)))) %>% 
              select(pool_id:removal_control, time_point, do_mg_l, po_umol_l:nh4_umol_l, heterotrophic_bacterioplankton_m_l:fi, prot, humic,allfDOM) %>%
              mutate(heterotrophic_bacterioplankton_m_l = heterotrophic_bacterioplankton_m_l *1000, # convert to per L
                     autotrophic_pico_eukaryotes_m_l =  autotrophic_pico_eukaryotes_m_l*1000,
                     synechoococcus_m_l = synechoococcus_m_l*1000
              )%>%
              pivot_longer(cols = do_mg_l:allfDOM)) %>%
  filter(foundation_spp != "Ocean")%>%
  mutate(value_norm = value - ocean_mean) %>%
  group_by(pool_id, before_after, removal_control,day_night, foundation_spp, sampling_group, name) %>%
  reframe(change = change_val(value_norm, time_point)) %>% # calculate the difference between start and end
  left_join(times)%>%# join with the times
  left_join(MetaData %>%
              clean_names %>%
              mutate(pool_id = as.character(pool_id)) %>%
              select(pool_id, before_after, surface_area,vol))%>% # add in the tide pool info
  mutate(rate_hr = change/diff_time,# difference in value per hour
         rate_m2_hr = rate_hr*vol/surface_area
  )

#### run two-way anovas
rate_norm %>%
  filter(day_night == "Day",
         foundation_spp != "Ocean",
         removal_control != "Ocean")%>%
  group_by(foundation_spp, pool_id,removal_control,name) %>%
  reframe(Rate_diff = rate_m2_hr[before_after == "After"] - rate_m2_hr[before_after == "Before"]) %>%
  ungroup() %>%
  filter(name %in% c("m_c","bix", "nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "fi") )%>%
  ungroup()%>%
  group_by(foundation_spp,name )%>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) {
                       lm(Rate_diff ~ removal_control, data = df)
                     })) %>%
  mutate(
    tidy = map(model, tidy),
    glance = map(model, glance)
  )%>%
  unnest(glance) %>%
  filter(p.value <=  0.05)


####

data_norm %>%
  group_by(foundation_spp, pool_id,removal_control,name, time_point) %>%
  reframe(value_diff = value_norm[month == "August (Upwelling)"] - value_norm[month == "July"]) %>%
  ungroup() %>%
  group_by(foundation_spp,name, time_point)%>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) {
                       lm(value_diff ~ removal_control, data = df)
                     })) %>%
  mutate(
    tidy = map(model, tidy),
    glance = map(model, glance)
  )%>%
  unnest(glance) %>%
  filter(p.value <=  0.05)

data_norm %>%
  group_by(foundation_spp, pool_id,removal_control,name, time_point) %>%
  reframe(value_diff = value_norm[month == "August (Upwelling)"] - value_norm[month == "July"]) %>%
  ungroup() %>%
  ggplot(aes(x = time_point, y = value_diff, color = removal_control))+
  geom_point()+
  facet_wrap(name~foundation_spp, scale = "free")

#### ONE SAMPLE T-TEST TO SEE OF THE VALUES ARE DIFFERENT FROM THE OCEAN i.e. not 0
data_norm %>%
  group_by(foundation_spp, month, removal_control,name, time_point) %>%
  nest() %>%
  mutate(model = map(data, 
                     ~t.test(.x$value_norm, mu = 0))) %>%
  mutate(
    tidy = map(model, tidy),
    glance = map(model, glance)
  )%>%
  unnest(glance) %>%
  filter(p.value <=  0.05,
         time_point == "start") %>%
  ggplot(aes(x = month, color = removal_control, y = estimate))+
  geom_point()+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1)+
  geom_hline(yintercept = 0)+
  facet_wrap(name~foundation_spp, scale = "free_y")


###### rates norm
Rates %>%
  group_by(foundation_spp, before_after, removal_control,name) %>%
  filter(
         foundation_spp != "Ocean")%>%
  nest() %>%
  mutate(model = map(data, 
                     ~t.test(.x$rate_hr, mu = 0))) %>%
  mutate(
    tidy = map(model, tidy),
    glance = map(model, glance)
  )%>%
  unnest(glance) %>%
 # filter(p.value <=  0.05) %>%
  ggplot(aes(x = before_after, color = removal_control, y = estimate))+
  geom_point()+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1)+
  geom_hline(yintercept = 0)+
  facet_wrap(name~foundation_spp, scale = "free_y")


###########

test<-Rates %>%
  filter(day_night == "Day",
         foundation_spp == "Phyllospadix",
         removal_control != "Ocean",
         name =="m_c") %>%
  mutate(before_after = factor(before_after, labels = c("Before","After")))

# fit<-brm(formula = rate_hr ~ before_after*removal_control,
#     data = test, warmup = 2000, iter = 10000 )

# summary(fit)
# plot(fit)
# # plot conditional effects for each predictor
# plot(conditional_effects(fit), ask = FALSE)

mc<-lmer(rate_hr ~ removal_control+before_after+(1|pool_id), data =Rates %>%
         filter(day_night == "Day",
                foundation_spp == "Phyllospadix",
                removal_control != "Ocean",
                name =="m_c" ) %>%
         mutate(removal_control = ifelse(before_after == "Before","Control",removal_control)))


anova(mc)
summary(mc)

## Make a plot of the benthic data

### order the pools to make it prettier
PoolOrder<-as_tibble(list(PoolOrder = c(1:32), 
                     PoolID = factor(c(4,1,20, 18, 5, 27,8, 29,
                                   21,26,6,3,2,28,22,7,13,31,14,19,30,9,25,
                                    12,32,15,11,10,17,16,23,24))))%>%
  mutate(PoolOrder2 = c(1:16,1:16))



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
  select(PoolID, Before_After, Removal_Control, Foundation_spp, Macroalgae = macroalgae , Diatoms, `Non-mussel Consumers` = consumers, CCA = allCCA, Mussels = AdjMusselCover, Surfgrass = AdjSurfgrassCover)%>%
  mutate(PoolID =  as.factor(PoolID)) %>%
  mutate(`Rock/Sand` = (100-(Macroalgae+`Non-mussel Consumers`+CCA+Mussels+Surfgrass)),
         Macroalgae = Macroalgae - Diatoms,
         PoolID = fct_reorder2(PoolID,`Rock/Sand`, Mussels)) %>% # the macroalgae includes diatom cover and I want to see the difference
  left_join(PoolOrder)%>%
  pivot_longer(cols = Macroalgae:`Rock/Sand`) %>%
  mutate(name = factor(name, levels = c("Surfgrass","Macroalgae","Mussels","Non-mussel Consumers", "Diatoms","CCA","Rock/Sand")),
         Foundation_spp = ifelse(Foundation_spp == "Mytilus","Mussel-dominated","Surfgrass-dominated"),
         Before_After = ifelse(Before_After == "Before","Before-Impact","After-Impact"),
         Before_After = factor(Before_After, levels = c("Before-Impact","After-Impact")))%>%
  ggplot(aes(x = PoolOrder2, y = value, fill = name))+
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cal_palette("tidepool", n = 7, type = "continuous"))+
  geom_vline(aes(xintercept = 8.5), linewidth = 1.25)+
  labs(x = "",
       y = "% Cover",
       fill = "")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16)
        )+
  facet_grid(Foundation_spp~Before_After, scale = "free_x")

ggsave(here("Output","BenthicCover.pdf"), width = 10, height = 8, device = cairo_pdf)

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
Day_rates_clean<-Day_rates %>%  
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
    
  )))

P_Nuts <- Day_rates_clean %>%
ggplot(aes(x = allproddom, y = rate_hr))+
  geom_vline(aes(xintercept  = 0))+
  geom_hline(aes(yintercept  = 0))+
  geom_point(aes(color =  month), alpha = 0.5)+
  geom_smooth(method = "lm",
              formula = 'y~poly(x,2)', color = "black",
              data = Day_rates_clean %>% # only plot the significant models
                filter(name %in% c("synechoococcus_m_l","po_umol_l","nn_umol_l","nh4_umol_l","do_mg_l")))+
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


## lmer models --  pool id to account for repeated measure as random effect
data_anova<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  filter(!name %in% c("fi","bix","hix","m_c") ) %>%
  select(pool_id, before_after,allproddom, name, rate_hr) %>%
  nest_by(name) %>%
  mutate(mod = list(lmer(rate_hr~poly(allproddom,2)+(1|pool_id), data = data)),
         modstat = list(broom::glance(mod)),
         res =  list(broom::tidy(mod)),
         ano = list(broom::tidy(anova(mod))))
 
data_anova %>%
  select(name, modstat, res, ano) %>%
  unnest(ano) %>%
  filter(term != "Residuals") %>%
  filter(p.value <= 0.05) # pull out just the significant values


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
       x = "heterotrophic bacteria (counts mL-1 hr-1)",
       shape = "Foundation Species")+
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

Day_rates_wide2<-Day_rates %>%  
  filter(foundation_spp != "Ocean") %>%
  select(-c(change:vol, rate_hr))%>%
  pivot_wider(values_from = rate_m2_hr, names_from = name)


##DO and Hetero
P_HDO<-Day_rates_wide2 %>%  
  mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                             removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                             removal_control == "Removal"& month != "July" ~ "Foundation spp. removed"))%>%
  mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation spp. removed")))%>%
  filter(foundation_spp != "Ocean")%>%
  mutate(foundation_spp = ifelse(foundation_spp == "Mytilus", "Mussel-dominated","Surfgrass-dominated"))%>%
  ggplot(aes(y = heterotrophic_bacterioplankton_m_l, x = do_mg_l))+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5, aes(color = month, shape = foundation_spp))+
  # geom_text(aes(x = 175, y = 0.1, label = "p = 0.036"))+
#  ylim(-.3,.3)+
  geom_smooth(method = "lm", color = "black",
              data = Day_rates_wide2 %>%  
                mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                                           removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                                           removal_control == "Removal"& month != "July" ~ "Foundation spp. removed"))%>%
                mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation spp. removed")))%>%
                filter(foundation_spp != "Ocean",
                       removal == "Unmanipulated"))+
  scale_color_manual(values = c("grey","grey8"))+
  labs(y = "Heterotrophic bacterial production <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
       x = "DO production <br> (mg m<sup>-2</sup> hr<sup>-1</sup>)",
       color = "Sampling Month",
       shape = " Foundation Species")+
  facet_wrap(~removal, scales = "free_x")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_markdown(size = 14),
        axis.title.y = element_markdown(size = 14),
        axis.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
  )
ggsave(here("Output","het_DO_regression.pdf"), width = 8, height = 4)


DO_het_mod<-lmer(data = Day_rates_wide2 %>%  
                 mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                                            removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                                            removal_control == "Removal"& month != "July" ~ "Foundation spp. removed"))%>%
                 mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation spp. removed")))%>%
                 filter(foundation_spp != "Ocean"), heterotrophic_bacterioplankton_m_l~do_mg_l*removal+(1|pool_id))

anova(DO_het_mod)
summary(DO_het_mod)

#### Het and Syn
##DO and Hetero
P_HPico<-Day_rates_wide %>%  
  mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                             removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                             removal_control == "Removal"& month != "July" ~ "Foundation spp. removed"))%>%
  mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation spp. removed")))%>%
  filter(foundation_spp != "Ocean")%>%
  ggplot(aes(x = heterotrophic_bacterioplankton_m_l, y = autotrophic_pico_eukaryotes_m_l))+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_point(alpha = 0.5,aes(color = month, shape = foundation_spp))+
  # geom_text(aes(x = 175, y = 0.1, label = "p = 0.036"))+
  #  ylim(-.3,.3)+
  geom_smooth(method = "lm", 
              data = Day_rates_wide %>%  
                mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                                           removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                                           removal_control == "Removal"& month != "July" ~ "Foundation spp. removed"))%>%
                mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation spp. removed")))%>%
                filter(foundation_spp != "Ocean",
                       removal == "Unmanipulated"), color = "black")+
  scale_color_manual(values = c("grey","grey8"))+
  labs(x = "Heterotrophic bacterial production <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
       y = "Auotrophic picoeukaryotes production <br> (# mL<sup>-1</sup> hr<sup>-1</sup>)",
       color = "")+
  facet_wrap(~removal, scales = "free_x")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_markdown(size = 14),
        axis.title.y = element_markdown(size = 14),
        axis.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14)
  )

DO_het_pico<-lm(data = Day_rates_wide %>%  
                 mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                                            removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                                            removal_control == "Removal"& month != "July" ~ "Foundation spp. removed"))%>%
                 mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation spp. removed")))%>%
                 filter(foundation_spp != "Ocean"), heterotrophic_bacterioplankton_m_l~autotrophic_pico_eukaryotes_m_l*removal)

anova(DO_het_pico)
summary(DO_het_pico)

P_HDO/P_HPico+plot_layout(guides = "collect")
ggsave(here("output","top_bottom_regress.pdf"), width = 8, height = 8, useDingbats = FALSE)


