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
  dplyr::select(PoolID, Foundation_spp, Day_Night, Before_After, Removal_Control, Group, Time_Point,
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
  group_by(before_after, time_point, day_night)%>%
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
