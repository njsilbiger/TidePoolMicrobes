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
library(scales)


### read in data #########
BenthicData<-read_csv(here("Data","Microbe_Clean","CommunityData.csv"))

MetaData<-read_csv(here("Data","Microbe_Clean","TidePoolDescriptions.csv"))

data_all<-read_csv(here("Data","Microbe_Clean","joinedData_edited.csv"))

# read in the BUETI 
#BT<-read_csv(here("Data","Microbe_Clean", "BUETI.csv"))
#mmol m-1 s-1
CT<-read_csv(here("Data","Microbe_Clean",  "CUTI.csv"))
#m2 s-1
# https://oceanview.pfeg.noaa.gov/products/upwelling/cutibeuti

#CUTI and BEUTI
#Jacox, M. G., C. A. Edwards, E. L. Hazen, and S. J. Bograd (2018) Coastal upwelling revisited: Ekman, Bakun, and improved upwelling indices for the U.S. west coast. Journal of Geophysical Research, doi:10.1029/2018JC014187. 
#Coastal Upwelling Transport Index

reverse2_trans <- function() {
  trans_new(
    "reverse2",
    function(x) -1 * as.numeric(x), # Force values to be numeric for Date objects
    function(x) -1 * as.numeric(x)
  )
}

CT %>%
  filter(latitude == 45) %>%
  mutate(date = as_date(date)) %>%
  ggplot(aes(x = date, y = CUTI))+
  geom_hline(aes(yintercept = 0), linetype = 2)+
  geom_ribbon(aes(xmin = ymd("2019-07-08"), 
                  xmax =ymd("2019-07-09") ), alpha = 0.2, fill = "red")+
  geom_ribbon(aes(xmin = ymd("2019-08-05"), 
                  xmax =ymd("2019-08-06") ), alpha = 0.2, fill = "red")+
  geom_point()+
  geom_point(data = CT %>%
               mutate(date = as_date(date)) %>%
               filter(latitude == 45,
                      date >= ymd("2019-07-08")&
                        date <= ymd("2019-07-09") |
                        date >= ymd("2019-08-05")&
                        date <= ymd("2019-08-06")
                        ),
             color = "red")+
  geom_path()+
  labs(x = "",
       y = "CUTI m<sup>2</sup> s<sup>-1</sup>")+
  coord_flip()+
  scale_x_continuous(trans = c("date", "reverse2"))+
  theme_bw()+
  theme(axis.title = element_markdown(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank())
 
ggsave(here("Output","CUTI.pdf"), height = 8, width = 6)



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

# make it long 

data_long<-data_all %>%
  ungroup()%>%
  filter(removal_control != "Removal") %>%
  filter(day_night == "Day")  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = prot+humic) %>%
  select(month, pool_id,day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot, humic,total_fDOM, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  group_by(month, name, foundation_spp,time_point, day_night, removal_control)%>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val= sd(value, na.rm = TRUE)/sqrt(n())) %>%
  mutate(time_point_clean = ifelse(time_point == "start","Early", "Late"))%>%
  mutate(nicenames = case_when(
     name == "autotrophic_pico_eukaryotes_m_l" ~"Autotrophic <br> (# mL<sup>-1</sup>)",
    name == "do_mg_l" ~ "DO <br> (mg L<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic <br> (# mL<sup>-1</sup>)",
      name == "synechoococcus_m_l" ~ "Synechoococcus <br> (# mL<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
    name == "po_umol_l" ~ "Phosphate <br> (&mu;mol L<sup>-1</sup>)",
    name == "nn_umol_l" ~ "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
      name == "humic" ~ "Humic-like <br> (Raman units)",
    name == "total_fDOM" ~"Total fDOM <br> (Raman units)",
      name == "prot" ~ "Proteinaceous <br> (Raman units)",
    name == "hix"~"HIX",
    name == "bix"~"BIX",
    name == "m_c"~"M:C")
  ) %>%
  mutate(nicenames = factor(nicenames, levels = c("DO <br> (mg L<sup>-1</sup>)",
                                                  "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Phosphate <br> (&mu;mol L<sup>-1</sup>)",
                                                    "Autotrophic <br> (# mL<sup>-1</sup>)",
                                                    "Synechoococcus <br> (# mL<sup>-1</sup>)",
                                                  "Heterotrophic <br> (# mL<sup>-1</sup>)",
                                                   "Humic-like <br> (Raman units)",
                                                   "Proteinaceous <br> (Raman units)",
                                                  "Total fDOM <br> (Raman units)",
                                                  "HIX","BIX","M:C"
                                                  
                                                  
  )))

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

### run models of rates with a BACI design

### Do this again with the real BACI design
mods_BACI<-Rates %>%
  filter(day_night == "Day",
         foundation_spp != "Ocean",
         removal_control != "Ocean")%>%
  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "m_c","bix","hix","do_mg_l") )%>%
  mutate(rate_m2_hr_sqrt = sign(rate_m2_hr)*sqrt(abs(rate_m2_hr)),
         rate_m2_hr_log = sign(rate_m2_hr)*log(abs(rate_m2_hr))
         )%>%
  group_by(foundation_spp, name)%>%
  mutate(#rate_diff_scale = as.numeric(scale(rate_m2_hr_sqrt, scale = TRUE,center = TRUE))
         rate_diff_scale = as.numeric(scale(rate_m2_hr_sqrt, scale = TRUE,center = TRUE))
         #  rate_hr_scale = as.numeric(scale(rate_m2_hr, scale = TRUE,center = TRUE))
  ) %>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) {
                       lmer(rate_diff_scale ~ removal_control*before_after +(1|pool_id),
                          data = df)
                     })) %>%
  mutate(
    tidy = map(model, function(x)tidy(x,effects="fixed")),
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
                                                  "&Delta;Heterotrophic Bacteria",
                                                  "&Delta;fDOM"
                                                  
  )))

# extract the coefs for the interaction term
r2<-mods_BACI%>%
  unnest(tidy)%>%
  filter(#p.value <=  0.05,
    #  effect == "fixed",
    term =="removal_controlRemoval:before_afterBefore",
    #  name !="heterotrophic_bacterioplankton_m_l"
  ) %>%
  mutate(alpha = ifelse(p.value<= 0.055,1, 0.6))

write_csv(mods_BACI%>%
            unnest(tidy), file = here("Output","Rates_models.csv"))

## Do the models with the avaerage concentration 
values<-data_all %>%
  ungroup()%>%
  filter(day_night == "Day",
         foundation_spp != "Ocean"
  )  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = humic+prot) %>%
  select(month,pool_id, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot,humic,m_c,total_fDOM, bix, hix,fi, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  group_by(foundation_spp, pool_id,removal_control,name, month) %>%
  summarise(mean_val = mean(value, na.rm = TRUE))%>%
  mutate(mean_val_sqrt = sqrt(mean_val),
         mean_val_log = log(mean_val))%>%
  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "do_mg_l","m_c","bix","hix") )%>%
  group_by(foundation_spp, name)%>%
#   mutate(value_scale = as.numeric(scale(mean_val_sqrt, scale = TRUE,center = TRUE))) %>%
  mutate(value_scale = as.numeric(scale(mean_val_sqrt, scale = TRUE,center = TRUE))) %>%
  ungroup()%>%
  mutate(nicenames = case_when(
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria",
    name == "nh4_umol_l" ~ "&Delta;Ammonium",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite",
    name == "bix"~"&Delta;BIX" ,
    name == "hix"~"&Delta;HIX",
    name == "m_c"~"&Delta;M:C",
    name == "total_fDOM"~"&Delta; Total fDOM"
  )
  )%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;Ammonium",
                                                  "&Delta;Nitrate+Nitrite",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic Bacteria",
                                                  "&Delta;fDOM"
                                                  
  )))

# run the models
mods2<-values %>%
  mutate(before_after = ifelse(month == "July", "Before","After"))%>%
   # mutate(value_scale_sqrt = sign(value_scale)*sqrt(abs(value_scale)))%>%
  group_by(foundation_spp,nicenames)%>%
  nest() %>%
  mutate(model = map(data,
                     function(df) {
                       lmer(value_scale~ removal_control*before_after+(1|pool_id), data = df) #log transformed log(abs(value_diff_scale))*sign(value_diff_scale)
                     })) %>%
  mutate(
    map(model, function(x)tidy(x,effects="fixed")),
    glance = map(model, glance)
  ) %>%
  rename(tidy = `map(model, function(x) tidy(x, effects = "fixed"))`)

# tidy the concentrations
conc2<-mods2%>%
  unnest(tidy)%>%
  filter(#p.value <=  0.05,
    #  effect == "fixed",
    term =="removal_controlRemoval:before_afterBefore",
    #  name !="heterotrophic_bacterioplankton_m_l"
  ) %>%
  mutate(alpha = ifelse(p.value<= 0.055,1, 0.6))

write_csv(mods2%>%
            unnest(tidy), file = here("Output","Concentration_models.csv"))

### make a plot of the interaction terms
con_int<-conc2 %>%
ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha, color = foundation_spp))+
  scale_color_manual(values = c("black","#34c230"))+
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

r_int<-r2 %>%
  ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha, color = foundation_spp))+
  scale_color_manual(values = c("black","#34c230"))+
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

BC_int<-r_int|con_int
ggsave(here("Output","BACIEffects_interaction.pdf"), width = 8, height = 8, device = cairo_pdf)

# To get the upwelling effect we cal only look at the unmanipulated pools 
## do a two-way ANOVA between month and foundation species to see effect of upwelling and species identity

Unmamipulated_mean<-data_all %>%
  ungroup()%>%
  filter(day_night == "Day",
         foundation_spp != "Ocean")  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = prot+humic) %>%
  select(month, pool_id,before_after, removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot,m_c, bix,hix, humic,total_fDOM, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "m_c","bix","hix","do_mg_l") )%>%
  group_by(foundation_spp, pool_id,removal_control, before_after,name) %>%
  summarise(mean_val = mean(value, na.rm = TRUE)) %>% # get mean for the low tide
  mutate(together = paste(removal_control, before_after),
         mean_val_sqrt = sqrt(mean_val))%>%
  filter(together != "Removal After"
         #  removal_control == "Control"
  )%>%
  group_by(name)%>%
  #mutate(value_diff_scale = as.numeric(scale(value_diff, scale = TRUE,center = TRUE))) %>%
  mutate(value_scale = as.numeric(scale(mean_val_sqrt, scale = TRUE,center = TRUE))) %>%
  ungroup()%>%
  mutate(nicenames = case_when(
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic Bacteria",
    name == "nh4_umol_l" ~ "Ammonium",
    name == "nn_umol_l" ~ "Nitrate+Nitrite",
    name == "bix"~"BIX" ,
    name == "hix"~"HIX",
    name == "m_c"~"M:C",
    name == "total_fDOM"~"Total fDOM"
  )
  )%>%
  mutate(nicenames = factor(nicenames, levels = c("Ammonium",
                                                  "Nitrate+Nitrite",
                                                  "BIX" ,
                                                  "M:C",
                                                  "HIX",
                                                  "Total fDOM",
                                                  "Heterotrophic Bacteria"
                                                  
  )))

# run the models
mods3<-Unmamipulated_mean %>%
  group_by(nicenames)%>%
  nest() %>%
  mutate(model = map(data,
                     function(df) {
                       lm(value_scale~ foundation_spp*before_after, data = df) #log transformed log(abs(value_diff_scale))*sign(value_diff_scale)
                     })) %>%
  mutate(
    tidy = map(model, function(x)tidy(x,effects="fixed")),
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
  mutate(alpha = ifelse(p.value <= 0.055, 1, 0.5))%>%
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


write_csv(mods3%>%
            unnest(tidy), here("Output","Concentration_controlonlymods.csv"))

#### do it again for the rates
mods<-Rates %>%
  filter(day_night == "Day",
         foundation_spp != "Ocean",
         removal_control != "Ocean")%>%
  mutate(together = paste(removal_control, before_after),
         rate_m2_hr_sqrt = sign(rate_m2_hr)*sqrt(abs(rate_m2_hr)))%>%
  filter(together != "Removal After"
         #  removal_control == "Control"
  )%>%
  filter(name %in% c("m_c","bix", "nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "hix","do_mg_l") )%>%
   group_by(name)%>%
  mutate(#rate_diff_scale = as.numeric(scale(rate_diff, scale = TRUE,center = TRUE)),
      rate_hr_scale = as.numeric(scale(rate_m2_hr_sqrt, scale = TRUE,center = TRUE))) %>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) {
                       lm(rate_hr_scale  ~ foundation_spp*before_after,
                          data = df)
                     })) %>%
  
   mutate(
    tidy = map(model, function(x)tidy(x,effects="fixed")),
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
  mutate(alpha = ifelse(p.value <= 0.055, 1, 0.5))%>%
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

write_csv(mods%>%
            unnest(tidy), here("Output","Rates_controlonlymods.csv"))


r1|conc1
ggsave(here("Output","ControlOnlyEffects_sqrt.pdf"), width = 8, height = 8, device = cairo_pdf)

Benthic_all<-BenthicData %>%
  select(pool_id = PoolID, Foundation_spp:Before_After, MusselCover, SurfgrassCover) %>%
  rename(foundation_spp = Foundation_spp, before_after = Before_After, removal_control = Removal_Control)%>%
  mutate(pool_id = as.character(pool_id),
         MusselCover = ifelse(is.na(MusselCover),0,MusselCover),
         SurfgrassCover = ifelse(is.na(SurfgrassCover),0,SurfgrassCover))%>%
  pivot_wider(names_from = before_after,
              values_from = c(MusselCover,SurfgrassCover)) %>%
  mutate(mussel_change = MusselCover_After - MusselCover_Before, # difference in cover from removal
         surfgrass_change = SurfgrassCover_After - SurfgrassCover_Before) 
 
Benthic_long<-Benthic_all %>%
  filter(foundation_spp == "Phyllospadix",
         removal_control == "Control") %>%
  select(pool_id:removal_control, cover_change_adj = SurfgrassCover_Before)%>%
  mutate(before_after = "Before") %>%
  bind_rows(
    Benthic_all %>%
      filter(foundation_spp == "Phyllospadix",
             removal_control == "Control") %>%
      select(pool_id:removal_control, cover_change_adj = SurfgrassCover_After)%>%
      mutate(before_after = "After") 
  ) %>%
  bind_rows(
    Benthic_all %>%
      filter(foundation_spp == "Phyllospadix",
             removal_control == "Removal") %>%
      select(pool_id:removal_control, cover_change_adj = SurfgrassCover_Before)%>%
      mutate(before_after = "Before")
  ) %>%
  bind_rows(
    Benthic_all %>%
      filter(foundation_spp == "Phyllospadix",
             removal_control == "Removal") %>%
      select(pool_id:removal_control, cover_change_adj = surfgrass_change)%>%
      mutate(before_after = "After")
  ) %>%
  
  bind_rows(Benthic_all %>%
            filter(foundation_spp == "Mytilus",
                   removal_control == "Control") %>%
            select(pool_id:removal_control, cover_change_adj = MusselCover_Before)%>%
            mutate(before_after = "Before")) %>%
        bind_rows(
              Benthic_all %>%
                filter(foundation_spp == "Mytilus",
                       removal_control == "Control") %>%
                select(pool_id:removal_control, cover_change_adj = MusselCover_After)%>%
                mutate(before_after = "After") 
            ) %>%
            bind_rows(
              Benthic_all %>%
                filter(foundation_spp == "Mytilus",
                       removal_control == "Removal") %>%
                select(pool_id:removal_control, cover_change_adj = MusselCover_Before)%>%
                mutate(before_after = "Before")
            ) %>%
            bind_rows(
              Benthic_all %>%
                filter(foundation_spp == "Mytilus",
                       removal_control == "Removal") %>%
                select(pool_id:removal_control, cover_change_adj = mussel_change)%>%
                mutate(before_after = "After")
            )

# add in the diatoms
Benthic_long <-Benthic_long %>% left_join(
BenthicData%>% 
  select(pool_id = PoolID, Foundation_spp:Before_After, Diatoms) %>%
  mutate(pool_id = as.character(pool_id))%>%
  rename(foundation_spp = Foundation_spp, before_after = Before_After, removal_control = Removal_Control)
)  

Rates %>% 
  filter(day_night == "Day",
         foundation_spp != "Ocean",
         foundation_spp == "Mytilus",
         name %in% c("bix","do_mg_l","heterotrophic_bacterioplankton_m_l","hix", "m_c", "nh4_umol_l", "nn_umol_l","total_fDOM"))%>%
  filter(name != "total_fDOM" | rate_m2_hr <0.1)%>%
  mutate(together = paste(before_after, removal_control),
         manipulated = ifelse(together == "After Removal","Manipulated", "Not Manipulated"))%>%
  left_join(Benthic_long) %>%
  ggplot(aes(x = cover_change_adj, y = rate_hr, color = before_after))+
  geom_point()+
  facet_wrap(name~manipulated, scale = "free")

values %>% 
  mutate(before_after = ifelse(month == "July","Before","After"))%>%
  left_join(Benthic_long) %>%
  filter(#day_night == "Day",
         foundation_spp != "Ocean",
         foundation_spp == "Phyllospadix",
         name %in% c("bix","do_mg_l","heterotrophic_bacterioplankton_m_l","hix", "m_c", "nh4_umol_l", "nn_umol_l","total_fDOM"))%>%
    mutate(together = paste(before_after, removal_control),
         manipulated = ifelse(together == "After Removal","Manipulated", "Not Manipulated"))%>%
  ggplot(aes(x =  cover_change_adj , y = log(mean_val)))+
  geom_point(aes( color = month))+
  geom_smooth(method = "lm")+
 # coord_trans(x = "log")+
  facet_wrap(name~manipulated, scale = "free")

Rates %>% 
  filter(name %in% c("m_c","bix", "nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "hix") )%>%
  left_join(Benthic_long) %>%
  filter(foundation_spp == "Phyllospadix")%>%
  mutate(together = paste(before_after, removal_control),
         manipulated = ifelse(together == "After Removal","Manipulated", "Not Manipulated"))%>%
  ggplot(aes(x =  cover_change_adj , y = rate_m2_hr))+
  geom_point(aes( color = before_after))+
  geom_smooth(method = "lm")+
  # coord_trans(x = "log")+
  facet_wrap(name~manipulated, scale = "free")


data_all %>%
  ungroup()%>%
  filter(day_night == "Day",
         #foundation_spp != "Ocean"
  )  %>%
  mutate(prot = tyrosine_like+tryptophan_like+ phenylalanine_like,
         humic = ultra_violet_humic_like+visible_humic_like+marine_humic_like,
         total_fDOM = humic+prot) %>%
  select(before_after, month,pool_id, day_night, time_point,removal_control, foundation_spp,do_mg_l,heterotrophic_bacterioplankton_m_l:autotrophic_pico_eukaryotes_m_l,prot,humic,m_c, bix, hix,fi, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_mg_l:nh4_umol_l) %>%
  group_by(foundation_spp, pool_id,removal_control,name, before_after) %>%
  summarise(mean_val = mean(value, na.rm = TRUE))%>%
  mutate(mean_val_sqrt = sqrt(mean_val))%>%
  filter(name %in% c("nn_umol_l","nh4_umol_l") )%>%
  mutate(together = paste(before_after, removal_control),
         manipulated = ifelse(together == "After Removal","Manipulated", "Not Manipulated")) %>%
  group_by(name, foundation_spp, manipulated, before_after) %>%
  summarise(means = mean(mean_val, na.rm = TRUE))
  

### test cover as a covariate
Long_all<-Rates %>% 
  filter(day_night == "Day",
         foundation_spp != "Ocean",
      #   foundation_spp == "Mytilus",
         name %in% c("bix","do_mg_l","heterotrophic_bacterioplankton_m_l","hix", "m_c", "nh4_umol_l", "nn_umol_l","total_fDOM"))%>%
  filter(name != "total_fDOM" | rate_m2_hr <0.1)%>%
  mutate(together = paste(before_after, removal_control),
         manipulated = ifelse(together == "After Removal","Manipulated", "Not Manipulated"))%>%
  left_join(Benthic_long) %>%
  mutate(rate_sqrt = sign(rate_m2_hr)*sqrt(abs(rate_m2_hr)),
         rate_log = sign(rate_m2_hr)*log(abs(rate_m2_hr)))

# make a square root function with negatives
sqrt_fcn<-function(x){sign(x)*sqrt(abs(x))}
# transform the axes with this function
tn <- trans_new("sqrt_fcn",
                function(x){sign(x)*sqrt(abs(x))},
                function(y){sign(y)*abs(y)^2}
               )
# get the average values for the plot
sum_rates<-Long_all %>%
  filter(name == "nn_umol_l",
         foundation_spp == "Mytilus") %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)"))) %>%
  group_by(month, removal_control) %>%
  summarise(mean_rate = mean(rate_m2_hr, na.rm = TRUE),
            se_rate = sd(rate_m2_hr, na.rm = TRUE)/sqrt(n()))

NN_mussel<-Long_all %>%
  filter(name == "nn_umol_l",
         foundation_spp == "Mytilus")  %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)"))) %>%
  ggplot()+
  geom_point(aes(x = month, y = rate_m2_hr, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_path(aes(x = month, y = rate_m2_hr, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_point(data = sum_rates,aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 3 )+
  geom_errorbar(data = sum_rates,aes(x = month, ymin = mean_rate - se_rate, ymax = mean_rate+se_rate, color = removal_control), width = 0.01, size = 1)+
  geom_path(data = sum_rates,aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 1 )+
  scale_y_continuous(transform = tn, 
                     breaks = c(-0.1, 0, 0.1, 0.25, 0.5))+
  labs(x = "",
       y = "&Delta;Nitrate+Nitrite <br>(&mu;mol m<sup>-2</sup> hr<sup>-1</sup>)",
       color = "")+
  scale_color_manual(values = c("#e02b35","#082a54"))+
  theme_bw()+
  theme(axis.title = element_markdown(size = 14))

sum_rates<-Long_all %>%
  filter(name == "heterotrophic_bacterioplankton_m_l",
         foundation_spp == "Phyllospadix") %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)"))) %>%
  group_by(month, removal_control) %>%
  summarise(mean_rate = mean(rate_m2_hr, na.rm = TRUE),
            se_rate = sd(rate_m2_hr, na.rm = TRUE)/sqrt(n()))

Het_Phylo<-Long_all %>%
  filter(name == "heterotrophic_bacterioplankton_m_l",
         foundation_spp == "Phyllospadix")  %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)"))) %>%
  ggplot()+
  geom_point(aes(x = month, y = rate_m2_hr, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_path(aes(x = month, y = rate_m2_hr, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_point(data = sum_rates,aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 3 )+
  geom_errorbar(data = sum_rates,aes(x = month, ymin = mean_rate - se_rate, ymax = mean_rate+se_rate, color = removal_control), width = 0.01, size = 1)+
  geom_path(data = sum_rates,aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 1 )+
  scale_y_continuous(transform = tn, breaks = c(-5000, -1000, 0,1000,5000, 10000, 20000))+
  labs(x = "",
       y = "&Delta;Heterotrophic Bacteria <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
       color = "")+
  scale_color_manual(values = c("#e02b35","#082a54"))+
  theme_bw()+
  theme(axis.title = element_markdown(size = 14))

# Same with significant values
#NH4
sum_rates<-values %>%
  filter(name%in% c("nh4_umol_l", "nn_umol_l","bix"),
         foundation_spp == "Mytilus") %>%
  group_by(name, month, removal_control) %>%
  summarise(mean_rate = mean(mean_val, na.rm = TRUE),
            se_rate = sd(mean_val, na.rm = TRUE)/sqrt(n()))

NH4_mussel<-values %>%
  filter(name == "nh4_umol_l",
         foundation_spp == "Mytilus")  %>%
  ggplot()+
  geom_point(aes(x = month, y = mean_val, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_path(aes(x = month, y = mean_val, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_point(data = sum_rates %>% filter(name == "nh4_umol_l"),aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 3 )+
  geom_errorbar(data = sum_rates %>% filter(name == "nh4_umol_l"),
                aes(x = month, ymin = mean_rate - se_rate, ymax = mean_rate+se_rate, color = removal_control), 
                width = 0.01, size = 1)+
  geom_path(data = sum_rates %>% filter(name == "nh4_umol_l"),
            aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 1 )+
  scale_y_continuous(transform = tn, breaks = c(2, 10, 20, 30))+
  scale_color_manual(values = c("#e02b35","#082a54"))+
  labs(x = "",
       y = "Ammonium <br>(&mu;mol L<sup>-1</sup>)",
       color = "")+
  theme_bw()+
  theme(axis.title = element_markdown(size = 14))

NN_mussel_v<-values %>%
  filter(name == "nn_umol_l",
         foundation_spp == "Mytilus")  %>%
  ggplot()+
  geom_point(aes(x = month, y = mean_val, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_path(aes(x = month, y = mean_val, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_point(data = sum_rates %>% filter(name == "nn_umol_l"),aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 3 )+
  geom_errorbar(data = sum_rates %>% filter(name == "nn_umol_l"),
                aes(x = month, ymin = mean_rate - se_rate, ymax = mean_rate+se_rate, color = removal_control), 
                width = 0.01, size = 1)+
  geom_path(data = sum_rates %>% filter(name == "nn_umol_l"),
            aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 1 )+
  scale_y_continuous(transform = tn, breaks = c(2, 10, 20, 30))+
  scale_color_manual(values = c("#e02b35","#082a54"))+
  labs(x = "",
       y = "Nitrate + Nitrite <br>(&mu;mol L<sup>-1</sup>)",
       color = "")+
  theme_bw()+
  theme(axis.title = element_markdown(size = 14))


BIX_mussel_v<-values %>%
  filter(name == "bix",
         foundation_spp == "Mytilus")  %>%
  ggplot()+
  geom_point(aes(x = month, y = mean_val, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_path(aes(x = month, y = mean_val, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_point(data = sum_rates %>% filter(name == "bix"),aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 3 )+
  geom_errorbar(data = sum_rates %>% filter(name == "bix"),
                aes(x = month, ymin = mean_rate - se_rate, ymax = mean_rate+se_rate, color = removal_control), 
                width = 0.01, size = 1)+
  geom_path(data = sum_rates %>% filter(name == "bix"),
            aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 1 )+
  scale_y_continuous(transform = tn)+
  scale_color_manual(values = c("#e02b35","#082a54"))+
  labs(x = "",
       y = "BIX",
       color = "")+
  theme_bw()+
  theme(axis.title = element_markdown(size = 14))

sum_rates<-values %>%
  filter(name%in% c("heterotrophic_bacterioplankton_m_l"),
         foundation_spp == "Phyllospadix") %>%
  group_by(name, month, removal_control) %>%
  summarise(mean_rate = mean(mean_val, na.rm = TRUE),
            se_rate = sd(mean_val, na.rm = TRUE)/sqrt(n()))

Het_Phylo_v<-values %>%
  filter(name == "heterotrophic_bacterioplankton_m_l",
         foundation_spp == "Phyllospadix")  %>% 
  ggplot()+
  geom_point(aes(x = month, y = mean_val, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_path(aes(x = month, y = mean_val, color = removal_control, group = pool_id), alpha = 0.2)+
  geom_point(data = sum_rates,aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 3 )+
  geom_errorbar(data = sum_rates,aes(x = month, ymin = mean_rate - se_rate, ymax = mean_rate+se_rate, color = removal_control), width = 0.01, size = 1)+
  geom_path(data = sum_rates,aes(x = month, y = mean_rate, color = removal_control, group = removal_control), size = 1 )+
  scale_y_continuous(transform = tn)+
  scale_color_manual(values = c("#e02b35","#082a54"))+
  labs(x = "",
       y = "Heterotrophic Bacteria <br>(# mL<sup>-1</sup>)",
       color = "")+
  theme_bw()+
  theme(axis.title = element_markdown(size = 14))

BC_int /((NN_mussel/ plot_spacer()/ plot_spacer()/Het_Phylo) | (NN_mussel_v/NH4_mussel/BIX_mussel_v/Het_Phylo_v))+plot_layout(guides = "collect",widths = c(1,1,1), heights = c(2, 4,4))

ggsave(here("Output","Composite_BACI.pdf"), width = 8, 
       height = 16, device = cairo_pdf)


#### Make same plot but with means instead of reaction norms

summary(lmer(rate_sqrt ~ removal_control*before_after + cover_change_adj +(1|pool_id), Long_all %>%
               filter(name == "nn_umol_l",
                      foundation_spp == "Mytilus")))

NN_rate_2<-Long_all %>%
  filter(name == "nn_umol_l",
         foundation_spp == "Mytilus") %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)"))) %>%
  group_by(pool_id, removal_control)%>%
  reframe(difference = rate_m2_hr[before_after == "After"] - rate_m2_hr[before_after == "Before"]) %>%
  ggplot(aes(x = removal_control, y = difference))+
  geom_point(alpha = 0.2)+
  stat_summary(size = 1)+
  labs(x = " ",
       y = "&Delta;Nitrate+Nitrite <br>(&mu;mol m<sup>-2</sup> hr<sup>-1</sup>)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_markdown(size = 14))

HBac_rate_2<-Long_all %>%
  filter(name == "heterotrophic_bacterioplankton_m_l",
         foundation_spp == "Phyllospadix") %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)"))) %>%
  group_by(pool_id, removal_control)%>%
  reframe(difference = rate_m2_hr[before_after == "After"] - rate_m2_hr[before_after == "Before"]) %>%
  ggplot(aes(x = removal_control, y = difference))+
  geom_point(alpha = 0.2, color = "#34c230")+
  stat_summary(size = 1, color = "#34c230")+
  labs(x = " ",
       y = "&Delta;Heterotrophic Bacteria <br> (# m<sup>-2</sup> hr<sup>-1</sup>)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_markdown(size = 14))

NN_V2<-values %>%
  filter(name == "nn_umol_l",
         foundation_spp == "Mytilus")%>%
  group_by(pool_id, removal_control)%>%
  reframe(difference = mean_val[month == "August (Upwelling)"] - mean_val[month == "July"]) %>%
  filter(difference <15)%>%
  ggplot(aes(x = removal_control, y = difference))+
  geom_point(alpha = 0.2)+
  stat_summary(size = 1)+
  labs(x = " ",
       y = "&Delta;N+N </sup> <br>(&mu;mol L<sup>-1</sup>)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_markdown(size = 14))

NH4_V2<-values %>%
  filter(name == "nh4_umol_l",
         foundation_spp == "Mytilus")%>%
  group_by(pool_id, removal_control)%>%
  reframe(difference = mean_val[month == "August (Upwelling)"] - mean_val[month == "July"]) %>%
  filter(difference > -20)%>%
  ggplot(aes(x = removal_control, y = difference))+
  geom_point(alpha = 0.2)+
  stat_summary(size = 1)+
  labs(x = " ",
       y = "&Delta;NH<sub>4</sub><sup>+</sup> <br>(&mu;mol L<sup>-1</sup>)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_markdown(size = 14))

BIX_V2<-values %>%
  filter(name == "bix",
         foundation_spp == "Mytilus")%>%
  group_by(pool_id, removal_control)%>%
  reframe(difference = mean_val[month == "August (Upwelling)"] - mean_val[month == "July"]) %>%
  filter(difference > -20)%>%
  ggplot(aes(x = removal_control, y = difference))+
  geom_point(alpha = 0.2)+
  stat_summary(size = 1)+
  labs(x = " ",
       y = "&Delta;BIX")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_markdown(size = 14))

Hbac_V2<-values %>%
  filter(name == "heterotrophic_bacterioplankton_m_l",
         foundation_spp == "Phyllospadix")%>%
  group_by(pool_id, removal_control)%>%
  reframe(difference = mean_val[month == "August (Upwelling)"] - mean_val[month == "July"]) %>%
  ggplot(aes(x = removal_control, y = difference))+
  geom_point(alpha = 0.2, color = "#34c230")+
  stat_summary(size = 1, color = "#34c230")+
   labs(x = " ",
       y = "&Delta;Heterotrophic Bacteria <br> (# mL<sup>-1</sup>)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_markdown(size = 14))

(BC_int&theme(strip.text = element_blank()))/((NN_rate_2/NN_V2/HBac_rate_2)|(NH4_V2/BIX_V2/Hbac_V2))+plot_layout(guides = "collect",widths = c(1,1,1), heights = c(3, 5,5))

ggsave(here("Output","Composite_BACI_means.pdf"), width = 7, 
       height = 13, device = cairo_pdf)

# Plot showing relationship between DO and HBac

Rates_wide<- Rates %>%  
  filter(day_night == "Day") %>%
  filter(foundation_spp != "Ocean") %>%
  select(-c(change:vol, rate_hr))%>%
  pivot_wider(values_from = rate_m2_hr, names_from = name) %>%  
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)")))


Rates_wide %>%
  mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                             removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                             removal_control == "Removal"& month != "July" ~ "Foundation species removed"))%>%
  mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation species removed")))%>%
  mutate(foundation_spp = ifelse(foundation_spp == "Mytilus", "Mussel-Dominated","Surfgrass-Dominated"))%>%
  ggplot(aes(y = heterotrophic_bacterioplankton_m_l, x = do_mg_l))+
  geom_hline(aes(yintercept  = 0), lty = 2)+
  geom_vline(aes(xintercept  = 0), lty = 2)+
  geom_point(aes(color = foundation_spp, shape = month))+
  # geom_text(aes(x = 175, y = 0.1, label = "p = 0.036"))+
  #  ylim(-.3,.3)+
  geom_smooth(method ="lm", color = "grey3", data = Rates_wide %>%  
                mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                                           removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                                           removal_control == "Removal"& month != "July" ~ "Foundation species removed"))%>%
                mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation species removed")))%>%
                filter(foundation_spp != "Ocean",
                       removal == "Unmanipulated"))+
  scale_color_manual(values = c("black","#34c230"))+
  scale_shape_manual(values = c(1,16))+
  labs(y = "Heterotrophic bacterial production <br> (# m<sup>-2</sup> hr<sup>-1</sup>)",
       x = "DO production <br> (mg m<sup>-2</sup> hr<sup>-1</sup>)",
       shape = "Sampling Month",
       color = " Foundation Species")+
  facet_wrap(~removal, scales = "free_x")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_markdown(size = 14),
        axis.title.y = element_markdown(size = 14),
        axis.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
  )
ggsave(here("Output","het_DO_regression.pdf"), width = 8, height = 4)


DO_het_mod<-lm(data = Rates_wide %>%  
                   mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                                              removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                                              removal_control == "Removal"& month != "July" ~ "Foundation spp. removed"))%>%
                   mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation spp. removed")))%>%
                   filter(foundation_spp != "Ocean"), heterotrophic_bacterioplankton_m_l~do_mg_l*removal)


anova(DO_het_mod)
summary(DO_het_mod)

# get the ocean temperatures
data_all %>%
  filter(day_night == "Day",
         removal_control == "Ocean") %>%
  select(before_after, temp_pool) %>%
  group_by(before_after)%>%
  summarise(mean_temp = mean(temp_pool, na.rm = TRUE),
            se_temp = sd(temp_pool)/sqrt(n()))

data_all %>%
  filter(day_night == "Day") %>%
  select(foundation_spp, before_after,removal_control, nh4_umol_l) %>%
  group_by(foundation_spp,before_after, removal_control)%>%
  summarise(mean_temp = mean(nh4_umol_l, na.rm = TRUE),
            se_temp = sd(nh4_umol_l)/sqrt(n()))


## explore this and run one-way t-tests
Long_all %>% 
  filter(name != "total_fDOM")%>%
  mutate(man = ifelse(manipulated == "Manipulated", "Foundation spp removed", "Unmanipulated"),
         foundation_spp = ifelse(foundation_spp == "Mytilus", "Mussels", "Surfgrass")) %>%
  mutate(nicenames = case_when(
    name == "do_mg_l" ~ "&Delta;DO <br> (mg)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria <br> (counts)",
    name == "nh4_umol_l" ~ "&Delta; Ammonium <br> (&mu;mol)",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite <br> (&mu;mol)",
    name == "bix"~"&Delta;BIX" ,
    name == "hix"~"&Delta;HIX",
    name == "m_c"~"&Delta;M:C"))%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;DO <br> (mg)",
                                                  "&Delta; Ammonium <br> (&mu;mol)",
                                                  "&Delta;Nitrate+Nitrite <br> (&mu;mol)",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic Bacteria <br> (counts)",
                                                  "&Delta;fDOM"
  )))%>%
  ggplot(aes(x= rate_m2_hr, y = before_after, shape = man))+
  geom_vline(xintercept = 0, color = "grey")+
  geom_point(alpha = 0.1)+
  stat_summary(size = 0.8)+
  scale_shape_manual(values = c(21,19))+
  labs(x = "Flux (m<sup>-2</sup> hr<sup>-1</sup>)", 
       y = "",
       shape = "")+
  ggh4x::facet_grid2(nicenames~foundation_spp, scales = "free_x", independent = "x")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside", 
        panel.grid.minor = element_blank(),
        strip.text = element_markdown(),
        axis.title.x = element_markdown(),
        legend.position = "bottom")
ggsave(filename = here("Output","RawFlux.png"), height = 12, width = 6)


## same for the values

# get the ocean values
ocean <-data_all %>%
  ungroup()%>%
  filter(day_night == "Day",
         removal_control == "Ocean") %>%
  select(before_after, do_mg_l, heterotrophic_bacterioplankton_m_l, nn_umol_l, nh4_umol_l, bix, hix, m_c) %>%
  group_by(before_after)%>%
  summarise_all(.funs = function(x){mean(x, na.rm = TRUE)}) %>%
  pivot_longer(cols = do_mg_l:m_c) %>%
  mutate(nicenames = case_when(
    name == "do_mg_l" ~ "&Delta;DO <br> (mg)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria <br> (counts)",
    name == "nh4_umol_l" ~ "&Delta; Ammonium <br> (&mu;mol)",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite <br> (&mu;mol)",
    name == "bix"~"&Delta;BIX" ,
    name == "hix"~"&Delta;HIX",
    name == "m_c"~"&Delta;M:C"))%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;DO <br> (mg)",
                                                  "&Delta; Ammonium <br> (&mu;mol)",
                                                  "&Delta;Nitrate+Nitrite <br> (&mu;mol)",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic Bacteria <br> (counts)",
                                                  "&Delta;fDOM"
  ))) %>%
  mutate(removal = "Ocean")


values %>%
  mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                             removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                             removal_control == "Removal"& month != "July" ~ "Foundation species removed"))%>%
  mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation species removed")))%>%
  filter(name != "total_fDOM")%>%
  mutate(foundation_spp = ifelse(foundation_spp == "Mytilus", "Mussels", "Surfgrass")) %>%
  mutate(before_after = ifelse(month == "July","Before","After"))%>%
  mutate(nicenames = case_when(
    name == "do_mg_l" ~ "&Delta;DO <br> (mg)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria <br> (counts)",
    name == "nh4_umol_l" ~ "&Delta; Ammonium <br> (&mu;mol)",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite <br> (&mu;mol)",
    name == "bix"~"&Delta;BIX" ,
    name == "hix"~"&Delta;HIX",
    name == "m_c"~"&Delta;M:C"))%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;DO <br> (mg)",
                                                  "&Delta; Ammonium <br> (&mu;mol)",
                                                  "&Delta;Nitrate+Nitrite <br> (&mu;mol)",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic Bacteria <br> (counts)",
                                                  "&Delta;fDOM"
  )))%>%
  ggplot(aes(x= mean_val, y = before_after, shape = removal))+
  #geom_vline(xintercept = 0, color = "grey")+
  geom_point(alpha = 0.1)+
  stat_summary(size = 0.8)+
  geom_point(data = ocean, aes(x= value, y = before_after), 
             color = "lightblue", size = 3)+
  scale_shape_manual(values = c(19,21,19))+
  labs(x = "Stock", 
       y = "",
       shape = "")+
  ggh4x::facet_grid2(nicenames~foundation_spp, scales = "free_x", independent = "x")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside", 
        panel.grid.minor = element_blank(),
        strip.text = element_markdown(),
        axis.title.x = element_markdown(),
        legend.position = "bottom")
ggsave(filename = here("Output","Rawvalue.png"), height = 12, width = 6)
