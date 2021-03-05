## IMOS plankton data products Indices 
## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)

## Created: Sept 2020
## Updated: 
## 1 Oct 2020 (Written to Git)
## 6th October 2020

suppressPackageStartupMessages({
  library(lubridate)
  library(lutz) # for the timezone calculations #remotes::install_github("jiho/castr")
  library(castr)
  library(vegan)
  library(data.table)  
  library(ncdf4) # devtools::install_github("mdsumner/ncdf4")
  library(tidyverse)
})

source("IMOS_Plankton_functions.R")
# source("../Satellite/fIMOS_MatchAltimetry.R")
# source("../Satellite/fIMOS_MatchMODIS.R")
# source("../Satellite/fIMOS_MatchGHRSST.R")

rawD <- "RawData"
outD <- "Output"

# uses mostly the same raw data from IMOS_PlanktonProducts_Create.R

# ensure we have all trips accounted for 
# note there are circumstances where a trip won't have a phyto and a zoo samples due to loss of sample etc.

NRSdat <- get_NRSTrips() %>% 
  select(-SampleDepth_m) %>% 
  distinct() # warning message on 'fast' versus 'accurate' method, in this case fast is better.

dNRSdat <- distinct(NRSdat, NRScode, .keep_all = TRUE) %>%  # Distinct rows for satellite, should be anyway
  rename(Date = SampleDateLocal) %>% 
  select(NRScode, Date, Latitude, Longitude)

# SST and Chlorophyll from CTD
CTD <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_CTD.csv"), na = "(null)",
                col_types = cols(PRES = col_double(), # columns start with  nulls so tidyverse annoyingly assigns col_logical()
                                 PAR = col_double(),
                                 SPEC_CNDC = col_double())) %>% 
  rename(NRScode = NRS_TRIP_CODE, SampleDepth_m = PRES_REL, CTDDensity_kgm3 = DENS, 
         CTDTemperature = TEMP, CTDPAR_umolm2s = PAR,
         CTDConductivity_sm = CNDC, CTDSpecificConductivity_Sm = SPEC_CNDC, CTDSalinity = PSAL, 
         CTDTurbidity_ntu = TURB, CTDChlF_mgm3 = CHLF) %>%
  filter(SampleDepth_m < 15) %>% # take average of top 10m as a surface value for SST and CHL, this is removing 17 casts as of nov 2020
  group_by(NRScode) %>% 
  summarise(CTD_SST_C = mean(CTDTemperature, na.rm = TRUE),
            CTDChlF_mgm3 = mean(CTDChlF_mgm3, na.rm = TRUE),
            .groups = "drop") %>%
  untibble()

CTD_MLD <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_CTD.csv"), na = "(null)",
                    col_types = cols(PRES = col_double(), # columns start with  nulls so tidyverse annoyingly assigns col_logical()
                                     PAR = col_double(),
                                     SPEC_CNDC = col_double())) %>% 
  rename(NRScode = NRS_TRIP_CODE, SampleDepth_m = PRES_REL, CTDDensity_kgm3 = DENS, 
         CTDTemperature = TEMP, CTDPAR_umolm2s = PAR,
         CTDConductivity_sm = CNDC, CTDSpecificConductivity_Sm = SPEC_CNDC, CTDSalinity = PSAL, 
         CTDTurbidity_ntu = TURB, CTDChlF_mgm3 = CHLF) %>% 
  select(NRScode, CTDDensity_kgm3, CTDTemperature, CTDSalinity, SampleDepth_m)

n = nrow(CTD_MLD %>% select(NRScode) %>% unique())
MLD <- data.frame(NRScode = NA, MLD = NA, type = NA)

# MLD by density using castr
for (i in 1:n) {
  dat <- CTD_MLD %>% select(NRScode) %>% unique() %>% mutate(NRScode = as.factor(NRScode))
  nrscode <- dat$NRScode[[i]] %>% droplevels()
  mldData <- CTD_MLD %>% filter(NRScode == nrscode) %>% arrange(SampleDepth_m)
  mld_dens <- mld(mldData$CTDDensity_kgm3, mldData$SampleDepth_m, n.smooth = 3, ref.depths = 9:11)
  MLD <- rbind(MLD, c(as.character(nrscode), mld_dens, "Dens"))
}

# MLD by T and S  
for (i in 1:n) {
  dat <- CTD_MLD %>% select(NRScode) %>% unique() %>% mutate(NRScode = as.factor(NRScode))
  nrscode <- "YON20120816"
  mldData <- CTD_MLD %>% filter(NRScode == nrscode) %>% arrange(SampleDepth_m)
  ref_T <- mldData %>% mutate(refd = abs(SampleDepth_m - 10),
                             rankrefd = ave(refd, FUN = . %>% order %>% order)) %>%
    filter(rankrefd == 1)
  refT <- ref_T$CTDTemperature - 0.4 # temp at 10m minus 0.4
  mldData <- mldData %>% filter(SampleDepth_m >= ref_T$SampleDepth_m)
  mld_t <- mldData %>% mutate(temp = abs(CTDTemperature - refT),
                              ranktemp = ave(temp, FUN = . %>% order %>% order)) %>%
  filter(ranktemp == 1)
  MLD_temp <- mld_t$SampleDepth_m 
  
  refS <- ref_T$CTDSalinity - 0.03 # temp at 10m minus 0.4
  mld_s <- mldData %>% mutate(temp = abs(CTDSalinity - refS),
                              ranksal = ave(temp, FUN = . %>% order %>% order)) %>%
    filter(ranksal == 1)
  MLD_sal <- mld_s$SampleDepth_m
  MLD <- rbind(MLD, c(as.character(nrscode), MLD_temp, "Temp"))
  MLD <- rbind(MLD, c(as.character(nrscode), MLD_sal, "Sal"))           
}

MLD_ARR <- read_csv("MLD_NRS_ARR_PvR.csv", na = "(null)") %>%
  mutate(Station = recode("North Stradbroke" = "NSI", 
                          "Port Hacking" = "PHB", 
                          "Kangaroo Island" = "KAI",
                          "Maria Island" = "MAI", 
                          "Rottnest Island" = "ROT",
                          "Darwin" = "DAR",
                          "Yongala" = "YON", Site),
         Dates = as.numeric(gsub("-", "", Date)),
         NRScode = paste0(Station, Dates),
         type = "ARR",
         NRSmatch = substr(NRScode, 1,nchar(NRScode)-2)) %>%
  select(NRScode, MLD, type, NRSmatch)

MLDtest <-  MLD %>% mutate(NRSmatch = substr(NRScode, 1,nchar(NRScode)-2)) %>%
  rbind(MLD_ARR) %>%
  mutate(Station = substr(NRScode, 0,3),
         Date = ymd(substr(NRScode, 4, 11)),
         MLD = as.numeric(MLD, na.rm = TRUE)) %>%
  filter(!Station %in% c("NIN", "ESP")) %>%
  mutate(Station = factor(Station, levels = c("DAR","YON", "NSI", "ROT", "PHB", "KAI", "MAI"))) %>% drop_na()

#x11(width = 6, height = 8)
mldplot <- ggplot(MLDtest) + geom_point(aes(x=Date, y=MLD, color = type), size = 0.5) +
  geom_line(aes(x=Date, y=MLD, color = type)) +
  facet_grid(Station~., scales = "free") +
  labs(y = "Mid Layer Depth (m)") +
  theme_bw() + theme(strip.background = element_blank())
mldplot

ggsave("mld_nrs.png", mldplot, dpi = 600)

MLDmon <- MLDtest %>% mutate(Mon = month(Date)) %>%
  group_by(type, Station, Mon) %>% summarise(MLD = mean(MLD, na.rm = TRUE), .groups = "drop")

#x11(width = 6, height = 8)
mldsplot <- ggplot(MLDmon) + geom_point(aes(x=Mon, y=MLD, color = type), size = 0.5) +
  geom_line(aes(x=Mon, y=MLD, color = type)) +
  facet_grid(Station~., scales = "free") +
  labs(y = "Mid Layer Depth (m)") +
  scale_x_continuous(breaks = seq(1,12,length.out = 12), label = c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  theme_bw() + theme(strip.background = element_blank())
mldsplot
ggsave("mld_nrs_mon.png", mldsplot, dpi = 600)

# look at correlations
x11(width = 8, height = 8)
MLDcorr <- MLDtest %>% select(type, MLD, NRScode) %>%
  mutate(NRScode = substr(NRScode, 1,nchar(NRScode)-2)) %>%
  pivot_wider(values_from = MLD, names_from = type, values_fn = mean) %>%
  mutate(Station = substr(NRScode, 0,3)) %>%
  mutate(Station = factor(Station, levels = c("DAR","YON", "NSI", "ROT", "PHB", "KAI", "MAI"))) %>% drop_na()

pairs(MLDcorr[,c(2:5)]) 

plots <- list()
counter <- 1

for (i in 1:7){
  stat <- as.character(MLDcorr$Station[[i]] %>% droplevels())
  dat <- MLDcorr %>% filter(Station == stat)
  pairs(dat[,c(2:5)]) 
  plots[[counter]] <- p
  counter <- counter + 1
}

plotsout <- plots[1:6]
plotsout
