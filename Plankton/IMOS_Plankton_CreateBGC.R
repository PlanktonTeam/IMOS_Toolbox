## IMOS BGC Combined Water Quality Parameters
## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)

## Created: Aug 2020
## Updated: 
## 24 Sept 2020 (Written to Git)
## 6th October 2020

suppressPackageStartupMessages({
  library(lutz)
  library(purrr)
  library(lubridate)
  library(data.table)
  library(tidyverse)
})

source("IMOS_Plankton_functions.R")

rawD <- "RawData"
outD <- "Output"

################################
## Bring in data for combined water quality
################################

# Each trip and depth combination for water quality parameters
# the number of rows in this table should equal that in comb, if not look out for duplicates and replicates
NRSTrips <- getNRSTrips() %>% select(-SampleType)

# you will get a warning about the fast method, this actually works better than the accurate method for this data set. 

# Hydrochemistry data 
Chemistry <- getChemistry()

# Zooplankton biomass
ZBiomass <-  getNRSTrips() %>% select(TripCode, Biomass_mgm3, Secchi_m) %>%
  mutate(SampleDepth_m = 'WC') 

# Pigments data
Pigments <- read_csv(paste0(rawD,.Platform$file.sep,"BGC_Pigments.csv"), na = "(null)") %>% 
  rename(TripCode = TRIP_CODE,
         SampleDepth_m = SAMPLEDEPTH_M) %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m)) %>% 
  filter(PIGMENTS_FLAG %in% c(0,1,2,5,8)) %>% # keep data flagged as good
  select(-PIGMENTS_FLAG) %>%
  untibble()

# Flow cytometry picoplankton data
Pico <- read_csv(paste0(rawD,.Platform$file.sep,"BGC_Picoplankton.csv"), na = "(null)") %>% 
  rename(TripCode = TRIP_CODE,
         SampleDepth_m = SAMPLEDEPTH_M, Prochlorococcus_cellsml = PROCHLOROCOCCUS_CELLSML, Synecochoccus_cellsml = SYNECOCHOCCUS_CELLSML, 
         Picoeukaryotes_cellsml = PICOEUKARYOTES_CELLSML) %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         Prochlorococcus_cellsml = ifelse(PROCHLOROCOCCUS_FLAG %in% c(3,4,9), NA, Prochlorococcus_cellsml), # remove bad data
         Synecochoccus_cellsml = ifelse(SYNECOCHOCCUS_FLAG %in% c(3,4,9), NA, Synecochoccus_cellsml),
         Picoeukaryotes_cellsml = ifelse(PICOEUKARYOTES_FLAG %in% c(3,4,9), NA, Picoeukaryotes_cellsml)) %>%
  group_by(TripCode, SampleDepth_m) %>% 
  summarise(Prochlorococcus_cellsml = mean(Prochlorococcus_cellsml, na.rm = TRUE), # mean of replicates
            Synecochoccus_cellsml = mean(Synecochoccus_cellsml, na.rm = TRUE),
            Picoeukaryotes_cellsml = mean(Picoeukaryotes_cellsml, na.rm = TRUE),
            .groups = "drop") %>% 
  untibble()

# Total suspended solid data
TSS <- read_csv(paste0(rawD,.Platform$file.sep,"BGC_TSS.csv"), na = "(null)") %>% 
  rename(TripCode = TRIP_CODE, SampleDepth_m = SAMPLEDEPTH_M, TSS_mgL = TSS_MGL, 
         InorganicFraction_mgL = INORGANICFRACTION_MGL, 
         OrganicFraction_mgL = ORGANICFRACTION_MGL, Secchi_m = SECCHIDEPTH_M) %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         TripCode = substring(TripCode,4),
         TSS_mg_L = ifelse(TSSFLAG %in% c(3,4,9), NA, TSS_mgL), # remove bad data
         InorganicFraction_mgL = ifelse(TSSFLAG %in% c(3,4,9), NA, InorganicFraction_mgL),
         OrganicFraction_mgL = ifelse(TSSFLAG %in% c(3,4,9), NA, OrganicFraction_mgL)) %>%
  group_by(TripCode, SampleDepth_m) %>% 
  summarise(TSS_mgL = mean(TSS_mgL, na.rm = TRUE), # mean of replicates
            InorganicFraction_mgL = mean(InorganicFraction_mgL, na.rm = TRUE),
            OrganicFraction_mgL = mean(OrganicFraction_mgL, na.rm = TRUE),
            .groups = "drop") %>% 
  drop_na(SampleDepth_m) %>%
  untibble()

# Secchi Disc        
Secchi <- read_csv(paste0(rawD,.Platform$file.sep,"BGC_TSS.csv"), na = "(null)") %>% 
  rename(TripCode = TRIP_CODE, SampleDepth_m = SAMPLEDEPTH_M, Secchi_m = SECCHIDEPTH_M) %>%
  select(TripCode, Secchi_m, SampleDepth_m) %>% 
  distinct() %>%
  mutate(SampleDepth_m = "WC",
         TripCode = substring(TripCode,4))

# CTD Cast Data
CTD <- getCTD() %>%
    mutate(SampleDepth_m = as.character(round(Depth_m, 0))) %>% 
    select(-c(Pressure_dbar)) %>%
    group_by(TripCode, SampleDepth_m) %>% summarise(CTDDensity_kgm3 = mean(WaterDensity_kgm3, na.rm = TRUE),
                                                   CTDTemperature = mean(Temperature_degC, na.rm = TRUE),
                                                   CTDConductivity_sm = mean(Conductivity_Sm, na.rm = TRUE),
                                                   CTDSalinity = mean(Salinity_psu, na.rm = TRUE),
                                                   CTDChlF_mgm3 = mean(Chla_mgm3, na.rm = TRUE),
                                                   CTDTurbidity_ntu = mean(Turbidity_NTU, na.rm = TRUE)) %>%
    untibble()

notrips <-  read_csv(paste0(rawD,.Platform$file.sep,"nrs_CTD.csv"), na = "(null)",
                     col_types = cols(PRES = col_double(), # columns start with nulls so tidyverse annoyingly assigns col_logical()
                                      PAR = col_double(),
                                      SPEC_CNDC = col_double())) %>% select(NRS_TRIP_CODE) %>% distinct()

# combine for all samples taken
Samples <- rbind(Chemistry %>% select(TripCode, SampleDepth_m),
                 Pico %>% select(TripCode, SampleDepth_m),
                 Pigments %>% select(TripCode, SampleDepth_m),
                 TSS %>% select(TripCode, SampleDepth_m),
                 ZBiomass %>% select(TripCode, SampleDepth_m)) %>%
  unique()

# Combined BGC data for each station at the sample depth
BGC <- Samples %>% left_join(NRSTrips,  by = c("TripCode")) %>% mutate(IMOSsampleCode = paste0('NRS',TripCode, '_', ifelse(SampleDepth_m == 'WC', 'WC', str_pad(SampleDepth_m, 3, side = "left", "0")))) %>%
  left_join(Chemistry, by = c("TripCode", "SampleDepth_m")) %>%
  left_join(Pico, by = c("TripCode", "SampleDepth_m")) %>%
  left_join(Pigments, by = c("TripCode", "SampleDepth_m")) %>%
  left_join(TSS, by = c("TripCode", "SampleDepth_m")) %>%
  left_join(CTD, by = c("TripCode", "SampleDepth_m")) 

# test table
# n should be 1, replicates or duplicate samples will have values > 1
test <- BGC %>% 
  group_by(TripCode, SampleDepth_m) %>% 
  summarise(n = n(),
            .groups = "drop")

# Check
max(test$n)

# save to github
fwrite(BGC, file = paste0(outD,.Platform$file.sep,"NRS_CombinedWaterQuality.csv"), row.names = FALSE)
