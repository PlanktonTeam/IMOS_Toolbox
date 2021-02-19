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
NRSTrips <- get_NRSTrips()

# you will get a warning about the fast method, this actually works better than the accurate method for this data set. 

# Hydrochemistry data 
Chemistry <- getChemistry()

# Zooplankton biomass
ZBiomass <-  read_csv(paste0(rawD,.Platform$file.sep,"nrs_biomass.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_CODE, Biomass_mgm3 = BIOMASS_MGM3) %>% 
  mutate(SampleDepth_m = "WC",
         NRScode = gsub('^.{3}|.{9}$', '', NRScode)) %>% 
  untibble()

# Pigments data
Pigments <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_pigments.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_TRIP_CODE,
         SampleDepth_m = SAMPLE_DEPTH_M) %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         NRScode = substring(NRScode,4)) %>% 
  filter(QC_FLAG %in% c(0,1,2,5,8)) %>% # keep data flagged as good
  select(-QC_FLAG) %>%
  untibble()

# Flow cytometry picoplankton data
Pico <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_picoplankton.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_TRIP_CODE,
         SampleDepth_m = SAMPLE_DEPTH_M, Prochlorococcus_cells_ml = PROCHLOROCOCCUS_CELLSPERML, Synecochoccus_cells_ml = SYNECOCHOCCUS_CELLSPERML, 
         Picoeukaryotes_cells_ml = PICOEUKARYOTES_CELLSPERML) %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         NRScode = substring(NRScode,4),
         Prochlorococcus_cells_ml = ifelse(PROCHLOROCOCCUS_FLAG %in% c(3,4,9), NA, Prochlorococcus_cells_ml), # remove bad data
         Synecochoccus_cells_ml = ifelse(SYNECOCHOCCUS_FLAG %in% c(3,4,9), NA, Synecochoccus_cells_ml),
         Picoeukaryotes_cells_ml = ifelse(PICOEUKARYOTES_FLAG %in% c(3,4,9), NA, Picoeukaryotes_cells_ml)) %>%
  group_by(NRScode, SampleDepth_m) %>% 
  summarise(Prochlorococcus_cells_ml = mean(Prochlorococcus_cells_ml, na.rm = TRUE), # mean of replicates
            Synecochoccus_cells_ml = mean(Synecochoccus_cells_ml, na.rm = TRUE),
            Picoeukaryotes_cells_ml = mean(Picoeukaryotes_cells_ml, na.rm = TRUE),
            .groups = "drop") %>% 
  untibble()

# Total suspended solid data
TSS <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_TSS.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_TRIP_CODE, SampleDepth_m = SAMPLE_DEPTH_M, TSS_mg_L = TSS_MG_PER_L, 
         InorganicFraction_mg_L = INORGANIC_FRACTION_MG_PER_L, 
         OrganicFraction_mg_L = ORGANIC_FRACTION_MG_PER_L, Secchi_m = SECCHI_DEPTH_M) %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         NRScode = substring(NRScode,4),
         TSS_mg_L = ifelse(TSS_FLAG %in% c(3,4,9), NA, TSS_mg_L), # remove bad data
         InorganicFraction_mg_L = ifelse(TSS_FLAG %in% c(3,4,9), NA, InorganicFraction_mg_L),
         OrganicFraction_mg_L = ifelse(TSS_FLAG %in% c(3,4,9), NA, OrganicFraction_mg_L)) %>%
  group_by(NRScode, SampleDepth_m) %>% 
  summarise(TSS_mg_L = mean(TSS_mg_L, na.rm = TRUE), # mean of replicates
            InorganicFraction_mg_L = mean(InorganicFraction_mg_L, na.rm = TRUE),
            OrganicFraction_mg_L = mean(OrganicFraction_mg_L, na.rm = TRUE),
            .groups = "drop") %>% 
  untibble()

# Secchi Disc        
Secchi <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_TSS.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_TRIP_CODE, SampleDepth_m = SAMPLE_DEPTH_M, Secchi_m = SECCHI_DEPTH_M) %>%
  select(NRScode, Secchi_m, SampleDepth_m) %>% 
  distinct() %>%
  mutate(SampleDepth_m = "WC",
         NRScode = substring(NRScode,4))

# CTD Cast Data
CTD <- getCTD() %>%
    mutate(SampleDepth_m = as.character(round(SampleDepth_m, 0))) %>% 
    select(-c(PRES, CTDSpecificConductivity_Sm)) %>%
    group_by(NRScode, SampleDepth_m) %>% summarise(CTDDensity_kgm3 = mean(CTDDensity_kgm3, na.rm = TRUE),
                                                   CTDTemperature = mean(CTDTemperature, na.rm = TRUE),
                                                   CTDPAR_umolm2s = mean(CTDPAR_umolm2s, na.rm = TRUE),
                                                   CTDConductivity_sm = mean(CTDConductivity_sm, na.rm = TRUE),
                                                   CTDSalinity = mean(CTDSalinity, na.rm = TRUE),
                                                   CTDChlF_mgm3 = mean(CTDChlF_mgm3, na.rm = TRUE),
                                                   CTDTurbidity_ntu = mean(CTDTurbidity_ntu, na.rm = TRUE)) %>%
    untibble()

notrips <-  read_csv(paste0(rawD,.Platform$file.sep,"nrs_CTD.csv"), na = "(null)",
                     col_types = cols(PRES = col_double(), # columns start with nulls so tidyverse annoyingly assigns col_logical()
                                      PAR = col_double(),
                                      SPEC_CNDC = col_double())) %>% select(NRS_TRIP_CODE) %>% distinct()


# Combined BGC data for each station at the sample depth
BGC <- NRSTrips %>% mutate(IMOSsampleCode = paste0('NRS',NRScode, '_', ifelse(SampleDepth_m == 'WC', 'WC', str_pad(SampleDepth_m, 3, side = "left", "0")))) %>%
  left_join(ZBiomass %>% 
              select(NRScode, SampleDepth_m, Biomass_mgm3), by = c("NRScode", "SampleDepth_m")) %>%
  left_join(Secchi,  by = c("NRScode", "SampleDepth_m")) %>%
  left_join(Chemistry, by = c("NRScode", "SampleDepth_m")) %>%
  left_join(Pico, by = c("NRScode", "SampleDepth_m")) %>%
  left_join(Pigments, by = c("NRScode", "SampleDepth_m")) %>%
  left_join(TSS, by = c("NRScode", "SampleDepth_m")) %>%
  left_join(CTD, by = c("NRScode", "SampleDepth_m")) 

# test table
# n should be 1, replicates or duplicate samples will have values > 1
test <- BGC %>% 
  group_by(NRScode, SampleDepth_m) %>% 
  summarise(n = n(),
            .groups = "drop")

# Check
max(test$n)

# save to github
fwrite(BGC, file = paste0(outD,.Platform$file.sep,"NRS_CombinedWaterQuality.csv"), row.names = FALSE)
