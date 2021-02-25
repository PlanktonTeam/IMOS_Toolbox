## Functions for operating
untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense

## Functions for bringing in data sets
## NRS 
get_NRSTrips <- function(){
  NRSTrips <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_trips.csv"), na = "(null)") %>% 
    rename(Station = STATION_NAME, Latitude = Y_COORD, Longitude = X_COORD, SampleDateLocal = TRIP_START_DATETIME_LOCAL, 
           NRScode = NRS_CODE, StationDepth_m = STATION_DEPTH, SampleDepth_m = SAMPLEDEPTH_M, stateCode = STATE_CODE,
           DaylightSavings = DAYLIGHT_SAVINGS_Y_N, UTCoffsetH = UTC_OFFSET_H) %>%
    # select(-SAMPLEDEPTH_M) %>%
    mutate(Year = year(SampleDateLocal),
           Month = month(SampleDateLocal),
           Day = day(SampleDateLocal),
           Time_24hr = str_sub(SampleDateLocal, -8, -1), # hms doesn"t seem to work on 00:00:00 times
           tz = tz_lookup_coords(Latitude, Longitude, method = "fast"),
           SampleDateUTC = with_tz(force_tzs(SampleDateLocal, tz, roll = TRUE), "UTC"),
           SampleDateLocal = as.character(SampleDateLocal)) %>% 
    select(-c(tz, DaylightSavings, UTCoffsetH, stateCode)) %>%
    distinct()
  return(NRSTrips)
}

# Bring in all NRS samples
getNRSSamples <- function(){
    NRSSamp <- read_csv(paste0(rawD,.Platform$file.sep,"PSampNRS.csv"), na = "(null)") %>% 
      rename(Sample = SAMPLE, Station = STATION, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateLocal = SAMPLEDATE, NRScode = NRS_CODE) %>%
      mutate(Year = year(SampleDateLocal),
             Month = month(SampleDateLocal),
             Day = day(SampleDateLocal),
             Time_24hr = str_sub(SampleDateLocal, -8, -1)) # hms doesn"t seem to work on 00:00:00 times
    return(NRSSamp)
}

# Bring in plankton data
getNRSPhytoData <- function(){
  NRSPdat <- read_csv(paste0(rawD,.Platform$file.sep,"NRS_phyto_raw.csv"), na = "(null)") %>%
    rename(Sample = SAMPLE, TaxonName = TAXON_NAME, TaxonGroup = TAXON_GROUP, Genus = GENUS, Species = SPECIES, 
           Cells_L = CELL_PER_LITRE, Biovolume_uM3_L = BIOVOLUME_UM3_PER_L)
  return(NRSPdat)
}

# Bring in Change Log
getNRSPhytoChangeLog <- function(){
  NRSPcl <- read_csv(paste0(rawD,.Platform$file.sep,"ChangeLogNRSP.csv"), na = "(null)") %>%
    rename(TaxonName = TAXON_NAME, StartDate = START_DATE, ParentName = PARENT_NAME)
  return(NRSPcl)
}

# Bring in all NRS zooplankton samples
# getNRSZooSamples <- function(){
#   NRSZsamp <- read_csv(paste0(rawD,.Platform$file.sep,"ZSampNRS.csv"), na = "(null)") %>% 
#     rename(Sample = SAMPLE, Station = STATION, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateLocal = SAMPLEDATE, NRScode = NRS_CODE) %>%
#     mutate(Year = year(SampleDateLocal),
#            Month = month(SampleDateLocal),
#            Day = day(SampleDateLocal),
#            Time_24hr = str_sub(SampleDateLocal, -8, -1)) # hms doesn"t seem to work on 00:00:00 times
#   return(NRSZsamp)
# }

# Bring in zooplankton  abundance data
getNRSZooData <- function(){
  NRSZdat <- read_csv(paste0(rawD,.Platform$file.sep,"NRS_zoop_raw.csv"), na = "(null)") %>%
    rename(Sample = SAMPLE, TaxonName = TAXON_NAME, Copepod = TAXON_GROUP, TaxonGroup = TAXON_GRP01, 
         Genus = GENUS, Species = SPECIES, ZAbund_m3 = TAXON_PER_M3)
  return(NRSZdat)
}

# Bring in zooplankton  abundance data
getNRSZooCount <- function(){
  NRSZcount <- read_csv(paste0(rawD,.Platform$file.sep,"NRS_zoop_count_raw.csv"), na = "(null)") %>%
    dplyr::rename("TaxonName" = "TAXON_NAME", "Copepod" = "TAXON_GROUP", 
                  "TaxonGroup" = "TAXON_GRP01", "NRScode" = "NRS_CODE",
                  "Genus" = "GENUS", "Species" = "SPECIES", "TaxonCount" = "TAXON_COUNT")
  return(NRSZcount)
}

# Bring in Change Log
getNRSZooChangeLog <- function(){
  NRSZcl <- read_csv(paste0(rawD,.Platform$file.sep,"ChangeLogNRSZ.csv"), na = "(null)") %>%
    rename(TaxonName = TAXON_NAME, StartDate = START_DATE, ParentName = PARENT_NAME)
  return(NRSZcl)
}

# Bring in Zooplankton biomass
getNRSZooBiomass <- function(){
  ZBiomass <-  read_csv(paste0(rawD,.Platform$file.sep,"nrs_biomass.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_CODE, Biomass_mgm3 = BIOMASS_MGM3) %>% 
  mutate(SampleDepth_m = "WC",
         NRScode = gsub('^.{3}|.{9}$', '', NRScode)) %>% 
  untibble()
  return(ZBiomass)
}

# Bring in copepod information table with sizes etc.
get_ZooInfo <- function(){
  ZInfo <- read_csv(paste0(rawD,.Platform$file.sep,"taxon_info.csv"), na = "(null)") %>% 
    dplyr::rename( "TaxonName" = "TAXON_NAME") %>% 
    untibble()
}


# Bring in chemistry data
getChemistry <- function(){
  chemistry <- read_csv(paste0(rawD,.Platform$file.sep,"chemistry.csv"), na = c("(null)", NaN)) %>% 
  rename(NRScode = NRS_TRIP_CODE,
         SampleDepth_m = SAMPLE_DEPTH_M, Silicate_umol_L = SILICATE_UMOL_PER_L, Nitrate_umol_L =  NITRATE_UMOL_PER_L,
         Phosphate_umol_L =  PHOSPHATE_UMOL_PER_L, Salinity = SALINITY, 
         Ammonium_umol_L =  AMMONIUM_UMOL_PER_L,
         Nitrite_umol_L =  NITRITE_UMOL_PER_L,
         TCO2_umol_kg =  TCO2_UMOL_PER_KG,
         TAlkalinity_umol_kg =  TALKALINITY_UMOL_PER_KG,
         Oxygen_umol_L =  OXYGEN_UMOL_PER_L) %>% 
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         NRScode = substring(NRScode,4),
         Silicate_umol_L = ifelse(SILICATE_FLAG %in% c(3,4,9), NA, Silicate_umol_L), # remove all data flagged as bad or probably bad
         Phosphate_umol_L = ifelse(PHOSPHATE_FLAG %in% c(3,4,9), NA, Phosphate_umol_L),
         Ammonium_umol_L = ifelse(AMMONIUM_FLAG %in% c(3,4,9), NA, Ammonium_umol_L),
         Nitrate_umol_L = ifelse(NITRATE_FLAG %in% c(3,4,9), NA, Nitrate_umol_L),
         Nitrite_umol_L = ifelse(NITRITE_FLAG %in% c(3,4,9), NA, Nitrite_umol_L),
         Oxygen_umol_L = ifelse(OXYGEN_FLAG %in% c(3,4,9), NA, Oxygen_umol_L),
         TCO2_umol_kg = ifelse(CARBON_FLAG %in% c(3,4,9), NA, TCO2_umol_kg),
         TAlkalinity_umol_kg = ifelse(ALKALINITY_FLAG %in% c(3,4,9), NA, TAlkalinity_umol_kg),
         Salinity = ifelse(SALINITY_FLAG %in% c(3,4,9), NA, Salinity)) %>%
  group_by(NRScode, SampleDepth_m) %>% 
  summarise(Silicate_umol_L = mean(Silicate_umol_L, na.rm = TRUE), # some replicated samples from error picking up PHB data, needs addressing in database
            Phosphate_umol_L = mean(Phosphate_umol_L, na.rm = TRUE),
            Ammonium_umol_L = mean(Ammonium_umol_L, na.rm = TRUE),
            Nitrate_umol_L = mean(Nitrate_umol_L, na.rm = TRUE),
            Nitrite_umol_L = mean(Nitrite_umol_L, na.rm = TRUE),
            Oxygen_umol_L = mean(Oxygen_umol_L, na.rm = TRUE),
            TCO2_umol_kg = mean(TCO2_umol_kg, na.rm = TRUE),
            TAlkalinity_umol_kg = mean(TAlkalinity_umol_kg, na.rm = TRUE),
            Salinity = mean(Salinity, na.rm = TRUE),
            .groups = "drop") %>% 
  ungroup() %>% 
  mutate_all(~ replace(., is.na(.), NA)) %>% 
  untibble()
  return(chemistry)
}

# get CTD data
getCTD <- function(){
  CTD <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_CTD.csv"), na = "(null)",
                col_types = cols(PRES = col_double(), # columns start with nulls so tidyverse annoyingly assigns col_logical()
                                 PAR = col_double(),
                                 SPEC_CNDC = col_double())) %>% 
  rename(NRScode = NRS_TRIP_CODE, SampleDepth_m = PRES_REL, CTDDensity_kgm3 = DENS, 
         CTDTemperature = TEMP, CTDPAR_umolm2s = PAR,
         CTDConductivity_sm = CNDC, CTDSpecificConductivity_Sm = SPEC_CNDC, 
         CTDSalinity = PSAL, CTDTurbidity_ntu = TURB, CTDChlF_mgm3 = CHLF) %>%
  untibble()
  
  CTD_remove <- CTD %>% group_by(NRScode) %>% summarise(n = n()) %>%
    filter(n == 1)
  CTD <- CTD %>% filter(!NRScode %in% CTD_remove$NRScode) # removing records where CTD length is one value
  
  return(CTD)
}

