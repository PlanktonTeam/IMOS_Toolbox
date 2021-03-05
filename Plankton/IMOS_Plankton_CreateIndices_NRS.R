## IMOS plankton data products Indices 
## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)

## Created: Sept 2020
## Updated: 
  ## 1 Oct 2020 (Written to Git)
  ## 6th October 2020

suppressPackageStartupMessages({
  library(lubridate)
  library(lutz)
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

NRSdat <- get_NRSTrips() %>% #ignore warning, 'fast' method does better here than 'accurate'
  select(-SampleDepth_m) %>% 
  distinct()

dNRSdat <- distinct(NRSdat, NRScode, .keep_all = TRUE) %>%  # Distinct rows for satellite, should be anyway
  rename(Date = SampleDateLocal) %>% 
  select(NRScode, Date, Latitude, Longitude)

# SST and Chlorophyll from CTD
CTD <- getCTD() %>%
  filter(SampleDepth_m < 15) %>% # take average of top 10m as a surface value for SST and CHL, this is removing 17 casts as of nov 2020
  group_by(NRScode) %>% 
  summarise(CTD_SST_C = mean(CTDTemperature, na.rm = TRUE),
            CTDChlF_mgm3 = mean(CTDChlF_mgm3, na.rm = TRUE),
            CTDSalinity = mean(CTDSalinity, na.rm = TRUE),
            .groups = "drop") %>%
  untibble()

# data set for calculating MLD
CTD_MLD <- getCTD() %>% 
  select(NRScode, CTDTemperature, CTDSalinity, CTDChlF_mgm3, SampleDepth_m)

n = nrow(CTD_MLD %>% select(NRScode) %>% unique())
MLD <- data.frame(NRScode = NA, MLD_temp = as.numeric(NA), MLD_sal = as.numeric(NA), DCM = as.numeric(NA))

# MLD by T and S (Ref: Condie & Dunn 2006)
# DCM from max f from CTD
for (i in 1:n) {
  dat <- CTD_MLD %>% select(NRScode) %>% unique() %>% mutate(NRScode = as.factor(NRScode))
  nrscode <- dat$NRScode[[i]] %>% droplevels()
  mldData <- CTD_MLD %>% filter(NRScode == nrscode) %>% arrange(SampleDepth_m)
    
  if (as.character(substr(nrscode, 0,3)) %in% c("DAR", "YON")){
    refDepth  <-  5
  }  
    
  if (!as.character(substr(nrscode, 0,3)) %in% c("DAR", "YON")){
    refDepth  <-  10
  }  
  
  ref_T <- mldData %>% mutate(refd = abs(SampleDepth_m - refDepth), # find depth nearest to 10 m
                              rankrefd = ave(refd, FUN = . %>% order %>% order)) %>%
    filter(rankrefd == 1)
  refT <- ref_T$CTDTemperature - 0.4 # temp at 10 m minus 0.4 deg C
  mldData <- mldData %>% filter(SampleDepth_m > ref_T$SampleDepth_m)
  mld_t <- mldData %>% mutate(temp = abs(CTDTemperature - refT),
                              ranktemp = ave(temp, FUN = . %>% order %>% order)) %>%
    filter(ranktemp == 1)
  MLD_temp <- mld_t$SampleDepth_m
  
  refS <- ref_T$CTDSalinity - 0.03 # temp at 10 m minus 0.4
  mld_s <- mldData %>% mutate(temp = abs(CTDSalinity - refS),
                              ranksal = ave(temp, FUN = . %>% order %>% order)) %>%
    filter(ranksal == 1)
  MLD_sal <- mld_s$SampleDepth_m
  
  dcm <- (mldData %>% filter(CTDChlF_mgm3 == max(CTDChlF_mgm3)))$SampleDepth_m 
  dcm[is_empty(dcm)] = NA
  
  MLD <- rbind(MLD, data.frame(NRScode = as.character(nrscode), MLD_temp = MLD_temp, MLD_sal = MLD_sal, DCM = dcm), stringsAsFactors = FALSE)  %>% drop_na(NRScode)    
}

# # Access satellite data for the sample dates using the IMOS_Toolbox
# 
# # If on Windows you will need to install a development 
# # version of ncdf4 which allows the use of OpenDAP
# if(.Platform$OS.type == "windows") {
#   warning("It looks like you are on a Windows PC - You will need to install a 
#   development version of ncdf4 which allows the use of OpenDAP. Please 
#   run devtools::install_github('mdsumner/ncdf4') to install or 
#   see 'https://github.com/mdsumner/ncdf4' for more information.")
# }
# 
# # Get GHRSST SST Data
# # Possible products to download are: 
# # dt_analysis, l2p_flags, quality_level, satellite_zenith_angle, sea_ice_fraction, sea_ice_fraction_dtime_from_sst, 
# # sea_surface_temperature, sea_surface_temperature_day_night, sses_bias, sses_count,sses_standard_deviation,
# # sst_count, sst_dtime, sst_mean, sst_standard_deviation, wind_speed, wind_speed_dtime_from_sst,
# res_temp <- "1d"
# res_spat <- 10 # Return the average of res_spat x res_spat pixels
# pr <- ("sea_surface_temperature")
# GHRSST <- fIMOS_MatchGHRSST(dNRSdat, pr, res_temp, res_spat)
# 
# # Get MODIS Data
# # Possible products
# # pr <- c("sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
# #         "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
# #         "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oc3", "chl_gsm", "K_490")
# 
# pr <- c("chl_oci")
# res_temp <- "1d"
# res_spat <- 10 # Return the average of res_spat x res_spat pixels
# MODIS <- fIMOS_MatchMODIS(dNRSdat, pr, res_temp, res_spat)
# 
# 
# # Get Altimetry (Gridded sea level anomaly, Gridded sea level, Surface geostrophic velocity)
# dNRSdat <- dNRSdat[1:3,]
# Alt <- fIMOS_MatchAltimetry(dNRSdat, res_spat)


# Nutrient data

Nuts <- getChemistry() %>% 
  group_by(NRScode) %>% 
  summarise(Silicate_umol_L = mean(Silicate_umol_L, na.rm = TRUE),
            Phosphate_umol_L = mean(Phosphate_umol_L, na.rm = TRUE),
            Ammonium_umol_L = mean(Ammonium_umol_L, na.rm = TRUE),
            Nitrate_umol_L = mean(Nitrate_umol_L, na.rm = TRUE),
            Nitrite_umol_L = mean(Nitrite_umol_L, na.rm = TRUE),
            Oxygen_umol_L = mean(Oxygen_umol_L, na.rm = TRUE),
            TCO2_umol_kg = mean(TCO2_umol_kg, na.rm = TRUE),
            TAlkalinity_umol_kg = mean(TAlkalinity_umol_kg, na.rm = TRUE),
            Salinity_psu = mean(Salinity, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate_all(~ replace(., is.na(.), NA)) %>% 
  untibble()

Pigments <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_pigments.csv"), na = "(null)") %>% 
  select(NRS_TRIP_CODE, SAMPLE_DEPTH_M, DV_CPHL_A_AND_CPHL_A) %>% 
  rename(NRScode = NRS_TRIP_CODE, SampleDepth_m = SAMPLE_DEPTH_M, Chla = DV_CPHL_A_AND_CPHL_A) %>%
  filter(SampleDepth_m <= 25) %>% # take average of top 10m as a surface value for SST and CHL
  # filter(SampleDepth_m == "WC") %>% 
  mutate(NRScode = str_replace(NRScode, "NRS", "")) %>% 
  group_by(NRScode) %>% 
  summarise(Chla_mgm3 = mean(Chla, na.rm = TRUE),
            .groups = "drop") %>%
  untibble()

# 
# dat <- left_join(Pigments, NRSdat, by = "NRScode")
# dat$SampleDepth_m <- as.numeric(str_replace_all(dat$SampleDepth_m, "WC", "150"))
# ggplot(data = dat, aes(x = Year, y = SampleDepth_m)) + geom_point()

# Total Zooplankton Abundance
ZooData <- getNRSSamples() %>% 
  left_join(getNRSZooData(), by = "Sample")

TZoo <- ZooData %>% 
  group_by(NRScode) %>% 
  summarise(ZoopAbundance_m3 = sum(ZAbund_m3, na.rm = TRUE),
            .groups = "drop")

TCope <- ZooData %>% 
  filter(Copepod == 'COPEPOD') %>% 
  group_by(NRScode, ) %>% 
  summarise(CopeAbundance_m3 = sum(ZAbund_m3, na.rm = TRUE),
            .groups = "drop")

# Total Zooplankton Biomass
ZBiomass <- getNRSSamples() %>% select(c(NRScode, Biomass_mgm3)) %>% 
  mutate(SampleDepth_m = "WC")

# Bring in copepod information table with sizes etc.
ZInfo <- get_ZooInfo() 

ACopeSize <- ZooData %>% 
  filter(Copepod == 'COPEPOD') %>%
  inner_join(ZInfo %>% select(SIZE_AVE_MM, TaxonName, DIET), by = "TaxonName") %>%
  mutate(abunSize = SIZE_AVE_MM * ZAbund_m3, 
         DIET = ifelse(DIET == 'CC', 'CC', 'CO')) %>%
  group_by(NRScode) %>% 
  summarise(AvgTotalLengthCopepod_mm = sum(abunSize, na.rm = TRUE)/sum(ZAbund_m3, na.rm = TRUE),
            .groups = "drop")

HCrat <- ZooData %>% 
  filter(Copepod == 'COPEPOD') %>%
  inner_join(ZInfo %>% select(TaxonName, DIET), by = "TaxonName") %>%
  mutate(DIET = ifelse(DIET == 'CC', 'CC', 'CO')) %>% 
  drop_na() %>%
  select(NRScode, DIET, ZAbund_m3) %>% 
  group_by(NRScode, DIET) %>% 
  summarise(sumdiet = sum(ZAbund_m3 , na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(values_from = sumdiet, names_from = DIET) %>%
  mutate(HerbivoreCarnivoreCopepodRatio = CO / (CO + CC)) %>% 
  untibble()

# Diversity, evenness etc.     

# Bring in plankton data
ZooCount <- getNRSSamples() %>% 
  left_join(getNRSZooCount(), by = "NRScode")

n <- ZooCount %>% 
  filter(Copepod == 'COPEPOD' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode) %>% 
  summarise(NoCopepodSpecies_Sample = n(), .groups = "drop")

ShannonCopepodDiversity <- ZooCount %>% 
  filter(Copepod == 'COPEPOD' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode, TaxonName) %>% 
  summarise(ZCount = sum(TaxonCount, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(values_from = ZCount, names_from = TaxonName, values_fill = 0) %>% 
  ungroup() %>%
  select(-NRScode) %>%
  diversity('shannon')

CopepodEvenness <- n %>% 
  cbind(ShannonCopepodDiversity) %>% 
  mutate(CopepodEvenness = ShannonCopepodDiversity / log(NoCopepodSpecies_Sample))

# Total Phyto abundance
PhytoData <- getNRSSamples() %>% 
  left_join(getNRSPhytoData(), by = "Sample") %>% 
  filter(TaxonGroup != 'Other')

# PhytoData <- PhytoData %>%
#   filter(str_detect(TaxonName, "Flagellate <10", negate = TRUE)) # Remove flagellates

PhytoC <- PhytoData %>% 
  select(NRScode, TaxonGroup, Cells_L, Biovolume_uM3_L) %>% 
  mutate(BV_Cell = Biovolume_uM3_L / Cells_L, # biovolume of one cell
         Carbon = ifelse(TaxonGroup == 'Dinoflagellate', 0.76*(BV_Cell)^0.819, # conversion to Carbon based on taxongroup and biovolume of cell
                         ifelse(TaxonGroup == 'Ciliate', 0.22*(BV_Cell)^0.939,
                                ifelse(TaxonGroup == 'Cyanobacteria', 0.2, 0.288*(BV_Cell)^0.811 ))),
         Carbon_L = Cells_L * Carbon) %>% # Carbon per litre
  group_by(NRScode) %>% 
  summarise(PhytoBiomassCarbon_pg_L = sum(Carbon_L),
            .groups = "drop")

TPhyto <-  PhytoData %>% 
  group_by(NRScode) %>% 
  summarise(AbundancePhyto_cells_L = sum(Cells_L, na.rm = TRUE),
            .groups = "drop")

DDrat <- PhytoData %>%
  filter(TaxonGroup %in% c('Centric diatom', "Pennate diatom", 'Dinoflagellate')) %>% 
  mutate(TaxonGroup = recode(TaxonGroup, 'Centric diatom' = 'Diatom', 'Pennate diatom' = 'Diatom')) %>%
  select(NRScode, TaxonGroup, Cells_L) %>% 
  group_by(NRScode, TaxonGroup) %>% 
  summarise(sumTG = sum(Cells_L, na.rm = TRUE),
            .groups = "drop") %>% 
  pivot_wider(values_from = sumTG, names_from = TaxonGroup) %>%
  mutate(DiatomDinoflagellateRatio = Diatom / (Diatom + Dinoflagellate)) %>% 
  untibble()

AvgCellVol <- PhytoData %>% 
  filter(!is.na(Biovolume_uM3_L)) %>% 
  group_by(NRScode) %>% 
  summarise(AvgCellVol_um3 = mean(sum(Biovolume_uM3_L)/sum(Cells_L)),
            .groups = "drop")

# Diversity (phyto, diatoms, dinos)
# stick to abundance data here or we lose all the data that Pru counted which we don't have counts for.

NP <- PhytoData %>% 
  filter(TaxonGroup != 'Other' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode) %>% 
  summarise(NoPhytoSpecies_Sample = n(),
            .groups = "drop")

ShannonPhytoDiversity <- PhytoData %>% 
  filter(TaxonGroup != 'Other' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode, TaxonName) %>% 
  summarise(Pdata = sum(Cells_L, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(values_from = Pdata, names_from = TaxonName, values_fill = 0) %>% 
  ungroup() %>%
  select(-NRScode) %>%
  diversity('shannon')

PhytoEven <- NP %>% 
  cbind(ShannonPhytoDiversity) %>% 
  mutate(PhytoEvenness = ShannonPhytoDiversity / log(NoPhytoSpecies_Sample))

NDia <- PhytoData %>% 
  filter(TaxonGroup %in% c('Centric diatom', 'Pennate diatom') & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode) %>% 
  summarise(NoDiatomSpecies_Sample = n(),
            .groups = "drop")

ShannonDiatomDiversity <- PhytoData %>% 
  filter(TaxonGroup %in% c('Centric diatom', 'Pennate diatom') & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode, TaxonName) %>% 
  summarise(Diadata = sum(Cells_L, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(values_from = Diadata, names_from = TaxonName, values_fill = 0) %>% 
  ungroup() %>%
  select(-NRScode) %>%
  diversity('shannon')

DiaEven <- NDia %>% 
  cbind(ShannonDiatomDiversity) %>% 
  mutate(DiatomEvenness = ShannonDiatomDiversity / log(NoDiatomSpecies_Sample))

NDino <-  PhytoData %>% 
  filter(TaxonGroup == 'Dinoflagellate' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode) %>% 
  summarise(NoDinoSpecies_Sample = n(),
            .groups = "drop")

ShannonDinoDiversity <- PhytoData %>% 
  filter(TaxonGroup  == 'Dinoflagellate' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode, TaxonName) %>% 
  summarise(Dinodata = sum(Cells_L, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(values_from = Dinodata, names_from = TaxonName, values_fill = 0) %>%
  ungroup() %>%
  select(-NRScode) %>%
  diversity('shannon')

DinoEven <- NDino %>% 
  cbind(ShannonDinoDiversity) %>% 
  mutate(DinoflagellateEvenness = ShannonDinoDiversity / log(NoDinoSpecies_Sample))

# make indices table (nrows must always equal nrows of Trips)
Indices <- NRSdat  %>%
  left_join(TZoo, by = ("NRScode")) %>%
  left_join(TCope, by = ("NRScode")) %>%
  left_join(ZBiomass %>% select(NRScode, Biomass_mgm3), by = ("NRScode")) %>%
  left_join(ACopeSize, by = ("NRScode")) %>%
  left_join(HCrat %>% select(-c('CO', 'CC')), ("NRScode")) %>%
  left_join(CopepodEvenness,  by = ("NRScode")) %>%
  left_join(PhytoC, by = ("NRScode")) %>%
  left_join(TPhyto, by = ("NRScode")) %>%
  left_join(DDrat %>% select(-c('Diatom', 'Dinoflagellate')), by = ("NRScode")) %>%
  left_join(AvgCellVol, by = ("NRScode")) %>%
  left_join(PhytoEven, by = ("NRScode")) %>%
  left_join(DiaEven, by = ("NRScode")) %>%
  left_join(DinoEven, by = ("NRScode")) %>%  
  left_join(CTD, by = ("NRScode")) %>%
  left_join(MLD, by = ("NRScode")) %>%
  left_join(Nuts, by = ("NRScode")) %>% 
  left_join(Pigments, by = ("NRScode"))


# %>% 
#   left_join(GHRSST %>% select(-c("Longitude", "Latitude", "Date")), by = ("NRScode")) %>% 
#   left_join(MODIS %>% select(-c("Longitude", "Latitude", "Date")), by = ("NRScode")) %>% 
#   left_join(Alt %>% select(c("NRScode", "GSLA", "GSL", "UCUR", "VCUR")), by = ("NRScode"))

fwrite(Indices, file = paste0(outD,.Platform$file.sep,"NRS_Indices.csv"), row.names = FALSE)

# test table
# n should be 1, replicates or duplicate samples will have values > 1
test <- Indices %>% group_by(NRScode) %>% summarise(n = n())


