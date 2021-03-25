## IMOS plankton data products Indices 
## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)

## Created: Sept 2020
## Updated: 11 Nov 2020

suppressPackageStartupMessages({
  library(lubridate)
  library(vegan)
  library(data.table)
  library(rgdal)
  library(rlang)
  library(tidyverse)
})

raw <- "https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/RawData/"
output <- "https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/Output/"
get_sat_data <- FALSE

source("IMOS_Plankton_functions.R")
# uses mostly the same raw data from IMOS_PlanktonProducts_Create.R

# ensure we have all trips accounted for 
# note there are circumstances where a trip won't have a phyto and a zoo samples due to loss of sample etc.

cprTrips <- read_csv(paste0(raw, "CPR_Samp.csv"), na = "(null)") %>% 
  rename(Sample = SAMPLE, Route = ROUTE, Region = REGION, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateUTC = SAMPLEDATEUTC, 
         SampleType = SAMPLETYPE, Biomass_mgm3 = BIOMASS_MGM3) %>%
  filter(!is.na(SampleType)) %>%
  mutate(Year = year(SampleDateUTC),
         Month = month(SampleDateUTC),
         Day = day(SampleDateUTC),
         Time_24hr = str_sub(SampleDateUTC, -8, -1), # hms doesn"t seem to work on 00:00:00 times
         SampleDateUTC = as.character(SampleDateUTC)) %>% 
  select(c(Sample, Latitude:Time_24hr, Region, Route, PCI, Biomass_mgm3))

mbr <-  readOGR("../General/Shapefiles/marine_regions")
mbr <- spTransform(mbr, CRS("+proj=longlat +datum=WGS84"))
## Read the segments.

segments <- cprTrips %>% 
  select(Longitude, Latitude) # file with columns named Longitude, Latitude
coordinates(segments) <- c("Longitude", "Latitude") ## Convert the segments table to a SpatialPolygonsDataFrame.
proj4string(segments) <- CRS("+proj=longlat +datum=WGS84")
segmbr <- over(segments, mbr) ## Perform the point in polygon overlay.
cprTrips$BioRegion <- segmbr$REGION


## Temporary fix for state-based waters near bioregions
# ggplot(data = cprTrips, aes(x = Longitude, y = Latitude, colour = BioRegion)) + geom_point()

cprTrips <- cprTrips %>% 
  mutate(BioRegion = case_when(Longitude >= mbr[17,]@bbox[1,1] &
                                 Longitude <= mbr[17,]@bbox[1,2] &
                                 Latitude >= mbr[17,]@bbox[2,1] &
                                 Latitude <= mbr[17,]@bbox[2,2] &
                                 are_na(BioRegion) == TRUE ~ "Coral Sea",
                               Longitude >= mbr[2,]@bbox[1,1] &
                                 Longitude <= mbr[2,]@bbox[1,2] &
                                 Latitude >= mbr[2,]@bbox[2,1] &
                                 Latitude <= mbr[2,]@bbox[2,2] &
                                 are_na(BioRegion) == TRUE ~ "Temperate East",
                               Longitude >= mbr[7,]@bbox[1,1] &
                                 Longitude <= mbr[7,]@bbox[1,2] &
                                 Latitude >= mbr[7,]@bbox[2,1] &
                                 Latitude <= mbr[7,]@bbox[2,2] &
                                 are_na(BioRegion) == TRUE ~ "South-west",
                               Longitude >= mbr[13,]@bbox[1,1] &
                                 Longitude <= mbr[13,]@bbox[1,2] &
                                 Latitude >= mbr[13,]@bbox[2,1] &
                                 Latitude <= mbr[13,]@bbox[2,2] &
                                 are_na(BioRegion) == TRUE ~ "South-east",
                               TRUE ~ BioRegion))

# ggplot(data = cprTrips, aes(x = Longitude, y = Latitude, colour = BioRegion)) + geom_point()

cprProps <- read_csv(paste0(raw, "CPR_SatData.csv"), na = "(null)") %>% 
  rename(Sample = SAMPLE, ChlorophyllSatellite_mgm3 = CHLA, WaterDepth_m = DEPTH_M)

satcpr <- cprTrips %>% 
  rename(Date = SampleDateUTC) 

if (get_sat_data == TRUE){
  # accessing the satelitte data from MODIS
  
  # Possible products
  # pr <- c("sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
  #         "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
  #         "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oc3", "chl_gsm", "K_490")
  
  pr <- c("sst", "chl_oc3")
  res_temp <- "1d"
  res_spat <- 10 # Return the average of res_spat x res_spat pixels
  
  # Get MODIS Data
  satcpr <- fIMOS_MatchMODIS(satcpr, pr, res_temp, res_spat)
  
  # Get Altimetry (Gridded sea level anomaly, Gridded sea level, Surface geostrophic velocity)
  satcpr <- fIMOS_MatchAltimetry(satcpr, res_spat)
} else{
  satcpr <- satcpr %>% 
    add_column(sst_1d = NA, chl_oc3_1d = NA)
}

# Total zoop abundance
zoodatacpr <-  cprZsamp %>% 
  left_join(cprZdat, by = "Sample")

TZoocpr <-  zoodatacpr %>% 
  group_by(Sample) %>% 
  summarise(ZoopAbundance_m3 = sum(ZAbun_m3, na.rm = TRUE))

TCopecpr <- zoodatacpr %>% 
  filter(Copepod == 'COPEPOD') %>% 
  group_by(Sample) %>% 
  summarise(CopeAbundance_m3 = sum(ZAbun_m3, na.rm = TRUE))

# Bring in copepod information table with sizes etc.
Zinfo <- get_ZooInfo() 

ACopeSizeCpr <- zoodatacpr %>% 
  filter(Copepod == 'COPEPOD') %>%
  inner_join(Zinfo %>% select(SIZE_AVE_MM, TaxonName, DIET), by = "TaxonName") %>%
  mutate(abunSize = SIZE_AVE_MM * ZAbun_m3, 
         DIET = ifelse(DIET == 'CC', 'CC', 'CO')) %>%
  group_by(Sample) %>% summarise(AvgTotalLengthCopepod_mm = sum(abunSize, na.rm = TRUE)/sum(ZAbun_m3, na.rm = TRUE))

HCratCpr <- zoodatacpr %>% 
  filter(Copepod == 'COPEPOD') %>%
  inner_join(Zinfo %>% select(TaxonName, DIET), by = "TaxonName") %>%
  mutate(DIET = ifelse(DIET == 'CC', 'CC', 'CO')) %>% drop_na() %>%
  select(Sample, DIET, ZAbun_m3) %>% 
  group_by(Sample, DIET) %>% summarise(sumdiet = sum(ZAbun_m3 , na.rm = TRUE)) %>% 
  pivot_wider(values_from = sumdiet, names_from = DIET) %>%
  mutate(HerbivoreCarnivoreCopepodRatio = CO / (CO + CC)) %>% untibble

# Diversity, evenness etc.     

# Bring in plankton data
CPRZcount <- read_csv(paste0(raw, "CPR_Zoop_CountRaw.csv"), na = "(null)") %>%
  rename(TaxonName = TAXON_NAME, Copepod = TAXON_GROUP, TaxonGroup = TAXON_GRP01, Sample = SAMPLE,
                Genus= GENUS, Species = SPECIES, TaxonCount = COUNTS, SampleVol_L = SAMPVOL_L)

zooCountCpr <-  cprZsamp %>% 
  left_join(CPRZcount, by = "Sample")

nCPR <-  zooCountCpr %>% 
  filter(Copepod == 'COPEPOD' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample) %>% 
  summarise(NoCopepodSpecies_Sample = n())

ShannonCopepodDiversityCPR <- zooCountCpr %>% 
  filter(Copepod == 'COPEPOD' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample, TaxonName) %>% summarise(ZCount = sum(TaxonCount, na.rm = TRUE)) %>%
  pivot_wider(values_from = ZCount, names_from = TaxonName, values_fill = 0) %>% ungroup() %>%
  select(-Sample) %>%
  diversity('shannon')

CopepodEvennessCPR <- nCPR %>% 
  cbind(ShannonCopepodDiversityCPR) %>% 
  mutate(CopepodEvenness = ShannonCopepodDiversityCPR / log(NoCopepodSpecies_Sample))

# Total Phyto abundance
phytodatacpr <- cprPsamp %>% 
  left_join(cprPdat, by = "Sample") %>% 
  filter(TaxonGroup != 'Other')

PhytoCcpr <- phytodatacpr %>% 
  select(Sample, TaxonGroup, PAbun_m3, BioVolume_um3m3) %>% 
  mutate(BV_Cell = BioVolume_um3m3 / PAbun_m3, # biovolume of one cell
         Carbon = ifelse(TaxonGroup == 'Dinoflagellate', 0.76*(BV_Cell)^0.819, # conversion to Carbon based on taxongroup and biovolume of cell
                         ifelse(TaxonGroup == 'Ciliate', 0.22*(BV_Cell)^0.939,
                                ifelse(TaxonGroup == 'Cyanobacteria', 0.2, 0.288*(BV_Cell)^0.811 ))),
         Carbon_m3 = PAbun_m3 * Carbon) %>% # Carbon per m3
  group_by(Sample) %>% 
  summarise(PhytoBiomassCarbon_pgm3 = sum(Carbon_m3))

TPhytoCpr <- phytodatacpr %>% 
  group_by(Sample) %>% 
  summarise(AbundancePhyto_cells_m3 = sum(PAbun_m3, na.rm = TRUE))

DDratcpr <- phytodatacpr %>%
  filter(TaxonGroup %in% c('Centric diatom', "Pennate diatom", 'Dinoflagellate')) %>% 
  mutate(TaxonGroup = recode(TaxonGroup, 'Centric diatom' = 'Diatom', 'Pennate diatom' = 'Diatom')) %>%
  select(Sample, TaxonGroup, PAbun_m3) %>% 
  group_by(Sample, TaxonGroup) %>% summarise(sumTG = sum(PAbun_m3, na.rm = TRUE)) %>% 
  pivot_wider(values_from = sumTG, names_from = TaxonGroup) %>%
  mutate(DiatomDinoflagellateRatio = Diatom / (Diatom + Dinoflagellate)) %>% 
  untibble

AvgCellVolcpr <- phytodatacpr %>% 
  filter(!is.na(BioVolume_um3m3)) %>% 
  group_by(Sample) %>% 
  summarise(AvgCellVol_um3 = mean(sum(BioVolume_um3m3)/sum(PAbun_m3)))

# Diversity (phyto, diatoms, dinos)
# stick to abundance data here as otherwise we have FOV counts

npcpr <-  phytodatacpr %>% 
  filter(TaxonGroup != 'Other' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample) %>% 
  summarise(NoPhytoSpecies_Sample = n())

ShannonPhytoDiversitycpr <- phytodatacpr %>% 
  filter(TaxonGroup != 'Other' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample, TaxonName) %>% 
  summarise(Pdata = sum(PAbun_m3, na.rm = TRUE)) %>%
  pivot_wider(values_from = Pdata, names_from = TaxonName, values_fill = 0) %>% 
  ungroup() %>%
  select(-Sample) %>%
  diversity('shannon')

PhytoEvencpr <- npcpr %>% 
  cbind(ShannonPhytoDiversitycpr) %>% 
  mutate(PhytoEvenness = ShannonPhytoDiversitycpr / log(NoPhytoSpecies_Sample))

ndiacpr <-  phytodatacpr %>% 
  filter(TaxonGroup %in% c('Centric diatom', 'Pennate diatom') & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample) %>% summarise(NoDiatomSpecies_Sample = n())

ShannonDiatomDiversitycpr <- phytodatacpr %>% 
  filter(TaxonGroup %in% c('Centric diatom', 'Pennate diatom') & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample, TaxonName) %>% 
  summarise(Diadata = sum(PAbun_m3, na.rm = TRUE)) %>%
  pivot_wider(values_from = Diadata, names_from = TaxonName, values_fill = 0) %>% 
  ungroup() %>%
  select(-Sample) %>%
  diversity('shannon')

DiaEvencpr <- ndiacpr %>% 
  cbind(ShannonDiatomDiversitycpr) %>% 
  mutate(DiatomEvenness = ShannonDiatomDiversitycpr / log(NoDiatomSpecies_Sample))

ndinocpr <- phytodatacpr %>% 
  filter(TaxonGroup == 'Dinoflagellate' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample) %>% summarise(NoDinoSpecies_Sample = n())

ShannonDinoDiversitycpr <- phytodatacpr %>% 
  filter(TaxonGroup  == 'Dinoflagellate' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample, TaxonName) %>% 
  summarise(Dinodata = sum(PAbun_m3, na.rm = TRUE)) %>%
  pivot_wider(values_from = Dinodata, names_from = TaxonName, values_fill = 0) %>% 
  ungroup() %>%
  select(-Sample) %>%
  diversity('shannon')

DinoEvencpr <- ndinocpr %>% 
  cbind(ShannonDinoDiversitycpr) %>% 
  mutate(DinoflagellateEvenness = ShannonDinoDiversitycpr / log(NoDinoSpecies_Sample))

# make indices table (nrows must always equal nrows of Trips)
IndicesCPR <-  cprTrips  %>% 
  left_join(TZoocpr, by = ("Sample")) %>%
  left_join(TCopecpr, by = ("Sample")) %>%
  left_join(ACopeSizeCpr, by = ("Sample")) %>%
  left_join(HCratCpr %>% select(-c('CO', 'CC')), by = ("Sample")) %>% 
  left_join(CopepodEvennessCPR,  by = ("Sample")) %>%
  left_join(PhytoCcpr, by = ("Sample")) %>%
  left_join(TPhytoCpr, by = ("Sample")) %>%
  left_join(DDratcpr %>% select(-c('Diatom', 'Dinoflagellate')), by = ("Sample")) %>%
  left_join(AvgCellVolcpr, by = ("Sample")) %>%
  left_join(PhytoEvencpr, by = ("Sample")) %>%
  left_join(DiaEvencpr, by = ("Sample")) %>%
  left_join(DinoEvencpr, by = ("Sample")) %>%   
#  left_join(satcpr %>% select(Sample, sst_1d, chl_oc3_1d), by = ("Sample")) %>%  #add once run , GSLA, GSL, UCUR, VCUR
  select(-Sample, -SampleType)

# make indices table (nrows must always equal nrows of Trips) - old one for IMOS
IMOSCPR_ind <-  cprTrips  %>%
  left_join(cprProps, by = ("Sample")) %>%
  left_join(TZoocpr, by = ("Sample")) %>%
  left_join(TCopecpr, by = ("Sample")) %>%
  left_join(ACopeSizeCpr, by = ("Sample")) %>%
  left_join(HCratCpr %>% select(-c('CO', 'CC')), by = ("Sample")) %>% 
  left_join(CopepodEvennessCPR,  by = ("Sample")) %>%
  left_join(PhytoCcpr, by = ("Sample")) %>%
  left_join(TPhytoCpr, by = ("Sample")) %>%
  left_join(DDratcpr %>% select(-c('Diatom', 'Dinoflagellate')), by = ("Sample")) %>%
  left_join(AvgCellVolcpr, by = ("Sample")) %>%
  left_join(PhytoEvencpr, by = ("Sample")) %>%
  left_join(DiaEvencpr, by = ("Sample")) %>%
  left_join(DinoEvencpr, by = ("Sample")) %>%   
  select(-Sample)

fwrite(IndicesCPR, file = paste0(output, "CPR_Indices.csv"), row.names = FALSE)

# test table
# n should be 1, replicates or duplicate samples will have values > 1
test <- IndicesCPR %>% 
  group_by(Latitude, Longitude, SampleDateUTC) %>% 
  summarise(n = n())

max(test$n)

