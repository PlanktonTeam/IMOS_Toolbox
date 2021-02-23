## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)
## IMOS plankton data products

## Created: May 2020
## Updated: 
## 21 July 2020 (Written to Git)
## 22 September 2020 (Updated data file structure)
## 13th November 2020 (Split CPR and NRS)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(data.table)
})

rawD <- "RawData"
outD <- "Output"

source("IMOS_Plankton_functions.R")

#### NRS Phytoplankton #### #################################################################################################################################

# Bring in all NRS phytoplankton samples, data and changelog
NRSPsamp <- getNRSPhytoSamples()
NRSPdat <- getNRSPhytoData()
NRSPcl <- getNRSPhytoChangeLog()

#### Raw Phytoplankton ####

NRSRawP1 <- left_join(NRSPsamp, NRSPdat, by = "Sample") %>% 
  select(-c(TaxonGroup, Genus, Species, Biovolume_uM3_L)) %>% 
  arrange(-desc(TaxonName)) 

NRSRawP <- NRSRawP1 %>% 
  pivot_wider(names_from = TaxonName, values_from = Cells_L, values_fill = list(Cells_L = 0)) %>% 
  arrange(desc(SampleDateLocal)) %>% 
  select(-Sample) %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(NRSRawP, file = paste0(outD,.Platform$file.sep,"NRS_phyto_raw_mat.csv"), row.names = FALSE)

#### Higher Trophic Groups Abund ####

NRSHTGP1 <- NRSPdat %>% 
  group_by(Sample, TaxonGroup) %>% 
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE), .groups = "drop") %>%
  filter(!TaxonGroup %in% c("Other","Coccolithophore", "Diatom","Protozoa")) 

NRSHTGP1 <- NRSPsamp %>% 
  left_join(NRSHTGP1, by = "Sample") %>% 
  mutate(TaxonGroup = ifelse(is.na(TaxonGroup), "Ciliate", TaxonGroup),
         Cells_L = ifelse(is.na(Cells_L), 0, Cells_L)) %>% 
  arrange(-desc(TaxonGroup))

NRSHTGP <-  NRSHTGP1 %>% 
  pivot_wider(names_from = TaxonGroup, values_from = Cells_L, values_fill = list(Cells_L = 0)) %>% 
  arrange(desc(SampleDateLocal)) %>% 
  select(-Sample)

fwrite(NRSHTGP, file = paste0(outD,.Platform$file.sep,"NRS_phyto_HTG_mat.csv"), row.names = FALSE)

#### Genus Abund ####

# Check genus are effected by change log
nrslg <- NRSPcl %>% 
  mutate(genus1 = word(TaxonName, 1),
         genus2 = word(ParentName, 1)) %>%
  mutate(same = ifelse(genus1==genus2, "yes", "no")) %>%
  filter(same == "no")# no changes at genera level

# for non change log species
NRSGenP1 <- NRSPdat %>% 
  filter(!TaxonName %in% levels(as.factor(nrslg$TaxonName))) %>% 
  group_by(Sample, Genus) %>% 
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE), .groups = "drop") %>% 
  drop_na(Genus) 

NRSGenP1 <- NRSPsamp %>% 
  left_join(NRSGenP1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Genus = ifelse(is.na(Genus), "Acanthoica", Genus),
         Cells_L = ifelse(is.na(Cells_L), 0, Cells_L)) %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(Cells_L = sum(Cells_L), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
NRSGenP2 <- NRSPdat %>% 
  filter(TaxonName %in% levels(as.factor(nrslg$TaxonName))) %>% 
  left_join(NRSPcl, by = "TaxonName") %>%
  mutate(Genus = as_factor(Genus)) %>% 
  drop_na(Genus) %>%
  group_by(Sample, StartDate, Genus) %>% 
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(NRSGenP2$Genus)) {
  Gen = NRSGenP2 %>% select(Genus) %>% unique() 
  Gen = as.character(Gen$Genus[i] %>% droplevels())
  Dates <- as.data.frame(NRSGenP2) %>% 
    filter(Genus == Gen) %>% 
    slice(1) %>% 
    droplevels()
  
  gen <- as.data.frame(NRSGenP2) %>% 
    filter(Genus == Gen) %>% 
    droplevels() 
  
  gen <- NRSPsamp %>% 
    left_join(gen, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Genus = replace(Genus, is.na(Genus), Dates$Genus),
           Cells_L = replace(Cells_L, StartDate>SampleDateLocal, -999), 
           Cells_L = replace(Cells_L, StartDate<SampleDateLocal & is.na(Cells_L), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
    summarise(Cells_L = sum(Cells_L), .groups = "drop") %>% 
    as.data.frame()     
  NRSGenP1 <- rbind(NRSGenP1, gen)
}

NRSGenP1 <- NRSGenP1 %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(Cells_L = max(Cells_L), .groups = "drop") %>% 
  arrange(-desc(Genus)) %>% 
  as.data.frame()

# select maximum value of duplicates, but leave -999 for all other occurrences as not regularly identified
NRSGenP <-  NRSGenP1 %>% 
  pivot_wider(names_from = Genus, values_from = Cells_L, values_fill = list(Cells_L = 0)) %>% 
  arrange(desc(SampleDateLocal))  %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(NRSGenP, file = paste0(outD,.Platform$file.sep,"NRS_phyto_genus_mat.csv"), row.names = FALSE)

#### Species Abund ####

# Check at what level we need change log
nrsls <- NRSPcl %>% 
  mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

NRSSpecP1 <- NRSPdat %>% 
  mutate(TaxonName = paste0(Genus, ' ', Species)) %>% # remove comments about with flagellates or cilliates etc.
  filter(!TaxonName %in% levels(as.factor(nrsls$TaxonName))
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("/", Species)) %>% 
  group_by(Sample, TaxonName) %>% 
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE), .groups = "drop")

NRSSpecP1 <- NRSPsamp %>% 
  left_join(NRSSpecP1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         TaxonName = ifelse(is.na(TaxonName), "Paralia sulcata", TaxonName),
         Cells_L = ifelse(is.na(Cells_L), 0, Cells_L)) %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(Cells_L = sum(Cells_L), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
NRSSpecP2 <- NRSPdat %>% 
  mutate(TaxonName = paste0(Genus, ' ', Species)) %>% # remove comments about with flagellates or cilliates etc.
  filter(TaxonName %in% levels(as.factor(nrsls$TaxonName))
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species)) %>% 
  left_join(NRSPcl, by = "TaxonName") %>%
  mutate(TaxonName = as_factor(TaxonName)) %>% 
  drop_na(TaxonName) %>%
  group_by(Sample, StartDate, TaxonName) %>% 
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(NRSSpecP2$TaxonName)) {
  Taxon = NRSSpecP2 %>% select(TaxonName) %>% unique() 
  Taxon = as.character(Taxon$TaxonName[i] %>% droplevels())
  Dates <- as.data.frame(NRSSpecP2) %>% 
    filter(TaxonName == Taxon) %>% 
    slice(1)  %>% 
    droplevels()
  
  spec <- as.data.frame(NRSSpecP2) %>% 
    filter(TaxonName == Taxon) %>% 
    droplevels() 
  
  spec <- NRSPsamp %>% 
    left_join(spec, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           TaxonName = replace(TaxonName, is.na(TaxonName), Dates$TaxonName),
           Cells_L = replace(Cells_L, StartDate>SampleDateLocal, -999), 
           Cells_L = replace(Cells_L, StartDate<SampleDateLocal & is.na(Cells_L), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, TaxonName) %>%
    summarise(Cells_L = sum(Cells_L), .groups = "drop") %>% 
    as.data.frame()     
  NRSSpecP1 <- rbind(NRSSpecP1, spec)
}

NRSSpecP1 <- NRSSpecP1 %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(Cells_L = max(Cells_L), .groups = "drop") %>% 
  arrange(-desc(TaxonName)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

NRSSpecP <-  NRSSpecP1 %>% 
  pivot_wider(names_from = TaxonName, values_from = Cells_L, values_fill = list(Cells_L = 0)) %>% 
  arrange(desc(SampleDateLocal))  %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(NRSSpecP, file = paste0(outD,.Platform$file.sep,"NRS_phyto_species_mat.csv"), row.names = FALSE)

###########################################################################################
#### Higher Trophic Groups BioV ####

NRSHTGPB1 <- NRSPdat %>% 
  group_by(Sample, TaxonGroup) %>% 
  summarise(BioV_um3_L = sum(Biovolume_uM3_L, na.rm = TRUE), .groups = "drop") %>%
  filter(!TaxonGroup %in% c("Other","Coccolithophore", "Diatom","Protozoa")) 

NRSHTGPB1 <- NRSPsamp %>% 
  left_join(NRSHTGPB1, by = "Sample") %>% 
  mutate(TaxonGroup = ifelse(is.na(TaxonGroup), "Ciliate", TaxonGroup),
         BioV_um3_L = ifelse(is.na(BioV_um3_L), 0, BioV_um3_L)) %>% 
  arrange(-desc(TaxonGroup))

NRSHTGPB <-  NRSHTGPB1 %>% 
  pivot_wider(names_from = TaxonGroup, values_from = BioV_um3_L, values_fill = list(BioV_um3_L = 0)) %>% 
  arrange(desc(SampleDateLocal)) %>% 
  select(-Sample)

fwrite(NRSHTGPB, file = paste0(outD,.Platform$file.sep,"NRS_phytoBioV_HTG_mat.csv"), row.names = FALSE)

#### Genus ####

# Check genus are effected by change log
nrslg <- NRSPcl %>% 
  mutate(genus1 = word(TaxonName, 1),
         genus2 = word(ParentName, 1)) %>%
  mutate(same = ifelse(genus1==genus2, "yes", "no")) %>%
  filter(same == "no")# no changes at genera level

# for non change log species
NRSGenPB1 <- NRSPdat %>% 
  filter(!TaxonName %in% levels(as.factor(nrslg$TaxonName))) %>% 
  group_by(Sample, Genus) %>% 
  summarise(BioV_um3_L = sum(Biovolume_uM3_L, na.rm = TRUE), .groups = "drop") %>% 
  drop_na(Genus) 

NRSGenPB1 <- NRSPsamp %>% 
  left_join(NRSGenPB1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Genus = ifelse(is.na(Genus), "Acanthoica", Genus),
         BioV_um3_L = ifelse(is.na(BioV_um3_L), 0, BioV_um3_L)) %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(BioV_um3_L = sum(BioV_um3_L), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
NRSGenPB2 <- NRSPdat %>% 
  filter(TaxonName %in% levels(as.factor(nrslg$TaxonName))) %>% 
  left_join(NRSPcl, by = "TaxonName") %>%
  mutate(Genus = as_factor(Genus)) %>% 
  drop_na(Genus) %>%
  group_by(Sample, StartDate, Genus) %>% 
  summarise(BioV_um3_L = sum(Biovolume_uM3_L, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(NRSGenPB2$Genus)) {
  Gen = NRSGenPB2 %>% select(Genus) %>% unique() 
  Gen = as.character(Gen$Genus[i] %>% droplevels())

  Dates <- as.data.frame(NRSGenPB2) %>% 
    filter(Genus == Gen) %>% 
    slice(1) %>% 
    droplevels()
  
  gen <- as.data.frame(NRSGenPB2) %>% 
    filter(Genus == Gen) %>% 
    droplevels() 
  
  gen <- NRSPsamp %>% 
    left_join(gen, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Genus = replace(Genus, is.na(Genus), Dates$Genus),
           BioV_um3_L = replace(BioV_um3_L, StartDate>SampleDateLocal, -999), 
           BioV_um3_L = replace(BioV_um3_L, StartDate<SampleDateLocal & is.na(BioV_um3_L), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
    summarise(BioV_um3_L = sum(BioV_um3_L), .groups = "drop") %>% 
    as.data.frame()     
  NRSGenPB1 <- rbind(NRSGenPB1, gen)
}

NRSGenPB1 <- NRSGenPB1 %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(BioV_um3_L = max(BioV_um3_L), .groups = "drop") %>% 
  arrange(-desc(Genus)) %>% 
  as.data.frame()

# select maximum value of duplicates, but leave -999 for all other occurrences as not regularly identified
NRSGenPB <-  NRSGenPB1 %>% 
  pivot_wider(names_from = Genus, values_from = BioV_um3_L, values_fill = list(BioV_um3_L = 0)) %>% 
  arrange(desc(SampleDateLocal))  %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(NRSGenPB, file = paste0(outD,.Platform$file.sep,"NRS_phytoBioV_genus_mat.csv"), row.names = FALSE)

#### Species ####

# Check at what level we need change log
nrsls <- NRSPcl %>% 
  mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

NRSPdat1 <- NRSPdat %>% 
  mutate(TaxonName = paste0(Genus, ' ', Species)) %>% # remove comments about with flagellates or ciliates etc.
  filter(!grepl("^(?!.*/~\\$)", TaxonName, perl = TRUE) &
                                Species == "spp." &
                               !TaxonName %in% levels(as.factor(nrsls$TaxonName)) )
NRSSpecPB1 <- NRSPdat %>% 
  filter(!TaxonName %in% levels(as.factor(nrsls$TaxonName)) &
         Species != "spp." & !is.na(Species) & !grepl("cf.", Species)) %>% 
  rbind(NRSPdat1) %>% 
  group_by(Sample, TaxonName) %>% 
  summarise(BioV_um3_L = sum(Biovolume_uM3_L, na.rm = TRUE), .groups = "drop")

NRSSpecPB1 <- NRSPsamp %>% 
  left_join(NRSSpecPB1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         TaxonName = ifelse(is.na(TaxonName), "Paralia sulcata", TaxonName),
         BioV_um3_L = ifelse(is.na(BioV_um3_L), 0, BioV_um3_L)) %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(BioV_um3_L = sum(BioV_um3_L), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
NRSPdat2 <- NRSPdat %>% 
  mutate(TaxonName = paste0(Genus, ' ', Species)) %>% # remove comments about with flagellates or ciliates etc.
  filter(!grepl("^(?!.*/~\\$)", TaxonName, perl = TRUE) &
                                 Species == "spp." &
                                 TaxonName %in% levels(as.factor(nrsls$TaxonName)))

NRSSpecPB2 <- NRSPdat %>% 
  filter(TaxonName %in% levels(as.factor(nrsls$TaxonName))
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species)) %>% 
  rbind(NRSPdat2) %>% 
  left_join(NRSPcl, by = "TaxonName") %>%
  mutate(TaxonName = as_factor(TaxonName)) %>% 
  drop_na(TaxonName) %>%
  group_by(Sample, StartDate, TaxonName) %>% 
  summarise(BioV_um3_L = sum(Cells_L, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(NRSSpecPB2$TaxonName)) {
  Taxon = NRSSpecPB2 %>% select(TaxonName) %>% unique() 
  Taxon = as.character(Taxon$TaxonName[i] %>% droplevels())

  Dates <- as.data.frame(NRSSpecPB2) %>% 
    filter(TaxonName == Taxon) %>% 
    slice(1)  %>% 
    droplevels()
  
  spec <- as.data.frame(NRSSpecPB2) %>% 
    filter(TaxonName == Taxon) %>% 
    droplevels() 
  
  spec <- NRSPsamp %>% 
    left_join(spec, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           TaxonName = replace(TaxonName, is.na(TaxonName), Dates$TaxonName),
           BioV_um3_L = replace(BioV_um3_L, StartDate>SampleDateLocal, -999), 
           BioV_um3_L = replace(BioV_um3_L, StartDate<SampleDateLocal & is.na(BioV_um3_L), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, TaxonName) %>%
    summarise(BioV_um3_L = sum(BioV_um3_L), .groups = "drop") %>% 
    as.data.frame()     
  NRSSpecPB1 <- rbind(NRSSpecPB1, spec)
}

NRSSpecPB1 <- NRSSpecPB1 %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(BioV_um3_L = max(BioV_um3_L), .groups = "drop") %>% 
  arrange(-desc(TaxonName)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

NRSSpecPB <-  NRSSpecPB1 %>% 
  pivot_wider(names_from = TaxonName, values_from = BioV_um3_L, values_fill = list(BioV_um3_L = 0)) %>% 
  arrange(desc(SampleDateLocal))  %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(NRSSpecPB, file = paste0(outD,.Platform$file.sep,"NRS_phytoBioV_species_mat.csv"), row.names = FALSE)

#### NRS Zooplankton #### #################################################################################################################################
# Bring in all NRS zooplankton samples, data and changelog
NRSZsamp <- getNRSZooSamples()
NRSZdat <- getNRSZooData()
NRSZcl <- getNRSZooChangeLog()

#### Raw Zooplankton ####

NRSRawZ1 <- left_join(NRSZsamp, NRSZdat, by = "Sample") %>% 
  select(-c(Sample, Copepod, TaxonGroup, Genus, Species)) %>% 
  arrange(-desc(TaxonName)) 

NRSRawZ <- NRSRawZ1 %>% 
  pivot_wider(names_from = TaxonName, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(NRSRawP, file = paste0(outD,.Platform$file.sep,"NRS_zoop_raw_mat.csv"), row.names = FALSE)

#### Higher Trophic Groups ####
nrsHTGZ1 <- NRSZdat %>% 
  group_by(Sample, TaxonGroup) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop") %>%
  filter(!TaxonGroup %in% c("Other")) 

nrsHTGZ1 <-  NRSZsamp %>% 
  left_join(nrsHTGZ1, by = "Sample") %>% 
  mutate(TaxonGroup = ifelse(is.na(TaxonGroup), "Copepod", TaxonGroup),
         ZAbund_m3 = ifelse(is.na(ZAbund_m3), 0, ZAbund_m3)) %>% arrange(-desc(TaxonGroup))

nrsHTGZ <-  nrsHTGZ1 %>% 
  pivot_wider(names_from = TaxonGroup, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) %>% 
  select(-Sample) %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(nrsHTGZ, file = paste0(outD,.Platform$file.sep,"NRS_zoop_HTG_mat.csv"), row.names = FALSE)

#### Genus ####

# Check genus are effected by change log
nrszlg <- NRSZcl %>% 
  mutate(genus1 = word(TaxonName, 1),
         genus2 = word(ParentName, 1)) %>%
  mutate(same = ifelse(genus1==genus2, "yes", "no")) %>%
  filter(same == "no")# no changes at genera level

# for non change log species
NRSGenZ1 <- NRSZdat %>% 
  filter(!TaxonName %in% levels(as.factor(nrszlg$TaxonName))) %>% 
  group_by(Sample, Genus) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop") %>% 
  drop_na(Genus) 

NRSGenZ1 <- NRSZsamp %>% left_join(NRSGenZ1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Genus = ifelse(is.na(Genus), "Acanthoica", word(Genus,1)), # bin subgenera together
         ZAbund_m3 = ifelse(is.na(ZAbund_m3), 0, ZAbund_m3))  %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
NRSGenP2 <- NRSZdat %>% 
  filter(TaxonName %in% levels(as.factor(nrszlg$TaxonName))) %>% 
  left_join(NRSZcl, by = "TaxonName") %>%
  mutate(Genus = as_factor(Genus)) %>% 
  drop_na(Genus) %>%
  group_by(Sample, StartDate, Genus) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(NRSGenP2$Genus)) {
  Gen <- NRSGenP2 %>% select(Genus) %>% unique()
  Gen = as.character(Gen$Genus[i] %>% droplevels()) 
  
  Dates <- as.data.frame(NRSGenP2) %>% 
    filter(Genus == Gen) %>% 
    slice(1) %>% 
    droplevels()
  
  gen <- as.data.frame(NRSGenP2) %>% 
    filter(Genus == Gen) %>% 
    droplevels() 
  
  gen <- NRSZsamp %>% 
    left_join(gen, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Genus = replace(Genus, is.na(Genus), Dates$Genus),
           ZAbund_m3 = replace(ZAbund_m3, StartDate>SampleDateLocal, -999), 
           ZAbund_m3 = replace(ZAbund_m3, StartDate<SampleDateLocal & is.na(ZAbund_m3), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
    summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
    as.data.frame()     
  NRSGenZ1 <- rbind(NRSGenZ1, gen)
}

NRSGenZ1 <- NRSGenZ1 %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(ZAbund_m3 = max(ZAbund_m3), .groups = "drop") %>% 
  arrange(-desc(Genus)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

NRSGenZ <-  NRSGenZ1 %>% 
  pivot_wider(names_from = Genus, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(NRSGenZ, file = paste0(outD,.Platform$file.sep,"NRS_zoop_genus_mat.csv"), row.names = FALSE)

#### Copepods ####

# Check at what level we need change log
nrsclc <- NRSZcl %>% 
  mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

NRSCop1 <- NRSZdat %>% 
  filter(!TaxonName %in% levels(as.factor(nrsclc$TaxonName)) & Copepod =="COPEPOD"
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  group_by(Sample, Species) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop")

NRSCop1 <- NRSZsamp %>% left_join(NRSCop1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Species = ifelse(is.na(Species), "Calanus Australis", Species), # avoids nulls in pivot
         ZAbund_m3 = ifelse(is.na(ZAbund_m3), 0, ZAbund_m3)) %>%  # avoids nulls in pivot
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
NRSCop2 <- NRSZdat %>% 
  filter(TaxonName %in% levels(as.factor(nrsclc$TaxonName)) & Copepod =="COPEPOD"
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  left_join(NRSZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>% 
  drop_na(Species) %>%
  group_by(Sample, StartDate, Species) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(NRSCop2$Species)) {
  Taxon <- NRSCop2 %>% select(Species) %>% unique()
  Taxon <- as.character(Taxon$Species[i] %>% droplevels())
  
  Dates <- as.data.frame(NRSCop2) %>% 
    filter(Species == Taxon) %>% 
    slice(1) %>% 
    droplevels()
  
  copes <- as.data.frame(NRSCop2) %>% 
    filter(Species == Taxon) %>% 
    droplevels() 
  
  copes <- NRSZsamp %>% 
    left_join(copes, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           ZAbund_m3 = replace(ZAbund_m3, StartDate>SampleDateLocal, -999), 
           ZAbund_m3 = replace(ZAbund_m3, StartDate<SampleDateLocal & is.na(ZAbund_m3), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
    summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
    as.data.frame()     
  NRSCop1 <- rbind(NRSCop1, copes)
}

NRSCop1 <- NRSCop1 %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbund_m3 = max(ZAbund_m3), .groups = "drop") %>% 
  arrange(-desc(Species)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

NRSCop <-  NRSCop1 %>% 
  pivot_wider(names_from = Species, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(NRSCop, file = paste0(outD,.Platform$file.sep,"NRS_zoop_copes_mat.csv"), row.names = FALSE)

#### Non-Copepods ####
# for non change log species

NRSnCop1 <- NRSZdat %>% 
  filter(!TaxonName %in% levels(as.factor(nrsclc$TaxonName)) & Copepod =="NON-COPEPOD"
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  group_by(Sample, Species) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop")

NRSnCop1 <- NRSZsamp %>% left_join(NRSnCop1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Species = ifelse(is.na(Species), "Calanus Australis", Species), # avoids nulls in pivot
         ZAbund_m3 = ifelse(is.na(ZAbund_m3), 0, ZAbund_m3)) %>%  # avoids nulls in pivot
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
NRSnCop2 <- NRSZdat %>% 
  filter(TaxonName %in% levels(as.factor(nrsclc$TaxonName)) & Copepod =="NON-COPEPOD"
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  left_join(NRSZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>%
  drop_na(Species) %>%
  group_by(Sample, StartDate, Species) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(NRSnCop2$Species)) {
  Taxon <- NRSnCop2 %>% select(Species) %>% unique()
  Taxon <- as.character(Taxon$Species[i] %>% droplevels())

  Dates <- as.data.frame(NRSnCop2) %>%
    filter(Species == Taxon) %>% 
    slice(1) %>% 
    droplevels()
  
  ncopes <- as.data.frame(NRSnCop2) %>% 
    filter(Species == Taxon) %>% 
    droplevels() 
  
  ncopes <- NRSZsamp %>% left_join(ncopes, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           ZAbund_m3 = replace(ZAbund_m3, StartDate>SampleDateLocal, -999), 
           ZAbund_m3 = replace(ZAbund_m3, StartDate<SampleDateLocal & is.na(ZAbund_m3), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
    summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
    as.data.frame()     
  NRSnCop1 <- rbind(NRSnCop1, ncopes)
}

NRSnCop1 <- NRSnCop1 %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbund_m3 = max(ZAbund_m3), .groups = "drop") %>% 
  arrange(-desc(Species)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurrences as not regularly identified

NRSnCop <- NRSnCop1 %>% 
  pivot_wider(names_from = Species, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) %>% 
  mutate(SampleDateLocal = as.character(SampleDateLocal))

fwrite(NRSnCop, file = paste0(outD,.Platform$file.sep,"NRS_zoop_noncopes_mat.csv"), row.names = FALSE)

# test <- read_csv(paste0(outD,.Platform$file.sep,"NRS_zoop_noncopes_mat.csv"))