# This is work in progress and the intellectual property of Jens Kattge and Susanne Tautenhahn.

# The routine is supposed to read TRY output, combine different TRY output files, clean taxonomy, 
# facilitate the exclusion of trait data along different aspects, incl. outliers, 
# and calculate species means and variations.  

#install.packages("tidyverse")
#install.packages("tseries")
#install.packages("hydroGOF")
#install.packages("epiDisplay")

library(epiDisplay)

library(tseries)
library(tidyverse)
#library(tidyr)
#library(dplyr)
#library(readr)

library(hydroGOF)

require(data.table)

library(ggplot2)
library(ggpubr)


#-------------------------------------------------
# read TRY data requests into data frames
#-------------------------------------------------

TRYdata1 <- fread("Input_Data/7571.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)

# check data frame: TRYdata1

dim(TRYdata1)
ls(TRYdata1) 
head(TRYdata1) 
#tail(TRYdata1) # avoid in case of large files
str(TRYdata1)

TRYdata2 <- fread("Input_Data/7956.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)

# ceck data frame: TRYdata2

dim(TRYdata2)
ls(TRYdata2) 
head(TRYdata2) 
#tail(TRYdata2) # avoid in case of large files
str(TRYdata2)

TRYdata3 <- fread("Input_Data/8200.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)

dim(TRYdata3)
ls(TRYdata3) 
head(TRYdata3) 
#tail(TRYdata3) # avoid in case of large files
str(TRYdata3)

# combine the two

TRYdata2 <- rbind(TRYdata2, TRYdata3)

dim(TRYdata2)
ls(TRYdata2) 
head(TRYdata2) 
#tail(TRYdata2) # avoid in case of large files
str(TRYdata2)


#-------------------------------------------------
# Read additional information
#-------------------------------------------------
# adding "genus" and "family" to enable BHPMF
# adding "plant growth form" to enable additional data selection

# Read additional Information from TRY: "species_level", "genus", "family" and "growthform" 
# Make sure, AccSpeciesID is unique, the file does not contain commas and all information on Genus and Family is complete (necessary for gap-filling)

AdditionalTRYdata.tmp <- fread("Input_Data/AccSpecies4DataPublication_unique.csv", header = T, sep = ",", dec = ".", quote = "", data.table = T)

dim(AdditionalTRYdata.tmp)
ls(AdditionalTRYdata.tmp) 
head(AdditionalTRYdata.tmp) 
tail(AdditionalTRYdata.tmp) 
str(AdditionalTRYdata.tmp)

AdditionalTRYdata <- AdditionalTRYdata.tmp %>%
	select(AccSpeciesID, SpeciesLevel, Genus, Family, GrowthForm)

dim(AdditionalTRYdata)
ls(AdditionalTRYdata) 
head(AdditionalTRYdata) 
tail(AdditionalTRYdata) 
str(AdditionalTRYdata)

rm(AdditionalTRYdata.tmp)

#-------------------------------------------------
# view traits and auxilliary Data
#-------------------------------------------------

# for small files:
#DataIDs <- TRYdata1 %>%
#  group_by(DataID, DataName, TraitID, TraitName) %>%
#  summarize(n = n()) %>%
#  arrange(DataName)%>%
#  write.csv(file = "Traits_auxilliary_Data_incl_Names.csv")

# for large requests: use IDs only, avoid names (link to names externally)
DataIDs <- TRYdata1 %>%
  group_by(DataID, TraitID) %>%
  summarize(n = n()) %>%
  arrange(TraitID)%>%
  write.csv(file = "Traits_auxilliary_Data.csv")

#-------------------------------------------------
# check number of trait records
#-------------------------------------------------

check_traits <- subset(TRYdata1, TraitID > 0, select=ObservationID:Comment)
dim(check_traits)
rm(check_traits)
check_trait_duplicates <- subset(TRYdata1, OrigObsDataID > 0, select=ObservationID:Comment)
dim(check_trait_duplicates)
rm(check_trait_duplicates)

gc()

#-------------------------------------------------
# For large requests: Reduce the size of the data frame by selecting only relevant columns and rows
#-------------------------------------------------

# select relevant columns 

work_data1 <- TRYdata1 %>%
	select(ObsDataID, ObservationID, AccSpeciesID, AccSpeciesName, ValueKindName, TraitID, TraitName, DataID, DataName, OriglName, OrigValueStr, OrigUnitStr, StdValue, UnitName, OrigObsDataID, ErrorRisk, Comment)

dim(work_data1)

# select relevant rows: all traits and only relevant auxilliary data

# 59 Latitude
# 60 Longitude
# 61 Altitude
# 6601 Sampling date
# 210 Leaf exposition
# 308 Experimental treatment
# 327 Exposition
# 413 Plant developmental status / plant age / maturity / plant life stage
# 1961 Health status of plants (vitality)

work_data1 <- subset(work_data1, TraitID > 0 | DataID == 59 | DataID == 60 | DataID == 61 | DataID == 6601 | DataID == 327 | DataID == 413 | DataID == 1961 | DataID == 210 | DataID == 308, select=ObsDataID:Comment)

dim(work_data1)


# do the same for the second TRY request

work_data2 <- TRYdata2 %>%
	select(ObsDataID, ObservationID, AccSpeciesID, AccSpeciesName, ValueKindName, TraitID, TraitName, DataID, DataName, OriglName, OrigValueStr, OrigUnitStr, StdValue, UnitName, OrigObsDataID, ErrorRisk, Comment)

work_data2 <- subset(work_data2, TraitID > 0 | DataID == 59 | DataID == 60 | DataID == 61 | DataID == 6601 | DataID == 327 | DataID == 413 | DataID == 1961 | DataID == 210 | DataID == 308, select=ObsDataID:Comment)

dim(work_data2)

#combine both:

work_data <- rbind(work_data1, work_data2)

dim(work_data)

nc <- ncol(work_data)

nc

#-------------------------------------------------
# write traits and auxilliary Data in csv file for visualization
#-------------------------------------------------

DataIDs <- work_data %>%
  group_by(DataID, DataName, TraitID, TraitName) %>%
  summarize(n = n()) %>%
  arrange(DataName)%>%
  write.csv(file = "Traits_auxilliary_Data.csv")

#-------------------------------------------------
# write traits Units in csv file for visualization
#-------------------------------------------------

TraitIDs <- work_data %>%
  group_by(TraitID, TraitName, UnitName) %>%
  summarize(n = n()) %>%
  arrange(TraitName)%>%
  write.csv(file = "Traits_Units.csv")

#-------------------------------------------------
# for small files, when a reduction of the size of the data frame is not necessary:
#-------------------------------------------------

# work_data <- TRYdata 

#-------------------------------------------------
# Add additional information
#-------------------------------------------------

work_data <- left_join(work_data, AdditionalTRYdata, by = "AccSpeciesID")

head(work_data)

work_data <- work_data[, c(1:4, (nc+1),(nc+2),(nc+3),(nc+4), 5:nc)]

#-------------------------------------------------
# check data.frame: work_data
#-------------------------------------------------

dim(work_data)
head(work_data)
tail(work_data)
str(work_data)

nc <- ncol(work_data)

nc

#--------------------------------------------------
# standardize species to format "genus species"
#--------------------------------------------------

# check for string lengths (number of words in AccSpeciesName)
species_names <- work_data %>%
  distinct(AccSpeciesName) %>%
  arrange(AccSpeciesName) %>%
  mutate(word_length = str_count(AccSpeciesName, "\\w+"))
  
# check for different cases in word_length > 2
species_names %>%
  filter(word_length == 5)
  
#  (species_names)

    # "var.", "subsp." - these can be reduced to "genus species"
    # species with authors (3 words) - these can be reduced to "genus species"
    # " x " - these are hybrids - pattern has to be preserved
    # "Abies borisii-regis" - species with compund species name - these have to be preserved (use blank space as word separator to keep these)

species_names <- species_names[complete.cases(species_names$AccSpeciesName), ]
species_names <- species_names %>%
  select(AccSpeciesName) 

# find hybrids and keep them, reduce the rest of the species names to "genus species"
species_names$Species <- ""


for (i in 1:nrow(species_names)) {
  if (str_detect(species_names$AccSpeciesName[i], " x ")) {
    species_names$Species[i] <- word(species_names$AccSpeciesName[i], 1, 3, sep = " ")
  } else {
    species_names$Species[i] <- word(species_names$AccSpeciesName[i], 1, 2, sep = " ")
  }
}
rm(i)

# combine species_names with work_data
work_data <- left_join(work_data, species_names, by = "AccSpeciesName")
names(work_data)
work_data <- work_data[, c(1:5, (nc+1), 6:nc)]
rm(species_names)

gc()

head(work_data)

dim(work_data)


#-------------------------------------------------------------
# save work_data_unfiltered as back-up
#-------------------------------------------------------------

work_data_unfiltered <- work_data

#-------------------------------------------------------------
# filter juveniles: exclude the whole Observation
#-------------------------------------------------------------
# check which DataID contains information about juveniles

        # 413 - Plant developmental status / plant age / maturity / plant life stage

# check different states of OriglValueStr
work_data %>%  
  filter(DataID == 413) %>%
  distinct(OriglName, OrigValueStr, OrigUnitStr, StdValue, Comment) %>%
  arrange(OrigValueStr) %>%
	write.csv(file = "Maturity.csv")

exclude <- work_data %>%
  filter(OrigValueStr == "immature"|
           OrigValueStr == "junvenil"|
           OrigValueStr == "juvenil (3 years)"|
           OrigValueStr == "juvenile"|
           OrigValueStr == "Juvenile"|
           OrigValueStr == "juvenile, 11-14 weeks"|
           OrigValueStr == "juvenile, 6 weeks"|
           OrigValueStr == "juveniles"|
           OrigValueStr == "sapling"|
           OrigValueStr == "Sapling (1 - 10 y)"|
           OrigValueStr == "saplings"|
           OrigValueStr == "seedling"|
           OrigValueStr == "Seedling"|
           OrigValueStr == "Seedling (0 - 1 y)"|
           OrigValueStr == "seedlings"|
           OrigValueStr == "seedlings, < 1/2 year"|
           OriglName == "Developmental stage" & OrigValueStr == "S"|
           OriglName == "Seedlings (True/False)" & OrigValueStr == "T"|
           OriglName == "Juvenile" & OrigValueStr == "Y")

#view(exclude)

exclude <- unique(exclude$ObservationID)
work_data$exclude <- work_data$ObservationID %in% exclude
work_data <- work_data[work_data$exclude == F, -(nc+2)]

rm(exclude)

dim(work_data)

#---------------------------
# filter: non-natural exposition: Exclude the whole Observation
#---------------------------

# check different states of OriglValueStr
work_data %>% 
  filter(DataID == 327 | DataID == 308 | DataID == 210) %>%
  distinct(OriglName, OrigValueStr, OrigUnitStr, StdValue, Comment) %>%
  arrange(OrigValueStr) %>%
	write.csv(file = "Exposition.csv")

exclude <- work_data %>%
  filter(OrigValueStr == "Climate Chamber"|
           OrigValueStr == "Climate chamber, non-limiting conditions, (cf. dataset reference)"|
           OrigValueStr == "climate chambers"|
           OrigValueStr == "controlled environment room"|
           OrigValueStr == "experimental"|
           OrigValueStr == "experimental treatment"|
           OrigValueStr == "GH"|
           OrigValueStr == "Glass house"|
           OrigValueStr == "Glasshouse"|
           OrigValueStr == "Glasshouse experiment"|
           OrigValueStr == "Greehouse"|
           OrigValueStr == "Green house"|
           OrigValueStr == "greenhouse"|
           OrigValueStr == "Greenhouse"|
           OrigValueStr == "Greenhouse plants"|
           OrigValueStr == "Greenhouse, grrowth container"|
           OrigValueStr == "Greenhouse, Indiana University"|
           OrigValueStr == "Greenhouse: highlight_highpH_competition"|
           OrigValueStr == "Greenhouse: highlight_highpH_nocompetition"|
           OrigValueStr == "Greenhouse: highlight_lowpH_competition"|
           OrigValueStr == "Greenhouse: highlight_lowpH_nocompetition"|
           OrigValueStr == "Greenhouse: lowleight_lowpH_competition"|
           OrigValueStr == "Greenhouse: lowlight_highpH_competition"|
           OrigValueStr == "Greenhouse: lowlight_highpH_nocompetition"|
           OrigValueStr == "Greenhouse: lowlight_lowpH_nocompetition"|
           OrigValueStr == "groth chamber"|
           OrigValueStr == "growth chamber"|
           OrigValueStr == "Growth chamber"|
           OrigValueStr == "Growth Chamber"|
           OrigValueStr == "Growth chamber, -N"|
           OrigValueStr == "Growth chamber, +N"|
           OrigValueStr == "growth chambers"|
           OrigValueStr == "growth_chamber"|
           OrigValueStr == "mesocosm"|
           OrigValueStr == "Shade - Natural environment"|
           OriglName == "Artificial growth conditions (G=greenhouse, C=growth chamber)" & OrigValueStr == "G"|
           OriglName == "growingCondition" & OrigValueStr == "GH"|
           OriglName == "Natural / Greenhouse" & OrigValueStr == "G")

exclude <- unique(exclude$ObservationID)
work_data$exclude <- work_data$ObservationID %in% exclude
work_data <- work_data[work_data$exclude == F, -(nc+2)]

rm(exclude)

dim(work_data)

#---------------------------
# filter: non healthy plants: Exclude the whole Observation
#---------------------------

work_data %>% 
  filter(DataID == 1961) %>%
  distinct(OriglName, OrigValueStr, OrigUnitStr, StdValue, Comment) %>%
  arrange(OrigValueStr) %>%
	write.csv(file = "Healthy.csv")

exclude <- work_data %>%
  filter(OrigValueStr == "Dead"|
          OriglName == "dummy" & OrigValueStr == "dummy")

exclude <- unique(exclude$ObservationID)
work_data$exclude <- work_data$ObservationID %in% exclude
work_data <- work_data[work_data$exclude == F, -(nc+2)]

rm(exclude)

dim(work_data)

#---------------------------
# filter duplicates:
#---------------------------

exclude <- work_data %>%
  filter(OrigObsDataID >0)

exclude <- unique(exclude$OrigObsDataID)
work_data$exclude <- work_data$OrigObsDataID %in% exclude
work_data <- work_data[work_data$exclude == F, -(nc+2)]

dim(work_data)

#---------------------------
# filter outliers
#---------------------------

# select all auxiliary data and only trait data with with ErrorRisk < 4 (< 3),

work_data <- subset(work_data, ErrorRisk < 4 | DataID == 59 | DataID == 60 | DataID == 61 | DataID == 6601 | DataID == 327 | DataID == 413 | DataID == 1961 | DataID == 210 | DataID == 308)

dim(work_data)

#-----------------------------------------
# filter non-representative Trait DataIDs
#-----------------------------------------

# DataID	DataName												TraitID	TraitName
# 2526	Leaf senescent carbon content per dry mass				13	Leaf carbon (C) content per leaf dry mass
# 1628	Leaf nitrogen content per dry mass (shaded leaves)		14	Leaf nitrogen (N) content per leaf dry mass
# 2527	Leaf senescent nitrogen (N) content per dry mass		14	Leaf nitrogen (N) content per leaf dry mass
# 2528	Leaf senescent phosphorus (P) content per dry mass		15	Leaf phosphorus (P) content per leaf dry mass
# 1083	Diameter at base										21	Stem diameter
# 1909	Stem diameter 10 cm above soil surface					21	Stem diameter
# 1908	Stem diameter at base (basal diameter)					21	Stem diameter
# 1084	Stem diameter at base of crown (seedlings)				21	Stem diameter
# 4152	Seed Curved Length										27	Seed length
# 4081	Leaf dry fraction										47	LDMC
# 542	Leaf nitrogen content per total area					50	Leaf nitrogen (N) content per leaf area
# 2534	Leaf senescent dry mass									55	Leaf dry mass (single leaf)
# 6594	SLA disc: mid-vein, petiole, rhachis excluded; shade	3086	SLA petiole, rhachis and midrib excluded
# 6592	SLA disc: mid-vein, petiole and rhachis excluded; sun	3086	SLA petiole, rhachis and midrib excluded
# 6493	SLA lamina less rhachis									3086	SLA petiole, rhachis and midrib excluded
# 2311	Height of seedlings										3106	Plant height vegetative
# 2081	Plant height of aquatic plants							3106	Plant height vegetative
# 6595	SLA lamina: petiole and rhachis excluded; shade			3115	SLA petiole excluded
# 6591	SLA leaf; shade											3116	SLA petiole included
# 4079	Whole plant SLA											3117	SLA undefined if petiole is in- or excluded

work_data <- work_data[!(work_data$DataID == "2526" | 
							work_data$DataID == "1628" |
							work_data$DataID == "2527" |
							work_data$DataID == "2528" |
							work_data$DataID == "1083" |
							work_data$DataID == "1909" |
							work_data$DataID == "1908" |
							work_data$DataID == "1084" |
							work_data$DataID == "4152" |
							work_data$DataID == "4081" |
							work_data$DataID == "542"  |
							work_data$DataID == "2534" |
							work_data$DataID == "6594" |
							work_data$DataID == "6592" |
							work_data$DataID == "6493" |
							work_data$DataID == "2311" |
							work_data$DataID == "2081" |
							work_data$DataID == "6595" |
							work_data$DataID == "6591" |
							work_data$DataID == "5089" |
							work_data$DataID == "5091" |
							work_data$DataID == "1923" |
							work_data$DataID == "1323" |
							work_data$DataID == "7149" |
							work_data$DataID == "4603" |
							work_data$DataID == "4605" |
							work_data$DataID == "4079" ),] 

dim(work_data)
#head(work_data)

#---------------------------
# Exclude ObsDataIDs with incorrect units
#---------------------------

exclude <- subset(work_data, 	(TraitID == 3116  & UnitName != "mm2 mg-1") |
								(TraitID == 144  & UnitName != "mm") |
								(TraitID == 3107  & UnitName != "m") |
								(TraitID == 21  & UnitName != "m") |
								(TraitID == 4  & UnitName != "g/cm3"))
dim(exclude)

exclude <- unique(exclude$ObsDataID)
work_data$exclude <- work_data$ObsDataID %in% exclude
work_data <- work_data[work_data$exclude == F, -(nc+2)]

rm(exclude)
dim(work_data)

#---------------------------
# delete plant height for trees < 10m, herbs > 10m, shrubs > 10m
#---------------------------

exclude <- subset(work_data, ((TraitID == 3106 | TraitID == 3107) & ((GrowthForm == "tree" & StdValue < 2) | (GrowthForm == "shrub" & StdValue >15) | (GrowthForm == "herb" & StdValue > 10))))

dim(exclude)
write.csv(exclude, file = "Height_excluded.csv")

exclude <- unique(exclude$ObsDataID)
work_data$exclude <- work_data$ObsDataID %in% exclude
work_data <- work_data[work_data$exclude == F, -(nc+2)]

rm(exclude)
dim(work_data)


#---------------------------------------------------------------------------------------
# select only trait values with numerical values
#---------------------------------------------------------------------------------------

check_traits <- subset(work_data, TraitID > 0, select=ObservationID:Comment)
dim(check_traits)

traits <- work_data %>%
  select(ObservationID, AccSpeciesID, AccSpeciesName, Species, Genus, Family, GrowthForm, TraitID, TraitName, StdValue, UnitName, ErrorRisk) %>%
  filter(complete.cases(TraitID) &      # exclude all entries with "" in TraitID
           complete.cases(StdValue))    # exclude potential categorical traits that don't have a StdValue, excludes also traits that have not yet been standardized in TRY

dim(traits)

#---------------------------------------------------------------------------------------
# select only species names of the form: genus and species epithet, remove NA (?)
#---------------------------------------------------------------------------------------

traits <- traits %>%
  filter(complete.cases(Species))

dim(traits)

#---------------------------------------------------------------------------------------
# transform data into long table format on traits, including latitude and longitude (additional auxiliary data can be added)
# ObservationID, AccSpeciesID, AccSpeciesName, Species, Genus, Family, GrowthForm, TraitID, TraitName, StdValue, Unit, Latitude, Longitude
#---------------------------------------------------------------------------------------
# !!! Here I use the preliminary output from filtered data (filtered by error risk, duplicates , maturity, experimental conditions, health)
# to prepare Data for calculation of species mean traits and gap-filling

Latitude <- work_data %>%
  select(ObservationID, DataName, StdValue) %>%
  filter(DataName == "Latitude") %>%
  rename(Latitude = StdValue) %>%
  select(ObservationID, Latitude) %>%
  distinct(ObservationID, .keep_all = T)

Longitude <- work_data %>%
  select(ObservationID, DataName, StdValue) %>%
  filter(DataName == "Longitude") %>%
  rename(Longitude = StdValue) %>%
  select(ObservationID, Longitude) %>%
  distinct(ObservationID, .keep_all = T)

traits_x_georef <- left_join(traits, Latitude, by = "ObservationID")
traits_x_georef <- left_join(traits_x_georef, Longitude, by = "ObservationID")

dim(traits_x_georef)
head(traits_x_georef)
tail(traits_x_georef)
str(traits_x_georef)

#---------------------------
# check trait ranges and frequency distributions
#---------------------------

# check 10 maximum values and 10 minimum values for lat, lon and alt 

#temp <- filter(traits_x_georef, DataID == "59")
#temp <- filter(traits_x_georef, DataID == "60")
#temp <- filter(traits_x_georef, DataID == "61")

#temp1 <- filter(temp, StdValue != "NA")
#temp2 <- arrange(temp1, StdValue)
#head(temp2, 10)
#tail(temp2, 10)

# check 10 maximum values and 10 minimum values for traits and plot histogram 

temp <- filter(work_data_unfiltered, TraitID == 3106 | TraitID == 3107)
temp1 <- filter(temp, StdValue != "NA")
temp2 <- arrange(temp1, StdValue)
head(temp2, 10)
tail(temp2, 10)
temp3_unfiltered <- data.frame(temp2$StdValue)
temp3_unfiltered$filter<-"unfiltered"

temp1 <- filter(traits_x_georef, TraitID == 3106 | TraitID == 3107)
temp2 <- arrange(temp1, StdValue)
head(temp2, 10)
tail(temp2, 10)
temp3 <- data.frame(temp2$StdValue)
temp3$filter<-"filtered"

temp4 <- rbind(temp3_unfiltered,temp3)

#head(temp4)

ggplot(temp4, aes(temp2.StdValue, fill = filter)) +
   geom_histogram(alpha = 0.5, position = 'identity', colour=("grey10")) +
   scale_x_log10() +
#   scale_x_continuous(limits = c(0,1)) +
   scale_fill_manual(values=c("grey10", "white"))

#---------------------------
# compare original and filteredl number of trait records
#---------------------------

check_traits_unfiltered <- subset(work_data_unfiltered, TraitID > 0, select=ObservationID:Comment)
dim(check_traits_unfiltered)
rm(check_traits_unfiltered)

check_traits_final <- subset(traits_x_georef, TraitID > 0, select=ObservationID:Longitude)
dim(check_traits_final)
rm(check_traits_final)

Final_number_of_AccSpecies <- distinct(traits_x_georef, AccSpeciesID)
dim(Final_number_of_AccSpecies)
rm(Final_number_of_AccSpecies)

gc()

#---------------------------------------------------------------------------------------
# write filtered trait data in long table format to csv  
#---------------------------------------------------------------------------------------

write.csv(traits_x_georef, file = "traits_x_georef_long_table.csv")

rm(traits, Latitude, Longitude)

#---------------------------------------------------------------------------------------
# transform data from long table format to wide table format: 
# using pivot_wider : https://tidyr.tidyverse.org/reference/pivot_wider.html
#---------------------------------------------------------------------------------------

# using pivot_wider on TraitID or TraitName can cause duplicate entries on ObservationID. In these cases the mean is calculated for StdValue and ErrorRisk. If we want to avoid this, we have to run Pivot on the DataIDs of the traits.

traits_x_georef %>%
  pivot_wider(names_from = c(TraitID, TraitName, UnitName), values_from = c(StdValue, ErrorRisk), values_fn = list(StdValue = mean, ErrorRisk = mean)) %>%
write.csv(file = "traits_x_georef_wide_table.csv")

#---------------------------------------------------------------------------------------
# Combine (sub)traits to broader traits
#---------------------------------------------------------------------------------------

# Plant height: 3106 3107 >> 18

traits_x_georef$TraitID[traits_x_georef$TraitID %in% 3106 ] <- 18
traits_x_georef$TraitID[traits_x_georef$TraitID %in% 3107 ] <- 18

# SLA

traits_x_georef$TraitID[traits_x_georef$TraitID %in% 3086 ] <- 11
traits_x_georef$TraitID[traits_x_georef$TraitID %in% 3115 ] <- 11
traits_x_georef$TraitID[traits_x_georef$TraitID %in% 3116 ] <- 11
traits_x_georef$TraitID[traits_x_georef$TraitID %in% 3117 ] <- 11

# Rooting depth

traits_x_georef$TraitID[traits_x_georef$TraitID %in% 1777 ] <- 6


# SRL

traits_x_georef$TraitID[traits_x_georef$TraitID %in% 2318 ] <- 1080
traits_x_georef$TraitID[traits_x_georef$TraitID %in%  614 ] <- 1080


# Leaf area

traits_x_georef$TraitID[traits_x_georef$TraitID %in%  3108 ] <- 3112
traits_x_georef$TraitID[traits_x_georef$TraitID %in%  3110 ] <- 3112

# Leaflet area

traits_x_georef$TraitID[traits_x_georef$TraitID %in%  3109 ] <- 3113
traits_x_georef$TraitID[traits_x_georef$TraitID %in%  3111 ] <- 3113


#---------------------------------------------------------------------------------------
#Show traits
#---------------------------------------------------------------------------------------

Traits_unique <- unique(traits_x_georef$TraitID)
Traits_unique <- sort(Traits_unique)
Traits_unique

#    4    6   11   13   14   15   18   21   26   27   46   47   50   55   78   95  138  144  145  146  163  169  223  224  237  281  282  289 1080 3112 3113 3114 3120

#---------------------------------------------------------------------------------------
#Number of species
#---------------------------------------------------------------------------------------

Species_unique <- unique(traits_x_georef$Species)
length(Species_unique)

#---------------------------------------------------------------------------------------
# aggregate at species level 
#---------------------------------------------------------------------------------------

head(traits_x_georef)
dim(traits_x_georef) 

detach(traits_x_georef)

Species.data <- aggregate(traits_x_georef$StdValue, by = list(Species=traits_x_georef$Species, TraitID=traits_x_georef$TraitID), FUN=c("count","mean","median","sd","se","min","max"), na.rm=TRUE)

head(Species.data)

write.csv(Species.data, file = "Species_mean_traits.csv")

Species_number_per_trait <- aggregate(Species.data$mean.traits_x_georef.StdValue, by = list(Species.data$TraitID), FUN=c("count"), na.rm=TRUE)

write.csv(Species_number_per_trait, file = "Number_of_species_per_trait.csv")

#---------------------------------------------------------------------------------------
# Transform to wide table:
#---------------------------------------------------------------------------------------

ls(traits_x_georef)

traits_x_georef_select <- select(traits_x_georef, AccSpeciesID, AccSpeciesName, Species, Genus, Family, GrowthForm, Latitude, Longitude, ObservationID, StdValue, TraitID)

ls(traits_x_georef_select)


wide_table <- traits_x_georef_select %>%
  pivot_wider(names_from = c(TraitID), values_from = c(StdValue), values_fn = list(StdValue = mean)) 

write.csv(wide_table, file = "wide_table.csv")

#---------------------------------------------------------------------------------------
# view trait histograms before and after filtering 
#---------------------------------------------------------------------------------------

for(i in c(4, 6, 11, 13, 14, 15, 18, 21, 26, 27, 46, 47, 50, 55, 78, 95, 138, 144, 145, 146, 163, 169, 223, 224, 237, 281, 282, 289, 1080, 3112, 3113, 3114, 3120)){

#for(i in c(1080)){

temp <- filter(work_data_unfiltered, TraitID == i)

if(i == 6){   temp <- filter(work_data_unfiltered, TraitID == 1777 | TraitID == 6 )}
if(i == 11){  temp <- filter(work_data_unfiltered, TraitID == 3086 | TraitID == 3115 | TraitID == 3116 | TraitID == 3117 | TraitID == 11)}
if(i == 18){  temp <- filter(work_data_unfiltered, TraitID == 3106 | TraitID == 3107 | TraitID == 18 )}
if(i == 1080){temp <- filter(work_data_unfiltered, TraitID == 2318 | TraitID == 614 | TraitID == 1080)}
if(i == 3112){temp <- filter(work_data_unfiltered, TraitID == 3108 | TraitID == 3110 | TraitID == 3112)}
if(i == 3113){temp <- filter(work_data_unfiltered, TraitID == 3109 | TraitID == 3111 | TraitID == 3113)}
if(i == 3114){temp <- filter(work_data_unfiltered, TraitID == 3109 | TraitID == 3111 | TraitID == 3113 | TraitID == 3108 | TraitID == 3110 | TraitID == 3112 | TraitID == 3114)}

temp1 <- filter(temp, StdValue != "NA")
temp2 <- arrange(temp1, StdValue)
temp3_unfiltered <- data.frame(temp2$StdValue)
temp3_unfiltered$filter<-"0_unfiltered"

temp1 <- filter(traits_x_georef, TraitID == i)
if(i == 3114){temp1 <- filter(traits_x_georef, TraitID == 3112 | TraitID == 3113 | TraitID == 3114)}
temp2 <- arrange(temp1, StdValue)
temp3 <- data.frame(temp2$StdValue)
temp3$filter<-"1_filtered"

temp4 <- rbind(temp3_unfiltered,temp3)

if(i == 4){
p4 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.05, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x= bquote('Stem specific density '~(g/cm^3)), y="Number of observations") +
  scale_x_continuous(breaks=c(0.0, 0.5, 1), labels=c("0.0", "0.5", "1.0"), limits=c(-0.1, 1.5)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p4)
}

if(i == 6){
p6 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.15, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x= bquote('Rooting depth '~(m)), y="") +
  scale_x_log10(breaks=c(0.01,0.1,1,10,100), limits=c(0.002,300), labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p6)
}

if(i == 11){
p11 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.09, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x= bquote("SLA" ~(mm^2 ~mg^-1)), y="") +
  scale_x_log10(breaks=c( 1, 10, 100), labels=c( "1", "10","100"), limits=c(0.3, 200)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p11)
}

if(i == 13){
p13 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=15, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf Cmass" ~(mg ~g^-1)), y="Number of observations") +
  scale_x_continuous(limits=c(200, 700)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p13)
}

if(i == 14){
p14 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.05, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf Nmass" ~(mg ~g^-1)), y="") +
  scale_x_log10(breaks=c( 10, 100), labels=c( "10","100"), limits=c(2, 100)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p14)
}

if(i == 15){
p15 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.075, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf Pmass" ~(mg ~g^-1)), y="Number of observations") +
  scale_x_log10(breaks=c(0.1, 1, 10), labels=c("0.1", "1", "10"), limits=c(0.05, 15)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p15)
}

if(i == 18){
p18 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.21, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Plant height" ~(m)), y="Number of observations") +
  scale_x_log10(breaks=c(0.01, 0.1, 1, 10, 100), labels=c("0.01", "0.1", "1", "10","100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p18)
}

if(i == 21){
p21 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.21, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Stem diameter" ~(m)), y="") +
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1, 10), labels=c("0.001","0.01", "0.1", "1", "10")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p21)
}

if(i == 26){
p26 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.35, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Seed mass" ~(mg)), y="") +
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1,10,100, 1000,10000,100000), labels=c("0.001","", "0.1","", "10","", "1000","","100,000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p26)
}

if(i == 27){
p27 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.1, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Seed length" ~(mm)), y="Number of observations") +
  scale_x_log10(breaks=c(0.1, 1,10,100), labels=c("0.1","1", "10","100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p27)
}

if(i == 46){
p46 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.1, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf thickness" ~(mm)), y="") +
  scale_x_log10(breaks=c(0.01, 0.1, 1, 10), labels=c("0.01","0.1","1","10")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p46)
}

if(i == 47){
p47 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=0.025, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("LDMC" ~(g ~g^-1)), y="") +
  scale_x_continuous(limits=c(0.,1)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p47)
}

if(i == 50){
p50 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.05, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf Narea" ~(g ~m^-2)), y="Number of observations") +
  scale_x_log10(breaks=c(0.1, 1, 10), labels=c("0.1", "1", "10"), limits=c(0.1, 15)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p50)
}

if(i == 55){
p55 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.2, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf dry mass" ~(mg)), y="") +
  scale_x_log10(breaks=c(0.1, 1, 10, 100, 1000, 10000), labels=c("0.1", "1", "10", "100", "", "100000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p55)
}


if(i == 78){
p78 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=1, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x="Leaf 15N", y="") +
  scale_x_continuous() +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p78)
}


if(i == 95){
p95 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=5, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Seed germination efficiency" ~("%")), y="Number of observations") +
  scale_x_continuous(breaks=c(0, 25, 50, 75, 100), labels=c("0", "25", "50", "75","100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p95)
}

if(i == 138){
p138 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.4, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Seed number per reproducton unit"), y="") +
  scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000, 1000000), labels=c("1", "", "100", "", "10000", "", "1000000", "")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p138)
}

if(i == 144){
p144 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.15, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf length" ~("cm")), y="") +
  scale_x_log10(breaks=c(1, 10, 100, 1000, 10000), labels=c("1", "10", "100", "1000", "10000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p144)
}

if(i == 145){
p145 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.2, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf width" ~("cm")), y="Number of observations") +
  scale_x_log10(breaks=c(0.01, 0.1, 1, 10, 100), labels=c("0.01", "0.1", "1", "10", "100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p145)
}

if(i == 146){
p146 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=3, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf C/N ratio" ~(g ~g^-1)), y="") +
  # scale_x_log10(breaks=c(0.01, 0.1, 1), labels=c("0.01", "0.1", "1"), limits=c(0.002, 2)) +
  scale_x_continuous(limits=c(0, 125)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p146)
}

if(i == 163){
p163 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.25, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf fresh mass" ~(g)), y="") +
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1, 10, 100), labels=c("0.001", "0.01", "0.1", "1", "10", "100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p163)
}

if(i == 169){
p169 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.2, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Stem conduit density" ~(mm^-2)), y="Number of observations") +
  scale_x_log10(breaks=c(1, 10, 100, 1000, 10000), labels=c("1", "10", "100", "1000", "10000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p169)
}

if(i == 223){
p223 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.125, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Chromosome number"), y="") +
  scale_x_log10(breaks=c(10, 50, 100, 500), labels=c("10", "50","100", "500")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p223)
}

if(i == 224){
p224 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.125, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Chromosome cDNA content" ~(pg)), y="") +
  scale_x_log10(breaks=c(0.1, 1, 10, 100), labels=c("0.1", "1", "10","100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p224)
}

if(i == 237){
p237 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.125, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Dispersal unit length" ~(mm)), y="Number of observations") +
  scale_x_log10(breaks=c(0.1, 1, 10, 100), labels=c("0.1", "1", "10","100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p237)
}

if(i == 281){
p281 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.125, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=~"Stem conduit diameter" ~(mu * m), y="") +
  scale_x_log10(breaks=c(0.1, 1, 10, 100), labels=c("0.1", "1", "10","100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p281)
}

if(i == 282){
p282 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.085, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=~"Stem conduit element length" ~(mu * m), y="") +
  scale_x_log10(breaks=c(100, 500, 1000), labels=c("100", "500","1000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p282)
}

if(i == 289){
p289 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.085, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=~"Wood fiber lengths" ~(mu * m), y="Number of observations") +
  scale_x_log10(breaks=c(100, 500, 1000, 5000, 10000), labels=c("100", "500","1000", "5000","10000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p289)
}

if(i == 1080){
p1080 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.15, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=~"Specific root length, SRL" ~(cm ~g^-1), y="") +
  scale_x_log10(breaks=c(100, 1000, 10000, 100000), labels=c("100", "1000", "10000", "100000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p1080)
}

if(i == 3120){
p3120 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.1, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf water content" ~(g ~g^-1)), y="") +
  scale_x_log10(breaks=c(1,10,100), labels=c("1", "10","100")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p3120)
}

if(i == 3112){
p3112 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.3, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf area" ~(mm^2)), y="Number of observations") +
  scale_x_log10(breaks=c(1,10,100, 1000,  10000, 100000, 1000000), labels=c("1", "10","100", "1000", "10000", "", "1000000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p3112)
}

if(i == 3113){
p3113 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.3, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaflet area" ~(mm^2)), y="") +
  scale_x_log10(breaks=c(1,10,100, 1000,  10000, 100000, 1000000), labels=c("1", "10","100", "1000", "10000", "", "1000000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p3113)
}

if(i == 3114){
p3114 <- ggplot(temp4, aes(temp2.StdValue, fill = filter)) + 
  geom_histogram(binwidth=.3, alpha=1.0, position="identity", colour=("grey10")) +
  labs(x=bquote("Leaf or leaflet area" ~(mm^2)), y="") +
  scale_x_log10(breaks=c(1,10,100, 1000,  10000, 100000, 1000000), labels=c("1", "10","100", "1000", "10000", "", "1000000")) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("white","grey50")) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "")
#print(p3114)
}

}

ggarrange(p4, p6, p11, p13, p14, p15, p18, p21, p26, p27, p46, p47, p50, p55, p78,
          labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
          ncol = 3, nrow = 5)

ggsave("TraitHistograms_observed1.pdf", width = 30, height = 40, units = "cm", dpi=300)

ggarrange(p95, p138, p144, p145, p146, p163, p169, p223, p224, p237, p281, p282, p289, p1080, p3120,
          labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
          ncol = 3, nrow = 5)

ggsave("TraitHistograms_observed2.pdf", width = 30, height = 40, units = "cm", dpi=300)

ggarrange(p3112, p3113, p3114, 
          labels=c("A", "B", "C"),
          ncol = 3, nrow = 5)

ggsave("TraitHistograms_observed3.pdf", width = 30, height = 40, units = "cm", dpi=300)

