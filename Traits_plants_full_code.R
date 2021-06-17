# This script can be used to reproduce and adapt the plant trait dataset (Dataset ID XXX) provided on Bexis.
# It is structured in 3 main parts.
#       - the first part takes as input the species names list and the plant abundance dataset (ID XXX) to create a table linking BE species
#     names to TRY names, including synonyms and aggregate species data. Abundance data is used for plants identified only at genus level,
#     for which an abundance-weighted average of all species from the same genus is used.
#       - the second part takes as input the raw TRY data (dataset ID XXX), and cleans it following a procedure that will be detailed below.
#       - the last part combines the two dataset to obtain a BE species x trait dataset (corresponding to dataset XXX) and the corresponding
#     CWM values.

# Author: Margot Neyret - Please get in touch if you have questions or find some bugs!

setwd(dir = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Bexis_upload/')

### Note: this code was build under R version 3.7, some packages (e.g. tpl) don't work under R 4.0

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### ******* PART 0: LOADING LIBRARIES AND DATA ******* ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

library(beepr)
library(data.table)
library(dplyr)
library(FD)
library(ggpmisc)
library(readr)
library(readxl)
library(reshape2)
library(stringr)
library(taxize)
library(Taxonstand)
library(tpl)
library(vegan)

### Species list.
# Species_raw corresponds to raw species names found in the exploratories species lists; Species_list_format names with homogenized formatting, and Species_list_corrected includes
# the individual species names included in aggregate species
Species_table = fread('Species_list.csv')

### TRY raw data
TRYdata11741 <-
  fread(
    "Try_request_11741.txt",
    header = T,
    sep = "\t",
    dec = ".",
    quote = "",
    data.table = T
  ) # Get dataset 27586 on Bexis
TRYdata = TRYdata11741


### TRY species list:
# For updates (e.g. adding new species) use the following to download from https://www.try-db.org/dnld/TryAccSpecies.txt
TryAccSpecies <- fread("TryAccSpecies.txt")
#list_try_species = TryAccSpecies$AccSpeciesName
# But when I checked the list wasn't up to date (i.e. did not correspond to species names as used in the main dataset)
# leading to inconsistencies, so just in case I use both
list_try_species = unique(c(
  unique(TRYdata$AccSpeciesName),
  TryAccSpecies$AccSpeciesName
))


### Plant abundances
#Abundances_grasslands = fread("200205_EP_species_diversity_GRL.txt") # Get synthesis abundance dataset 25626 (or new version) on Bexis
# Species_info_grasslands = fread("200205_EP_species_info_GRL.txt") # Get synthesis info dataset 24608 (or new version) on Bexis

# New
Abundances_grasslands = fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/210112_EP_species_diversity_GRL_BEXIS.txt")
Species_info_grasslands = fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/210112_EP_species_info_GRL_BEXIS.txt") # Get synthesis info dataset 24608 (or new version) on Bexis

Abundances_grasslands = Abundances_grasslands[Species %in% Species_info_grasslands[Group_broad == 'Plant', unique(Species)], ]

Abundances_forests <-
  fread("25886.txt") # Get dataset on Bexis: 25886_Vegetation Records for 151 Forest EPs, 2009 - 2018_1.5.4


### Below-ground traits: merge 2 existing datasets
Root_traits_all = fread("26587.txt") # Get dataset 26587 on Bexis
Root_traits_myco = fread("26546.txt")[, c('species', 'col')] # Get dataset 26587 on Bexis
Root_traits_myco[, scientificName := recode(
  species,
  Veronica_teucrium = 'Veronica_austriaca',
  Tripleurospermum_perforatum = 'Tripleurospermum_inodorum',
  Hieracium_pilosella = 'Pilosella_officinarum',
  Leontodon_autumnalis = 'Scorzoneroides_autumnalis',
  Anthemis_tinctoria = 'Cota_tinctoria'
)]
Root_traits_myco[, c(
  'traitName',
  'traitID',
  'traitValue',
  'traitUnit',
  "scientificNameStd",
  "traitNameStd",
  "traitValueStd",
  "traitUnitStd",
  "measurementID"
) :=
  list(
    'Mycorrhizal_inf_int',
    '1030',
    col,
    '(unitless)',
    gsub('_', '', 'scientificName'),
    'Mycorrhizal infection intensity',
    col,
    '%',
    NA
  )]
Root_traits_myco = merge(Root_traits_myco,
                         unique(Root_traits_all[, c(
                           'scientificName',
                           'warnings',
                           'taxonID',
                           'taxonRank',
                           'kingdom',
                           'phylum',
                           'class',
                           'order',
                           'family'
                         )]),
                         by = 'scientificName',
                         all.x = T)

Root_traits = rbind(Root_traits_all, Root_traits_myco[, .SD, .SDcols = colnames(Root_traits_myco)[!(colnames(Root_traits_myco) %in% c('species', 'col'))]])

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### ******* PART 1: CREATING THE FULL SPECIES TABLE ******* ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

#  We first need to correct some species names to fit with TRY standards
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Neottia nidus avis', 'Neottia nidus-avis', x)
})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Impatiens noli tangere', 'Impatiens noli-tangere', x)
})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Buphtalmum salicifolium', 'Buphthalmum salicifolium', x)
})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Silene flos.cuculi', 'Silene flos-cuculi', x)
})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Athyrium filix femina', 'Athyrium filix-femina', x)
})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Capsella bursa.pastoris', 'Capsella bursa-pastoris', x)
})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Dryopteris filix mas', 'Dryopteris filix-mas', x)
})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Mycelis muralis', 'Lactuca muralis', x)
})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Persicaria lapathifolium', 'Persicaria lapathifolia', x)
})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub("Valerianella officinalis", 'Valeriana officinalis', x)
})] # I guess that was probably a typo
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub("heliscopia", 'helioscopia', x)
})] # I guess that was probably a typo
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub("Alisma plantago-aquatica",
       'Alisma plantago-aquatica subsp. orientale',
       x)
})] # Only subsp. present in TRY, better than nothing
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub(
    "Achillea distans",
    '"Achillea distans subsp. stricta___Achillea distans subsp. tanacetifolia"',
    x
  )
})] # Only subsp. present in TRY, better than nothing
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x) {
  gsub('Taraxacum sect ruderalia', 'Taraxacum officinale', x)
})]

Species_table[, Species_list_corrected := gsub('[\"]*', '', Species_list_corrected)]

# Remove unidentified species from the list
unidentifiable_species = c(
  "unknown deciduous",
  "Tree seedling",
  "Lamiaceae sp",
  "Asteraceae sp",
  "Caryophyllaceae sp",
  "Baumkeimling sp",
  'Orchidaceae sp',
  'Brassicaceae sp'
)
Species_table = Species_table[!(Species_list_format %in% unidentifiable_species), ]
Species_table = Species_table[Species_list_corrected != '', ]

#### The next step is to match species names to the TRY standard (mostly based on The Plant List) ####
# Two options:
#     1. If the current species names is found in the Try species list, then we keep only this accepted name. Otherwise,
#   we look for potential synonyms and retain only those present in the TRY species list.
#     2. (DEFAULT) We look for synonyms for all species names, regardless of whether the accepted name is in TRY. This options allows
#   to get data that might be lost otherwise (i.e. if TRY has more data for a synonym than for the accepted name.)
# To change this default, change the function name from 'get_all_synonyms_and_accepted2' to "get_all_synonyms_and_accepted1' on line ~134.

# Option 1: take synonyms only if the accepted name is not in TRY
get_all_synonyms_and_accepted1 = function(x, match_list) {
  x = unlist(strsplit(as.character(x), '___'))
  names_ids = sapply(x, function(y) {
    print(y)
    tpl_res = tpl.get(y, return.synonyms = T)
    id = NA
    spe_names = NA
    if (y %in% list_try_species) {
      spe_names = y
      id = tpl_res$all.entries$id
    }
    if (grepl('sp$', y)) {
      spe_names = y
      id = NA
    }
    else {
      if (tpl_res$all.entries$name %in% list_try_species) {
        spe_names = tpl_res$all.entries$name
        id = tpl_res$all.entries$id
      }
      else {
        synonyms = tpl_res$synonyms[tpl_res$synonyms$confidence.level %in% c('M', "H") &
                                      tpl_res$synonyms$name.synonym %in% list_try_species, ]$name.synonym
        spe_names = unique(synonyms)
        id = unique(tpl_res$synonyms[tpl_res$synonyms$confidence.level %in% c('M', "H") &
                                       tpl_res$synonyms$name.synonym %in% list_try_species, ]$id)
      }
    }
    return(list('spe_names' = as.character(spe_names), 'ids' = id))
  })
  names = unlist(names_ids['spe_names', ])
  names = unique(names[!is.na(names) & names != 'NA'])
  names = paste(names, collapse = '___')
  
  ids = names_ids['ids', ]
  ids = unique(ids[!is.na(ids) & ids != 'NA'])
  ids = paste(ids, collapse = '___')
  return(list('Synonyms' = names, 'IDs' = ids))
}

# Option 2: get ALL synonyms, even if the accepted name is in TRY
get_all_synonyms_and_accepted2 = function(x, match_list) {
  x = unlist(strsplit(as.character(x), '___'))
  names_ids = sapply(x, function(y) {
    tpl_res = tpl.get(y, return.synonyms = T)
    spe_names = c(y,
                  tpl_res$all.entries$name,
                  tpl_res$synonyms[tpl_res$synonyms$confidence.level %in% c('M', "H") &
                                     tpl_res$synonyms$name.synonym %in% list_try_species, ]$name.synonym)
    id = c(tpl_res$all.entries$id,
           tpl_res$synonyms[tpl_res$synonyms$confidence.level %in% c('M', "H") &
                              tpl_res$synonyms$name.synonym %in% list_try_species, ]$id)
    
    return(list('spe_names' = as.character(spe_names), 'ids' = id))
  })
  names = unlist(names_ids['spe_names', ])
  names = unique(names[!is.na(names) & names != 'NA'])
  names = paste(names, collapse = '___')
  
  ids = names_ids['ids', ]
  ids = unique(ids[!is.na(ids) & ids != 'NA'])
  ids = paste(ids, collapse = '___')
  return(list('Synonyms' = names, 'IDs' = ids))
}

Species_table$rowID = 1:nrow(Species_table)
Species_table[, c('Synonyms', 'IDs') := sapply(Species_list_corrected, function(x) {
  get_all_synonyms_and_accepted2(x, list_try_species)
}), by = rowID] #%***% Change get_all_synonyms_and_accepted2 to get_all_synonyms_and_accepted1 here for sensitivity analysis

beep() # This will warn you when the request is finished.


Species_table_2 = copy(Species_table)
Species_table_2$All_names = unlist(Species_table$Synonyms)


#### Some species are unresolved at the genus ####
# When other species have been identified in the same land-use type (forest or grasslands), we will attribute the average
# of the traits of those species, weighted by their relative abundance
# When NO other species have been identified in the same LU type, but some have in the other, we will use the values in the
# other LU
# Otherwise, we'll go through each case independently

# First I add a column, just for these problematic species, to differentiate between forest and grassland estimates
# and I duplicate the corresponding rows
Species_table_2$Forest_or_grassland_estimate = as.character(Species_table_2$Forest_or_grassland_estimate)
Genus_sp = Species_table_2[(
  grepl(" sp$", Species_list_format) |
    grepl("Rhinanthus$", Species_list_format) |
    grepl("Rhinanthus aggr.", Species_list_format)
) & !grepl('Inula', Species_list_format), ]
Genus_sp[, Forest_or_grassland_estimate := 'Forest']
Species_table_2[(
  grepl(" sp$", Species_list_format) |
    grepl("Rhinanthus$", Species_list_format) |
    grepl("Rhinanthus aggr.", Species_list_format)
) &
  !grepl('Inula', Species_list_format), Forest_or_grassland_estimate := 'Grassland']
Species_table_2 = rbind(Species_table_2,
                        Genus_sp)

# Then I fill the new column with a character string, containing all species from the same genus with their relative cover
# and a list of species IDs
for (s in sort(unique(Species_table_2[(
  grepl(" sp$", Species_list_format) |
  grepl("Rhinanthus$", Species_list_format) |
  grepl("Rhinanthus aggr.", Species_list_format)
) & !grepl('Inula', Species_list_format), Species_list_format]))) {
  print(s)
  sp = gsub(' sp', '', s)
  sp = gsub(' aggr.', '', sp)
  
  # In grasslands: fetch all rows w the same species, calculate relative abundance
  df_grass = Abundances_grasslands[grepl(sp, Species) &
                                     !grepl("sp$", Species) &
                                     !grepl("aggr", Species) & !grepl("sp_cf$", Species), ]
  DF_grass = df_grass[, list(cover = sum(value, na.rm = T)), by = Species][, tot := round(100 *
                                                                                            cover / sum(cover, na.rm = T))]
  if (DF_grass[, sum(cover)] == 0) {
    DF_grass$tot = 100
  }
  DF_grass[, Species_format := sapply(Species, function(x) {
    unique(Species_table_2[Species_raw %in% c(x, gsub(' ', '_', x), gsub('_', ' ', x)), Species_list_format])
  })]
  grass_sp = DF_grass[, Sp_prop := paste(Species_format, '=', tot, '%', sep = ''), by = .I][, paste(Sp_prop, collapse = '___')]
  
  grass_ids = Species_table_2[Species_list_format %in% DF_grass$Species_format, paste(unique(IDs), collapse = '___')]
  
  # Same in forests
  df_for = Abundances_forests[grepl(sp, Species) &
                                !grepl("sp$", Species) &
                                !grepl("aggr", Species) & !grepl("sp_cf$", Species), ]
  DF_for = df_for[, list(cover = sum(Cover, na.rm = T)), by = Species][, tot := round(100 *
                                                                                        cover / sum(cover, na.rm = T))]
  if (DF_for[, sum(cover)] == 0) {
    DF_for$tot = 100
  }
  DF_for[, Species_format := sapply(Species, function(x) {
    unique(Species_table_2[Species_raw %in% c(x, gsub(' ', '_', x), gsub('_', ' ', x)), Species_list_format])
  })]
  for_sp = DF_for[, Sp_prop := paste(Species_format, '=', tot, '%', sep = ''), by = .I][, paste(Sp_prop, collapse = '___')]
  for_ids = Species_table_2[Species_list_format %in% DF_for$Species_format, paste(unique(IDs), collapse = '___')]
  
  # If no species found in one of the LU type, replace by the other
  if (grass_sp == '' & for_sp != '') {
    grass_sp = for_sp
    grass_ids = for_ids
  }
  if (grass_sp != '' & for_sp == '') {
    for_sp = grass_sp
    for_ids = grass_ids
  }
  
  # If no species is found in any of the abundance datasets, we take the species found in the other dataset
  if (grass_sp == '' & for_sp == '') {
    all_sp = Species_table_2[grepl(sp, Species_raw), All_names]
    species = unique(unlist(strsplit(paste(
      all_sp, collapse = '___'
    ), '___')))
    species = unlist(sapply(species, function(x) {
      unique(Species_table_2[Species_list_format == x, Species_list_format])
    }))
    species = unique(species[species != ''])
    species = species[!grepl("sp$", species) &
                        !grepl("aggr", species)]
    species = paste(species,
                    '=',
                    round(100 / length(species)),
                    '%',
                    sep = '',
                    collapse = '___')
    
    all_ids = Species_table_2[grepl(sp, Species_raw), IDs]
    ids = unique(unlist(strsplit(
      paste(all_ids, collapse = '___'), '___'
    )))
    ids = ids[ids != "logical(0)" & ids != "character(0)"]
    ids = paste(ids, collapse = '___')
    
    grass_ids = ids
    for_ids = ids
    grass_sp = species
    for_sp = species
  }
  
  Species_table_2[Forest_or_grassland_estimate == 'Forest' &
                    Species_list_format == s, c('All_names', 'IDs') := list(for_sp, for_ids)]
  Species_table_2[Forest_or_grassland_estimate == 'Grassland' &
                    Species_list_format == s,  c('All_names', 'IDs') := list(grass_sp, grass_ids)]
}
beep('mario')

Species_table_2 = Species_table_2[!is.na(Synonyms) & Synonyms != '', ]
Species_table_2 = Species_table_2[All_names != '=Inf%', ]

#fwrite(Species_table_2, file = 'Species_table_complete.csv')



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### ******* PART 2: Cleaning TRY raw data ******* ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# This code cleans the raw output from TRY into a trait x species matrix usable in the exploratories. #
# It is based on a script provided by J. Kattge & S. Tautenhahn, adapted and extended by M. Neyret

# It goes through the following steps:
# 1. Check and format species names
# 2. Identify (and exclude) duplicates
# 3. Identify (and exclude) non-mature plants
# 4. Identify (and exclude) non-healthy plants
# 5. Identify (and exclude) non-standard explosition (i.e. shade or experiments)
# 6. Filter non-standard trait measurements
# 7. Identify (and exclude) measurements with non-standard units
# 8. Identify (and exclude) measurements from different climate zones
# 9. Average by contribution and exclude outliers
# 10. Get average value for each trait/species
# 11. Export data


# Note: For all steps excluding data, we make a list of all observations to exclude
# They will actually be excluded only at steps 10-11

# *************************************** #
#### 1. Check and format species names ####
# *************************************** #
# For now, all species should already be matched with exploratories species, so this part is not needed.
# For future updates it is also possible to rematch species to strict genus + species denomination by
# uncommenting the following lines.

# Species names should be only two words (genus + species)
# TRYdata[, temp.word_length := str_count(AccSpeciesName, "\\w+")] # So we split by blank space, etc. to check

# TRYdata[temp.word_length > 2, unique(AccSpeciesName)] # Check for the different cases
# TRYdata[, new.AccSpeciesName := AccSpeciesName] # Create new columns
# TRYdata[temp.word_length > 2, new.AccSpeciesName := word(AccSpeciesName, 1, 2, sep = " ")] # Keep only 2 first words of the name

# TRYdata <- TRYdata[, .SD, .SDcols = colnames(TRYdata)[!grepl("temp.", colnames(TRYdata))]] # Removes temporary columns
TRYdata <-
  fread(
    "Try_request_11741.txt",
    header = T,
    sep = "\t",
    dec = ".",
    quote = "",
    data.table = T
  )

# ****************************************** #
#### 2. Identify (and exclude) duplicates ####
# ****************************************** #
#### OUTDATED First we check that no ObsDataID is duplicated
##length(TRYdata$ObsDataID) == length(unique(TRYdata$ObsDataID)) # -> OK
## Then there is apparently entries that were entered multiple times
## obs_exclude_duplicate <- unique(TRYdata[DataName == "Duplicate data", ObservationID])
## So we remove all the entries with corresponding ObservationID
## TRYdata = TRYdata[!(ObservationID %in% obs_exclude_duplicate),]

#### NEW DUPLICATE MANAGEMENT METHOD
# The column OrigObsDataID indicates which measurements are duplicated; however some appear only once.

# we list all the OrigObsDataID (supposed to give reference to original data)
origObsDataID = TRYdata[, unique(OrigObsDataID)]
# we keep only those that are actually present in the rest of the data
origObsDataID = origObsDataID[origObsDataID %in% TRYdata$ObsDataID]
# We exclude this data
TRYdata = TRYdata[!(OrigObsDataID %in% origObsDataID), ]

### Old method ###
#all_duplicates = TRYdata[!is.na(StdValue) & TraitName != '', list(average = mean(StdValue), n = length(StdValue)), by = c( 'AccSpeciesName', 'TraitName','Dataset')]
#all_duplicates[, duplicates := duplicated(.SD), .SDcols = c('AccSpeciesName', 'TraitName', 'average', 'n')]
#all_duplicates[, duplicates := duplicated(.SD, fromLast = FALSE), .SDcols = c('AccSpeciesName', 'TraitName', 'average', 'n')]

#all_duplicates = all_duplicates[duplicates == TRUE,]
#all_duplicates = all_duplicates[n >1,]

#for (i in 1:nrow(all_duplicates)){
#  duplicate_info = (all_duplicates[i,])
#  TRYdata[AccSpeciesName == duplicate_info$AccSpeciesName &
#            TraitName == duplicate_info$TraitName &
#            Dataset == duplicate_info$Dataset, duplicate := TRUE]
#}

#TRYdata = TRYdata[is.na(duplicate),]


# ************************************ #
#### 3. Identify  non-mature plants ####
# ************************************ #
TRYdata[DataID %in% c(413), unique(OrigValueStr)] # Check different levels and complete if needed below ("dummy")

TRYdata[, temp.juvenile := OrigValueStr == "dummy" |
          OrigValueStr == "immature" |
          OrigValueStr == "junvenil" |
          OrigValueStr == "juvenil (3 years)" |
          OrigValueStr == "juvenile" |
          OrigValueStr == "Juvenile" |
          OrigValueStr == "juvenile, 11-14 weeks" |
          OrigValueStr == "juvenile, 6 weeks" |
          OrigValueStr == "juveniles" |
          OrigValueStr == "sapling" |
          OrigValueStr == "Sapling (1 - 10 y)" |
          OrigValueStr == "saplings" |
          OrigValueStr == "seedling" |
          OrigValueStr == "Seedling" |
          OrigValueStr == "Seedling (0 - 1 y)" |
          OrigValueStr == "seedlings" |
          OrigValueStr == "seedlings, < 1/2 year" |
          OriglName == "Developmental stage" & OrigValueStr == "S" |
          OriglName == "Seedlings (True/False)" &
          OrigValueStr == "T" |
          OriglName == "Juvenile" & OrigValueStr == "Y"]


obs_exclude_juvenile <-
  unique(TRYdata[temp.juvenile == TRUE, ObservationID])

TRYdata <-
  TRYdata[, .SD, .SDcols = colnames(TRYdata)[!grepl("temp.", colnames(TRYdata))]] # Removes temporary columns


# *********************************************** #
#### 4. Identify non-healthy plants and organs ####
# *********************************************** #
TRYdata[DataID %in% c(3057, 1961), unique(OrigValueStr)] # Check different levels and complete below in 'dummy' if needed
TRYdata[DataID %in% c(3057, 1961), temp.unhealthy := OrigValueStr %in% c("Dead", "dummy")]

obs_exclude_dead <-
  unique(TRYdata[temp.unhealthy == TRUE, ObservationID])

TRYdata <-
  TRYdata[, .SD, .SDcols = colnames(TRYdata)[!grepl("temp.", colnames(TRYdata))]] # Removes temporary columns


# ***************************************** #
#### 5. Identify non-standard exposition ####
# ***************************************** #
TRYdata[, temp.exposition := OrigValueStr %in% c(
  "growth-chamber",
  "growth-chamber",
  "Climate chamber",
  "Climate Chamber",
  "Climate chamber, non-limiting conditions, (cf. dataset reference)",
  "climate chambers",
  "controlled environment room",
  "experimental",
  "experimental treatment",
  "GH",
  "Glass house",
  "Glasshouse",
  "Glasshouse experiment",
  "Greehouse",
  "Green house",
  "greenhouse",
  "Greenhouse",
  "Greenhouse plants",
  "Controlled climate chamber",
  "climate chamber",
  "open-top chamber",
  "open-sided growth chamber",
  "Greenhouse, grrowth container",
  "Greenhouse, Indiana University",
  "Greenhouse: highlight_highpH_competition",
  "Greenhouse: highlight_highpH_nocompetition",
  "Greenhouse: highlight_lowpH_competition",
  "Greenhouse: highlight_lowpH_nocompetition",
  "Greenhouse: lowleight_lowpH_competition",
  "Greenhouse: lowlight_highpH_competition",
  "Greenhouse: lowlight_highpH_nocompetition",
  "Greenhouse: lowlight_lowpH_nocompetition",
  "groth chamber",
  "growth chamber",
  "Growth chamber",
  "Growth Chamber",
  "Growth chamber, -N",
  "Growth chamber, +N",
  "growth chambers",
  "growth_chamber",
  "mesocosm",
  "Shade - Natural environment"
) |
  OriglName == "Artificial growth conditions (G=greenhouse, C=growth chamber)" &
  OrigValueStr == "G" |
  OriglName == "growingCondition" & OrigValueStr == "GH" |
  OriglName == "Natural / Greenhouse" & OrigValueStr == "G"]

TRYdata[DataID %in% c(327, 308, 210) &
          !temp.exposition, unique(OrigValueStr)] # Check different levels and complete above in 'dummy' if needed

obs_exclude_exposition <-
  unique(TRYdata[temp.exposition == TRUE, ObservationID])

TRYdata <-
  TRYdata[, .SD, .SDcols = colnames(TRYdata)[!grepl("temp.", colnames(TRYdata))]] # Removes temporary columns

# *********************************************** #
#### 6. Filter non-standard trait measurements ####
# *********************************************** #
# Check the different modalities of DataName for each trait
tapply(TRYdata[TraitName != "" & DataID != "",]$DataName,
       TRYdata[TraitName != "" & DataID != "",]$TraitName,
       unique)

# We check trait by trait and remove the Data which are not standard
remove_DataName <- c()

# Trait: `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included`
# DataName: "SLA: petiole  included"    "SLA: petiole included (2)" "SLA: petiole included (1)"
## -> all should be ok

# Trait:  `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded`
# DataName: "SLA: undefined if petiole in- or excluded"  "SLA: undefined if petiole in- or excluded (1)" "Leaf mass per area (LMA)"
## -> all should be ok

# Trait: `Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)`
# DataName: "Leaf dry matter content per leaf water-saturated mass (LDMC)"
#           "Leaf water content per leaf water-saturated mass (LWC, 1-LDMC)"  ** !!!! ???? Remove **
#           "LDMC min"       ** Remove **
#           "LDMC max"       ** Remove **
#           "Leaf dry fraction"
#           "Leaf dry matter concentration  (NIRS based)"
#           "Leaf dry matter content (min)"
#           "Leaf dry matter content (max)"
#           "Leaf lamina dry matter content (LDMC)"

remove_DataName <-
  c(
    remove_DataName,
    "Leaf water content per leaf water-saturated mass (LWC, 1-LDMC)",
    "LDMC min",
    "LDMC max"
  )

# Trait: $`Leaf nitrogen (N) content per leaf dry mass`
# DataName: "Leaf nitrogen content per dry mass (Nmass)"
#           "N amount at 15 N measurement"
#           "Leaf senescent nitrogen (N) content per dry mass" ** Remove **
#           "Leaf nitrogen concentration  (NIRS based)"

remove_DataName <-
  c(remove_DataName,
    "Leaf senescent nitrogen (N) content per dry mass")

# Trait: `Leaf phosphorus (P) content per leaf dry mass`
# DataName: "Leaf phosphorus content per dry mass (Pmass)"
#           "Leaf senescent phosphorus (P) content per dry mass" ** Remove **

remove_DataName <-
  c(remove_DataName,
    "Leaf senescent phosphorus (P) content per dry mass")

# Trait: `Mycorrhiza type`
# DataName: "Mycorrhizal type"
#           "Original term for mycorrhizal type given by Selivanov"
#           "Mycorrhiza: Arbuscular mycorrhizal (AM) fungi"
#           "Mycorrhiza: Ectomycorrhizal (ECM) fungi"
#           "Mycorrhiza: No Mycorrhizal (NM) fungi"
#           "Mycorrhiza: Ericoid Mycorrhiza (ER) fungi"
#           "Mycorrhiza: Orchid mycorrhizal (ORM) fungi"
#           "Mycorrhiza: Arbutoid mycorrhizal (ARB) fungi"
#           "Mycorrhiza: Presence of any non-AM mycorrhizal fungus (binary)"
#           "Mycorrhiza: Presence of any  alternative to AM mycorrhizal fungus (binary)"
#           "Mycorrhiza: Likelihood the species was in the 'Stable AM Loss' state inferred under our best HRM-model (SI Figure 1, Extended Methods)" ** Remove **
#           "Mycorrhiza: Likelihood the species was in the 'Labile' state inferred under our best HRM-model (SI Figure 1, Extended Methods)"         ** Remove **
#           "Mycorrhiza: Likelihood the species was in the 'Stable AM' state inferred under our best HRM-model (SI Figure 1, Extended Methods)"      ** Remove **
#           "Mycorrhiza: Likelihood the species had lost AM interactions under our best HRM-model"                                                   ** Remove **
#           "Mycorrhiza: Likelihood the species had retained AM interactions under our best HRM-model"                                               ** Remove **
#           "Mycorrhiza: Inferred binary AM state under our best HRM-model (i.e. \"Yes\" if AM_retained_likelihood > 0.5)"                           ** Remove **
#           "Mycorrhiza type according to Maherali"

remove_DataName <- c(
  remove_DataName,
  "Mycorrhiza: Likelihood the species was in the 'Stable AM Loss' state inferred under our best HRM-model (SI Figure 1, Extended Methods)",
  "Mycorrhiza: Likelihood the species was in the 'Labile' state inferred under our best HRM-model (SI Figure 1, Extended Methods)",
  "Mycorrhiza: Likelihood the species was in the 'Stable AM' state inferred under our best HRM-model (SI Figure 1, Extended Methods)",
  "Mycorrhiza: Likelihood the species had lost AM interactions under our best HRM-model",
  "Mycorrhiza: Likelihood the species had retained AM interactions under our best HRM-model",
  "Mycorrhiza: Inferred binary AM state under our best HRM-model (i.e. \"Yes\" if AM_retained_likelihood > 0.5)"
)

# Trait: Plant height vegetative`
#           "Plant height vegetative"                                           "Vegetative plant height, not elongated"
#           "Flowering plant height, heighest leaf elongated" ** Remove **       "Flowering plant height, heighest leaf not elongated"
#           "Average tree height"                                                "Height at 20 Years"
#           "Maximum plant height" ** Remove **                                  "Plant height observed"
#           "Plant height (unspecified if vegetative or reproductive)"           "Plant height of aquatic plants"
#           "Height of seedlings"                                                "Maximum height"
#           "Maximum height min"  ** Remove **                                   "Maximum height max"   ** Remove **
#           "Maximum height extreme"    ** Remove **                              "Plant height from flora"
#           "Plant height vegetative min."                             "          Plant height vegetative max."
#           "Vegetative plant height elongated"   ** Remove **

remove_DataName <- c(
  remove_DataName,
  "Flowering plant height, heighest leaf elongated",
  "Maximum plant height",
  "Maximum height min",
  "Maximum height extreme",
  "Maximum height max",
  "Vegetative plant height elongated"
)

# Trait: Seed dry mass`
#           "Seed dry mass"                               "Seed mass"
#           "Seed mass original value: mean"              "Seed mass min"  ** Remove **
#           "Seed mass max"    ** Remove **               "Seed mass original value: min" ** Remove **
#           "Seed mass original value: max"  ** Remove **

remove_DataName <- c(
  remove_DataName,
  "Seed mass min",
  "Seed mass max",
  "Seed mass original value: min",
  "Seed mass original value: max"
)
TRYdata[TraitName == "Seed dry mass", table(DataName)]


# Trait: Stem specific density (SSD) or wood density (stem dry mass per stem fresh volume)`
#           "Wood density; stem specific density; wood specific gravity"
# --> OK

# Trait: Mycorrhizal infection intensity
#           "Intensity of mycorrhizal infection"
#           "Root length colonisation by vesicular-arbucular mycorrhizal fungi"

# --> OK

TRYdata[TraitName == 'Mycorrhizal infection intensity', table(DataName)]
remove_DataName <- c(remove_DataName,
                     "Root length colonisation by vesicular-arbucular mycorrhizal fungi")


ObsDataID_exclude_nonstandard <-
  unique(TRYdata[DataName %in% remove_DataName, ObsDataID])



### I'm also renaming mycorrhiza into homogeneous names
TRYdata[TraitName == "Mycorrhiza type" &
          !(ObsDataID %in% ObsDataID_exclude_nonstandard),
        StdValue_myco := ifelse(
          OrigValueStr %in% c(
            "ECTO",
            "ecto",
            "E.ch.ect.",
            "Ecto",
            "EC",
            "ectomycorrhizal",
            "Mycorrhiza: Ectomycorrhizal (ECM) fungi",
            "EM",
            "ec?",
            "Ectomycorrhiza"
          ),
          "-Ectomycorrhiza",
          ifelse(
            OrigValueStr %in% c(
              "Arbutoid",
              "AbtM",
              "Ecto arbut.",
              "E.t.ect.arb.",
              "Mycorrhiza: Arbutoid mycorrhizal (ARB) fungi"
            ),
            "-Arbutoid mycorrhiza",
            ifelse(
              OrigValueStr %in% c(
                "VAM",
                "AMNM",
                "NM/AM",
                "?va",
                "VA",
                "Mycorrhiza: Arbuscular mycorrhizal (AM) fungi",
                "AM",
                "arbuscular",
                "Gram-AM",
                "Forb-AM",
                "Ph.th.end."
              ),
              "-Arbuscular mycorrhiza",
              # ifelse(OrigValueStr %in% c("AMNM", "NM/AM"), "-Weak arbuscular mycorrhiza",
              ifelse(
                OrigValueStr %in% c("AM + EM", "EC/AM"),
                "-Ecto and/or arbuscular mycorrhiza",
                ifelse(
                  OrigValueStr %in% c(
                    "ericoid",
                    "Mycorrhiza: Ericoid Mycorrhiza (ER) fungi",
                    "Ericoid",
                    "E.t.ect.er.",
                    "ErM",
                    "EM",
                    "Ecto er.",
                    "Endo er.",
                    "ER",
                    "ERICOID",
                    "pyroloid"
                  ),
                  "-Ericoid mycorrhiza",
                  ifelse(
                    OrigValueStr %in% c(
                      "Orchid",
                      "OrM",
                      "orchid",
                      "E.t.end.",
                      "Mycorrhiza: Orchid mycorrhizal (ORM) fungi"
                    ),
                    "-Orchid mycorrhiza",
                    ifelse(
                      OrigValueStr %in% c(
                        "DS",
                        "Ps.end.",
                        "Ectendo",
                        "Endo",
                        "ecto, ectendo",
                        "Endo-unidentified",
                        "Ecto, endo, VA"
                      ),
                      "-Other or undefined endomycorrhiza",
                      ifelse(
                        OrigValueStr %in% c(
                          "NM",
                          0,
                          "Non",
                          "no",
                          "absent",
                          "Mycorrhiza: No Mycorrhizal (NM) fungi",
                          "N"
                        ),
                        "-Non mycorrhizal",
                        NA
                      )
                    )
                  )
                )
              )
              # )
            )
          )
        )]



# ******************************************************** #
#### 7.  Identify measurements with non-standard units ####
# ******************************************************** #
# tapply(
#  TRYdata[TraitName != "" & DataID != "", ]$UnitName,
#  TRYdata[TraitName != "" & DataID != "", ]$TraitName,
#  unique
# )

# -> We can use the StdValue column, which already standardize values, for all traits (!!! to check if new traits are added)
# This solves this issue


# *********************************************************** #
#### 8. Identify measurements from different climate zones ####
# *********************************************************** #
## We currently decided to keep all datasets in the analysis, whatever their geographic origin is.
## You can change this by uncommenting the obs_exclude_geography line in the following session

histogram(TRYdata[TRYdata$DataID == 59, StdValue, by = Dataset]$StdValue, breaks = 100) # Latitude: bw 45 et 55
histogram(TRYdata[TRYdata$DataID == 60, StdValue, by = Dataset]$StdValue, breaks = 100) # Longitude: -5-30 (Europe) or 5-25 (Central Europe)
histogram(TRYdata[TRYdata$DataID == 61, StdValue, by = Dataset]$StdValue, breaks = 100) # Altitude: < 1500

data_exclude_lat = TRYdata[TRYdata$DataID == 59 &
                             (StdValue < 45 | StdValue > 55), ]$ObservationID
data_exclude_lon = TRYdata[TRYdata$DataID == 60 &
                             (StdValue < 5 | StdValue > 25), ]$ObservationID
data_exclude_alt = TRYdata[TRYdata$DataID == 61 &
                             (StdValue > 1500) , ]$ObservationID

obs_exclude_geography = unique(c(data_exclude_lat, data_exclude_lon, data_exclude_alt))

# ******************************************************* #
####  9. Average by contributor, then exclude outliers ####
# ******************************************************* #

# We first need to remove all previous identified sources of bias to effectively identify "outliers".
#TRYdata_all_data = copy(TRYdata)
#TRYdata = copy(TRYdata_all_data)
TRYdata <- TRYdata[!(
  ObservationID %in%
    c(#  obs_exclude_geography, # %***% Uncomment to restrict geographc origin
      # obs_exclude_exposition,  # %***% Comment to take data from all experimental conditions
      obs_exclude_dead,
      obs_exclude_juvenile) |
    ObsDataID %in% ObsDataID_exclude_nonstandard
),]

### We then average data for each contributor, to avoid too high weights from one dataset to bias the data

## 1. For quantitative traits -> take the mean of all estimates, removing min/max values when calculating the mean
# For most traits, the dominant value type is trait mean or single measurements, and we have enough data to keep only those and remove
# StdValues that correspond to minimum/maxmum measurements in the initial datasets.
# However, this is not the case for Myco infection intensity, for which maximum values are provided in most cases, so we keep these values
TRYdata[TraitName == 'Mycorrhizal infection intensity', ValueKindName := tolower(ValueKindName)]

TRYdata_contributors <-
  TRYdata[TraitName != "Mycorrhiza type" &
            !(ValueKindName %in% c('Minimum', 'Maximum', 'Low', 'High')),
          list(StdValue = mean(StdValue, na.rm = T)), by = c(
            "LastName",
            "FirstName",
            "DatasetID",
            "Dataset",
            "SpeciesName",
            "AccSpeciesID",
            "AccSpeciesName",
            "TraitID",
            "TraitName",
            "UnitName",
            "Reference"
          )]

# To calculate the range we do not exclude min/max
TRYdata_contributors_min_max <-
  TRYdata[TraitName != "Mycorrhiza type",
          list(Min = min(StdValue, na.rm = T),
               Max = max(StdValue, na.rm = T)),
          by = c(
            "LastName",
            "FirstName",
            "DatasetID",
            "Dataset",
            "SpeciesName",
            "AccSpeciesID",
            "AccSpeciesName",
            "TraitID",
            "TraitName",
            "UnitName",
            "Reference"
          )]
TRYdata_contributors_min_max[grepl('Inf', as.character(Min)), Min := NA]
TRYdata_contributors_min_max[grepl('Inf', as.character(Max)), Max := NA]

TRYdata_contributors_min_max[, Range := ifelse(Min != Max, paste(round(as.numeric(Min), 2), round(as.numeric(Max), 2), sep = ' '), round(as.numeric(Min), 2))][, Range := gsub("Inf -Inf", NA, Range)]
TRYdata_contributors = merge(TRYdata_contributors, TRYdata_contributors_min_max)

## 2. For mycorrhizal type traits
# Need to correct one binary database first
TRYdata[Dataset == "Mycorrhizal Association Database" &
          TraitName == "Mycorrhiza type", OrigValueStr := ifelse(OrigValueStr == "Yes", DataName, NA)]

# Function to choose the most common type, if it's twice as common as the 2nd most common, otherwise describe as 'Mixed'
choose_myco <- function(data) {
  data = data[!is.na(data) & data != 'NA']
  x <- table(data)
  if (length(x) == 0) {
    value = NA
    prop = NA
  }
  else {
    # if ("-Weak arbuscular mycorrhiza" %in% names(x) & "-Arbuscular mycorrhiza" %in% names(x)) {
    #    x["-Arbuscular mycorrhiza"] <- x["-Arbuscular mycorrhiza"] + x["-Weak arbuscular mycorrhiza"]
    #    x <- x[names(x) != "-Weak arbuscular mycorrhiza"]
    #  }
    if (length(x) == 1) {
      value = names(x)
      prop = '100 %'
    }
    else {
      X <- sort(x, decreasing = TRUE)
      if (X[1] >= 2 * X[2]) {
        value = names(X[1])
        prop = paste(round(X[1] / sum(X), 2) * 100, '%')
      }
      else {
        value = 'Mixed'
        prop = 'NA'
      }
    }
  }
  return(list(value = value, proportion = prop))
}

TRYdata_contributors_myco = TRYdata[TraitName == "Mycorrhiza type", list(
  StdValue = as.character(choose_myco(StdValue_myco)$value),
  Range = paste(unique(StdValue_myco), collapse = ' '),
  Min = NA,
  Max = NA
), by = c(
  "LastName",
  "FirstName",
  "DatasetID",
  "Dataset",
  "SpeciesName",
  "AccSpeciesID",
  "AccSpeciesName",
  "TraitID",
  "TraitName",
  "UnitName",
  "Reference"
)]

TRYdata_contributors = rbind(TRYdata_contributors, TRYdata_contributors_myco)

## We also want to have an idea of how many measurements were made.
TRYdata[, Replicates2 := Replicates]
# If Replicates = '' but StdValue or StdValue_myco non NA, we consider that there is 1 replicate
TRYdata[!(is.na(StdValue) &
            is.na(StdValue_myco)) & Replicates == '', Replicates2 := 1]
# If ValueKindName is single value, we consider that it's one replicate.
TRYdata[ValueKindName == 'Single', Replicates2 := 1]
# If Replicate == Originat str value, it probably means that it's an error -> replace by 1.
TRYdata[Replicates == OrigValueStr, Replicates2 := 1]
# If StdValue = NA, we can't count it -> replace by 0
TRYdata[is.na(StdValue) & is.na(StdValue_myco), Replicates2 := 0]
# There is also a category "Number of replicates" in DataName column, which we need to re-merge into the main dataset to use it
Replicates_per_Observation = TRYdata[DataName == 'Number of replicates', list(Rep = unique(OrigValueStr)), by = c('ObservationID')]

TRYdata[TraitName != '' , Replicates2 := sapply(ObservationID, function(x) {
  if (length(Replicates_per_Observation[ObservationID == x, ]$Rep) > 0) {
    Replicates_per_Observation[ObservationID == x, ]$Rep
  }
  else
    (1)
})]
TRYdata[, Replicates2 := as.numeric(Replicates2)]

TRYdata_contributors_n_measurements <-
  TRYdata[!(ValueKindName %in% c('Minimum', 'Maximum', 'Low', 'High')),
          list(n_replicates = sum(Replicates2)),
          by = c(
            "LastName",
            "FirstName",
            "DatasetID",
            "Dataset",
            "SpeciesName",
            "AccSpeciesID",
            "AccSpeciesName",
            "TraitID",
            "TraitName",
            "UnitName",
            "Reference"
          )]


### Merge into 1 dataset
TRYdata_contributors_all <- merge(TRYdata_contributors,
                                  TRYdata_contributors_n_measurements)

# identify outliers
is.outlier <- function(x) {
  print(x)
  outliers <- boxplot.stats(x)$out
  is.outlier <- x %in% outliers
  return(is.outlier)
}

TRYdata_contributors_all[TraitName != "Mycorrhiza type", is.outlier := is.outlier(as.numeric(StdValue)), by = c("TraitName", "AccSpeciesName")]
TRYdata_contributors_all[TraitName == "Mycorrhiza type", is.outlier := FALSE, by = c("TraitName", "AccSpeciesName")]

TRYdata_contributors_all = TRYdata_contributors_all[TraitName != '', ]
TRYdata_contributors_all = TRYdata_contributors_all[!(TraitName != "Mycorrhiza type" &
                                                        UnitName == ''), ]

# *********************************************************** #
#### 10. Finally, get average value for each trait/species ####
# *********************************************************** #

### Number of datasets per value
numberDatasets = TRYdata_contributors_all[is.outlier == FALSE,
                                          list(
                                            numberDatasets = length(unique(Dataset)),
                                            References = paste(unique(Reference), collapse = '___')
                                          ),
                                          by = c("TraitName", "AccSpeciesName", 'AccSpeciesID')]

TRYdata_contributors_all[AccSpeciesName %in% All_names &
                           TraitName == trait, ]

###  For quantitative traits ###
Final_data <-
  TRYdata_contributors_all[TraitName != "Mycorrhiza type" &
                             !is.outlier,
                           list(
                             StdValue = as.character(mean(as.numeric(StdValue), na.rm = T)),
                             traitSd =  as.character(sd(as.numeric(StdValue), na.rm = T)),
                             numberReplicates = sum(n_replicates, na.rm = T),
                             Min = min(Min, na.rm = T),
                             Max = max(Max, na.rm = T)
                           ),
                           by = c("TraitName", "AccSpeciesName", "UnitName", 'AccSpeciesID')]
Final_data[traitSd == 0, traitSd := NA]
Final_data[, Range := ifelse(Min != Max, paste(round(as.numeric(Min), 2), round(as.numeric(Max), 2), sep = ' '), round(as.numeric(Min), 2))][, Range := gsub("Inf -Inf", NA, Range)]

### We, again, need to aggregate for mycorrhiza
Final_data_myco = TRYdata_contributors_all[TraitName == "Mycorrhiza type",
                                           list(
                                             StdValue = as.character(choose_myco(StdValue)$value),
                                             traitSd = NA,
                                             numberReplicates = sum(n_replicates),
                                             Min = NA,
                                             Max = NA,
                                             Range = sapply(Range, function(x) {
                                               paste(unique(unlist(strsplit(x, ' '))), collapse = " ")
                                             })
                                           ),
                                           by = c("TraitName", "AccSpeciesName", "UnitName",  'AccSpeciesID')]

Final_data <- rbind(Final_data,
                    Final_data_myco)

Final_data = merge(Final_data, numberDatasets)

Final_data[Final_data$traitSd == 0, traitSd := NA]

Final_data[TraitName == 'Mycorrhiza type' &
             StdValue == "Unknown", StdValue := NA]
Final_data = Final_data[TraitName != '', ]

Final_data$Range = as.character(Final_data$Range)

Final_data2 = Final_data[, list(
  traitName = TraitName,
  scientificName = AccSpeciesName,
  taxonID = AccSpeciesID,
  traitValue = StdValue,
  traitUnit = UnitName,
  traitValueStd = StdValue,
  traitSd = traitSd,
  numberReplicates = numberReplicates,
  Min = Min,
  Max = Max,
  Range = Range,
  numberDatasets = numberDatasets ,
  References = References
)]


# ********************************************************** #
####       11. Export dataset, raw and by species         ####
# ********************************************************** #
#fwrite(Final_data2, "TRY_plant_traits.csv")
#fwrite(Final_data2, "TRY_plant_traits_exposition.csv")


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### ******* PART 3: Trait to CWM values ******* ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# This script takes as input the abundances and species-level traits of plant species
# and outputs a CWM matrix averaged for all considered years.

# This is the matching table between Exploratories and TRY data
Species_table <- data.table(read.csv("Species_table_complete.csv"))
Species_table = Species_table[, c(
  'Species_raw',
  "Species_list_format",
  "Species_list_corrected",
  "Synonyms",
  "All_names",
  'Forest_or_grassland_estimate'
)]

# This is the AG trait data
AG_traits <- fread("TRY_plant_traits.csv", header = TRUE)
AG_traits = AG_traits[traitName != '', c(
  'traitName',
  'scientificName',
  'traitUnit',
  'traitValue',
  'traitSd',
  'traitValueStd',
  "numberReplicates",
  "Min",
  "Max",
  "Range",
  "numberDatasets",
  "References"
)]
AG_traits[, traitName := gsub(' ', '_', traitName)]
AG_traits[, traitName := gsub('[,/()]', '', traitName)]
# 0 replicate does not make sense when a StdValue is provided. A 0 is present when we don't know how many replicates were
# included initially. We replace it by 1.
AG_traits[!is.na(traitValue) &
            numberReplicates == 0, numberReplicates := 1]
AG_traits[traitSd == '', traitSd := NA]

# This is the BG trait data. We don't use seed mass as we have data from TRY
colnames(AG_traits)[!(colnames(AG_traits) %in% colnames(Root_traits))]
Root_traits[, scientificName := gsub('_', ' ', scientificName)]
Root_traits[, c(
  'traitSd',
  "numberReplicates",
  "Min",
  "Max",
  "Range",
  "numberDatasets",
  "References",
  'traitValue',
  'traitValueStd'
) :=
  list(
    NA,
    NA,
    NA,
    NA,
    NA,
    1,
    "Bexis dataset 26587",
    as.character(traitValue),
    as.character(traitValueStd)
  )]

All_traits = rbind(AG_traits[, c(
  'traitName',
  'scientificName',
  'traitUnit',
  'traitValue',
  'traitSd',
  'traitValueStd',
  "numberReplicates",
  "Min",
  "Max",
  "Range",
  "numberDatasets",
  "References"
)],
Root_traits[traitName != "Seed_mass", c(
  'traitName',
  'scientificName',
  'traitUnit',
  'traitValue',
  'traitSd',
  'traitValueStd',
  "numberReplicates",
  "Min",
  "Max",
  "Range",
  "numberDatasets",
  "References"
)])

search_species_table <-
  function(species_raw,
           trait_name,
           all_names,
           Trait_Table) {
    All_names <- unlist(strsplit(as.character(all_names), "___"))
    # This function takes as input a "raw" species name, as found in the BE, and associates to it the different species (e.g. for aggregate species)
    Trait_data <-
      Trait_Table[scientificName %in% All_names &
                    traitName == as.character(trait_name),]
    trait_value <- NA
    trait_value2 <- NA
    unit_value <- NA
    sd_value_w <- NA
    sd_type <- NA
    trait_min <- NA
    trait_max <- NA
    trait_range <- NA
    numberDatasets = NA
    numberReplicates = NA
    References = NA
    
    if (nrow(Trait_data) > 0)  {
      References = ifelse(
        length(Trait_data$References) == 1,
        Trait_data$References,
        paste(unique(unlist(
          strsplit(Trait_data$References, '___', useBytes = TRUE)
        )), collapse = '___')
      )
      
      numberDatasets = sum(as.numeric(Trait_data$numberDatasets))
      numberReplicates = sum(as.numeric(Trait_data$numberReplicates))
      Trait_data = Trait_data[!is.na(traitValue) &
                                !is.nan(traitValue) & traitValue != 'NaN', ]
      if (!any(is.na(as.numeric(Trait_data$traitValue)))) {
        # If trait quantitative
        #print(Trait_data)
        trait_value <-
          as.character(mean(as.numeric(Trait_data$traitValue), na.rm = T))
        trait_value2 <-
          as.character(weighted.mean(
            as.numeric(Trait_data$traitValue),
            as.numeric(Trait_data$numberDatasets),
            na.rm = T
          ))
        
        if (is.na(trait_value2)) {
          trait_value2 = trait_value
        }
        
        if (!all(is.na(Trait_data$Min))) {
          trait_min = min(as.numeric(Trait_data$Min), na.rm = T)
          trait_max = max(as.numeric(Trait_data$Max), na.rm = T)
        }
        if (all(is.na(Trait_data$Min))) {
          trait_min = min(as.numeric(Trait_data$traitValue), na.rm = T)
          trait_max = max(as.numeric(Trait_data$traitValue), na.rm = T)
        }
        
        trait_range = paste(signif(trait_min, 2), signif(trait_max, 2), collapse = ' ')
        
        if (length(Trait_data$traitValue) > 1) {
          #traits_not_na = Trait_data[!is.na(StdValue) & StdValue != 'NA' ,]
          sd_value_w = sd(Trait_data$traitValue, na.rm = T)
          
          # Here I chose to calculate sd as simply the sd of average values per species.
          # To be more mathematically correct one should calculate it as the weighted sd as sqrt(SOMME[a * var(A)]),
          # with a the relative weight equal to the number of contribution divided by total number
          # I did not because in some cases we don't have any sd value for a given species
          
          #if (nrow(traits_not_na)> 0){
          #sd = as.numeric(traits_not_na$traitSd)
          #  sd = as.numeric(traits_not_na$traitSd)
          
          #  n = as.numeric(traits_not_na$numberDatasets)
          #  means = as.numeric(Trait_data$StdValue)
          #  sd_value_w <- as.character(sqrt(sum(n*(sd^2 + (means-as.numeric(trait_value2))^2),  na.rm = T)/sum(n, na.rm = T)))
          #   sd(as.numeric(traits_not_na$StdValue), na.rm = T)
          #  }
          #  if (nrow(traits_not_na)== 0){
          #
          #  }
          sd_type = "pooled variance"
        }
        if (length(Trait_data$traitValue) == 1) {
          sd_type = "single species variance"
          sd_value_w = Trait_data$traitSd
          trait_min <- Trait_data$Min
          trait_max <- Trait_data$Max
          trait_range <- gsub('___', ' ', Trait_data$Range)
        }
      }
      
      else {
        # trait qualitative
        trait_value <-
          paste(unique(Trait_data$traitValue[!is.na(Trait_data$traitValue) &
                                               Trait_data$traitValue != 'NA']), collapse = "___")
        trait_value2 = trait_value
        sd_value_w = NA
        trait_min <- NA
        trait_max <- NA
        trait_range <- NA
        
      }
      unit_value <- unique(Trait_data$traitUnit)
      unit_value = unit_value[!is.na(unit_value)]
      unit_value = unit_value[!unit_value %in% c('NA', '')]
      unit_value <- paste(unique(unit_value), collapse = "___")
    }
    res = c(
      traitValue = trait_value,
      traitValue2 = trait_value2,
      traitUnit = unit_value,
      traitSd = sd_value_w,
      Min = trait_min,
      Max = trait_max,
      Range = trait_range,
      References = References,
      numberDatasets = numberDatasets,
      numberReplicates = numberReplicates
    )
    # print(res)
    return(as.list(res))
  }

# Note: traitValue is the trait value averaged over all synonyms/species names
# traitValue2 is the average weighted by the respective number of contributors which provided data on the specific name.

# *************************************************************************************************** #
#### 1. Match trait data with species names used in Exploratories based on species matching table  ####
# *************************************************************************************************** #

trait_names <- unique(All_traits$traitName)
Traits_dummy <-
  data.table(matrix(
    ncol = length(trait_names),
    nrow = nrow(Species_table)
  ))
Traits_dummy[, (colnames(Traits_dummy)) := lapply(.SD, as.character), .SDcols = colnames(Traits_dummy)]
names(Traits_dummy) <- as.character(trait_names)
Species_table2 <-
  cbind(Species_table[, .SD, .SD = colnames(Species_table)[!(colnames(Species_table) %in% c('IDs', 'rowID'))]], Traits_dummy)

# Melt the datatable
Species_table_melt <- melt.data.table(
  Species_table2,
  id.vars = c(
    "Species_raw",
    "Species_list_format",
    "Species_list_corrected",
    "Forest_or_grassland_estimate",
    "Synonyms",
    "All_names"
  ),
  value.name = "traitValue",
  variable.name = "traitName"
)

Species_table_melt$traitValue <-
  as.character(Species_table_melt$traitValue)
Species_table_melt[, traitUnit := character()]
# Remove duplicates if there is any
Species_table_melt <- unique(Species_table_melt)


# Here we fill each trait value with the average of the synonyms or corresponding species
Species_table_melt[, c(
  "traitValue",
  'traitValue2',
  "traitUnit",
  'traitSd',
  'Min',
  'Max',
  'Range',
  'References',
  'numberDatasets',
  'numberReplicates'
) := search_species_table(Species_raw, traitName, All_names, All_traits), by = c("Species_raw", "traitName")]
Species_table_melt2 = copy(Species_table_melt)
# For species only identified at the genus level, we have to take into account the relative proportion
# of each species to calculate weighted mean and variance for each.
# %***% Change traitValue2 to traitValue in the following lines if you wish to use non-weighted data
for (s in unique(Species_table_melt2[grepl("%", All_names), All_names])) {
  print(s)
  list_spe <- data.table(Raw = unlist(strsplit(s, "___")))
  list_spe$Species_list_format <- sapply(list_spe$Raw, function(x) {
    unlist(strsplit(x, "="))[1]
  })
  list_spe$Weight <-
    as.numeric(gsub("%", "", sapply(list_spe$Raw, function(x) {
      unlist(strsplit(x, "="))[2]
    })))
  list_spe = list_spe[Species_list_format != 'NA' & Weight > 0, ]
  
  for (trait in trait_names) {
    trait_value = NA
    unit = NA
    trait_value = NA
    trait_sd = NA
    trait_min = NA
    trait_max = NA
    trait_range = NA
    
    print(trait)
    DF <-
      Species_table_melt2[Species_list_format %in% c(list_spe$Species_list_format,
                                                     gsub(' ', '_', list_spe$Species_list_format)) &
                            traitName == trait,]
    DF = unique(DF[, .SD, .SDcols = colnames(DF)[!(colnames(DF) %in% c('Species_raw'))]])
    DF_with_weight <-
      merge(list_spe[, c("Species_list_format", "Weight")], DF)
    
    if (sum(!is.na(as.numeric(DF_with_weight$traitValue))) > 0) {
      # Check if at least one numeric value
      trait_value <-
        weighted.mean(as.numeric(DF_with_weight$traitValue2),
                      DF_with_weight$Weight,
                      na.rm = T)
      
      # Calculation of weighted variance based on https://stats.stackexchange.com/questions/51442/weighted-variance-one-more-time
      
      if (length(unique(DF_with_weight$traitValue2)) == 1) {
        trait_sd = unique(DF_with_weight$traitSd)
      }
      if (length(unique(DF_with_weight$traitValue2)) > 1) {
        v1 = sum(DF_with_weight$Weight)
        v2 =  sum(DF_with_weight$Weight ^ 2)
        trait_var = v1 / (v1 ^ 2 - v2) * sum(DF_with_weight$Weight * (as.numeric(DF_with_weight$traitValue2) - trait_value) ^
                                               2)
        trait_sd = sqrt(trait_var)
      }
      
      if (!all(is.na(DF_with_weight$Min))) {
        trait_min = min(as.numeric(DF_with_weight$Min), na.rm = T)
        trait_max = min(as.numeric(DF_with_weight$Max), na.rm = T)
      }
      if (all(is.na(DF_with_weight$Min))) {
        trait_min = min(as.numeric(DF_with_weight$traitValue2), na.rm = T)
        trait_max = min(as.numeric(DF_with_weight$traitValue2), na.rm = T)
      }
      
      trait_range = paste(signif(trait_min, 2), signif(trait_max, 2), sep = ' ')
      numberDatasets = sum(as.numeric(DF_with_weight$numberDatasets))
      numberReplicates = sum(as.numeric(DF_with_weight$numberReplicates))
      References = ifelse(
        length(DF_with_weight$References) == 1,
        DF_with_weight$References,
        paste(unique(unlist(
          strsplit(
            paste(DF_with_weight$References[!is.na(DF_with_weight$References)], collapse = '___'),
            '___',
            useBytes = TRUE
          )
        )), collapse = '___')
      )
    }
    else {
      values = DF_with_weight$traitValue2
      values = values[!is.na(values) &
                        values != 'NA' & values != '']
      trait_value <- paste(unique(values), collapse = "___")
    }
    unit = unique(DF_with_weight$traitUnit)
    if (length(unit) > 1) {
      unit = unit[!is.na(unit)]
    }
    Species_table_melt2[All_names == s & traitName == trait,
                        c("traitValue",
                          "traitUnit",
                          "traitValue2",
                          "traitSd",
                          "Min",
                          "Max",
                          "Range") :=
                          list(trait_value,
                               unit,
                               trait_value,
                               trait_sd,
                               trait_min,
                               trait_max,
                               trait_range)]
    Species_table_melt2[All_names == s &
                          traitName == trait, ]$References = References
    Species_table_melt2[All_names == s &
                          traitName == trait, ]$numberDatasets = numberDatasets
    Species_table_melt2[All_names == s &
                          traitName == trait, ]$numberReplicates = numberReplicates
  }
}

Species_table_melt2[traitName == 'Mycorrhiza_type', traitValue := ifelse(grepl('Mixed', traitValue), 'Mixed', traitValue)]
Species_table_melt2[traitName == 'Mycorrhiza_type', traitValue := gsub('[___]', '', traitValue)]
Species_table_melt2[traitName == 'Mycorrhiza_type', traitValue := gsub('-', '', traitValue)]
Species_table_melt2[traitName == 'Mycorrhiza_type', traitValue := gsub('NA', '', traitValue)]
Species_table_melt2[traitName == 'Mycorrhiza_type' &
                      (is.na(traitValue) | traitValue == ''), traitValue := NA]
Species_table_melt2[traitName == 'Mycorrhiza_type' &
                      traitName == 'Mycorrhiza_type' %in% c("-Non mycorrhizal-Arbuscular mycorrhiza",
                                                            "-Arbuscular mycorrhiza-Non mycorrhizal"), traitValue := 'Mixed']
Species_table_melt2[traitName == 'Mycorrhiza_type', traitValue2 := traitValue]
Species_table_melt2[traitValue %in% c('', 'NA', 'NaN', "NA___NaN"), traitValue := NA]
Species_table_melt2[traitValue2 %in% c('', 'NA', 'NaN', "NA___NaN"), traitValue := NA]
Species_table_melt2[traitSd == '', traitSd := NA]
Species_table_melt2 = Species_table_melt2[!is.na(traitValue)]
Species_table_melt2[Range == '', Range := NA]

# Last pass to make sure that everything fits with species names in the Abundance datasets:
# we add both spellings
Species_table_melt2[, Species_list_format := gsub(' ', '_', Species_list_format)]

Species_duplicates = copy(Species_table_melt2[Species_list_format %in% c(
  'Pinus_mugo',
  'Euphrasia_sp',
  'Trifolium_montanum',
  'Erigeron_annuus',
  'Medicago_varia',
  'Bromus_hordeaceus_aggr_incl_B_commutatus'
),])
Species_duplicates[, Species_list_format := recode(
  Species_list_format,
  'Pinus_mugo' = 'Pinus_cf_mugo',
  'Euphrasia_sp' = 'Euphrasia_sp_cf',
  'Trifolium_montanum' = 'Trifolium_cf_montanum',
  'Erigeron_annuus' = 'Erigernon_annuus',
  'Medicago_varia' = "Medicago_x_varia",
  'Bromus_hordeaceus_aggr_incl_B_commutatus' = 'Bromus_hordeaceus_aggr.incl_B_commutatus'
)]
Species_table_melt2 = rbind(Species_table_melt2, Species_duplicates)
Species_table_melt2$verbatimTraitName = Species_table_melt2$traitName
Species_table_melt2[, traitName := recode(
  verbatimTraitName,
  'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included' = 'SLA_with_petiole',
  'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded' = 'SLA_all',
  'Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC' = 'LDMC',
  'Leaf_nitrogen_N_content_per_leaf_dry_mass' = 'LeafN',
  'Leaf_phosphorus_P_content_per_leaf_dry_mass' = 'LeafP',
  'Mycorrhiza_type' = 'Myco_type',
  'Mycorrhizal_infection_intensity' = 'Myco_intensity',
  'Plant_height_vegetative' = 'Height',
  'Seed_dry_mass' = 'Seed_mass',
  'Stem_specific_density_SSD_or_wood_density_stem_dry_mass_per_stem_fresh_volume' = 'SSD'
)]

Species_table_melt2[, traitID := recode(
  traitName,
  'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included' = '3116',
  'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded' = '3117',
  'Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC' = '47',
  'Leaf_nitrogen_N_content_per_leaf_dry_mass' = '14',
  'Leaf_phosphorus_P_content_per_leaf_dry_mass' = '15',
  'Mycorrhiza_type' = '7',
  'Mycorrhizal_infection_intensity' = '1030',
  'Plant_height_vegetative' = '3106',
  'Seed_dry_mass' = '26',
  'Stem_specific_density_SSD_or_wood_density_stem_dry_mass_per_stem_fresh_volume' = '4',
  'Specific_root_length' = '1080',
  'Root_weight_ratio' = '9',
  'Rooting_depth' = '6',
  'Fine_roots_diameter' = '896',
  'Root_tissue_density' = '82'
)]
Species_table_melt2[, traitUnit := recode(
  traitName,
  'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included' = 'mm2/mg',
  'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded' = 'mm2/mg1 ',
  'Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC' = ' g/g',
  'Leaf_nitrogen_N_content_per_leaf_dry_mass' = 'mg/g',
  'Leaf_phosphorus_P_content_per_leaf_dry_mass' = 'mg/g',
  'Mycorrhiza_type' = '',
  'Mycorrhizal_infection_intensity' = '%',
  'Plant_height_vegetative' = 'm',
  'Seed_dry_mass' = 'mg',
  'Stem_specific_density_SSD_or_wood_density_stem_dry_mass_per_stem_fresh_volume' = 'g/cm3',
  'Specific_root_length' = 'm/g',
  'Root_weight_ratio' = 'unitless',
  'Rooting_depth' = 'cm',
  'Fine_roots_diameter' = 'mm',
  'Root_tissue_density' = 'mg/cm3'
)]
Species_table_melt2[Forest_or_grassland_estimate == 'Grassland', Comment := 'Use this estimate in GRASSLANDS']
Species_table_melt2[Forest_or_grassland_estimate == 'Forest', Comment := 'Use this estimate in FORESTS']
Species_table_melt2 = Species_table_melt2[!is.na(traitValue2), ]

# Need to add taxonID and scientificNameverbatim to fit the standard
addID = function(All_names, traitID) {
  All_names = gsub('=', '', All_names)
  All_names = gsub('%', '', All_names)
  All_names = gsub('[0-9]*', '', All_names)
  all_names = unlist(strsplit(All_names, '___'))
  data = TRYdata[(AccSpeciesName %in% all_names |
                    SpeciesName %in% all_names) &
                   !is.na(StdValue) &
                   traitID == traitID, list(
                     verbatimScientificName = unique(AccSpeciesName),
                     taxonID = unique(AccSpeciesID)
                   )]
  return(list(
    taxonID = paste(data$taxonID, collapse = '___'),
    verbatimScientificName = paste(data$verbatimScientificName, collapse = '___')
  ))
}

Species_table_melt2[, c('taxonID', 'verbatimScientificName') := addID(All_names, traitID), by = All_names]


fwrite(Species_table_melt2, "Plant_traits_final_details_all.csv")
fwrite(Species_table_melt2[!(
  traitName %in% c(
    'Specific_root_length',
    'Root_weight_ratio',
    'Rooting_depth',
    'Fine_roots_diameter',
    'Root_tissue_density'
  )
),],
"Plant_traits_final_details_above_ground.csv")

# Make it clean
Species_trait_data_clean = unique(Species_table_melt2[, list(
  scientificName = Species_list_format,
  #   verbatimScientificName,
  #   taxonID,
  traitName = traitName,
  # verbatimTraitName = traitName,
  traitID ,
  traitValue ,
  traitSd ,
  traitRange = Range,
  traitUnit,
  Comment = Comment,
  numberDatasets ,
  numberReplicates ,
  References
)])

Species_trait_data_clean[traitRange == '', traitRange := 'NA']
Species_trait_data_clean[numberDatasets == '', numberDatasets := 'NA']
Species_trait_data_clean[numberReplicates == '', numberReplicates := 'NA']

unique(Species_trait_data_clean$traitID)

fwrite(Species_trait_data_clean,
       "MAIN_Plant_traits_final_all.csv",
       sep = ";")
fwrite(Species_trait_data_clean[!(
  traitName %in% c(
    'Specific_root_length',
    'Root_weight_ratio',
    'Rooting_depth',
    'Fine_roots_diameter',
    'Root_tissue_density'
  )
),],
"MAIN_Plant_traits_final_all_TRY.csv",
sep = ";")




# *********************************** #
#### 2. Match with abundance data  ####
# *********************************** #
#Species_trait_data_clean = fread( "MAIN_Plant_traits_final_all.csv", sep = ";")

Species_trait_data_clean = fread("MAIN_Plant_traits_final_all.csv", sep = ";")

# Individual PCAs
Data_cast = dcast.data.table(unique(Species_trait_data_clean[Comment != "Use this estimate in FORESTS", ]),
                             scientificName ~ traitName,
                             value.var = 'traitValue')
Data_cast[, c(
  "Fine_roots_diameter",
  "Height",
  "LDMC" ,
  "LeafN",
  "LeafP",
  "Myco_intensity",
  "Root_tissue_density" ,
  "Root_weight_ratio"  ,
  "Rooting_depth"   ,
  "SLA_all"    ,
  "SLA_with_petiole"   ,
  "SSD"     ,
  "Seed_mass"     ,
  "Specific_root_length"
) := lapply(.SD, as.numeric), .SDcols = c(
  "Fine_roots_diameter",
  "Height",
  "LDMC" ,
  "LeafN",
  "LeafP",
  "Myco_intensity",
  "Root_tissue_density" ,
  "Root_weight_ratio"  ,
  "Rooting_depth"   ,
  "SLA_all"    ,
  "SLA_with_petiole"   ,
  "SSD"     ,
  "Seed_mass"     ,
  "Specific_root_length"
)]
Data_cast[Myco_type == 'Arbuscular mycorrhizaNon mycorrhizal', Myco_type := 'Non mycorrhizalArbuscular mycorrhiza']
Data_cast[Myco_type %in% c(
  'Mixed',
  'Non mycorrhizalArbuscular mycorrhiza',
  'Arbuscular mycorrhizaOther or undefined endomycorrhiza'
), Myco_type := 'Mixed']
Data_cast[, Myco_type := factor(Myco_type)]
pca_AG = dudi.pca(mice::complete(mice(Data_cast[, c("Height",
                                                    "LDMC" ,
                                                    "LeafN",
                                                    "LeafP",
                                                    "SLA_with_petiole" ,
                                                    "SSD"  ,
                                                    "Seed_mass")])))
draw_dudi_mix(
  Data = pca_AG,
  c(1, 3),
  save = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Ind_level/PlantsAG13.pdf',
  new_var_names = c("Height"   ,        "LDMC"    ,         "LeafN"        ,    "LeafP"      ,      "SLA" , "SSD"  ,  "Seed_mass")
)

pca_BG = dudi.mix(Data_cast[complete.cases(Data_cast[, c(
  "Fine_roots_diameter",
  "Myco_intensity",
  "Myco_type",
  "Root_tissue_density" ,
  "Root_weight_ratio"  ,
  "Rooting_depth"   ,
  "Specific_root_length"
)])
, c(
  "Fine_roots_diameter",
  "Myco_intensity",
  "Myco_type",
  "Root_tissue_density" ,
  "Root_weight_ratio"  ,
  "Rooting_depth"   ,
  "Specific_root_length"
)],
scannf = FALSE, nf = 3)
draw_dudi_mix(
  Data = pca_BG,
  c(1, 3),
  save = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Ind_level/PlantsBG13.pdf',
  new_var_names = c(
    "Fine_roots_diameter" ,
    "Myco_intensity",
    "AM",
    "MixedM",
    "NoM",
    "OtherM",
    "Root_tissue_density"  ,
    "Root_weight_ratio",
    "Rooting_depth" ,
    "Specific_root_length"
  )
)

### Homogenise species name ###
Abundances = rbind(Abundances_grasslands[, list(
  Plot = Plot,
  Species = Species,
  Year = Year,
  Cover = value
)],
Abundances_forests[#Layer %in% c('H', 'S') %***% Uncomment and complete if you want to restrict the calculation of CWM in forests to some layers
  , list(Cover = sum(Cover)), by = list('Plot' = Useful_EPPlotID, Year, Species)])

# Check if all species are in the trait dataset
unique(Abundances[!(Species %in% Species_trait_data_clean$scientificName), Species])
# -> ok, other species we couldn't find in TRY or not identified enough

# ************************************ #
#### 3. Calculate coverage and CWM  ####
# ************************************ #
# Functions
library(FD)

# In some cases we might want to exclude tree saplings from grassland CWM measures, so here is a list of some tree species found in grasslands
Species_trait_data_clean[traitName == 'Height' &
                           as.numeric(traitValue) > 5, unique(scientificName)]
list_trees = c(
  "Acer_campestre",
  "Acer_platanoides" ,
  "Acer_pseudoplatanus",
  "Acer_sp",
  "Alnus_glutinosa",
  "Betula_pendula",
  "Carpinus_betulus",
  "Crataegus_monogyna",
  "Crataegus_sp" ,
  "Fagus_sylvatica",
  "Fraxinus_excelsior"  ,
  "Malus_sylvestris" ,
  "Picea_abies" ,
  "Pinus_sylvestris",
  "Populus_sp" ,
  "Populus_tremula" ,
  "Prunus_avium"   ,
  "Prunus_sp" ,
  "Quercus_petraea"   ,
  "Quercus_robur" ,
  "Quercus_rubra"  ,
  "Quercus_sp"    ,
  "Sorbus_aria" ,
  "Sorbus_aucuparia" ,
  "Sorbus_torminalis",
  "Tilia_platyphyllos"  ,
  "Tilia_sp"       ,
  "Ulmus_minor"      ,
  "Crataegus_laevigata"    ,
  "Abies_alba"     ,
  "Aesculus_hippocastanum",
  "Alnus_incana" ,
  "Larix_decidua" ,
  "Malus_sp"   ,
  "Prunus_padus" ,
  "Prunus_serotina" ,
  "Pseudotsuga_menziesii" ,
  "Pyrus_communis",
  "Pyrus_pyraster" ,
  "Robinia_pseudoacacia",
  "Salix_caprea" ,
  "Salix_cinerea",
  "Salix_sp"  ,
  "Taxus_baccata",
  "Tilia_cordata"  ,
  "Ulmus_glabra" ,
  "Ulmus_laevis"  ,
  "Alnus_sp"    ,
  "Betula_sp"    ,
  "Larix_sp"  ,
  "Populus_nigra"  ,
  "Ulmus_sp" ,
  "Carya_ovata"  ,
  "Pinus_mugo"   ,
  "Pinus_cf_mugo" ,
  'Tree_seedling'
)
list_bushes = c("Hedera_helix" ,
                "Juniperus_communis",
                'Clematis_vitalba',
                'Sambucus_nigra')

# This function reformats abundance and trait datasets (e.g. by selecting species with abundance > 0
# and pooling functionally redundant species)
# These "functionally redundant" species are species for which all traits estimates are the same (for instance, Acer sp.
# and Acer_campestre if Acer sp. estimates are based on Acer_campestre only)
# Then applies the dbFD function
wrapper_dbFD = function(Ab_data0,
                        Trait_data0,
                        trait_names,
                        exclude_species = '',
                        threshold = 0.8,
                        method_threshold = 'per_plot') {
 
  Trait_data0$traitValue = as.numeric(Trait_data0$traitValue)
  
   # Both Abundance and Trait data must be in the long format
  Ab_data = copy(Ab_data0)
  Trait_data = copy(Trait_data0)
  setDT(Trait_data)
  setDT(Ab_data)
  
  ### 0. We first subset the species based on abundance, to get rid of rare/weird species that break the dbFD function
  # First method: we keep most abundant species OVERALL, to reach an overall abundance of x%
  if (!(method_threshold %in% c('overall', 'per_plot'))) {
    print("method_threshold must be 'overall' or 'per_plot'")
    stop()
  }
  Abb = Ab_data[Cover > 0,]
  if (method_threshold == 'overall') {
    Abb = Abb[,  sum(Cover, na.rm = T), by = c('Species')][order(V1, decreasing = T), ]
    Abb[, prop := V1 / sum(V1)]
    Abb[, cumprop := cumsum(prop)]
    most_common = Abb[cumprop < threshold, unique(Species)]
  }
  # Second method: we keep most abundant species within each plot, to reach an overall abundance of at least x% in every plot
  if (method_threshold == 'per_plot') {
    Abb = Abb[,  sum(Cover, na.rm = T), by = c('Species', 'Plot', 'Year')][order(Plot, V1, decreasing = T), ]
    Abb[, prop := V1 / sum(V1), by = c('Plot',  'Year')]
    Abb[, cumprop := cumsum(prop), by = c('Plot', 'Year')]
    most_common = Abb[cumprop < threshold, unique(Species)]
  }
  
  
  ### 1. Calculation of CWM and FD ###
  # First we take the common species between the abundance and trait data
  common_sp = intersect(intersect(Ab_data$Species, Trait_data$scientificName),
                        most_common)
  # Reformat datasets
  Trait_data_cast = dcast.data.table(Trait_data[scientificName %in% common_sp, ], scientificName ~ traitName, value.var = 'traitValue')
  Trait_data_cast[, (trait_names[trait_names != "Mycorrhiza_type"]) := lapply(.SD, as.numeric), .SDcols = trait_names[trait_names != "Mycorrhiza_type"]]
  Ab_data_cast = dcast.data.table(Ab_data[Species %in% common_sp, ],
                                  Plot + Year ~ Species,
                                  value.var = 'Cover',
                                  fill = 0)
  # We make sure that both datasets data had no rows (Traits) or rows or columns (Abundance) with only NAs
  Trait_data = Trait_data_cast[Trait_data_cast[,!all(is.na(.SD)), .SDcols = trait_names, by = 1:nrow(Trait_data_cast)]$V1, ]
  Ab_data_cast = Ab_data_cast[(rowSums(Ab_data_cast[,-c(1, 2)], na.rm = T) > 0) , ]
  species_not_0 = colSums(Ab_data_cast[,-c(1, 2)], na.rm = T) > 0
  species_not_0 = names(species_not_0[species_not_0 == TRUE])
  Ab_data_cast = Ab_data_cast[, .SD, .SDcols = colnames(Ab_data_cast)[colnames(Ab_data_cast) %in% c('Plot', 'Year', species_not_0)]]
  # Then we re-check that all species are present at least once
  common_sp1 = intersect(Trait_data$scientificName, colnames(Ab_data_cast))
  Trait_data_cast = Trait_data_cast[scientificName %in% common_sp1, ]
  Ab_data_cast = Ab_data_cast[, .SD, .SDcols = colnames(Ab_data_cast)[colnames(Ab_data_cast) %in% c('Plot', 'Year', common_sp1)]]
  # We have to get rid of species that have only one trait data or it crashes the dbFD function
  one_trait_species = Trait_data[Trait_data[, rowSums(!is.na(.SD)) < 2, .SDcols = trait_names], scientificName]
  common_sp2 = common_sp1[!(common_sp1 %in% one_trait_species)]
  # We can also exclude some additional species if needed
  common_sp_FD = common_sp2[!(common_sp2 %in% exclude_species)]
  
  # One of the issues here is that some species are functionally redundant, ie they have same trait values.
  # So we have to identify these species and aggregate them in the Abundance matrix
  Trait_data_cast[, (trait_names) := lapply(.SD, function(x) {
    signif(x, 5)
  }), .SDcols = trait_names]
  Trait_data_cast[, duplicated := duplicated(.SD), .SDcols = trait_names]
  # We identify duplicate rows and pool the corresponding rows in the Abundance dataset
  duplicated_species = unlist(Trait_data_cast[duplicated == TRUE, scientificName]) # Duplicate species ID
  for (dup_sp in duplicated_species) {
    trait_values = unlist(Trait_data_cast[scientificName == dup_sp, .SD, .SDcols = trait_names])
    trait_values = sort(trait_values, na.last = TRUE)
    if (length(trait_values[!is.na(trait_values)]) == 1)
      next # If there is only one trait value, the species will be removed at the next step
    l = length(trait_values)
    same_traits_species = Trait_data_cast[(get(names(trait_values)[1]) == trait_values[1] |
                                             is.na(trait_values[1]))  &
                                            (get(names(trait_values)[2]) == trait_values[2] |
                                               is.na(trait_values[2])) &
                                            (get(names(trait_values)[3]) == trait_values[3]  |
                                               is.na(trait_values[3])) &
                                            (get(names(trait_values)[4]) == trait_values[4]  |
                                               is.na(trait_values[4])) &
                                            (get(names(trait_values)[5]) == trait_values[5] |
                                               is.na(trait_values[5])), unique(scientificName)]
    # same_traits_species = Trait_data_cast[Trait_data_cast[, Reduce(`&`, lapply(.SD, `%in%`, trait_values[!is.na(trait_values)])), .SDcols = trait_names] , unique(scientificName)]
    same_traits_species = same_traits_species[same_traits_species != dup_sp]
    Ab_data_cast[, (dup_sp) := rowSums(.SD), .SDcols = c(dup_sp, same_traits_species)]
    Ab_data_cast = Ab_data_cast[, .SD, .SDcols = colnames(Ab_data_cast)[!(colnames(Ab_data_cast) %in% same_traits_species)]]
    Trait_data_cast = Trait_data_cast[!(scientificName %in% same_traits_species), ]
  }
  Ab_data_cast = data.frame(Ab_data_cast)
  rownames(Ab_data_cast) = paste(Ab_data_cast$Plot, Ab_data_cast$Year, sep = '_')
  Ab_data_cast = Ab_data_cast[,-c(1, 2)]
  
  Trait_data_cast = data.frame(Trait_data_cast)
  rownames(Trait_data_cast) = Trait_data_cast$scientificName
  Trait_data_cast = Trait_data_cast[, trait_names]
  Trait_data_cast = Trait_data_cast[order(rownames(Trait_data_cast)), ]
  Ab_data_cast = Ab_data_cast[, order(colnames(Ab_data_cast))]
  
  # We calculate CWM on all species, but FD on the restricted species list that we just created
  CWM_output = functcomp(x = as.matrix(Trait_data_cast), a = as.matrix(Ab_data_cast))
  
  Ab_data_cast = Ab_data_cast[, intersect(colnames(Ab_data_cast), common_sp_FD)][colSums(Ab_data_cast[, intersect(colnames(Ab_data_cast), common_sp_FD)]) >
                                                                                   0, ]
  Ab_data_cast = Ab_data_cast[rowSums(Ab_data_cast, na.rm = T) > 0, ]
  
  dbFD_output = dbFD(
    x = as.matrix(Trait_data_cast[colnames(Ab_data_cast), ]),
    a = as.matrix(Ab_data_cast),
    corr = 'lingoes',
    stand.FRic = T,
    calc.CWM = F,
    m = 5
  )
  
  ### 2. Calculation of coverage (i.e. how many individuals are actually included in the calculation) ###
  # I have to start from scratch as the total abundance must be with all species, not on the restricted list
  #  I'm cheating a bit: I replace all non-na traits values by 1 and NAs by ~0
  # Then calculating CWM on these traits provides the approximate coverage for each trait
  Ab_data = copy(Ab_data0)
  Trait_data = copy(Trait_data0)
  Trait_data_cast_cov = dcast.data.table(Trait_data, scientificName ~ traitName, value.var = 'traitValue')
  Trait_data_cast_cov[, (trait_names[trait_names != "Mycorrhiza_type"]) := lapply(.SD, as.numeric), .SDcols = trait_names[trait_names != "Mycorrhiza_type"]]
  Ab_data_cast_cov = dcast.data.table(Ab_data,
                                      Plot + Year ~ Species,
                                      value.var = 'Cover',
                                      fill = 0)
  not_0_sp = names(colSums(Ab_data_cast_cov[,-c(1, 2)], na.rm = T)[colSums(Ab_data_cast_cov[,-c(1, 2)], na.rm = T) > 0])
  Ab_data_cast_cov = data.frame(Ab_data_cast_cov[, .SD, .SDcols = c('Plot', 'Year', not_0_sp)])
  rownames(Ab_data_cast_cov) = paste(Ab_data_cast_cov$Plot, Ab_data_cast_cov$Year, sep = '_')
  Ab_data_cast_cov = Ab_data_cast_cov[,-c(1, 2)]
  # We add all species that are included in Ab_data_cast but not Trait_data_cast as NAs in Trait_data_cast
  no_trait_species = not_0_sp[!not_0_sp %in% Trait_data_cast_cov$scientificName]
  No_trait_species = data.table(scientificName = no_trait_species)
  No_trait_species[, (colnames(Trait_data_cast_cov)[-1]) := NA]
  Trait_data_cast_cov2 = data.frame(rbind(Trait_data_cast_cov, No_trait_species))
  rownames(Trait_data_cast_cov2) = Trait_data_cast_cov2$scientificName
  Trait_data_cast_cov2 = Trait_data_cast_cov2[not_0_sp, trait_names]
  Trait_data_cast_cov2 = Trait_data_cast_cov2[order(rownames(Trait_data_cast_cov2)), ]
  Ab_data_cast_cov = Ab_data_cast_cov[, order(colnames(Ab_data_cast_cov))]
  
  # All NAs are 0, all the others are 1
  Trait_data_cast_cov2[!is.na(Trait_data_cast_cov2)] = 1.000000000000000000000000001
  Trait_data_cast_cov2[is.na(Trait_data_cast_cov2)] =  -0.000000000000000000000000001
  # First for CWM - we keep only the species that were counted in the CWM calculation, for the other traits are set to 0
  Trait_data_cast_cov2[!(rownames(Trait_data_cast_cov2) %in% rownames(Trait_data_cast)), ] = -0.000000000000000000000000001
  coverage_cwm = functcomp(a = as.matrix(Ab_data_cast_cov),
                           x = as.matrix(Trait_data_cast_cov2))
  # Then for FD - we keep only the species that were counted in the FD calculation, for the other traits are set to 0
  Trait_data_cast_cov2[!(rownames(Trait_data_cast_cov2) %in% rownames(Trait_data_cast)), ] = -0.000000000000000000000000001
  Trait_data_cast_cov2[!(rownames(Trait_data_cast_cov2) %in% common_sp_FD), ] = -0.000000000000000000000000001
  
  #Trait_data_cast_cov2[!(rownames(Trait_data_cast_cov2) %in% most_common) & abs(rowSums(Trait_data_cast_cov2)) > 0.000000000000000000000001,]
  coverage_fd = functcomp(a = as.matrix(Ab_data_cast_cov),
                          x = as.matrix(Trait_data_cast_cov2))
  output = list(
    CWM = data.frame(CWM_output),
    FD = data.frame(
      FRic = dbFD_output$FRic,
      FEve = dbFD_output$FEve,
      FDiv = dbFD_output$FDiv,
      FDis = dbFD_output$FDis,
      RaoQ = dbFD_output$RaoQ
    ),
    Coverage_CWM = data.frame(coverage_cwm),
    Coverage_FD = data.frame(coverage_fd)
  )
  return(output)
}

## This is kept here for debugging purpose: in case some distances are found to be NAs, run this
## on the Trait_data_cast (after excluding necessary species) to check which species pairs cause issues
dat = Trait_data_cast[intersect(rownames(Trait_data_cast), common_sp_FD), ]
for (i in 1:(nrow(Trait_data_cast[intersect(rownames(Trait_data_cast), common_sp_FD), ]) -
             1)) {
  for (j in (i + 1):nrow(Trait_data_cast[intersect(rownames(Trait_data_cast), common_sp_FD), ])) {
    gwd = cluster::daisy(Trait_data_cast[intersect(rownames(Trait_data_cast), common_sp_FD), ][c(i, j), ], "gower")
    if (is.na(gwd) | is.nan(gwd) | (gwd == 0)) {
      print(gwd)
    }
  }
}


# In grasslands,  Above-ground
CWM_grasslands_AG_no_tree = wrapper_dbFD(
  Trait_data0 =  unique(Species_trait_data_clean[!grepl('FOREST', Comment), .SD, .SDcols = colnames(Species_trait_data_clean)[colnames(Species_trait_data_clean) != 'traitSd']]),
  Ab_data0 = Abundances[grepl('G', Plot)
                        &
                          !(Species %in% c(list_trees,  "Juniperus_communis")) #%***% Comment here if you want to keep trees in the calculation of the CWM
                        , ],
  trait_names = c(
    "SLA_all",
    "SLA_with_petiole",
    "LDMC",
    "LeafN",
    "LeafP",
    "Height",
    "Seed_mass",
    "SSD"
  ),
  threshold = 0.99,
  'per_plot'
)


# In forests, Above-ground
# Note: I have to remove 'Potentilla_sterilis' here as it crashes the function (that's because it shares no non-na trait
# with another species, 'Poa_chaixii' - might be corrected in the future).
# This species has an overall relative abundance of 0.086 % -> probably negligible
CWM_forests_AG = wrapper_dbFD(
  Trait_data0 = unique(Species_trait_data_clean[!grepl('GRASS', Comment), .SD, .SDcols = colnames(Species_trait_data_clean)[colnames(Species_trait_data_clean) != 'traitSd']]),
  Ab_data0 = Abundances[grepl('W', Plot),],
  exclude_species = 'Potentilla_sterilis',
  trait_names = c(
    "SLA_all",
    "SLA_with_petiole",
    "LDMC",
    "LeafN",
    "LeafP",
    "Height",
    "Seed_mass",
    "SSD"
  ),
  threshold = 1,#0.99,
  'per_plot'
)
# Some plots without FD data? --> Mostly Fagus forests with only 1 or 2 species

# In grasslands, Below-ground
# Note: I have to remove 'Persicaria_lapathifolium', 'Juncus_inflexus', 'Juncus_effusus' and 'Geranium_dissectum' here as they crashes the function (that's because
# they have only two available traits)
# These species have an overall relative abundance of 0.0001 %, 0.0046 %, 0.059 %  and 0.065 %-> probably negligible
CWM_grasslands_BG_no_tree = wrapper_dbFD(
  Trait_data0 = unique(Species_trait_data_clean[!grepl('FOREST', Comment), .SD, .SDcols = colnames(Species_trait_data_clean)[colnames(Species_trait_data_clean) != 'traitSd']]),
  Ab_data0 = Abundances[grepl('G', Plot)
                        &
                          !(Species %in% c(list_trees,  "Juniperus_communis")) #%***% Comment here if you want to keep trees in the calculation of the CWM
                        , ],
  exclude_species = c(
    'Persicaria_lapathifolium',
    'Juncus_inflexus',
    'Juncus_effusus',
    'Geranium_dissectum',
    'Chaerophyllum_aureum',
    'Cruciata_laevipes'
  ),
  trait_names = c(
    "Specific_root_length",
    "Root_weight_ratio",
    "Rooting_depth",
    "Fine_roots_diameter",
    "Root_tissue_density",
    'Mycorrhizal_inf_int'
  ),
  threshold = 0.99,
  'per_plot'
)

# In forests, Below-ground
CWM_forests_BG = wrapper_dbFD(
  Trait_data0 = unique(Species_trait_data_clean[!grepl('GRASS', Comment), .SD, .SDcols = colnames(Species_trait_data_clean)[colnames(Species_trait_data_clean) != 'traitSd']]),
  Ab_data0 = Abundances[grepl('W', Plot),],
  exclude_species = c(
    'Persicaria_lapathifolium',
    'Juncus_inflexus',
    'Juncus_sp',
    'Juncus_effusus',
    'Geranium_dissectum',
    'Chaerophyllum_aureum',
    'Cruciata_laevipes'
  ),
  trait_names = c(
    "Specific_root_length",
    "Root_weight_ratio",
    "Rooting_depth",
    "Fine_roots_diameter",
    "Root_tissue_density",
    'Mycorrhizal_inf_int'
  ),
  threshold = 0.99,
  'per_plot'
)
#### Export data ###
reformat_data = function(data, index = 'CWM') {
  if (!(index %in% c('CWM', 'FD'))) {
    print("Index should be 'CWM' or 'FD'")
    stop()
  }
  if (index == 'CWM') {
    print("calculating CWM")
    trait_data = data.table(data$CWM)
    cov_data = data.table(data$Coverage_CWM)
    trait_data[, c('Plot', 'Year') := list(substr(rownames(data$CWM), 1, 5), substr(rownames(data$CWM), 7, 11))]
    cov_data[, c('Plot', 'Year') := list(substr(rownames(data$Coverage_CWM), 1, 5), substr(rownames(data$Coverage_CWM), 7, 11))]
  }
  if (index == 'FD')  {
    print("calculating FD")
    trait_data = data.table(data$FD)
    cov_data = data.table(data$Coverage_FD)
    cov_data[, c('Plot', 'Year') := list(substr(rownames(data$Coverage_FD), 1, 5), substr(rownames(data$Coverage_FD), 7, 11))]
    trait_data[, c('Plot', 'Year') := list(substr(rownames(data$FD), 1, 5), substr(rownames(data$FD), 7, 11))]
    
  }
  
  if (index == 'CWM') {
    trait_data_melt = melt.data.table(
      trait_data,
      id.vars = c('Plot', 'Year'),
      value.name = index,
      variable.name = "traitName"
    )
    trait_data_melt[, 'CWM' := ifelse(CWM>1000, signif(CWM, 4),signif(CWM, 3))]
    cov_data_melt = melt.data.table(
      cov_data,
      id.vars = c('Plot', 'Year'),
      value.name = 'Coverage',
      variable.name = "traitName"
    )
    cov_data_melt[, Coverage := signif(Coverage, 3)]
  }
  
  if (index == 'FD') {
    trait_data_melt = trait_data[, lapply(.SD, function(x) {
      signif(x, 3)
    }), .SDcols = colnames(trait_data)[!(colnames(trait_data) %in% c('Plot', 'Year'))], by = c('Plot', 'Year')]
    cov_data_melt = cov_data[, list(Coverage_min = signif(min(.SD), 3)), by = c('Plot', 'Year')]
    print(cov_data_melt)
    print(trait_data_melt)
    
  }
  
  print("Merging")
  
  cwm_cov = merge.data.table(trait_data_melt, cov_data_melt)
  print(cwm_cov)
  print("Adding columns")
  
  if (index == 'CWM') {
    cwm_cov[, traitName := recode(traitName, 'SLA_all' = 'SLA')]
    cwm_cov[, verbatimTraitName := recode(
      traitName,
      'SLA_with_petiole' = 'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included' ,
      'SLA' = 'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded'  ,
      'LDMC' = 'Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC'  ,
      'LeafN' = 'Leaf_nitrogen_N_content_per_leaf_dry_mass'  ,
      'LeafP' =  'Leaf_phosphorus_P_content_per_leaf_dry_mass'  ,
      'Myco_type' = 'Main_mycorrhiza_type' ,
      'Myco_intensity' = 'Mycorrhizal_infection_intensity'  ,
      'Mycorrhizal_inf_int' = 'Mycorrhizal_infection_intensity'  ,
      'Height' = 'Plant_height_vegetative' ,
      'Seed_mass' = 'Seed_dry_mass',
      'SSD' = 'Stem_specific_density_SSD_or_wood_density_stem_dry_mass_per_stem_fresh_volume'
    )]
    
    cwm_cov[, traitID := recode(
      verbatimTraitName,
      'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included' = '3116',
      'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded' = '3117',
      'Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC' = '47',
      'Leaf_nitrogen_N_content_per_leaf_dry_mass' = '14',
      'Leaf_phosphorus_P_content_per_leaf_dry_mass' = '15',
      'Main_mycorrhiza_type' = '7',
      'Mycorrhizal_infection_intensity' = '1030',
      'Plant_height_vegetative' = '3106',
      'Seed_dry_mass' = '26',
      'Stem_specific_density_SSD_or_wood_density_stem_dry_mass_per_stem_fresh_volume' = '4',
      'Specific_root_length' = '1080',
      'Root_weight_ratio' = '9',
      'Rooting_depth' = '6',
      'Fine_roots_diameter' = '896',
      'Root_tissue_density' = '82'
    )]
    
    cwm_cov[, traitUnit := recode(
      verbatimTraitName,
      'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included' = 'mm2/mg',
      'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded' = 'mm2/mg',
      'Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC' = 'g/g',
      'Leaf_nitrogen_N_content_per_leaf_dry_mass' = 'mg/g',
      'Leaf_phosphorus_P_content_per_leaf_dry_mass' = 'mg/g',
      'Main_mycorrhiza_type' = '',
      'Mycorrhizal_infection_intensity' = '%/100',
      'Plant_height_vegetative' = 'm',
      'Seed_dry_mass' = 'mg',
      'Stem_specific_density_SSD_or_wood_density_stem_dry_mass_per_stem_fresh_volume' = 'g/cm3',
      'Specific_root_length' = 'm/g',
      'Root_weight_ratio' = 'unitless',
      'Rooting_depth' = 'cm',
      'Fine_roots_diameter' = 'mm',
      'Root_tissue_density' = 'mg/cm3'
    )]
  }
  print('Exporting')
  
  cwm_cov$EP_PlotID = gsub('0([1-9])', '\\1', cwm_cov$Plot)
  
  # Excluding values with < 10% coverage
  cwm_cov[Coverage < 0.00000001, Coverage := 0]
  cwm_cov[Coverage < 0.01, CWM := 'NA']

  return(cwm_cov)
}

# Without trees in grasslands
all_CWM_AG_no_tree = rbind(reformat_data(CWM_grasslands_AG_no_tree),
                           reformat_data(CWM_forests_AG))
fwrite(all_CWM_AG_no_tree, "CWM_Plants_AG_no_tree.csv")

all_CWM_BG_no_tree = rbind(reformat_data(CWM_grasslands_BG_no_tree),
                           reformat_data(CWM_forests_BG))
fwrite(all_CWM_BG_no_tree, "CWM_Plants_BG_no_tree.csv")


all_FD_AG_no_tree = rbind(
  reformat_data(CWM_grasslands_AG_no_tree, 'FD'),
  reformat_data(CWM_forests_AG, 'FD')
)
fwrite(all_FD_AG_no_tree, "FD_Plants_AG_no_tree.csv")

all_FD_BG_no_tree = rbind(
  reformat_data(CWM_grasslands_BG_no_tree, 'FD'),
  reformat_data(CWM_forests_BG, 'FD')
)
fwrite(all_FD_BG_no_tree, "all_FD_BG_no_tree.csv")

all_CWM_no_tree = rbind(all_CWM_AG_no_tree, all_CWM_BG_no_tree)
all_CWM_no_tree = all_CWM_no_tree[grepl('G', Plot), ]
fwrite(all_CWM_no_tree, "all_CWM_no_tree.csv")

all_FD_no_tree = cbind(rbind(all_FD_AG_no_tree, all_FD_BG_no_tree),
                       rep(
                         c('Above-ground traits', 'Below-ground traits'),
                         c(nrow(all_FD_AG_no_tree), nrow(all_FD_BG_no_tree))
                       ))
fwrite(all_FD_no_tree, "all_FD_no_tree.csv")


CWM_test = reformat_data(CWM_grasslands_for_Antonios)
FD_test = reformat_data(CWM_grasslands_for_Antonios, 'FD')

fwrite(FD_test, "all_FD_no_tree_for_Antonios_NEW.csv")

# Quick check
