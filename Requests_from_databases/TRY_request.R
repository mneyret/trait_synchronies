############ Create, match species list to TRY #######
library(Taxonstand)
library(readxl)
library(readr)
library(tpl) 
library(flora)
library(dplyr)
library(beepr)
library(data.table)


### These are all the different dataset...
BEFUp_list <- read_excel("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/SP_LIST_EP_BEFUp.xlsx")
Forest_veg_list <- read_delim("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/25886_Vegetation Records for 151 Forest EPs, 2009 - 2018_1.5.4/25886.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)
Forest_stand_list <- read_delim("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/18269_Stand composition, abundance on all forest EPs, 2008 – 2014_1.5.14/18269.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)
Forest_stand_list2 <- read_delim("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/22907_Stand composition based on 2nd forest inventory (abundance, basal area and crown projection area) on all forest EPs, 2014 – 2018_1.4.6/22907.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)
Grassland_list <- read_delim("~/Desktop/Research/Senckenberg/Project_Landscape_MF/Landscape_composition/Data/Raw_data/Pour Margot/190411_EP_species_info_GRL_140519.txt", 
                           ";", escape_double = FALSE, trim_ws = TRUE)
Dataset_clean <- fread("~/Desktop/Research/Senckenberg/Project_Landscape_MF/Landscape_composition/Data/Raw_data/New_synthesis_data/Dataset_clean.txt")
Abundances = Dataset_clean[Group_broad == 'Plant',]
New_sp_Ralph <- read_excel("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/EPs_Relevee_2008_2016_2019_SpName_Ralph_Bolliger.xlsx")
Aggregates_species<- data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/Aggregates_species_from_Ralph.xlsx"))
Aggregates_species_melt = melt.data.table(Aggregates_species, id.vars = c('SpList_Margot', 'Aggregate'))
Aggregates_species_melt = Aggregates_species_melt[!(value %in% c('0', '1', NA, 'Species we probably saw and aggregat because are two similar vegetatif','Species we probably saw andgat because are two similar vegetatif', "Species we probably saw and aggregated because are two similar vegetatif"))]

# I also add additionnal species based on Ralph species matching if need later
additionnal_species = unique(c(New_sp_Ralph$OldVersion_2008_2016_Names_To_Reassign))
additionnal_species = additionnal_species[additionnal_species != 'Rumex_x_pratensis']
additionnal_species[!(additionnal_species %in% Species_table$Species_list_format)]
additionnal_species = gsub("_", ' ', additionnal_species)
additionnal_species[additionnal_species == "Persicaria vivipara (Bistorta vivipara in the dataset of 2016)"] = "Persicaria vivipara"
additionnal_species[additionnal_species == "Cerastium nutans (Cerastium glutinosum)"] = "Cerastium nutans"
additionnal_species = c(additionnal_species, "Cerastium glutinosum",
                        'Primula veris',
                        'Ononis repens',
                        'Rubus corylifolius',
                        'Viola riviniana',
                        'Bromus commutatus',
                        'Rhinanthus aggr.',
                        'Corydalis cava', 'Corydalis intermedia', 'Corydalis pumila',  'Corydalis solida', 'Pseudofumaria lutea',
                        'Lycopodium clavatum', 'Lycopodium × oellgaardii','Lycopodium tristachyum', 'Lycopodium × zeilleri', 'Lycopodium alpinum L.',
                        "Lycopodium annotinum L.", 'Lycopodium clavatum L.', 'Lycopodium complanatum L.', 'Lycopodium × issleri (Rouy) Domin')

# Notes
# These genus do not have any identified species, and so were not included in the TRY request. Maybe do a new request by looking most common species?
# Corydalis:
# Lycopodium: 

# ...  which I combine into one species list
Species_list = unique(c(BEFUp_list$SPECIES, additionnal_species, Forest_veg_list$Species, Forest_stand_list$Species,  Forest_stand_list2$Species, New_sp_Ralph$NewVersion_Dataset_26106_original_name, Grassland_list[Grassland_list$Group_broad == 'Plant',]$Species ))
Species_table = data.frame(Species_raw = Species_list)
setDT(Species_table)

Species_table[, Species_list_format := gsub('_', ' ', Species_raw)]


### This is the list of TRY species
TryAccSpecies <- read_delim("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/TryAccSpecies.txt", '\t')

TryAccSpecies$AccSpeciesName[grepl('Allium', TryAccSpecies$AccSpeciesName)]



### Try to match both lists. For this:
#      - Consider that aggregate species can be attributed the main species' traits
Species_table[, Species_list_format := gsub("aggr$", "aggr.", Species_list_format)]
Species_table[, Species_list_format := gsub("agg$", "aggr.", Species_list_format)]
Species_table[, Species_list_format := gsub("agg.$", "aggr.", Species_list_format)]
Species_table[, Species_list_format := gsub("aggr..", "aggr.", Species_list_format)]

#      - Homogenise "sp" and "spec"
Species_table[, Species_list_format :=  gsub(" spec", " sp", Species_list_format)]
Species_table[, Species_list_format :=  gsub(" sp\\.", " sp", Species_list_format)]

#      - Get rid of cf, except for Allium which is represented by a list of multiple species
Species_table[!grepl('Allium', Species_list_format), Species_list_format :=  gsub(" cf", "", Species_list_format)]

# Get rid of weird symbols
Species_table[, Species_list_format :=  gsub("×", "x", Species_list_format)]

Species_table[grepl('Bromus hordeaceus', Species_list_format) | grepl('Bromus_hordeaceus', Species_list_format), Species_list_format := "Bromus hordeaceus aggr incl B commutatus"]
Species_table[, Species_list_corrected := Species_list_format]

Aggregates_species_agg = tapply(Aggregates_species_melt$value, Aggregates_species_melt$SpList_Margot, function(x){paste(x, collapse = '___')}, simplify = TRUE)
for (aggr.sp in unique(Aggregates_species_melt[!is.na(SpList_Margot),]$SpList_Margot)){
  Species_table[Species_list_format == aggr.sp, Species_list_corrected := Aggregates_species_agg[[aggr.sp]]]
}

for (i in 1:nrow(Species_table)){
  if (is.null(Species_table$Species_list_corrected[i]) | is.na(Species_table$Species_list_corrected[i])){
    Species_table$Species_list_corrected[i] =  Species_table$Species_list_format[i]
  }}

#      - Correct some species names
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Neottia nidus avis', 'Neottia nidus-avis', x)})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Impatiens noli tangere', 'Impatiens noli-tangere', x)})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Buphtalmum salicifolium', 'Buphthalmum salicifolium', x)})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Silene flos.cuculi', 'Silene flos-cuculi', x)})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Athyrium filix femina', 'Athyrium filix-femina', x)})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Capsella bursa.pastoris', 'Capsella bursa-pastoris', x)})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Dryopteris filix mas', 'Dryopteris filix-mas', x)})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Mycelis muralis', 'Lactuca muralis', x)})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Persicaria lapathifolium', 'Persicaria lapathifolia', x)})]
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub("Valerianella officinalis", 'Valeriana officinalis', x)})] # I guess that was probably a typo 
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub("heliscopia", 'helioscopia', x)})] # I guess that was probably a typo 
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub("Alisma plantago-aquatica", 'Alisma plantago-aquatica subsp. orientale', x)})] # Only subsp. present in TRY, better than nothing
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub("Achillea distans", '"Achillea distans subsp. stricta___Achillea distans subsp. tanacetifolia"', x)})] # Only subsp. present in TRY, better than nothing
Species_table[, Species_list_corrected := sapply(Species_list_corrected, function(x){gsub('Taraxacum sect ruderalia', 'Taraxacum officinale', x)})]


Species_table[Species_list_format == 'Inula sp', 
              Species_list_corrected := "Inula conyzae___Inula helenium___Inula salicina___Inula hirta___Inula helvetica___Inula germanica___Inula britannica"] # I also add Inulas to the lot, taken all species from Rothmahler Flora


Species_table[, Species_list_corrected := gsub('[\"]*','', Species_list_corrected)]

# Remove unidentified species from the list
unidentifiable_species = c("unknown deciduous", "Tree seedling", "Lamiaceae sp", "Asteraceae sp", "Caryophyllaceae sp","Baumkeimling sp", 'Orchidaceae sp', 'Brassicaceae sp')
Species_table = Species_table[!(Species_list_format %in% unidentifiable_species),]


##### Get all synonyms from the plant list, but keep only those that are in TRY list
# Version 1: take synonyms only if the accepted name is not in TRY
get_all_synonyms_and_accepted1 = function(x, match_list){
  x = unlist(strsplit(as.character(x), '___'))
  names_ids = sapply(x, function(y){
    print(y)
    tpl_res = tpl.get(y, return.synonyms = T)
    id = NA
    spe_names = NA
    if (y %in% TryAccSpecies$AccSpeciesName){
      spe_names = y
      id = tpl_res$all.entries$id
    }
    if (grepl('sp$', y )){
      spe_names = y
      id = NA
    }
    else {
      if (tpl_res$all.entries$name %in% TryAccSpecies$AccSpeciesName){
        spe_names = tpl_res$all.entries$name
        id = tpl_res$all.entries$id
      }
      else {
        synonyms = tpl_res$synonyms[
          tpl_res$synonyms$confidence.level %in% c('M',"H") & tpl_res$synonyms$name.synonym %in% TryAccSpecies$AccSpeciesName,]$name.synonym
        spe_names = unique(synonyms)
        id = unique(tpl_res$synonyms[tpl_res$synonyms$confidence.level %in% c('M',"H") & tpl_res$synonyms$name.synonym %in% TryAccSpecies$AccSpeciesName,]$id)
      }  
    }
    return(list('spe_names' = as.character(spe_names), 'ids' = id))
  })
  names = unlist(names_ids['spe_names',])
  names = unique(names[!is.na(names) & names != 'NA'])
  names = paste(names, collapse = '___')
  
  ids = names_ids['ids',]
  ids = unique(ids[!is.na(ids) & ids != 'NA'])
  ids = paste(ids, collapse = '___')
  return(list('Synonyms' = names, 'IDs' = ids))}


# Version 2: get ALL synonyms, even if the acepted name is in TRY
get_all_synonyms_and_accepted2 = function(x, match_list){
  x = unlist(strsplit(as.character(x), '___'))
  names_ids = sapply(x, function(y){
    tpl_res = tpl.get(y, return.synonyms = T)
    spe_names = c(y,
                  tpl_res$all.entries$name,
                  tpl_res$synonyms[tpl_res$synonyms$confidence.level %in% c('M',"H") & tpl_res$synonyms$name.synonym %in% TryAccSpecies$AccSpeciesName,]$name.synonym)
    id = c(tpl_res$all.entries$id,
           tpl_res$synonyms[tpl_res$synonyms$confidence.level %in% c('M',"H") & tpl_res$synonyms$name.synonym %in% TryAccSpecies$AccSpeciesName,]$id)
    
    return(list('spe_names' = as.character(spe_names), 'ids' = id))
  })
  names = unlist(names_ids['spe_names',])
  names = unique(names[!is.na(names) & names != 'NA'])
  names = paste(names, collapse = '___')
  
  ids = names_ids['ids',]
  ids = unique(ids[!is.na(ids) & ids != 'NA'])
  ids = paste(ids, collapse = '___')
  return(list('Synonyms' = names, 'IDs' = ids))}

#paste(Species_table$Species_list_corrected, unlist(Species_table$Synonyms))
Species_table$rowID = 1:nrow(Species_table)
Species_table[, c('Synonyms', 'IDs') := sapply(Species_list_corrected, function(x){get_all_synonyms_and_accepted1(x, TryAccSpecies$AccSpeciesName)}), by = rowID]


Species_table_2 = copy(Species_table)

Species_table_2$All_names = unlist(Species_table$Synonyms)

# We save it to avoid re-reconnecting to TPL
#fwrite(Species_table_2, '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/Species_table_temp.csv')
#Species_table_2 = fread('/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/Species_table_temp.csv')

beep()

#### Some species could not be found ####
# Carex_ovalis: unsure which is the right synonym
# Ribes uva-crispa: none of the synonyms can be found in TRY
# Cerastium glutinosum: unsure which synonym
# Alchemilla flavicoma: none of the synonyms can be found in TRY
# Alchemilla gaillardiana: none of the synonyms can be found in TRY

# Need to check:  Taraxacum sect ruderalia

#### Most remaining species are unresolved at the genus or family level ####
# When other species have been identified in the same LU type, we will attribute the average of the traits of those species, weighted by their relative abundance
# When NO other species have been identified in the same LU type, but some have in the other, we will use the values in the other LU
# Otherwise, we'll go through each case independently

# First I add a column, just for these problematic species, to differentiate between forest and grassland estimates
# and I duplicate the corresponding rows
Species_table_2$Forest_or_grassland_estimate = as.character(Species_table_2$Forest_or_grassland_estimate)
Genus_sp = Species_table_2[grepl(" sp$", Species_list_format) & !grepl('Inula', Species_list_format),]
Genus_sp[, Forest_or_grassland_estimate := 'Forest']
Species_table_2[grepl(" sp$", Species_list_format) & !grepl('Inula', Species_list_format), Forest_or_grassland_estimate := 'Grassland']
Species_table_2 = rbind(Species_table_2,
                       Genus_sp)
setDT(Forest_veg_list)

# Then I fill the new column with a character string, containing all species from the same genus with their relative cover
# and a list of species IDs
for (s in sort(unique(Species_table_2[grepl("sp$", Species_list_corrected), Species_list_format]))){
  print(s)
  sp = gsub(' sp', '', s)

  # In grasslands: fetch all rows w the same species, calculate relative abundance
  df_grass = Abundances[grepl(sp, Species) & !grepl("sp$", Species),]
  DF_grass = df_grass[, list(cover = sum(value, na.rm = T)), by = Species][, tot := round(100*cover/sum(cover, na.rm = T))]
  if (DF_grass[, sum(cover)] == 0){DF_grass$tot = 100}
  DF_grass[, Species_format := sapply(Species, function(x){unique(Species_table_2[Species_raw %in% c(x, gsub(' ', '_', x), gsub('_', ' ', x)), Species_list_format])})]
  grass_sp = DF_grass[, Sp_prop := paste(Species_format, '=', tot, '%', sep = ''), by = .I][, paste(Sp_prop, collapse = '___')]

  grass_ids = Species_table_2[Species_list_format %in% DF_grass$Species_format, paste(unique(IDs), collapse = '___')]
  
  # Same in forests
  df_for = Forest_veg_list[grepl(sp, Species) & !grepl("sp$", Species),]
  DF_for = df_for[, list(cover = sum(Cover, na.rm = T)), by = Species][, tot := round(100*cover/sum(cover, na.rm = T))]
  if (DF_for[, sum(cover)] == 0){DF_for$tot = 100}
  DF_for[, Species_format := sapply(Species, function(x){unique(Species_table_2[Species_raw %in% c(x, gsub(' ', '_', x), gsub('_', ' ', x)), Species_list_format])})]
  for_sp = DF_for[, Sp_prop := paste(Species_format, '=', tot, '%', sep = ''), by = .I][, paste(Sp_prop, collapse = '___')]
  for_ids = Species_table_2[Species_list_format %in% DF_for$Species_format, paste(unique(IDs), collapse = '___')]
  
  # If no species found in one of the LU type, replace by the other
  if (grass_sp == '' & for_sp != ''){
    grass_sp = for_sp
    grass_ids = for_ids
    }
  if (grass_sp != '' & for_sp == ''){
    for_sp = grass_sp
    for_ids = grass_ids}
  
  # If no species is found in any of the abundance datasets, we take the species found in the other dataset
  if (grass_sp == '' & for_sp == ''){
    all_sp = Species_table_2[grepl(sp, Species_raw), All_names]
    species = unique(unlist(strsplit(paste(all_sp, collapse = '___'), '___')))
    species = unlist(sapply(species, function(x){unique(Species_table_2[Species_list_format == x, Species_list_format])}))
    species = unique(species[species != ''])
    species = paste(species, '=', round(100 / length(species)), '%', sep = '', collapse = '___')
    
    all_ids = Species_table_2[grepl(sp, Species_raw), IDs]
    ids = unique(unlist(strsplit(paste(all_ids, collapse = '___'), '___')))
    ids = ids[ids != "logical(0)" & ids != "character(0)" ]
    ids = paste(ids, collapse = '___')
    
    grass_ids = ids
    for_ids = ids
    grass_sp = species
    for_sp = species
  }
  
  Species_table_2[Forest_or_grassland_estimate == 'Forest' & Species_list_format == s, c('All_names', 'IDs') := list(for_sp, for_ids) ]
  Species_table_2[Forest_or_grassland_estimate == 'Grassland' & Species_list_format == s,  c('All_names', 'IDs') := list(grass_sp, grass_ids)]
}

Species_table_2 = Species_table_2[!is.na(Synonyms) & Synonyms!= '',]
fwrite(Species_table_2, file = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/Species_table_no_synonyms.csv')
#Species_table_2 = read.csv('/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/Species_table.csv')
beep('mario')


##### Quick check for potential second TRY request
# Species that were already requested
request1 = fread('/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/Try_request_list.csv')
# New species, incl all synonyms
new_w_syn = fread( file = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request_10498/Final_data/Species_table22.csv')
all_sp_w_syn = unique(unlist(strsplit(paste(new_w_syn$All_names, collapse = '___'), '___')))
all_sp_w_syn = gsub('[0-9%=]*', '', all_sp_w_syn)
# New species, not incl synonyms
new_wo_syn = fread( file = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request_10498/Final_data/Species_table.csv')
all_sp_wo_syn = unique(unlist(strsplit(paste(new_wo_syn$All_names, collapse = '___'), '___')))
all_sp_wo_syn = gsub('[0-9%=]*', '', all_sp_wo_syn)

all_sp_w_syn = all_sp_w_syn[(all_sp_w_syn %in% TryAccSpecies$AccSpeciesName)]
all_sp_wo_syn = all_sp_wo_syn[(all_sp_wo_syn %in% TryAccSpecies$AccSpeciesName)]

unique(all_sp_w_syn[!(all_sp_w_syn %in% request1$AccSpeciesName)])
unique(all_sp_wo_syn[!(all_sp_wo_syn %in% request1$AccSpeciesName)])x

unique(all_sp_w_syn[!(all_sp_w_syn %in% all_sp_wo_syn)])
length(unique(all_sp_w_syn))
length(unique(all_sp_wo_syn))
length(unique(request1$AccSpeciesName))

all_sp = unique(c(all_sp_w_syn, all_sp_wo_syn, request1$AccSpeciesName ))

write.csv(TryAccSpecies[AccSpeciesName %in% all_sp, ], '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Species_lists/Try_request_list_SECOND_REQUESR.csv')
paste(sort(unique(TryAccSpecies[AccSpeciesName %in% all_sp, AccSpeciesID])), collapse = ', ')
