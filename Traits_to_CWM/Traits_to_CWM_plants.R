# This script takes as input the abundances and species-level traits of PLANT species
# and outputs a CWM matrix averaged for all considered years.

library(data.table)
library(reshape2)
library(vegan)
library(taxize)
library(FD)
library(ggpmisc)

setwd("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis")

# This is the plant abundance dataset
#Abundances = fread("Data/Raw_data/Abundances/Dataset_clean.txt")
Abundances = Abundances[Group_broad == "Plant",]      

# This is the matching table between Exploratories and TRY data
Species_table <- 
  data.table(read.csv("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/Species_table22.csv"))
Species_table = Species_table[, c('Species_raw', "Species_list_format", "Species_list_corrected", "Synonyms", "All_names", 'Forest_or_grassland_estimate' )]

# This is the AG trait data
AG_traits <- fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/TRY_plant_traits.csv", header = TRUE)

AG_traits = AG_traits[TraitName != '', c('TraitName', 'new.AccSpeciesName', 'UnitName', 'StdValue', 'Std_std', "Reps", "Min", "Max", "Range", "N_contrib", "References")]

# This is the BG trait data
BG_traits_raw = fread("Data/Raw_data/Plants/Below-ground_From_Joana/22326.txt")
BG_traits = melt(BG_traits_raw, id.vars = 'scientificName', variable.name = 'TraitName', value.name = 'StdValue')
BG_traits[, new.AccSpeciesName := gsub('[_]', ' ', scientificName)]
BG_traits[, new.AccSpeciesName := gsub('Lychnis flos-cuculi', 'Silene flos-cuculi', new.AccSpeciesName)]
BG_traits[, new.AccSpeciesName := gsub('Taraxacum sect. Ruderalia', 'Taraxacum sect ruderalia', new.AccSpeciesName)]


All_traits = rbind(AG_traits, BG_traits, fill=TRUE)
All_traits[, TraitName := gsub(' ', '_', TraitName)]
All_traits[, TraitName := gsub('[,/()]', '', TraitName)]

search_species_table <- function(species_raw, trait_name, all_names, Trait_Table) {
  All_names <- unlist(strsplit(as.character(all_names), "___"))
 # print(trait_name)
#  print(All_names)
  # This function takes as input a "raw" species name, as found in the BE, and associates to it the different species (e.g. for aggregate species)
  Trait_data <- Trait_Table[new.AccSpeciesName %in% All_names & TraitName == as.character(trait_name), ]
  if (nrow(Trait_data) == 0) {
    trait_value <- NA
    trait_value2 <- NA
    unit_value <- NA
    sd_value_w <- NA
    sd_type <- NA
    trait_min <- NA
    trait_max <- NA
    trait_range <- NA
    N_contrib = NA
    Reps = NA
    References = NA
  }
  else {
    References = ifelse(length(Trait_data$References) == 1, Trait_data$References,
                        paste(unique(unlist(strsplit(Trait_data$References, '___', useBytes = TRUE))), collapse = '___'))
    
    N_contrib = sum(as.numeric(Trait_data$N_contrib))
    Reps = sum(as.numeric(Trait_data$Reps))
    if (!is.na(as.numeric(Trait_data$StdValue))) { # If trait quantitative
      #print(Trait_data)
      trait_value <- as.character(mean(as.numeric(Trait_data$StdValue), na.rm = T))
      trait_value2 <- as.character(weighted.mean(as.numeric(Trait_data$StdValue), as.numeric(Trait_data$N_contrib), na.rm = T))
      
      if (is.na(trait_value2)){
        trait_value2 = trait_value
      }
      
      trait_min = min(as.numeric(Trait_data$Min), na.rm = T)
      trait_max = max(as.numeric(Trait_data$Max), na.rm = T)
      trait_range = paste(round(trait_min, 2), round(trait_max, 2), collapse = '___')
      
      if (length(Trait_data$StdValue) > 1){
      traits_not_na =   Trait_data[!is.na(Std_std) & Std_std != 'NA' ,]
        
      # Here I calculate the weighted sd as sqrt(SOMME[a * var(A)]), with a the relative weight equal to the number of contribution divided by total number
      sd = as.numeric(traits_not_na$Std_std)
      n = as.numeric(traits_not_na$N_contrib)
      means = as.numeric(Trait_data$StdValue)
      sd_value_w <- as.character(sqrt(sum(n*(sd^2 + (means-as.numeric(trait_value2))^2))/sum(n)))
      sd_type = "pooled variance"
      }
      else {
        sd_type = "single species variance"
        sd_value_w = Trait_data$Std_std
        trait_min <- Trait_data$Min
        trait_max <- Trait_data$Max
        trait_range <- Trait_data$Range
      }}
      
    else { # trait qualitative
      trait_value <- paste(unique(Trait_data$StdValue[!is.na(Trait_data$StdValue) & Trait_data$StdValue != 'NA']), collapse = "___")
      trait_value2 = trait_value
      sd_value_w = NA
      trait_min <- NA
      trait_max <- NA
      trait_range <- NA
      
      }
    unit_value <- unique(Trait_data$UnitName)
    unit_value = unit_value[!is.na(unit_value)]
    unit_value = unit_value[! unit_value %in% c('NA', '')]
    unit_value <- paste(unique(unit_value), collapse = "___")
  }
  res = c(StdValue = trait_value, StdValue2 = trait_value2, UnitName = unit_value, Std_std = sd_value_w,
             Min = trait_min, Max = trait_max, Range = trait_range, 
          References = References, N_contrib = N_contrib, Reps = Reps)
 # print(res)
  return(as.list(res))
}

# Add columns for traits
trait_names <- unique(All_traits$TraitName)
Traits_dummy <- data.table(matrix(ncol = length(trait_names), nrow = nrow(Species_table)))
Traits_dummy[, (colnames(Traits_dummy)) := lapply(.SD, as.character), .SDcols = colnames(Traits_dummy)]
names(Traits_dummy) <- as.character(trait_names)
Species_table2 <- cbind(Species_table[, .SD, .SD = colnames(Species_table)[!(colnames(Species_table) %in% c('IDs', 'rowID'))]], Traits_dummy)

# Melt the datatable
Species_table_melt <- melt.data.table(Species_table2, id.vars = c(
  "Species_raw", "Species_list_format", "Species_list_corrected", "Forest_or_grassland_estimate",
  "Synonyms",  "All_names"
), value.name = "StdValue", variable.name = "TraitName")

Species_table_melt$StdValue <- as.character(Species_table_melt$StdValue)
Species_table_melt[, UnitName := character()]
# Remove duplicates if there is any
Species_table_melt <- unique(Species_table_melt)


# Here we fill each trait value with the average of the synonyms or corresponding species
Species_table_melt[, c("StdValue",'StdValue2', "UnitName", 'Std_std', 'Min', 'Max', 'Range', 'References', 'N_contrib', 'Reps') := search_species_table(Species_raw, TraitName, All_names, All_traits), by = c("Species_raw", "TraitName")]

Species_table_melt2 = copy(Species_table_melt)

# For species only identified at the genus level, we have to take into account the relative proportion of each species.
for (s in unique(Species_table_melt2[grepl("%", All_names), All_names])) {
  list_spe <- data.table(Raw = unlist(strsplit(s, "___")))
  list_spe$Species_list_format <- sapply(list_spe$Raw, function(x) {
    unlist(strsplit(x, "="))[1]
  })
  list_spe$Weight <- as.numeric(gsub("%", "", sapply(list_spe$Raw, function(x) {
    unlist(strsplit(x, "="))[2]
  })))
  list_spe = list_spe[Species_list_format != 'NA' & Weight >0,]
  
  for (trait in trait_names) {

    DF <- Species_table_melt2[Species_list_format %in% list_spe$Species_list_format & TraitName == trait, ]
    DF_with_weight <- merge(list_spe[, c("Species_list_format", "Weight")], DF)
    if (trait == 'ALLO'){
      print(DF_with_weight)
    }
   #if(grepl('Acer', s)){ print(DF_with_weight)}
    if (sum(!is.na(as.numeric(DF_with_weight$StdValue))) > 0) { # Check if at least one numeric value
      trait_value <- weighted.mean(as.numeric(DF_with_weight$StdValue2), DF_with_weight$Weight, na.rm = T)
      
      if (trait == 'ALLO'){
        print(trait_value)
      }
      # Calculation of weighted variance based on https://stats.stackexchange.com/questions/51442/weighted-variance-one-more-time
      v1 = sum(DF_with_weight$Weight)
      v2 =  sum(DF_with_weight$Weight^2)
      trait_var = v1/(v1^2-v2)*sum(DF_with_weight$Weight*(as.numeric(DF_with_weight$StdValue2) - trait_value)^2)
      trait_sd = sqrt(trait_var)
      
      trait_min = min(as.numeric(DF_with_weight$Min), na.rm = T)
      trait_max = min(as.numeric(DF_with_weight$Max), na.rm = T)
      trait_range = paste(round(trait_min, 2), round(trait_max, 2), sep = '___')
      N_contrib = sum(as.numeric(DF_with_weight$N_contrib))
      Reps = sum(as.numeric(DF_with_weight$Reps))
      References = ifelse(length(DF_with_weight$References) == 1, DF_with_weight$References,
                          paste(unique(unlist(strsplit(paste(DF_with_weight$References[!is.na(DF_with_weight$References)], collapse = '___'), '___', useBytes = TRUE))), collapse = '___'))
      }
    else {
      values = DF_with_weight$StdValue
      values = values[!is.na(values) & values != 'NA' & values != '']
      trait_value <- paste(unique(values), collapse = "___")
    }
    unit = unique(DF_with_weight$UnitName)
    if(length(unit) >1){unit = unit[!is.na(unit)]}
    Species_table_melt2[All_names == s & TraitName == trait, 
                       c("StdValue","UnitName", "StdValue2", "Std_std",  "Min",   "Max",   "Range") :=
                    list(trait_value, unit, trait_value, trait_sd,   trait_min, trait_max, trait_range)]
    Species_table_melt2[All_names == s & TraitName == trait,]$References = References
    Species_table_melt2[All_names == s & TraitName == trait,]$N_contrib = N_contrib
    Species_table_melt2[All_names == s & TraitName == trait,]$Reps = Reps

    
    }
}

Species_table_melt2[TraitName == 'Mycorrhiza_type', StdValue2 := ifelse(grepl('Mixed', StdValue), 'Mixed', StdValue)]
Species_table_melt2[TraitName == 'Mycorrhiza_type', StdValue2 := gsub('[___]', '', StdValue)]
Species_table_melt2[TraitName == 'Mycorrhiza_type', StdValue2 := gsub('NA', '', StdValue)]
Species_table_melt2[TraitName == 'Mycorrhiza_type' & (is.na(StdValue) | StdValue == ''), StdValue2 := NA]

Species_table_melt2[TraitName == 'Mycorrhiza_type' & StdValue%in% c("-Non mycorrhizal-Arbuscular mycorrhiza", 
                                                                   "-Arbuscular mycorrhiza-Non mycorrhizal"), StdValue2 := 'Mixed']
Species_table_melt2[StdValue %in% c('', 'NA', 'NaN', "NA___NaN"), StdValue := NA]

fwrite(Species_table_melt2, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/Plant_traits_final_details.csv")

# Only the pretty one (for grasslands)
Species_trait_data_clean = Species_table_melt2[Forest_or_grassland_estimate != 'Forest', list(scientificNameStd = Species_list_format,
                                                                                              traitNameStd = TraitName,
                                                                                              traitValueStd = StdValue2,
                                                                                              traitUnitStd = UnitName)]
fwrite(Species_trait_data_clean, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/Plant_traits_final.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
###### Match with Abundance data #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### a. Homogenise species name ####

# Abundance dataset
Abundances$Species = gsub('_', ' ', Abundances$Species)
Abundances$Species = gsub('Pinus cf mugo', 'Pinus mugo', Abundances$Species)
Abundances$Species = gsub('Euphrasia sp cf', 'Euphrasia sp', Abundances$Species)
Abundances$Species = gsub('Trifolium cf montanum', 'Trifolium montanum', Abundances$Species)
Abundances$Species = gsub('Erigernon', 'Erigeron', Abundances$Species)
Abundances$Species = gsub('Rhinanthus aggr.', 'Rhinanthus', Abundances$Species)
#Abundances$Species = gsub('Trifolium montanum', 'Trifolium cf montanum', Abundances$Species)
Abundances$Species = gsub('Bromus hordeaceus aggr.incl B commutatus', 'Bromus hordeaceus aggr incl B commutatus', Abundances$Species)

# Check if all species are in AG dataset
unique(Abundances[!(Species %in% Species_table_melt2$Species_list_format), Species])
# -> ok

#### AG traits ####
Plant_traits = unique(Species_table_melt2[Forest_or_grassland_estimate != 'Forest',-1])
#Plant_traits[TraitName == "Plant_height_vegetative", c('TraitName', 'StdValue') := list('log_Plant_height', log(as.numeric(StdValue)))]
#Plant_traits[TraitName == "Seed -=dry mass", c('TraitName', 'StdValue') := list('log Seed mass', log(as.numeric(StdValue)))]
Plant_traits_log = Plant_traits[TraitName %in% c("Plant_height_vegetative","Seed_dry_mass", "Plant height vegetative","Seed dry mass"),]
Plant_traits_log[TraitName %in% c("Plant_height_vegetative","Plant height vegetative"), c('TraitName', 'StdValue2') := list('Log_height_vegetative', log(as.numeric(StdValue2)))]
Plant_traits_log[TraitName %in% c("Seed_dry_mass", "Seed dry mass"), c('TraitName', 'StdValue2') := list('log_seed_mass', log(as.numeric(StdValue2)))]

Plant_traits = rbind(Plant_traits,Plant_traits_log)

Plant_traits[, TraitName := gsub(' ', '_', TraitName)]
Trait_names = as.character(unique(Plant_traits$TraitName))

ggplot(Plant_traits[TraitName != "Mycorrhiza type",list( StdValue = as.vector(scale(as.numeric(StdValue2)))), by =  TraitName], 
       aes(StdValue)) + geom_histogram() + facet_wrap(~ TraitName)
Plant_traits_cast = dcast.data.table(Plant_traits, Species_list_format ~ TraitName, value.var = 'StdValue2')
Plant_traits_cast[, (Trait_names[Trait_names != "Mycorrhiza_type"]) := lapply(.SD, as.numeric), .SDcols = Trait_names[Trait_names != "Mycorrhiza_type"]]

# Coverage high to very high for most plots/traits
CC = check_coverage(Plant_traits_cast, Abundances, Trait_names, trait_taxo = 'Species_list_format', 'Species')
CC[, lapply(.SD, function(x){list(min = min(x), median = mean(x), max = max(x))}), .SDcols = Trait_names]
Plants_CWM = my_cwm(Plant_traits_cast, Abundances, Trait_names, trait_taxo = 'Species_list_format', 'Species')

write.csv(Plants_CWM, "Data/CWM_data/CWM_Plants.csv")

ggplot(Plant_traits_cast, aes(x = Mycorrhiza_type, y = Mycorrhizal_infection_intensity)) + geom_violin() #+ 
  #geom_jitter(width = 1)


#### Sensitivity analyses #####
library(ggpmisc)

# Raw trait data
AG_traits <- fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/TRY_plant_traits.csv", header = TRUE)
AG_traits_shade =   fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/TRY_plant_traits_with_shade.csv")
AG_traits_geographic = fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/TRY_plant_traits_onlyCEurope.csv")

All_traits_values = merge.data.table(unique(AG_traits[, list(TraitName, new.AccSpeciesName, value_baseline = StdValue)]),
                                     merge.data.table(unique(AG_traits_shade[, list(TraitName, new.AccSpeciesName, value_growing = StdValue)]),
                                                       unique(AG_traits_geographic[, list(TraitName, new.AccSpeciesName, value_geogr = StdValue)]), by = c('TraitName', 'new.AccSpeciesName')),
                                     by = c('TraitName', 'new.AccSpeciesName'))
All_traits_values[, c('value_baseline', 'value_growing', 'value_geogr') := lapply(.SD, as.numeric), .SDcols = c('value_baseline', 'value_growing', 'value_geogr')]
All_traits_values[TraitName == 'Seed dry mass', c('TraitName', 'value_baseline', 'value_growing', 'value_geogr') := list('log Seed Dry Mass',log(value_baseline), log(value_growing), log(value_geogr)) ]
All_traits_values[, c('value_baseline', 'value_growing', 'value_geogr') := as.data.table(scale(.SD)), .SDcols = c('value_baseline', 'value_growing', 'value_geogr'), by = TraitName]
All_traits_values = All_traits_values[TraitName != 'Mycorrhiza type']
All_traits_values[, TraitName := gsub(' ', '_', TraitName)]
All_traits_values[, TraitName := gsub('[,)(/]', '', TraitName)]

All_traits_values[, Trait_simple := dplyr::recode(gsub(' ', '_', TraitName), 'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included'  = 'LMA_incl_pet',                    
                                        'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded'  = 'LMA',
                                        'Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC'        = 'LDMC',                           
                                        'Leaf_nitrogen_N_content_per_leaf_dry_mass'      = 'Ncontent',                                                     
                                        'Leaf_phosphorus_P_content_per_leaf_dry_mass'       = 'Pcontent',                                                  
                                        'Mycorrhizal_infection_intensity'   = 'Myco_infection_intensity',                                                                  
                                        'Plant_height_vegetative'       = 'Height',                                                                      
                                        'Seed_dry_mass'            = 'Seed mass',                                                                           
                                        'Stem_specific_density_SSD_or_wood_density_stem_dry_mass_per_stem_fresh_volume'  = 'SSD')]


baseline_growing_persp =ggplot(All_traits_values, aes(x= value_baseline, y = value_growing)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y~x) +
  facet_wrap(~Trait_simple) +  geom_abline(slope = 1, intercept = 0, lwd = 0.5, lty = 8) + theme_bw() +
  stat_poly_eq(formula = x~y, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, coef.digits = 2) + xlab('Species trait value based on main dataset') + ylab('Species trait value based on dataset incl. all growing conditions')
ggsave(baseline_growing_persp, file =  '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Metadata/baseline_growing_persp.pdf', width = 9, height = 9)

baseline_geogr_persp = ggplot(All_traits_values, aes(x= value_baseline, y = value_geogr)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y~x) +
  facet_wrap(~Trait_simple) +  geom_abline(slope = 1, intercept = 0, lwd = 0.5, lty = 8) + theme_bw() +
  stat_poly_eq(formula = x~y, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, coef.digits = 2)  + xlab('Species trait value based on main dataset') + ylab('Species trait value based on dataset restricted to C Europe')

ggsave(baseline_geogr_persp, file =  '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Metadata/baseline_geogr_persp.pdf', width = 9, height = 9)


# CWM data
Plants_CWM = fread("Data/CWM_data/CWM_Plants.csv")
Plants_CWM_shade = fread("Data/CWM_data/CWM_Plants_shade.csv")
Plants_CWM_geographic = fread("Data/CWM_data/CWM_Plants_geographic.csv")
Plants_CWM_nosynonyms = fread("Data/CWM_data/CWM_Plants_nosynonyms.csv")

Plants_CWM_melt = melt(Plants_CWM[, -1], id.vars = c('Plot', 'Year'), variable.name = 'Trait', value.name = 'CWM_baseline')
Plants_CWM_shade_melt = melt(Plants_CWM_shade[, -1], id.vars = c('Plot', 'Year'), variable.name = 'Trait', value.name = 'CWM_growing')
Plants_CWM_geographic_melt = melt(Plants_CWM_geographic[, -1], id.vars = c('Plot', 'Year'), variable.name = 'Trait', value.name = 'CWM_geogr')
Plants_CWM_nosynonyms = melt(Plants_CWM_nosynonyms[, -1], id.vars = c('Plot', 'Year'), variable.name = 'Trait', value.name = 'CWM_synonyms')

All_cwm = merge.data.table(merge.data.table(Plants_CWM_melt, Plants_CWM_nosynonyms, by = c('Plot', 'Year', 'Trait')),
                           merge.data.table(Plants_CWM_shade_melt, Plants_CWM_geographic_melt, by = c('Plot', 'Year', 'Trait')),
                           by = c('Plot', 'Year', 'Trait'))

All_cwm = All_cwm[!grepl('Mycorrhiza_type', Trait),]
All_cwm[, c('CWM_baseline', 'CWM_growing', 'CWM_geogr', 'CWM_synonyms') :=  as.data.table(scale(.SD)), .SDcols = c('CWM_baseline', 'CWM_growing', 'CWM_geogr', 'CWM_synonyms'), by = 'Trait']
All_cwm[, Trait_simple := dplyr::recode(Trait, 'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included'  = 'LMA_incl_pet',                    
                                 'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded'  = 'LMA',
                                 'Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC'        = 'LDMC',                           
                                 'Leaf_nitrogen_N_content_per_leaf_dry_mass'      = 'Ncontent',                                                     
                                 'Leaf_phosphorus_P_content_per_leaf_dry_mass'       = 'Pcontent',                                                  
                                 'Mycorrhizal_infection_intensity'   = 'Myco_infection_intensity',                                                                  
                                 'Plant_height_vegetative'       = 'Height',                                                                      
                                 'Seed_dry_mass'            = 'Seed mass',                                                                           
                                 'Stem_specific_density_SSD_or_wood_density_stem_dry_mass_per_stem_fresh_volume'  = 'SSD')]

baseline_shade = ggplot(All_cwm, aes(x= CWM_baseline, y = CWM_growing)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y~x) + theme_bw() +
  facet_wrap(~Trait_simple) + geom_abline(slope = 1, intercept = 0, lwd = 0.5, lty = 8) + 
  stat_poly_eq(formula = x~y, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, coef.digits = 2) + xlab('CWM based on main dataset') + ylab('CWM based on dataset incl. all growing conditions')
ggsave(baseline_shade, file =  '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Metadata/baseline_shade.pdf', width = 9, height = 9)

baseline_syno = ggplot(All_cwm, aes(x= CWM_baseline, y = CWM_synonyms)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y~x) +theme_bw() +
  facet_wrap(~Trait_simple) + geom_abline(slope = 1, intercept = 0, lwd = 0.5, lty = 8) + 
  stat_poly_eq(formula = x~y, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, coef.digits = 2)  + xlab('CWM based on main dataset') + ylab('CWM based on dataset without synonym species')
ggsave(baseline_syno, file =  '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Metadata/baseline_syno.pdf', width = 9, height = 9)

baseline_geogr = ggplot(All_cwm, aes(x= CWM_baseline, y = CWM_geogr)) + geom_point() + geom_smooth(method = 'lm') +
  facet_wrap(~Trait_simple) + geom_abline(slope = 1, intercept = 0, lwd = 0.5, lty = 8) + theme_bw() +
  stat_poly_eq(formula = x~y, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, coef.digits = 2)  + xlab('CWM based on main dataset') + ylab('CWM based on dataset restricted to Central Europe')

ggsave(baseline_geogr, file =  '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Metadata/baseline_geogr.pdf', width = 9, height = 9)

