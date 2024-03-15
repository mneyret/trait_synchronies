# This is the script to merge multitrophic CWM trait data used in Neyret, M., Le Provost, G., Boesing, A.L. et al. A slow-fast trait continuum at the whole community level in relation to land-use intensification. Nat Commun 15, 1251 (2024).
# It calls separate scripts for each functional guild. Please check the 1_Get_data file to get the raw data.

library(ade4) # v 1.7-22
library(betapart) # v 1.6
library(car) # v 3.1-1
library(data.table) # v 1.14.8
library(dplyr) # v 1.1.0
library(factoextra) # v 1.0.7
library(FD) # v 1.0-12.1
library(ggcorrplot) # v 0.1.4
library(ggpmisc) # v 0.5.2
library(Hmisc) # v 5.0-1
library(mice) # v 3.15.0 
library(plyr) # v 1.1.0
library(readr) # v 2.1.4
library(readxl) # v 1.4.2  
library(reshape2) # v 1.4.4 
library(stringr) # v 1.5.0
library(taxize) # v 0.9.100  
library(Taxonstand) # v 2.4
library(textclean) # v 0.9.3
library(traitdataform) # v 0.6.8 # https://github.com/EcologicalTraitData/traitdataform
library(vegan) # v 2.6-4 

wdir = "trait_synchronies/"
setwd(wdir)

# Get the main abundance dataset (most CWM will be based on abundances from this dataset)
## Raw diversity

allsp <- fread("Data/Abundance_data/27707_2_Dataset/27707_2_data.csv")
allsp$Species = gsub('_$', '', allsp$Species ) # Remove _ if last character
## Species information
fgs <- fread("Data/Abundance_data/27706_3_Dataset/27706_3_data.csv")
fgs$Species = gsub(' $', '', fgs$Species )
fgs$Species = gsub(' ', '_', fgs$Species )
Abundance_all <- merge.data.table(allsp, fgs, by ="Species", all.x=TRUE)
Abundance_all[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1, 3), '0', substr(Plot, 4, 4), sep = ''))]

### Run CWM analyses for all individual groups ####
source('Functions.R')
rerun_CWM = FALSE
if (rerun_CWM == TRUE){
#source('Traits_to_CWM/Traits_to_CWM_Plants.R')
source('Traits_to_CWM/Traits_to_CWM_Birds.R')
source('Traits_to_CWM/Traits_to_CWM_Protists.R') 
  
source('Traits_to_CWM/Traits_to_CWM_Bacterias.R') 
  
source('Traits_to_CWM/Traits_to_CWM_Bats.R') 
source('Traits_to_CWM/Traits_to_CWM_Protists.R') 
source('Traits_to_CWM/Traits_to_CWM_Arthropods.R') 
#source('Traits_to_CWM/Traits_to_CWM_Butterflies_Moths.R')
}

### Reimport and merge all datasets ####
for (weighted in c(TRUE#, FALSE
                   )){
  all_data = data.table()
for (group in c('birds', 'bats', 'plants', 'arthropods_below_omni_carni', 'arthropods_below_herb', 'arthropods_above_carni', 'arthropods_above_herb', 
                'butterflies',  'coll', 'mites', 'microbes', 'protists_sec_cons', 'protists_bact','protists')){

  if (weighted == TRUE){  data = fread(paste('Data/CWM_data/CWM_', group, '.csv', sep = ''))}
  if (weighted == FALSE){  data = fread(paste('Data/CWM_data/CWM_', group, '_noweight.csv', sep = ''))}
  
  data[, Group := group]
  data[, traitName2 := paste(group, traitName, sep ="_")]
  all_data = rbind(all_data, data)
}


### Add additional, community-level traits for fungi ####
## FB ratio from corresponding datasets
microb_soil_prop_2011 <- fread("Data/Trait_data/20250_3_Dataset/20250_3_data.csv") # https://www.bexis.uni-jena.de/ddm/data/Showdata/20250
microb_soil_prop_2014 <- fread("Data/Trait_data/20251_Dataset/20251.txt") # https://www.bexis.uni-jena.de/ddm/data/Showdata/20251

microb_soil_prop = rbindlist(list(microb_soil_prop_2011[, c('Year', 'EP_Plot_ID', 'fungi_bacteria')], 
                                  microb_soil_prop_2014[, c('Year', 'EP_Plot_ID', 'fungi_bacteria')]))
microb_soil_prop[, Plot := ifelse(nchar(EP_Plot_ID) == 5, EP_Plot_ID, paste(substr(EP_Plot_ID, 1, 3), '0', substr(EP_Plot_ID, 4, 4), sep = ''))]

Abundances_fungi = Abundance_all[Group_broad == "soilfungi",]
Prop_fun_group = Abundances_fungi[, list(value = sum(value)), by = c('Plot', 'Fun_group_broad', 'Year')]
Prop_fun_group = dcast.data.table(Prop_fun_group, Year+ Plot~Fun_group_broad, value.var = 'value')
Prop_fun_group_final = Prop_fun_group[, c('total', 'Year') := list(rowSums(.SD), as.numeric(Year)), .SD = c('soilfungi.decomposer', 'soilfungi.other', 'soilfungi.pathotroph', 'soilfungi.symbiont')]
Prop_fun_group_final[, c('soilfungi.decomposer', 'soilfungi.other', 'soilfungi.pathotroph', 'soilfungi.symbiont') := lapply(.SD, function(x){x/total}), .SD = c('soilfungi.decomposer', 'soilfungi.other', 'soilfungi.pathotroph', 'soilfungi.symbiont')]
Prop_fun_group_final[, traitCoverage := 1 - soilfungi.other]
microb_soil_prop2 = merge(microb_soil_prop[, -'EP_Plot_ID'], Prop_fun_group_final[, list(Year,  Plot, soilfungi.decomposer, soilfungi.pathotroph, soilfungi.symbiont, traitCoverage)], by = c('Plot', 'Year'), all = T)

microb_melt = melt.data.table(microb_soil_prop2, id.vars = c('Plot', 'Year', 'traitCoverage'), variable.name = "traitName", value.var = 'traitValue', value.name = 'traitValue')
microb_melt[traitName == 'fungi_bacteria', c('traitName', 'traitCoverage', 'traitUnit', 'traitDescription', 'traitDataRef', 'TraitDataID', 'AbundanceDataID', 'Group', 'traitName2') :=
                                            list('FB_ratio', NA, 'NA', 'Ratio of fungi:bacteria based on PLFA','NA' , 'Bexis ID 20250, 20251', 'NA', 'microbes', 'microbes_FB_ratio')]
microb_melt[traitName == 'soilfungi.decomposer', c('traitName', 'traitUnit', 'traitDescription', 'traitDataRef', 'TraitDataID', 'AbundanceDataID', 'Group', 'traitName2') :=
              list('Decomposer', '%', 'Proportion of decomposer fungi','NA' , 'NA', 'Bexis ID 26470, 26472', 'microbes', 'microbes_decomposer')]
microb_melt[traitName == 'soilfungi.pathotroph', c('traitName', 'traitUnit', 'traitDescription', 'traitDataRef', 'TraitDataID', 'AbundanceDataID', 'Group', 'traitName2') :=
              list('Pathotrophs', '%', 'Proportion of pathogenic fungi','NA' , 'NA', 'Bexis ID 26470, 26472', 'microbes', 'microbes_pathogens')]
microb_melt[traitName == 'soilfungi.symbiont', c('traitName', 'traitUnit', 'traitDescription', 'traitDataRef', 'TraitDataID', 'AbundanceDataID', 'Group', 'traitName2') :=
              list('Symbionts', '%', 'Proportion of symbiont fungi','NA' , 'NA', 'Bexis ID 26470, 26472', 'microbes', 'microbes_symbiont')]


# Add to main dataset
all_data = rbind(all_data, microb_melt)


# Keep only slow-fast traits
#traits_use=c( # birds
#       "birds_logBody_Mass",        "birds_logIncub_time",       "birds_logLongevity" ,            "birds_logNumber_offspring",    'birds_GenLength',          
# 
#  # Bats
#  "bats_logBody_mass"  ,         "bats_Lifespan"                , 'bats_Number_offspring',             
#  
#  # Plants 
#  "plants_Fine_roots_diameter", "plants_LDMC", "plants_LeafN", "plants_LeafP", "plants_Root_tissue_density", "plants_SLA_all", "plants_Seed_mass",
#             # Most arthropods
#  "arthropods_below_omni_carni_Dispersal_ability", "arthropods_below_omni_carni_logBody_Size", 
#"arthropods_below_herb_Dispersal_ability", "arthropods_below_herb_logBody_Size", 
#"arthropods_below_herb_Feeding_generalism", "arthropods_above_carni_Dispersal_ability", 
#"arthropods_above_carni_logBody_Size", "arthropods_above_herb_Dispersal_ability", 
#"arthropods_above_herb_logBody_Size", "arthropods_above_herb_Feeding_generalism", 
#"arthropods_above_herb_Generations", 
#              # Butterflies
#"butterflies_logFlight", 'butterflies_logSize',
#"butterflies_Wintering_stage", "butterflies_Generalism_use", 
#"butterflies_Voltinism_use", 
#              # Collembola and mites 
#"coll_logSize", "coll_Gen_per_year", "coll_Depth_preference",  "coll_Repro_sex", 
#"mites_logMass", "mites_Feeding_spec", 
#"mites_Habitat_spec", "mites_Repro_sex", 
#"mites_DaystoAdult", 
#              # Bacteria and fungi
# "microbes_logVolume", "microbes_OC_ratio", "microbes_FB_ratio", "microbes_pathogens", 'microbes_Genome_size',
#              # Protists
#"protists_sec_cons_Size", "protists_bact_Size",  "protists_nutrition_code_primary_cons"
#)


# Clean up
all_data[Group == 'microbes', c('Group','traitName2') := list("Bact_fun", gsub('microbes', 'bact_fungi',traitName2))]

all_data[, traitCoverage := as.character(round(traitCoverage, 2))]
all_data[traitName %in% c('FB_ratio', 'Symbionts', 'Decomposer', 'Pathotrophs'), traitCoverage := 'NA']
all_data[traitName %in% c('FB_ratio', 'Symbionts', 'Decomposer', 'Pathotrophs'), traitDescription := paste(traitDescription, '(trait coverage not defined because community-level measure)')]


# Clean up dataset
all_data[, traitValue := as.character(round(traitValue, 6))]
all_data[is.na(traitValue), traitValue := 'NA']
all_data[is.na(TraitDataID) | TraitDataID == '' | TraitDataID == "NA", TraitDataID := 'No trait data from bexis']
all_data[is.na(traitDataRef)| traitDataRef == ""|  traitDataRef == "NA", traitDataRef := 'No external data source']
all_data[, traitDataRef := gsub('\n', '', traitDataRef)]
all_data[traitUnit == '' |traitUnit == 'NA' | traitUnit == 'unitless' , traitUnit := 'no unit']
all_data[traitUnit == 'µm', traitUnit := 'microm'] 
all_data = all_data[grepl('G', Plot),] # Exclude forest plots (protists)
all_data[,EP_PlotID := gsub('G0', 'G', Plot)] # match BExis format

datasets = '20067, 31368, 26587, 27586,  31122,  21228,  24468, 24426,  20250, 20251, 21446, 21447, 21448, 21449, 24690, 25306 , 27707, 19849, 19850 , 27707, 27386 , 27707, 21969, 27406 , 27707, 24468, 24426,  26470, 26472'

# Save

if (weighted == TRUE){      fwrite(all_data, "Data/CWM_data/All_CWM_data.csv", sep = "\t")}

if (weighted == FALSE){     fwrite(all_data, "Data/CWM_data/All_CWM_data_noweight.csv", sep = ";")}

}


