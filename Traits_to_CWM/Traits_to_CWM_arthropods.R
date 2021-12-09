# This script takes as input the abundances and species-level traits of arthropod species
# and outputs a CWM matrix averaged for all considered years.
library(data.table)
library(reshape2)
library(vegan)
library(taxize)
library(FD)
library(readr)
library(textclean)
library(readxl)
library(dplyr)
library(traitdataform)
library(mice)

setwd("/Users/Margot/Desktop/Research/Senckenberg/Data/")

###### Functions ######
# Standardise Name in datatable
get_gbif = function(x){
  if(is.na(x)){
    return(NA)
  } else {
    if (grepl('sp[ec]*[\\.]*$',x)){
      x_genus = gsub('sp[ec]*[\\.]*$','', x)
      x_genus = gsub(' ', '', x_genus)
      genus = try(get_gbif_taxonomy(x_genus)$scientificName)
      if (class(genus) == 'try-error') {return(paste(x_genus, 'sp.'))}
      return(paste(genus, 'sp.'))
    }
    else{
  names = get_gbif_taxonomy(unique(x))
  if(names$warnings =="No matching species concept! "){
    return(names$verbatimScientificName)
  } else{
    return(names$scientificName)
}}}
}


###### Input data ####### 
### Traits
#published in Goessner 2015
Arthropod_traits_Goessner = setDT(read_delim("Traits/Arthropods/ArthropodSpeciesTraits.csv", ";", escape_double = FALSE, trim_ws = TRUE))
Arthropod_traits_Goessner[, gbifName := get_gbif(SpeciesID), by = SpeciesID]

### New version
Arthropod_traits_grasslands_Raw = fread('Traits/Arthropods/31122_6_Dataset/31122_6_data.csv')
Arthropod_traits_grasslands_Raw = unique(Arthropod_traits_grasslands_Raw)
Arthropod_traits_grasslands_Raw[, gbifName := get_gbif(SpeciesID), by = SpeciesID]
Arthropod_traits_grasslands_Raw[, Dispersal_ability := as.numeric(gsub(',', '.', Dispersal_ability))]

### Additional traits
# Generation time for herbivore
Voltinism_traits <- data.table(read_excel("Traits/Arthropods/Voltinism_herbivore.xlsx"))
Voltinism_traits[, gbifName := get_gbif(Species_raw), by = Species_raw]
Voltinism_traits[, Generations := as.numeric(Generations)]
Voltinism_traits = Voltinism_traits[!is.na(gbifName),]
# Traits from Birkhofer et al. 2017 (Araneae)
#Birkhofer_araneae = data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Arthropods/Birkhofer_2017/Traits_Communities.xlsx", sheet = 'Araneae'))
#Birkhofer_araneae[, new_scientificName := get_gbif(gsub('\\.', ' ', Site)), by = Site]
#Birkhofer_araneae$new_scientificName %in% Arthropod_traits_grasslands$new_scientificName

### Abundances
Abundances_all = fread("Abundances/Dataset_clean.txt")
Abundances_all[, length(unique(Plot)),by = c('Group_broad', 'Year') ]
Abundances_2008 = Abundances_all[Group_broad %in% c("Araneae" , "Coleoptera", "Hemiptera" ,"Orthoptera"),] 
Abundances_temporal = fread('Abundances/21969_4_Dataset/21969_4_data.csv')
Abundances_temporal = Abundances_temporal[!is.na(Species),]
Abundances_temporal[, Plot := ifelse(nchar(PlotID) == 5, PlotID, paste(substr(PlotID, 1, 3), '0', substr(PlotID, 4, 4), sep = ''))]
Abundances_temporal_by_plot = Abundances_temporal[, list(value = sum(NumberAdults, na.rm = T)), by = list('Plot' = Plot, 'Species' = Species, 'Year' = CollectionYear, 'Order' = Order)]
Abundances_temporal_by_plot[, gbifName := get_gbif(Species), by = Species]


###### Data cleaning and standardisation ####### 
# One species wasn't found by the function get_gbif
Abundances_temporal_by_plot[Species == 'Pachybrachius fracticollis' , new_scientificName := 'Pachymerus fracticollis']
Arthropod_traits_grasslands[SpeciesID == 'Pachybrachius fracticollis' , new_scientificName := 'Pachymerus fracticollis']


# Merge with other trait dataset
Arthropod_traits_grasslands = merge(Arthropod_traits_grasslands_Raw, Arthropod_traits_Goessner[, list(Feeding_guild_Goessner = Feeding_guild, Feeding_specialization, gbifName)], all.x = T)
Arthropod_traits_grasslands = merge(Arthropod_traits_grasslands, Voltinism_traits[, list(Generations, gbifName)], all.x = T)
Arthropod_traits_grasslands = unique(Arthropod_traits_grasslands)

Arthropod_traits_grasslands[, new_scientificName := ifelse(SpeciesID %in% Abundances_temporal_by_plot$Species, SpeciesID, gbifName), by = SpeciesID]
Abundances_temporal_by_plot[, new_scientificName := ifelse(Species %in% Arthropod_traits_grasslands_Raw$SpeciesID, Species, gbifName), by = Species]


#%%%%%%%%%%%%%%%%%%%%%#
#### Recode traits ####
#%%%%%%%%%%%%%%%%%%%%%#

# Carnivores are secondary consumers and the rest are primary consumers
Arthropod_traits_grasslands[, Feeding_guild_simple := dplyr::recode(Feeding_guild, 
                                                              "h"        = 'primary',
                                                              "f"        = 'primary', 
                                                              "m"        = 'primary',
                                                              'd'        = 'primary',
                                                              "c"       = 'secondary',
                                                              )]
# To classify omnivores, let's look more precisely at the Goessner dataset.
# Species which are carnivores at at least one of their lifestages are classified as secondary consumers, the other as primary consumers
Arthropod_traits_grasslands[Feeding_guild == 'o' & grepl('c', Feeding_guild_Goessner) & Feeding_guild_Goessner != 'h-(c)', Feeding_guild_simple := 'secondary']
Arthropod_traits_grasslands[Feeding_guild == 'o' & Feeding_guild_simple == 'o', Feeding_guild_simple := 'primary']

# Simplyfy stratum_use into below-ground (soil, ground) and above-ground
Arthropod_traits_grasslands[, Stratum_use_simple := dplyr::recode(Stratum_use, 
                                                                    "g" = "below",
                                                                    "h" = "above",
                                                                    "i" = "above",
                                                                    "s" = "below",
                                                                    "t" = "above",
                                                                    "u" = "above",
                                                                    "w" = "NA")
]


# Also using the old dataset
Arthropod_traits_grasslands[, Feeding_generalism := as.numeric(factor(Feeding_specialization, levels = c('m', 'o', 'p')))]


##### Add taxa identified only at genus level ####
# Some in abundances are not ID at species level. Let's see if we can interpolate the species-level data to genus level.
un_id_species = unique(Abundances_temporal_by_plot[!(Species %in% Arthropod_traits_grasslands$SpeciesID), list(Species, new_scientificName)])
un_id_species[, genus := gsub(' sp.', '', new_scientificName)]
un_id_species[, Genus := gsub(' sp.', '', new_scientificName)]
un_id_species[, SpeciesID := Species]

# Order, Family
un_id_species[, c( 'Order', 'Family') := list(unique( unique(Arthropod_traits_grasslands[Genus == genus, Order]),
                                                          unique(Arthropod_traits_grasslands[Genus == genus, Family]))), by = genus]
# Numeric traits
un_id_species[, c('Mean_BodySize', 'Dispersal_ability', 'Generations') := list(as.numeric(ifelse( Arthropod_traits_grasslands[Genus == genus, mean(Mean_BodySize, na.rm = T)] > 2*Arthropod_traits_grasslands[Genus == genus, sd(Mean_BodySize, na.rm = T)],
                                                                                  Arthropod_traits_grasslands[Genus == genus, mean(Mean_BodySize, na.rm = T)],
                                                                                  NA)),
                                                                           as.numeric(ifelse(Arthropod_traits_grasslands[Genus == genus, mean(Dispersal_ability, na.rm = T)] > 2* Arthropod_traits_grasslands[Genus == genus, sd(Dispersal_ability, na.rm = T)],
                                                                                  Arthropod_traits_grasslands[Genus == genus, mean(Dispersal_ability, na.rm = T)],
                                                                                  NA)),
                                                                           as.numeric(ifelse(Arthropod_traits_grasslands[Genus == genus, mean(Generations, na.rm = T)] > 2* Arthropod_traits_grasslands[Genus == genus, sd(Generations, na.rm = T)],
                                                                                  Arthropod_traits_grasslands[Genus == genus, mean(Generations, na.rm = T)],
                                                                                  NA))), by = genus]
# Qualitative traits traits
un_id_species[, Stratum_use_simple := as.character(ifelse(
  Arthropod_traits_grasslands[Genus == genus, max(table(Stratum_use_simple[!is.na(Stratum_use_simple)]), na.rm = T)>0.9*sum(table(Stratum_use_simple[!is.na(Stratum_use_simple)]), na.rm = T)],
  Arthropod_traits_grasslands[Genus == genus, names(sort(table(Stratum_use_simple[!is.na(Stratum_use_simple)]), decreasing = T))[1]],
                                                                 NA)), by = genus]
un_id_species[, Feeding_generalism := as.character(ifelse(
  Arthropod_traits_grasslands[Genus == genus, max(table(Feeding_generalism[!is.na(Feeding_generalism)]), na.rm = T)>0.9*sum(table(Feeding_generalism[!is.na(Feeding_generalism)]), na.rm = T)],
  Arthropod_traits_grasslands[Genus == genus, names(sort(table(Feeding_generalism[!is.na(Feeding_generalism)]), decreasing = T))[1]],
  NA)), by = genus]
un_id_species[, Feeding_guild_simple := as.character(ifelse(
  Arthropod_traits_grasslands[Genus == genus, max(table(Feeding_guild_simple[!is.na(Feeding_guild_simple)]), na.rm = T)>0.9*sum(table(Feeding_guild_simple[!is.na(Feeding_guild_simple)]), na.rm = T)],
  Arthropod_traits_grasslands[Genus == genus, names(sort(table(Feeding_guild_simple[!is.na(Feeding_guild_simple)]), decreasing = T))[1]],
  NA)), by = genus]

un_id_species = un_id_species[new_scientificName != " sp." ]

# Merge with initial dataset
Arthropod_traits = unique(rbind(Arthropod_traits_grasslands[, list(SpeciesID, Order, Family, Genus, Species, Mean_BodySize, Mean_BodySize, Dispersal_ability, Feeding_guild_simple, Stratum_use_simple, Feeding_generalism, Generations, new_scientificName)],
                                              un_id_species[, list(SpeciesID, Order, Family, Genus, Species, Mean_BodySize, Mean_BodySize, Dispersal_ability, Feeding_guild_simple, Stratum_use_simple, Feeding_generalism, Generations, new_scientificName)]))
Arthropod_traits = unique(Arthropod_traits)
Arthropod_traits = Arthropod_traits[!is.na(new_scientificName),]

Arthropod_traits$Feeding_generalism = as.numeric(Arthropod_traits$Feeding_generalism)
###### Check trait distribution
Arthropod_traits[, logBody_Size := log(Mean_BodySize)]


### Calculate CWM
# Some species from abundance are not included in trait so need to adapt the number reported:
# AG, herb: + 1 Auchenorrhyncha spec., Typhlocybinae spec., Stenodemini, Deltocephalinae spec."

# Herbivore, AG
traits_h_AG = Arthropod_traits[Stratum_use_simple == 'above'  & 
                                          Feeding_guild_simple == 'primary' &
                                          new_scientificName %in% Abundances_temporal_by_plot$new_scientificName, .SD, 
                               .SDcols = c('Feeding_generalism', 'Dispersal_ability', 'logBody_Size', 'Generations')]
traits_h_AG[, lapply(.SD, function(x){length(x[!is.na(x)])})]

cor.test(traits_h_AG$logBody_Size, traits_h_AG$Feeding_generalism)
pca_species = dudi.pca(traits_h_AG[!is.na(Generations) & !is.na(Feeding_generalism) & !is.na(Dispersal_ability),], scannf = FALSE, nf = 2)
fviz_pca(pca_species)

Ab_h_AG = Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' , SpeciesID] 
                                      | new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary', new_scientificName]
                                      | Species %in% c("Auchenorrhyncha spec.", "Typhlocybinae spec." ,  "Stenodemini spec.", "Deltocephalinae spec."), 
                                      list(value = sum(value), Year = 2011), by = c('new_scientificName','Species', 'Plot')]
Coverage_Arthropods_above_herb = check_coverage(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' ,], 
                                                Ab_h_AG,
                                                c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Generations'), 'new_scientificName', 'new_scientificName')
Coverage_Arthropods_above_herb[, -c(1,2)][, lapply(.SD, min)]


CWM_Arthropods_above_herb = my_cwm(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' ,], 
                                   Abundances_temporal_by_plot, c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Generations'), 'new_scientificName', 'new_scientificName')
test_pca = dudi.pca(complete(mice(CWM_Arthropods_above_herb[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability','Feeding_generalism', 'logBody_Size',  'Generations')][, -1])), scannf = FALSE, nf = 2)
fviz_pca(test_pca)


# Carnivores, AG ---> Need to take axis 2
trait_c_AG = Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName,lapply(.SD, as.numeric), .SDcols = c('Dispersal_ability', 'logBody_Size')]

trait_c_AG[, lapply(.SD, function(x){length(x[!is.na(x)])})]

pca_c_AG = dudi.pca(trait_c_AG[complete.cases(trait_c_AG),], scannf = FALSE, nf = 2)
fviz_pca(pca_c_AG)
cor.test(trait_c_AG$Dispersal_ability, trait_c_AG$Active_hunt)


Coverage_Arthropods_above_carni = check_coverage(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' ,], 
                                                Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' , SpeciesID] |
                                                                              new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' , new_scientificName], 
                                                                            list(value = sum(value), Year = 2011), by = c('new_scientificName', 'Plot')],
                                                c( 'Dispersal_ability',"logBody_Size"), 'new_scientificName', 'new_scientificName')

CWM_Arthropods_above_omni_carni = my_cwm(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple  == 'secondary',], 
                                         Abundances_temporal_by_plot, c('Dispersal_ability',"logBody_Size"), 'new_scientificName', 'new_scientificName')


pca_c_AG = dudi.pca(complete(mice(CWM_Arthropods_above_omni_carni[, lapply(.SD, mean, na.rm = T), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size')][, -1])), , scannf = FALSE, nf = 2)
fviz_pca(pca_c_AG)

# Herbivores, BG
trait_h_BG = Arthropod_traits[Stratum_use_simple == 'below' & Feeding_guild_simple  == 'primary'  &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName,lapply(.SD, as.numeric), 
                              .SDcols = c('Dispersal_ability', 'logBody_Size','Feeding_generalism')]
trait_h_BG[, lapply(.SD, function(x){length(x[!is.na(x)])})]


pca_h_BG = dudi.pca(trait_h_BG[complete.cases(trait_h_BG),], scannf = FALSE, nf = 2)
fviz_pca(pca_h_BG)
cor.test(trait_h_BG$logBody_Size, trait_h_BG$Dispersal_ability)

Coverage_Arthropods_below_herb = check_coverage(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' ,], 
                                                Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' , SpeciesID] |
                                                                              new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' , new_scientificName], 
                                                                            list(value = sum(value), Year = 2011), by = c('new_scientificName', 'Plot')],
                                                c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism'), 'new_scientificName', 'new_scientificName')
CWM_Arthropods_below_herb = my_cwm(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' ,], 
                                   Abundances_temporal_by_plot, c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Generations'), 'new_scientificName', 'new_scientificName')


# Carnivores, BG
trait_c_BG = Arthropod_traits[Stratum_use_simple == 'below' & Feeding_guild_simple == 'secondary' &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName,lapply(.SD, as.numeric), .SDcols = c('Dispersal_ability', 'logBody_Size','Feeding_generalism')]

trait_c_BG[, lapply(.SD, function(x){length(x[!is.na(x)])})]

pca_c_BG = dudi.pca(trait_c_BG[complete.cases(trait_c_BG),], scannf = FALSE, nf = 2)
fviz_pca(pca_c_BG)
cor.test(trait_c_BG$logBody_Size, trait_c_BG$Dispersal_ability)

Coverage_Arthropods_below_carni = check_coverage(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' ,], 
                                                Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' , SpeciesID] |
                                                                              new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' , new_scientificName], 
                                                                            list(value = sum(value), Year = 2011), by = c('new_scientificName', 'Plot')],
                                                c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism'), 'new_scientificName', 'new_scientificName')
CWM_Arthropods_below_omni_carni  = my_cwm(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' ,], 
                                   Abundances_temporal_by_plot, c( 'Dispersal_ability',"logBody_Size"), 'SpeciesID', 'Species')

test_pca = dudi.pca(mice::complete(mice(CWM_Arthropods_below_omni_carni[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size')][, -1])), scannf = FALSE, nf = 2)
fviz_pca(test_pca)

# Save datasets
write.csv(CWM_Arthropods_above_herb[,list( "Ah_Dispersal" = Dispersal_ability, 
                                         "Ah_BodySize" = logBody_Size,
                                         "Ah_Generalism" = Feeding_generalism,
                                         "Ah_Generations" = Generations,
                                         "Plot" = Plot,
                                         "Year" = Year )], 

          "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_above_herb.csv")
 
write.csv(CWM_Arthropods_below_herb[,list( "Ah_b_Dispersal" = Dispersal_ability, 
                                           "Ah_b_BodySize" = logBody_Size,
                                           "Ah_b_Generalism" = Feeding_generalism,
                                           "Plot" = Plot,
                                           "Year" = Year )], 
          
          "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_below_herb.csv")


write.csv(CWM_Arthropods_above_omni_carni[,list( "Aoc_Dispersal" = Dispersal_ability, 
                                           "Aoc_BodySize" = logBody_Size,
                                           "Plot" = Plot,
                                           "Year" = Year )], 
          
          "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_above_omni_carni.csv")



write.csv(CWM_Arthropods_below_omni_carni[,list(
                                            "Aoc_b_Dispersal" = Dispersal_ability, 
                                            "Aoc_b_BodySize" = logBody_Size,
                                            "Plot" = Plot,
                                            "Year" = Year )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Arthropods_below_omni_carni.csv")

######  !!! MORRRRRE TRRAAAAAITS #####
### Collembola
coll_traits = data.table(read_excel(('Traits/Collembola_Mites/Coll_Traits_NEW.xlsx')))
coll_traits[, c( 'Size', 'Size_Adult', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales') := lapply(.SD, as.numeric),
            .SDcols = c( 'Size','Size_Adult',  'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales')]

# Check distribution
coll_traits[, lapply(.SD, function(x){hist(log(x))}), .SDcols = c( 'Size', 'Size_Adult', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales')]
coll_traits[, c( "Size.log"):= lapply(.SD, log), .SDcols = c( "Size_Adult")]
coll_traits[, Depth_preference := ifelse(Vertical_preference == 'Euedaphic', 2,
                                  ifelse( Vertical_preference == 'Epiedaphic',0,
                                  1))]
coll_traits[, Repro_sex := as.numeric(Reproduction_Ruslan == 'bisex')]
coll_traits[, Gen_per_year := ifelse(Phenology == 'multivoltine', 3, 
                                     ifelse(Phenology == 'bivoltine', 2, 
                                            NA))]
coll_traits[, Species := gsub(' ', '_', Species)]
coll_traits[, Species := gsub('__', '_', Species)]


pca_species = dudi.pca(coll_traits[complete.cases(coll_traits[,c( 'Size.log', 'Gen_per_year', 'Depth_preference', 'Repro_sex')]),c( 'Size.log', 'Gen_per_year', 'Depth_preference', 'Repro_sex')], scannf = FALSE, nf = 3)
pca_coll_species= fviz_pca_biplot(pca_species, geom = c("point"), repel = T, axes = c(1,2), title = 'Collembola')
ggsave(pca_coll_species,file= '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Species_PCA_Collembola.pdf', width = 5, height = 5)    



# Plot AEG43: only one individual, undefined, was found so we remove it
Abundances_coll = Abundances_all[Group_broad =='Collembola' & Plot != 'AEG43',]

# Check number of species for each trait
coll_traits[Species %in% Abundances_coll$Species , lapply(.SD, function(x){length(x[!is.na(x) & x != 'NA'])})]


Coll_coverage = check_coverage(coll_traits, Abundances_coll, c( 'Size.log', 'Gen_per_year', 'Depth_preference', 'Repro_sex', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'#, 'Vertical_preference'
), 'Species', 'Species')
Coll_cwm_all = my_cwm(coll_traits, Abundances_coll, c( 'Size.log', 'Gen_per_year', 'Depth_preference', 'Repro_sex', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'#, 'Vertical_preference'
                                                       ), 'Species', 'Species')

Coll_cwm_all[is.nan(Gen_per_year), Gen_per_year := mean(Gen_per_year)]

pca_coll = dudi.pca(complete(mice(Coll_cwm_all[ ,c( 'Size.log', 'Gen_per_year', 'Depth_preference', 'Repro_sex_1')])), scannf = FALSE, nf = 3)
fviz_pca_biplot(pca_coll, geom = c("point"), repel = T, axes = c(1,2), title = 'Collembola')

cor.test(pca_coll$l1$RS1, env_data_lui[Plot %in% Coll_cwm_all$Plot,]$LUI)


write.csv(Coll_cwm_all[, list( "col_Size" = Size.log,
                               "col_Depth" = Depth_preference,
                               "col_Gen_per_Year" = Gen_per_year,
                               "col_Sex" = Repro_sex_1,
                               "Plot" = Plot ,
                               "Year" = Year     )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Coll_all.csv")

### Mites
mites_traits = fread('Traits/Collembola_Mites/2019_Mite_fauna_Traits.csv')
mites_traits[, Species := V1]

# Correct some species names
mites_traits[, Species := recode(Species,
                                # Carabodes_forsslundi
                                 'Ceratozetella_sellnicki' = 'Ceratozetes_sellnicki',
                                 'Nanhermannia_nanus' = 'Nanhermannia_nana',
                                 #'Nothrus_borussicus'
                                 #'Peloptulus_montanus'
                                 'Phthiracarus_globosus' = 'Phthiracarus_globulus',
                                 'Poecilochthonius_italicus' = 'Poeliochthonius_italicus',
                                 'Poecilochthonius_spiciger' = 'Poeliochthonius_spiciger',
                                 'Scutovertex_sp.' = 'Scutovertex_sp',
                                 )]


mites_traits[, Habitat_spec := as.numeric(recode(Vertical_distribution.Morphotype,
                                     "non-spec" = '1',
                                     "soil" = '2',
                                     "surface" = '3',
                                     "litter" = '4'))]
mites_traits[, Feeding_spec := as.numeric(recode(Feeding,
                                                "omnivorous" = '1',
                                                "herbifungivorous" = '2',
                                                "herbivorous" = '3',
                                                "fungivorous" = '4'))]
mites_traits[, Repro_sex := as.numeric(ifelse(Reproduction == 'parthenogenetic', 0,
                                   ifelse( Reproduction == 'sex', 1,
                                           NA)))]

mites_traits[, Mass.log := log(Mass)]

pca_species = dudi.pca(mites_traits[,c( 'Mass.log', 'Feeding_spec', 'Habitat_spec', 'Repro_sex', 'DaystoAdult')], scannf = FALSE, nf = 3)
pca_mites_species = fviz_pca_biplot(pca_species, geom = c("point"), repel = T, axes = c(1,2), title = 'Mites')
ggsave(pca_mites_species,file= '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Species_PCA_Collembola.pdf', width = 5, height = 5)    




Abundances_mites = Abundances_all[Group_fine=='Oribatida', ]
Abundances_mites[!(Species %in% mites_traits$Species) & value > 0, unique(Species)]

# For species id to genus level, we take the average or value of other species from the same genus if it's consistent (mean > 2sd)
# Carabodes_sp: only one species
#mites_traits[grepl('Carabodes', Species),]
mites_traits = rbind(mites_traits, mites_traits[grepl('Carabodes', Species), list(V1 = 'Carabodes_sp', Vertical_distribution.Morphotype, Surface,Litter,  Soil, Feeding, Reproduction, Mass, size100_500, 
                                                                                  size500_900, size900,  ash,  CuticuleCalcium, CuticuleNo, DaystoAdult, Species = 'Carabodes_sp' , Habitat_spec,Feeding_spec, Repro_sex, Mass.log)])

# Carabodes_sp: only one species
#mites_traits[grepl('Nothrus', Species),]
mites_traits = rbind(mites_traits, mites_traits[grepl('Nothrus', Species), list(V1= 'Nothrus_sp', 
                                                                                Vertical_distribution.Morphotype = unique(Vertical_distribution.Morphotype), 
                                                                                Surface = unique(Surface),
                                                                                Litter = unique(Litter), 
                                                                                Soil= unique(Soil),
                                                                                Feeding= unique(Feeding), 
                                                                                Reproduction= unique(Feeding),
                                                                                Mass = mean(Mass), 
                                                                                size100_500= mean(size100_500), 
                                                                                size500_900= mean(size500_900), 
                                                                                size900= mean(size900), 
                                                                                ash = unique(ash), 
                                                                                CuticuleCalcium= unique(CuticuleCalcium), 
                                                                                CuticuleNo= unique(CuticuleNo),
                                                                                DaystoAdult= unique(DaystoAdult), 
                                                                                Species = 'Nothrus_sp' ,
                                                                                Habitat_spec= unique(Habitat_spec),
                                                                                Feeding_spec= unique(Feeding_spec),
                                                                                Repro_sex= unique(Repro_sex), 
                                                                                Mass.log= mean(Mass.log))])

#Pelotpulus_sp
mites_traits = rbind(mites_traits, mites_traits[grepl('Peloptulus', Species), list(V1 = 'Pelotpulus_sp', Vertical_distribution.Morphotype, Surface,Litter,  Soil, Feeding, Reproduction, Mass, size100_500, 
                                                                                   size500_900, size900,  ash,  CuticuleCalcium, CuticuleNo, DaystoAdult, Species = 'Pelotpulus_sp', Habitat_spec,Feeding_spec, Repro_sex, Mass.log)])

#Trichoribates_sp
mites_traits = rbind(mites_traits, mites_traits[grepl('Trichoribates', Species), list(V1 = 'Trichoribates_sp', Vertical_distribution.Morphotype, Surface,Litter,  Soil, Feeding, Reproduction, Mass, size100_500, 
                                                                                   size500_900, size900,  ash,  CuticuleCalcium, CuticuleNo, DaystoAdult, Species = 'Trichoribates_sp', Habitat_spec,Feeding_spec, Repro_sex, Mass.log)])



mites_traits[Species %in% Abundances_mites$Species, lapply(.SD, function(x){length(x[!is.na(x)])}),]

### Calculate coverage and CWM

Mites_CC = check_coverage(mites_traits, Abundances_mites, c( 'Mass.log', 'Feeding_spec', 'Habitat_spec', 'Repro_sex', 'DaystoAdult'),
                       'Species', 'Species')
Mites_CC[, average_coverage := mean(c(Mass.log,Feeding_spec, Habitat_spec, Repro_sex , DaystoAdult)), by = 1:nrow(Mites_CC)]
Mites_CC = Mites_CC[average_coverage > 0.2,]

Mites_CC[, -c(1, 2)][, lapply(.SD, min)]
Mites_CC[, -c(1, 2)][, lapply(.SD, median)]
Mites_CC[, -c(1, 2)][, lapply(.SD, max)]


Mites_cwm_all = my_cwm(mites_traits, Abundances_mites, c( 'Mass.log', 'Feeding_spec', 'Habitat_spec', 'Repro_sex', 'DaystoAdult'),
                       'Species', 'Species')

Mites_cwm_all = Mites_cwm_all[Plot %in% Mites_CC[average_coverage >= 0.2, Plot],]
 

pca_mites = dudi.pca(complete(mice(Mites_cwm_all[ ,c( 'Mass.log', 'Feeding_spec', 'Habitat_spec','DaystoAdult', 'Repro_sex_1')])), scannf = FALSE, nf = 3)
fviz_pca_biplot(pca_mites, geom = c("point"), repel = T, axes = c(1,2), title = 'Mites')

cor.test(pca_mites$l1$RS1, env_data_lui[Plot %in% Mites_cwm_all$Plot,]$LUI)


write.csv(Mites_cwm_all[, list( "mites_Hab_spec" = Habitat_spec,
                                "mites_Sex" = Repro_sex_1,
                                "mites_Mass" = Mass.log,
                                "mites_Feed_spec" = Feeding_spec,
                                "mites_DaysAdult" = DaystoAdult,
                                "Plot" = Plot,
                                "Year" = Year )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Mites_all.csv")


#### Myriapoda ####
# Traits from Birkhofer et al. 2017 (Araneae)
Birkhofer_myriapoda = data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Arthropods/Birkhofer_2017/Traits_Communities.xlsx", sheet = 'Chilopoda'))
Birkhofer_myriapoda[, new_scientificName := get_gbif(gsub('\\.', ' ', Site)), by = Site]
Abundance_myr = Abundances_all[Group_broad == 'Myriapoda' & Trophic_level == 'secondary.consumer.myriapod',]
Abundance_myr[, new_scientificName := get_gbif(Species), by = Species]

Abundance_myr[ , unique(new_scientificName[!new_scientificName %in% Birkhofer_myriapoda[ , unique(new_scientificName)]])]
Birkhofer_myriapoda[ , unique(new_scientificName)]


Birkhofer_myriapoda[, Size := weighted.mean(c(1, 2, 3), .SD), .SDcols = c('6-19mm', '19-32mm', '>32mm'), by = new_scientificName]
Birkhofer_myriapoda[, Mobility := weighted.mean(c(0, 1), .SD), .SDcols = c('mobility.ow', 'mobility.high'), by = new_scientificName]

Abundance_myr[, sum(value[new_scientificName %in% Birkhofer_myriapoda$new_scientificName])/sum(value), by = Plot][, range(V1, na.rm = T)]



#### Ants/Hymenoptera: looked at BETSI and Global Ant Database, no data for the species we need
# get_gbif_taxonomy(Abundances_all[Group_broad == 'Formicidae', unique(Species)])

#Abundances_all[Group_broad %in% c('Formicidae', 'Hymenoptera'), sum(value), by = Species][order(V1, decreasing = T),]

# Campanotus ligniperdus ,  Formica clara          ,  Formica fusca           , Formica pratensis     
# Formica sanguinea      ,  Lasius alienus         ,  Lasius myops            , Lasius paralienus     
# Lasius umbratus        ,  Myrmecina graminicola  ,  Myrmica curvithorax     , Myrmica lobicornis    
# Myrmica lonae          ,  Myrmica sabuleti       ,  Myrmica schenki         , Myrmica specioides    
# Tapinoma erractium     ,  Tapinoma subboreale    ,  Temnothorax unifasciatum, Tetramorium caespitum



### Ants ####
# Traits
#Heuss2019 <- data.table(read_excel("Traits/Ants/Heuss2019.xlsx", 
#                                   col_types = c("text", "numeric", "numeric", 
#                                                 "numeric", "numeric", "numeric", 
#                                                 "numeric", "numeric", "numeric", 
#                                                 "numeric", "numeric", "numeric")))
#Heuss2019[, Species := gsub(' ', '_', Species)]
#
#Heuss2019[Species == 'Campanotus_ligniperda', Species := 'Campanotus_ligniperdus']
#trait_ants = c('Nectar','Plant', 'Zoopha','Tropho','WL', 'Dom',    'CS',  'nQ',  'nN', 'CFT', 'Strata_forage')
#
#ant_abundance = fread('/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Ants/23986_2_Dataset/23986_2_data.csv')
#ant_abundance_melt = melt.data.table(ant_abundance, id.vars = c('Species', 'Plot'), measure.vars = 'Presence_absence')
#ant_abundance_melt[, Plot := ifelse(nchar(Plot) == 5, Plot, paste0(substr(Plot, 1, 3), '0', substr(Plot, 4, 4)))]
#ant_abundance_melt[, Year := '2014_15']
#ant_abundance_melt[, Species := gsub(' ', '_', Species)]
#
#
#CWM_ants = my_cwm(Traits0 = Heuss2019, ant_abundance_melt, trait_names = trait_ants,
#                  trait_taxo = 'Species', abundance_taxo = 'Species')
#
#CWM_ants_above = my_cwm(Traits0 = Heuss2019[Strata_forage>0,], ant_abundance_melt, trait_names = trait_ants,
#                  trait_taxo = 'Species', abundance_taxo = 'Species')
#
#CWM_ants_below = my_cwm(Traits0 = Heuss2019[Strata_forage<0,], ant_abundance_melt, trait_names = trait_ants,
#                        trait_taxo = 'Species', abundance_taxo = 'Species')
#
#fwrite(CWM_ants[, list(
#  'ant_zoopha' = Zoopha,
#  "ant_nectar" = Nectar      ,
#  "ant_plant" = Plant    ,
#  "ant_tropho" = Tropho       ,
#  "ant_length" = WL     ,
#  "ant_dom" = Dom_1       ,
#  "ant_colsize" = CS        ,
#  "ant_polygyny" = nQ         ,
#  "ant_polydomy" = nN   ,
#  "ant_colfound" = CFT ,
#  "ant_forageabove" = Strata_forage  ,
#   Plot,
#  Year = 2014
#)], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_ants.csv")
#
#fwrite(CWM_ants_above[, list(
#  'ant_a_zoopha' = Zoopha,
#  "ant_a_nectar" = Nectar      ,
#  "ant_a_plant" = Plant    ,
#  "ant_a_tropho" = Tropho       ,
#  "ant_a_length" = WL     ,
#  "ant_a_dom" = Dom_1       ,
#  "ant_a_colsize" = CS        ,
#  "ant_a_polygyny" = nQ         ,
#  "ant_a_polydomy" = nN   ,
#  "ant_a_colfound" = CFT ,
#  "ant_a_forageabove" = Strata_forage  ,
#  Plot,
#  Year = 2014
#)], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_ants_above.csv")
#
#fwrite(CWM_ants_below[, list(
#  'ant_b_zoopha' = Zoopha,
#  "ant_b_nectar" = Nectar      ,
#  "ant_b_plant" = Plant    ,
#  "ant_b_tropho" = Tropho       ,
#  "ant_b_length" = WL     ,
#  "ant_b_dom" = Dom_1       ,
#  "ant_b_colsize" = CS        ,
#  "ant_b_polygyny" = nQ         ,
# # "ant_polydomy" = nN   ,
#  "ant_b_colfound" = CFT ,
#  "ant_b_forageabove" = Strata_forage  ,
#  Plot,
#  Year = 2014
#)], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_ants_below.csv")



