
# This script demonstrate the analysis conducted for the manuscript "A fast-slow trait continuum at the level of entire communities" by Neyret et al. 

# Author: Margot Neyret - Please get in touch if you have questions.

# This script takes as input the abundances and species-level traits of arthropod species
# and outputs a CWM matrix averaged for all considered years.


##############################################
#### ***** HEMIPTERA COLEOPTERA ETC ***** ####
##############################################

# ********************************** #
#### 1. Load and merge trait data ####
# ********************************** #

## Traits
# Published in Goessner 2015 - it has additional traits compared to the new version
Arthropod_traits_Goessner = setDT(read_delim("Data/Trait_data/ArthropodSpeciesTraits.csv", ";", escape_double = FALSE, trim_ws = TRUE)) # https://datadryad.org/stash/dataset/doi:10.5061/dryad.53ds2
Arthropod_traits_Goessner[, gbifName := get_gbif(SpeciesID), by = SpeciesID]

# New version - most complete trait dataset
Arthropod_traits_grasslands_Raw = fread('Data/Trait_data/31122_6_Dataset/31122_6_data.csv') # https://www.bexis.uni-jena.de/ddm/data/Showdata/31122
Arthropod_traits_grasslands_Raw = unique(Arthropod_traits_grasslands_Raw)
Arthropod_traits_grasslands_Raw[, gbifName := get_gbif(SpeciesID), by = SpeciesID]
Arthropod_traits_grasslands_Raw[, Dispersal_ability := as.numeric(gsub(',', '.', Dispersal_ability))]

# Generation time for herbivore (compiled from different datasets)
Voltinism_traits <- data.table(read_excel("Data/Additional_data/Voltinism_herbivore.xlsx")) # Available in Additional_data folder
Voltinism_traits[, gbifName := get_gbif(Species_raw), by = Species_raw]
Voltinism_traits[, Generations := as.numeric(Generations)]
Voltinism_traits = Voltinism_traits[!is.na(gbifName),]

### Abundances
Abundances_temporal = fread('Data/Abundance_data/21969_4_Dataset/21969_4_data.csv') #https://www.bexis.uni-jena.de/ddm/data/Showdata/21969
Abundances_temporal = Abundances_temporal[!is.na(Species),]
Abundances_temporal[, Plot := ifelse(nchar(PlotID) == 5, PlotID, paste(substr(PlotID, 1, 3), '0', substr(PlotID, 4, 4), sep = ''))]
Abundances_temporal_by_plot = Abundances_temporal[, list(value = sum(NumberAdults, na.rm = T)), by = list('Plot' = Plot, 'Species' = Species, 'Year' = CollectionYear, 'Order' = Order)]
Abundances_temporal_by_plot[, gbifName := get_gbif(Species), by = Species]

###### Data cleaning and standardization ####### 
# One species can't be found by the function get_gbif for some reason
Abundances_temporal_by_plot[Species == 'Pachybrachius fracticollis' , new_scientificName := 'Pachymerus fracticollis']
Arthropod_traits_grasslands_Raw[SpeciesID == 'Pachybrachius fracticollis' , new_scientificName := 'Pachymerus fracticollis']

# Merge trait datasets
Arthropod_traits_grasslands = merge(Arthropod_traits_grasslands_Raw, Arthropod_traits_Goessner[, list(Feeding_guild_Goessner = Feeding_guild, Feeding_specialization, gbifName)], all.x = T)
Arthropod_traits_grasslands = merge(Arthropod_traits_grasslands, Voltinism_traits[, list(Generations, gbifName)], all.x = T)
Arthropod_traits_grasslands = unique(Arthropod_traits_grasslands)

# Merging primarily on the names originally provided in the datasets (which fit in most cases).
# If no fit is found, I use the names provided by the gbif function.
# This is done because in some cases the gbif function merges synonyms which are considered separately in the trait databases
Arthropod_traits_grasslands[, new_scientificName := ifelse(SpeciesID %in% Abundances_temporal_by_plot$Species, SpeciesID, gbifName), by = SpeciesID]
Abundances_temporal_by_plot[, new_scientificName := ifelse(Species %in% Arthropod_traits_grasslands_Raw$SpeciesID, Species, gbifName), by = Species]


# ************************** #
#### 2.Recode trait data ####
# ************************* #
# Carnivores are secondary consumers and the rest are primary consumers
Arthropod_traits_grasslands[, Feeding_guild_simple := dplyr::recode(Feeding_guild, 
                                                              "h"        = 'primary',
                                                              "f"        = 'primary', 
                                                              "m"        = 'primary',
                                                              'd'        = 'primary',
                                                              "c"       = 'secondary',
                                                              )]

# We classify species which are carnivores at at least one of their life stages as secondary consumers, the other as primary consumers. 
Arthropod_traits_grasslands[Feeding_guild == 'o' & grepl('c', Feeding_guild_Goessner) & Feeding_guild_Goessner != 'h-(c)', Feeding_guild_simple := 'secondary']
Arthropod_traits_grasslands[Feeding_guild == 'o' & Feeding_guild_simple == 'o', Feeding_guild_simple := 'primary']

# Simplify stratum_use into below-ground (soil, ground) and above-ground
Arthropod_traits_grasslands[, Stratum_use_simple := dplyr::recode(Stratum_use, 
                                                                    "g" = "below",
                                                                    "h" = "above",
                                                                    "i" = "above",
                                                                    "s" = "below",
                                                                    "t" = "above",
                                                                    "u" = "NA",
                                                                    "w" = "NA")]

# We transform feeding specialization into a numerical variable
Arthropod_traits_grasslands[, Feeding_generalism := as.numeric(factor(Feeding_specialization, levels = c('m', 'o', 'p')))]

##### Add taxa identified only at genus level ####
# Some taxa in abundances are not ID at species level. We can partly interpolate the species-level data to genus level.
un_id_species = unique(Abundances_temporal_by_plot[!(Species %in% Arthropod_traits_grasslands$SpeciesID), list(Species, new_scientificName)])
un_id_species[, genus := gsub(' sp.', '', new_scientificName)]
un_id_species[, Genus := gsub(' sp.', '', new_scientificName)]
un_id_species[, SpeciesID := Species]

# Order, Family
un_id_species[, c( 'Order', 'Family') := list(unique( unique(Arthropod_traits_grasslands[Genus == genus, Order]),
                                                          unique(Arthropod_traits_grasslands[Genus == genus, Family]))), by = genus]

# Numeric traits: we attribute genus-level values to unidentified species if all data is consistent, i.e. mean > 2*sd - otherwise NA
un_id_species[, c('Mean_BodySize', 'Dispersal_ability', 'Generations') := list(as.numeric(ifelse( Arthropod_traits_grasslands[Genus == genus, mean(Mean_BodySize, na.rm = T)] > 2*Arthropod_traits_grasslands[Genus == genus, sd(Mean_BodySize, na.rm = T)],
                                                                                  Arthropod_traits_grasslands[Genus == genus, mean(Mean_BodySize, na.rm = T)],
                                                                                  NA)),
                                                                           as.numeric(ifelse(Arthropod_traits_grasslands[Genus == genus, mean(Dispersal_ability, na.rm = T)] > 2* Arthropod_traits_grasslands[Genus == genus, sd(Dispersal_ability, na.rm = T)],
                                                                                  Arthropod_traits_grasslands[Genus == genus, mean(Dispersal_ability, na.rm = T)],
                                                                                  NA)),
                                                                           as.numeric(ifelse(Arthropod_traits_grasslands[Genus == genus, mean(Generations, na.rm = T)] > 2* Arthropod_traits_grasslands[Genus == genus, sd(Generations, na.rm = T)],
                                                                                  Arthropod_traits_grasslands[Genus == genus, mean(Generations, na.rm = T)],
                                                                                  NA))), by = genus]

# Qualitative traits traits: we attribute genus-level values to unidentified species if all data is consistent, i.e. 90% species share the same trait value
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
un_id_species[, Stratum_use := as.character(ifelse(
  Arthropod_traits_grasslands[Genus == genus, max(table(Stratum_use[!is.na(Stratum_use)]), na.rm = T)>0.9*sum(table(Stratum_use[!is.na(Stratum_use)]), na.rm = T)],
  Arthropod_traits_grasslands[Genus == genus, names(sort(table(Stratum_use[!is.na(Stratum_use)]), decreasing = T))[1]],
  NA)), by = genus]
un_id_species[, Feeding_guild := as.character(ifelse(
  Arthropod_traits_grasslands[Genus == genus, max(table(Feeding_guild[!is.na(Feeding_guild)]), na.rm = T)>0.9*sum(table(Feeding_guild[!is.na(Feeding_guild)]), na.rm = T)],
  Arthropod_traits_grasslands[Genus == genus, names(sort(table(Feeding_guild[!is.na(Feeding_guild)]), decreasing = T))[1]],
  NA)), by = genus]


un_id_species = un_id_species[new_scientificName != " sp." ]

# Merge with initial dataset
Arthropod_traits = unique(rbind(Arthropod_traits_grasslands[, list(SpeciesID, Order, Family, Genus, Species, Mean_BodySize, Dispersal_ability, Feeding_guild, Stratum_use, Feeding_guild_simple, Stratum_use_simple, Feeding_generalism, Generations, new_scientificName)],
                                              un_id_species[, list(SpeciesID, Order, Family, Genus, Species, Mean_BodySize, Dispersal_ability, Feeding_guild, Stratum_use, Feeding_guild_simple, Stratum_use_simple, Feeding_generalism, Generations, new_scientificName)]))
Arthropod_traits = unique(Arthropod_traits)
Arthropod_traits = Arthropod_traits[!is.na(new_scientificName),]

Arthropod_traits$Feeding_generalism = as.numeric(Arthropod_traits$Feeding_generalism)

# Save species-matched trait data
traitUnits = c("unitless","mm","unitless",'unitless', 'unitless',"unitless")
traitDescription = c("Dispersal ability ranging from 0 to 1; definition differs across taxonomic groups (see Goessner et al. 2016)",
                     "Mean body length",
                     "Main stratum used during life cycle  'g' (ground), '' (herb layer), 's' (soil), 't' (shrub and tree layer), 'u' (unspecific) or 'w' (water)",
                     "Feeding guild ('c' (carnivore), 'h' (herbivore), 'o' (omnivore), 'd' (detritivore), 'm' (mycetophagous) or 'f' (fungivore))",
                     "Feeding generalism index coded as 1 = monophage (feeding on 1 genus), 2 = oligophage (on one plant lineage), polyphages (on more than one plant lineage) ",
                     "Number of generations per year"
)
traitDataID = c("Bexis ID 31122","Bexis ID 31122", "Bexis ID 31122","Bexis ID 31122","Bexis ID 31122",'NA')

traitRef = c("DOI: 10.1038/sdata.2015.13","DOI: 10.1038/sdata.2015.13","DOI: 10.1038/sdata.2015.13","DOI: 10.1038/sdata.2015.13","DOI: 10.1038/sdata.2015.13",
             'Bakewell, A.T., Davis, K.E., Freckleton, R.P., Isaac, N.J.B., Mayhew, P.J., 2020. Comparing Life Histories across Taxonomic Groups in Multiple Dimensions: How Mammal-Like Are Insects? The American Naturalist 195, 70–81. https://doi.org/10.1086/706195; http://michentsoc.org/gle-pdfs/vol26no2.pdf, https://andrewsforest.oregonstate.edu/sites/default/files/lter/pubs/pdf/pub1895.pdf, https://wiki.pestinfo.org/wiki/Brassicogethes_aeneus; Neff, F., M.C. Resch, A. Marty, J.D. Rolley, M. Schütz, A.C. Risch, and M.M. Gossner. 2020. Long-term restoration success of insect herbivore communities in semi-natural grasslands: a functional approach. Ecological Applications; Nickel, H. 2003. The leafhoppers and planthoppers of Germany (Hemiptera, Auchenorrhyncha): Goecke & Evers, Sofia - Moscow / Keltern; ; Saulich, A., Musolin, D., 2021. Seasonal Development of Plant Bugs (Heteroptera, Miridae): Subfamily Mirinae, Tribe Stenodemini. Entomological Review 101, 147–161. https://doi.org/10.1134/S0013873821020019; Wipfli, Mark S.; Wedberg, John L.; Hogg, David B.; and Syverud, Thomas D. 1989. "Insect Pests Associated With Birdsfoot Trefoil, Lotus Corniculatus, in Wisconsin," The Great Lakes Entomologist, vol 22 (1)')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c('Dispersal_ability', 'Mean_BodySize','Stratum_use', 'Feeding_guild',  'Feeding_generalism', 'Generations'  )

Arthropod_traits_melt = melt.data.table(unique(Arthropod_traits[, .SD, .SDcols = c('Dispersal_ability', 'Mean_BodySize','Stratum_use', 'Feeding_guild',  'Feeding_generalism', 'Generations', 'new_scientificName' )]), id.vars = c( 'new_scientificName'), variable.name = 'traitName', value.name = 'traitValue')
Arthropod_traits_info = add_info(Arthropod_traits_melt[traitName %in% c('Dispersal_ability', 'Mean_BodySize','Stratum_use', 'Feeding_guild',  'Feeding_generalism', 'Generations' )], traitRefs = traitRef, traitDataIDs = traitDataID, traitDescriptions = traitDescription, traitUnits = traitUnits, traitsOnly = TRUE)

fwrite(Arthropod_traits_info, "Data/Temporary_data/Arthropod_traits.csv")


###### Check trait distribution
Arthropod_traits[, logBody_Size := log(Mean_BodySize)]

# ************************** #
#### 3. Species-level PCA ####
# ************************** #

### Species-level PCAs, coverage and CWM per trophic level
# Some species from abundance are not included in trait so need to adapt the number reported:
# AG, herb: + 1 Auchenorrhyncha spec., Typhlocybinae spec., Stenodemini, Deltocephalinae spec."

# Herbivore, AG
traits_h_AG = Arthropod_traits[Stratum_use_simple == 'above'  & 
                                          Feeding_guild_simple == 'primary' &
                                          new_scientificName %in% Abundances_temporal_by_plot$new_scientificName, .SD, 
                               .SDcols = c('Feeding_generalism', 'Dispersal_ability', 'logBody_Size', 'Generations')]

pca_h_AG_sp = dudi.pca(mice::complete(mice(traits_h_AG[, list(Feeding_generalism, Dispersal_ability,Body_logSize = logBody_Size, Generations)])), scannf = FALSE, nf = 2)
gg_h_AG = fviz_pca(pca_h_AG_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                  col.ind = "steelblue",
                  fill.ind = "white",
                  col.var = "black")
ggsave(gg_h_AG, file = 'Results/species_pca_arth_h_AG.pdf',  width = 6, height = 5)


# Carnivores, AG ---> Here the fast-slow axis is axis 2
trait_c_AG = Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName,lapply(.SD, as.numeric), .SDcols = c('Dispersal_ability', 'logBody_Size')]

trait_c_AG[, lapply(.SD, function(x){length(x[!is.na(x)])})]

pca_c_AG_sp = dudi.pca(mice::complete(mice(trait_c_AG[, list( Dispersal_ability, Body_logSize = logBody_Size)])), scannf = FALSE, nf = 2)
gg_c_AG = fviz_pca(pca_c_AG_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                   col.ind = "steelblue",
                   fill.ind = "white",
                   col.var = "black")
ggsave(gg_c_AG, file ='Results/species_pca_arth_c_AG.pdf', width = 6, height = 5)

# Herbivores, BG
trait_h_BG = Arthropod_traits[Stratum_use_simple == 'below' & Feeding_guild_simple  == 'primary'  &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName,lapply(.SD, as.numeric), 
                              .SDcols = c('Dispersal_ability', 'logBody_Size','Feeding_generalism')]
trait_h_BG[, lapply(.SD, function(x){length(x[!is.na(x)])})]

pca_h_BG_sp = dudi.pca(mice::complete(mice(trait_h_BG[, list(Feeding_generalism, Dispersal_ability, Body_logSize = logBody_Size)])), scannf = FALSE, nf = 2)
gg_h_BG_sp = fviz_pca(pca_h_BG_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                   col.ind = "steelblue",
                   fill.ind = "white",
                   col.var = "black")
ggsave(gg_h_BG_sp, file ='Results/species_pca_arth_h_BG.pdf', width = 6, height = 5)

# Carnivores, BG
trait_c_BG = Arthropod_traits[Stratum_use_simple == 'below' & Feeding_guild_simple == 'secondary' &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName,lapply(.SD, as.numeric), .SDcols = c('Dispersal_ability', 'logBody_Size','Feeding_generalism')]

pca_c_BG_sp = dudi.pca(mice::complete(mice(trait_c_BG[, list(Dispersal_ability, Body_logSize = logBody_Size)])), scannf = FALSE, nf = 2)
gg_c_BG_sp = fviz_pca(pca_c_BG_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                      col.ind = "steelblue",
                      fill.ind = "white",
                      col.var = "black")
ggsave(gg_c_BG_sp, file ='Results/species_pca_arth_c_BG.pdf', width = 6, height = 5)


# ********************************** #
#### 4. Community-weighted traits ####
# ********************************** #

# Herbivore, AG
Coverage_Arthropods_above_herb = check_coverage(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' ,], 
                                                Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' , SpeciesID] |
                                                                              new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' , new_scientificName]], c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Generations'), 'new_scientificName', 'new_scientificName')

CWM_Arthropods_above_herb = my_cwm(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' ,], 
                                   Abundances_temporal_by_plot, c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Generations'), 'new_scientificName', 'new_scientificName')

# Carnivores, AG ---> Here the fast-slow axis is axis 2
Coverage_Arthropods_above_omni_carni = check_coverage(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' ,], 
                                                 Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' , SpeciesID] |
                                                                               new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' , new_scientificName], 
                                                                             list(value, Year), by = c('new_scientificName', 'Plot')],
                                                 c( 'Dispersal_ability',"logBody_Size"), 'new_scientificName', 'new_scientificName')

CWM_Arthropods_above_omni_carni = my_cwm(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple  == 'secondary',], 
                                         Abundances_temporal_by_plot, c('Dispersal_ability',"logBody_Size"), 'new_scientificName', 'new_scientificName')

# Herbivores, BG
Coverage_Arthropods_below_herb = check_coverage(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' ,], 
                                                Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' , SpeciesID] |
                                                                              new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' , new_scientificName], 
                                                                            list(value, Year), by = c('new_scientificName', 'Plot')],
                                                c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism'), 'new_scientificName', 'new_scientificName')
CWM_Arthropods_below_herb = my_cwm(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' ,], 
                                   Abundances_temporal_by_plot, c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Generations'), 'new_scientificName', 'new_scientificName')


# Carnivores, BG
Coverage_Arthropods_below_carni = check_coverage(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' ,], 
                                                 Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' , SpeciesID] |
                                                                               new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' , new_scientificName], 
                                                                             list(value, Year), by = c('new_scientificName', 'Plot')],
                                                 c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism'), 'new_scientificName', 'new_scientificName')
CWM_Arthropods_below_omni_carni  = my_cwm(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' ,], 
                                          Abundances_temporal_by_plot, c( 'Dispersal_ability',"logBody_Size"), 'SpeciesID', 'Species')


# Melt and merge
CWM_CC_Arthropods_below_herb= merge.data.table(melt.data.table(CWM_Arthropods_below_herb[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                               melt.data.table(Coverage_Arthropods_below_herb[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))

CWM_CC_Arthropods_above_carni = merge.data.table(melt.data.table(CWM_Arthropods_above_omni_carni[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                                melt.data.table(Coverage_Arthropods_above_omni_carni[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


CWM_CC_Arthropods_above_herb = merge.data.table(melt.data.table(CWM_Arthropods_above_herb[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                               melt.data.table(Coverage_Arthropods_above_herb[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))

CWM_CC_Arthropods_below_omni_carni = merge.data.table(melt.data.table(CWM_Arthropods_below_omni_carni[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                       melt.data.table(Coverage_Arthropods_below_carni[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


# Add info
traitUnits = c("unitless","mm","unitless","unitless")
traitDescription = c("Dispersal ability ranging from 0 to 1; definition differs across taxonomic groups (see Goessner et al. 2016)",
                     "Mean body length",
                     "Feeding generalism index coded as 1 = monophage (feeding on 1 genus), 2 = oligophage (on one plant lineage), polyphages (on more than one plant lineage) ",
                     "Number of generations per year"
                     )
traitDataID = c("Bexis ID 31122","Bexis ID 31122","Bexis ID 31122",'')

traitRef = c("DOI: 10.1038/sdata.2015.13","DOI: 10.1038/sdata.2015.13","DOI: 10.1038/sdata.2015.13",
             'Bakewell, A.T., Davis, K.E., Freckleton, R.P., Isaac, N.J.B., Mayhew, P.J., 2020. Comparing Life Histories across Taxonomic Groups in Multiple Dimensions: How Mammal-Like Are Insects? The American Naturalist 195, 70–81. https://doi.org/10.1086/706195; http://michentsoc.org/gle-pdfs/vol26no2.pdf, https://andrewsforest.oregonstate.edu/sites/default/files/lter/pubs/pdf/pub1895.pdf, https://wiki.pestinfo.org/wiki/Brassicogethes_aeneus; Neff, F., M.C. Resch, A. Marty, J.D. Rolley, M. Schütz, A.C. Risch, and M.M. Gossner. 2020. Long-term restoration success of insect herbivore communities in semi-natural grasslands: a functional approach. Ecological Applications; Nickel, H. 2003. The leafhoppers and planthoppers of Germany (Hemiptera, Auchenorrhyncha): Goecke & Evers, Sofia - Moscow / Keltern; ; Saulich, A., Musolin, D., 2021. Seasonal Development of Plant Bugs (Heteroptera, Miridae): Subfamily Mirinae, Tribe Stenodemini. Entomological Review 101, 147–161. https://doi.org/10.1134/S0013873821020019; Wipfli, Mark S.; Wedberg, John L.; Hogg, David B.; and Syverud, Thomas D. 1989. "Insect Pests Associated With Birdsfoot Trefoil, Lotus Corniculatus, in Wisconsin," The Great Lakes Entomologist, vol 22 (1)')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c('Dispersal_ability', 'logBody_Size',  'Feeding_generalism', 'Generations'  )

CWM_CC_Arthropods_above_herb = add_info(CWM_CC_Arthropods_above_herb, traitRef, traitDataID, traitDescription, traitUnits, c('21969'))
CWM_CC_Arthropods_above_carni = add_info(CWM_CC_Arthropods_above_carni, traitRef, traitDataID, traitDescription, traitUnits, c('21969'))
CWM_CC_Arthropods_below_herb = add_info(CWM_CC_Arthropods_below_herb, traitRef, traitDataID, traitDescription, traitUnits, c('21969'))
CWM_CC_Arthropods_below_omni_carni = add_info(CWM_CC_Arthropods_below_omni_carni, traitRef, traitDataID, traitDescription, traitUnits, c('21969'))

fwrite(CWM_CC_Arthropods_above_herb, "Data/CWM_data/CWM_arthropods_above_herb.csv")
fwrite(CWM_CC_Arthropods_above_carni, "Data/CWM_data/CWM_arthropods_above_carni.csv")
fwrite(CWM_CC_Arthropods_below_herb, "Data/CWM_data/CWM_arthropods_below_herb.csv")
fwrite(CWM_CC_Arthropods_below_omni_carni, "Data/CWM_data/CWM_arthropods_below_omni_carni.csv")


# ************************************** #
#### 5. Non-weighted community traits ####
# ************************************** #
Abundances_temporal_by_plot_presence_absence = Abundances_temporal_by_plot[, list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Species , Order, gbifName, new_scientificName)]
Abundances_temporal_by_plot_presence_absence[value>1, value := 1]

# Herbivore, AG
Coverage_Arthropods_above_herb_noweight = check_coverage(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' ,], 
                                                Abundances_temporal_by_plot_presence_absence[Species %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' , SpeciesID] |
                                                                              new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' , new_scientificName]], c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Generations'), 'new_scientificName', 'new_scientificName')


CWM_Arthropods_above_herb_noweight = my_cwm(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' ,], 
                                   Abundances_temporal_by_plot_presence_absence, c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Generations'), 'new_scientificName', 'new_scientificName')

# Carnivores, AG 
Coverage_Arthropods_above_omni_carni_noweight = check_coverage(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' ,], 
                                                      Abundances_temporal_by_plot_presence_absence[Species %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' , SpeciesID] |
                                                                                    new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' , new_scientificName], 
                                                                                  list(value, Year), by = c('new_scientificName', 'Plot')],
                                                      c( 'Dispersal_ability',"logBody_Size"), 'new_scientificName', 'new_scientificName')

CWM_Arthropods_above_omni_carni_noweight = my_cwm(Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple  == 'secondary',], 
                                         Abundances_temporal_by_plot_presence_absence, c('Dispersal_ability',"logBody_Size"), 'new_scientificName', 'new_scientificName')

# Herbivores, BG
Coverage_Arthropods_below_herb_noweight = check_coverage(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' ,], 
                                                Abundances_temporal_by_plot_presence_absence[Species %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' , SpeciesID] |
                                                                              new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' , new_scientificName], 
                                                                            list(value, Year), by = c('new_scientificName', 'Plot')],
                                                c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism'), 'new_scientificName', 'new_scientificName')
CWM_Arthropods_below_herb_noweight = my_cwm(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' ,], 
                                   Abundances_temporal_by_plot_presence_absence, c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Generations'), 'new_scientificName', 'new_scientificName')


# Carnivores, BG
Coverage_Arthropods_below_carni_noweight = check_coverage(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' ,], 
                                                 Abundances_temporal_by_plot_presence_absence[Species %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' , SpeciesID] |
                                                                               new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' , new_scientificName], 
                                                                             list(value, Year), by = c('new_scientificName', 'Plot')],
                                                 c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism'), 'new_scientificName', 'new_scientificName')
CWM_Arthropods_below_omni_carni_noweight  = my_cwm(Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' ,], 
                                          Abundances_temporal_by_plot_presence_absence, c( 'Dispersal_ability',"logBody_Size"), 'SpeciesID', 'Species')




# Melt and merge
CWM_CC_Arthropods_below_herb_noweight = merge.data.table(melt.data.table(CWM_Arthropods_below_herb_noweight[,Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                               melt.data.table(Coverage_Arthropods_below_herb_noweight[,Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))

CWM_CC_Arthropods_above_carni_noweight = merge.data.table(melt.data.table(CWM_Arthropods_above_omni_carni_noweight[, Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                                 melt.data.table(Coverage_Arthropods_above_omni_carni_noweight[, Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


CWM_CC_Arthropods_above_herb_noweight = merge.data.table(melt.data.table(CWM_Arthropods_above_herb_noweight[, Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                                melt.data.table(Coverage_Arthropods_above_herb_noweight[, Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))

CWM_CC_Arthropods_below_omni_carni_noweight = merge.data.table(melt.data.table(CWM_Arthropods_below_omni_carni_noweight[, Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                                      melt.data.table(Coverage_Arthropods_below_carni_noweight[, Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


# Add info

CWM_CC_Arthropods_above_herb_noweight = add_info(CWM_CC_Arthropods_above_herb_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('21969'))
CWM_CC_Arthropods_above_carni_noweight = add_info(CWM_CC_Arthropods_above_carni_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('21969'))
CWM_CC_Arthropods_below_herb_noweight = add_info(CWM_CC_Arthropods_below_herb_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('21969'))
CWM_CC_Arthropods_below_omni_carni_noweight = add_info(CWM_CC_Arthropods_below_omni_carni_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('21969'))

fwrite(CWM_CC_Arthropods_above_herb_noweight, "Data/CWM_data/CWM_arthropods_above_herb_noweight.csv")
fwrite(CWM_CC_Arthropods_above_carni_noweight, "Data/CWM_data/CWM_arthropods_above_carni_noweight.csv")
fwrite(CWM_CC_Arthropods_below_herb_noweight, "Data/CWM_data/CWM_arthropods_below_herb_noweight.csv")
fwrite(CWM_CC_Arthropods_below_omni_carni_noweight, "Data/CWM_data/CWM_arthropods_below_omni_carni_noweight.csv")

# ***************************** #
#### 5. Turnover (Table S6) #### 
# ***************************** #

data_lui <- fread("Data/Environment_function_data/LUI_standardized_global.txt") # from https://www.bexis.uni-jena.de/lui/LUICalculation/index; new components, standardised, global, all regions, all years
data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

comm.testaH_AG = dcast(Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' , SpeciesID] |
                              new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'primary' , new_scientificName], 
                            list(value = sum(value, na.rm = T)), by = c('new_scientificName', 'Plot')], Plot ~ new_scientificName, fill = 0)
comm.testaH_BG = dcast(Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' , SpeciesID] |
                                                     new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'primary' , new_scientificName], 
                                                   list(value = sum(value, na.rm = T)), by = c('new_scientificName', 'Plot')], Plot ~ new_scientificName, fill = 0)
comm.testaC_AG = dcast(Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' , SpeciesID] |
                                                     new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'above'  &  Feeding_guild_simple == 'secondary' , new_scientificName], 
                                                   list(value = sum(value, na.rm = T)), by = c('new_scientificName', 'Plot')], Plot ~ new_scientificName, fill = 0)
comm.testaC_BG = dcast(Abundances_temporal_by_plot[Species %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' , SpeciesID] |
                                                     new_scientificName %in% Arthropod_traits[Stratum_use_simple == 'below'  &  Feeding_guild_simple == 'secondary' , new_scientificName], 
                                                   list(value = sum(value, na.rm = T)), by = c('new_scientificName', 'Plot')], Plot ~ new_scientificName, fill = 0)

rownames(comm.testaH_AG)= comm.testaH_AG$Plot
rownames(comm.testaH_BG)= comm.testaH_BG$Plot
rownames(comm.testaC_AG)= comm.testaC_AG$Plot
rownames(comm.testaC_BG)= comm.testaC_BG$Plot

comm.testaH_AG = comm.testaH_AG[,-1]
comm.testaH_BG = comm.testaH_BG[,-1]
comm.testaC_AG = comm.testaC_AG[,-1]
comm.testaC_BG = comm.testaC_BG[,-1]

beta.multi.abund(comm.testaH_AG)
beta.multi.abund(comm.testaH_BG)
beta.multi.abund(comm.testaC_AG)
beta.multi.abund(comm.testaC_BG)


commaH_AG_min_max = matrix(c(colSums(comm.testaH_AG[min_lui_plots,]),colSums(comm.testaH_AG[max_lui_plots,])), nrow = 2)
commaH_BG_min_max = matrix(c(colSums(comm.testaH_BG[min_lui_plots,], na.rm = T),colSums(comm.testaH_BG[max_lui_plots,], na.rm = T)), nrow = 2)
commaC_AG_min_max = matrix(c(colSums(comm.testaC_AG[min_lui_plots,]),colSums(comm.testaC_AG[max_lui_plots,])), nrow = 2)
commaC_BG_min_max = matrix(c(colSums(comm.testaC_BG[min_lui_plots,]),colSums(comm.testaC_BG[max_lui_plots,])), nrow = 2)

beta.multi.abund(commaH_AG_min_max)
beta.multi.abund(commaH_BG_min_max)
beta.multi.abund(commaC_AG_min_max)
beta.multi.abund(commaC_BG_min_max)


################################
#### ***** COLLEMBOLA ***** ####
################################

# ********************************** #
#### 1. Load and merge trait data ####
# ********************************** #
coll_traits = data.table(read_excel(('Data/Trait_data/Coll_Traits_NEW.xlsx')))
coll_traits[, c( 'Size', 'Size_Adult', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales') := lapply(.SD, as.numeric),
            .SDcols = c( 'Size','Size_Adult',  'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales')]

# Check distribution
coll_traits[, lapply(.SD, function(x){hist(log(x))}), .SDcols = c( 'Size', 'Size_Adult', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales')]
coll_traits[, c( "logSize"):= lapply(.SD, log), .SDcols = c( "Size_Adult")]
coll_traits[, Depth_preference := ifelse(Vertical_preference == 'Euedaphic', 2,
                                  ifelse( Vertical_preference == 'Epiedaphic',0,
                                  1))]
coll_traits[, Repro_sex := as.numeric(Reproduction_Ruslan == 'bisex')]
coll_traits[, Gen_per_year := ifelse(Phenology == 'multivoltine', 3, 
                                     ifelse(Phenology == 'bivoltine', 2, 
                                            NA))]
coll_traits[, Species := gsub(' ', '_', Species)]
coll_traits[, Species := gsub('__', '_', Species)]

# Save species-matched trait data
# Add info
# logSize         Gen_per_year     Depth_preference Ocelli           Repro_sex        Pigment          Furca            PAO             
# PSO              Asp              Scales   
traitUnits = c("mm (log-transformed)","unitless","unitless","unitless","unitless", 'unitless', 'unitless','unitless', 'unitless', 'unitless', 'unitless')
traitDescription = c("Body size",
                     "Number of generations per year (coded as bivoltine = 2, multivoltine = 3",
                     "Depth preference (coded as Euedaphic = 2, hemiedaphic = 1, Epiedaphic=0)",
                     "Number of ocelli",
                     "Reproduction type (bisexual/parthenogenetic)",
                     "Presence of pigments",
                     "Presence of a furca",
                     "Presence of a Post Antennal Organ",
                     "Presence of pseudocelli",
                     "Presence of anal spines",
                     "Presence of scales")
traitDataID = c("NA","NA","NA","NA", "NA",'NA','NA','NA','NA','NA', 'NA')

traitRef = c("Saifutdinov, Zaytsev, John, Baulechner & Wolters","Saifutdinov, John, Baulechner & Wolters","Saifutdinov, John, Baulechner & Wolters","Saifutdinov, John, Baulechner & Wolters", "Polierer, Scheu, 2016; Chernova et al., 2010; Saifutdinov et a., 2018",'Saifutdinov, John, Baulechner & Wolters','Saifutdinov, John, Baulechner & Wolters','Saifutdinov, John, Baulechner & Wolters','Saifutdinov, John, Baulechner & Wolters','Saifutdinov, John, Baulechner & Wolters', 'Saifutdinov, John, Baulechner & Wolters')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c('Size_Adult', 'Phenology', 'Vertical_preference', 'Ocelli', 'Reproduction_Ruslan','Pigment','Furca', 'PAO','PSO','Asp', 'Scales')

coll_traits_melt = melt.data.table(coll_traits[, list(Size = Size_Adult, Phenology, Vertical_preference, Ocelli, Reproduction = Reproduction_Ruslan,Pigment,Furca, PAO,PSO,Asp, Scales)], variable.name ='traitName', value.var = 'traitValue')
coll_traits_info = add_info(coll_traits_melt, traitRef, traitDataID, traitDescription, traitUnits, c('27406 synthesised in 27707'))

fwrite(coll_traits_info, "Data/Temporary_data/Coll_traits.csv")


# ************************** #
#### 2. Species-level PCA ####
# ************************** #
pca_coll_sp = dudi.pca(mice::complete(mice(coll_traits[, c( 'logSize', 'Gen_per_year', 'Depth_preference', 'Repro_sex')])), scannf = FALSE, nf = 2)
gg_coll_sp = fviz_pca(pca_coll_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                      col.ind = "steelblue",
                      fill.ind = "white",
                      col.var = "black")
ggsave(gg_coll_sp, file = 'Results/species_pca_coll.pdf', width = 6, height = 5)


pca_species = dudi.pca(coll_traits[complete.cases(coll_traits[,c( 'logSize', 'Gen_per_year', 'Depth_preference', 'Repro_sex')]),c( 'logSize', 'Gen_per_year', 'Depth_preference', 'Repro_sex')], scannf = FALSE, nf = 3)
pca_coll_species= fviz_pca_biplot(pca_species, geom = c("point"), repel = T, axes = c(1,2), title = 'Collembola')
ggsave(pca_coll_species,file= 'Results/Species_PCA_Collembola.pdf', width = 5, height = 5)    

# Plot AEG43: only one individual, undefined, was found so we remove the plot
Abundances_coll = Abundance_all[Group_broad =='Collembola' & Plot != 'AEG43',]

# Check number of species for each trait
coll_traits[Species %in% Abundances_coll$Species , lapply(.SD, function(x){length(x[!is.na(x) & x != 'NA'])})]


# ********************************** #
#### 3. Community-weighted traits ####
# ********************************** #
CC_coll = check_coverage(coll_traits, Abundances_coll, c( 'logSize', 'Gen_per_year', 'Depth_preference', 'Repro_sex', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'#, 'Vertical_preference'
), 'Species', 'Species')
CWM_coll = my_cwm(coll_traits, Abundances_coll, c( 'logSize', 'Gen_per_year', 'Depth_preference', 'Repro_sex', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'#, 'Vertical_preference'
                                                       ), 'Species', 'Species')

# Some traits are categorical  so we change the name for merging (e.g. Repro_sex_1 means proportion with sexual repro, Repro_sex_0 means proportion with no sexual repro, we keep Repro_sex as Repro_sex_1)
CWM_coll[, c('Repro_sex','Pigment','Furca','PAO','PSO', 'Asp', 'Scales') := 
            list(Repro_sex_1, Pigment_1,Furca_1,PAO_1,PSO_1, Asp_1, Scales_1)] 

# Melt and merge
CWM_CC_coll= merge.data.table(melt.data.table(CWM_coll[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                melt.data.table(CC_coll[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


# Add info
traitUnits = c("mm (log-transformed)","unitless","unitless","unitless","%", '%', '%','%', '%', '%', '%')
traitDescription = c("Body size",
                     "Number of generations per year (coded as bivoltine = 2, multivoltine = 3",
                     "Depth preference (coded as Euedaphic = 2, hemiedaphic = 1, Epiedaphic=0)",
                     "Number of ocelli",
                     "% of sexual reproduction",
                     "% displaying pigments",
                     "% displaying a furca",
                     "% displaying Post Antennal Organ",
                     "% displaying pseudocelli",
                     "Presence of anal spines",
                     "% displaying scales")
traitDataID = c("NA","NA","NA","NA", "NA",'NA','NA','NA','NA','NA', 'NA')

traitRef = c("Saifutdinov, Zaytsev, John, Baulechner & Wolters","Saifutdinov, John, Baulechner & Wolters","Saifutdinov, John, Baulechner & Wolters","Saifutdinov, John, Baulechner & Wolters", "Polierer, Scheu, 2016; Chernova et al., 2010; Saifutdinov et a., 2018",'Saifutdinov, John, Baulechner & Wolters','Saifutdinov, John, Baulechner & Wolters','Saifutdinov, John, Baulechner & Wolters','Saifutdinov, John, Baulechner & Wolters','Saifutdinov, John, Baulechner & Wolters', 'Saifutdinov, John, Baulechner & Wolters')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c('logSize', 'Gen_per_year', 'Depth_preference', 'Ocelli', 'Repro_sex','Pigment','Furca', 'PAO','PSO','Asp', 'Scales')

CWM_CC_coll = add_info(CWM_CC_coll, traitRef, traitDataID, traitDescription, traitUnits, c('27406 synthesised in 27707'))

fwrite(CWM_CC_coll, "Data/CWM_data/CWM_coll.csv")


# ************************************** #
#### 4. Non-weighted community traits ####
# ************************************** #
Abundances_coll_presence_absence = Abundances_coll[, list(value = sum(value, na.rm = T), Year = 2019), by = list(Plot, Species)]
Abundances_coll_presence_absence[value>1, value := 1]

CC_coll_noweight = check_coverage(coll_traits, Abundances_coll_presence_absence, c( 'logSize', 'Gen_per_year', 'Depth_preference', 'Repro_sex', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'#, 'Vertical_preference'
), 'Species', 'Species')
CWM_coll_noweight = my_cwm(coll_traits, Abundances_coll_presence_absence, c( 'logSize', 'Gen_per_year', 'Depth_preference', 'Repro_sex', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'#, 'Vertical_preference'
), 'Species', 'Species')

# Some traits are categorical  so we change the name for merging (e.g. Repro_sex_1 means proportion with sexual repro, Repro_sex_0 means proportion with no sexual repro, we keep Repro_sex as Repro_sex_1)
CWM_coll_noweight[, c('Repro_sex','Pigment','Furca','PAO','PSO', 'Asp', 'Scales') := 
           list(Repro_sex_1, Pigment_1,Furca_1,PAO_1,PSO_1, Asp_1, Scales_1)] 

# Melt and merge
CWM_CC_coll_noweight= merge.data.table(melt.data.table(CWM_coll_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                              melt.data.table(CC_coll_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))

CWM_CC_coll_noweight = add_info(CWM_CC_coll_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('27406 synthesised in 27707'))

fwrite(CWM_CC_coll_noweight, "Data/CWM_data/CWM_coll_noweight.csv")

# ***************************** #
#### 5. Turnover (Table S6) #### 
# ***************************** #
comm.test_coll = dcast(Abundances_coll,  Plot~Species, value.var = 'value', fill = 0)
rownames(comm.test_coll)= comm.test_coll$Plot
comm.test_coll = comm.test_coll[,-1]

beta.multi.abund(comm.test_coll)

comm_coll_min_max = matrix(c(colSums(comm.test_coll[min_lui_plots,]),colSums(comm.test_coll[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_coll_min_max)


########################
#### ORIBATID MITES ####
########################

# ******************* #
#### 1. Load data ####
# ****************** #

mites_traits = fread('Data/Trait_data/2019_Mite_fauna_Traits.csv')
mites_traits[, Species := V1]
Abundances_mites = Abundance_all[Group_fine=='Oribatida', ]
Abundances_mites = Abundances_mites[Species != 'Gamasidae',] # These are not Oribatids

# Correct some species names manually
mites_traits[, Species := dplyr::recode(Species,
                                 'Ceratozetella_sellnicki' = 'Ceratozetes_sellnicki',
                                 'Nanhermannia_nanus' = 'Nanhermannia_nana',
                                 'Phthiracarus_globosus' = 'Phthiracarus_globulus',
                                 'Poecilochthonius_italicus' = 'Poeliochthonius_italicus',
                                 'Poecilochthonius_spiciger' = 'Poeliochthonius_spiciger',
                                 'Scutovertex_sp.' = 'Scutovertex_sp'
                                 )]

mites_traits[, Habitat_spec := as.numeric(dplyr::recode(Vertical_distribution.Morphotype,
                                     "non-spec" = '1',
                                     "soil" = '2',
                                     "surface" = '3',
                                     "litter" = '4'))]
mites_traits[, Feeding_spec := as.numeric(dplyr::recode(Feeding,
                                                "omnivorous" = '1',
                                                "herbifungivorous" = '2',
                                                "herbivorous" = '3',
                                                "fungivorous" = '4'))]
mites_traits[, Repro_sex := as.numeric(ifelse(Reproduction == 'parthenogenetic', 0,
                                   ifelse( Reproduction == 'sex', 1,
                                           NA)))]

mites_traits[, logMass := log(Mass)]


# For species id to genus level, we take the average or value of other species from the same genus if it's consistent (mean > 2sd)
# Carabodes_sp: only one species
mites_traits = rbind(mites_traits, mites_traits[grepl('Carabodes', Species), list(V1 = 'Carabodes_sp', Vertical_distribution.Morphotype, Surface,Litter,  Soil, Feeding, Reproduction, Mass, size100_500, 
                                                                                  size500_900, size900,  ash,  CuticuleCalcium, CuticuleNo, DaystoAdult, Species = 'Carabodes_sp' , Habitat_spec,Feeding_spec, Repro_sex, logMass)])

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
                                                                                logMass= mean(logMass))])

#Pelotpulus_sp
mites_traits = rbind(mites_traits, mites_traits[grepl('Peloptulus', Species), list(V1 = 'Pelotpulus_sp', Vertical_distribution.Morphotype, Surface,Litter,  Soil, Feeding, Reproduction, Mass, size100_500, 
                                                                                   size500_900, size900,  ash,  CuticuleCalcium, CuticuleNo, DaystoAdult, Species = 'Pelotpulus_sp', Habitat_spec,Feeding_spec, Repro_sex, logMass)])

#Trichoribates_sp
mites_traits = rbind(mites_traits, mites_traits[grepl('Trichoribates', Species), list(V1 = 'Trichoribates_sp', Vertical_distribution.Morphotype, Surface,Litter,  Soil, Feeding, Reproduction, Mass, size100_500, 
                                                                                   size500_900, size900,  ash,  CuticuleCalcium, CuticuleNo, DaystoAdult, Species = 'Trichoribates_sp', Habitat_spec,Feeding_spec, Repro_sex, logMass)])

# Save species-matched trait data
traitUnits = c("mg","unitless","unitless","unitless","unitless","unitless","day")
traitDescription = c("Body mass",
                     "Feeding niche",
                     "Habitat",
                     "Feeding specialisation (coded based on expert knowledge as: omnivorous = 1, herbifungivorous = 2, herbivorous = 3, fungivorous = 4)",
                     "Habitat specialisation (coded based on expert knowledge as: non-specialist = 1, soil-dwelling = 2, surface-dwelling = 3, litter-dwelling = 4)",
                     "Reproduction: sexual or parthogenetic",
                     "Number of days to reach maturity")
traitDataID = c("NA","NA","NA","NA","NA","NA", "NA")

traitRef = c("Saifutdinov, Zaytsev, John, Baulechner & Wolters","Saifutdinov, Zaytsev, John, Baulechner & Wolters","Saifutdinov, Zaytsev, John, Baulechner & Wolters","Saifutdinov, Zaytsev, John, Baulechner & Wolters","Saifutdinov, Zaytsev, John, Baulechner & Wolters",
             "Saifutdinov, Zaytsev, John, Baulechner & Wolters","Krivolutsky, 1995; Bellido A. 1970. Canadian Journal of Zoology 68(10):2221-2229; Grishina, 1991; Luxton, 1981")

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c( 'Mass', 'Feeding', 'Vertical_distribution.Morphotype','Feeding_spec', 'Habitat_spec' , 'Repro_sex', 'DaystoAdult')

mites_traits_melt = melt.data.table(mites_traits[, list(Species, Mass, Feeding, Habitat = Vertical_distribution.Morphotype,Feeding_spec, Habitat_spec, Repro_sex, DaystoAdult)], id.var = 'Species', value.var = 'traitValue', variable.name = 'traitName')
mites_traits_info = add_info(mites_traits_melt, traitRef, traitDataID, traitDescription, traitUnits, c('27406 synthesised in 27707'))

fwrite(mites_traits_info, "Data/Temporary_data/Mites_traits.csv")


# ************************** #
#### 2. Species-level PCA ####
# ************************** #
pca_mites_sp = dudi.pca(mice::complete(mice(mites_traits[, c( 'logMass', 'Repro_sex', 'Feeding_spec', 'Habitat_spec')])), scannf = FALSE, nf = 2)
gg_mites_sp = fviz_pca(pca_mites_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                       col.ind = "steelblue",
                       fill.ind = "white",
                       col.var = "black")
ggsave(gg_mites_sp, file = '/Users/Margot/Desktop/Research/Senckenberg/Documents/Papers/Traits/Figures/species_pca_mites.pdf', width = 6, height = 5)


# ********************************** #
#### 3. Community-weighted traits ####
# ********************************** #
# Check species coverage
mites_traits[Species %in% Abundances_mites$Species, lapply(.SD, function(x){length(x[!is.na(x)])}),]


CC_Mites = check_coverage(mites_traits, Abundances_mites, c( 'logMass', 'Feeding_spec', 'Habitat_spec', 'Repro_sex', 'DaystoAdult'),
                       'Species', 'Species')

CWM_Mites = my_cwm(mites_traits, Abundances_mites, c( 'logMass', 'Feeding_spec', 'Habitat_spec', 'Repro_sex', 'DaystoAdult'),
                       'Species', 'Species')

# Remove some plots with too little data
CC_Mites[, average_coverage := mean(c(logMass,Feeding_spec, Habitat_spec, Repro_sex , DaystoAdult)), by = 1:nrow(CC_Mites)]
CWM_Mites = CWM_Mites[Plot %in% CC_Mites[average_coverage >= 0.5, Plot],]

# Repro_sex is a categorical variable so we change the name for merging
CWM_Mites[, Repro_sex := Repro_sex_1] # Repro_sex_1 means sexual reproduction

# Melt and merge
CWM_CC_mites = merge.data.table(melt.data.table(CWM_Mites[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                       melt.data.table(CC_Mites[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


# Add info
traitUnits = c("mg (log-transformed)","","","%","day")
traitDescription = c("Body mass",
                     "Feeding specialisation (coded based on expert knowledge as: omnivorous = 1, herbifungivorous = 2, herbivorous = 3, fungivorous = 4)",
                     "Habitat specialisation (coded based on expert knowledge as: non-specialist = 1, soil-dwelling = 2, surface-dwelling = 3, litter-dwelling = 4)",
                     "Proportion of species with sexual reproduction",
                     "Number of days to reach maturity")
traitDataID = c("","","","", "")

traitRef = c("Saifutdinov, Zaytsev, John, Baulechner & Wolters","Saifutdinov, Zaytsev, John, Baulechner & Wolters","Saifutdinov, Zaytsev, John, Baulechner & Wolters",
             "Saifutdinov, Zaytsev, John, Baulechner & Wolters","Krivolutsky, 1995; Bellido A. 1970. Canadian Journal of Zoology 68(10):2221-2229; Grishina, 1991; Luxton, 1981")

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c( 'logMass', 'Feeding_spec', 'Habitat_spec', 'Repro_sex', 'DaystoAdult')

CWM_CC_mites = add_info(CWM_CC_mites, traitRef, traitDataID, traitDescription, traitUnits, c('27406 synthesised in 27707'))

fwrite(CWM_CC_mites, "Data/CWM_data/CWM_mites.csv")


# ************************************** #
#### 4. Non-weighted community traits ####
# ************************************** #

Abundances_mites_presence_absence = Abundances_mites[, list(value = sum(value, na.rm = T), Year = 2019), by = list(Plot, Species)]
Abundances_mites_presence_absence[value>1, value := 1]

CC_Mites_noweight = check_coverage(mites_traits, Abundances_mites_presence_absence, c( 'logMass', 'Feeding_spec', 'Habitat_spec', 'Repro_sex', 'DaystoAdult'),
                          'Species', 'Species')

CWM_Mites_noweight = my_cwm(mites_traits, Abundances_mites_presence_absence, c( 'logMass', 'Feeding_spec', 'Habitat_spec', 'Repro_sex', 'DaystoAdult'),
                   'Species', 'Species')

# Remove some plots with too little data
CC_Mites_noweight[, average_coverage := mean(c(logMass,Feeding_spec, Habitat_spec, Repro_sex , DaystoAdult)), by = 1:nrow(CC_Mites_noweight)]
CWM_Mites_noweight = CWM_Mites_noweight[Plot %in% CC_Mites_noweight[average_coverage >= 0.2, Plot],]

# Repro_sex is a categorical variable so we change the name for merging
CWM_Mites_noweight[, Repro_sex := Repro_sex_1] # Repro_sex_1 means sexual reproduction

# Melt and merge
CWM_CC_mites_noweight = merge.data.table(melt.data.table(CWM_Mites_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                melt.data.table(CC_Mites_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


# Add info
CWM_CC_mites_noweight = add_info(CWM_CC_mites_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('27406 synthesised in 27707'))

fwrite(CWM_CC_mites_noweight, "Data/CWM_data/CWM_mites_noweight.csv")

# ***************************** #
#### 5. Turnover (Table S6) #### 
# ***************************** #

comm.test_mites = dcast(Abundances_mites,  Plot~Species, value.var = 'value', fill = 0)
rownames(comm.test_mites)= comm.test_mites$Plot
comm.test_mites = comm.test_mites[,-1]

beta.multi.abund(comm.test_mites)

comm_mites_min_max = matrix(c(colSums(comm.test_mites[min_lui_plots,]),colSums(comm.test_mites[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_mites_min_max)


