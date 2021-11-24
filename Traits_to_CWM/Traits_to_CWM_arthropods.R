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

####  Compare trait datasets #### 
#published in Goessner 2015
Arthropod_traits_Goessner = setDT(read_delim("Traits/Arthropods/ArthropodSpeciesTraits.csv", ";", escape_double = FALSE, trim_ws = TRUE))
Arthropod_traits_Goessner[, new_scientificName := get_gbif(SpeciesID), by = SpeciesID]

## New version updated by Nadja Simons
Arthropod_traits_grasslands = fread('Traits/Arthropods/NEW_Arthropod_traits_0ct2020/Arthropod_traits_grasslands_std.csv')
Arthropod_traits_grasslands[scientificName == 'Sphaeroderma rubidum', order := "Coleoptera"]
Arthropod_species_info_grasslands = fread('Traits/Arthropods//NEW_Arthropod_traits_0ct2020/taxa_grasslands.csv')
Arthropod_species_info_grasslands[scientificName == 'Sphaeroderma rubidum', order := "Coleoptera"]

Arthropod_traits_grasslands[, new_scientificName := get_gbif(scientificName), by = scientificName]

## Additional traits



# Traits for herbivores
Neff2020_traits = data.table(read.csv("~/Desktop/Research/Senckenberg/Data/Traits/Arthropods/Neff_et_al_2020.Ecological_Applications.csv"))
Neff2020_traits[, new_scientificName := get_gbif(Species), by = Species]

# Traits from Birkhofer et al. 2017 (Araneae)
Birkhofer_araneae = data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Arthropods/Birkhofer_2017/Traits_Communities.xlsx", sheet = 'Araneae'))
Birkhofer_araneae[, new_scientificName := get_gbif(gsub('\\.', ' ', Site)), by = Site]
Birkhofer_araneae$new_scientificName %in% Arthropod_traits_grasslands$new_scientificName


####### Abundances ########

Abundances_all = fread("Abundances/Dataset_clean.txt")
Abundances_all[, length(unique(Plot)),by = c('Group_broad', 'Year') ]
Abundances_2008 = Abundances_all[Group_broad %in% c("Araneae" , "Coleoptera", "Hemiptera" ,"Orthoptera"),] 

Abundances_temporal = fread('Abundances/21969_4_Dataset/21969_4_data.csv')
Abundances_temporal = Abundances_temporal[!is.na(Species),]
Abundances_temporal[, Plot := ifelse(nchar(PlotID) == 5, PlotID, paste(substr(PlotID, 1, 3), '0', substr(PlotID, 4, 4), sep = ''))]
Abundances_temporal_by_plot = Abundances_temporal[, list(value = sum(NumberAdults, na.rm = T)), by = list('Plot' = Plot, 'Species' = Species, 'Year' = CollectionYear, 'Order' = Order)]
Abundances_temporal_by_plot[, new_scientificName := get_gbif(Species), by = Species]


# Species names to check
#[1] "Acanephodus_onopordi"       "Acanthodelphax_spinosa"     "Aphrodes_bicincta"          "Asiorestia_ferruginea"     
#[5] "Asiorestia_transversa"      "Furcipus_rectirostris"      "Hyledelphax_elegantula"     "Hypera_zoilus"             
#[9] "Kybos_smaragdula"           "Laodelphax_striatella"      "Megadelphax_sordidula"      "Polydrusus_sericeus"       
#[13] "Ribautodelphax_albostriata" "Ribautodelphax_collina"     "Xanthodelphax_flaveola"    "Xanthodelphax_straminea"  
#
#Abundances_temporal[, sum(NumberAdults[Species %in% Arthropod_traits_grasslands$scientificName]/sum(NumberAdults, na.rm = T), na.rm = T), by = PlotID][, median(V1)]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Recode traits ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Make sure 1 name = 1 trait value
Arthropod_traits_grasslands[, verbatimTraitValue := gsub(',', '.', verbatimTraitValue) ]
Arthropod_traits_grasslands_unique = Arthropod_traits_grasslands[, list(verbatimTraitValue = ifelse(traitName %in% c('Dispersal_ability', 'Mean_BodySize'), 
                                                                                                as.character(mean(as.numeric(verbatimTraitValue), na.rm = T)),
                                                                                                paste(unique(verbatimTraitValue), sep = '-'))), by = list(traitName, new_scientificName, family, order)]

Arthropod_traits_cast = dcast.data.table(Arthropod_traits_grasslands_unique, new_scientificName + family + order ~ traitName, value.var = 'verbatimTraitValue')
Arthropod_traits_cast$Dispersal_ability = as.numeric(gsub(',', '.', Arthropod_traits_cast$Dispersal_ability))
Arthropod_traits_cast[, Stratum_use_simple := dplyr::recode(Stratum_use, "g"          = 's',
                                                            "g-(h)"          = 's',
                                                            "g-t"            = 's',
                                                            "g-h"            = 's',
                                                            "g-h-(t)"        = 's',
                                                            "s-t"            = 's', ###?
                                                            "s-h"            = 's', ###?
                                                            "g-(h)-(t)"      = 's',
                                                            "s"              = 's',
                                                            "g-h-t"         = 's', ##?
                                                            "(s)-g"          = 's',
                                                            "(s)-g-(h)-(t)" = 's',
                                                            "s-g" = 's',
                                                                   "t"             = 't',
                                                                  "h"              = 'h',
                                                                  "i"              = 'NA',
                                                                  "h-t"           = 'h',
                                                                  "(h)-t"         = 't',
                                                                  "h-(t)"         = 'h',
                                                                  "(g)-h"          = 'h',
                                                               #   ""              = 'h',
                                                                  "(s)-(g)-h-(t)"  = 'h',
                                                                  "w-h"            = 'h',
                                                                  "(g)-h-(t)"     = 'h',
                                                                   "w"             = 'NA',
                                                                   "(g)-(h)-t"      = 't',
                                                                   "(s)-(g)-h"      = 'h'
                                                                   )] 
Arthropod_traits_cast[, Feeding_guild_simple := dplyr::recode(Feeding_guild, 
                                                              "h"        = 'primary',
                                                              "c"        = 'secondary',
                                                              "c-h"      = 'omnivore_secondary',
                                                              "c-d"      = 'omnivore_secondary',
                                                              "f"        = 'primary',
                                                              "d-h"      = 'omnivore_primary',
                                                              "d"        = 'primary',
                                                              "ce"       = 'secondary', # carnivore extraintestinal?
                                                              "hsp"      = 'primary',# Herbibore, sucking?
                                                              "h-(c)"    = 'primary',
                                                              "d-f"      = 'primary',
                                                              "h-hcf"    = 'primary',# Herbibore, feeding on buds and flowers?
                                                               "c-d-h"   = 'omnivore_secondary',
                                                               "hs"      = 'primary', #Herbivore, ??
                                                               "f-h"     = 'omnivore_primary',
                                                               "c-(h)"   = 'omnivore_secondary',
                                                               "hss"     = 'primary', #Herbivore, seed-feeding
                                                               "m"       = 'omnivore_primary', # Detritivore, mold, fungi
                                                               "c-h-hc"  = 'omnivore_primary', # Ophonus and amara, mostly herbivore according to Dennis?
                                                               "cs-hs"  = 'omnivore_secondary',#?
                                                               "dx-hcf"  = 'omnivore_primary', # Detritivore, mold, fungi not sure if other insects
                                                               "hss-(hs)" = 'primary', #Herbivore, ??
                                                               "c-m"    = 'omnivore_primary',# cambivore??
                                                               "cs"      = 'omnivore_secondary',#?
                                                               "hsm"    = 'primary',
                                                               "c-d-dd"  = 'omnivore_primary',#mostly detritivore (dung) but can feed on small arthropods
                                                               "c-f"     = 'omnivore_secondary',
                                                               "c-hcf"   = 'omnivore_secondary',
                                                               "ca"      = 'omnivore_secondary',
                                                               "d-do"    = 'omnivore_primary',
                                                               "hs-(cs)" = 'omnivore_primary',
                                                               "ca-m"    = 'omnivore_secondary',
                                                               "cs-(hs)"  = 'omnivore_secondary',
                                                               "cs-hsf-(hs)" = 'omnivore_secondary',
                                                               "cx-dx-hcf"  = 'omnivore_secondary',
                                                               "d-m"     = 'omnivore_primary',
                                                               "dc-do-h" = 'omnivore_secondary',
                                                               "do"       = 'omnivore_primary',#detritovire/dung
                                                               "h_r"     = 'primary',
                                                               "h-hc"   = 'primary',#???
                                                               "h-hcf-m" = 'primary',#???
                                                               "hcf"      = 'primary',#???
                                                               "hcf-m"    = 'omnivore_secondary',# maybe???
                                                               "hsf"     = 'primary',#???
                                                               "hsf-(cs)"  = 'primary',#???
                                                               "hsf-(hs)"  = 'primary'#???
                                                              )]
 
Arthropod_traits_cast[, Stratum_use_numeric := mgsub(Stratum_use, c('s', 'g', 'h', 't', 'w', 'i'), c(0, 1, 2, 3, NA, NA))] 
Arthropod_traits_cast[, Stratum_use_simple_numeric := sapply(Stratum_use_numeric, function(x){
  xlist = unlist(strsplit(x, '-'))
  xlist_numeric =  as.numeric(gsub('[\\(\\)]', '', xlist))
  weights = ifelse(grepl('[\\(\\)]', xlist), 0.5, 1)
  return(weighted.mean(xlist_numeric, weights))
})] 

# Also using the old dataset
Arthropod_traits_Goessner$Feeding_generalism <- c(1,2,3)[match(Arthropod_traits_Goessner$Feeding_specialization, c("m", "o", "p"))]

Arthropod_traits = merge(Arthropod_traits_cast, Arthropod_traits_Goessner[, list( Body_SizeG = Body_Size, Feeding_generalism, Feeding_mode, new_scientificName)], by = 'new_scientificName', all = T)
Arthropod_traits$Mean_BodySize = as.numeric(Arthropod_traits$Mean_BodySize )

# Add Neff 2020
Arthropod_traits = merge.data.table(Arthropod_traits, Neff2020_traits[, list( Generations, Hibernation, new_scientificName)], by = 'new_scientificName', all = T)
Arthropod_traits[, Hibernation_late := as.numeric(recode(Hibernation, 'e' = 1,
                                                                            "e/n" = 1.5,
                                                                            "n" = 2,
                                                                            "a" = 3 ))]


# Fill-in Aranae hunting traits from Birkhofer et al. (2017).
Arthropod_traits = merge.data.table(Arthropod_traits, Birkhofer_araneae[, list(Active_hunt = weighted.mean(c(1, 0), c(cursorial.hunting, web.hunting))), by = new_scientificName], by = 'new_scientificName', all = T)
Arthropod_traits$logBody_Size = log(Arthropod_traits$Mean_BodySize)




###### Species names
#Abundances$species_check = Abundances$Species
#Arthropod_traits[, Chewer := Feeding_mode == 'c']
#Arthropod_traits[, Extrainst := Feeding_mode == 'e']
#Arthropod_traits[, Suckers := Feeding_mode == 's']

###### Check trait distribution
Arthropod_traits[, sapply(.SD, function(x){hist(x)}), .SDcols = c("Mean_BodySize", "Dispersal_ability", "Feeding_generalism")]


### Calculate CWM
### All together
Arthropod_traits[, Primary := Feeding_guild_simple == 'primary',]
traits_all =Arthropod_traits[new_scientificName %in% Abundances_temporal_by_plot$new_scientificName, .SD, .SDcols = c('Feeding_generalism', 'Dispersal_ability', 'logBody_Size')]
pca_species = dudi.pca(traits_all[complete.cases(traits_all),], scannf = FALSE, nf = 2)
fviz_pca(pca_species)

CWM_Arthropods_all = my_cwm(Arthropod_traits, 
                                   Abundances_temporal_by_plot, c('Feeding_generalism', 'Dispersal_ability', 'logBody_Size', 'Primary'), 'Species', 'Species')
test_pca = dudi.pca(complete(mice(CWM_Arthropods_all[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size', 'Feeding_generalism', 'Primary_TRUE')][, -1])), scannf = FALSE, nf = 2)
fviz_pca(test_pca)

test = merge(CWM_Arthropods_all[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size', 'Feeding_generalism', 'Primary_TRUE')],
             env_data_lui, by = 'Plot')
cor.test(test_pca$l1$RS1, env_data_lui$LUI)
cor.test(test$LUI, test$Dispersal_ability)
cor.test(test$LUI, test$Feeding_generalism)
cor.test(test$LUI, test$logBody_Size)


# Herbivore, AG
traits_h_AG = Arthropod_traits[Stratum_use_simple != 's'  & Feeding_guild_simple %in% c('omnivore_primary', 'primary') &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName, .SD, 
                               .SDcols = c('Feeding_generalism', 'Dispersal_ability', 'logBody_Size', 'Generations')]
cor.test(traits_h_AG$logBody_Size, traits_h_AG$Feeding_generalism)
pca_species = dudi.pca(traits_h_AG[!is.na(Generations),], scannf = FALSE, nf = 2)
fviz_pca(pca_species)

Ab_h_AG = Abundances_temporal_by_plot[new_scientificName %in% Arthropod_traits[Stratum_use_simple != 's'  & Feeding_guild_simple %in% c('omnivore_primary', 'primary'),new_scientificName], .N, by = c('Plot', 'Year')]
CWM_Arthropods_above_herb = my_cwm(Arthropod_traits[Stratum_use_simple != 's'  & Feeding_guild_simple %in% c('omnivore_primary', 'primary'),], 
                                   Abundances_temporal_by_plot, c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Feeding_mode', 'Generations'), 'new_scientificName', 'new_scientificName')
test_pca = dudi.pca(complete(mice(CWM_Arthropods_above_herb[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size',  'Generations')][, -1])), scannf = FALSE, nf = 2)
fviz_pca(test_pca)

test = merge(CWM_Arthropods_above_herb[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size', 'Feeding_generalism', 'Generations', 'Hibernation_late')],
             env_data_lui, by = 'Plot')
cor.test(test_pca$l1$RS1, env_data_lui$LUI)
cor.test(test$LUI, test$Dispersal_ability)
cor.test(test$LUI, test$Feeding_generalism)
cor.test(test$LUI, test$logBody_Size)

# Carnivores, AG ---> Need to take axis 2
trait_c_AG = Arthropod_traits[Stratum_use_simple != 's' & Feeding_guild_simple  %in% c('omnivore_secondary', 'secondary') &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName,lapply(.SD, as.numeric), .SDcols = c('Dispersal_ability', 'logBody_Size', 'Active_hunt')]
pca_c_AG = dudi.pca(trait_c_AG[complete.cases(trait_c_AG),], scannf = FALSE, nf = 2)
fviz_pca(pca_c_AG)
cor.test(trait_c_AG$Dispersal_ability, trait_c_AG$Active_hunt)

CWM_Arthropods_above_omni_carni = my_cwm(Arthropod_traits[Stratum_use_simple != 's' & Feeding_guild_simple  %in% c('omnivore_secondary', 'secondary'),], 
                                         Abundances_temporal_by_plot, c('Dispersal_ability',"logBody_Size"), 'new_scientificName', 'new_scientificName')
pca_c_AG = dudi.pca(complete(mice(CWM_Arthropods_above_omni_carni[, lapply(.SD, mean, na.rm = T), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size')][, -1])), , scannf = FALSE, nf = 2)
fviz_pca(pca_c_AG)
test = merge(CWM_Arthropods_above_omni_carni[, lapply(.SD, mean, na.rm = T), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size')],
             env_data_lui, by = 'Plot')

cor.test(test$LUI, test$Dispersal_ability)
cor.test(test$LUI, pca_c_AG$l1$RS2)
cor.test(test$LUI, scale(test$Dispersal_ability) +scale(-test$logBody_Size))


# Herbivores, BG
trait_h_BG = Arthropod_traits[Stratum_use_simple == 's' & Feeding_guild_simple  %in% c('omnivore_primary', 'primary') &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName,lapply(.SD, as.numeric), 
                              .SDcols = c('Dispersal_ability', 'logBody_Size','Feeding_generalism')]
pca_h_BG = dudi.pca(trait_h_BG[complete.cases(trait_h_BG),], scannf = FALSE, nf = 2)
fviz_pca(pca_h_BG)
cor.test(trait_h_BG$logBody_Size, trait_h_BG$Dispersal_ability)



CC_Arthropods_below_herb = check_coverage(Arthropod_traits[Stratum_use_simple == 's'  & Feeding_guild_simple %in% c('omnivore_primary', 'primary'),], 
                                          Abundances_temporal_by_plot[new_scientificName %in% Arthropod_traits[Stratum_use_simple == 's'  & Feeding_guild_simple %in% c('omnivore_primary', 'primary'),new_scientificName],], 
       c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Feeding_mode'), 'new_scientificName', 'new_scientificName')

CWM_Arthropods_below_herb = my_cwm(Arthropod_traits[Stratum_use_simple == 's'  & Feeding_guild_simple %in% c('omnivore_primary', 'primary'),], Abundances_temporal_by_plot, 
                                   c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Feeding_mode'), 'new_scientificName', 'new_scientificName')
test_pca = dudi.pca(mice::complete(mice(CWM_Arthropods_below_herb[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size')][, -1])), scannf = FALSE, nf = 2)
fviz_pca(test_pca)
test = merge(CWM_Arthropods_below_herb[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size', 'Feeding_generalism')],
             env_data_lui, by = 'Plot')
cor.test(test$Dispersal_ability, test$logBody_Size)
cor.test(test$LUI, c(scale(test$logBody_Size)))
cor.test(test[Plot %in% CWM_Arthropods_below_herb$Plot,]$LUI, test_pca$l1$RS1)

# Carnivores, BG
trait_c_BG = Arthropod_traits[Stratum_use_simple == 's' & Feeding_guild_simple  %in% c('omnivore_secondary', 'secondary') &
                                new_scientificName %in% Abundances_temporal_by_plot$new_scientificName,lapply(.SD, as.numeric), .SDcols = c('Dispersal_ability', 'logBody_Size','Feeding_generalism')]
pca_c_BG = dudi.pca(trait_c_BG[complete.cases(trait_c_BG),], scannf = FALSE, nf = 2)
fviz_pca(pca_c_BG)
cor.test(trait_c_BG$logBody_Size, trait_c_BG$Dispersal_ability)

CWM_Arthropods_below_omni_carni = my_cwm(Arthropod_traits[Stratum_use_simple == 's' & Feeding_guild_simple  %in% c('omnivore_secondary', 'secondary'),], Abundances_temporal_by_plot, c( 'Dispersal_ability',"logBody_Size", 'Feeding_generalism', 'Feeding_mode'), 'new_scientificName', 'new_scientificName')
test_pca = dudi.pca(mice::complete(mice(CWM_Arthropods_below_omni_carni[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size')][, -1])), scannf = FALSE, nf = 2)
fviz_pca(test_pca)
test = merge(CWM_Arthropods_below_omni_carni[, lapply(.SD, mean), by = Plot, .SDcols = c('Dispersal_ability', 'logBody_Size')],
             env_data_lui, by = 'Plot')

cor.test(test$LUI, test$Dispersal_ability)
cor.test(test$LUI, test_pca$l1$RS1)
cor.test(test$Dispersal_ability, test$logBody_Size)
cor.test(test$LUI, -scale01(test$logBody_Size) + scale01( test$Dispersal_ability))

# Save datasets
write.csv(CWM_Arthropods_above_herb[,list( "Ah_Dispersal" = Dispersal_ability, 
                                         "Ah_BodySize" = logBody_Size,
                                        # "Ah_StratumHerb" = Stratum_use_simple_h/(Stratum_use_simple_h +Stratum_use_simple_NA +Stratum_use_simple_t),
                                         "Ah_Generalism" = Feeding_generalism,
                                         "Ah_Generations" = Generations,
                                      #   "Ah_chewers" = Feeding_mode_c,
                                        # "Ah_suckers" = Feeding_mode_s,
                                         #"Ah_abundance" = arthro_herb_abundance,
                                         #"Ah_total_mass" = arthro_herb_abundance * Mean_BodySize,
                                         "Plot" = Plot,
                                         "Year" = Year )], 

          "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_above_herb.csv")
 
write.csv(CWM_Arthropods_below_herb[,list( "Ah_b_Dispersal" = Dispersal_ability, 
                                           "Ah_b_BodySize" = logBody_Size,
                                           # "Ah_StratumHerb" = Stratum_use_simple_h/(Stratum_use_simple_h +Stratum_use_simple_NA +Stratum_use_simple_t),
                                           "Ah_b_Generalism" = Feeding_generalism,
                                        #   "Ah_Generations" = Generations,
                                           #   "Ah_chewers" = Feeding_mode_c,
                                           # "Ah_suckers" = Feeding_mode_s,
                                           #"Ah_abundance" = arthro_herb_abundance,
                                           #"Ah_total_mass" = arthro_herb_abundance * Mean_BodySize,
                                           "Plot" = Plot,
                                           "Year" = Year )], 
          
          "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_below_herb.csv")


write.csv(CWM_Arthropods_above_omni_carni[,list( "Aoc_Dispersal" = Dispersal_ability, 
                                           "Aoc_BodySize" = logBody_Size,
                                     #      "Aoc_StratumHerb" = Stratum_use_simple_h/(Stratum_use_simple_h  +Stratum_use_simple_t),
                                          # "Aoc_generalism" = Feeding_generalism,
                                       #    "Aoc_chewers" = Feeding_mode_c,
                                      #     "Aoc_suckers" =  Feeding_mode_s,
                                      #     "Aoc_extraint" = Feeding_mode_e,
                                           #"Ah_abundance" = arthro_herb_abundance,
                                           #"Ah_total_mass" = arthro_herb_abundance * Mean_BodySize,
                                           "Plot" = Plot,
                                           "Year" = Year )], 
          
          "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_above_omni_carni.csv")



write.csv(CWM_Arthropods_below_omni_carni[,list(
                                            "Aoc_b_Dispersal" = Dispersal_ability, 
                                            "Aoc_b_BodySize" = logBody_Size,
                                        #    "Aoc_b_extraint" = Feeding_mode_e,
                                            #"Aoc_b_chewers" = Feeding_mode_c,
                                            #"Aoc_b_suckers" = Feeding_mode_s,
                                            # "Aoc_b_generalism" = Feeding_generalism,
                                            # "Ac_abundance" = arthro_carni_abundance,
                                            # "Ac_total_mass" = arthro_carni_abundance * Mean_BodySize,
                                            "Plot" = Plot,
                                            "Year" = Year )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Arthropods_below_omni_carni.csv")

#write.csv(CWM_Arthropods_above_carni, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_above_carni.csv")
#write.csv(CWM_Arthropods_below_herb, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Arthropods_below_herb.csv")
#write.csv(CWM_Arthropods_below_carni, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_below_carni.csv")

#write.table(FD_Arthropods, file = "Data/CWM_data/FD_Arthropods.csv")
#write.csv(FD_Arthropods_ground, file = "Data/CWM_data/FD_Arthropods_ground.csv")
#write.csv(FD_Arthropods_above, file = "Data/CWM_data/FD_Arthropods_above.csv")


######  !!! MORRRRRE TRRAAAAAITS #####
#pollinators = fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/FlosDocs/data/pollinator_traits/18987.txt")

#### Heteroptera ####
# # ==> Included in the Big Dataset
# heteroptera_traits = data.table(read_excel("Traits/Arthropods/arthropod_traits_c_westphal/Heteroptera\ Traits\ -\ Sagrario\ Gwen\ Christoph.xlsx"))
# 
# # Check distribution
# heteroptera_traits[, lapply(.SD, hist), .SDcols = c( "Size.med")]
# heteroptera_traits[, c(  "logSize"):= lapply(.SD, log), .SDcols = c(  "Size.med")]
# 
# 
# heteroptera_traits$Species[!(heteroptera_traits$Species%in% Abundances_all$Species)]
# Abundances_all[Group_broad == 'Hemiptera',sort(unique(Species[(Species %in% heteroptera_traits$Species)]))]
# heteroptera_traits$Species[!(heteroptera_traits$Species %in% gsub(' ','_',Arthropod_traits_grasslands$scientificName))]
# #"Acalypta_marginata"
# #"Amblytylus_brevicollis" 
# #"Berytinus_hirticornis" 
# #"Berytinus_crassipes"      
# # "Capsus_pilifer"            
# #"Geocoris_grylloides"       
# #"Ischnocoris_hemipterus"    
# #"Kleidocerys_resedae"      
# # "Lygus_punctatus"           
# #"Nabis_ericetorum"          
# #"Peritrechus_geniculatus"   
# #"Pithanus_maerkelli"       
# # "Podops_inuncta"            
# #"Saldula_fucicola"          
# #"Saldula_orthochila" 
# # Taphropeltus_contractus"  
# #"Tropistethus_holosericeus"
# 
# Abundances_heteroptera = Abundances_all[Species %in% heteroptera_traits$Species,]
# heteroptera_traits[, c('Gen.per.year', 'Specialism') := list(as.numeric(dplyr::recode(Generations, 'Univoltine' = 1, Bivoltine = 2)),
#                                                              as.numeric(dplyr::recode(Feeding_niche...6, 'Poly' = 0, Oligo = 1, Mono = 2)))]
# 
# # Species PCA
# pca_heteroptera = dudi.pca(heteroptera_traits[, c('Gen.per.year', 'Specialism',  'Size.med','First.Month')], scannf = FALSE, nf = 2)
# fviz_pca_biplot(pca_heteroptera)
# 
# 
# heteroptera_cwm = my_cwm(heteroptera_traits, Abundances_heteroptera, c( 'Size.med', 'First.Month' ,'Gen.per.year', 'Specialism', 'ActivityRange', 'logSize'),
#                    'Species', 'Species')
# pca_heteroptera = dudi.pca(heteroptera_cwm[complete.cases(heteroptera_cwm), c('Gen.per.year', 'Specialism',  'Size.med','ActivityRange', 'First.Month')], scannf = FALSE, nf = 2)
# fviz_pca_biplot(pca_heteroptera)
# 
# write.csv(heteroptera_cwm, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_heteroptera.csv")



#### Flower visitors ####
# # ===> Mostly included in butterflies dataset

# flower_visitors_traits = data.table(read_excel("Traits/Arthropods/arthropod_traits_c_westphal/Flower\ visitor\ Traits\ -\ Gwen\ Heike.xlsx", sheet = 2))
# flower_visitors_traits[, Species := gsub(' ', '_', species) ]
# flower_visitors_traits[, size := as.numeric(size) ]
# 
# # Check distribution
# flower_visitors_traits[, lapply(.SD, function(x){hist(log(x))}), .SDcols = c( "size", "hibernation", "migration", "eggs", "phagie_H", "flight","density",  "eggmaturation", "distribution")]
# flower_visitors_traits[, c( "logSize", "logHibernation", "logMigration", "logEggs", "logPhagie_H", "logFlight","logDensity"):= lapply(.SD, log), .SDcols = c( "size", "hibernation", "migration", "eggs", "phagie_H", "flight","density")]
# 
# 
# Abundances_visitors = Abundances_all[Species %in% flower_visitors_traits$Species,]
# 
# flower_visitors_traits[,c( 'size', 'hibernation' ,'migration', 'eggs', 'phagie_H', 'flight', 'strategy', 'generation', 'eggmaturation') := lapply(.SD, as.numeric), 
#               .SDcols = c( 'size', 'hibernation' ,'migration', 'eggs', 'phagie_H', 'flight', 'strategy', 'generation', 'eggmaturation')]
# flower_visitors_traits_cwm = my_cwm(flower_visitors_traits, Abundances_visitors, c( "logSize", "logHibernation", "logMigration",  "logDensity", 'phagie_H', 'flight', 'strategy', 'generation', 'eggmaturation'),
#                          'Species', 'Species')
# pca_visitors= dudi.pca(flower_visitors_traits_cwm[complete.cases(flower_visitors_traits_cwm),-c('Plot', 'Year')], scannf = FALSE, nf = 2)
# fviz_pca_biplot(pca_visitors)
# write.csv(flower_visitors_traits_cwm, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_flower_visitors.csv")


### Collembola
coll_traits = data.table(read_excel(('Traits/Collembola_Mites/Coll_Traits_NEW.xlsx')))

#BE_species = intersect(coll_traits$Species, Abundances_all[Group_broad =='Collembola',Species])
#soilbiostore = fread('Traits/Collembola_Mites/SoilBioStoreSpecies.csv')
#soilbiostore$Species = gsub(' ', '_', soilbiostore$Species )
#
#coll_full = merge(coll_traits, soilbiostore[, .SD, .SDcols = c('Mode of reproduction'    , 'Phenology', 'Species')], by = 'Species', all.x = T)
##write_csv(coll_full, 'Traits/Collembola_Mites/Coll_Traits_2019+SoilBioStore.xls')
#no_data_spe = BE_species[!(BE_species %in% soilbiostore$Species)]
#no_data_gbif = data.table(get_gbif_taxonomy(no_data_spe))
#no_data_stand = no_data_gbif[, ifelse(is.na(scientificName),verbatimScientificName, gsub(' ', '_', scientificName) )]

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
                                     ifelse(Phenology == 'bivoltine', 2, NA))]
coll_traits[, Species := gsub(' ', '_', Species)]
pca_species = dudi.pca(coll_traits[complete.cases(coll_traits[,c( 'Size.log', 'Gen_per_year', 'Depth_preference', 'Repro_sex')]),c( 'Size.log', 'Gen_per_year', 'Depth_preference', 'Repro_sex')], scannf = FALSE, nf = 3)
pca_coll_species= fviz_pca_biplot(pca_species, geom = c("point"), repel = T, axes = c(1,2), title = 'Collembola')
ggsave(pca_coll_species,file= '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Species_PCA_Collembola.pdf', width = 5, height = 5)    

Abundances_coll = Abundances_all[Group_broad =='Collembola',]

Coll_cwm_all = my_cwm(coll_traits, Abundances_coll, c( 'Size.log', 'Gen_per_year', 'Depth_preference', 'Repro_sex', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'#, 'Vertical_preference'
                                                       ), 'Species', 'Species')

pca_coll = dudi.pca(complete(mice(Coll_cwm_all[ ,c( 'Size.log', 'Gen_per_year', 'Depth_preference', 'Repro_sex_1')])), scannf = FALSE, nf = 3)
fviz_pca_biplot(pca_coll, geom = c("point"), repel = T, axes = c(1,2), title = 'Collembola')

cor.test(pca_coll$l1$RS1, env_data_lui[Plot %in% Coll_cwm_all$Plot,]$LUI)


write.csv(Coll_cwm_all[, list( "col_Size" = Size.log,
                               "col_Depth" = Depth_preference,
                               "col_Gen_per_Year" = Gen_per_year,
                               "col_Sex" = Repro_sex,
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
mites_traits[!(Species %in% Abundances_mites$Species), unique(Species)]


Mites_cwm_all = my_cwm(mites_traits, Abundances_mites, c( 'Mass.log', 'Feeding_spec', 'Habitat_spec', 'Repro_sex', 'DaystoAdult'),
                       'Species', 'Species')
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
#### Mollusca ####

#### Acari ####

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



