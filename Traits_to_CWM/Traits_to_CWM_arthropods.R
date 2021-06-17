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

setwd("/Users/Margot/Desktop/Research/Senckenberg/Data/")

Arthropod_traits0 = setDT(read_delim("Traits/Arthropods/ArthropodSpeciesTraits.csv", 
                             ";", escape_double = FALSE, trim_ws = TRUE))
## New version
Arthropod_traits_grasslands = fread('Traits/Arthropods/NEW_Arthropod_traits_0ct2020/Arthropod_traits_grasslands_std.csv')
Arthropod_traits_grasslands[scientificName == 'Sphaeroderma rubidum', order := "Coleoptera"]
Arthropod_species_info_grasslands = fread('Traits/Arthropods//NEW_Arthropod_traits_0ct2020/taxa_grasslands.csv')
Arthropod_species_info_grasslands[scientificName == 'Sphaeroderma rubidum', order := "Coleoptera"]

Abundances_all = fread("Abundances/Dataset_clean.txt")
Abundances_all[, length(unique(Plot)),by = c('Group_broad', 'Year') ]
Abundances = Abundances_all[Group_broad %in% c("Araneae" , "Coleoptera", "Hemiptera" ,"Orthoptera"),] 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Recode traits ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

Arthropod_traits_cast = dcast.data.table(unique(Arthropod_traits_grasslands), verbatimScientificName ~ traitName, value.var = 'verbatimTraitValue')
Arthropod_traits_cast$Dispersal_ability = as.numeric(gsub(',', '.', Arthropod_traits_cast$Dispersal_ability))
Arthropod_traits_cast[, Stratum_use_simple := dplyr::recode(Stratum_use, "g"          = 's',
                                                                  "t"             = 't',
                                                                  "h"              = 'h',
                                                                  "i"              = 'NA',
                                                                  "h-t"           = 'h',
                                                                  "g-(h)"          = 's',
                                                                  "g-t"            = 's',
                                                                  "g-h"            = 'h',
                                                                  "(h)-t"         = 't',
                                                                  "h-(t)"         = 'h',
                                                                  "s-h"            = 'h',
                                                                  "(g)-h"          = 'h',
                                                               #   ""              = 'h',
                                                                  "g-h-(t)"        = 's',
                                                                  "s-t"            = 's',
                                                                  "(s)-(g)-h-(t)"  = 'h',
                                                                  "g-(h)-(t)"      = 's',
                                                                  "w-h"            = 'h',
                                                                  "s"              = 's',
                                                                  "(g)-h-(t)"     = 'h',
                                                                   "w"             = 'NA',
                                                                   "(g)-(h)-t"      = 't',
                                                                   "g-h-t"         = 'h',
                                                                   "(s)-g"          = 's',
                                                                   "(s)-g-(h)-(t)" = 's',
                                                                   "(s)-(g)-h"      = 'h',
                                                                   "s-g" = 's'
                                                                   )] 
Arthropod_traits_cast[, Stratum_use_numeric := mgsub(Stratum_use, c('s', 'g', 'h', 't', 'w', 'i'), c(0, 1, 2, 3, NA, NA))] 
Arthropod_traits_cast[, Stratum_use_simple_numeric := sapply(Stratum_use_numeric, function(x){
  xlist = unlist(strsplit(x, '-'))
  xlist_numeric =  as.numeric(gsub('[\\(\\)]', '', xlist))
  weights = ifelse(grepl('[\\(\\)]', xlist), 0.5, 1)
  return(weighted.mean(xlist_numeric, weights))
})] 
Arthropod_traits_cast$Species = Arthropod_traits_cast$verbatimScientificName

# Also using the old dataset
Arthropod_traits0$Feeding_generalism <- c(1,2,3)[match(Arthropod_traits0$Feeding_specialization, c("m", "o", "p"))]
Arthropod_traits0$trophic_level <- ifelse(Arthropod_traits0$Feeding_guild_short %in% c("c"), 'secondary_cons', 
                                          ifelse(Arthropod_traits0$Feeding_guild_short %in% c('o'),'omnivores', 'primary_cons'))
Arthropod_traits0$Species = Arthropod_traits0$SpeciesID

Arthropod_traits_cast2 = merge(Arthropod_traits_cast, Arthropod_traits0[, c('trophic_level', 'Feeding_generalism', 'Feeding_mode','Species')], by = 'Species', all = T)
Arthropod_traits_cast2$Mean_BodySize = as.numeric(Arthropod_traits_cast2$Mean_BodySize )
Arthropod_traits_cast2[, Species := gsub(' ', '_', Species)]

#### Check correspondance between 2 datasets ###
missing_trait_species = unique(Abundances$Species[!(Abundances$Species %in% Arthropod_traits_cast2$Species)])

length(missing_trait_species)/ length(unique(Abundances$Species))
# --> Only 2.7% species absent
Abundances[, sum(value[Species %in% missing_trait_species], na.rm = T)/ sum(value, na.rm = T)]
# --> correponding to 1% individuals

###### Species names
Abundances$species_check = Abundances$Species
Arthropod_traits_cast2$Species = gsub(' ', '_', Arthropod_traits_cast2$Species)
Arthropod_traits = Arthropod_traits_cast2

###### Check trait distribution
Arthropod_traits[, sapply(.SD, function(x){hist(x)}), .SDcols = c("Mean_BodySize", "Dispersal_ability", "Feeding_generalism")]
Arthropod_traits$logBody_Size = log(Arthropod_traits$Mean_BodySize)
traits <- c("logBody_Size", 'Mean_BodySize', "Dispersal_ability", "Feeding_generalist", 'Feeding_mode')#"Feeding_chewers", "Feeding_suckers")



arthro_herb_abundances = Abundances[Species %in% Arthropod_traits[trophic_level == 'primary_cons', Species], list(arthro_herb_abundance = sum(value, na.rm = T)), by = c('Plot')]
arthro_carni_abundances = Abundances[Species %in% Arthropod_traits[trophic_level == 'secondary_cons', Species], list(arthro_carni_abundance = sum(value, na.rm = T)), by = c('Plot')]

# Calcualte CWM
CWM_Arthropods_all_herb = my_cwm(Arthropod_traits[trophic_level == 'primary_cons',], Abundances, c( 'Dispersal_ability', 'Mean_BodySize','logBody_Size', 'Stratum_use','Stratum_use_simple', 'Stratum_use_simple_numeric', 'Feeding_generalism', 'Feeding_mode'),'Species', 'Species')
CWM_Arthropods_all_carni = my_cwm(Arthropod_traits[trophic_level == 'secondary_cons',], Abundances, c( 'Dispersal_ability', 'Mean_BodySize',"logBody_Size",'Stratum_use_simple','Stratum_use_simple_numeric', 'Feeding_generalism', 'Feeding_mode'),'Species', 'Species')

CWM_Arthropods_all_herb = merge.data.table(CWM_Arthropods_all_herb, arthro_herb_abundances)
CWM_Arthropods_all_carni = merge.data.table(CWM_Arthropods_all_carni, arthro_carni_abundances)

CWM_Arthropods_above_herb = my_cwm(Arthropod_traits[Stratum_use_simple != 's'  & trophic_level == 'primary_cons',], Abundances, c( 'Dispersal_ability', 'Mean_BodySize','Stratum_use_simple',"logBody_Size", 'Feeding_generalism', 'Feeding_mode'), 'Species', 'Species')
CWM_Arthropods_above_carni = my_cwm(Arthropod_traits[Stratum_use_simple != 's' & trophic_level == 'secondary_cons',], Abundances, c( 'Dispersal_ability', 'Mean_BodySize','Stratum_use_simple',"logBody_Size", 'Feeding_generalism', 'Feeding_mode'), 'Species', 'Species')
CWM_Arthropods_above_omni = my_cwm(Arthropod_traits[Stratum_use_simple != 's' & trophic_level == 'omnivores',], Abundances, c( 'Dispersal_ability','Stratum_use_simple', 'Mean_BodySize',"logBody_Size", 'Feeding_generalism', 'Feeding_mode'), 'Species', 'Species')

CWM_Arthropods_below_herb = my_cwm(Arthropod_traits[Stratum_use_simple == 's'  & trophic_level == 'primary_cons',], Abundances, c( 'Dispersal_ability', 'Mean_BodySize',"logBody_Size", 'Feeding_generalism', 'Feeding_mode'), 'Species', 'Species')
CWM_Arthropods_below_carni = my_cwm(Arthropod_traits[Stratum_use_simple == 's' & trophic_level == 'secondary_cons',], Abundances, c( 'Dispersal_ability', 'Mean_BodySize',"logBody_Size", 'Feeding_mode'), 'Species', 'Species')
CWM_Arthropods_below_omni = my_cwm(Arthropod_traits[Stratum_use_simple == 's' & trophic_level == 'omnivores',], Abundances, c( 'Dispersal_ability', 'Mean_BodySize',"logBody_Size", 'Feeding_generalism', 'Feeding_mode'), 'Species', 'Species')

# I try to add to the herbivore dataset the proportion of herbivores per year
Prop_herbivores = merge(Abundances[Species %in% Arthropod_traits[Stratum_use_simple !='g' & trophic_level == 'primary_cons', gsub(' ', '_', Species)], list(ab_herb = sum(value)), by = Plot],
  Abundances[Species %in%  gsub(' ', '_', Arthropod_traits$Species), list(ab_tot = sum(value)), by = Plot])
Prop_herbivores[, prop_herbivores := ab_herb/ab_tot]

Arthropod_traits_cast[Stratum_use_simple !='g', Species]
                      
CWM_Arthropods_above_herb = merge(CWM_Arthropods_above_herb, Prop_herbivores[, c('Plot', 'prop_herbivores')], by = 'Plot')

# Save datasets
write.csv(CWM_Arthropods_above_herb[,list( "Ah_Dispersal" = Dispersal_ability, 
                                         "Ah_BodySize" = logBody_Size,
                                         "Ah_StratumHerb" = Stratum_use_simple_h/(Stratum_use_simple_h +Stratum_use_simple_NA +Stratum_use_simple_t),
                                         "Ah_generalism" = Feeding_generalism,
                                         "Ah_chewers" = Feeding_mode_c,
                                         #"Ah_abundance" = arthro_herb_abundance,
                                         #"Ah_total_mass" = arthro_herb_abundance * Mean_BodySize,
                                         "Plot" = Plot,
                                         "Year" = Year )], 

          "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_above_herb.csv")
write.csv(CWM_Arthropods_above_omni[,list( "Ao_Dispersal" = Dispersal_ability, 
                                           "Ao_BodySize" = logBody_Size,
                                           "Ao_StratumHerb" = Stratum_use_simple_h/(Stratum_use_simple_h +Stratum_use_simple_NA  +Stratum_use_simple_t),
                                           "Ao_generalism" = Feeding_generalism,
                                           "Ao_chewers" = Feeding_mode_c,
                                           #"Ah_abundance" = arthro_herb_abundance,
                                           #"Ah_total_mass" = arthro_herb_abundance * Mean_BodySize,
                                           "Plot" = Plot,
                                           "Year" = Year )], 
          
          "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_above_omni.csv")


write.csv(CWM_Arthropods_above_carni[,list( "Ac_Dispersal" = Dispersal_ability, 
                                          "Ac_BodySize" = logBody_Size,
                                          "Ac_extraint" = Feeding_mode_e,
                                          "Ac_chewers" = Feeding_mode_c,
                                        # "Ac_abundance" = arthro_carni_abundance,
                                        # "Ac_total_mass" = arthro_carni_abundance * Mean_BodySize,
                                          "Plot" = Plot,
                                          "Year" = Year )], 
          "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_above_carni.csv")
write.csv(CWM_Arthropods_below_carni[,list( "Ac_b_Dispersal" = Dispersal_ability, 
                                            "Ac_b_BodySize" = logBody_Size,
                                            "Ac_b_extraint" = Feeding_mode_e,
                                            "Ac_b_chewers" = Feeding_mode_c,
                                            # "Ac_abundance" = arthro_carni_abundance,
                                            # "Ac_total_mass" = arthro_carni_abundance * Mean_BodySize,
                                            "Plot" = Plot,
                                            "Year" = Year )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Arthropods_below_carni.csv")
#write.csv(CWM_Arthropods_above_carni, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_above_carni.csv")
#write.csv(CWM_Arthropods_below_herb, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Arthropods_below_herb.csv")
#write.csv(CWM_Arthropods_below_carni, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data//CWM_Arthropods_below_carni.csv")

#write.table(FD_Arthropods, file = "Data/CWM_data/FD_Arthropods.csv")
#write.csv(FD_Arthropods_ground, file = "Data/CWM_data/FD_Arthropods_ground.csv")
#write.csv(FD_Arthropods_above, file = "Data/CWM_data/FD_Arthropods_above.csv")


######  !!! MORRRRRE TRRAAAAAITS #####
#pollinators = fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/FlosDocs/data/pollinator_traits/18987.txt")

### Butterflies ####
butterflies_traits = fread("Traits/Arthropods/arthropod_traits_c_westphal/butterfly_traits.csv")
butterflies_traits[Species == 'Leptidea_sinapis/reali', Species := 'Leptidea_sinapis_reali']
butterflies_traits[Species == 'Melitaea_athalia-Komplex', Species := 'Melitaea_athalia_Komplex']
butterflies_traits[Species == 'Colias_crocea_',Species := 'Colias_crocea']
butterflies_traits[Species == 'Aglais_io',Species := 'Inachis_io']
butterflies_traits[Species == 'Aglais_urticae',Species := 'Nymphalis_urticae']
butterflies_traits[Species == 'Aricia_agestis',Species := 'Polyommatus_agestis']
#butterflies_traits[Species == 'Aricia_eumedon',Species := 'Nymphalis urticae']
butterflies_traits[Species == 'Cyanaris_semiargus',Species := 'Polyommatus_semiargus']
butterflies_traits[Species == 'Phengaris_arion',Species := 'Maculinea_arion']

# Check distribution
butterflies_traits[, lapply(.SD, hist), .SDcols = c(c( "Phagie"    ,   "Size"       ,  "Migration",    "Fl.period" ,   "Egg.no"  , "Gen.per.year"))]

butterflies_traits[, c(  "logSize",   "logEgg.no", 'Specialism'):= lapply(.SD, log), .SDcols = c(  "Size",   "Egg.no", 'Phagie')]
butterflies_traits[, log_tot_offspring := log(Egg.no*Gen.per.year)]
# Species-level PCA
pca_butterflies = dudi.mix(butterflies_traits[,c( "Specialism"    , 'log_tot_offspring',  "logSize"       ,  "Migration",    "Fl.period" ,   "logEgg.no"  , "Gen.per.year")],  scannf = FALSE, nf = 3)
draw_dudi_mix(pca_butterflies, c(1, 2))

Abundances_butterflies = Abundances_all[Species %in% butterflies_traits$Species,]
but_abundance = Abundances_butterflies[, list(but_abundance = sum(value, na.rm = T)), by = Plot]

butterfly_cwm = my_cwm(butterflies_traits, Abundances_butterflies, c( "Specialism"    , 'log_tot_offspring',   "logSize"       ,"Size"       ,  "Migration",    "Fl.period" ,   "Egg.no",   "logEgg.no"  , "Gen.per.year"),
       'Species', 'Species')
butterfly_cwm = merge.data.table(but_abundance, butterfly_cwm, by = 'Plot')



pca_butterflies = dudi.pca(butterfly_cwm[Plot != 'AEG35', c( "Specialism" , 'log_tot_offspring'  ,   "logSize"       ,  "Migration",    "Fl.period" #,   "logEgg.no"  , "Gen.per.year"
                                                             )],
                           scannf = FALSE, nf = 3)
fviz_pca_biplot(pca_butterflies)
fviz_pca_biplot(pca_butterflies, axes = c(1,3))



write.csv(butterfly_cwm[ , list("but_Specialism" = Specialism      ,
                         "but_Offspring" = log_tot_offspring,
                         "but_Size" = logSize,
                         "but_Mig" = Migration,
                         "but_FlPeriod" = Fl.period,
                         "but_Eggs"= logEgg.no,
                         "but_GenYear" = Gen.per.year,
                         "but_abundance" = but_ab,
                         "but_total_mass" = but_ab*Size,
                         "Plot" = Plot   ,
                         "Year" = Year), ], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Butterflies.csv")



#### Moths ####
moths_traits = data.table(read_excel("Traits/Arthropods/arthropod_traits_c_westphal/Moths\ Traits\ -\ Sagrario.xlsx"))

# Check distribution
moths_traits[, lapply(.SD, hist), .SDcols = c( "Fl.period", "Specialisation", "Gen.per.year", "Wingspan.med")]
moths_traits[, c(  "logSize"):= lapply(.SD, log), .SDcols = c(  "Wingspan.med")]
moths_traits[, c("Fl.period", "Fl.period_corr"):= list(Fl.period*Gen.per.year, Fl.period)]

Abundances_moths = Abundances_all[Species %in% moths_traits$Species,]
moths_cwm = my_cwm(moths_traits, Abundances_moths, c( 'Fl.period','Fl.period_corr', 'logSize', 'Specialisation' ,'Gen.per.year', 'Wingspan.med'),
                       'Species', 'Species')
pca_moths= dudi.pca(moths_cwm[complete.cases(moths_cwm), c('Fl.period', 'Specialisation', 'Gen.per.year', 'logSize')], scannf = FALSE, nf = 2)
fviz_pca_biplot(pca_moths)
write.csv(moths_cwm[, list("mot_FlPeriod" = Fl.period,   
                           "mot_FlPeriod_corr" = Fl.period_corr, 
                           "mot_Size" = logSize, 
                           "mot_Spec" = Specialisation,
                           "mot_GenYear" = Gen.per.year,
                           'Year' = 2008,
                           'Plot' = Plot)], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Moths.csv")

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
coll_traits = fread('Traits/Collembola_Mites/Coll_Traits_2019.csv')
coll_traits[, c( 'Size', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales') := lapply(.SD, as.numeric),
            .SDcols = c( 'Size', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales')]

# Check distribution
coll_traits[, lapply(.SD, function(x){hist(log(x))}), .SDcols = c( 'Size', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales')]
coll_traits[, c( "logSize"):= lapply(.SD, log), .SDcols = c( "Size")]



pca_species = dudi.pca(coll_traits[,c( 'logSize', 'Size', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales')])
pca_coll_species= fviz_pca_biplot(pca_species, habillage = coll_traits$Vertical_preference)
                                                                   

Abundances_coll = Abundances_all[Group_broad =='Collembola',]
Coll_cwm_epi = my_cwm(coll_traits[Vertical_preference != 'Euedaphic',], Abundances_coll, c( 'logSize', 'Size', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'),
                                    'Species', 'Species')
Coll_cwm_eue = my_cwm(coll_traits[Vertical_preference == 'Euedaphic',], Abundances_coll, c( 'logSize','Size', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'),
                      'Species', 'Species')
Coll_cwm_all = my_cwm(coll_traits, Abundances_coll, c( 'logSize','Size', 'Ocelli', 'Pigment', 'Furca', 'PAO', 'PSO', 'Asp', 'Scales'#, 'Vertical_preference'
                                                       ),
                      'Species', 'Species')


pca_coll_all= dudi.pca(Coll_cwm_all[complete.cases(Coll_cwm_all),c('Size',    'Ocelli', 'PAO_1',  'Asp_1',    'Scales_1',
                                                                                               'Pigment_1', 'Furca_1')], scannf = FALSE, nf = 3)
fviz_pca_biplot(pca_coll_all, habillage = coll_traits$Vertical_preference)
pca_coll_eue = dudi.pca(Coll_cwm_eue[complete.cases(Coll_cwm_eue),c('Size',    'Ocelli', 'PAO_1',  'Asp_1',    'Scales_1',
                                                                    'Pigment_1', 'Furca_1')], scannf = FALSE, nf = 3)
fviz_pca_biplot(pca_coll_eue)

write.csv(Coll_cwm_all[, list( "col_Size" = logSize,
                               "col_Ocelli" = Ocelli ,
                               "col_Pigments" = Pigment_1,
                               "col_Furca" = Furca_1,
                               "col_PAO" = PAO_1   ,
                               "col_PSO" = PSO_1 ,
                               "col_Asp" = Asp_1 ,
                               "col_Scales" = Scales_1,
                               "Plot" = Plot ,
                               "Year" = Year     )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Coll_all.csv")
write.csv(Coll_cwm_epi[, list( "colEpi_Size" = logSize,
                               "colEpi_Ocelli" = Ocelli ,
                               "colEpi_Pigments" = Pigment_1,
                               "colEpi_Furca" = Furca_1,
                               "colEpi_PAO" = PAO_1   ,
                               "colEpi_Asp" = Asp_1 ,
                               "colEpi_Scales" = Scales_1,
                               "Plot" = Plot ,
                               "Year" = Year     )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Coll_epi.csv")
write.csv(Coll_cwm_eue[, list( "colEu_Size" = logSize,
                               "colEu_Furca" = Furca_1,
                               "colEu_PAO" = PAO_1   ,
                               "colEu_PSO" = PSO_1 ,
                               "colEu_Asp" = Asp_1 ,
                               "colEu_Scales" = Scales_1,
                               "Plot" = Plot ,
                               "Year" = Year     )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Coll_eue.csv")




### Mites
mites_traits = fread('Traits/Collembola_Mites/2019_Mite_fauna_Traits.csv')
mites_traits[, Size := weighted.mean(c(300,700,1000), c(size100_500, size500_900, size900)), by = 1:nrow(mites_traits)]
mites_traits[, Species := V1]
# Check distribution
mites_traits[, lapply(.SD, function(x){hist(log(x))}), .SDcols = c( 'Mass', 'DaystoAdult')]


Abundances_mites = Abundances_all[Group_broad =='Acari',]
Mites_cwm_all = my_cwm(mites_traits, Abundances_mites, c( 'Surface', 'Litter', 'Soil', 'Reproduction', 'Mass', 'Size', 'CuticuleCalcium', 'DaystoAdult'),
                      'Species', 'Species')
pca_mites_all= dudi.pca(Mites_cwm_primary[complete.cases(Mites_cwm_primary),c(# 'Surface_1', 'Litter_1', 'Soil_1', 
  'Reproduction_sex', 'Mass', 'Size', 'DaystoAdult')], scannf = FALSE, nf = 3)
fviz_pca_biplot(pca_mites_all,
                c(1,2))#, habillage = substr(Mites_cwm_all$Plot,1,1)‹)


Mites_cwm_primary = my_cwm(mites_traits[Feeding != 'omnivorous',], Abundances_mites, c( 'Surface', 'Litter', 'Soil', 'Reproduction', 'Mass', 'Size', 'CuticuleCalcium', 'DaystoAdult'),
                       'Species', 'Species')
pca_mites_all= dudi.pca(Mites_cwm_all[complete.cases(Mites_cwm_all),c(# 'Surface_1', 'Litter_1', 'Soil_1', 
  'Reproduction_parthenogenetic', 'Mass', 'Size',  'CuticuleCalcium_1', 'DaystoAdult')], scannf = FALSE, nf = 3)



Mites_cwm_soil = my_cwm(mites_traits[Vertical_distribution.Morphotype == 'soil',], Abundances_mites, c( 'Reproduction', 'Mass', 'Size', 'CuticuleCalcium', 'DaystoAdult'),
                       'Species', 'Species')
Mites_cwm_litter_surface = my_cwm(mites_traits[Vertical_distribution.Morphotype != 'soil',], Abundances_mites, c( 'Reproduction', 'Mass', 'Size', 'CuticuleCalcium', 'DaystoAdult'),
                        'Species', 'Species')
Mites_cwm_surface= my_cwm(mites_traits[Surface == 1,], Abundances_mites, c( 'Reproduction', 'Mass', 'Size', 'CuticuleCalcium', 'DaystoAdult'),
                         'Species', 'Species')

write.csv(Mites_cwm_primary[, list(             
                                "mites_Surface" = Surface_1,
                                "mites_Soil" = Soil_1,
                                "mites_Sex" = Reproduction_sex,
                                "mites_Mass" = Mass,
                                "mites_Size" = Size            ,
                               # "mites_Cuticule" = CuticuleCalcium_1,
                                "mites_DaysAdult" = DaystoAdult,
                                "Plot" = Plot,
                                "Year" = Year )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Mites_primary.csv")
write.csv(Mites_cwm_soil[, list(
  "mitesSo_Sex" = Reproduction_sex,
  "mitesSo_Mass" = Mass,
  "mitesSo_Size" = Size            ,
  "mitesSo_DaysAdult" = DaystoAdult,
  "Plot" = Plot,
  "Year" = Year )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Mites_soil.csv")
write.csv(Mites_cwm_litter_surface[, list(             
  "mitesL_Sex" = Reproduction_sex,
  "mitesL_Mass" = Mass,
  "mitesL_Size" = Size            ,
  "mitesL_Cuticule" = CuticuleCalcium_1,
  "mitesL_DaysAdult" = DaystoAdult,
  "Plot" = Plot,
  "Year" = Year )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Mites_litter_surface.csv")
write.csv(Mites_cwm_surface[, list(             
  "mitesSu_Sex" = Reproduction_sex,
  "mitesSu_Mass" = Mass,
  "mitesSu_Size" = Size            ,
  "mitesSu_Cuticule" = CuticuleCalcium_1,
  "mitesSu_DaysAdult" = DaystoAdult,
  "Plot" = Plot,
  "Year" = Year )], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Mites_surface.csv")


#### Mollusca ####

#### Acari ####

#### Myriapoda ####
# Looked in BETSI, missing the most abundance species :(
Abundances_all[Group_broad %in% c('Myriapoda'), sum(value), by = Species][order(V1, decreasing = T),]

"Myriapoda"      "Neuroptera"     "Dictyoptera"    "Dermaptera"    
"Protists"       "Opiliones"      "soilfungi"   

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
Heuss2019 <- data.table(read_excel("Traits/Ants/Heuss2019.xlsx", 
                                   col_types = c("text", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric")))
Heuss2019[, Species := gsub(' ', '_', Species)]

Heuss2019[Species == 'Campanotus_ligniperda', Species := 'Campanotus_ligniperdus']
trait_ants = c('Nectar','Plant', 'Zoopha','Tropho','WL', 'Dom',    'CS',  'nQ',  'nN', 'CFT', 'Strata_forage')

ant_abundance = fread('/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Ants/23986_2_Dataset/23986_2_data.csv')
ant_abundance_melt = melt.data.table(ant_abundance, id.vars = c('Species', 'Plot'), measure.vars = 'Presence_absence')
ant_abundance_melt[, Plot := ifelse(nchar(Plot) == 5, Plot, paste0(substr(Plot, 1, 3), '0', substr(Plot, 4, 4)))]
ant_abundance_melt[, Year := '2014_15']
ant_abundance_melt[, Species := gsub(' ', '_', Species)]


CWM_ants = my_cwm(Traits0 = Heuss2019, ant_abundance_melt, trait_names = trait_ants,
                  trait_taxo = 'Species', abundance_taxo = 'Species')

CWM_ants_above = my_cwm(Traits0 = Heuss2019[Strata_forage>0,], ant_abundance_melt, trait_names = trait_ants,
                  trait_taxo = 'Species', abundance_taxo = 'Species')

CWM_ants_below = my_cwm(Traits0 = Heuss2019[Strata_forage<0,], ant_abundance_melt, trait_names = trait_ants,
                        trait_taxo = 'Species', abundance_taxo = 'Species')

fwrite(CWM_ants[, list(
  'ant_zoopha' = Zoopha,
  "ant_nectar" = Nectar      ,
  "ant_plant" = Plant    ,
  "ant_tropho" = Tropho       ,
  "ant_length" = WL     ,
  "ant_dom" = Dom_1       ,
  "ant_colsize" = CS        ,
  "ant_polygyny" = nQ         ,
  "ant_polydomy" = nN   ,
  "ant_colfound" = CFT ,
  "ant_forageabove" = Strata_forage  ,
   Plot,
  Year = 2014
)], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_ants.csv")

fwrite(CWM_ants_above[, list(
  'ant_a_zoopha' = Zoopha,
  "ant_a_nectar" = Nectar      ,
  "ant_a_plant" = Plant    ,
  "ant_a_tropho" = Tropho       ,
  "ant_a_length" = WL     ,
  "ant_a_dom" = Dom_1       ,
  "ant_a_colsize" = CS        ,
  "ant_a_polygyny" = nQ         ,
  "ant_a_polydomy" = nN   ,
  "ant_a_colfound" = CFT ,
  "ant_a_forageabove" = Strata_forage  ,
  Plot,
  Year = 2014
)], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_ants_above.csv")

fwrite(CWM_ants_below[, list(
  'ant_b_zoopha' = Zoopha,
  "ant_b_nectar" = Nectar      ,
  "ant_b_plant" = Plant    ,
  "ant_b_tropho" = Tropho       ,
  "ant_b_length" = WL     ,
  "ant_b_dom" = Dom_1       ,
  "ant_b_colsize" = CS        ,
  "ant_b_polygyny" = nQ         ,
 # "ant_polydomy" = nN   ,
  "ant_b_colfound" = CFT ,
  "ant_b_forageabove" = Strata_forage  ,
  Plot,
  Year = 2014
)], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_ants_below.csv")



