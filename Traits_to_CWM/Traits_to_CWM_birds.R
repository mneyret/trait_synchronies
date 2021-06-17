library(data.table)
library(reshape2)
library(vegan)
library(taxize)
library(FD)
library(readr)
library(textclean)
library(readxl)
library(factoextra)
library(dplyr)
library(mice)

setwd("/Users/Margot/Desktop/Research/Senckenberg/Data/")
Abundances_all = fread("Abundances/Dataset_clean.txt")
Abundances = Abundances_all[Group_broad %in% c("Birds"),] 

bird_traits_1 = fread('/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Birds/20067_2_Dataset/20067_2_data.csv')
bird_traits_1[, Species := paste(Genus, gsub(' ', '_',species), sep = '_')]
library(traitdataform)
bird_traits_1[, Species_gbif := lapply(Species, function(x){
             tax = get_gbif_taxonomy(x)
             correct = unlist(gsub(' ', '_', tax$scientificName))
             return(correct)}), by = Species]
bird_traits_1[Species == 'Columba_livia_', Species_gbif := 'Columba_livia']
Abundances[Species == 'Columba_livia_', Species_gbif := 'Columba_livia']

Abundances[, Species_gbif := lapply(Species, function(x){
  tax = get_gbif_taxonomy(x)
  correct = unlist(gsub(' ', '_', tax$scientificName))
  return(correct)}), by = Species]

bird_traits_1[, c('maximum_broods_per_year', 'maximum_age_EURING','wingbeat_frequency','wing_area') :=
              lapply(.SD, as.numeric), .SDcols = c('maximum_broods_per_year', 'maximum_age_EURING','wingbeat_frequency','wing_area')]

bird_traits = data.table(complete(mice(bird_traits_1)))
colnames(bird_traits)
bird_traits[, hist(tail_length_max)]

bird_traits[, c('log_wing_len',
                'log_body_len',
                'log_tail_len',
                'log_bill_len',
                'log_tarsus_len',
                'log_body_m',
                'log_Nestling_st',
                'log_wing_s',
                'log_wing_a') := lapply(.SD, log),
            .SDcols = c('wing_length_max',
                        'body_length_max',
                        'tail_length_max',
                        'bill_length_max',
                        'tarsus_length_max',
                        'body_mass_max',
                        'Nestling_stage_max',
                        'wing_span_max',
                        'wing_area')]

bird_traits[, wing_body_ratio := wing_length_max/body_length_max]
bird_traits[, tot_offspring := as.numeric(clutch_size_max)*as.numeric(maximum_broods_per_year)]
bird_traits[, fun_group := ifelse(Functional_Group %in% c('carnivore', 'omnivore'), 'carnivore',
                                  ifelse(Functional_Group == 'insectivore','insectivore',
                                         'herbivore'))]

trait_selection = c('log_wing_len',
                  'log_body_len',
                  'log_tail_len',
                  'log_bill_len',
                  'log_tarsus_len',
                  'log_body_m',
                  'log_Nestling_st',
                  'log_wing_s',
                  'clutch_size_max',
                  "maximum_broods_per_year" ,
                  "incubation_time_max",
                  'tot_offspring',
                  "maximum_age_EURING",
                  'wing_body_ratio'
                  )


### Species-level PCA
pca_birds = dudi.pca(bird_traits[, .SD, .SDcols = trait_selection], scannf = FALSE, nf = 3)
fviz_pca(pca_birds)


### Community-level

# Herbivores
CWM_birds_herb = my_cwm(bird_traits[fun_group == 'herbivore',], Abundances, trait_selection,'Species_gbif', 'Species_gbif')

# Insectivores
CWM_birds_insect = my_cwm(bird_traits[fun_group == 'insectivore',], Abundances, trait_selection,'Species_gbif', 'Species_gbif')

# Omnivores and predators
CWM_birds_pred = my_cwm(bird_traits[fun_group == 'carnivore',], Abundances, trait_selection,'Species_gbif', 'Species_gbif')

fwrite(CWM_birds_pred[, list( "Plot" = Plot               ,
                              "Year" = Year,
                                 "Bp_WingLen" = log_wing_len       ,
                                 "Bp_BodyLen" = log_body_len      ,
                                 "Bp_TailLen" = log_tail_len      ,
                                 "Bp_BillLen" = log_bill_len      ,
                                 "Bp_TarsusLen" = log_tarsus_len  ,
                                 "Bp_BodyMass" = log_body_m     ,
                                 "Bp_Nestling" = log_Nestling_st,
                                 "Bp_WingS" = log_wing_s,
                                 "Bp_ClutchS" = clutch_size_max,
                                # ppB_BroodMax" = maximum_broods_per_year,
                                 "Bp_Incub" = incubation_time_max,
                                # "Bp_TOffsprings" = tot_offspring  ,
                                 "Bp_AgeMax" = maximum_age_EURING,
                                 "Bp_WingBody" = wing_body_ratio
                                )
                          ], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Birds_pred.csv")
fwrite(CWM_birds_herb[, list( "Plot" = Plot               ,
                               "Year" = Year,
                                                     "Bh_WingLen" = log_wing_len       ,
                                                     "Bh_BodyLen" = log_body_len      ,
                                                     "Bh_TailLen" = log_tail_len      ,
                                                     "Bh_BillLen" = log_bill_len      ,
                                                     "Bh_TarsusLen" = log_tarsus_len  ,
                                                     "Bh_BodyMass" = log_body_m     ,
                                                     "Bh_Nestling" = log_Nestling_st,
                                                     "Bh_WingS" = log_wing_s,
                                                     "Bh_ClutchS" = clutch_size_max,
                                                     "Bh_BroodMax" = maximum_broods_per_year,
                                                     "Bh_Incub" = incubation_time_max,
                                                   #  "Bh_TOffsprings" = tot_offspring  ,
                                                     "Bh_AgeMax" = maximum_age_EURING,
                                                     "Bh_WingBody" = wing_body_ratio
)
], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Birds_herb.csv")
fwrite(CWM_birds_insect[, list( "Plot" = Plot               ,
                         "Year" = Year,
                             "Bi_WingLen" = log_wing_len       ,
                             "Bi_BodyLen" = log_body_len      ,
                             "Bi_TailLen" = log_tail_len      ,
                             "Bi_BillLen" = log_bill_len      ,
                             "Bi_TarsusLen" = log_tarsus_len  ,
                             "Bi_BodyMass" = log_body_m     ,
                             "Bi_Nestling" = log_Nestling_st,
                             "Bi_WingS" = log_wing_s,
                             "Bi_ClutchS" = clutch_size_max,
                             "Bi_BroodMax" = maximum_broods_per_year,
                             "Bi_Incub" = incubation_time_max,
                         #    "Bi_TOffsprings" = tot_offspring  ,
                             "Bi_AgeMax" = maximum_age_EURING,
                             "Bi_WingBody" = wing_body_ratio#,
                             #"B_prop_herb" = fun_group_herbivore
)
], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_birds_insect.csv")

