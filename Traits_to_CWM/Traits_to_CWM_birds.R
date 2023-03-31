# This script takes as input the abundances and species-level traits of bird species found 
# in the Exploratories grasslands and outputs a matched trait dataset, a CWM matrix for all considered years and a species-level PCA.



Abundance_birds = Abundance_all[Group_broad == "Birds",]

# Check weird species?
summed = dcast(Abundance_birds[, list(value =sum(value)), by = c('Species', 'Plot')], Plot~Species)
fviz_pca(PCA(summed[-1]))

# ********************************** #
#### 1. Load and merge trait data ####
# ********************************** #

# Birds trait compiled within the Exploratories
bird_traits = fread('Data/Trait_data/210819_explo_birdtraits_fromCat.csv') # Exploratories: 31368
bird_traits[grepl('Delichon_urbica', species_latin ),  species_latin:= 'Delichon_urbicum']

# Pigot 2020 for strategies and trophic levels
Traits_Pigot2020 <- data.table(read_excel("Data/Trait_data/Traits_Pigot2020.xlsx")) # https://www.nature.com/articles/s41559-019-1070-4#Sec24 supplementary dataset 1

# Avonet 
Avonet <- fread("Data/Trait_data/AVONET1_BirdLife.csv") # https://figshare.com/articles/dataset/AVONET_morphological_ecological_and_geographical_data_for_all_birds_Tobias_et_al_2021_Ecology_Letters_/16586228

# Bird et al 2020 for generation time
birds_gen_time = data.table(read_excel('Data/Trait_data/Bird_et_al_2020.xlsx')) #https://conbio.onlinelibrary.wiley.com/doi/full/10.1111/cobi.13486 Table S1
birds_gen_time[, scientific_name := `Scientific name`]

# Species matching based on gbif
bird_traits$species_latin[!(bird_traits$species_latin %in% gsub(' ', '_', birds_gen_time$scientific_name))]

birds_gen_time = rbind(birds_gen_time,
                       rbind(birds_gen_time[grepl('Corvus corone', scientific_name),][, scientific_name := 'Corvus_corone_cornix'],
                             birds_gen_time[grepl('Corvus corone', scientific_name),][, scientific_name := 'Corvus_corone_corone']))

birds_gen_time[grepl('Parus', scientific_name), unique(scientific_name)]

birds_gen_time[grepl('Linaria cannabina',     birds_gen_time$scientific_name),  scientific_name:= 'Carduelis_cannabina']
birds_gen_time[grepl('Spinus spinus',         birds_gen_time$scientific_name),  scientific_name:= 'Carduelis_spinus']
birds_gen_time[grepl('Chloris chloris',       birds_gen_time$scientific_name),  scientific_name:= 'Carduelis_chloris']
birds_gen_time[grepl('Lophophanes cristatus', birds_gen_time$scientific_name),  scientific_name:= 'Parus cristatus']
birds_gen_time[grepl('Poecile palustris',     birds_gen_time$scientific_name),  scientific_name:= 'Parus_palustris']
birds_gen_time[grepl('Regulus ignicapilla',   birds_gen_time$scientific_name),  scientific_name:= 'Regulus ignicapillus']
birds_gen_time[grepl('Sylvia curruca',        birds_gen_time$scientific_name),  scientific_name:= 'Sylvia_curucca']
birds_gen_time[grepl('Dryobates minor',       birds_gen_time$scientific_name),  scientific_name:= 'Dendrocopos minor']
birds_gen_time[grepl('Saxicola rubicola',     birds_gen_time$scientific_name),  scientific_name:= 'Motacilla rubicola']
birds_gen_time[grepl('Leiopicus_medius',     birds_gen_time$scientific_name),  scientific_name := 'Dendrocopos_medius']
birds_gen_time[grepl('Poecile_montanus',      birds_gen_time$scientific_name),  scientific_name := 'Parus_montanus']

birds_gen_time[, species_latin := gsub(' ', '_', scientific_name) ]#
birds_gen_time[, Max_longevity := `Maximum longevity` ]#

# Species matching based on gbif for Pigot et al. 2020
Traits_Pigot2020 = rbind(Traits_Pigot2020,
                         rbind(Traits_Pigot2020[grepl('Corvus_corone', Traits_Pigot2020$Binomial),][, Binomial := 'Corvus_corone_cornix'],
                               Traits_Pigot2020[grepl('Corvus_corone', Traits_Pigot2020$Binomial),][, Binomial := 'Corvus_corone_corone']))

Traits_Pigot2020[grepl('Parus_caeruleus', Traits_Pigot2020$Binomial),     Binomial := 'Cyanistes_caeruleus']
Traits_Pigot2020[grepl('Miliaria_calandra', Traits_Pigot2020$Binomial),   Binomial := 'Emberiza_calandra']
Traits_Pigot2020[grepl('Parus_ater', Traits_Pigot2020$Binomial),          Binomial := 'Periparus_ater']
Traits_Pigot2020[grepl('Regulus_ignicapilla', Traits_Pigot2020$Binomial), Binomial := "Regulus_ignicapillus"]
Traits_Pigot2020[grepl('Saxicola_torquatus', Traits_Pigot2020$Binomial),  Binomial := "Saxicola_rubicola" ]
Traits_Pigot2020[grepl('Sylvia_curruca', Traits_Pigot2020$Binomial),      Binomial := "Sylvia_curucca"  ]
Traits_Pigot2020[grepl('Leiopicus_medius',      Traits_Pigot2020$Binomial),  Binomial := 'Dendrocopos_medius']
Traits_Pigot2020[grepl('Poecile_montanus',      Traits_Pigot2020$Binomial),  Binomial := 'Parus_montanus']
Traits_Pigot2020[, species_latin := Binomial ]

# Species matching based on gbif for Avonet
Avonet = rbind(Avonet, rbind(Avonet[grepl('Corvus corone', Avonet$Species1),][, Species1 := 'Corvus_corone_cornix'],
                             Avonet[grepl('Corvus corone', Avonet$Species1),][, Species1 := 'Corvus_corone_corone']))

Avonet[grepl('Parus caeruleus', Avonet$Species1),     Species1 := 'Cyanistes_caeruleus']
Avonet[grepl('Miliaria calandra', Avonet$Species1),   Species1 := 'Emberiza_calandra']
Avonet[grepl('Parus ater', Avonet$Species1),          Species1 := 'Periparus_ater']
Avonet[grepl('Regulus ignicapilla', Avonet$Species1), Species1 := "Regulus_ignicapillus"]
Avonet[grepl('Saxicola torquatus', Avonet$Species1),  Species1 := "Saxicola_rubicola" ]
Avonet[grepl('Sylvia curruca', Avonet$Species1),      Species1 := "Sylvia_curucca"  ]
Avonet[grepl('Leiopicus medius',      Avonet$Species1),  Species1:= 'Dendrocopos_medius']
Avonet[grepl('Poecile montanus',      Avonet$Species1),  Species1:= 'Parus_montanus']
Avonet[grepl('Linaria cannabina',      Avonet$Species1),  Species1:= 'Carduelis_cannabina']
Avonet[grepl('Chloris chloris',        Avonet$Species1),  Species1:= 'Carduelis_chloris']
Avonet[grepl('Spinus spinus',         Avonet$Species1),  Species1:= 'Carduelis_spinus']
Avonet[grepl('Dryobates minor',       Avonet$Species1),  Species1:= 'Dendrocopos minor']
Avonet[grepl('Lophophanes cristatus',  Avonet$Species1),  Species1:= 'Parus cristatus']
Avonet[grepl('Poecile palustris',      Avonet$Species1),  Species1:= 'Parus_palustris']
Avonet[, species_latin := gsub(' ', '_', Species1) ]

#### Check unmatched species ####
bird_traits[!(species_latin %in% Traits_Pigot2020$species_latin), unique(species_latin)]
bird_traits[!(species_latin %in% Avonet$species_latin), unique(species_latin)]


### Merge
bird_traits = merge.data.table(bird_traits, Traits_Pigot2020[, .SD, .SDcols = c('TrophicLevel', 'TrophicNiche','ForagingNiche','species_latin')],
                               by = 'species_latin')
bird_traits = merge.data.table(bird_traits, birds_gen_time[, .SD, .SDcols = c('Max_longevity', 'GenLength','species_latin')],
                               by = 'species_latin', all.x = T)
bird_traits = merge.data.table(bird_traits, Avonet[, .SD, .SDcols = c('Primary.Lifestyle','species_latin')],
                               by = 'species_latin', all.x = T)


## Check trophic levels, especially for omnivores. Omnivorous birds that also eat insects or or animals are classified as tertiary consumers
bird_traits[, table(TrophicLevel, TrophicNiche)]

bird_traits[Functional_Group == 'carnivore' , trophic_level  := 'carnivore']
bird_traits[Functional_Group == 'insectivore',  trophic_level  := 'insectivore']
bird_traits[Functional_Group %in% c('granivore',   'herbivore'), trophic_level  := 'herbivore']
bird_traits[Functional_Group == 'omnivore' & TrophicLevel == 'Carnivore', trophic_level  := 'carnivore']            
bird_traits[Functional_Group == 'omnivore' & TrophicLevel == 'Herbivore', trophic_level  := 'herbivore']            
bird_traits[Functional_Group == 'omnivore' & TrophicLevel == 'Omnivore', trophic_level  := 'insectivore']            


# Transform traits
bird_traits[, c('logWing_len',
                'logBody_len',
                'logTail_len',
                'logBill_len',
                'logTarsus_len',
                'logBody_Mass',
                'logNestling_st',
                'logWing_s',
                'logWing_a',
                'logClutch',
                'logIncub_time',
                'logAge',
                "logLongevity",
                'logBrood'#,
               # 'log_GenLength'
                )        :=   lapply(.SD, log),
            .SDcols = c('wing_length_max',
                        'body_length_max',
                        'tail_length_max',
                        'bill_length_max',
                        'tarsus_length_max',
                        'body_mass_max',
                        'Nestling_stage_max',
                        'wing_span_max',
                        'wing_area',
                        'clutch_size_max',
                        'incubation_time_max',
                        'maximum_age_EURING',
                        "Max_longevity",
                        'maximum_broods_per_year')]

bird_traits[, Number_offspring := as.numeric(clutch_size_max)*as.numeric(maximum_broods_per_year)]
bird_traits[, logNumber_offspring := log(Number_offspring)]


# All traits
all_traits = c('logWing_len',
               'logBody_len',
               'logTail_len',
               'logBill_len',
               'logTarsus_len',
               'logBody_Mass',
               'logNestling_st',
               'logWing_s',
               'logWing_a',
               'logClutch',
               'logIncub_time',
               'logAge',
               "logLongevity",
               'logBrood', 
               'logNumber_offspring',
               'GenLength')


# We will actually use only a subset of traits
trait_selection =  c( 'GenLength',
                      "logIncub_time",
                      'logNumber_offspring',
                      "logLongevity",
                      'logBody_Mass'
)


# Add info
traitUnits = c("mm (log-transformed)","mm (log-transformed)","mm (log-transformed)","mm (log-transformed)","mm (log-transformed)","g (log-transformed)",'day (log-transformed)',"mm (log-transformed)","mm^2 (log-transformed)",'(log-transformed)','day (log-transformed)','year (log-transformed)','year (log-transformed)','(log-transformed)','(log-transformed)', 'year', 'Hz', 'unitless', 'unitless')
traitDescription = c("Maximum wing length, measured from bow to tip (typically flattened)",
                     "Maximum body length of live bird from bill tip to longest tail feather",
                     "Maximum length of tail, measured after Svensson",
                     "Maximum length of bill from tip to front of cranium",
                     "Maximum tarsus length; von Blotzheim, U.N.G.; Bauer, K.M. Handbuch Der V?gel Mitteleuropas; Aula: Wiesbaden, Germany, 1998.",
                     "Maximum body mass from the literature",
                     'Maximum number of days a nestling is on the nest, typically from hatching to flying',
                     "Maximum wingspan of flattened and stretched wings from one wing tip to the other",
                     "Maximum area of the wing covered if stretched fully",
                     'Maximum number of eggs in clutch',
                     'Maximum number of days from egg laying to hatching',
                     'Maximum age observed of banded birds (in years, rounded to the next year, if any month was specified)',
                     'Maximum longevity observed from both wild and captive animals',
                     'Maximum number of broods observed per species so far',
                     'Maximum number of offspring per year (calculated manually as number of broods * clutch size)',
                     'Generation length',
                     "Wingbeat frequency",
                     'Main foraging niche',
                     'Lifestyle')
traitDataID = c("Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","Bexis ID 20067, 31368","NA","Bexis ID 20067, 31368",'NA', 'NA', "Bexis ID 20067, 31368", 'NA', 'NA')

traitRef = c("https://www.mdpi.com/2306-5729/2/2/12/htm; Svensson, L. Identification Guide to European Passerines; British Trust for Ornithology: Stockholm, Sweden, 1992.;
             von Blotzheim, U.N.G.; Bauer, K.M. Handbuch Der V?gel Mitteleuropas; Aula: Wiesbaden, Germany, 1998.","https://www.mdpi.com/2306-5729/2/2/12/htm","https://www.mdpi.com/2306-5729/2/2/12/htm",
             "https://www.mdpi.com/2306-5729/2/2/12/htm; Svensson, L. Identification Guide to European Passerines; British Trust for Ornithology: Stockholm, Sweden, 1992.;
             von Blotzheim, U.N.G.; Bauer, K.M. Handbuch Der V?gel Mitteleuropas; Aula: Wiesbaden, Germany, 1998.",
             "https://www.mdpi.com/2306-5729/2/2/12/htm","https://www.mdpi.com/2306-5729/2/2/12/htm",'https://www.mdpi.com/2306-5729/2/2/12/htm',"https://www.mdpi.com/2306-5729/2/2/12/htm","https://www.mdpi.com/2306-5729/2/2/12/htm", 'https://www.mdpi.com/2306-5729/2/2/12/htm','https://www.mdpi.com/2306-5729/2/2/12/htm','https://www.mdpi.com/2306-5729/2/2/12/htm','https://doi.org/10.1111/cobi.13486',
             'https://www.mdpi.com/2306-5729/2/2/12/htm', 'Calculated manually', ' https://doi.org/10.1111/cobi.13486', 'https://www.mdpi.com/2306-5729/2/2/12/htm',  'https://doi.org/10.1038/s41559-019-1070-4', 'https://doi.org/10.1038/s41559-019-1070-4')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c('wing_length_max',
                                                                                       'body_length_max',
                                                                                       'tail_length_max',
                                                                                       'bill_length_max',
                                                                                       'tarsus_length_max',
                                                                                       'body_mass_max',
                                                                                       'Nestling_stage_max',
                                                                                       'wing_span_max',
                                                                                       'wing_area',
                                                                                       'clutch_size_max',
                                                                                       'incubation_time_max',
                                                                                       'maximum_age_EURING',
                                                                                       "Max_longevity",
                                                                                       'maximum_broods_per_year',
                                                                                       'Number_offspring',
                                                                                       'GenLength',
                                                                                       'wingbeat_frequency',
                                                                                       'ForagingNiche',
                                                                                       'Primary.Lifestyle')

bird_traits_melt = melt.data.table(unique(bird_traits[, .SD, .SDcols = c('wing_length_max',
                                                                         'body_length_max',
                                                                         'tail_length_max',
                                                                         'bill_length_max',
                                                                         'tarsus_length_max',
                                                                         'body_mass_max',
                                                                         'Nestling_stage_max',
                                                                         'wing_span_max',
                                                                         'wing_area',
                                                                         'clutch_size_max',
                                                                         'incubation_time_max',
                                                                         'maximum_age_EURING',
                                                                         "Max_longevity",
                                                                         'maximum_broods_per_year',
                                                                         'Number_offspring',
                                                                         'GenLength',
                                                                         'wingbeat_frequency',
                                                                         'ForagingNiche',
                                                                         'Primary.Lifestyle','species_latin')]), id.vars = c( 'species_latin'), , variable.name = 'traitName', value.name = 'traitValue')
bird_traits_info = add_info(bird_traits_melt, traitRefs = traitRef, traitDataIDs = traitDataID, traitDescriptions = traitDescription, traitUnits = traitUnits, traitsOnly = TRUE)
fwrite(bird_traits_info, "Data/Temporary_data/Bird_traits.csv")


# ************************** #
#### 2. Species-level PCA ####
# ************************** #

# Insectivores and Carnivore only
# Without imputation
pca_birds = dudi.pca(bird_traits[trophic_level %in% c('insectivore','carnivore'),][complete.cases(bird_traits[trophic_level %in% c('insectivore','carnivore'), ..trait_selection]), ..trait_selection], scannf = FALSE, nf = 3)
pca_birds= fviz_pca_biplot(pca_birds, geom = c("point"), repel = T, axes = c(1,2), title = 'Insectivore birds')
# With imputation
pca_birds_sp = dudi.pca(mice::complete(mice(bird_traits[trophic_level %in% c('insectivore','carnivore'),list(GenLength, Incub.log =  logIncub_time,    NOffspring.log = logNumber_offspring, Longevity.log = logLongevity, Body_len.log =logBody_len)])), scannf = FALSE, nf = 2)
gg_birds_sp = fviz_pca(pca_birds_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                       col.ind = "steelblue",
                       fill.ind = "white",
                       col.var = "black")
ggsave(gg_birds_sp, file =  'Results/species_pca_birds.pdf', width = 6, height = 5)




# ********************************** #
#### 3. Community-weighted traits ####
# ********************************** #
# Coverage
CC_birds_insect = check_coverage(bird_traits[trophic_level %in% c('insectivore','carnivore'),], 
                                 Abundance_birds[Species %in% bird_traits[trophic_level %in% c('insectivore','carnivore'),species_latin],], 
                                 all_traits,'species_latin', 'Species')

#### Community weighted mean ####
CWM_birds_insect = my_cwm(bird_traits[trophic_level %in% c('insectivore','carnivore'),], 
                                      Abundance_birds[Species %in% bird_traits[trophic_level %in% c('insectivore','carnivore'),species_latin],],
                                      all_traits,'species_latin', 'Species')


# Melt and merge
CWM_CC_birds_insect = merge.data.table(melt.data.table(CWM_birds_insect[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                             melt.data.table(CC_birds_insect[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


# Add info
traitUnits = c("mm (log-transformed)","mm (log-transformed)","mm (log-transformed)","mm (log-transformed)","mm (log-transformed)","g (log-transformed)",'day (log-transformed)',"mm (log-transformed)","mm^2 (log-transformed)",'(log-transformed)','day (log-transformed)','year (log-transformed)','year (log-transformed)','(log-transformed)','(log-transformed)', 'year')
traitDescription = c("Maximum wing length, measured from bow to tip (typically flattened)",
                     "Maximum body length of live bird from bill tip to longest tail feather",
                     "Maximum length of tail, measured after Svensson",
                     "Maximum length of bill from tip to front of cranium",
                     "Maximum tarsus length; von Blotzheim, U.N.G.; Bauer, K.M. Handbuch Der V?gel Mitteleuropas; Aula: Wiesbaden, Germany, 1998.",
                     "Maximum body mass from the literature",
                     'Maximum number of days a nestling is on the nest, typically from hatching to flying',
                     "Maximum wingspan of flattened and stretched wings from one wing tip to the other",
                     "Maximum area of the wing covered if stretched fully",
                     'Maximum number of eggs in clutch',
                     'Maximum number of days from egg laying to hatching',
                     'Maximum age observed of banded birds (in years, rounded to the next year, if any month was specified)',
                     'Maximum longevity observed from both wild and captive animals',
                     'Maximum number of broods observed per species so far',
                     'Maximum number of offspring per year (calculated manually as number of broods * clutch size)',
                     'Generation length')
traitDataID = c("Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","Bexis ID 20067, x","","Bexis ID 20067, x",'', '')

traitRef = c("https://www.mdpi.com/2306-5729/2/2/12/htm; Svensson, L. Identification Guide to European Passerines; British Trust for Ornithology: Stockholm, Sweden, 1992.;
             von Blotzheim, U.N.G.; Bauer, K.M. Handbuch Der V?gel Mitteleuropas; Aula: Wiesbaden, Germany, 1998.","https://www.mdpi.com/2306-5729/2/2/12/htm","https://www.mdpi.com/2306-5729/2/2/12/htm",
             "https://www.mdpi.com/2306-5729/2/2/12/htm; Svensson, L. Identification Guide to European Passerines; British Trust for Ornithology: Stockholm, Sweden, 1992.;
             von Blotzheim, U.N.G.; Bauer, K.M. Handbuch Der V?gel Mitteleuropas; Aula: Wiesbaden, Germany, 1998.",
             "https://www.mdpi.com/2306-5729/2/2/12/htm","https://www.mdpi.com/2306-5729/2/2/12/htm",'https://www.mdpi.com/2306-5729/2/2/12/htm',"https://www.mdpi.com/2306-5729/2/2/12/htm","https://www.mdpi.com/2306-5729/2/2/12/htm", 'https://www.mdpi.com/2306-5729/2/2/12/htm','https://www.mdpi.com/2306-5729/2/2/12/htm','https://www.mdpi.com/2306-5729/2/2/12/htm','https://doi.org/10.1111/cobi.13486',
             'https://www.mdpi.com/2306-5729/2/2/12/htm', 'Calculated manually', '')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = all_traits


CWM_CC_birds_insect = add_info(CWM_CC_birds_insect, traitRef, traitDataID, traitDescription, traitUnits, c('21446, 21447, 21448, 21449, 24690, 25306 synthesised in 27707'))

fwrite(CWM_CC_birds_insect, "Data/CWM_data/CWM_birds.csv")

# fwrite(CWM_birds_insect[, list( "Plot" = Plot               ,
#                                "Year" = Year,
#                                "Bi_Size" = log_body_len     ,
#                                "Bi_Incub" = logIncub_time,
#                                "Bi_TOffsprings" = log_offspring  ,
#                                "Bi_AgeMax" = log_longevity,
#                                "Bi_GenLength" = GenLength)
#
#], paste(cwm_path, "CWM_birds_insect.csv", sep = ''))


# ************************************** #
#### 4. Non-weighted community traits ####
# ************************************** #

# We also calculate unweighted traits (i.e. average trait of all species present in a plot) to check the ffect of 
# turnover v. changes in abundance

Abundances_birds_presence_absence = Abundance_birds[, list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Species)]
Abundances_birds_presence_absence[value>1, value := 1]

CC_birds_insect_noweight   = check_coverage(bird_traits[trophic_level %in% c('insectivore','carnivore'),], 
                                           Abundances_birds_presence_absence[Species %in% bird_traits[trophic_level %in% c('insectivore','carnivore'),species_latin],], 
                                 all_traits,'species_latin', 'Species')
CWM_birds_insect_noweight  = my_cwm(bird_traits[trophic_level %in% c('insectivore','carnivore'),], Abundances_birds_presence_absence, c(all_traits),'species_latin', 'Species')


# Melt and merge
CWM_CC_birds_insect_noweight = merge.data.table(melt.data.table(CWM_birds_insect_noweight[, Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                       melt.data.table(CC_birds_insect_noweight[, Year := NA], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))

CWM_CC_birds_insect_noweight = add_info(CWM_CC_birds_insect_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('21446, 21447, 21448, 21449, 24690, 25306 synthesised in 27707'))

fwrite(CWM_CC_birds_insect_noweight, "Data/CWM_data/CWM_birds_noweight.csv")

#fwrite(CWM_birds_insect_noweight[, list( "Plot" = Plot               ,
#                                "Year" = Year,
#                                "Bi_Size" = log_body_len     ,
#                                "Bi_Incub" = logIncub_time,
#                                "Bi_TOffsprings" = log_offspring  ,
#                                "Bi_AgeMax" = log_longevity,
#                                "Bi_GenLength" = GenLength)
                        
#], paste(cwm_path, "CWM_birds_insect_noweight.csv", sep = ''))


# ************************************ #
#### Turnover accross LUI gradient ####
# ************************************ #

# A more direct way to check for turnover

data_lui <- fread("Data/Environment_function_data/LUI_standardized_global.txt") # from https://www.bexis.uni-jena.de/lui/LUICalculation/index; new components, standardised, global, all regions, all years
data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

comm.test = data.frame(dcast(Abundance_birds[Species %in% bird_traits[trophic_level %in% c('insectivore','carnivore'), species_latin], list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Species)],  Plot~Species, value.var = 'value', fill = 0))
rownames(comm.test)= comm.test$Plot
comm.test = comm.test[,-1]

beta.multi.abund(comm.test)

comm_min_max = matrix(c(colSums(comm.test[min_lui_plots,]),colSums(comm.test[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_min_max)
