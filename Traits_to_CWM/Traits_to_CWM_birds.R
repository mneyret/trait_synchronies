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
# Abundances
## Raw diversity
allsp <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/210112_EP_species_diversity_GRL_BEXIS.txt")
allsp$Species = gsub('_$', '', allsp$Species ) # Remove _ if last character
## Species information
fgs <- fread("//Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/210112_EP_species_info_GRL_BEXIS.txt")
fgs$Species = gsub(' $', '', fgs$Species )
fgs$Species = gsub(' ', '_', fgs$Species )
Abundance_all <- merge.data.table(allsp, fgs, by ="Species", all.x=TRUE)
Abundance_all[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1, 3), '0', substr(Plot, 4, 4), sep = ''))]

# Explo birds traits
bird_traits = fread('/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Birds/210819_explo_birdtraits_fromCat.csv')
bird_traits[grepl('Delichon_urbica', species_latin ),  species_latin:= 'Delichon_urbicum']

# Pigot 2020 for strategies and trophic levels
Traits_Pigot2020 <- data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Birds/Traits_Pigot2020.xlsx"))

# Avonet for strategies and trophic levels
Avonet <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Birds/AVONET/ELEData/TraitData/AVONET1_BirdLife.csv")

# Bird et al 2020 for generation time
birds_gen_time = data.table(read_excel('/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Birds/Bird_et_al_2020.xlsx'))
birds_gen_time[, scientific_name := `Scientific name`]

# Species matching based on gbif for McMahon et al. 2020
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
Avonet = rbind(Avonet,
                         rbind(Avonet[grepl('Corvus corone', Avonet$Species1),][, Species1 := 'Corvus_corone_cornix'],
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


#### Merge
bird_traits = merge.data.table(bird_traits, Traits_Pigot2020[, .SD, .SDcols = c('TrophicLevel', 'TrophicNiche','ForagingNiche','species_latin')],
by = 'species_latin')
bird_traits = merge.data.table(bird_traits, birds_gen_time[, .SD, .SDcols = c('Max_longevity', 'GenLength','species_latin')],
                               by = 'species_latin', all.x = T)
bird_traits = merge.data.table(bird_traits, Avonet[, .SD, .SDcols = c('Primary.Lifestyle','species_latin')],
                               by = 'species_latin', all.x = T)


## Check trophic levels, especially for omnivores
bird_traits[, table(TrophicLevel, TrophicNiche)]

bird_traits[Functional_Group == 'carnivore' , trophic_level  := 'carnivore']
bird_traits[Functional_Group == 'insectivore',  trophic_level  := 'insectivore']
bird_traits[Functional_Group %in% c('granivore',   'herbivore'), trophic_level  := 'herbivore']
bird_traits[Functional_Group == 'omnivore' & TrophicLevel == 'Carnivore', trophic_level  := 'carnivore']            
bird_traits[Functional_Group == 'omnivore' & TrophicLevel == 'Herbivore', trophic_level  := 'herbivore']            
bird_traits[Functional_Group == 'omnivore' & TrophicLevel == 'Omnivore', trophic_level  := 'insectivore']            


bird_traits[, c('log_wing_len',
                'log_body_len',
                'log_tail_len',
                'log_bill_len',
                'log_tarsus_len',
                'log_body_m',
                'log_Nestling_st',
                'log_wing_s',
                'log_wing_a',
                'log_clutch',
                'log_incub',
                'logAge',
                "log_longevity",
                'log_brood') := lapply(.SD, log),
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

bird_traits[, wing_body_ratio := wing_length_max/body_length_max]
bird_traits[, log_offspring := log(as.numeric(clutch_size_max)*as.numeric(maximum_broods_per_year))]

trait_pols = c( 'GenLength',
                  "log_incub",
                  'log_offspring',
                  "log_longevity"
                  )
trait_morpho =c(
                'log_body_len'
)
trait_selection = c(trait_pols, trait_morpho)
traits_main = c(trait_pols[trait_pols != 'tot_offspring'], 'log_body_m', 'wing_body_ratio')

### Compare species with traits and all species
Abundances[(Species %in% bird_traits[trophic_level %in% c('insectivore','carnivore'),species_latin]), unique(Species)]
Abundances[(Species %in% bird_traits[!(trophic_level %in% c('insectivore','carnivore')),species_latin]) & value >0, unique(Plot)]




### Species-level PCA

# Insectivores and Carnivore
pca_birds_insect = dudi.pca(bird_traits[trophic_level %in% c('insectivore','carnivore'),][complete.cases(bird_traits[trophic_level %in% c('insectivore','carnivore'), ..trait_selection]), ..trait_selection], scannf = FALSE, nf = 3)
pca_birds_insect_species= fviz_pca_biplot(pca_birds_insect, geom = c("point"), repel = T, axes = c(1,2), title = 'Insectivore birds')
ggsave(pca_birds_insect_species,file= '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Species_PCA_Birds_insect.pdf', width = 5, height = 5)    

pca_birds_sp = dudi.pca(mice::complete(mice(bird_traits[trophic_level %in% c('insectivore','carnivore'),list(GenLength, Incub.log =  log_incub,    NOffspring.log = log_offspring, Longevity.log = log_longevity, Body_len.log =log_body_len)])), scannf = FALSE, nf = 2)
gg_birds_sp = fviz_pca(pca_birds_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                       col.ind = "steelblue",
                       fill.ind = "white",
                       col.var = "black")
ggsave(gg_birds_sp, file = '/Users/Margot/Desktop/Research/Senckenberg/Documents/Papers/Traits/Figures/species_pca_birds.pdf', width = 6, height = 5)


# Coverage
CC_birds_insect = check_coverage(bird_traits[trophic_level %in% c('insectivore','carnivore'),], 
                                 Abundances[Species %in% bird_traits[trophic_level %in% c('insectivore','carnivore'),species_latin],], 
                                 c(trait_selection),'species_latin', 'Species')

# Community weighted mean
CWM_birds_insect = my_cwm(bird_traits[trophic_level %in% c('insectivore','carnivore'),], Abundances, c(trait_selection),'species_latin', 'Species')

# Look at PCAs
pca_ins  = dudi.pca(CWM_birds_insect[complete.cases(CWM_birds_insect[,.SD, .SDcols = c(trait_selection)]), .SD, .SDcols = c(trait_selection)], scannf = FALSE, nf = 3)

pca = pca_ins
cwm = CWM_birds_insect
tot_pca = fviz_pca(pca, geom.ind = 'text')
tot_pca
quanti.coord <- supcol(pca, data.frame(scale(env_data_lui[Plot %in% cwm$Plot, c( 'LUI')])))$cosup * pca$eig[1:3]^2
tot_pca12 <-fviz_add(tot_pca, quanti.coord, axes = c(1, 2), "arrow", color = "blue", linetype = "solid", repel = T,
                     addlabel = T
)
tot_pca12
cor.test(pca_herb$l1$RS1, env_data_lui[Plot %in% CWM_birds_herb$Plot,]$LUI)

tot_pca12

### Save data

fwrite(CWM_birds_insect[, list( "Plot" = Plot               ,
                                "Year" = Year,
                                "Bi_Size" = log_body_len     ,
                                "Bi_Incub" = log_incub,
                                "Bi_TOffsprings" = log_offspring  ,
                                "Bi_AgeMax" = log_longevity,
                                "Bi_GenLength" = GenLength)

], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_birds_insect.csv")

# Non-weighted community traits

Abundances_birds_presence_absence = Abundances[, list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Species)]
Abundances_birds_presence_absence[value>1, value := 1]

CWM_birds_insect_noweight  = my_cwm(bird_traits[trophic_level %in% c('insectivore','carnivore'),], Abundances_birds_presence_absence, c(trait_selection),'species_latin', 'Species')

fwrite(CWM_birds_insect_noweight[, list( "Plot" = Plot               ,
                                "Year" = Year,
                                "Bi_Size" = log_body_len     ,
                                "Bi_Incub" = log_incub,
                                "Bi_TOffsprings" = log_offspring  ,
                                "Bi_AgeMax" = log_longevity,
                                "Bi_GenLength" = GenLength)
                        
], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_birds_insect_noweight.csv")


### Check turnover
data_lui <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Environment/LUI_input_data/LUI_standardized_global.txt")
data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

library(betapart)
comm.test = dcast(Abundances[Species %in% bird_traits[trophic_level %in% c('insectivore','carnivore'), species_latin], list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Species)],  Plot~Species, value.var = 'value', fill = 0)
rownames(comm.test)= comm.test$Plot
comm.test = comm.test[,-1]

beta.multi.abund(comm.test)

comm_min_max = matrix(c(colSums(comm.test[min_lui_plots,]),colSums(comm.test[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_min_max)
